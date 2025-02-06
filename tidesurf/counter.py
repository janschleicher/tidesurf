# cython: embedsignature=True, annotation_typing=False

"""Module for counting UMIs with reads mapping to transcripts."""

import logging
from bisect import bisect
from itertools import zip_longest
from typing import Dict, List, Literal, Optional, Tuple

import cython
import numpy as np
import polars as pl
import pysam
from scipy.sparse import csr_matrix, lil_matrix
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

from tidesurf.read import ReadType, SpliceType
from tidesurf.transcript import Exon, GenomicFeature, Intron, Strand, TranscriptIndex

log = logging.getLogger(__name__)


@cython.ccall
def _get_splice_type(read_type: int) -> int:
    """Return the corresponding SpliceType for a ReadType."""
    if read_type == ReadType.INTRON:
        return SpliceType.UNSPLICED
    elif read_type == ReadType.EXON_EXON or read_type == ReadType.EXON:
        return SpliceType.SPLICED
    else:
        return SpliceType.AMBIGUOUS


@cython.cclass
class UMICounter:
    """
    Counter for unique molecular identifiers (UMIs) with reads mapping
    to transcripts in single-cell RNA-seq data.

    :param transcript_index: Transcript index.
    :param orientation: Orientation in which reads map to transcripts.
        Either "sense" or "antisense".
    :param min_intron_overlap: Minimum overlap of reads with introns
        required to consider them intronic.
    :param cbc_chunk_size: Number of cells per chunk for demultiplexing
        UMIs. Increasing the chunk size will result in higher memory
        requirements.
    :param multi_mapped_reads: Whether to count multi-mapped reads.
    """

    __slots__ = [
        "transcript_index",
        "orientation",
        "MIN_INTRON_OVERLAP",
        "CBC_CHUNK_SIZE",
        "multi_mapped_reads",
    ]
    transcript_index: TranscriptIndex
    """Transcript index for extraction of overlapping transcripts."""
    orientation = cython.declare(str, visibility="readonly")
    """Orientation in which reads map to transcripts."""
    MIN_INTRON_OVERLAP = cython.declare(cython.int, visibility="readonly")
    """Minimum overlap of reads with introns required to consider them intronic."""
    CBC_CHUNK_SIZE = cython.declare(cython.int, visibility="readonly")
    """Number of cells per chunk for demultiplexing UMIs."""
    multi_mapped_reads = cython.declare(cython.bint, visibility="readonly")
    """Whether to count multi-mapped reads."""

    def __init__(
        self,
        transcript_index: TranscriptIndex,
        orientation: Literal["sense", "antisense"],
        min_intron_overlap: int = 5,
        cbc_chunk_size: int = 10,
        multi_mapped_reads: bool = False,
    ) -> None:
        self.transcript_index = transcript_index
        self.orientation = orientation
        self.MIN_INTRON_OVERLAP = min_intron_overlap
        self.CBC_CHUNK_SIZE = cbc_chunk_size
        self.multi_mapped_reads = multi_mapped_reads

    @cython.ccall
    def count(
        self,
        bam_file: str,
        filter_cells: bool = False,
        whitelist: Optional[str] = None,
        num_umis: Optional[int] = None,
    ) -> Tuple[np.ndarray, np.ndarray, Dict[str, csr_matrix]]:
        """
        Count UMIs with reads mapping to transcripts.

        :param bam_file: Path to BAM file.
        :param filter_cells: Whether to filter cells.
        :param whitelist: If `filter_cells` is True: path to cell
            barcode whitelist file. Mutually exclusive with `num_umis`.
        :param num_umis: If `filter_cells` is True: set to an integer to
            only keep cells with at least that many UMIs. Mutually
            exclusive with `whitelist`.
        :return: cells (array of shape (n_cells,)), genes (array of
            shape (n_genes,)), counts (sparse matrix of shape
            (n_cells, n_genes)).
        """
        wl = {}
        if filter_cells:
            if not whitelist and not num_umis:
                raise ValueError(
                    "Either whitelist or num_umis must be provided when filter_cells==True."
                )
            elif whitelist and num_umis:
                raise ValueError(
                    "Whitelist and num_umis are mutually exclusive arguments."
                )
            elif whitelist:
                log.info(f"Reading whitelist from {whitelist}.")
                wl = set(
                    pl.read_csv(whitelist, has_header=False)[:, 0].str.strip_chars()
                )

        aln_file = pysam.AlignmentFile(bam_file, mode="r")
        total_reads = 0
        for idx_stats in aln_file.get_index_statistics():
            total_reads += idx_stats.total

        with logging_redirect_tqdm():
            results = {}
            log.info("Processing reads from BAM file.")
            skipped_reads = {"unmapped": 0, "no or multimapped transcripts": 0}
            if filter_cells and whitelist:
                skipped_reads["whitelist"] = 0
            for bam_read in tqdm(
                aln_file, total=total_reads, desc="Processing BAM file", unit=" reads"
            ):
                if (
                    bam_read.is_unmapped
                    or bam_read.mapping_quality
                    != 255  # discard reads with mapping quality != 255
                    or not bam_read.has_tag("CB")
                    or not bam_read.has_tag("UB")
                ):
                    skipped_reads["unmapped"] += 1
                    continue
                if filter_cells and whitelist:
                    if bam_read.get_tag("CB") not in wl:
                        skipped_reads["whitelist"] += 1
                        continue
                res = self._process_read(bam_read)
                if res is not None:
                    cbc, results_list = res
                    if cbc in results:
                        results[cbc].extend(results_list)
                    else:
                        results[cbc] = results_list
                else:
                    skipped_reads["no or multimapped transcripts"] += 1
        log.info(
            f"Skipped {', '.join([f'{n_reads:,} reads ({reason})' for reason, n_reads in skipped_reads.items()])}."
        )

        # Deduplicate cell barcodes and UMIs.
        counts_dfs = []
        log.info("Determining splice types and deduplicating UMIs.")
        with logging_redirect_tqdm(), pl.StringCache():
            for chunk in tqdm(
                zip_longest(*[iter(results.items())] * self.CBC_CHUNK_SIZE),
                total=len(results) / self.CBC_CHUNK_SIZE,
                desc="Deduplicating UMIs",
                unit=" CBCs",
                unit_scale=self.CBC_CHUNK_SIZE,
            ):
                df = (
                    pl.concat(
                        [
                            pl.DataFrame(
                                res_lst,
                                schema={
                                    "umi": pl.Categorical,
                                    "gene": pl.String,
                                    "read_type": pl.UInt8,
                                    "weight": pl.Float32,
                                },
                                strict=False,
                                orient="row",
                            ).with_columns(pl.lit(cbc, dtype=pl.String).alias("cbc"))
                            for cbc, res_lst in [
                                item for item in chunk if item is not None
                            ]
                        ]
                    )
                    .group_by("cbc", "umi", "gene", "read_type")
                    .agg(
                        pl.col("weight").sum()
                    )  # Count ReadTypes per umi/gene combination
                    .with_columns(
                        pl.col("read_type")
                        .replace(old=int(ReadType.EXON_EXON), new=int(ReadType.EXON))
                        .alias("read_type_")
                    )
                    .select(
                        pl.exclude("weight"),
                        (pl.sum("weight").over("cbc", "umi", "gene")).alias("total"),
                        (pl.sum("weight").over("cbc", "umi", "gene", "read_type_")),
                    )
                    .select(
                        pl.all(),
                        (pl.col("weight") / pl.col("total")).alias("percentage"),
                    )
                    .filter(  # Remove read types with low counts and percentage (exonic types together)
                        ~(
                            ((pl.col("weight") < 2) & (pl.col("percentage") < 0.1))
                            | (pl.col("percentage") < 0.1)
                        )
                    )
                    .group_by("cbc", "umi", "gene")
                    # Keep the first ReadType, order: INTRON, EXON_EXON, AMBIGUOUS, EXON
                    .agg(pl.min("read_type"), pl.max("total"))
                    # Remove UMIs that are only supported by multimapped reads
                    .filter(pl.col("total") >= 1)
                    .with_columns(
                        pl.col("read_type")
                        .map_elements(  # Map ReadType to SpliceType
                            _get_splice_type,
                            return_dtype=pl.UInt8,
                        )
                        .alias("splice_type")
                    )
                    .drop("read_type")
                )

                # Keep the gene with the highest read support
                df = (
                    df.group_by("cbc", "umi")
                    .agg(pl.col("gene"), pl.col("total"), pl.col("splice_type"))
                    .with_columns(
                        (
                            pl.when(pl.col("total").list.len() > 1)
                            .then(
                                pl.col("total").map_batches(
                                    _argmax_vec, return_dtype=pl.Int16
                                )
                            )
                            .otherwise(pl.lit(0, dtype=pl.Int16))
                        ).alias("idx")
                    )
                    # Ties for maximal read support (represented by -1)
                    # are discarded
                    .filter(pl.col("idx") >= 0)
                    .with_columns(
                        pl.col("gene").list.get(pl.col("idx")),
                        pl.col("splice_type").list.get(pl.col("idx")),
                    )
                    .group_by("cbc", "gene", "splice_type")
                    .len()
                )
                counts_dfs.append(df)

        log.info("Aggregating counts from cell chunks.")
        # Concatenate the cell-wise count DataFrames
        results_df = pl.concat(counts_dfs)

        cells = np.asarray(sorted(results_df["cbc"].unique()))
        genes = np.asarray(sorted(results_df["gene"].unique()))
        n_cells = cells.shape[0]
        n_genes = genes.shape[0]

        # Map cells and genes to integer indices
        cbc_map = {cbc: i for i, cbc in enumerate(cells)}
        gene_map = {gene: i for i, gene in enumerate(genes)}

        results_df = results_df.with_columns(
            pl.col("cbc").replace_strict(cbc_map).name.suffix("_idx"),
            pl.col("gene").replace_strict(gene_map).name.suffix("_idx"),
        )

        assert n_cells == results_df["cbc_idx"].max() + 1
        assert n_genes == results_df["gene_idx"].max() + 1

        # Construct sparse matrices
        counts = {
            key: lil_matrix((n_cells, n_genes), dtype=np.int32)
            for key in [SpliceType.SPLICED, SpliceType.UNSPLICED, SpliceType.AMBIGUOUS]
        }
        for splice_type, mat in counts.items():
            df_ = results_df.filter(pl.col("splice_type") == int(splice_type))
            idx = df_.select("cbc_idx", "gene_idx").to_numpy()
            mat[idx[:, 0], idx[:, 1]] = np.asarray(df_["len"])

        counts = {splice_type.name.lower(): mat for splice_type, mat in counts.items()}

        if filter_cells and num_umis:
            log.info(f"Filtering cells with at least {num_umis} UMIs.")
            idx = (
                counts["spliced"].sum(axis=1).A1
                + counts["unspliced"].sum(axis=1).A1
                + counts["unspliced"].sum(axis=1).A1
            ) >= num_umis
            cells = cells[idx]
            counts = {key: value[idx] for key, value in counts.items()}

        return (
            cells,
            genes,
            {key: csr_matrix(val) for key, val in counts.items()},
        )

    @cython.ccall
    def _process_read(
        self, read: pysam.AlignedSegment
    ) -> Optional[Tuple[str, List[Tuple[str, str, int, float]]]]:
        """
        Process a single read.

        :param read: The read to process.
        :return: cell barcode, list of UMI, gene name, and read type.
        """
        cbc = str(read.get_tag("CB"))
        umi = str(read.get_tag("UB"))
        chromosome = read.reference_name
        strand = Strand("+") if read.is_forward else Strand("-")
        start = read.reference_start
        end = read.reference_end - 1  # pysam reference_end is exclusive
        length = read.infer_read_length()

        if self.orientation == "antisense":
            strand = strand.antisense()

        overlapping_transcripts = self.transcript_index.get_overlapping_transcripts(
            chromosome=chromosome,
            strand=str(strand),
            start=start,
            end=end,
        )

        # Only keep transcripts with minimum overlap of 50% of the read length.
        min_overlap = length // 2
        overlapping_transcripts = [
            t
            for t in overlapping_transcripts
            if t.overlaps(
                chromosome=chromosome,
                strand=str(strand),
                start=start,
                end=end,
                min_overlap=min_overlap,
            )
        ]

        if not overlapping_transcripts:
            return None

        # Determine length of read without soft-clipped bases and count
        # inserted bases (present in read, but not in reference)
        clipped_length = length
        insertion_length = 0
        for cigar_op, n_bases in read.cigartuples:
            if cigar_op == pysam.CSOFT_CLIP:
                clipped_length -= n_bases
            elif cigar_op == pysam.CINS:
                insertion_length += n_bases

        # For each gene, determine the type of read alignment
        read_types_per_gene = {
            trans.gene_name: set() for trans in overlapping_transcripts
        }
        for trans in overlapping_transcripts:
            # Loop over exons and introns
            total_exon_overlap = 0
            total_intron_overlap = 0
            n_exons = 0
            left_idx = max(
                bisect(trans.regions, start, key=_genomic_feature_sort_key) - 1, 0
            )
            for region in trans.regions[left_idx:]:
                if region.start > end:
                    break
                if isinstance(region, Exon):
                    exon_overlap = read.get_overlap(region.start, region.end + 1)
                    total_exon_overlap += exon_overlap
                    if exon_overlap > 0:
                        n_exons += 1
                elif isinstance(region, Intron):
                    total_intron_overlap += read.get_overlap(
                        region.start, region.end + 1
                    )
                else:
                    raise ValueError("Unknown region type.")

            # Assign read alignment region for this transcript to exonic if
            # at most MIN_INTRON_OVERLAP - 1 bases do not overlap with exons
            if (
                clipped_length - total_exon_overlap - insertion_length
                < self.MIN_INTRON_OVERLAP
            ):
                # More than one exon: exon-exon junction
                if n_exons > 1:
                    read_types_per_gene[trans.gene_name].add(ReadType.EXON_EXON)
                elif n_exons == 1:
                    read_types_per_gene[trans.gene_name].add(ReadType.EXON)
                else:
                    raise ValueError("Exon overlap without exons.")
            # Special case: if read overlaps with only first exon and the
            # region before or with only last exon and the region after
            elif (
                left_idx == 0
                and start < trans.regions[left_idx].start
                and end <= trans.regions[left_idx].end
            ) or (
                left_idx == len(trans.regions) - 1
                and end > trans.regions[left_idx].end
                and start >= trans.regions[left_idx].start
            ):
                read_types_per_gene[trans.gene_name].add(ReadType.EXON)
            elif total_intron_overlap >= self.MIN_INTRON_OVERLAP:
                read_types_per_gene[trans.gene_name].add(ReadType.INTRON)

        # Determine ReadType for each mapped gene
        processed_reads = []
        n_genes = len(read_types_per_gene)
        if n_genes > 1 and not self.multi_mapped_reads:
            return None
        for gene_name, read_types in read_types_per_gene.items():
            if not read_types:
                continue
            # Return all genes with their ReadTypes and corresponding weight
            if ReadType.EXON_EXON in read_types:
                read_type = ReadType.EXON_EXON
            elif len(read_types) == 1:
                read_type = read_types.pop()
            else:
                read_type = ReadType.AMBIGUOUS_READ

            processed_reads.append((umi, gene_name, int(read_type), 1.0 / n_genes))
        if not processed_reads:
            return None
        return cbc, processed_reads

    @property
    def transcript_index(self) -> TranscriptIndex:
        return self.transcript_index


@cython.ccall
def _genomic_feature_sort_key(gen_feat: GenomicFeature):
    return gen_feat.start


@cython.ccall
def _argmax(lst: np.ndarray) -> np.int64:
    _, indices, value_counts = np.unique(lst, return_index=True, return_counts=True)
    if value_counts[-1] > 1:
        return -1
    else:
        return indices[-1]


_argmax_vec = np.vectorize(_argmax)
