import numpy as np
import polars as pl
from bisect import bisect
import pysam
from scipy.sparse import csr_matrix, lil_matrix
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
from tidesurf.transcript import TranscriptIndex, Strand
from enum import Enum
from typing import Literal, Tuple, Optional, Dict, List
import logging

log = logging.getLogger(__name__)


class SpliceType(Enum):
    """
    Enum for read/UMI splice types.
    """

    UNSPLICED = 0
    AMBIGUOUS = 1
    SPLICED = 2

    def __int__(self):
        return self.value


class ReadType(Enum):
    """
    Enum for read alignment types.
    """

    INTRON = 0
    EXON_EXON = 1
    AMBIGUOUS = 2
    EXON = 3

    def __int__(self):
        return self.value

    def get_splice_type(self):
        if self == ReadType.INTRON:
            return SpliceType.UNSPLICED
        elif self == ReadType.EXON_EXON or self == ReadType.EXON:
            return SpliceType.SPLICED
        else:
            return SpliceType.AMBIGUOUS


class UMICounter:
    """
    Counter for unique molecular identifiers (UMIs) with reads mapping
    to transcripts in single-cell RNA-seq data.

    :param transcript_index: Transcript index.
    :param orientation: Orientation in which reads map to transcripts.
    Either "sense" or "antisense".
    :param multi_mapped: Whether to count multi-mapped reads.
    """

    def __init__(
        self,
        transcript_index: TranscriptIndex,
        orientation: Literal["sense", "antisense"],
        multi_mapped: bool = False,
    ) -> None:
        self.transcript_index = transcript_index
        self.orientation = orientation
        self.multi_mapped = multi_mapped

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
        only keep cells with at least that many UMIs. Mutually exclusive
        with `whitelist`.
        :return: cells (array of shape (n_cells,)), genes (array of
        shape (n_genes,)), counts (sparse matrix of shape (n_cells, n_genes)).
        """
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
                whitelist = set(
                    pl.read_csv(whitelist, has_header=False)[:, 0].str.strip_chars()
                )

        aln_file = pysam.AlignmentFile(bam_file)
        total_reads = 0
        for idx_stats in aln_file.get_index_statistics():
            total_reads += idx_stats.total

        with logging_redirect_tqdm():
            results = []
            log.info("Processing reads from BAM file.")
            skipped_reads = 0
            for bam_read in tqdm(
                aln_file, total=total_reads, desc="Processing BAM file", unit=" reads"
            ):
                if (
                    bam_read.is_unmapped
                    or bam_read.mapping_quality
                    != 255  # discard reads with mapping quality < 255
                    or not bam_read.has_tag("CB")
                    or not bam_read.has_tag("UB")
                ):
                    skipped_reads += 1
                    continue
                if filter_cells and whitelist:
                    if bam_read.get_tag("CB") not in whitelist:
                        continue
                res = self._process_read(bam_read)
                if res is not None:
                    results.append(res)
                else:
                    skipped_reads += 1
        log.info(f"Skipped {skipped_reads:,} reads.")

        # Deduplicate cell barcodes and UMIs.
        log.info("Deduplicating cell barcodes and UMIs.")
        results = (
            pl.DataFrame(
                results,
                schema={"cbc": str, "umi": str, "gene": str, "read_type": pl.UInt8},
                strict=False,
                orient="row",
            )
            .group_by("cbc", "umi", "gene", "read_type")
            .len(name="count")  # Count ReadTypes per cbc/umi/gene combination
            .select(
                pl.all(), (pl.sum("count").over("cbc", "umi", "gene")).alias("total")
            )
            .select(pl.all(), (pl.col("count") / pl.col("total")).alias("percentage"))
            .filter(  # Remove read types with low counts or percentage
                ~(
                    ((pl.col("count") < 2) & (pl.col("percentage") < 0.1))
                    | (pl.col("percentage") < 0.1)
                )
            )
            .group_by("cbc", "umi", "gene")
            # Keep the first ReadType, order: INTRON, EXON_EXON, AMBIGUOUS, EXON
            .agg(pl.min("read_type"), pl.max("total"))
            .with_columns(
                pl.col("read_type")
                .map_elements(  # Map ReadType to SpliceType
                    lambda x: int(ReadType(x).get_splice_type()), return_dtype=pl.UInt8
                )
                .alias("splice_type")
            )
            .drop("read_type")
        )

        # Resolve multi-mapped UMIs.
        def _argmax(lst: List[int]) -> int:
            _, indices, value_counts = np.unique(
                lst, return_index=True, return_counts=True
            )
            if value_counts[-1] > 1:
                return -1
            else:
                return indices[-1]

        _argmax_vec = np.vectorize(_argmax)

        if not self.multi_mapped:
            # Keep the gene with the highest read support
            results = (
                results.group_by("cbc", "umi")
                .agg(pl.col("gene"), pl.col("total"), pl.min("splice_type"))
                .with_columns(
                    (
                        pl.when(pl.col("total").list.len() > 1)
                        .then(
                            pl.col("total").map_batches(
                                _argmax_vec, return_dtype=pl.Int8
                            )
                        )
                        .otherwise(pl.lit(0, dtype=pl.Int8))
                    ).alias("idx")
                )
                # Ties for maximal read support (represented by -1)
                # are discarded
                .filter(pl.col("idx") >= 0)
                .with_columns(
                    pl.col("gene").list.get(pl.col("idx")),
                )
                .drop("total", "idx")
            )
        else:
            results = results.drop("total")

        # Do the rest of the counting.
        log.info("Counting UMIs.")
        cells = results["cbc"].unique().sort().to_numpy()
        genes = results["gene"].unique().sort().to_numpy()
        counts_dict = self.__count(
            results.to_numpy(),
            cells,
            genes,
        )
        counts = {
            key: lil_matrix((cells.shape[0], genes.shape[0]), dtype=np.int32)
            for key in ["spliced", "unspliced", "ambiguous"]
        }
        for splice_type, counts_dict_splice in counts_dict.items():
            idx = np.asarray(list(counts_dict_splice.keys()))
            counts[splice_type.name.lower()][idx[:, 0], idx[:, 1]] = np.asarray(
                list(counts_dict_splice.values())
            )

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

    def _process_read(
        self, read: pysam.AlignedSegment
    ) -> Optional[Tuple[str, str, str, int]]:
        """
        Process a single read.

        :param read: The read to process.
        :return: cell barcode, UMI, gene name, and read type.
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
        # TODO: Does this make sense?
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

        read_types = set()
        for trans in overlapping_transcripts:
            # Loop over exons
            total_exon_overlap = 0
            n_exons = 0
            left_idx = max(bisect(trans.exons, start, key=lambda x: x.start) - 1, 0)
            for exon in trans.exons[left_idx:]:
                if exon.start > end:
                    break
                exon_overlap = read.get_overlap(exon.start, exon.end + 1)
                total_exon_overlap += exon_overlap
                if exon_overlap > 0:
                    n_exons += 1

            # Assign read alignment region for this transcript to exonic
            # if at most 5 bases do not overlap with exons
            if length - total_exon_overlap <= 5:
                # More than one exon: exon-exon junction
                # TODO: Should I check for Ns in cigar string?
                if n_exons > 1:
                    read_types.add(ReadType.EXON_EXON)
                elif n_exons == 1:
                    read_types.add(ReadType.EXON)
                else:
                    raise ValueError("Exon overlap without exons.")
            # Special case: if read overlaps with only first exon and the
            # region before or with only last exon and the region after
            elif (
                (left_idx == 0 and start < trans.exons[left_idx].start)
                or (
                    left_idx == len(trans.exons) - 1 and end > trans.exons[left_idx].end
                )
            ) and total_exon_overlap == read.get_overlap(
                trans.exons[left_idx].start, trans.exons[left_idx].end + 1
            ):
                read_types.add(ReadType.EXON)
            else:
                read_types.add(ReadType.INTRON)

        # Get gene names.
        gene_names = {t.gene_name for t in overlapping_transcripts}
        if not self.multi_mapped and len(gene_names) > 1:
            return None
        elif len(gene_names) == 1:
            if ReadType.EXON_EXON in read_types:
                read_type = ReadType.EXON_EXON
            elif len(read_types) == 1:
                read_type = read_types.pop()
            else:
                read_type = ReadType.AMBIGUOUS

            return (
                cbc,
                umi,
                gene_names.pop(),
                int(read_type),
            )
        elif self.multi_mapped:
            raise NotImplementedError("Multi-mapped reads not yet implemented.")

    @staticmethod
    def __count(
        arr: np.ndarray, cells: np.ndarray, genes: np.ndarray
    ) -> Dict[SpliceType, Dict[Tuple[int, int], int]]:
        cbc_map = {cbc: idx for idx, cbc in enumerate(cells)}
        gene_map = {gene: idx for idx, gene in enumerate(genes)}
        counts_dict = {
            SpliceType.SPLICED.value: {},
            SpliceType.UNSPLICED.value: {},
            SpliceType.AMBIGUOUS.value: {},
        }
        for row in tqdm(arr, desc="Counting s/u/a UMIs", unit=" UMIs"):
            cbc, gene, splice_type = (
                cbc_map[row[0]],
                gene_map[row[2]],
                row[3],
            )
            if (cbc, gene) not in counts_dict[splice_type]:
                counts_dict[splice_type][cbc, gene] = 0
            counts_dict[splice_type][cbc, gene] += 1
        return {SpliceType(key): val for key, val in counts_dict.items()}
