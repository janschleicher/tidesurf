import numpy as np
import polars as pl
import pysam
from scipy.sparse import csr_matrix, lil_matrix
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
from tidesurf.transcript import TranscriptIndex, Strand
from enum import Enum
from typing import Literal, Tuple, Optional, Dict
import logging

log = logging.getLogger(__name__)


class SpliceType(Enum):
    """
    Enum for read/UMI splice types.
    """

    UNSPLICED = 0
    AMBIGUOUS = 1
    SPLICED = 2

    def __lt__(self, other):
        return self.value < other.value

    def __int__(self):
        return self.value


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
                    continue
                if filter_cells and whitelist:
                    if bam_read.get_tag("CB") not in whitelist:
                        continue
                res = self._process_read(bam_read)
                if res is not None:
                    results.append(res)

        # Deduplicate cell barcodes and UMIs.
        log.info("Deduplicating cell barcodes and UMIs.")
        results = (
            pl.DataFrame(
                results,
                schema={"cbc": str, "umi": str, "gene": str, "splice_type": int},
                strict=False,
                orient="row",
            )
            .sort(by="splice_type")
            .unique(
                subset=[
                    "cbc",
                    "umi",
                    "gene",
                ],
                keep="first",
            )
            # Will keep the first splice type, unspliced before
            # ambiguous before spliced
        )

        # Remove multi-mapped UMIs.
        if not self.multi_mapped:
            results = results.unique(subset=["cbc", "umi"], keep="none")

        # Do the rest of the counting.
        log.info("Counting UMIs.")
        cells = results["cbc"].unique().sort().to_numpy()
        genes = results["gene"].unique().sort().to_numpy()
        counts_dict = _count(
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
    ) -> Optional[Tuple[str, str, str, SpliceType]]:
        """
        Process a single read.

        :param read: The read to process.
        :return: cell barcode, UMI, and gene name.
        """
        cbc = read.get_tag("CB")
        umi = read.get_tag("UB")
        chromosome = (read.reference_name,)  # [0],
        chromosome = chromosome[0]
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
        if not overlapping_transcripts:
            return None

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

        splice_types = set()
        for trans in overlapping_transcripts:
            # Loop over exons
            total_exon_overlap = 0
            for exon in trans.exons:
                total_exon_overlap += read.get_overlap(exon.start, exon.end + 1)

            # Assign splice type for this transcript to spliced if at most 5 bases do not overlap with exons
            if length - total_exon_overlap <= 5:
                splice_types.add(SpliceType.SPLICED)
            else:
                splice_types.add(SpliceType.UNSPLICED)

        # Get gene names.
        gene_names = {t.gene_name for t in overlapping_transcripts}
        if not self.multi_mapped and len(gene_names) > 1:
            return None
        elif len(gene_names) == 1:
            return (
                cbc,
                umi,
                gene_names.pop(),
                int(splice_types.pop())
                if len(splice_types) == 1
                else int(SpliceType.AMBIGUOUS),
            )
        elif self.multi_mapped:
            raise NotImplementedError("Multi-mapped reads not yet implemented.")


def _count(arr, cells, genes):
    cbc_map = {cbc: idx for idx, cbc in enumerate(cells)}
    gene_map = {gene: idx for idx, gene in enumerate(genes)}
    counts_dict = {
        SpliceType.SPLICED: {},
        SpliceType.UNSPLICED: {},
        SpliceType.AMBIGUOUS: {},
    }
    for line in tqdm(arr, desc="Counting s/u/a UMIs", unit=" UMIs"):
        cbc, gene, splice_type = (
            cbc_map[line[0]],
            gene_map[line[2]],
            SpliceType(line[3]),
        )
        if (cbc, gene) not in counts_dict[splice_type]:
            counts_dict[splice_type][cbc, gene] = 0
        counts_dict[splice_type][cbc, gene] += 1
    return counts_dict
