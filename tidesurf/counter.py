import numpy as np
import pandas as pd
import pysam
from scipy.sparse import csr_matrix
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
from tidesurf.transcript import TranscriptIndex, Strand
from tidesurf.read import Read
from typing import Literal, Tuple, Optional
import logging

log = logging.getLogger(__name__)


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

    def count(self, bam_file: str) -> Tuple[np.ndarray, np.ndarray, csr_matrix]:
        """
        Count UMIs with reads mapping to transcripts.

        :param bam_file: Path to BAM file.
        :return: cells (array of shape (n_cells,)), genes (array of
        shape (n_genes,)), counts (sparse matrix of shape (n_cells, n_genes)).
        """
        aln_file = pysam.AlignmentFile(bam_file)
        total_reads = 0
        for idx_stats in aln_file.get_index_statistics():
            total_reads += idx_stats.total

        with logging_redirect_tqdm():
            results = []
            log.info("Processing reads from BAM file.")
            for bam_read in tqdm(
                aln_file, total=total_reads, desc="Reading BAM file", unit="reads"
            ):
                if (
                    bam_read.is_unmapped
                    or bam_read.mapping_quality
                    != 255  # discard reads with mapping quality < 255
                    or not bam_read.has_tag("CB")
                    or not bam_read.has_tag("UB")
                ):
                    continue
                read = Read(
                    bam_read.get_tag("CB"),
                    bam_read.get_tag("UB"),
                    bam_read.reference_name,
                    Strand("+") if bam_read.is_forward else Strand("-"),
                    bam_read.reference_start,
                    bam_read.reference_end - 1,  # pysam reference_end is exclusive
                    bam_read.infer_read_length(),
                )
                res = self._process_read(read)
                if res is not None:
                    results.append(res)

        # Deduplicate cell barcodes and UMIs.
        log.info("Deduplicating cell barcodes and UMIs.")
        results = pd.DataFrame(results, columns=["cbc", "umi", "gene"]).drop_duplicates(
            keep="first"
        )

        # Remove multi-mapped UMIs.
        if not self.multi_mapped:
            results = results.drop_duplicates(subset=["cbc", "umi"], keep=False)

        # Do the rest of the counting.
        log.info("Counting UMIs.")
        cells = results["cbc"].astype("category").cat.categories.values.astype("<U32")
        genes = results["gene"].astype("category").cat.categories.values.astype("<U32")
        counts_dict = _count(
            results.values.astype("<U32"),
            cells,
            genes,
        )
        idx = np.asarray(list(counts_dict.keys()))
        counts = csr_matrix(
            (np.asarray(list(counts_dict.values())), (idx[:, 0], idx[:, 1]))
        )

        return (
            cells,
            genes,
            counts,
        )

    def _process_read(self, read: Read) -> Optional[Tuple[str, str, str]]:
        """
        Process a single read.

        :param read: The read to process.
        :return: cell barcode, UMI, and gene name.
        """
        strand = read.strand if self.orientation == "sense" else read.strand.antisense()
        overlapping_transcripts = self.transcript_index.get_overlapping_transcripts(
            chromosome=read.chromosome,
            strand=str(strand),
            start=read.start,
            end=read.end,
        )
        if not overlapping_transcripts:
            return None

        # Only keep transcripts with minimum overlap of 50% of the read length.
        # TODO: Does this make sense?
        min_overlap = read.length // 2
        overlapping_transcripts = [
            t
            for t in overlapping_transcripts
            if t.overlaps(
                chromosome=read.chromosome,
                strand=str(strand),
                start=read.start,
                end=read.end,
                min_overlap=min_overlap,
            )
        ]

        # Get gene names.
        gene_names = {t.gene_name for t in overlapping_transcripts}
        if not self.multi_mapped and len(gene_names) > 1:
            return None
        elif len(gene_names) == 1:
            return read.cbc, read.umi, gene_names.pop()
        elif self.multi_mapped:
            raise NotImplementedError("Multi-mapped reads not yet implemented.")


def _count(arr, cells, genes):
    cbc_map = {cbc: idx for idx, cbc in enumerate(cells)}
    gene_map = {gene: idx for idx, gene in enumerate(genes)}
    counts_dict = {}
    for line in arr:
        cbc, gene = cbc_map[line[0]], gene_map[line[2]]
        if (cbc, gene) not in counts_dict:
            counts_dict[cbc, gene] = 0
        counts_dict[cbc, gene] += 1
    return counts_dict
