import numpy as np
import pandas as pd
import pysam
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
from joblib import Parallel, delayed
from tidesurf import TranscriptIndex, Strand
from typing import Literal, Tuple, Optional


class UMICounter:
    """
    Counter for unique molecular identifiers (UMIs) with reads mapping
    to transcripts in single-cell RNA-seq data.

    :param transcript_index: Transcript index.
    :param orientation: Orientation in which reads map to transcripts.
    Either "sense" or "antisense".
    :param multi_mapped: Whether to count multi-mapped reads.
    :param threads: Number of threads to use.
    """

    def __init__(
        self,
        transcript_index: TranscriptIndex,
        orientation: Literal["sense", "antisense"],
        multi_mapped: bool = False,
        threads: int = 1,
    ) -> None:
        self.transcript_index = transcript_index
        self.orientation = orientation
        self.multi_mapped = multi_mapped
        self.threads = threads

    def count(self, bam_file: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Count UMIs with reads mapping to transcripts.

        :param bam_file: Path to BAM file.
        :return: cells (array of shape (n_cells,)), genes (array of
        shape (n_genes,)), counts (array of shape (n_cells, n_genes)).
        """
        aln_file = pysam.AlignmentFile(bam_file)
        total_reads = 0
        for idx_stats in aln_file.get_index_statistics():
            total_reads += idx_stats.total

        with logging_redirect_tqdm():
            results = list(
                tqdm(
                    Parallel(
                        n_jobs=self.threads, return_as="generator", backend="threading"
                    )(delayed(self._process_read)(read) for read in aln_file),
                    total=total_reads,
                    desc="Processing reads",
                    unit="reads",
                )
            )
        # Deduplicate cell barcodes and UMIs.
        results = (
            pd.DataFrame(results, columns=["cbc", "umi", "gene"])
            .dropna()
            .drop_duplicates(keep="first")
        )

        # Remove multi-mapped UMIs.
        if not self.multi_mapped:
            results = results.drop_duplicates(subset=["cbc", "umi"], keep=False)

        # Do the rest of the counting.
        results = (
            results.astype("category")
            .groupby(["cbc", "gene"], observed=False)
            .size()
            .unstack()
        )

        return (
            results.index.values.astype(str),
            results.columns.values.astype(str),
            results.values.astype(int),
        )

    def _process_read(
        self, read: pysam.AlignedSegment
    ) -> Optional[Tuple[str, str, str]]:
        """
        Process a single read.

        :param read: The read to process.
        :return: cell barcode, UMI, and gene name.
        """
        # Get cell barcode and UMI.
        if (
            read.is_unmapped
            or read.mapping_quality != 255  # discard reads with mapping quality < 255
            or read.get_tag("NH") != 1  # discard multimapped reads
            or not read.has_tag("CB")
            or not read.has_tag("UB")
        ):
            return None
        cbc = read.get_tag("CB")
        umi = read.get_tag("UB")

        # Find overlapping transcripts in the right orientation.
        strand = Strand("+") if read.is_forward else Strand("-")
        if self.orientation == "antisense":
            strand = strand.antisense()
        chromosome = read.reference_name
        start = read.reference_start
        end = read.reference_end - 1  # pysam reference_end is exclusive
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
        min_overlap = read.infer_read_length() // 2
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

        # Get gene names.
        gene_names = {t.gene_name for t in overlapping_transcripts}
        if not self.multi_mapped and len(gene_names) > 1:
            return None
        elif len(gene_names) == 1:
            return cbc, umi, gene_names.pop()
        elif self.multi_mapped:
            raise NotImplementedError("Multi-mapped reads not yet implemented.")
