"""Module for working with genomic features and GTF files."""

import logging
from bisect import bisect
from typing import Dict, List, Optional, Set, Tuple, Union

import cython
from cython.cimports.tidesurf.enums import Strand
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

log = logging.getLogger(__name__)

ALLOWED_BIOTYPES = {
    "protein_coding",
    "lncRNA",
    "antisense",
    "IG_C_gene",
    "IG_D_gene",
    "IG_J_gene",
    "IG_LV_gene",
    "IG_V_gene",
    "IG_V_pseudogene",
    "IG_J_pseudogene",
    "IG_C_pseudogene",
    "TR_C_gene",
    "TR_D_gene",
    "TR_J_gene",
    "TR_V_gene",
    "TR_V_pseudogene",
    "TR_J_pseudogene",
}


@cython.cclass
class GenomicFeature:
    """
    A genomic feature on a particular strand on a chromosome. Identified
    by a gene ID, gene name, transcript ID, and transcript name.

    Parameters
    ----------
    gene_id: str
        ID of the corresponding gene.
    gene_name: str
        Name of the corresponding gene.
    transcript_id: str
        ID of the corresponding transcript.
    transcript_name: str
        Name of the corresponding transcript.
    chromosome: str
        Chromosome on which the feature is located.
    strand: Strand
        Strand on which the feature is located.
    start: int
        Genomic start position of the feature (0-based).
    end: int
        Genomic end position of the feature (0-based
    """

    def __init__(
        self,
        gene_id: str,
        gene_name: str,
        transcript_id: str,
        transcript_name: str,
        chromosome: str,
        strand: Strand,
        start: int,
        end: int,
    ) -> None:
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.transcript_id = transcript_id
        self.transcript_name = transcript_name
        self.chromosome = chromosome
        self.strand = strand
        self.start = start
        self.end = end

    @cython.embedsignature(False)
    def overlaps(
        self,
        chromosome: str,
        strand: Strand,
        start: int,
        end: int,
        min_overlap: int = 1,
    ) -> bool:
        """
        overlaps(chromosome: str, strand: Strand, start: int, end: int, min_overlap: int = 1) -> bool

        Check if the feature overlaps with a given region.

        Parameters
        ----------
        chromosome
            Chromosome of interest.
        strand
            Strand of interest.
        start
            Genomic start position of region.
        end
            Genomic end position of region.
        min_overlap
            Minimum number of overlapping bases (default: `1`).

        Returns
        -------
        bool
            Whether the feature overlaps with the region by at least
            `min_overlap` bases.
        """
        if self.chromosome != chromosome or self.strand != strand:
            return False
        assert start <= end, "Start position must be less than or equal to end position"
        return (
            self.start <= end
            and self.end >= start
            and min(self.end - start + 1, end - self.start + 1) >= min_overlap
        )

    def __lt__(self, other) -> bool:
        if self.chromosome != other.chromosome or self.strand != other.strand:
            raise ValueError("Cannot compare features on different chromosomes/strands")
        if self.start != other.start:
            return self.start < other.start
        else:
            return self.end < other.end

    def __gt__(self, other) -> bool:
        if self.chromosome != other.chromosome or self.strand != other.strand:
            raise ValueError("Cannot compare features on different chromosomes/strands")
        if self.start != other.start:
            return self.start > other.start
        else:
            return self.end > other.end

    def __eq__(self, other) -> bool:
        return (
            self.gene_id == other.gene_id
            and self.transcript_id == other.transcript_id
            and self.chromosome == other.chromosome
            and self.strand == other.strand
            and self.start == other.start
            and self.end == other.end
        )

    def __hash__(self) -> int:
        return hash(
            (
                self.gene_id,
                self.transcript_id,
                self.chromosome,
                self.strand,
                self.start,
                self.end,
            )
        )

    def __repr__(self) -> str:
        return (
            f"<GenomicFeature {self.chromosome}:{self.start:,}-"
            f"{self.end:,} on '{self.strand}' strand at {hex(id(self))}>"
        )


@cython.cclass
class Exon(GenomicFeature):
    """
    An exon of a transcript. Identified by an exon ID and exon number.

    Parameters
    ----------
    gene_id: str
        ID of the corresponding gene.
    gene_name: str
        Name of the corresponding gene.
    transcript_id: str
        ID of the corresponding transcript.
    transcript_name: str
        Name of the corresponding transcript.
    chromosome: str
        Chromosome on which the exon is located.
    strand: Strand
        Strand on which the exon is located.
    start: int
        Genomic start position of the exon (0-based).
    end: int
        Genomic end position of the exon (0-based).
    exon_id: str
        ID of the exon.
    exon_number: int
        Number of the exon in the transcript.
    """

    def __init__(
        self,
        gene_id: str,
        gene_name: str,
        transcript_id: str,
        transcript_name: str,
        chromosome: str,
        strand: Strand,
        start: int,
        end: int,
        exon_id: str,
        exon_number: int,
    ) -> None:
        super(Exon, self).__init__(
            gene_id=gene_id,
            gene_name=gene_name,
            transcript_id=transcript_id,
            transcript_name=transcript_name,
            chromosome=chromosome,
            strand=strand,
            start=start,
            end=end,
        )
        self.exon_id = exon_id
        self.exon_number = exon_number

    def __repr__(self) -> str:
        return (
            f"<Exon {self.exon_id}, No. {self.exon_number} for transcript "
            f"{self.transcript_id} {self.chromosome}:{self.start:,}-"
            f"{self.end:,} on '{self.strand}' strand at {hex(id(self))}>"
        )


@cython.cclass
class Intron(GenomicFeature):
    """
    An intron of a transcript.

    Parameters
    ----------
    gene_id: str
        ID of the corresponding gene.
    gene_name: str
        Name of the corresponding gene.
    transcript_id: str
        ID of the corresponding transcript.
    transcript_name: str
        Name of the corresponding transcript.
    chromosome: str
        Chromosome on which the exon is located.
    strand: Strand
        Strand on which the exon is located.
    start: int
        Genomic start position of the exon (0-based).
    end: int
        Genomic end position of the exon (0-based).
    """

    def __init__(
        self,
        gene_id: str,
        gene_name: str,
        transcript_id: str,
        transcript_name: str,
        chromosome: str,
        strand: Strand,
        start: int,
        end: int,
    ) -> None:
        super(Intron, self).__init__(
            gene_id=gene_id,
            gene_name=gene_name,
            transcript_id=transcript_id,
            transcript_name=transcript_name,
            chromosome=chromosome,
            strand=strand,
            start=start,
            end=end,
        )


@cython.cclass
class Transcript(GenomicFeature):
    """
    A transcript. Contains a list of exons and introns.

    Parameters
    ----------
    gene_id: str
        ID of the corresponding gene.
    gene_name: str
        Name of the corresponding gene.
    transcript_id: str
        ID of the corresponding transcript.
    transcript_name: str
        Name of the corresponding transcript.
    chromosome: str
        Chromosome on which the exon is located.
    strand: Strand
        Strand on which the exon is located.
    start: int
        Genomic start position of the exon (0-based).
    end: int
        Genomic end position of the exon (0-based).
    regions: List[Union[Exon, Intron]]
        List of exons and introns in the transcript. A
        :class:`Transcript` object can be initialized without regions
        (default: `None`), in which case they should be added later. If
        only exons are added, introns can be inserted with
        :meth:`~tidesurf.transcript.Transcript.sort_regions`.
    """

    def __init__(
        self,
        gene_id: str,
        gene_name: str,
        transcript_id: str,
        transcript_name: str,
        chromosome: str,
        strand: Strand,
        start: int,
        end: int,
        regions: Optional[List[Union[Exon, Intron]]] = None,
    ) -> None:
        super(Transcript, self).__init__(
            gene_id=gene_id,
            gene_name=gene_name,
            transcript_id=transcript_id,
            transcript_name=transcript_name,
            chromosome=chromosome,
            strand=strand,
            start=start,
            end=end,
        )
        if regions is None:
            self.regions = []
        else:
            self.regions = regions

    @cython.embedsignature(False)
    def add_exon(self, exon: Exon):
        """
        add_exon(exon: Exon)

        Add an exon to the transcript.

        Parameters
        ----------
        exon:
            Exon to add.
        """
        if exon not in self.regions:
            self.regions.append(exon)

    @cython.embedsignature(False)
    def sort_regions(self):
        """
        sort_regions()

        Sort regions by start position and insert introns.
        """
        self.regions = sorted(set(self.regions))
        all_regions = []
        for i, exon in enumerate(self.regions[:-1]):
            # Don't have to check the first region as it is always an exon
            if isinstance(self.regions[i + 1], Intron):
                log.warning("Intron found in regions. Skipping intron insertion.")
                return
            all_regions.append(exon)
            all_regions.append(
                Intron(
                    self.gene_id,
                    self.gene_name,
                    self.transcript_id,
                    self.transcript_name,
                    self.chromosome,
                    self.strand,
                    exon.end + 1,
                    self.regions[i + 1].start - 1,
                )
            )
        all_regions.append(self.regions[-1])
        self.regions = all_regions

    def __eq__(self, other) -> bool:
        return super.__eq__(self, other) and self.regions == other.regions

    def __repr__(self) -> str:
        return (
            f"<Transcript {self.transcript_id} {self.chromosome}:{self.start:,}-"
            f"{self.end:,} on '{self.strand}' strand at {hex(id(self)) } "
            f"containing {len(self.regions)} exons/introns>"
        )

    def __hash__(self) -> int:
        return hash(self.transcript_id)


@cython.cclass
class GTFLine:
    """
    A line from a GTF file, corresponding to particular genomic feature.

    Parameters
    ----------
    chromosome: str
        Chromosome of the feature.
    source: str
        Source of the feature.
    feature: str
        Type of feature.
    start: int
        Genomic start position of feature (0-based).
    end: int
        Genomic end position of feature (0-based).
    score: str
        Feature score.
    strand: Strand
        Strand of the feature.
    frame: str
        Frame of the feature.
    attributes: Dict[str, str]
        Additional attributes of the feature.
    """

    def __init__(
        self,
        chromosome: str,
        source: str,
        feature: str,
        start: int,
        end: int,
        score: str,
        strand: Strand,
        frame: str,
        attributes: Dict[str, str],
    ) -> None:
        self.chromosome = chromosome
        self.source = source
        self.feature = feature
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.frame = frame
        self.attributes = attributes

    def __lt__(self, other) -> bool:
        if self.chromosome != other.chromosome:
            return self.chromosome < other.chromosome
        elif self.strand != other.strand:
            return self.strand < other.strand
        elif self.start != other.start:
            return self.start < other.start
        # Make sure that transcripts come before exons
        else:
            feature_order = {"transcript": 0, "exon": 1}
            return feature_order[self.feature] < feature_order[other.feature]

    def __gt__(self, other) -> bool:
        if self.chromosome != other.chromosome:
            return self.chromosome > other.chromosome
        elif self.strand != other.strand:
            return self.strand > other.strand
        elif self.start != other.start:
            return self.start > other.start
        # Make sure that transcripts come before exons
        else:
            feature_order = {"transcript": 0, "exon": 1}
            return feature_order[self.feature] > feature_order[other.feature]


@cython.cclass
class TranscriptIndex:
    """
    An index of transcripts from a GTF file. Allows for quick retrieval
    of transcripts on a particular chromosome and strand.

    Parameters
    ----------
    gtf_file: str
        Path to GTF file.
    """

    def __init__(self, gtf_file: str) -> None:
        self.transcripts = {}
        self.transcripts_by_region = {}
        self.read_gtf(gtf_file)

    def _add_transcripts(
        self,
        chrom_transcript_dict: Dict[str, Transcript],
        start_end_positions: List[Tuple[int, int, Transcript]],
        curr_chrom: str,
        curr_strand: int,
    ):
        self.transcripts.update(chrom_transcript_dict)
        start_end_positions.sort()
        regions = [(0, set())]
        for pos, is_end, trans in start_end_positions:
            if is_end == 0:
                # Multiple transcripts starting at the same position
                if regions[-1][0] == pos:
                    regions[-1] = (pos, regions[-1][1] | {trans})
                else:
                    regions.append((pos, regions[-1][1] | {trans}))
            else:
                # Multiple transcripts ending at the same position
                if regions[-1][0] == pos + 1:
                    regions[-1] = (pos + 1, regions[-1][1] - {trans})
                else:
                    regions.append((pos + 1, regions[-1][1] - {trans}))
        self.transcripts_by_region[curr_chrom, curr_strand] = regions

    @cython.embedsignature(False)
    def read_gtf(self, gtf_file: str):
        """
        read_gtf(gtf_file: str)

        Read a GTF file and construct an index of transcripts.

        Parameters
        ----------
        gtf_file
            Path to GTF file.
        """
        lines = []

        # Read the GTF file
        with logging_redirect_tqdm(), open(gtf_file, "r") as gtf:
            for line in tqdm(gtf, desc="Reading GTF file", unit=" lines"):
                # Skip header lines and comments
                if line.startswith("#"):
                    continue

                # Parse the GTF line
                (
                    curr_chrom,
                    source,
                    feature,
                    start_str,
                    end_str,
                    score,
                    curr_strand_str,
                    frame,
                    attributes_str,
                ) = line.strip().split("\t")

                curr_strand = Strand.PLUS if curr_strand_str == "+" else Strand.MINUS

                # Only keep exons and transcripts
                if feature not in ["exon", "transcript"]:
                    continue

                start, end = int(start_str), int(end_str)
                attributes = {
                    key: value.strip('"')
                    for attr in attributes_str.split("; ")
                    for key, value in [attr.split(" ")]
                }
                # Don't include transcripts with non-allowed biotypes
                # such as nonsense_mediated_decay
                if attributes["transcript_type"] not in ALLOWED_BIOTYPES:
                    continue
                gtf_line = GTFLine(
                    chromosome="chrM" if curr_chrom == "chrMT" else curr_chrom,
                    source=source,
                    feature=feature,
                    start=start - 1,  # Convert to 0-based
                    end=end - 1,  # Convert to 0-based
                    score=score,
                    strand=curr_strand,
                    frame=frame,
                    attributes=attributes,
                )
                lines.append(gtf_line)
        lines.sort()

        # Construct index from lines
        chrom_transcript_dict = {}
        start_end_positions = []
        curr_chrom, curr_strand = None, None
        with logging_redirect_tqdm():
            for line in tqdm(
                lines, desc="Adding transcripts to index", unit=" GTF lines"
            ):
                if line.chromosome != curr_chrom or line.strand != curr_strand:
                    # Going to a new chromosome-strand pair:
                    # add the previous one to transcripts_by_regions
                    if curr_chrom is not None and curr_strand is not None:
                        if (
                            curr_chrom,
                            curr_strand,
                        ) in self.transcripts_by_region.keys():
                            raise ValueError("GTF file was not sorted properly.")
                        self._add_transcripts(
                            chrom_transcript_dict=chrom_transcript_dict,
                            start_end_positions=start_end_positions,
                            curr_chrom=curr_chrom,
                            curr_strand=curr_strand,
                        )
                    curr_chrom = line.chromosome
                    curr_strand = line.strand
                    start_end_positions = []
                    chrom_transcript_dict = {}
                # Add new transcript to dictionary
                if line.feature == "transcript":
                    if (
                        line.attributes["transcript_id"]
                        not in chrom_transcript_dict.keys()
                    ):
                        transcript = Transcript(
                            gene_id=line.attributes["gene_id"],
                            gene_name=line.attributes["gene_name"],
                            transcript_id=line.attributes["transcript_id"],
                            transcript_name=line.attributes["transcript_name"],
                            chromosome=line.chromosome,
                            strand=line.strand,
                            start=line.start,
                            end=line.end,
                        )
                        chrom_transcript_dict[line.attributes["transcript_id"]] = (
                            transcript
                        )
                        start_end_positions.append((transcript.start, 0, transcript))
                        start_end_positions.append((transcript.end, 1, transcript))
                # Add new exon to corresponding transcript
                elif line.feature == "exon":
                    exon = Exon(
                        gene_id=line.attributes["gene_id"],
                        gene_name=line.attributes["gene_name"],
                        transcript_id=line.attributes["transcript_id"],
                        transcript_name=line.attributes["transcript_name"],
                        chromosome=line.chromosome,
                        strand=line.strand,
                        start=line.start,
                        end=line.end,
                        exon_id=line.attributes["exon_id"],
                        exon_number=int(line.attributes["exon_number"]),
                    )
                    chrom_transcript_dict[line.attributes["transcript_id"]].add_exon(
                        exon
                    )

        # Add last chromosome-strand pair
        if curr_chrom is not None and curr_strand is not None:
            if (curr_chrom, curr_strand) in self.transcripts_by_region.keys():
                raise ValueError("GTF file was not sorted properly.")
            self._add_transcripts(
                chrom_transcript_dict=chrom_transcript_dict,
                start_end_positions=start_end_positions,
                curr_chrom=curr_chrom,
                curr_strand=curr_strand,
            )

        # Sort exons in all transcripts
        for transcript in self.transcripts.values():
            transcript.sort_regions()

    def get_transcript(self, transcript_id: str) -> Optional[Transcript]:
        """
        Get a transcript by its ID.

        Parameters
        ----------
        transcript_id
            Transcript ID.

        Returns
        -------
        Optional[Transcript]
            Transcript object with the given ID if it is in the index.
        """
        if transcript_id in self.transcripts.keys():
            return self.transcripts[transcript_id]
        return None

    @cython.embedsignature(False)
    def get_overlapping_transcripts(
        self,
        chromosome: str,
        strand: Strand,
        start: int,
        end: int,
    ) -> List[Transcript]:
        """
        get_overlapping_transcripts(chromosome: str, strand: Strand, start: int, end: int) -> List[Transcript]

        Get transcripts that overlap with a given region.

        Parameters
        ----------
        chromosome
            Chromosome of interest.
        strand
            Strand of interest.
        start
            Genomic start position of the region.
        end
            Genomic end position of the region.

        Returns
        -------
        List[Transcript]
            List of transcripts that overlap with the region.
        """
        assert (
            start <= end
        ), "Start position must be less than or equal to end position."
        # strand = Strand(strand)
        if (chromosome, strand) not in self.transcripts_by_region.keys():
            return []
        overlapping_transcripts = set()

        # Find region of query start position
        left_idx = (
            bisect(
                self.transcripts_by_region[chromosome, strand],
                start,
                key=_bisect_sort_key,
            )
            - 1
        )
        # Find region of query end position
        right_idx = (
            bisect(
                self.transcripts_by_region[chromosome, strand],
                end,
                key=_bisect_sort_key,
            )
            - 1
        )

        # Get all overlapping transcripts
        if left_idx == right_idx:
            overlapping_transcripts.update(
                self.transcripts_by_region[chromosome, strand][left_idx][1]
            )
        else:
            for i in range(left_idx, right_idx + 1):
                overlapping_transcripts.update(
                    self.transcripts_by_region[chromosome, strand][i][1]
                )

        return sorted(overlapping_transcripts)


@cython.cfunc
@cython.inline
def _bisect_sort_key(x: Tuple[int, Set[Transcript]]) -> int:
    return x[0]
