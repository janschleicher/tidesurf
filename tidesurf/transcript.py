from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional


class GenomicFeature:
    __slots__ = [
        "gene_id",
        "gene_name",
        "transcript_id",
        "transcript_name",
        "chromosome",
        "strand",
        "start",
        "end",
        "exons",
        "introns",
    ]

    def __init__(
        self,
        gene_id: str,
        gene_name: str,
        transcript_id: str,
        transcript_name: str,
        chromosome: str,
        strand: str,
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

    def __lt__(self, other) -> bool:
        assert (
            self.chromosome == other.chromosome
        ), "Cannot compare features on different chromosomes"
        return self.start < other.start

    def __gt__(self, other) -> bool:
        assert (
            self.chromosome == other.chromosome
        ), "Cannot compare features on different chromosomes"
        return self.start > other.start

    def __repr__(self) -> str:
        return f"<GenomicFeature {self.chromosome}:{self.start}-{self.end} on '{self.strand}' strand at {hex(id(self))}>"


class Exon(GenomicFeature):
    __slots__ = ["exon_id", "exon_number"]

    def __init__(
        self,
        gene_id: str,
        gene_name: str,
        transcript_id: str,
        transcript_name: str,
        chromosome: str,
        strand: str,
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

    def __eq__(self, other) -> bool:
        return self.exon_id == other.exon_id

    def __repr__(self) -> str:
        return (
            f"<Exon {self.exon_id}, No. {self.exon_number} for transcript "
            f"{self.transcript_id} {self.chromosome}:{self.start}-"
            f"{self.end} on '{self.strand}' strand at {hex(id(self))}>"
        )


class Transcript(GenomicFeature):
    __slots__ = ["exons"]

    def __init__(
        self,
        gene_id: str,
        gene_name: str,
        transcript_id: str,
        transcript_name: str,
        chromosome: str,
        strand: str,
        start: int,
        end: int,
        exons: Optional[List[Exon]] = None,
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
        if exons is None:
            self.exons = []
        else:
            self.exons = exons

    def add_exon(self, exon: Exon) -> None:
        if exon not in self.exons:
            self.exons.append(exon)

    def sort_exons(self) -> None:
        self.exons.sort()

    def __eq__(self, other) -> bool:
        return self.transcript_id == other.transcript_id

    def __repr__(self) -> str:
        return (
            f"<Transcript {self.transcript_id} {self.chromosome}:{self.start}-"
            f"{self.end} on '{self.strand}' strand at {hex(id(self))}>"
        )


@dataclass
class GTFLine:
    chromosome: str
    source: str
    feature: str
    start: int
    end: int
    score: str
    strand: str
    frame: str
    attributes: Dict[str, str]

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


class TranscriptIndex:
    __slots__ = ["transcripts"]
    transcripts: Dict[Tuple[str, str], List[Transcript]]

    def __init__(self, gtf_file: str) -> None:
        self.transcripts = {}
        self.read_gtf(gtf_file)

    def read_gtf(self, gtf_file: str) -> None:
        lines = []

        # Read the GTF file
        with open(gtf_file, "r") as gtf:
            for line in gtf:
                # Skip header lines and comments
                if line.startswith("#"):
                    continue

                # Parse the GTF line
                (
                    curr_chrom,
                    source,
                    feature,
                    start,
                    end,
                    score,
                    curr_strand,
                    frame,
                    attributes_str,
                ) = line.strip().split("\t")

                # Only keep exons and transcripts
                if feature not in ["exon", "transcript"]:
                    continue

                start, end = int(start), int(end)
                attributes = {
                    key: value.strip('"')
                    for attr in attributes_str.split("; ")
                    for key, value in [attr.split(" ")]
                }
                gtf_line = GTFLine(
                    chromosome=curr_chrom,
                    source=source,
                    feature=feature,
                    start=start,
                    end=end,
                    score=score,
                    strand=curr_strand,
                    frame=frame,
                    attributes=attributes,
                )
                lines.append(gtf_line)
        lines.sort()

        # Construct index from lines
        chrom_transcript_dict = {}
        curr_chrom, curr_strand = None, None
        for line in lines:
            if line.chromosome != curr_chrom or line.strand != curr_strand:
                # Going to a new chromosome-strand pair:
                # add the previous one to transcripts
                if curr_chrom is not None and curr_strand is not None:
                    if (curr_chrom, curr_strand) in self.transcripts.keys():
                        raise ValueError("GTF file was not sorted properly.")
                    self.transcripts[curr_chrom, curr_strand] = sorted(
                        chrom_transcript_dict.values()
                    )
                    for transcript in self.transcripts[curr_chrom, curr_strand]:
                        transcript.sort_exons()
                curr_chrom = line.chromosome
                curr_strand = line.strand
                chrom_transcript_dict = {}
            # Add new transcript to dictionary
            if line.feature == "transcript":
                if line.attributes["transcript_id"] not in chrom_transcript_dict.keys():
                    chrom_transcript_dict[line.attributes["transcript_id"]] = (
                        Transcript(
                            gene_id=line.attributes["gene_id"],
                            gene_name=line.attributes["gene_name"],
                            transcript_id=line.attributes["transcript_id"],
                            transcript_name=line.attributes["transcript_name"],
                            chromosome=line.chromosome,
                            strand=line.strand,
                            start=line.start,
                            end=line.end,
                        )
                    )
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
                chrom_transcript_dict[line.attributes["transcript_id"]].add_exon(exon)

        # Add last chromosome-strand pair
        if curr_chrom is not None and curr_strand is not None:
            if (curr_chrom, curr_strand) in self.transcripts.keys():
                raise ValueError("GTF file was not sorted properly.")
            self.transcripts[curr_chrom, curr_strand] = sorted(
                chrom_transcript_dict.values()
            )
            for transcript in self.transcripts[curr_chrom, curr_strand]:
                transcript.sort_exons()

    def get_overlapping_transcripts(
        self, chromosome: str, start: int, end: int
    ) -> List[Transcript]:
        # TODO: Implement binary search
        pass
