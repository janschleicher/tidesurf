from typing import List


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

    def __repr__(self) -> str:
        return f"<GenomicFeature {self.chromosome}:{self.start}-{self.end} on '{self.strand}' strand at {hex(id(self))}>"


class Exon(GenomicFeature):
    __slots__ = ["exon_number"]

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
        self.exon_number = exon_number

    def __repr__(self) -> str:
        return (
            f"<Exon No. {self.exon_number} for transcript "
            f"{self.transcript_id} {self.chromosome}:{self.start}-"
            f"{self.end} on '{self.strand}' strand at {hex(id(self))}>"
        )


class Intron(GenomicFeature):
    __slots__ = ["intron_number"]

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
        intron_number: int,
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
        self.intron_number = intron_number

    def __repr__(self) -> str:
        return (
            f"<Intron No. {self.intron_number} for transcript "
            f"{self.transcript_id} {self.chromosome}:{self.start}-"
            f"{self.end} on '{self.strand}' strand at {hex(id(self))}>"
        )


class Transcript(GenomicFeature):
    __slots__ = ["exons", "introns"]

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
        exons: List[Exon],
        introns: List[Intron],
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
        self.exons = exons
        self.introns = introns

    def __repr__(self) -> str:
        return (
            f"<Transcript {self.transcript_id} {self.chromosome}:{self.start}-"
            f"{self.end} on '{self.strand}' strand at {hex(id(self))}>"
        )
