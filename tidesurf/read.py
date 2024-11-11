from dataclasses import dataclass
from tidesurf.transcript import Strand


@dataclass
class Read:
    cbc: str
    umi: str
    chromosome: str
    strand: Strand
    start: int
    end: int
    length: int

    def __repr__(self) -> str:
        return (
            f"<Read with CBC = {self.cbc}, UMI = {self.umi} mapped to"
            f" {self.chromosome}{self.strand}:{self.start:,}-"
            f"{self.end:,} at {hex(id(self))}>"
        )
