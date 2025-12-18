from importlib.metadata import PackageNotFoundError, version

from tidesurf.counter import UMICounter
from tidesurf.enums import ReadType, SpliceType, Strand
from tidesurf.transcript import (
    Exon,
    GenomicFeature,
    Intron,
    Transcript,
    TranscriptIndex,
)

try:
    __version__ = version("tidesurf")
except PackageNotFoundError:
    pass

__all__ = [
    "Exon",
    "GenomicFeature",
    "Intron",
    "ReadType",
    "SpliceType",
    "Strand",
    "Transcript",
    "TranscriptIndex",
    "UMICounter",
]
