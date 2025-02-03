from importlib.metadata import PackageNotFoundError, version

from tidesurf.counter import UMICounter
from tidesurf.transcript import (
    Exon,
    GenomicFeature,
    Strand,
    Transcript,
    TranscriptIndex,
)

try:
    __version__ = version("tidesurf")
except PackageNotFoundError:
    pass

__all__ = [
    "Strand",
    "GenomicFeature",
    "Exon",
    "Transcript",
    "TranscriptIndex",
    "UMICounter",
]
