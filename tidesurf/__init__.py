from importlib.metadata import version, PackageNotFoundError
from tidesurf.transcript import (
    Strand,
    GenomicFeature,
    Exon,
    Transcript,
    TranscriptIndex,
)
from tidesurf.counter import UMICounter

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
