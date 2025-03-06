"""Module containing enums for tidesurf."""

cimport cython

@cython.embedsignature.format("python")
cpdef Strand antisense(Strand strand):
    """Return the antisense of a given strand."""
    if strand == Strand.PLUS:
        return Strand.MINUS
    elif strand == Strand.MINUS:
        return Strand.PLUS
    else:
        raise ValueError("Invalid strand")
