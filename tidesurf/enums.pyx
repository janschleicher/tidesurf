"""Module containing enums for tidesurf."""

cpdef Strand antisense(Strand strand):
    """Return the antisense of a given strand."""
    if strand == Strand.PLUS:
        return Strand.MINUS
    elif strand == Strand.MINUS:
        return Strand.PLUS
    else:
        raise ValueError("Invalid strand")
