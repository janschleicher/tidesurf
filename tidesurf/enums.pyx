cpdef Strand antisense(Strand strand):
    if strand == Strand.PLUS:
        return Strand.MINUS
    elif strand == Strand.MINUS:
        return Strand.PLUS
    else:
        raise ValueError("Invalid strand")
