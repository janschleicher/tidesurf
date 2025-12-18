cpdef enum Strand:
    """Enum for strand information."""
    PLUS = 0
    MINUS = 1


cpdef enum ReadType:
    """Enum for read alignment types."""
    INTRON = 0
    EXON_EXON = 1
    AMBIGUOUS_READ = 2
    EXON = 3


cpdef enum SpliceType:
    """Enum for read/UMI splice types."""
    UNSPLICED = 0
    AMBIGUOUS = 1
    SPLICED = 2

cpdef Strand antisense(Strand strand)