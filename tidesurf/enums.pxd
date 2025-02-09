cpdef enum Strand:
    PLUS = 0
    MINUS = 1


cpdef enum ReadType:
    INTRON = 0
    EXON_EXON = 1
    AMBIGUOUS_READ = 2
    EXON = 3


cpdef enum SpliceType:
    UNSPLICED = 0
    AMBIGUOUS = 1
    SPLICED = 2

cpdef Strand antisense(Strand strand)