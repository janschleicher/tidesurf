cpdef enum SpliceType:
    UNSPLICED = 0
    AMBIGUOUS = 1
    SPLICED = 2


cpdef enum ReadType:
    INTRON = 0
    EXON_EXON = 1
    AMBIGUOUS_READ = 2
    EXON = 3