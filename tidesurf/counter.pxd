from tidesurf.transcript cimport TranscriptIndex
from pysam.libcalignedsegment cimport AlignedSegment

cdef class UMICounter:
    cdef readonly TranscriptIndex transcript_index
    cdef readonly str orientation
    cdef readonly int MIN_INTRON_OVERLAP
    cdef readonly bint multi_mapped_reads

    cpdef tuple count(
        self,
        str bam_file,
        bint filter_cells,
        str whitelist=*,
        int num_umis=*,
    )
    cpdef tuple _process_read(self, AlignedSegment read)