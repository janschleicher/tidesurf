import cython
from tidesurf.transcript cimport GenomicFeature, TranscriptIndex
from pysam.libcalignedsegment cimport AlignedSegment

@cython.final
cdef class UMICounter:
    cdef readonly TranscriptIndex transcript_index
    """Transcript index for extraction of overlapping transcripts."""
    cdef readonly str orientation
    """Orientation in which reads map to transcripts."""
    cdef readonly int MIN_INTRON_OVERLAP
    """Minimum overlap of reads with introns required to consider them intronic."""
    cdef readonly bint multi_mapped_reads
    """Whether to count multi-mapped reads."""
    
    cpdef tuple count(
        self,
        str bam_file,
        bint filter_cells=*,
        str whitelist=*,
        int num_umis=*,
    )
    cpdef inline tuple _process_read(self, AlignedSegment read)
