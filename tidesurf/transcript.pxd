from tidesurf.enums cimport Strand

cdef class GenomicFeature:
    cdef readonly str gene_id
    cdef readonly str gene_name
    cdef readonly str transcript_id
    cdef readonly str transcript_name
    cdef readonly str chromosome
    cdef readonly Strand strand
    cdef readonly int start
    cdef readonly int end

    cpdef bint overlaps(
        self,
        str chromosome,
        Strand strand,
        int start,
        int end,
        int min_overlap=*,
    )


cdef class Exon(GenomicFeature):
    cdef readonly exon_id
    cdef readonly exon_number

cdef class Intron(GenomicFeature):
    pass


cdef class Transcript(GenomicFeature):
    cdef readonly list regions

    cpdef void add_exon(self, Exon exon)
    cpdef void sort_regions(self)


cdef class TranscriptIndex:
    cdef readonly dict transcripts
    cdef readonly dict transcripts_by_region

    cpdef void _add_transcripts(
        self,
        dict chrom_transcript_dict,
        list start_end_positions,
        str curr_chrom,
        int curr_strand,
    )
    cpdef void read_gtf(self, str gtf_file)
    cpdef list get_overlapping_transcripts(
        self,
        str chromosome,
        Strand strand,
        int start,
        int end,
    )

cdef int _bisect_sort_key(tuple x)
