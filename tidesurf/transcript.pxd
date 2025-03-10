from tidesurf.enums cimport Strand

cdef class GenomicFeature:
    cdef readonly str gene_id
    """ID of the corresponding gene."""
    cdef readonly str gene_name
    """Name of the corresponding gene."""
    cdef readonly str transcript_id
    """ID of the corresponding transcript."""
    cdef readonly str transcript_name
    """Name of the corresponding transcript."""
    cdef readonly str chromosome
    """Chromosome on which the feature is located."""
    cdef readonly Strand strand
    """Strand on which the feature is located."""
    cdef readonly int start
    """Genomic start position of the feature (0-based)."""
    cdef readonly int end
    """Genomic end position of the feature (0-based)."""

    cpdef bint overlaps(
        self,
        str chromosome,
        Strand strand,
        int start,
        int end,
        int min_overlap=*,
    )


cdef class Exon(GenomicFeature):
    cdef readonly str exon_id
    """ID of the exon."""
    cdef readonly int exon_number
    """Number of the exon in the transcript."""

cdef class Intron(GenomicFeature):
    pass


cdef class Transcript(GenomicFeature):
    cdef readonly list regions
    """List of exons and introns in the transcript."""

    cpdef void add_exon(self, Exon exon)
    cpdef void sort_regions(self)


cdef class GTFLine:
    cdef readonly str chromosome
    """Chromosome of the feature."""
    cdef readonly str source
    """Source of the feature."""
    cdef readonly str feature
    """Type of feature."""
    cdef readonly int start
    """Genomic start position of the feature (0-based)."""
    cdef readonly int end
    """Genomic end position of the feature (0-based)."""
    cdef readonly str score
    """Feature score."""
    cdef readonly Strand strand
    """Strand of the feature."""
    cdef readonly str frame
    """Frame of the feature."""
    cdef readonly dict attributes
    """Additional attributes of the feature."""



cdef class TranscriptIndex:
    cdef readonly dict transcripts
    """Dictionary of transcripts by ID."""
    cdef readonly dict transcripts_by_region
    """Dictionary of transcript intervals by chromosome and strand."""

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
