import pytest

from tidesurf.enums import Strand, antisense
from tidesurf.transcript import (
    Exon,
    GenomicFeature,
    Intron,
    Transcript,
    TranscriptIndex,
)

TEST_GTF_FILE = "test_data/genes.gtf"


def test_strand():
    strand_1 = Strand.PLUS
    assert antisense(strand_1) == Strand.MINUS

    strand_2 = Strand.MINUS
    assert antisense(strand_2) == Strand.PLUS

    assert strand_1 < strand_2
    assert strand_2 > strand_1


def test_genomic_feature():
    gen_feat_1 = GenomicFeature(
        "",
        "",
        "",
        "",
        "chr1",
        Strand.PLUS,
        9_000_000,
        9_001_000,
    )
    gen_feat_2 = GenomicFeature(
        "",
        "",
        "",
        "",
        "chr1",
        Strand.PLUS,
        9_000_500,
        9_001_000,
    )
    gen_feat_3 = GenomicFeature(
        "",
        "",
        "",
        "",
        "chr1",
        Strand.PLUS,
        9_000_000,
        9_001_500,
    )
    gen_feat_4 = GenomicFeature(
        "",
        "",
        "",
        "",
        "chr2",
        Strand.PLUS,
        9_000_000,
        9_001_500,
    )
    assert gen_feat_1 < gen_feat_2
    assert gen_feat_1 < gen_feat_3
    assert gen_feat_3 < gen_feat_2

    assert gen_feat_2 > gen_feat_1
    assert gen_feat_3 > gen_feat_1
    assert gen_feat_2 > gen_feat_3

    with pytest.raises(ValueError):
        gen_feat_1 < gen_feat_4
    with pytest.raises(ValueError):
        gen_feat_4 > gen_feat_1

    assert str(gen_feat_1).startswith(
        "<GenomicFeature chr1:9,000,000-9,001,000 on '0' strand at 0x"
    )

    start, end = 8_999_900, 9_000_100
    assert gen_feat_1.overlaps(
        "chr1", Strand.PLUS, start, end
    ), f"Genomic feature {gen_feat_1} should overlap region chr1+ {start:,}-{end:,}."
    assert not gen_feat_2.overlaps(
        "chr1", Strand.PLUS, start, end
    ), f"Genomic feature {gen_feat_2} should not overlap region chr1+ {start:,}-{end:,}."
    assert not gen_feat_1.overlaps(
        "chr2", Strand.PLUS, start, end
    ), f"Genomic feature {gen_feat_1} should not overlap region chr2+ {start:,}-{end:,}."
    assert not gen_feat_1.overlaps(
        "chr1", Strand.MINUS, start, end
    ), f"Genomic feature {gen_feat_1} should not overlap region chr1- {start:,}-{end:,}."

    # Vary the min_overlap parameter
    end = 9_000_000
    chrom, strand = "chr1", Strand.PLUS
    assert gen_feat_1.overlaps(
        chrom, strand, start, end
    ), f"Genomic feature {gen_feat_1} should overlap region {chrom}{strand} {start:,}-{end:,} by >= 1 base."
    assert not gen_feat_1.overlaps(
        chrom, strand, start, end, min_overlap=2
    ), f"Genomic feature {gen_feat_1} should not overlap region {chrom}{strand} {start:,}-{end:,} by >= 2 bases."
    end = 9_000_001
    assert gen_feat_1.overlaps(
        chrom, strand, start, end, min_overlap=2
    ), f"Genomic feature {gen_feat_1} should overlap region {chrom}{strand} {start:,}-{end:,} by >= 2 bases."
    end = 9_000_008
    assert gen_feat_1.overlaps(
        chrom, strand, start, end, min_overlap=2
    ), f"Genomic feature {gen_feat_1} should overlap region {chrom}{strand} {start:,}-{end:,} by >= 2 bases."
    assert not gen_feat_1.overlaps(
        chrom, strand, start, end, min_overlap=10
    ), f"Genomic feature {gen_feat_1} should not overlap region {chrom}{strand} {start:,}-{end:,} by >= 10 bases."
    assert gen_feat_1.overlaps(
        chrom, strand, start, 9_000_009, min_overlap=10
    ), f"Genomic feature {gen_feat_1} should overlap region {chrom}{strand} {start:,}-9,000,009 by >= 10 bases."

    start, end = 9_001_000, 9_001_100
    assert gen_feat_1.overlaps(
        chrom, strand, start, end
    ), f"Genomic feature {gen_feat_1} should overlap region {chrom}{strand} {start:,}-{end:,} by >= 1 base."
    assert not gen_feat_1.overlaps(
        chrom, strand, start, end, min_overlap=2
    ), f"Genomic feature {gen_feat_1} should not overlap region {chrom}{strand} {start:,}-{end:,} by >= 2 bases."
    start = 9_000_999
    assert gen_feat_1.overlaps(
        chrom, strand, start, end, min_overlap=2
    ), f"Genomic feature {gen_feat_1} should overlap region {chrom}{strand} {start:,}-{end:,} by >= 2 bases."
    start = 9_000_992
    assert gen_feat_1.overlaps(
        chrom, strand, start, end, min_overlap=2
    ), f"Genomic feature {gen_feat_1} should overlap region {chrom}{strand} {start:,}-{end:,} by >= 2 bases."
    assert not gen_feat_1.overlaps(
        chrom, strand, start, end, min_overlap=10
    ), f"Genomic feature {gen_feat_1} should not overlap region {chrom}{strand} {start:,}-{end:,} by >= 10 bases."
    assert gen_feat_1.overlaps(
        chrom, strand, 9_000_991, end, min_overlap=10
    ), f"Genomic feature {gen_feat_1} should overlap region {chrom}{strand} 9,000,991-{end:,} by >= 10 bases."


def test_transcript():
    exon_1 = Exon(
        "ENSMUSG00000025911",
        "Adhfe1",
        "ENSMUST00000144177",
        "Adhfe1-203",
        "chr1",
        Strand.PLUS,
        9_547_948,
        9_548_151,
        "ENSMUSE00000754750",
        1,
    )
    assert str(exon_1).startswith(
        "<Exon ENSMUSE00000754750, No. 1 for transcript "
        "ENSMUST00000144177 chr1:9,547,948-9,548,151 on '0' strand at 0x"
    )
    exon_2 = Exon(
        "ENSMUSG00000025911",
        "Adhfe1",
        "ENSMUST00000144177",
        "Adhfe1-203",
        "chr1",
        Strand.PLUS,
        9_549_907,
        9_549_944,
        "ENSMUSE00001146417",
        2,
    )
    intron_1 = Intron(
        "ENSMUSG00000025911",
        "Adhfe1",
        "ENSMUST00000144177",
        "Adhfe1-203",
        "chr1",
        Strand.PLUS,
        9_548_152,
        9_549_906,
    )

    transcript_1 = Transcript(
        "ENSMUSG00000025911",
        "Adhfe1",
        "ENSMUST00000144177",
        "Adhfe1-203",
        "chr1",
        Strand.PLUS,
        9_547_948,
        9_577_970,
        [],
    )
    assert str(transcript_1).startswith(
        "<Transcript ENSMUST00000144177 chr1:9,547,948-9,577,970 on '0' strand at 0x"
    )
    assert transcript_1.regions == [], "Exons should be empty list."
    # Test adding exons
    transcript_1.add_exon(exon_1)
    transcript_1.add_exon(exon_2)
    assert transcript_1.regions == [exon_1, exon_2], "Exons were not added correctly."
    # Exon that is already present should not be added again
    transcript_1.add_exon(exon_1)
    assert transcript_1.regions == [
        exon_1,
        exon_2,
    ], "Exon that was already present was not handled correctly."
    transcript_1.sort_regions()
    assert transcript_1.regions == [
        exon_1,
        intron_1,
        exon_2,
    ], "Exon sorting and intron insertion failed."

    # Check that introns are not inserted again
    transcript_1.sort_regions()
    assert transcript_1.regions == [
        exon_1,
        intron_1,
        exon_2,
    ], "Second sorting failed and/or inserted further introns."

    transcript_2 = Transcript(
        "ENSMUSG00000025911",
        "Adhfe1",
        "ENSMUST00000144177",
        "Adhfe1-203",
        "chr1",
        Strand.PLUS,
        9_547_948,
        9_577_970,
        [exon_1, exon_2],
    )
    assert (
        transcript_1 != transcript_2
    ), "Transcripts should differ in their exons/introns."

    transcript_2.sort_regions()
    assert transcript_1 == transcript_2, "Transcripts should be equal after sorting."


def test_transcript_index():
    transcript_idx = TranscriptIndex(gtf_file=TEST_GTF_FILE)

    # Check that the index is constructed as expected
    for (
        chromosome,
        strand,
    ), regions in transcript_idx.transcripts_by_region.items():
        assert regions == sorted(
            regions, key=lambda x: x[0]
        ), f"Regions for chromosome {chromosome} {strand} are not sorted"
        for pos, transcripts in regions:
            for transcript in transcripts:
                assert (
                    transcript.chromosome == chromosome
                ), f"Chromosome mismatch: {transcript.chromosome} != {chromosome}"
                assert (
                    transcript.strand == strand
                ), f"Strand mismatch: {transcript.strand} != {strand}"
                assert (
                    transcript.start <= pos <= transcript.end
                ), f"Transcript {transcript.transcript_id} does not overlap region starting at {pos}"

    # Test that no transcripts are returned when the chromosome/strand
    # is not in the index
    chromosome, strand, start, end = "chr22", Strand.PLUS, 1_234_567, 1_234_667
    assert not transcript_idx.get_overlapping_transcripts(
        chromosome=chromosome, strand=strand, start=start, end=end
    )

    # Test that no transcripts are returned when there is no overlap
    chromosome, strand, start, end = "chr1", Strand.PLUS, 8_000_000, 8_000_150
    assert not transcript_idx.get_overlapping_transcripts(
        chromosome=chromosome, strand=strand, start=start, end=end
    )
    start, end = 9_908_450, 9_908_600
    assert not transcript_idx.get_overlapping_transcripts(
        chromosome=chromosome, strand=strand, start=start, end=end
    )

    # Test that the correct overlapping transcripts are returned
    start, end = 9_886_058, 9_886_361
    overlapping_transcripts = transcript_idx.get_overlapping_transcripts(
        chromosome=chromosome, strand=strand, start=start, end=end
    )
    for transcript in overlapping_transcripts:
        assert transcript.overlaps(chromosome, strand, start, end), (
            f"Transcript {transcript.transcript_id} does not overlap "
            f"region {chromosome} {strand} {start} {end}"
        )
        assert transcript.gene_name == "Sgk3"
    assert (
        len(overlapping_transcripts) == 4
    ), f"Expected 4 transcripts, found {len(overlapping_transcripts)}"
    transcript_ids = [x.transcript_id for x in overlapping_transcripts]
    assert transcript_ids == [
        "ENSMUST00000166384",
        "ENSMUST00000168907",
        "ENSMUST00000171265",
        "ENSMUST00000097826",
    ]
    for i, trans_id in enumerate(transcript_ids):
        assert transcript_idx.get_transcript(trans_id) == overlapping_transcripts[i]

    # Test that the correct overlapping transcripts are returned
    start, end = 9_848_250, 9_848_450
    overlapping_transcripts = transcript_idx.get_overlapping_transcripts(
        chromosome=chromosome, strand=strand, start=start, end=end
    )
    for transcript in overlapping_transcripts:
        assert transcript.overlaps(chromosome, strand, start, end), (
            f"Transcript {transcript.transcript_id} does not overlap "
            f"region {chromosome} {strand} {start} {end}"
        )
    assert (
        len(overlapping_transcripts) == 5
    ), f"Expected 5 transcripts, found {len(overlapping_transcripts)}"
    transcript_ids = [x.transcript_id for x in overlapping_transcripts]
    assert transcript_ids == [
        "ENSMUST00000166384",
        "ENSMUST00000168907",
        "ENSMUST00000171265",
        "ENSMUST00000188830",
        "ENSMUST00000097826",
    ]

    # Test case where one transcript overlaps by only its last base
    start, end = 9_899_852, 9_899_951
    overlapping_transcripts = transcript_idx.get_overlapping_transcripts(
        chromosome=chromosome, strand=strand, start=start, end=end
    )
    transcript_ids = [x.transcript_id for x in overlapping_transcripts]
    assert transcript_ids == [
        "ENSMUST00000166384",
        "ENSMUST00000168907",
        "ENSMUST00000171265",
        "ENSMUST00000097826",  # This transcript overlaps by only its last base
    ]
    # Shift by one base to the right: the last transcript should not be included
    start, end = 9_899_853, 9_899_952
    overlapping_transcripts = transcript_idx.get_overlapping_transcripts(
        chromosome=chromosome, strand=strand, start=start, end=end
    )
    transcript_ids = [x.transcript_id for x in overlapping_transcripts]
    assert transcript_ids == [
        "ENSMUST00000166384",
        "ENSMUST00000168907",
        "ENSMUST00000171265",
    ]

    # Try to get a transcript that does not exist
    assert transcript_idx.get_transcript("ENSMUST00000000000") is None
