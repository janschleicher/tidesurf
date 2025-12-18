from typing import List, Tuple

import pysam
import pytest

from tidesurf import TranscriptIndex, UMICounter
from tidesurf.enums import ReadType

TEST_GTF_FILE = "test_data/genes.gtf"
TRANSCRIPT_INDEX = TranscriptIndex(TEST_GTF_FILE)


def mock_read(
    barcode: str,
    identifier: str,
    reference_name: str,
    reference_start: int,
    cigar_tuples: List[Tuple[int, int]],
    mapping_quality: int = 255,
    is_unmapped: bool = False,
    is_reverse: bool = False,
    read_length: int = 100,
) -> pysam.AlignedSegment:
    """
    Create a mock read for testing purposes.
    :param barcode: Cell barcode.
    :param identifier: Unique molecular identifier.
    :param reference_name: Name of the reference sequence.
    :param reference_start: Position of the first aligned base in the
    reference sequence.
    :param reference_end: Position of the base after the last aligned
    base in the reference sequence.
    :param cigar_tuples: Cigar alignment tuples.
    :param mapping_quality: Mapping quality, integer between 0 and 255.
    :param is_unmapped: Whether the read is unmapped.
    :param is_reverse: Whether the read is mapped to the reverse strand.
    :param read_length: Length of the read.
    :return: Read object.
    """
    read = pysam.AlignedSegment(
        pysam.AlignmentHeader.from_dict(
            {"SQ": [{"SN": reference_name, "LN": 195_154_279}]}
        )
    )
    read.reference_name = reference_name
    read.reference_start = reference_start
    read.query_sequence = "A" * read_length
    read.cigar = cigar_tuples
    read.mapping_quality = mapping_quality
    read.is_unmapped = is_unmapped
    read.is_reverse = is_reverse
    read.set_tags([("CB", barcode), ("UB", identifier)])
    return read


READS = [
    (
        "ACTG",
        "TCGA",
        "chr1",
        9_865_240,
        [(pysam.CMATCH, 100)],
        False,
        [("Sgk3", ReadType.EXON)],
    ),
    (
        "ACTG",
        "TCGA",
        "chr1",
        9_865_340,
        [
            (pysam.CMATCH, 8),
            (pysam.CREF_SKIP, 3034),
            (pysam.CMATCH, 84),
            (pysam.CREF_SKIP, 3163),
            (pysam.CMATCH, 8),
        ],
        False,
        [("Sgk3", ReadType.EXON_EXON)],
    ),
    (
        "ACTG",
        "ATAT",
        "chr1",
        9_865_335,
        [(pysam.CMATCH, 13), (pysam.CREF_SKIP, 3034), (pysam.CMATCH, 87)],
        False,
        [("Sgk3", ReadType.EXON_EXON)],
    ),
    (
        "ACTG",
        "ATAT",
        "chr1",
        9_865_335,
        [(pysam.CMATCH, 100)],
        False,
        [("Sgk3", ReadType.INTRON)],
    ),
    (
        "ACTG",
        "ATAT",
        "chr1",
        9_866_123,
        [(pysam.CMATCH, 100)],
        False,
        [("Sgk3", ReadType.INTRON)],
    ),
    (
        "ACTG",
        "GTAT",
        "chr1",
        9_671_600,
        [(pysam.CMATCH, 100)],
        True,
        [("Mybl1", ReadType.AMBIGUOUS_READ)],
    ),
    (
        "ACTG",
        "GTAT",
        "chr1",
        9_672_550,
        [(pysam.CMATCH, 100)],
        True,
        [("Mybl1", ReadType.EXON)],
    ),
    (
        "ACTG",
        "GTAT",
        "chr1",
        9_672_600,
        [(pysam.CMATCH, 56), (pysam.CREF_SKIP, 444), (pysam.CMATCH, 44)],
        True,
        [("Mybl1", ReadType.EXON_EXON)],
    ),
    (
        "ACTG",
        "CACA",
        "chr1",
        9_672_550,
        [(pysam.CMATCH, 100)],
        True,
        [("Mybl1", ReadType.EXON)],
    ),
    (
        "ACTG",
        "CACA",
        "chr1",
        9_672_600,
        [(pysam.CMATCH, 100)],
        True,
        [("Mybl1", ReadType.INTRON)],
    ),
    (
        "ACTG",
        "TATA",
        "chr1",
        9_656_123,
        [(pysam.CMATCH, 100)],
        False,
        [("", None)],
    ),
    (
        "ACTG",
        "TATA",
        "chr1",
        9_672_600,
        [(pysam.CMATCH, 100)],
        False,
        [("", None)],
    ),
    (
        "ACTG",
        "TACT",
        "chrM",
        7_895,
        [(pysam.CMATCH, 100)],
        False,
        [("mt-Atp8", ReadType.EXON), ("mt-Atp6", ReadType.EXON)],
    ),
]


@pytest.mark.parametrize("multi_mapped_reads", [False, True])
def test_read_processing(multi_mapped_reads: bool) -> None:
    counter = UMICounter(
        TRANSCRIPT_INDEX, orientation="sense", multi_mapped_reads=multi_mapped_reads
    )
    for cbc, umi, ref_name, ref_start, cigar, is_reversed, expected_result in READS:
        read = mock_read(cbc, umi, ref_name, ref_start, cigar, is_reverse=is_reversed)
        res = counter._process_read(read)
        if not res:
            assert (expected_result == [("", None)]) or (
                not multi_mapped_reads and len(expected_result) > 1
            ), "Read was filtered out when it should not be."
        else:
            cbc_pred, res_list = res
            assert len(res_list) == len(
                expected_result
            ), f"Expected {len(expected_result)} gene(s), got {len(res_list)}."
            assert cbc_pred == cbc, "Cell barcode mismatch."
            for (
                umi_pred,
                gene_name_pred,
                read_type_pred,
                weight_pred,
            ) in res_list:
                assert umi_pred == umi, "Unique molecular identifier mismatch."
                assert (
                    gene_name_pred,
                    ReadType(read_type_pred),
                ) in expected_result, "Gene name or read type mismatch."
                assert weight_pred == 1.0 / len(expected_result), "Weight mismatch."
