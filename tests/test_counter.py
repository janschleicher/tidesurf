import anndata as ad
import numpy as np
import pandas as pd
import pytest

from tidesurf import TranscriptIndex, UMICounter


@pytest.mark.parametrize(
    "filter_cells, whitelist, num_umis, five_prime",
    [
        (False, None, -1, True),
        (
            True,
            "test_data/test_dir_count/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
            -1,
            True,
        ),
        (True, "test_data/whitelist.tsv", -1, True),
        (True, None, 10, True),
        (False, None, -1, False),
        (
            True,
            "test_data/test_dir_count_3p/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
            -1,
            False,
        ),
        (True, "test_data/whitelist_3p.tsv", -1, False),
        (True, None, 10, False),
    ],
)
def test_counter(
    filter_cells: bool, whitelist: str, num_umis: int, five_prime: bool
) -> None:
    t_idx = TranscriptIndex("test_data/genes.gtf")
    if five_prime:
        orientation = "antisense"
    else:
        orientation = "sense"
    counter = UMICounter(transcript_index=t_idx, orientation=orientation)
    cells, genes, counts = counter.count(
        bam_file=f"test_data/{'test_dir_count' if five_prime else 'test_dir_count_3p'}/outs/possorted_genome_bam.bam",
        filter_cells=filter_cells,
        whitelist=whitelist,
        num_umis=num_umis,
    )
    x_ts = (counts["spliced"] + counts["unspliced"] + counts["ambiguous"]).toarray()
    if five_prime:
        adata_cr = ad.read_h5ad("test_data/adata_cr_out.h5ad")
    else:
        adata_cr = ad.read_h5ad("test_data/adata_cr_out_3p.h5ad")
    x_cr = adata_cr[cells, genes].X.toarray()

    assert np.allclose(
        x_cr, x_ts, atol=5, rtol=0.05
    ), "Discrepancy between tidesurf and cellranger outputs is too big."

    for gene in adata_cr.var_names:
        assert (
            gene in genes or adata_cr[:, gene].X.sum() <= 1
        ), f"Gene {gene} with total count > 1 is missing in tidesurf output."

    # Mitochondrial genes should not have any unspliced or ambiguous counts
    assert (
        np.sum(
            counts["unspliced"].toarray()[:, pd.Series(genes).str.contains("(?i)^MT-")]
        )
        == 0
    ), "Mitochondrial genes do not have unspliced counts."
    assert (
        np.sum(
            counts["ambiguous"].toarray()[:, pd.Series(genes).str.contains("(?i)^MT-")]
        )
        == 0
    ), "Mitochondrial genes do not have ambiguous counts."


def test_counter_exceptions():
    t_idx = TranscriptIndex("test_data/genes.gtf")
    counter = UMICounter(transcript_index=t_idx, orientation="antisense")
    with pytest.raises(
        ValueError, match="Whitelist and num_umis are mutually exclusive arguments."
    ):
        counter.count(
            bam_file="test_data/test_dir_count/outs/possorted_genome_bam.bam",
            filter_cells=True,
            whitelist="test_data/whitelist.tsv",
            num_umis=10,
        )

    with pytest.raises(
        ValueError,
        match="Either whitelist or num_umis must be provided when filter_cells==True.",
    ):
        counter.count(
            bam_file="test_data/test_dir_count/outs/possorted_genome_bam.bam",
            filter_cells=True,
            whitelist=None,
            num_umis=-1,
        )
