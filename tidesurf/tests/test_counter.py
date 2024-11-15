from tidesurf import UMICounter, TranscriptIndex
import numpy as np
import anndata as ad
import pytest


@pytest.mark.parametrize(
    "filter_cells, whitelist, num_umis",
    [
        (False, None, None),
        (
            True,
            "test_data/test_dir_count/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
            None,
        ),
        (True, "test_data/whitelist.tsv", None),
        (True, None, 10),
    ],
)
def test_counter(filter_cells: bool, whitelist: str, num_umis: int) -> None:
    t_idx = TranscriptIndex("test_data/genes.gtf")
    counter = UMICounter(
        transcript_index=t_idx, orientation="antisense", multi_mapped=False
    )
    cells, genes, counts = counter.count(
        bam_file="test_data/test_dir_count/outs/possorted_genome_bam.bam",
        filter_cells=filter_cells,
        whitelist=whitelist,
        num_umis=num_umis,
    )
    x_ts = (counts["spliced"] + counts["unspliced"] + counts["ambiguous"]).toarray()
    adata_cr = ad.read_h5ad("test_data/adata_cr_out.h5ad")
    x_cr = adata_cr[cells, genes].X.toarray()

    assert np.allclose(
        x_cr, x_ts, atol=5, rtol=0.05
    ), "Discrepancy between tidesurf and cellranger outputs is too big."

    for gene in adata_cr.var_names:
        assert (
            gene in genes or adata_cr[:, gene].X.sum() == 0
        ), f"Gene {gene} with nonzero count is missing in tidesurf output."


def test_counter_exceptions():
    t_idx = TranscriptIndex("test_data/genes.gtf")
    counter = UMICounter(
        transcript_index=t_idx, orientation="antisense", multi_mapped=False
    )
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
            num_umis=None,
        )
