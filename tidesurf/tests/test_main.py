import os
import shutil
import numpy as np
import anndata as ad
import pytest

TEST_OUT_5P = "test_data/adata_cr_out.h5ad"


@pytest.mark.parametrize(
    "sample_dir, gtf_file, orientation, test_out",
    [
        ("test_data/test_dir_count", "test_data/genes.gtf", "antisense", TEST_OUT_5P),
        ("test_data/test_dir_multi", "test_data/genes.gtf", "antisense", TEST_OUT_5P),
    ],
)
@pytest.mark.parametrize("multi_mapped", [False])
def test_main(
    sample_dir: str, gtf_file: str, orientation: str, multi_mapped: bool, test_out: str
):
    os.system(
        f"tidesurf -o test_out --orientation {orientation} {'-m' if multi_mapped else ''} -t 4 {sample_dir} {gtf_file}"
    )
    adata_cr = ad.read_h5ad(TEST_OUT_5P)
    adata_ts = ad.read_h5ad(
        "test_out/tidesurf.h5ad"
        if sample_dir.endswith("count")
        else "test_out/tidesurf_sample_1.h5ad"
    )
    x_cr = adata_cr[adata_ts.obs_names, adata_ts.var_names].X.toarray()
    x_ts = adata_ts.X.toarray()

    assert np.allclose(
        x_cr, x_ts, atol=5, rtol=0.05
    ), "Discrepancy between tidesurf and cellranger outputs is too big."

    for gene in adata_cr.var_names:
        assert (
            gene in adata_ts.var_names or adata_cr[:, gene].X.sum() == 0
        ), f"Gene {gene} with nonzero count is missing in tidesurf output."

    shutil.rmtree("test_out")
