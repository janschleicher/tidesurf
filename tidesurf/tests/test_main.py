import os
import shutil
import numpy as np
import anndata as ad
import pytest
from typing import Optional

TEST_OUT_5P = "test_data/adata_cr_out.h5ad"


@pytest.mark.parametrize(
    "sample_dir, gtf_file, orientation, test_out",
    [
        ("test_data/test_dir_count", "test_data/genes.gtf", "antisense", TEST_OUT_5P),
        ("test_data/test_dir_multi", "test_data/genes.gtf", "antisense", TEST_OUT_5P),
    ],
)
@pytest.mark.parametrize("multi_mapped", [False])
@pytest.mark.parametrize(
    "filter_cells, whitelist, num_umis",
    [
        (False, None, None),
        (True, "cellranger", None),
        (True, "test_data/whitelist.tsv", None),
        (True, None, 10),
        (True, "cellranger", 10),
    ],
)
def test_main(
    sample_dir: str,
    gtf_file: str,
    orientation: str,
    multi_mapped: bool,
    filter_cells: bool,
    whitelist: Optional[str],
    num_umis: Optional[int],
    test_out: str,
):
    os.system(
        f"tidesurf -o test_out --orientation {orientation} "
        f"{'-m ' if multi_mapped else ''}"
        f"{'--filter_cells ' if filter_cells else ''}"
        f"{f'--whitelist {whitelist} ' if whitelist else ''}"
        f"{f'--num_umis {num_umis} ' if num_umis else ''}"
        f"{sample_dir} {gtf_file}"
    )
    adata_cr = ad.read_h5ad(TEST_OUT_5P)
    if whitelist and num_umis:
        assert not os.path.exists("test_out"), (
            "No output should be generated with both whitelist and "
            "num_umis present (mutually exclusive arguments)."
        )
        return
    adata_ts = ad.read_h5ad(
        "test_out/tidesurf.h5ad"
        if sample_dir.endswith("count")
        else "test_out/tidesurf_sample_1.h5ad"
    )
    if num_umis:
        assert np.all(adata_ts.X.sum(axis=1) >= num_umis), "Cells with too few UMIs."
    if filter_cells:
        assert (
            set(adata_ts.obs_names) - set(adata_cr.obs_names) == set()
        ), "Cells found with tidesurf that are not in Cell Ranger output."
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
