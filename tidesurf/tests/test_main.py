import os
import shutil
from typing import Optional, Tuple

import anndata as ad
import numpy as np
import pytest

from tidesurf.main import main

OUT_DIR = "test_out"
TEST_GTF_FILE = "test_data/genes.gtf"
TEST_OUT_CR_5P = "test_data/adata_cr_out.h5ad"
TEST_OUT_CR_3P = "test_data/adata_cr_out_3p.h5ad"
TEST_OUT_TS_NO_FILTER = "test_data/tidesurf_out/tidesurf_no_filter.h5ad"
TEST_OUT_TS_FILTER_CR = "test_data/tidesurf_out/tidesurf_filter_cr.h5ad"
TEST_OUT_TS_FILTER_UMI = "test_data/tidesurf_out/tidesurf_filter_umi.h5ad"


def make_cmd(
    sample_dir: str,
    orientation: str,
    multi_mapped_reads: bool,
    no_filter_cells: bool,
    whitelist: Optional[str],
    num_umis: int,
) -> Tuple[str, str]:
    if orientation == "sense" and whitelist:
        whitelist = whitelist.replace("whitelist", "whitelist_3p")
    cmd = (
        f"tidesurf -o {OUT_DIR} --orientation {orientation} "
        f"{'--no_filter_cells ' if no_filter_cells else ''}"
        f"{f'--whitelist {whitelist} ' if whitelist else ''}"
        f"{f'--num_umis {num_umis} ' if num_umis != -1 else ''}"
        f"{'--multi_mapped_reads ' if multi_mapped_reads else ''}"
        f"{sample_dir} {TEST_GTF_FILE}"
    )

    return whitelist, cmd


def check_output(
    sample_dir: str,
    orientation: str,
    multi_mapped_reads: bool,
    no_filter_cells: bool,
    whitelist: Optional[str],
    num_umis: int,
    test_out_cr: str,
    test_out_ts: str,
):
    adata_cr = ad.read_h5ad(test_out_cr)
    if whitelist and num_umis != -1:
        assert not os.path.exists(OUT_DIR), (
            "No output should be generated with both whitelist and "
            "num_umis present (mutually exclusive arguments)."
        )
        return
    adata_ts = ad.read_h5ad(
        f"{OUT_DIR}/tidesurf.h5ad"
        if "count" in sample_dir
        else f"{OUT_DIR}/tidesurf_sample_1.h5ad"
    )

    # Compare with expected output
    if multi_mapped_reads:
        test_out_ts = test_out_ts.replace(".h5ad", "_mm.h5ad")
    if orientation == "sense":
        test_out_ts = test_out_ts.replace(".h5ad", "_3p.h5ad")
    adata_ts_true = ad.read_h5ad(test_out_ts)
    assert adata_ts_true.shape == adata_ts.shape, "Output shape mismatch."
    assert np.all(adata_ts_true.obs == adata_ts.obs), "Output obs mismatch."
    assert np.all(adata_ts_true.var == adata_ts.var), "Output var mismatch."
    assert np.all(
        adata_ts_true.X.toarray() == adata_ts.X.toarray()
    ), "Output X mismatch."
    for layer in adata_ts.layers.keys():
        assert np.all(
            adata_ts_true.layers[layer].toarray() == adata_ts.layers[layer].toarray()
        ), f"Output layer {layer} mismatch."

    # Check correct filtering
    if num_umis:
        assert np.all(adata_ts.X.sum(axis=1) >= num_umis), "Cells with too few UMIs."
    if not no_filter_cells:
        assert (
            set(adata_ts.obs_names) - set(adata_cr.obs_names) == set()
        ), "Cells found with tidesurf that are not in Cell Ranger output."

    # Compare with Cell Ranger output
    x_cr = adata_cr[adata_ts.obs_names, adata_ts.var_names].X.toarray()
    x_ts = adata_ts.X.toarray()

    assert np.allclose(
        x_cr, x_ts, atol=5, rtol=0.05
    ), "Discrepancy between tidesurf and cellranger outputs is too big."

    for gene in adata_cr.var_names:
        assert (
            gene in adata_ts.var_names or adata_cr[:, gene].X.sum() <= 1
        ), f"Gene {gene} with total count > 1 is missing in tidesurf output."

    # Make sure mitochondrial genes do not have unspliced or ambiguous counts
    assert (
        np.sum(
            adata_ts[:, adata_ts.var_names.str.contains("(?i)^MT-")].layers["unspliced"]
        )
        == 0
    ), "Mitochondrial genes do not have unspliced counts."

    assert (
        np.sum(
            adata_ts[:, adata_ts.var_names.str.contains("(?i)^MT-")].layers["ambiguous"]
        )
        == 0
    ), "Mitochondrial genes do not have ambiguous counts."

    shutil.rmtree(OUT_DIR)


@pytest.mark.parametrize(
    "sample_dir, orientation, test_out_cr",
    [
        (
            "test_data/test_dir_count",
            "antisense",
            TEST_OUT_CR_5P,
        ),
        (
            "test_data/test_dir_multi",
            "antisense",
            TEST_OUT_CR_5P,
        ),
        ("test_data/test_dir_count_3p", "sense", TEST_OUT_CR_3P),
    ],
)
@pytest.mark.parametrize("multi_mapped_reads", [False, True])
@pytest.mark.parametrize(
    "no_filter_cells, whitelist, num_umis, test_out_ts",
    [
        (True, None, -1, TEST_OUT_TS_NO_FILTER),
        (False, None, -1, TEST_OUT_TS_FILTER_CR),
        (False, "cellranger", -1, TEST_OUT_TS_FILTER_CR),
        (False, "test_data/whitelist.tsv", -1, TEST_OUT_TS_FILTER_CR),
        (False, None, 10, TEST_OUT_TS_FILTER_UMI),
        (False, "cellranger", 10, None),
    ],
)
def test_tidesurf(
    sample_dir: str,
    orientation: str,
    multi_mapped_reads: bool,
    no_filter_cells: bool,
    whitelist: Optional[str],
    num_umis: int,
    test_out_cr: str,
    test_out_ts: str,
):
    whitelist, cmd = make_cmd(
        sample_dir,
        orientation,
        multi_mapped_reads,
        no_filter_cells,
        whitelist,
        num_umis,
    )
    os.system(cmd)
    check_output(
        sample_dir,
        orientation,
        multi_mapped_reads,
        no_filter_cells,
        whitelist,
        num_umis,
        test_out_cr,
        test_out_ts,
    )


@pytest.mark.parametrize(
    "sample_dir, orientation, test_out_cr",
    [
        (
            "test_data/test_dir_count",
            "antisense",
            TEST_OUT_CR_5P,
        ),
        (
            "test_data/test_dir_multi",
            "antisense",
            TEST_OUT_CR_5P,
        ),
        ("test_data/test_dir_count_3p", "sense", TEST_OUT_CR_3P),
    ],
)
@pytest.mark.parametrize("multi_mapped_reads", [False, True])
@pytest.mark.parametrize(
    "no_filter_cells, whitelist, num_umis, test_out_ts",
    [
        (True, None, -1, TEST_OUT_TS_NO_FILTER),
        (False, None, -1, TEST_OUT_TS_FILTER_CR),
        (False, "cellranger", -1, TEST_OUT_TS_FILTER_CR),
        (False, "test_data/whitelist.tsv", -1, TEST_OUT_TS_FILTER_CR),
        (False, None, 10, TEST_OUT_TS_FILTER_UMI),
        (False, "cellranger", 10, None),
    ],
)
def test_main(
    sample_dir: str,
    orientation: str,
    multi_mapped_reads: bool,
    no_filter_cells: bool,
    whitelist: Optional[str],
    num_umis: int,
    test_out_cr: str,
    test_out_ts: str,
):
    whitelist, cmd = make_cmd(
        sample_dir,
        orientation,
        multi_mapped_reads,
        no_filter_cells,
        whitelist,
        num_umis,
    )

    arg_list = cmd.split(" ")[1:]

    if whitelist and num_umis != -1:
        with pytest.raises(SystemExit):
            main(arg_list)
        return
    main(arg_list)

    check_output(
        sample_dir,
        orientation,
        multi_mapped_reads,
        no_filter_cells,
        whitelist,
        num_umis,
        test_out_cr,
        test_out_ts,
    )
