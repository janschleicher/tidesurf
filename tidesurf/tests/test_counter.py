from tidesurf import UMICounter, TranscriptIndex
import numpy as np
import anndata as ad


def test_counter():
    t_idx = TranscriptIndex("test_data/genes.gtf")
    counter = UMICounter(
        transcript_index=t_idx, orientation="antisense", multi_mapped=False, threads=4
    )
    cells, genes, x_ts = counter.count(
        "test_data/test_dir_count/outs/possorted_genome_bam.bam"
    )
    adata_cr = ad.read_h5ad("test_data/adata_cr_out.h5ad")
    x_cr = adata_cr[cells, genes].X.toarray()

    assert np.allclose(
        x_cr, x_ts, atol=5, rtol=0.05
    ), "Discrepancy between tidesurf and cellranger outputs is too big."

    for gene in adata_cr.var_names:
        assert (
            gene in genes or adata_cr[:, gene].X.sum() == 0
        ), f"Gene {gene} with nonzero count is missing in tidesurf output."
