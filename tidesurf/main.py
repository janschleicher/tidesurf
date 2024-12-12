import argparse
import os
import glob
import re
from datetime import datetime
from pathlib import Path
import anndata as ad
import tidesurf
from tidesurf.transcript import TranscriptIndex
from tidesurf.counter import UMICounter
from typing import Literal, Optional
import logging

log = logging.getLogger(__name__)


def run(
    sample_dir: str,
    gtf_file: str,
    output: str,
    orientation: Literal["sense", "antisense"] = "sense",
    filter_cells: bool = False,
    whitelist: Optional[str] = None,
    num_umis: Optional[int] = None,
) -> None:
    log.info("Building transcript index.")
    t_idx = TranscriptIndex(gtf_file)
    cr_pipeline = "count"
    # Try cellranger count output
    bam_files = [f"{sample_dir}/outs/possorted_genome_bam.bam"]
    sample_ids = [""]
    if not os.path.isfile(bam_files[0]):
        cr_pipeline = "multi"
        # Try cellranger multi output
        bam_files = glob.glob(
            f"{sample_dir}/outs/per_sample_outs/*/count/sample_alignments.bam"
        )
        if not bam_files:
            log.error(f"No Cell Ranger BAM files found in directory {sample_dir}.")
            raise FileNotFoundError(
                f"No Cell Ranger BAM files found in directory {sample_dir}."
            )
        sample_ids = [
            re.search(r"outs/per_sample_outs/(.*)/count", f).group(1) for f in bam_files
        ]

    counter = UMICounter(transcript_index=t_idx, orientation=orientation)
    log.info(
        f"Counting reads mapped to transcripts in {counter.orientation} orientation."
    )

    os.makedirs(output, exist_ok=True)

    for bam_file, sample_id in zip(bam_files, sample_ids):
        log.info(f"Processing {bam_file}.")
        if whitelist == "cellranger":
            if cr_pipeline == "count":
                wl = glob.glob(
                    f"{sample_dir}/outs/filtered_feature_bc_matrix/barcodes.*"
                )
            else:
                wl = glob.glob(
                    f"{sample_dir}/outs/per_sample_outs/{sample_id}/count/sample_filtered_feature_bc_matrix/barcodes.*"
                )
            if not wl:
                log.error("No whitelist found in Cell Ranger output.")
                return
            else:
                wl = wl[0]
        else:
            wl = whitelist
        cells, genes, counts = counter.count(
            bam_file=bam_file,
            filter_cells=filter_cells,
            whitelist=wl,
            num_umis=num_umis,
        )
        log.info("Writing output.")
        counts_matrix = counts["spliced"] + counts["unspliced"] + counts["ambiguous"]
        adata = ad.AnnData(X=counts_matrix, layers=counts)
        adata.obs_names = cells
        adata.var_names = genes
        f_name = "tidesurf.h5ad" if not sample_id else f"tidesurf_{sample_id}.h5ad"
        adata.write_h5ad(Path(os.path.join(output, f_name)))


def main() -> None:
    start_time = datetime.now()
    parser = argparse.ArgumentParser(
        description="Program: tidesurf (Tool for IDentification and "
        "Enumeration of Spliced and Unspliced Read Fragments)\n"
        f"Version: {tidesurf.__version__}",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s {tidesurf.__version__}",
    )
    parser.add_argument(
        "--orientation",
        type=str,
        default="sense",
        choices=["sense", "antisense"],
        help="Orientation of reads with respect to transcripts. For 10x"
        " Genomics, use 'sense' for three prime and 'antisense' for five prime.",
    )
    parser.add_argument(
        "-o", "--output", type=str, default="tidesurf_out", help="Output directory."
    )
    parser.add_argument(
        "--filter_cells",
        action="store_true",
        help="Filter cells based on a whitelist.",
    )
    arg_group = parser.add_mutually_exclusive_group()
    arg_group.add_argument(
        "--whitelist",
        type=str,
        help="Whitelist for cell filtering. Set to 'cellranger' to use "
        "barcodes in the sample directory. Alternatively, provide a "
        "path to a whitelist.",
    )
    arg_group.add_argument(
        "--num_umis",
        type=int,
        help="Minimum number of UMIs for filtering a cell.",
    )
    parser.add_argument(
        "sample_dir",
        metavar="SAMPLE_DIR",
        help="Sample directory containing Cell Ranger output.",
    )
    parser.add_argument(
        "gtf_file", metavar="GTF_FILE", help="GTF file with transcript information."
    )
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s][%(module)s][%(levelname)s] - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.FileHandler(
                filename=os.path.join(args.output, "tidesurf.log"),
                mode="w",
            ),
            logging.StreamHandler(),
        ],
    )

    log.info(f"Running tidesurf {tidesurf.__version__}.")
    log.info(f"Processing sample directory: {args.sample_dir}")
    run(
        sample_dir=args.sample_dir,
        gtf_file=args.gtf_file,
        output=args.output,
        orientation=args.orientation,
        filter_cells=args.filter_cells,
        whitelist=args.whitelist,
        num_umis=args.num_umis,
    )
    end_time = datetime.now()
    log.info(f"Finished in {end_time - start_time}.")


if __name__ == "__main__":
    main()
