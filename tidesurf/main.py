import argparse
import os
import glob
import re
from pathlib import Path
from scipy.sparse import csr_matrix
import anndata as ad
import tidesurf
from typing import Literal
import logging

log = logging.getLogger(__name__)


def run(
    sample_dir: str,
    gtf_file: str,
    output: str,
    orientation: Literal["sense", "antisense"] = "sense",
    multi_mapped: bool = False,
    threads: int = 1,
) -> None:
    log.info("Building transcript index.")
    t_idx = tidesurf.TranscriptIndex(gtf_file)
    # Try cellranger count output
    bam_files = [f"{sample_dir}/outs/possorted_genome_bam.bam"]
    sample_ids = [""]
    if not os.path.isfile(bam_files[0]):
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

    counter = tidesurf.UMICounter(
        transcript_index=t_idx,
        orientation=orientation,
        multi_mapped=multi_mapped,
        threads=threads,
    )
    log.info(
        f"Counting reads mapped to transcripts in {counter.orientation} orientation."
    )

    os.makedirs(output, exist_ok=True)

    for bam_file, sample_id in zip(bam_files, sample_ids):
        log.info(f"Processing {bam_file}.")
        cells, genes, counts = counter.count(bam_file=bam_file)
        log.info("Writing output.")
        adata = ad.AnnData(X=csr_matrix(counts))
        adata.obs_names = cells
        adata.var_names = genes
        f_name = "tidesurf.h5ad" if not sample_id else f"tidesurf_{sample_id}.h5ad"
        adata.write_h5ad(Path(os.path.join(output, f_name)))


def main() -> None:
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
        "-t", "--threads", type=int, default=1, help="Number of threads"
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
        "-m",
        "--multi_mapped",
        action="store_true",
        help="Include multi-mapped reads (not recommended).",
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

    log.info(
        f"Running tidesurf {tidesurf.__version__} with {args.threads} thread{'s' if args.threads != 1 else ''}."
    )
    log.info(f"Processing sample directory: {args.sample_dir}")
    run(
        sample_dir=args.sample_dir,
        gtf_file=args.gtf_file,
        output=args.output,
        orientation=args.orientation,
        multi_mapped=args.multi_mapped,
        threads=args.threads,
    )
