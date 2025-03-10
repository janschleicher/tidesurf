import argparse
import glob
import logging
import os
import re
from datetime import datetime
from pathlib import Path
from typing import List, Literal, Optional

import anndata as ad
from cython.cimports.tidesurf.counter import UMICounter
from cython.cimports.tidesurf.transcript import TranscriptIndex

import tidesurf

log = logging.getLogger(__name__)


def parse_args(arg_list: Optional[List[str]] = None) -> argparse.Namespace:
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
        " Genomics, use 'sense' for three prime and 'antisense' for "
        "five prime.",
    )
    parser.add_argument(
        "-o", "--output", type=str, default="tidesurf_out", help="Output directory."
    )
    parser.add_argument(
        "--no_filter_cells",
        action="store_true",
        help="Do not filter cells.",
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
        "--min_intron_overlap",
        type=int,
        default=5,
        help="Minimum number of bases that a read must overlap with an "
        "intron to be considered intronic.",
    )
    parser.add_argument(
        "--multi_mapped_reads",
        action="store_true",
        help="Take reads mapping to multiple genes into account "
        "(default: reads mapping to more than one gene are discarded).",
    )
    parser.add_argument(
        "sample_dir",
        metavar="SAMPLE_DIR",
        help="Sample directory containing Cell Ranger output.",
    )
    parser.add_argument(
        "gtf_file", metavar="GTF_FILE", help="GTF file with transcript information."
    )
    return parser.parse_args(arg_list)


def run(
    sample_dir: str,
    gtf_file: str,
    output: str,
    orientation: Literal["sense", "antisense"] = "sense",
    filter_cells: bool = False,
    whitelist: Optional[str] = None,
    num_umis: Optional[int] = None,
    min_intron_overlap: int = 5,
    multi_mapped_reads: bool = False,
) -> None:
    """
    Run tidesurf on a 10x sample directory.
    :param sample_dir: 10x Cell Ranger count/multi output directory.
    :param gtf_file: Path to GTF file with transcript annotations.
    :param output: Path to output directory.
    :param orientation: Orientation in which reads map to transcripts.
    :param filter_cells: Whether to filter cells.
    :param whitelist: If `filter_cells` is True: path to cell
        barcode whitelist file. Set to 'cellranger' to use barcodes in
        the sample directory. Mutually exclusive with `num_umis`.
    :param num_umis: If `filter_cells` is True: set to an integer to
        only keep cells with at least that many UMIs. Mutually exclusive
        with `whitelist`.
    :param min_intron_overlap: Minimum overlap of reads with introns
        required to consider them intronic.
    :param multi_mapped_reads: Whether to count multi-mapped reads
    :return:
    """
    log.info("Building transcript index.")
    t_idx = TranscriptIndex(gtf_file)
    cr_pipeline = "count"
    # Try cellranger count output
    bam_files = [os.path.join(sample_dir, "outs/possorted_genome_bam.bam")]
    sample_ids = [""]
    if not os.path.isfile(bam_files[0]):
        cr_pipeline = "multi"
        # Try cellranger multi output
        bam_files = glob.glob(
            os.path.join(
                sample_dir, "outs/per_sample_outs/*/count/sample_alignments.bam"
            )
        )
        if not bam_files:
            log.error(f"No Cell Ranger BAM files found in directory {sample_dir}.")
            raise FileNotFoundError(
                f"No Cell Ranger BAM files found in directory {sample_dir}."
            )
        sample_ids = [
            re.search(r"outs/per_sample_outs/(.*)/count", f).group(1) for f in bam_files
        ]

    counter = UMICounter(
        transcript_index=t_idx,
        orientation=orientation,
        min_intron_overlap=min_intron_overlap,
        multi_mapped_reads=multi_mapped_reads,
    )
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
        if num_umis is None:
            num_umis = -1
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


def main(arg_list: Optional[List[str]] = None) -> None:
    start_time = datetime.now()

    args = parse_args(arg_list)

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

    # Default behavior for filtering: use cellranger whitelist
    if not args.no_filter_cells and not args.whitelist and not args.num_umis:
        args.whitelist = "cellranger"

    log.info(f"Running tidesurf {tidesurf.__version__}.")
    log.info(f"Processing sample directory: {args.sample_dir}")
    run(
        sample_dir=args.sample_dir,
        gtf_file=args.gtf_file,
        output=args.output,
        orientation=args.orientation,
        filter_cells=not args.no_filter_cells,
        whitelist=args.whitelist,
        num_umis=args.num_umis,
        min_intron_overlap=args.min_intron_overlap,
        multi_mapped_reads=args.multi_mapped_reads,
    )
    end_time = datetime.now()
    log.info(f"Finished in {end_time - start_time}.")


if __name__ == "__main__":
    main()
