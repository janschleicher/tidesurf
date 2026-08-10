import argparse
import glob
import logging
import os
import re
from datetime import datetime
from itertools import combinations
from pathlib import Path
from typing import List, Literal, Optional, Union

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
    parser.add_argument(
        "--bam_path",
        type=str,
        nargs="+",
        default=None,
        help="Explicit path to one or more BAM files. The sample "
        "directory will be ignored if this is given. If this argument "
        "is used, the positional arguments must be separated from it "
        "by another argument, by ' -- ', or they must precede it.",
    )
    arg_group = parser.add_mutually_exclusive_group()
    arg_group.add_argument(
        "--whitelist",
        type=str,
        nargs="+",
        help="Whitelist for cell filtering. Set to 'cellranger' to use "
        "barcodes in the sample directory. Alternatively, provide a "
        "path to a whitelist. If multiple BAM files are passed to "
        "'bam_path', one whitelist can be passed per BAM file. If this "
        "argument is used, the positional arguments must be separated "
        "from it by another argument, by ' -- ', or they must precede it.",
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
        "--export_umi_tables",
        action="store_true",
        help="Export tables with splice type for UMIs.",
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
    bam_path: Optional[Union[str, List[str]]] = None,
    whitelist: Optional[Union[str, List[str]]] = None,
    num_umis: Optional[int] = None,
    min_intron_overlap: int = 5,
    multi_mapped_reads: bool = False,
    export_umi_tables: bool = False,
) -> None:
    """
    Run tidesurf on a 10x sample directory.

    Parameters
    ----------
    sample_dir
        10x Cell Ranger count/multi output directory.
    gtf_file
        Path to GTF file with transcript annotations.
    output
        Path to output directory.
    orientation
        Orientation in which reads map to transcripts.
    filter_cells
        Whether to filter cells.
    bam_path
        Explicit path to BAM file. The sample directory will be ignored
        if an explicit path is given. Can also be a list with the paths
        to multiple BAM files. If this is used, 'sample_dir' still has
        to be passed, but can be any string, including ' '.
    whitelist
        If `filter_cells` is True: path to cell barcode whitelist file(s).
        Set to 'cellranger' to use barcodes in the sample directory.
        Mutually exclusive with `num_umis`.
    num_umis
        If `filter_cells` is True: set to an integer to only keep cells
        with at least that many UMIs. Mutually exclusive with `whitelist`.
    min_intron_overlap
        Minimum overlap of reads with introns required to consider them
        intronic.
    multi_mapped_reads
        Whether to count multi-mapped reads
    export_umi_tables
        Whether to export intermediate UMI tables. If True, for each
        cell, a table containing information about each UMI is saved
        to the output directory before (multi) and after (single)
        resolving UMIs mapping to multiple genes, as
        umi_table_[multi|single]_gene_<CELL_BARCODE>.parquet, respectively.
    """
    log.info("Building transcript index.")
    t_idx = TranscriptIndex(gtf_file)

    # Use explicitly provided BAM files if available
    if bam_path is not None:
        if whitelist == "cellranger":
            raise ValueError(
                "Cannot use 'cellranger' for whitelist with explicit BAM file paths."
            )
        log.info(
            f"Using explicitly specified BAM files and ignoring sample directory ({sample_dir})"
        )
        cr_pipeline = "unknown"
        if isinstance(bam_path, str):
            bam_files = [bam_path]
        else:
            bam_files = bam_path
        sample_ids = [f_path.split(".bam")[0].split("/")[-1] for f_path in bam_files]
    else:
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
            if bam_files:
                sample_ids = [
                    re.search(r"outs/per_sample_outs/(.*)/count", f) for f in bam_files
                ]
                sample_ids = [s_id.group(1) for s_id in sample_ids if s_id is not None]
            else:
                cr_pipeline = "unknown"
                log.info(
                    f"Sample directory {sample_dir} does not have Cell Ranger output "
                    f"structure; looking for any BAM files in the directory."
                )
                bam_files = glob.glob(
                    os.path.join(sample_dir, "**/*.bam"), recursive=True
                )
                if not bam_files:
                    raise FileNotFoundError(
                        f"No BAM files found in directory {sample_dir}."
                    )
                sample_ids = [
                    f_path.split(".bam")[0].split("/")[-1] for f_path in bam_files
                ]

    if len(sample_ids) > 1:
        # Append integer to sample IDs if not unique
        for sample_1, sample_2 in combinations(sample_ids, 2):
            if sample_1 == sample_2:
                sample_ids = [f"{s_id}_{i}" for i, s_id in enumerate(sample_ids)]
                break

    if isinstance(whitelist, list):
        if len(whitelist) == 1:
            whitelist = whitelist[0]
        elif len(whitelist) != len(bam_files):
            log.error(
                f"Mismatch between number of BAM files ({len(bam_files)}) and whitelists ({len(whitelist)})."
            )
            raise ValueError(
                f"Mismatch between number of BAM files ({len(bam_files)}) and whitelists ({len(whitelist)})."
            )

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

    for i, (bam_file, sample_id) in enumerate(zip(bam_files, sample_ids)):
        log.info(f"Processing {bam_file}.")
        if whitelist == "cellranger":
            if cr_pipeline == "count":
                wl = glob.glob(
                    f"{sample_dir}/outs/filtered_feature_bc_matrix/barcodes.*"
                )
            elif cr_pipeline == "multi":
                wl = glob.glob(
                    f"{sample_dir}/outs/per_sample_outs/{sample_id}/count/sample_filtered_feature_bc_matrix/barcodes.*"
                )
            else:
                wl = glob.glob(f"{sample_dir}/**/barcodes.*", recursive=True)
                log.info(wl)
            if not wl:
                log.error("No whitelist found in Cell Ranger output.")
                return
            else:
                wl = wl[0]
        elif isinstance(whitelist, list):
            # The case where it's just one whitelist is handled above
            wl = whitelist[i]
        else:
            wl = whitelist
        if num_umis is None:
            num_umis = -1
        cells, genes, counts = counter.count(
            bam_file=bam_file,
            filter_cells=filter_cells,
            whitelist=wl,
            num_umis=num_umis,
            umi_table_dir=output if export_umi_tables else None,
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
    if (isinstance(args.whitelist, list) and args.whitelist[0] == "cellranger") or (
        not args.no_filter_cells and not args.whitelist and not args.num_umis
    ):
        args.whitelist = "cellranger"

    log.info(f"Running tidesurf {tidesurf.__version__}.")
    log.info(f"Processing sample directory: {args.sample_dir}")
    run(
        sample_dir=args.sample_dir,
        gtf_file=args.gtf_file,
        output=args.output,
        orientation=args.orientation,
        filter_cells=not args.no_filter_cells,
        bam_path=args.bam_path,
        whitelist=args.whitelist,
        num_umis=args.num_umis,
        min_intron_overlap=args.min_intron_overlap,
        multi_mapped_reads=args.multi_mapped_reads,
        export_umi_tables=args.export_umi_tables,
    )
    end_time = datetime.now()
    log.info(f"Finished in {end_time - start_time}.")


if __name__ == "__main__":
    main()
