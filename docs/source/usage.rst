Command line usage
==================

After installation, **tidesurf** can be run from the command line.
The following usage information is displayed when running the program with the ``-h`` or ``--help`` flag:

.. code-block:: console

    usage: tidesurf [-h] [-v] [--orientation {sense,antisense}] [-o OUTPUT]
                [--no_filter_cells] [--bam_path BAM_PATH [BAM_PATH ...]]
                [--whitelist WHITELIST [WHITELIST ...] | --num_umis NUM_UMIS]
                [--min_intron_overlap MIN_INTRON_OVERLAP]
                [--multi_mapped_reads] [--export_umi_tables]
                SAMPLE_DIR GTF_FILE

    Program: tidesurf (Tool for IDentification and Enumeration of Spliced and Unspliced Read Fragments)
    Version: 0.3.2.dev2

    positional arguments:
      SAMPLE_DIR            Sample directory containing Cell Ranger output.
      GTF_FILE              GTF file with transcript information.

    options:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit
      --orientation {sense,antisense}
                            Orientation of reads with respect to transcripts. For
                            10x Genomics, use 'sense' for three prime and
                            'antisense' for five prime.
      -o OUTPUT, --output OUTPUT
                            Output directory.
      --no_filter_cells     Do not filter cells.
      --bam_path BAM_PATH [BAM_PATH ...]
                            Explicit path to one or more BAM files. The sample
                            directory will be ignored if this is given. If this
                            argument is used, the positional arguments must be
                            separated from it by another argument, by ' -- ', or
                            they must precede it.
      --whitelist WHITELIST [WHITELIST ...]
                            Whitelist for cell filtering. Set to 'cellranger' to
                            use barcodes in the sample directory. Alternatively,
                            provide a path to a whitelist. If multiple BAM files
                            are passed to 'bam_path', one whitelist can be passed
                            per BAM file. If this argument is used, the positional
                            arguments must be separated from it by another
                            argument, by ' -- ', or they must precede it.
      --num_umis NUM_UMIS   Minimum number of UMIs for filtering a cell.
      --min_intron_overlap MIN_INTRON_OVERLAP
                            Minimum number of bases that a read must overlap with
                            an intron to be considered intronic.
      --multi_mapped_reads  Take reads mapping to multiple genes into account
                            (default: reads mapping to more than one gene are
                            discarded).
      --export_umi_tables   Export tables with splice type for UMIs.