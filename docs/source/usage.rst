Command line usage
==================

After installation, **tidesurf** can be run from the command line.
The following usage information is displayed when running the program with the ``-h`` or ``--help`` flag:

.. code-block:: console

    usage: tidesurf [-h] [-v] [--orientation {sense,antisense}] [-o OUTPUT]
                [--no_filter_cells]
                [--whitelist WHITELIST | --num_umis NUM_UMIS]
                [--min_intron_overlap MIN_INTRON_OVERLAP]
                [--multi_mapped_reads]
                SAMPLE_DIR GTF_FILE

    Program: tidesurf (Tool for IDentification and Enumeration of Spliced and Unspliced Read Fragments)
    Version: 0.2.0

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
      --whitelist WHITELIST
                            Whitelist for cell filtering. Set to 'cellranger' to
                            use barcodes in the sample directory. Alternatively,
                            provide a path to a whitelist.
      --num_umis NUM_UMIS   Minimum number of UMIs for filtering a cell.
      --min_intron_overlap MIN_INTRON_OVERLAP
                            Minimum number of bases that a read must overlap with
                            an intron to be considered intronic.
      --multi_mapped_reads  Take reads mapping to multiple genes into account
                            (default: reads mapping to more than one gene are
                            discarded).