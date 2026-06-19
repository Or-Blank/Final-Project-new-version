"""
main.py
===========
Command-line interface for BCR Lineage Tracer.

"""

import argparse
import sys

from .pipeline import run
from .gui import launch_gui


def main():
    parser = argparse.ArgumentParser(
        prog="bcr_lineage_tracer",
        description=(
            "BCR Lineage Tracer — reconstruct B cell receptor clonal lineage\n"
            "trees from single-cell VDJ sequencing data (.xlsx)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--input", "-i",
        help="Path to the input Excel workbook (.xlsx).",
    )
    parser.add_argument(
        "--output-dir", "-o",
        default="bcr_lineage_output",
        help="Directory to write tree images and mutation_table.xlsx  "
             "[default: bcr_lineage_output]",
    )
    parser.add_argument(
        "--clone-id",
        default=None,
        help="Process only this clone_id (useful for testing / debugging).",
    )
    parser.add_argument(
        "--collapse-threshold",
        type=float, default=1e-6,
        help="Branch-length threshold below which internal nodes are collapsed "
             "into polytomies  [default: 1e-6]",
    )
    parser.add_argument(
        "--no-isotype-refine",
        action="store_true",
        help="Disable the CSR-aware NNI refinement pass.",
    )
    parser.add_argument(
        "--max-clones",
        type=int, default=None,
        help="Process at most N clones (useful for quick checks).",
    )
    parser.add_argument(
        "--gui",
        action="store_true",
        help="Launch the graphical user interface.",
    )

    args = parser.parse_args()

    # No --input provided → fall back to GUI
    if args.gui or not args.input:
        launch_gui()
        return

    run(
        input_path=args.input,
        output_dir=args.output_dir,
        clone_id=args.clone_id,
        collapse_threshold=args.collapse_threshold,
        refine_isotypes=not args.no_isotype_refine,
        max_clones=args.max_clones,
    )


if __name__ == "__main__":
    main()
