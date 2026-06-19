"""
constants.py
============
Biological and visual constants used across the BCR Lineage Tracer package.

"""

from typing import Dict

# ---------------------------------------------------------------------------
# Isotype class-switch recombination (CSR) ordering
# ---------------------------------------------------------------------------
# Rank along the IGH locus (5' → 3').  Lower rank = "earlier" in CSR.
# A valid parent → child edge requires:
#   child_rank >= parent_rank
# i.e. you can switch *forward* (IgM → IgG) but never backward (IgG → IgM).

# Aliases cover both IMGT-style names (IGHG1) and short names (IGG1) and
# the mouse-specific subclasses (IGHG2B, IGHG2C) present in the heavy-only
# dataset.  Unknown isotypes get rank 3 (neutral mid-point) so they never
# trigger a spurious CSR violation.

ISOTYPE_ORDER: Dict[str, int] = {
    # IgM / IgD  —  naïve isotypes, always upstream
    "IGHM": 0, "IGM": 0,
    "IGHD": 0, "IGD": 0,
    # IgG sub-classes (human order: G3 → G1 → G2 → G4)
    "IGHG3": 1, "IGG3": 1,
    "IGHG1": 2, "IGG1": 2,
    "IGHA1": 3, "IGA1": 3,   # IgA1 sits between G1 and G2 on the locus
    "IGHG2": 4, "IGG2": 4,
    # Mouse-specific IgG2 sub-classes (present in All_clones dataset)
    "IGHG2B": 4,
    "IGHG2C": 4,
    "IGHG4": 5, "IGG4": 5,
    "IGHE":  6, "IGE":  6,
    "IGHA2": 7, "IGA2": 7,
}

# ---------------------------------------------------------------------------
# Sequence gap / placeholder characters
# ---------------------------------------------------------------------------
# IMGT-coordinate alignments use '.' for framework insertions relative to
# the reference; '-' is used by BioPython pairwise alignments.  Both must
# be stripped before computing Hamming distances and before Fitch parsimony.

GAP_CHARS = (".", "-")

# ---------------------------------------------------------------------------
# Matplotlib marker cycle for cluster_annotated (cell-type) encoding
# ---------------------------------------------------------------------------
# Cycles through this list when assigning node shapes in the tree plot.
# Using distinct geometric shapes makes cell types distinguishable even on
# a black-and-white printout.

MARKER_CYCLE = ["o", "s", "^", "D", "v", "P", "X", "*", "h", "<", ">"]
