"""
visualization.py
================
plot_tree() — renders a BCR clonal lineage tree as a rectangular cladogram.

Node labels
-----------
Observed leaf nodes are shown as  seq1, seq2, seq3 …  (assigned in the order
leaves appear top-to-bottom in the tree).  Internal inferred-ancestor nodes
are labelled  anc1, anc2 …  The germline root always shows "Germline".

The mapping  {original_cell_id → short_label}  is returned as a third element
of the tuple so the pipeline can add it to the Excel output.

X-axis
------
The cumulative mutation-distance axis now shows numeric tick marks so readers
can read off exact distances from the figure.

Visual design
-------------
- No axes frame (all four spines hidden).
- Legend anchored in figure coordinates so it never overlaps tree labels.
- Axes occupy the left 68 % of the figure; the remaining 32 % is label + legend space.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from Bio.Phylo.BaseTree import Clade, Tree

from .constants import MARKER_CYCLE


# ── layout ────────────────────────────────────────────────────────────────────

def _layout(tree: Tree) -> Tuple[Dict[Clade, float], Dict[Clade, float]]:
    """Return x (cumulative branch length) and y (leaf slot) position dicts."""
    x_pos: Dict[Clade, float] = {}

    def _assign_x(cl: Clade, x: float) -> None:
        x_pos[cl] = x
        for child in cl.clades:
            _assign_x(child, x + (child.branch_length or 0.0))

    _assign_x(tree.root, 0.0)

    y_pos: Dict[Clade, float] = {}
    for i, leaf in enumerate(tree.get_terminals()):
        y_pos[leaf] = float(i)

    def _assign_y(cl: Clade) -> float:
        if cl in y_pos:
            return y_pos[cl]
        child_ys = [_assign_y(c) for c in cl.clades]
        y = sum(child_ys) / len(child_ys)
        y_pos[cl] = y
        return y

    _assign_y(tree.root)
    return x_pos, y_pos


def _build_color_map(tree: Tree, color_by: str) -> Dict[str, object]:
    vals = sorted({str(getattr(c, color_by, "?")) for c in tree.find_clades()})
    cmap = plt.get_cmap("tab10")
    color_map: Dict[str, object] = {}
    for i, v in enumerate(vals):
        if v == "Germline":
            color_map[v] = "black"
        elif v == "Ancestral":
            color_map[v] = "lightgrey"
        else:
            color_map[v] = cmap(i % 10)
    return color_map


def _build_shape_map(tree: Tree, shape_by: str) -> Dict[str, str]:
    vals = sorted({str(getattr(c, shape_by, "?")) for c in tree.find_clades()})
    return {v: MARKER_CYCLE[i % len(MARKER_CYCLE)] for i, v in enumerate(vals)}


def _build_label_map(tree: Tree) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Assign short display labels to every clade.

    Rules
    -----
    - Germline root            → "Germline"
    - Observed leaf nodes      → "seq1", "seq2", … (top-to-bottom order)
    - Internal ancestral nodes → "anc1", "anc2", … (pre-order traversal)

    Returns
    -------
    display_map  : {clade.name → short_label}   used for tree drawing
    cell_id_map  : {original_cell_id → short_label}  returned to pipeline
                   (only contains observed cells, not Germline / ancestral)
    """
    display_map: Dict[str, str]  = {}
    cell_id_map: Dict[str, str]  = {}

    # Germline first
    if tree.root.name:
        display_map[tree.root.name] = "Germline"

    # Observed leaves  — numbered in top-to-bottom order
    seq_counter = 1
    for leaf in tree.get_terminals():
        if getattr(leaf, "is_germline", False):
            display_map[leaf.name] = "Germline"
        else:
            label = f"seq{seq_counter}"
            seq_counter += 1
            display_map[leaf.name] = label
            if leaf.name:
                cell_id_map[leaf.name] = label

    # Internal nodes (inferred ancestors) — numbered in pre-order
    anc_counter = 1
    for cl in tree.find_clades(order="preorder"):
        if cl.is_terminal():
            continue
        if getattr(cl, "is_germline", False):
            if cl.name:
                display_map[cl.name] = "Germline"
        else:
            if cl.name and cl.name not in display_map:
                display_map[cl.name] = f"anc{anc_counter}"
                anc_counter += 1

    return display_map, cell_id_map


# ── main public function ──────────────────────────────────────────────────────

def plot_tree(
    tree: Tree,
    color_by: str = "isotype",
    shape_by: Optional[str] = "cluster_annotated",
    title: str = "",
    output_path: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
    fig: Optional[plt.Figure] = None,
) -> Tuple[plt.Figure, plt.Axes, Dict[str, str]]:
    """Draw a rectangular cladogram for `tree`.

    Parameters
    ----------
    tree         : annotated Bio.Phylo Tree from LineageTracer.build()
    color_by     : clade attribute used for node fill colour
    shape_by     : clade attribute used for node marker shape (None = circles)
    title        : figure title
    output_path  : if given, save PNG here (150 dpi, tight bounding box)
    ax / fig     : optional pre-existing axes to draw into

    Returns
    -------
    (fig, ax, cell_id_map)
        cell_id_map : {original_cell_id → seq_label}  for Excel annotation
    """
    x_pos, y_pos   = _layout(tree)
    display_map, cell_id_map = _build_label_map(tree)
    n_leaves        = len(tree.get_terminals())

    # ── figure sizing ─────────────────────────────────────────────────────
    if ax is None:
        n_color = len({str(getattr(c, color_by, "?"))
                       for c in tree.find_clades()})
        n_shape = (len({str(getattr(c, shape_by, "?"))
                        for c in tree.find_clades()})
                   if shape_by else 0)
        n_legend_rows = n_color + n_shape

        legend_panel = max(3.5, n_legend_rows * 0.22)
        tree_panel   = 10.0
        fig_width    = tree_panel + legend_panel
        fig_height   = max(3.0, 0.35 * n_leaves)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Axes occupy the left 68 % of figure width.
    # The remaining 32 % is reserved for seq labels and the legend.
    fig.subplots_adjust(left=0.02, right=0.68, top=0.96, bottom=0.06)

    # ── draw branches ────────────────────────────────────────────────────
    for cl in tree.find_clades():
        if not cl.clades:
            continue
        px = x_pos[cl]
        child_ys = [y_pos[c] for c in cl.clades]
        ax.plot([px, px], [min(child_ys), max(child_ys)],
                color="grey", lw=1, zorder=1)
        for child in cl.clades:
            ax.plot([px, x_pos[child]], [y_pos[child], y_pos[child]],
                    color="grey", lw=1, zorder=1)

    # ── colour / shape maps ───────────────────────────────────────────────
    color_map = _build_color_map(tree, color_by)
    shape_map = _build_shape_map(tree, shape_by) if shape_by else {}

    # ── draw nodes and labels ─────────────────────────────────────────────
    for cl in tree.find_clades():
        x, y    = x_pos[cl], y_pos[cl]
        cv      = str(getattr(cl, color_by, "?"))
        is_germ = getattr(cl, "is_germline", False)

        marker = (
            shape_map.get(str(getattr(cl, shape_by, "?")), "o")
            if shape_by
            else ("s" if is_germ else "o")
        )

        ax.scatter(
            x, y,
            s=160 if is_germ else 70,
            marker=marker,
            facecolor=color_map.get(cv, "grey"),
            edgecolor="black",
            linewidth=0.8,
            zorder=3,
        )

        # Show short label for all leaf nodes and the germline root.
        # Internal ancestral nodes are NOT labelled to keep the tree clean.
        short = display_map.get(cl.name or "", "")
        if (cl.is_terminal() or is_germ) and short:
            ax.text(x, y, f"  {short}",
                    va="center", ha="left", fontsize=7, zorder=4)

    # ── legend ────────────────────────────────────────────────────────────
    handles: List[Line2D] = [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=color_map[v],
               markeredgecolor="black", markersize=8,
               label=f"{color_by}: {v}")
        for v in color_map
    ]
    if shape_by:
        for v, mk in shape_map.items():
            handles.append(
                Line2D([0], [0], marker=mk, color="w",
                       markerfacecolor="lightgrey",
                       markeredgecolor="black", markersize=8,
                       label=f"{shape_by}: {v}")
            )

    ax.legend(
        handles=handles,
        bbox_to_anchor=(0.70, 0.97),
        bbox_transform=fig.transFigure,
        loc="upper left",
        fontsize=6,
        frameon=True,
        framealpha=0.9,
        edgecolor="#cccccc",
        borderaxespad=0,
    )

    # ── axes cosmetics ────────────────────────────────────────────────────
    # Remove the frame (all four spines hidden).
    for spine in ax.spines.values():
        spine.set_visible(False)

    # X-axis: show numeric tick marks for mutation distance.
    # Use AutoLocator to pick sensible intervals, then format the numbers
    # with enough decimal places that they are readable at typical branch
    # length scales (1e-4 to 1e-1).
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=False))
    ax.ticklabel_format(style="plain", axis="x")
    ax.tick_params(axis="x", labelsize=7, length=4, color="#555555")

    # Keep only a faint horizontal grid at x-tick positions to help
    # readers trace values across the wide figure.
    ax.xaxis.grid(True, linestyle=":", linewidth=0.5,
                  color="#dddddd", zorder=0)
    ax.set_axisbelow(True)

    ax.set_xlabel("Cumulative mutation distance from germline",
                  labelpad=8, fontsize=8)
    ax.set_yticks([])
    ax.tick_params(left=False)
    ax.set_title(title, fontsize=9, pad=10)

    # ── save ──────────────────────────────────────────────────────────────
    if output_path and fig:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)

    return fig, ax, cell_id_map