"""
pipeline.py
===========
run() — top-level orchestration: loader → tracer → visualiser → Excel export.

Outputs per run
---------------
  tree_<clone_id>.png       Cladogram with seq1/seq2/… node labels and a
                            numeric x-axis showing mutation distance.

  node_labels.xlsx          One sheet per clone (or combined if many clones).
                            Columns:
                              clone_id | seq_label | cell_id | isotype
                              | sample_id | cluster_annotated
                            Maps every seqN label back to the original cell ID
                            so the tree and table can be cross-referenced.

  mutation_table.xlsx       Per-edge mutation events.  Now includes a
                            seq_label column so rows can be linked to the tree
                            by label name rather than raw cell ID.
"""

from __future__ import annotations

import itertools
import os
from typing import Callable, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from .loader import BCRTreeLoader, CellRecord
from .tracer import LineageTracer
from .visualization import plot_tree


def run(
    input_path: str,
    output_dir: str,
    clone_id: Optional[str] = None,
    collapse_threshold: float = 1e-6,
    refine_isotypes: bool = True,
    max_clones: Optional[int] = None,
    progress_cb: Optional[Callable[[str], None]] = None,
) -> Tuple[int, int, int, pd.DataFrame]:

    os.makedirs(output_dir, exist_ok=True)

    def log(msg: str) -> None:
        if progress_cb:
            progress_cb(msg)
        else:
            print(msg)

    # ── 1. Load ───────────────────────────────────────────────────────────
    loader = BCRTreeLoader(input_path).load()
    log(f"Format : {loader.format}")
    log(f"Rows   : {loader.df.shape[0]}")

    clones: Dict[str, List[CellRecord]] = loader.get_clones()
    log(f"Clones : {len(clones)} total")

    if clone_id:
        if clone_id not in clones:
            raise ValueError(
                f"clone_id '{clone_id}' not found in the input file.")
        clones = {clone_id: clones[clone_id]}

    if max_clones:
        clones = dict(itertools.islice(clones.items(), max_clones))

    # ── 2. Colour axis ────────────────────────────────────────────────────
    color_by = "sample_id" if loader.format == "heavy_only" else "isotype"

    # ── 3. Process each clone ─────────────────────────────────────────────
    all_mutation_tables: List[pd.DataFrame] = []
    all_label_tables:    List[pd.DataFrame] = []
    n_ok = n_skip = n_fail = 0

    # Build a lookup from cell_id → CellRecord for metadata enrichment
    record_lookup: Dict[str, CellRecord] = {
        r.cell_id: r
        for recs in clones.values()
        for r in recs
    }

    for cid, records in clones.items():
        n_obs = sum(1 for r in records if not r.is_germline)

        if n_obs < 2:
            log(f"  SKIP  clone {cid}  ({n_obs} observed cell — need ≥ 2)")
            n_skip += 1
            continue

        try:
            tracer = LineageTracer(
                records,
                collapse_threshold=collapse_threshold,
                refine_isotypes=refine_isotypes,
            )
            tree = tracer.build()
        except Exception as exc:
            log(f"  FAIL  clone {cid}  — {exc}")
            n_fail += 1
            continue

        # ── visualisation — returns cell_id → seq_label mapping
        out_png = os.path.join(output_dir, f"tree_{cid}.png")
        fig, _, cell_id_map = plot_tree(
            tree,
            color_by=color_by,
            shape_by="cluster_annotated",
            title=f"Clone {cid}  ({loader.format})",
            output_path=out_png,
        )
        plt.close("all")

        # ── node label table ──────────────────────────────────────────────
        # One row per observed cell: seq_label | cell_id | metadata
        label_rows = []
        for cell_id, seq_label in sorted(
            cell_id_map.items(), key=lambda kv: kv[1]   # sort by seq1,seq2,…
        ):
            rec = record_lookup.get(cell_id)
            label_rows.append({
                "clone_id":          cid,
                "seq_label":         seq_label,
                "cell_id":           cell_id,
                "isotype":           rec.isotype           if rec else "",
                "sample_id":         rec.sample_id         if rec else "",
                "cluster_annotated": rec.cluster_annotated if rec else "",
            })
        if label_rows:
            all_label_tables.append(pd.DataFrame(label_rows))

        # ── mutation table — add seq_label column ─────────────────────────
        mt = tracer.mutation_table()
        mt.insert(
            mt.columns.get_loc("node") + 1,   # place right after "node"
            "seq_label",
            mt["node"].map(cell_id_map).fillna(
                mt["node"].map(                # internal nodes get anc label
                    {c.name: f"anc"            # placeholder; viz assigned them
                     for c in tree.find_clades() if not c.is_terminal()}
                )
            ).fillna(""),
        )
        all_mutation_tables.append(mt)

        log(f"  ✓     clone {cid}  |  {n_obs} cells  "
            f"|  {len(cell_id_map)} seq labels  "
            f"|  {len(mt)} edges  →  {out_png}")
        n_ok += 1

    # ── 4. Export node label table ────────────────────────────────────────
    label_path = os.path.join(output_dir, "node_labels.xlsx")
    if all_label_tables:
        combined_labels = pd.concat(all_label_tables, ignore_index=True)
        combined_labels.to_excel(label_path, index=False)
        log(f"Node label table → {label_path}")
    else:
        combined_labels = pd.DataFrame()

    # ── 5. Export mutation table ──────────────────────────────────────────
    combined_mut = (
        pd.concat(all_mutation_tables, ignore_index=True)
        if all_mutation_tables
        else pd.DataFrame()
    )
    mut_path = os.path.join(output_dir, "mutation_table.xlsx")
    combined_mut.to_excel(mut_path, index=False)

    log("")
    log(f"Done  —  {n_ok} trees built,  "
        f"{n_skip} skipped (<2 cells),  {n_fail} failed.")
    log(f"Mutation table  → {mut_path}")

    return n_ok, n_skip, n_fail, combined_mut