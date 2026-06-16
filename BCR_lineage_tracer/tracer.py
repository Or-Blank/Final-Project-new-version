"""
tracer.py
=========
LineageTracer  —  builds a germline-rooted B cell clonal lineage tree for a
single clone and reports per-node mutation events.

Algorithm overview
------------------
1.  Pairwise distance matrix
    Alignment-based Hamming distance (BioPython PairwiseAligner, global
    mode) normalised by alignment length.  Handles unequal sequence lengths
    that arise when the H+L concatenation differs by a few bases between
    cells.

2.  Neighbor-Joining tree
    BioPython DistanceTreeConstructor.nj().  The NJ algorithm is O(n³) but
    runs in < 1 s for the clone sizes present in these datasets (≤ 181 cells).

3.  Germline rooting
    The synthetic germline node is passed as the outgroup.  The redundant
    bifurcating root introduced by root_with_outgroup() is collapsed so the
    germline itself becomes the tree root (one node, not two).

4.  Polytomy collapsing
    Internal branches below `collapse_threshold` (default 1e-6) are
    collapsed into polytomies.  BCR sequences are short and low-divergence;
    forcing bifurcation would imply resolution that the data do not support.

5.  Isotype-aware NNI refinement  (optional, default on)
    A lightweight local-search pass performs Nearest-Neighbour Interchange
    moves that reduce the count of parent → child isotype transitions
    violating the irreversibility of Class Switch Recombination (CSR).
    This mirrors the conceptual goal of the TRIBAL algorithm without
    requiring a full probabilistic CSR model.
    Only moves that reduce (or tie with) total branch length are accepted,
    so the sequence-based topology is perturbed as little as possible.

6.  Ancestral sequence reconstruction  (Fitch parsimony)
    Each internal NJ node receives a reconstructed sequence via Fitch's
    maximum-parsimony algorithm.  This makes the mutation table meaningful
    for every edge in the tree, not just terminal branches.
    When sequences within a clone differ in length (uncommon — usually a
    sign of raw VDJ sequences without IMGT coordinate anchoring), the
    nearest-leaf sequence is used as a fallback.

7.  Mutation table
    For every non-germline node, nucleotide substitutions and the
    resulting amino-acid changes (relative to the parent node's sequence)
    are enumerated and returned as a pandas DataFrame.
"""

from __future__ import annotations

import warnings
from typing import Dict, List, Optional, Tuple

import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.Phylo.BaseTree import Clade, Tree
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio.Seq import Seq

from .constants import GAP_CHARS, ISOTYPE_ORDER
from .loader import CellRecord


# ── helpers ──────────────────────────────────────────────────────────────────

def _isotype_rank(iso) -> int:
    """Map an isotype string to its CSR position (lower = more naïve)."""
    if iso is None or iso in ("Germline", "Ancestral"):
        return -1   # always valid as a parent
    return ISOTYPE_ORDER.get(str(iso).strip().upper(), 3)


def _ungap(seq: str) -> str:
    for ch in GAP_CHARS:
        seq = seq.replace(ch, "")
    return seq


# ── LineageTracer ──────────────────────────────────────────────────────────────

class LineageTracer:
    """Build and annotate a lineage tree for one BCR clone.

    Parameters
    ----------
    records             : list of CellRecord (must include exactly one germline)
    collapse_threshold  : branch-length cutoff for polytomy collapsing
    refine_isotypes     : whether to run the CSR-aware NNI pass
    max_refine_iter     : safety cap on NNI iterations

    Example
    -------
    >>> tracer = LineageTracer(clone_records)
    >>> tree   = tracer.build()
    >>> df     = tracer.mutation_table()
    """

    def __init__(
        self,
        records: List[CellRecord],
        collapse_threshold: float = 1e-6,
        refine_isotypes: bool = True,
        max_refine_iter: int = 100,
    ):
        if not records:
            raise ValueError("records list is empty.")
        self.records = records
        self.by_id: Dict[str, CellRecord] = {r.cell_id: r for r in records}
        self.collapse_threshold = collapse_threshold
        self.refine_isotypes = refine_isotypes
        self.max_refine_iter = max_refine_iter

        germs = [r for r in records if r.is_germline]
        if not germs:
            raise ValueError("No germline record found in this clone.")
        self.germline_id: str = germs[0].cell_id

        self.tree: Optional[Tree] = None

        # BioPython global aligner — scores measure dissimilarity:
        # match=0, mismatch/gap=-1; used only when sequences differ in length
        self._aligner = PairwiseAligner()
        self._aligner.mode = "global"
        self._aligner.match_score = 0
        self._aligner.mismatch_score = -1
        self._aligner.open_gap_score = -1
        self._aligner.extend_gap_score = -1

    # ── public API ────────────────────────────────────────────────────────────

    def build(self) -> Tree:
        """Run the full pipeline and return the annotated Bio.Phylo Tree."""
        dm   = self._distance_matrix()
        tree = DistanceTreeConstructor().nj(dm)

        self._annotate(tree)
        self._root_on_germline(tree)
        self._collapse_polytomies(tree)
        if self.refine_isotypes:
            self._isotype_refine(tree)
        self._reconstruct_ancestors(tree)

        self.tree = tree
        return tree

    def mutation_table(self) -> pd.DataFrame:
        """Return a DataFrame with one row per non-germline tree edge.

        Columns
        -------
        clone_id, node, parent, is_observed, isotype, sample_id,
        cluster_annotated, branch_length,
        n_nt, nt_changes,   # nucleotide-level changes vs parent
        n_aa, aa_changes    # amino-acid-level changes vs parent
        """
        if self.tree is None:
            raise RuntimeError("Call build() before mutation_table()")

        clone_id = self.records[0].clone_id
        rows = []
        for cl in self.tree.find_clades(order="preorder"):
            if getattr(cl, "is_germline", False):
                continue
            par = self._parent(cl)
            if par is None:
                continue
            cs = getattr(cl,  "sequence", "")
            ps = getattr(par, "sequence", "")
            if not cs or not ps:
                continue
            nt, aa = self._diff(ps, cs)
            rows.append({
                "clone_id":          clone_id,
                "node":              cl.name,
                "parent":            par.name,
                "is_observed":       cl.name in self.by_id and not self.by_id[cl.name].is_germline,
                "isotype":           getattr(cl, "isotype", "?"),
                "sample_id":         getattr(cl, "sample_id", "?"),
                "cluster_annotated": getattr(cl, "cluster_annotated", "?"),
                "branch_length":     cl.branch_length,
                "number_nucleotides_changes":              len(nt),
                "nucleotides_changes":        ";".join(nt) if nt else "",
                "number_amino_acid_changes":              len(aa),
                "amino_acid_changes":        ";".join(aa) if aa else "",
            })
        return pd.DataFrame(rows)

    # ── distance matrix ────────────────────────────────────────────────────

    def _seq_dist(self, a: str, b: str) -> float:
        if not a or not b:
            return 1.0
        if len(a) == len(b):
            n = len(a)
            d = sum(1 for x, y in zip(a, b)
                    if x != y and x != "N" and y != "N")
            return d / n if n else 0.0
        # Unequal lengths: align first
        aln = self._aligner.align(a, b)[0]
        ra, rb = str(aln[0]), str(aln[1])
        n = len(ra)
        d = sum(1 for x, y in zip(ra, rb)
                if x != y
                and x not in GAP_CHARS and y not in GAP_CHARS
                and x != "N" and y != "N")
        return d / n if n else 0.0

    def _distance_matrix(self) -> DistanceMatrix:
        names = [r.cell_id for r in self.records]
        n = len(names)
        mat = []
        for i in range(n):
            row = []
            for j in range(i + 1):
                row.append(
                    0.0 if i == j
                    else round(self._seq_dist(
                        self.records[i].sequence,
                        self.records[j].sequence), 6)
                )
            mat.append(row)
        return DistanceMatrix(names=names, matrix=mat)

    # ── tree construction & manipulation ──────────────────────────────────

    def _annotate(self, tree: Tree):
        """Attach metadata from CellRecord objects to each clade."""
        for cl in tree.find_clades():
            rec = self.by_id.get(cl.name) if cl.name else None
            if rec:
                cl.isotype           = rec.isotype
                cl.sample_id         = rec.sample_id
                cl.cluster_annotated = rec.cluster_annotated
                cl.is_germline       = rec.is_germline
            else:
                # Internal NJ nodes have no observed cell → mark as ancestral
                cl.isotype           = "Ancestral"
                cl.sample_id         = "Ancestral"
                cl.cluster_annotated = "Ancestral"
                cl.is_germline       = False

    def _root_on_germline(self, tree: Tree):
        """Re-root tree so the germline is the single root node."""
        gc = next(
            (c for c in tree.find_clades() if c.name == self.germline_id), None
        )
        if gc is None:
            warnings.warn(
                "Germline clade not found; tree left on NJ midpoint root.",
                UserWarning, stacklevel=2,
            )
            return

        tree.root_with_outgroup(gc)

        # root_with_outgroup() creates a new bifurcating root with two
        # children: gc (branch_length ≈ 0) and the rest.  Collapse that
        # extra node so the germline itself is the root.
        nr = tree.root
        if gc in nr.clades and len(nr.clades) == 2:
            other = next(c for c in nr.clades if c is not gc)
            extra = gc.branch_length or 0.0
            gc.branch_length = 0.0
            other.branch_length = (other.branch_length or 0.0) + extra
            gc.clades = [other]
            tree.root = gc
        else:
            tree.root = gc

        # Ensure root carries germline metadata
        tree.root.isotype           = "Germline"
        tree.root.sample_id         = "Germline"
        tree.root.cluster_annotated = "Germline"
        tree.root.is_germline       = True

    def _clamp_branch_lengths(self, tree: Tree):
        """Clamp negative branch lengths to zero.
        BioPython's NJ implementation can produce small negative branch
        lengths due to floating-point rounding.  These are biologically
        meaningless and must be zeroed before collapse or display."""
        for cl in tree.find_clades():
            if cl.branch_length is not None and cl.branch_length < 0:
                cl.branch_length = 0.0

    def _collapse_polytomies(self, tree: Tree):
        """Collapse very-short internal branches into polytomies."""
        self._clamp_branch_lengths(tree)
        t = self.collapse_threshold
        tree.collapse_all(
            lambda c: (
                c is not tree.root
                and not c.is_terminal()
                and c.branch_length is not None
                and c.branch_length <= t
            )
        )

    # ── isotype-aware NNI ─────────────────────────────────────────────────

    def _count_violations(self, tree: Tree) -> int:
        """Count parent→child edges that violate CSR directionality."""
        v = 0
        for cl in tree.find_clades():
            for ch in cl.clades:
                if _isotype_rank(getattr(ch, "isotype", None)) < \
                   _isotype_rank(getattr(cl, "isotype", None)):
                    v += 1
        return v

    def _isotype_refine(self, tree: Tree):
        """Greedy NNI pass: swap grandchild ↔ sibling to reduce CSR violations."""
        for _ in range(self.max_refine_iter):
            cur_v = self._count_violations(tree)
            cur_l = tree.total_branch_length()
            best_gain, best_move = 0, None

            for node in [c for c in tree.find_clades() if c.clades]:
                for child in list(node.clades):
                    if not child.clades:
                        continue
                    for gc in list(child.clades):
                        for sib in list(node.clades):
                            if sib is child:
                                continue
                            self._swap(node, sib, child, gc)
                            nv = self._count_violations(tree)
                            nl = tree.total_branch_length()
                            gain = cur_v - nv
                            if gain > best_gain or (
                                gain == best_gain and gain > 0 and nl < cur_l
                            ):
                                best_gain = gain
                                best_move = (node, sib, child, gc)
                            self._swap(node, gc, child, sib)   # revert

            if best_move is None or best_gain <= 0:
                break
            self._swap(*best_move)

    @staticmethod
    def _swap(node: Clade, sib: Clade, child: Clade, gc: Clade):
        """Swap `sib` (child of `node`) with `gc` (grandchild via `child`)."""
        node.clades.remove(sib)
        child.clades.remove(gc)
        node.clades.append(gc)
        child.clades.append(sib)

    # ── Fitch ancestral sequence reconstruction ───────────────────────────

    def _reconstruct_ancestors(self, tree: Tree):
        """Assign a .sequence attribute to every clade.

        Observed leaves keep their actual sequence.
        Internal nodes get reconstructed ancestral sequences via Fitch
        maximum-parsimony.  The germline root retains its reference sequence.
        Falls back to nearest-leaf sequence when lengths differ within a clone.
        """
        seqs: Dict[str, str] = {r.cell_id: r.sequence for r in self.records}
        lengths = {len(s) for s in seqs.values()}

        if len(lengths) == 1:
            self._fitch(tree, seqs, lengths.pop())
        else:
            warnings.warn(
                f"Clone {self.records[0].clone_id}: sequences have "
                f"differing lengths {sorted(lengths)}. "
                "Falling back to nearest-leaf ancestral sequence assignment.",
                UserWarning, stacklevel=2,
            )
            for cl in tree.find_clades():
                if cl.is_terminal() and cl.name in seqs:
                    cl.sequence = seqs[cl.name]
                else:
                    leaves = cl.get_terminals()
                    cl.sequence = next(
                        (seqs[l.name] for l in leaves if l.name in seqs), ""
                    )

    def _fitch(self, tree: Tree, seqs: Dict[str, str], L: int):
        """Classic two-pass Fitch parsimony for ancestral sequence inference."""
        state: Dict[Clade, List[set]] = {}

        # Bottom-up: compute Fitch state sets
        for cl in tree.find_clades(order="postorder"):
            if cl.is_terminal():
                s = seqs.get(cl.name, "N" * L)
                state[cl] = [{ch} for ch in s]
            else:
                ch_sets = [state[c] for c in cl.clades]
                pos = []
                for i in range(L):
                    sets_i = [cs[i] for cs in ch_sets]
                    inter  = set.intersection(*sets_i)
                    pos.append(inter if inter else set.union(*sets_i))
                state[cl] = pos

        # Top-down: resolve each position preferring the parent's chosen base
        root = tree.root
        if getattr(root, "is_germline", False) and root.name in seqs:
            root.sequence = seqs[root.name]
        else:
            root.sequence = "".join(sorted(s)[0] for s in state[root])

        def assign(cl: Clade, par_seq: str):
            if not hasattr(cl, "sequence"):
                if cl.is_terminal() and cl.name in seqs:
                    cl.sequence = seqs[cl.name]
                else:
                    chars = []
                    for i in range(L):
                        pset = state[cl][i]
                        pc   = par_seq[i]
                        chars.append(pc if pc in pset else sorted(pset)[0])
                    cl.sequence = "".join(chars)
            for c in cl.clades:
                assign(c, cl.sequence)

        for c in root.clades:
            assign(c, root.sequence)

    # ── mutation diff ─────────────────────────────────────────────────────

    def _diff(self, parent_seq: str, child_seq: str
              ) -> Tuple[List[str], List[str]]:
        """Return (nt_changes, aa_changes) between parent and child sequences."""
        if len(parent_seq) == len(child_seq):
            ra, rb = parent_seq, child_seq
        else:
            aln = self._aligner.align(parent_seq, child_seq)[0]
            ra, rb = str(aln[0]), str(aln[1])

        nt_changes = [
            f"{pa}{i+1}{pb}"
            for i, (pa, pb) in enumerate(zip(ra, rb))
            if pa != pb
            and pa not in GAP_CHARS and pb not in GAP_CHARS
            and pa != "N" and pb != "N"
        ]

        # Amino-acid level: translate ungapped sequences in frame 1
        au, bu = _ungap(ra), _ungap(rb)
        at = au[: len(au) - len(au) % 3]
        bt = bu[: len(bu) - len(bu) % 3]
        try:
            aa_a = str(Seq(at).translate()) if at else ""
            aa_b = str(Seq(bt).translate()) if bt else ""
        except Exception:
            aa_a = aa_b = ""

        aa_changes = [
            f"{pa}{i+1}{pb}"
            for i, (pa, pb) in enumerate(zip(aa_a, aa_b))
            if pa != pb
        ]
        if len(aa_a) != len(aa_b):
            aa_changes.append(f"indel(len {len(aa_a)}->{len(aa_b)})")

        return nt_changes, aa_changes

    # ── tree helpers ───────────────────────────────────────────────────────

    def _parent(self, target: Clade) -> Optional[Clade]:
        if not self.tree:
            return None
        path = self.tree.get_path(target)
        if not path:
            return None
        return self.tree.root if len(path) == 1 else path[-2]
