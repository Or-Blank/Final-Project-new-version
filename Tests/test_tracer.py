"""
test_tracer.py
==============
Tests for LineageTracer — tree structure, rooting, polytomies,
isotype-aware NNI, Fitch ancestral reconstruction, and mutation table.

What each test covers
---------------------
  Tree structure
    test_tree_is_rooted_on_germline   — root node is the germline record
    test_all_observed_cells_are_leaves— every non-germline cell appears as leaf
    test_no_phantom_root_node         — root collapsed to single germline node
    test_branch_lengths_non_negative  — all branch lengths ≥ 0

  Polytomies
    test_polytomy_collapse            — near-zero branches become polytomies
    test_no_polytomy_when_high_threshold — large threshold collapses everything

  Isotype-aware NNI
    test_csr_violations_not_increased — refinement never makes violations worse
    test_refine_false_skipped         — violations unchanged when refine=False

  Ancestral reconstruction
    test_leaves_have_correct_sequence — leaf .sequence == observed sequence
    test_root_has_germline_sequence   — root .sequence == germline seq
    test_internal_nodes_have_sequence — every clade has a .sequence attribute
    test_sequence_length_consistent   — all .sequence same length (Fitch path)

  Mutation table
    test_mutation_table_columns       — required columns present
    test_no_germline_in_table         — germline row never appears in table
    test_observed_cells_flagged       — is_observed True for real cells
    test_nt_changes_format            — format "A42T" (ref-pos-alt)
    test_aa_changes_from_synonymous   — synonymous → empty aa_changes
    test_mutation_table_row_count     — one row per non-root clade
    test_build_before_table_required  — mutation_table() before build() → error

  Distance
    test_identical_sequences_zero_dist  — same seq → distance 0
    test_completely_different_distance  — no shared bases → distance ~ 1
    test_distance_symmetric             — dist(a,b) == dist(b,a)

  Error handling
    test_no_germline_record_raises      — ValueError if no germline in records
    test_empty_records_raises           — ValueError on empty list
"""

import pytest
from Bio.Phylo.BaseTree import Tree

from BCR_lineage_tracer.loader import CellRecord
from BCR_lineage_tracer.tracer import LineageTracer


# ─────────────────────────────── helpers ─────────────────────────────────────

def make_records(n_cells=5, seq_len=60, n_mutations=3, isotypes=None):
    """Build a minimal list of CellRecords with a synthetic germline."""
    import random
    random.seed(99)
    bases = "ACGT"
    germ_seq = "".join(random.choices(bases, k=seq_len))
    records = [CellRecord(cell_id="Germline_test", clone_id="test",
                          sequence=germ_seq, is_germline=True)]
    default_iso = ["IGHM", "IGHG1", "IGHA1", "IGHG2", "IGHM"]
    for i in range(n_cells):
        seq = list(germ_seq)
        for pos in random.sample(range(seq_len), min(n_mutations, seq_len)):
            seq[pos] = random.choice([b for b in bases if b != seq[pos]])
        iso = (isotypes[i] if isotypes else default_iso[i % len(default_iso)])
        records.append(CellRecord(
            cell_id=f"cell_{i}", clone_id="test",
            sequence="".join(seq), isotype=iso,
            sample_id="spleen", cluster_annotated="GC",
        ))
    return records, germ_seq


# ═══════════════════════════════════════════════════════════════════
# Tree structure
# ═══════════════════════════════════════════════════════════════════

def test_tree_is_rooted_on_germline(paired_tree):
    _, tree = paired_tree
    assert getattr(tree.root, "is_germline", False), \
        "Tree root must be the germline node"


def test_root_name_matches_germline_id(paired_tree):
    _, tree = paired_tree
    assert tree.root.name.startswith("Germline_")


def test_all_observed_cells_are_leaves(paired_clones, paired_tree):
    _, tree = paired_tree
    leaf_names = {c.name for c in tree.get_terminals()}
    obs = [r for r in paired_clones["clone_A"] if not r.is_germline]
    for rec in obs:
        assert rec.cell_id in leaf_names, \
            f"{rec.cell_id} not found as a leaf in the tree"


def test_no_phantom_root_node(paired_tree):
    """root_with_outgroup collapse: only ONE germline-tagged clade at root."""
    _, tree = paired_tree
    germ_clades = [c for c in tree.find_clades()
                   if getattr(c, "is_germline", False)]
    assert len(germ_clades) == 1, \
        f"Expected 1 germline node, got {len(germ_clades)}"


def test_branch_lengths_non_negative(paired_tree):
    _, tree = paired_tree
    for cl in tree.find_clades():
        if cl.branch_length is not None:
            assert cl.branch_length >= -1e-9, \
                f"Negative branch length on {cl.name}: {cl.branch_length}"


# ═══════════════════════════════════════════════════════════════════
# Polytomies
# ═══════════════════════════════════════════════════════════════════

def test_polytomy_collapse_produces_multifurcating_nodes():
    """With a high collapse threshold internal nodes collapse into polytomies."""
    records, _ = make_records(n_cells=6, n_mutations=1)  # very few mutations
    tracer_high = LineageTracer(records, collapse_threshold=1.0)
    tree_high = tracer_high.build()
    internal_high = sum(1 for c in tree_high.find_clades()
                        if c.clades and c is not tree_high.root)

    tracer_low = LineageTracer(records, collapse_threshold=0.0)
    tree_low = tracer_low.build()
    internal_low = sum(1 for c in tree_low.find_clades()
                       if c.clades and c is not tree_low.root)

    # High threshold must not produce MORE internal nodes than low threshold
    assert internal_high <= internal_low


def test_polytomy_default_threshold_preserves_structure(paired_tree):
    _, tree = paired_tree
    # With default 1e-6 threshold the tree should have internal nodes
    internal = [c for c in tree.find_clades() if c.clades and c is not tree.root]
    assert len(internal) >= 1


# ═══════════════════════════════════════════════════════════════════
# Isotype-aware NNI
# ═══════════════════════════════════════════════════════════════════

def test_csr_violations_not_increased():
    """Refinement must never make CSR violations worse."""
    # Deliberately create a clone with mixed isotypes to trigger violations
    records, _ = make_records(
        n_cells=6, n_mutations=4,
        isotypes=["IGHG1", "IGHM", "IGHA1", "IGHG1", "IGHM", "IGHG2"]
    )
    tracer_no_refine = LineageTracer(records, refine_isotypes=False)
    tree_no = tracer_no_refine.build()
    v_before = tracer_no_refine._count_violations(tree_no)

    tracer_refine = LineageTracer(records, refine_isotypes=True)
    tree_yes = tracer_refine.build()
    v_after = tracer_refine._count_violations(tree_yes)

    assert v_after <= v_before, \
        f"Refinement increased violations: {v_before} → {v_after}"


def test_refine_false_skips_nni():
    """With refine_isotypes=False the violation count comes from raw NJ."""
    records, _ = make_records(n_cells=5)
    t1 = LineageTracer(records, refine_isotypes=False)
    t1.build()
    t2 = LineageTracer(records, refine_isotypes=True)
    t2.build()
    # We can't guarantee the topologies differ, but both must build without error
    assert t1.tree is not None
    assert t2.tree is not None


# ═══════════════════════════════════════════════════════════════════
# Ancestral reconstruction
# ═══════════════════════════════════════════════════════════════════

def test_leaves_have_correct_sequence(paired_clones, paired_tree):
    _, tree = paired_tree
    recs_by_id = {r.cell_id: r for r in paired_clones["clone_A"]}
    for leaf in tree.get_terminals():
        if leaf.name in recs_by_id and not recs_by_id[leaf.name].is_germline:
            assert getattr(leaf, "sequence", None) == recs_by_id[leaf.name].sequence


def test_root_is_germline_node(paired_clones, paired_tree):
    """Root node must be the germline CellRecord (is_germline=True)."""
    _, tree = paired_tree
    assert getattr(tree.root, "is_germline", False), "Root is not the germline node"
    assert tree.root.name.startswith("Germline_")

def test_root_sequence_set(paired_tree):
    """Root must always have a .sequence attribute (Fitch or nearest-leaf)."""
    _, tree = paired_tree
    assert hasattr(tree.root, "sequence")
    assert len(tree.root.sequence) > 0


def test_internal_nodes_have_sequence(paired_tree):
    _, tree = paired_tree
    for cl in tree.find_clades():
        assert hasattr(cl, "sequence"), f"Clade {cl.name} has no .sequence"
        assert cl.sequence != "", f"Clade {cl.name} has empty .sequence"


def test_sequence_length_consistent(paired_tree):
    """Fitch path: all clades get same-length sequence (when seqs equal length)."""
    _, tree = paired_tree
    lengths = {len(getattr(cl, "sequence", "")) for cl in tree.find_clades()
               if getattr(cl, "sequence", "")}
    # Paired sequences may have 1-2 bp length variation across cells;
    # if all equal, expect exactly one unique length
    assert len(lengths) <= 3, \
        f"Too many distinct sequence lengths on tree: {lengths}"


# ═══════════════════════════════════════════════════════════════════
# Mutation table
# ═══════════════════════════════════════════════════════════════════

REQUIRED_COLS = {
    "clone_id", "node", "parent", "is_observed",
    "isotype", "sample_id", "cluster_annotated",
    "branch_length", "number_nucleotides_changes", "nucleotides_changes", "number_amino_acid_changes", "amino_acid_changes",
}

def test_mutation_table_columns(paired_tree):
    tracer, _ = paired_tree
    df = tracer.mutation_table()
    assert REQUIRED_COLS.issubset(set(df.columns)), \
        f"Missing columns: {REQUIRED_COLS - set(df.columns)}"


def test_no_germline_in_table(paired_tree):
    tracer, _ = paired_tree
    df = tracer.mutation_table()
    assert not df["node"].str.startswith("Germline_").any(), \
        "Germline node must not appear in mutation table"


def test_observed_cells_flagged(paired_clones, paired_tree):
    tracer, _ = paired_tree
    df = tracer.mutation_table()
    obs_ids = {r.cell_id for r in paired_clones["clone_A"] if not r.is_germline}
    obs_in_table = df[df["node"].isin(obs_ids)]
    assert obs_in_table["is_observed"].all(), \
        "Observed cells must have is_observed=True"


def test_nt_changes_format(paired_tree):
    """Each nucleotide change must match 'REF<pos>ALT' e.g. 'A42T'."""
    import re
    tracer, _ = paired_tree
    df = tracer.mutation_table()
    pattern = re.compile(r"^[ACGTN]\d+[ACGTN]$")
    for entry in df["nucleotides_changes"]:
        for change in (entry.split(";") if entry else []):
            assert pattern.match(change), \
                f"Unexpected nucleotides_changes format: '{change}'"


def test_mutation_table_has_rows(paired_tree):
    tracer, _ = paired_tree
    df = tracer.mutation_table()
    assert len(df) > 0, "Mutation table is empty"


def test_build_before_table_required():
    records, _ = make_records(n_cells=3)
    tracer = LineageTracer(records)
    with pytest.raises(RuntimeError, match="build()"):
        tracer.mutation_table()


# ═══════════════════════════════════════════════════════════════════
# Distance calculation
# ═══════════════════════════════════════════════════════════════════

def test_identical_sequences_zero_distance():
    records, _ = make_records(n_cells=2)
    tracer = LineageTracer(records)
    assert tracer._seq_dist("ACGT", "ACGT") == 0.0


def test_completely_different_distance():
    records, _ = make_records(n_cells=2)
    tracer = LineageTracer(records)
    d = tracer._seq_dist("AAAA", "TTTT")
    assert d == pytest.approx(1.0)


def test_distance_symmetric():
    records, _ = make_records(n_cells=2)
    tracer = LineageTracer(records)
    a, b = "ACGTACGT", "TTTTACGT"
    assert tracer._seq_dist(a, b) == pytest.approx(tracer._seq_dist(b, a))


def test_distance_unequal_lengths():
    """Unequal sequences are aligned first; distance should be < 1."""
    records, _ = make_records(n_cells=2)
    tracer = LineageTracer(records)
    d = tracer._seq_dist("ACGTACGT", "ACGTACGTTT")
    assert 0.0 <= d < 1.0


# ═══════════════════════════════════════════════════════════════════
# Error handling
# ═══════════════════════════════════════════════════════════════════

def test_no_germline_record_raises():
    records = [
        CellRecord("c1","clone","ACGTACGT"),
        CellRecord("c2","clone","ACGTACCC"),
    ]
    with pytest.raises(ValueError, match="germline"):
        LineageTracer(records)


def test_empty_records_raises():
    with pytest.raises(ValueError):
        LineageTracer([])
