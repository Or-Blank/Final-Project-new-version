"""
test_visualization.py
=====================
Tests for plot_tree() — checks that figures are produced without errors
and that the legend / axes reflect the correct metadata.

What each test covers
---------------------
  test_plot_returns_fig_ax       — returns (Figure, Axes) tuple
  test_output_png_written        — file written when output_path given
  test_color_by_isotype          — legend contains isotype labels
  test_color_by_sample_id        — legend contains tissue labels
  test_no_shape_by               — shape_by=None doesn't crash
  test_empty_shape_by_none       — shape_by=None still labels x-axis
  test_title_in_figure           — title string appears on axes
  test_all_leaves_have_text      — leaf labels rendered (text artists > 0)
"""

import os
import pytest
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from BCR_lineage_tracer.visualization import plot_tree


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


def test_plot_returns_fig_ax(paired_tree):
    _, tree = paired_tree
    result = plot_tree(tree, title="test")
    assert len(result) == 3
    fig, ax, _ = result
    import matplotlib.figure, matplotlib.axes
    assert isinstance(fig, matplotlib.figure.Figure)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_output_png_written(paired_tree, tmp_path):
    _, tree = paired_tree
    out = str(tmp_path / "tree.png")
    plot_tree(tree, output_path=out)
    assert os.path.exists(out)
    assert os.path.getsize(out) > 1000   # non-trivial PNG


def test_color_by_isotype(paired_tree):
    _, tree = paired_tree
    fig, ax, _ = plot_tree(tree, color_by="isotype")
    legend_texts = [t.get_text() for t in ax.get_legend().get_texts()]
    assert any("isotype" in t for t in legend_texts)


def test_color_by_sample_id(heavy_tree):
    _, tree = heavy_tree
    fig, ax, _ = plot_tree(tree, color_by="sample_id")
    legend_texts = [t.get_text() for t in ax.get_legend().get_texts()]
    assert any("sample_id" in t for t in legend_texts)


def test_no_shape_by_does_not_crash(paired_tree):
    _, tree = paired_tree
    fig, ax, _ = plot_tree(tree, shape_by=None)
    assert fig is not None


def test_title_in_figure(paired_tree):
    _, tree = paired_tree
    fig, ax, _ = plot_tree(tree, title="My Clone ABC")
    assert "My Clone ABC" in ax.get_title()


def test_x_axis_label(paired_tree):
    _, tree = paired_tree
    _, ax, _ = plot_tree(tree)
    assert "mutation" in ax.get_xlabel().lower()


def test_leaf_text_artists_present(paired_tree):
    _, tree = paired_tree
    fig, ax, _ = plot_tree(tree)
    texts = [t.get_text() for t in ax.texts]
    assert len(texts) > 0, "No text labels found on tree plot"


def test_germline_included_in_plot(paired_tree):
    _, tree = paired_tree
    fig, ax, _ = plot_tree(tree)
    labels = [t.get_text().strip() for t in ax.texts]
    assert any("Germline" in l for l in labels)
