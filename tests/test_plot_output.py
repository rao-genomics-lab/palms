"""Where a figure goes, decided in one place.

Before ``utils/plot_output``, five modules built ``<data_path>/plots`` by hand,
four save policies coexisted, and the ``savefig`` path a ``plot:*`` node recorded
was a bare relative guess that matched none of them. These tests pin the
properties that replaced all of that.

Pure stdlib plus matplotlib's Agg canvas: no Qt, no napari, no dataset.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from palms.utils import plot_output as po  # noqa: E402


# ── formats ──────────────────────────────────────────────────────────────────

def test_the_default_writes_both_a_png_and_a_pdf():
    """Issue #35's literal ask, and the state app.py seeds."""
    assert po.plot_formats({}) == ["png", "pdf"]


def test_an_explicit_choice_is_honoured():
    assert po.plot_formats({"plot_formats": ["svg"]}) == ["svg"]
    assert po.plot_formats({"plot_formats": ["png", "pdf"]}) == ["png", "pdf"]


def test_the_old_scalar_setting_still_reads():
    """``plot_format: "svg"`` is what a session saved before this change."""
    assert po.plot_formats({"plot_formats": "svg"}) == ["svg"]


@pytest.mark.parametrize("state", [
    {},                                  # missing
    {"plot_formats": []},                # emptied
    {"plot_formats": ["tiff", "eps"]},   # nothing we can write
])
def test_an_unusable_setting_falls_back_rather_than_writing_nothing(state):
    """A plot the user asked for and cannot find afterwards is the failure this
    module exists to remove — so never resolve to an empty format list, and
    never hand ``savefig`` an extension that would raise inside a worker."""
    assert po.plot_formats(state) == ["png", "pdf"]


def test_primary_format_is_what_a_single_name_would_be():
    assert po.primary_format({"plot_formats": ["pdf", "png"]}) == "pdf"


# ── names and paths ──────────────────────────────────────────────────────────

@pytest.mark.parametrize("raw,expected", [
    ("leiden_r1.0", "leiden_r1.0"),
    ("leiden (res 1.0)", "leiden_res_1.0"),
    ("a/b\\c", "a_b_c"),
    ("  ", "plot"),          # never an empty filename
    ("", "plot"),
])
def test_safe_stem_makes_a_clustering_key_safe_to_put_in_a_filename(raw, expected):
    assert po.safe_stem(raw) == expected


def test_plots_dir_does_not_create_unless_asked(tmp_path):
    """Read-only by default: a preview provider asks where a figure *would* go,
    and drawing a preview pane must not touch the disk."""
    directory = po.plots_dir(tmp_path)
    assert directory == tmp_path / "plots"
    assert not directory.exists()

    assert po.plots_dir(tmp_path, create=True).exists()


def test_save_paths_gives_one_file_per_format(tmp_path):
    paths = po.save_paths(tmp_path, "dotplot_leiden", state={})
    assert [p.name for p in paths] == ["dotplot_leiden.png", "dotplot_leiden.pdf"]
    assert all(p.parent == tmp_path / "plots" for p in paths)


def test_batch_dir_keeps_a_volcano_run_under_plots(tmp_path):
    assert po.batch_dir(tmp_path, "volcano_leiden r1.0") == (
        tmp_path / "plots" / "volcano_leiden_r1.0")


# ── saving ───────────────────────────────────────────────────────────────────

def _figure():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    return fig


def test_save_figure_writes_every_path_and_creates_the_directory(tmp_path):
    fig = _figure()
    paths = po.save_paths(tmp_path, "demo", state={})
    written = po.save_figure(fig, paths)

    assert len(written) == 2
    for path in paths:
        assert path.exists() and path.stat().st_size > 0
    # The directory did not exist beforehand; save_figure is what made it.
    assert (tmp_path / "plots").is_dir()


# ── what gets recorded ───────────────────────────────────────────────────────

def test_recorded_paths_are_relative_to_the_dataset(tmp_path):
    """An exported notebook replays with ``data_path`` bound and its cwd there,
    so an absolute path would pin the notebook to one machine."""
    paths = po.save_paths(tmp_path, "demo", state={})
    assert po.recorded_paths(tmp_path, paths) == [
        "plots/demo.png", "plots/demo.pdf"]


def test_a_path_outside_the_dataset_stays_absolute(tmp_path):
    """A volcano batch the user redirected elsewhere did not move into the
    dataset; recording it as if it had would name a file that is not there."""
    outside = tmp_path.parent / "elsewhere" / "v.png"
    assert po.recorded_paths(tmp_path, [outside]) == [str(outside)]


def test_recorded_save_code_creates_the_directory_before_writing(tmp_path):
    """``savefig`` does not make its parent, and a replayed notebook has no
    ``plots/`` until something does — so the recorded cell has to."""
    code = po.recorded_save_code(["plots/demo.png", "plots/demo.pdf"])
    assert code.count("mkdir(parents=True, exist_ok=True)") == 2

    fig = _figure()
    namespace = {"fig": fig, "Path": Path}
    import os
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        exec(compile(code, "<recorded>", "exec"), namespace)  # noqa: S102
    finally:
        os.chdir(cwd)
    assert (tmp_path / "plots" / "demo.png").exists()
    assert (tmp_path / "plots" / "demo.pdf").exists()


def test_recorded_save_code_can_name_a_different_figure_expression():
    """The rank-genes dotplot binds ``dotplot``, not ``fig``."""
    code = po.recorded_save_code(["plots/d.png"], fig_expr="dotplot")
    assert "dotplot.savefig(" in code
