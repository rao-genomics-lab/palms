"""Renaming a dataset without breaking its recorded provenance.

Nothing stores a dataset *name*, so the viewer opens a renamed dataset fine and
its content-hashed cache stays fresh. What breaks is the set of absolute paths
recorded inside the provenance graph and the session attrs — and it does not
break quietly: ``app.py`` re-emits the preamble for the current ``data_path`` on
every launch, and ``ProvGraph.upsert`` then flags every transitive descendant
stale. So the first launch after a hand-rolled ``mv`` marks the whole notebook
⚠. ``palms-rename-dataset`` exists to repair the graph before that happens.

The properties that matter, and that each test below pins:

* every recorded reference to the old directory is rewritten, in all four
  places it can hide (graph sidecar, session attr, session path attrs, CNV
  sidecars);
* a path pointing *outside* the dataset directory is left byte-identical — that
  file did not move;
* ``--dry-run`` writes nothing at all;
* unsafe situations refuse rather than half-apply.

Run headless:  QT_QPA_PLATFORM=offscreen pytest tests/test_rename_dataset.py
"""
from __future__ import annotations

import hashlib
import json
import sys
import warnings
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

pytest.importorskip("spatialdata")

from palms.scripts import rename_dataset as rd  # noqa: E402
from palms.scripts.rename_dataset import (  # noqa: E402
    PreflightError, infer_old_path, main, sub_text, sub_value,
)


@pytest.fixture(autouse=True)
def _quiet():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


# ── fixture: a dataset directory with every path-bearing artifact ────────────

def _graph_items(data_path: Path, outside: Path) -> list[dict]:
    """A graph shaped like a real one: a preamble plus a clustering read_csv.

    Those are the only two node kinds that carry an absolute path in practice —
    verified against a real store before this tool was written. The ``he:load``
    node carries a path *outside* the dataset, which must survive untouched.
    """
    return [
        {"id": "preamble", "code": f'data_path = Path(r"{data_path}")',
         "deps": [], "kind": "setup", "label": "Setup & data loading", "seq": 1},
        {"id": "clustering:leiden_r1.0",
         "code": (f'clust_df = pd.read_csv(r"{data_path}/analysis/clustering/'
                  f'leiden_r1.0/clusters.csv", index_col=0)'),
         "deps": ["preamble"], "kind": "artifact", "seq": 2},
        {"id": "he:load", "code": f'he = imread(r"{outside}/hande.ome.tif")',
         "deps": ["preamble"], "kind": "terminal", "seq": 3},
    ]


@pytest.fixture
def dataset(tmp_path, make_table):
    """A dataset directory: real zarr store, session attrs, sidecars, outputs."""
    np = pytest.importorskip("numpy")
    import zarr
    from spatialdata import SpatialData
    from spatialdata.models import Labels2DModel

    data_path = tmp_path / "old_name"
    data_path.mkdir()
    outside = tmp_path / "elsewhere"
    outside.mkdir()
    (data_path / "experiment.xenium").write_text(json.dumps({"pixel_size": 0.2125}))

    cache = data_path / "sdata_cached.zarr"
    labels = Labels2DModel.parse(np.arange(16, dtype=np.int32).reshape(4, 4))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        SpatialData(labels={"lab": labels},
                    tables={"table": make_table("OLD")}).write(cache)

    items = _graph_items(data_path, outside)
    session = zarr.open_group(str(cache / "viewer_session"), mode="a",
                              use_consolidated=False)
    session.attrs["prov_graph"] = items
    session.attrs["he_path"] = f"{data_path}/hande.ome.tif"
    session.attrs["arms_geojson_path"] = f"{outside}/tiles.geojson"
    session.attrs["arms_he_path"] = None
    session.attrs["roi_count"] = 2
    from palms.utils.zarr_safe import consolidate
    consolidate(cache)

    sidecars = data_path / "viewer_cache"
    sidecars.mkdir()
    (sidecars / "prov_graph.json").write_text(json.dumps(items))
    (sidecars / "cnv_copykat_result.json").write_text(json.dumps({
        "backend": "copykat",
        "adata_cnv_path": f"{data_path}/viewer_cache/adata_cnv_cache_copykat.h5ad",
        "heatmap_png": f"{data_path}/plots/cnv_heatmap_copykat_x.png",
        "n_cells": 100,
    }))

    (data_path / "analysis.py").write_text("# stale\n")
    from palms.utils.notebook_export import write_notebook
    write_notebook([("code", "# stale")], data_path / "analysis_notebook.ipynb")

    return {"data_path": data_path, "outside": outside, "cache": cache,
            "items": items, "tmp": tmp_path}


def _read_session_attrs(cache: Path) -> dict:
    from palms.utils.session import _read_prev_attrs
    return _read_prev_attrs(cache)


def _digest(path: Path) -> dict[str, str]:
    """Content hash of every file under *path*, for proving nothing was written."""
    return {
        str(p.relative_to(path)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(path.rglob("*")) if p.is_file()
    }


# ── the substitution rule ────────────────────────────────────────────────────

def test_a_path_prefix_is_rewritten_but_a_longer_name_is_not():
    # /d/foo must not match inside /d/foobar or /d/foo.bak: those are different
    # directories that did not move.
    assert sub_text('Path(r"/d/foo")', "/d/foo", "/d/bar") == ('Path(r"/d/bar")', 1)
    assert sub_text('r"/d/foo/a.csv"', "/d/foo", "/d/bar") == ('r"/d/bar/a.csv"', 1)
    assert sub_text('"/d/foobar/x"', "/d/foo", "/d/bar") == ('"/d/foobar/x"', 0)
    assert sub_text('"/d/foo.bak/x"', "/d/foo", "/d/bar") == ('"/d/foo.bak/x"', 0)


def test_sub_value_recurses_and_leaves_non_strings_alone():
    payload = {"a": "/d/foo/x", "b": 3, "c": [{"d": "/d/foo"}, None]}
    fixed, n = sub_value(payload, "/d/foo", "/d/bar")
    assert n == 2
    assert fixed == {"a": "/d/bar/x", "b": 3, "c": [{"d": "/d/bar"}, None]}


def test_the_old_path_is_inferred_from_the_preamble(dataset):
    assert infer_old_path(dataset["items"]) == str(dataset["data_path"])
    assert infer_old_path([{"id": "normalize", "code": "sc.pp.log1p(adata)"}]) is None


# ── the move ─────────────────────────────────────────────────────────────────

def test_rename_rewrites_every_recorded_reference(dataset):
    data_path, outside = dataset["data_path"], dataset["outside"]
    new_path = dataset["tmp"] / "new_name"

    assert main([str(data_path), "new_name"]) == 0

    assert not data_path.exists()
    assert (new_path / "experiment.xenium").exists()

    sidecar = json.loads((new_path / "viewer_cache" / "prov_graph.json").read_text())
    by_id = {n["id"]: n for n in sidecar}
    assert by_id["preamble"]["code"] == f'data_path = Path(r"{new_path}")'
    assert str(new_path) in by_id["clustering:leiden_r1.0"]["code"]
    assert str(data_path) not in json.dumps(sidecar)

    attrs = _read_session_attrs(new_path / "sdata_cached.zarr")
    assert str(data_path) not in json.dumps(attrs["prov_graph"])
    assert attrs["he_path"] == f"{new_path}/hande.ome.tif"
    assert attrs["roi_count"] == 2          # untouched attrs survive the swap

    cnv = json.loads((new_path / "viewer_cache" / "cnv_copykat_result.json").read_text())
    assert cnv["adata_cnv_path"].startswith(str(new_path))
    assert cnv["n_cells"] == 100

    # The one path that did not move stays exactly as recorded.
    assert by_id["he:load"]["code"] == f'he = imread(r"{outside}/hande.ome.tif")'
    assert attrs["arms_geojson_path"] == f"{outside}/tiles.geojson"


def test_a_full_destination_path_moves_rather_than_renames(dataset):
    dest = dataset["tmp"] / "somewhere" / "else"
    dest.parent.mkdir()
    assert main([str(dataset["data_path"]), str(dest)]) == 0
    assert (dest / "experiment.xenium").exists()
    items = json.loads((dest / "viewer_cache" / "prov_graph.json").read_text())
    assert infer_old_path(items) == str(dest)


def test_derived_outputs_are_regenerated_from_the_repaired_graph(dataset):
    new_path = dataset["tmp"] / "new_name"
    assert main([str(dataset["data_path"]), "new_name"]) == 0

    code = (new_path / "analysis.py").read_text()
    assert f'data_path = Path(r"{new_path}")' in code
    assert "# stale" not in code

    from palms.utils.notebook_export import read_notebook
    sources = "\n".join(src for _kind, src in read_notebook(new_path / NB))
    assert str(new_path) in sources
    assert str(dataset["data_path"]) not in sources


NB = "analysis_notebook.ipynb"


def test_derived_outputs_are_not_created_when_absent(dataset):
    (dataset["data_path"] / "analysis.py").unlink()
    (dataset["data_path"] / NB).unlink()
    new_path = dataset["tmp"] / "new_name"
    assert main([str(dataset["data_path"]), "new_name"]) == 0
    assert not (new_path / "analysis.py").exists()
    assert not (new_path / NB).exists()


def test_derived_outputs_are_replaced_not_written_in_place(dataset):
    """Every write this tool makes must swap an inode, never truncate one.

    Found by probing the real thing against a ``cp -al`` snapshot: the two
    derived files were the only writes not going through ``atomic_json``, so
    regenerating them wrote *through* the hardlink and silently edited the
    snapshot as well. The property is cheap to state and pins the whole class:
    a hardlink made before the rename must still hold the old bytes after it.
    """
    import os
    data_path, tmp = dataset["data_path"], dataset["tmp"]
    snapshot = tmp / "snapshot"
    snapshot.mkdir()
    for name in ("analysis.py", NB):
        os.link(data_path / name, snapshot / name)
    before = {name: (snapshot / name).read_bytes() for name in ("analysis.py", NB)}

    assert main([str(data_path), "new_name"]) == 0

    for name in ("analysis.py", NB):
        assert (snapshot / name).read_bytes() == before[name], (
            f"{name} was written in place, mutating a hardlinked snapshot"
        )
    # ...and the real target really did change.
    assert (tmp / "new_name" / "analysis.py").read_bytes() != before["analysis.py"]


# ── repair mode ──────────────────────────────────────────────────────────────

def test_repair_reports_a_repair_rather_than_a_move(dataset, capsys):
    import os
    new_path = dataset["tmp"] / "moved_by_hand"
    os.rename(dataset["data_path"], new_path)
    assert main([str(new_path), "--repair"]) == 0
    out = capsys.readouterr().out
    assert "Repaired:" in out and "Moved" not in out
    assert "recorded path was:" in out


def test_repair_fixes_a_dataset_moved_by_hand(dataset):
    import os
    data_path = dataset["data_path"]
    new_path = dataset["tmp"] / "moved_by_hand"
    os.rename(data_path, new_path)          # what a user would do

    assert main([str(new_path), "--repair"]) == 0

    items = json.loads((new_path / "viewer_cache" / "prov_graph.json").read_text())
    assert infer_old_path(items) == str(new_path)


def test_repair_accepts_an_explicit_from_path(dataset):
    import os
    new_path = dataset["tmp"] / "moved_by_hand"
    os.rename(dataset["data_path"], new_path)
    assert main([str(new_path), "--repair", "--from", str(dataset["data_path"])]) == 0
    items = json.loads((new_path / "viewer_cache" / "prov_graph.json").read_text())
    assert infer_old_path(items) == str(new_path)


def test_repair_is_a_no_op_when_the_path_already_matches(dataset, capsys):
    assert main([str(dataset["data_path"]), "--repair"]) == 0
    assert "Nothing to do" in capsys.readouterr().out


def test_repair_without_a_preamble_refuses_rather_than_guessing(dataset, capsys):
    sidecar = dataset["data_path"] / "viewer_cache" / "prov_graph.json"
    sidecar.write_text(json.dumps([{"id": "normalize", "code": "sc.pp.log1p(adata)"}]))
    import os
    new_path = dataset["tmp"] / "moved"
    os.rename(dataset["data_path"], new_path)
    # The session attr still has a preamble, so blank that too.
    from palms.utils.zarr_safe import safe_group_update
    with safe_group_update(new_path / "sdata_cached.zarr", "viewer_session") as (g, _):
        g.attrs["prov_graph"] = []
    assert main([str(new_path), "--repair"]) == 1
    assert "could not infer" in capsys.readouterr().err


# ── dry run ──────────────────────────────────────────────────────────────────

def test_dry_run_writes_nothing(dataset, capsys):
    data_path = dataset["data_path"]
    before = _digest(dataset["tmp"])
    assert main([str(data_path), "new_name", "--dry-run"]) == 0
    assert _digest(dataset["tmp"]) == before
    assert data_path.exists()
    out = capsys.readouterr().out
    assert "Would move" in out and "Would rewrite" in out
    assert "Would regenerate analysis.py" in out


# ── refusals ─────────────────────────────────────────────────────────────────

def test_an_existing_destination_is_refused(dataset, capsys):
    (dataset["tmp"] / "taken").mkdir()
    assert main([str(dataset["data_path"]), "taken"]) == 1
    assert "already exists" in capsys.readouterr().err
    assert dataset["data_path"].exists()


def test_a_running_copykat_job_is_refused(dataset, capsys):
    plots = dataset["data_path"] / "plots"
    plots.mkdir()
    (plots / "copykat_RUNNING.txt").write_text("12345")
    assert main([str(dataset["data_path"]), "new_name"]) == 1
    assert "CopyKAT" in capsys.readouterr().err
    assert dataset["data_path"].exists()


def test_a_non_dataset_directory_is_refused(tmp_path, capsys):
    (tmp_path / "random").mkdir()
    assert main([str(tmp_path / "random"), "new_name"]) == 1
    assert "does not look like a Xenium dataset" in capsys.readouterr().err


def test_a_cache_only_export_is_accepted(tmp_path, dataset):
    # A Crop Dataset export has no raw Xenium files; its zarr *is* the data.
    (dataset["data_path"] / "experiment.xenium").unlink()
    assert rd.is_dataset_dir(dataset["data_path"])
    assert main([str(dataset["data_path"]), "cropped"]) == 0


def test_unreadable_root_metadata_is_refused(dataset, capsys):
    (dataset["cache"] / "zarr.json").write_text("{ not json")
    assert main([str(dataset["data_path"]), "new_name"]) == 1
    assert "unreadable" in capsys.readouterr().err
    assert dataset["data_path"].exists()


def test_a_cross_device_move_refuses_instead_of_copying(dataset, monkeypatch, capsys):
    import errno as _errno
    import os as _os

    def _boom(src, dst):
        raise OSError(_errno.EXDEV, "Invalid cross-device link")

    monkeypatch.setattr(_os, "rename", _boom)
    assert main([str(dataset["data_path"]), "new_name"]) == 1
    err = capsys.readouterr().err
    assert "different filesystems" in err and "--repair" in err
    assert dataset["data_path"].exists()


def test_naming_both_a_destination_and_repair_is_an_error(dataset, capsys):
    assert main([str(dataset["data_path"]), "new_name", "--repair"]) == 2
    assert "not both" in capsys.readouterr().err


def test_neither_a_destination_nor_repair_is_an_error(dataset, capsys):
    assert main([str(dataset["data_path"])]) == 2
    assert "--repair" in capsys.readouterr().err


# ── source guard ─────────────────────────────────────────────────────────────

def test_the_store_is_only_reached_through_the_safe_write_helpers():
    """A rename must not become another way to write the zarr directly.

    ``safe_group_update`` and ``atomic_json`` are the sanctioned paths; a bare
    ``zarr.open_group`` write, a ``shutil.move`` of the store, or an ``rmtree``
    would each reintroduce a loss window this codebase spent real effort
    closing. Reading the session attrs goes through ``session._read_prev_attrs``
    for the same reason.

    Parsed rather than grepped: this module *documents* why it does not call
    ``shutil.move``, and a text search cannot tell an explanation from a call.
    """
    import ast

    forbidden = {
        "shutil.move", "shutil.rmtree", "shutil.copytree", "os.remove",
        "os.unlink", "zarr.open", "zarr.open_group", "zarr.consolidate_metadata",
    }
    tree = ast.parse(Path(rd.__file__).read_text())
    called = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            called.add(ast.unparse(node.func))
    assert not (called & forbidden), f"forbidden call(s): {sorted(called & forbidden)}"
    # Path.unlink()/rmdir() on any receiver: this tool never deletes.
    assert not {c for c in called if c.split(".")[-1] in ("unlink", "rmdir", "rmtree")}

    source = Path(rd.__file__).read_text()
    assert "safe_group_update" in source
    assert "atomic_json" in source
