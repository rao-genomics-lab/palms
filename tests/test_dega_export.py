"""The DegaFile export wrapper: where it writes, and what it refuses to write to.

``utils/dega_export`` exists for one reason — ``celldega.pre.main`` decompresses
several hundred MB *into the directory it is pointed at*, and pointing it at raw
10x output is therefore not an option. Everything measured here is some form of
that promise: the raw output is read through symlinks and comes back
byte-for-byte unchanged, and every write lands under ``viewer_cache/``.

**Nothing here requires celldega.** The wrapper's own logic — path derivation,
the symlink farm, the extraction that exists because ``gzip -dk`` refuses a
symlink — is pure stdlib and is tested for real. The one function that does
reach for celldega is tested against a stub module, which is stronger than
skipping: it pins the arguments the wrapper passes, which is exactly the thing
that would drift when the pin moves. The only genuinely skipped test is the one
that asserts what the *absence* of celldega looks like.
"""

from __future__ import annotations

import ast
import gzip
import importlib.util
import os
import sys
import tarfile
import zipfile
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from palms.utils import dega_export as de  # noqa: E402

CELLDEGA_INSTALLED = importlib.util.find_spec("celldega") is not None


# ── fixtures ─────────────────────────────────────────────────────────────────

@pytest.fixture
def dataset(tmp_path):
    """A dataset directory shaped like the parts of a Xenium bundle that matter.

    Only the four archives ``celldega.pre`` unpacks plus enough neighbours to
    tell "raw output" from "ours" — no real formats, since nothing here decodes
    anything but the compression.
    """
    data_path = tmp_path / "Xenium_V1_test_outs"
    data_path.mkdir()
    (data_path / "experiment.xenium").write_text('{"run_name": "test"}')
    (data_path / "transcripts.parquet").write_bytes(b"PAR1" * 32)
    (data_path / "morphology_focus").mkdir()
    (data_path / "morphology_focus" / "focus_0000.ome.tif").write_bytes(b"II*\0")

    # The four archives, each really compressed.
    with gzip.open(data_path / "cells.csv.gz", "wb") as fh:
        fh.write(b"cell_id,x,y\n1,0.0,0.0\n")
    with zipfile.ZipFile(data_path / "cells.zarr.zip", "w") as zf:
        zf.writestr(".zgroup", '{"zarr_format": 2}')
    for name in ("analysis", "cell_feature_matrix"):
        payload = tmp_path / name
        payload.mkdir()
        (payload / "contents.txt").write_text(name)
        with tarfile.open(data_path / f"{name}.tar.gz", "w:gz") as tf:
            tf.add(payload, arcname=name)

    # The viewer's own directories, which must not be linked into the farm.
    (data_path / "viewer_cache").mkdir()
    (data_path / "plots").mkdir()
    (data_path / "transcript_cache").mkdir()
    (data_path / "sdata_cached.zarr").mkdir()
    return data_path


def _snapshot(root: Path) -> dict:
    """Every file under *root*, by relative path, with its size and mtime."""
    out = {}
    for path in sorted(root.rglob("*")):
        if path.is_dir():
            continue
        stat = path.lstat()
        out[str(path.relative_to(root))] = (stat.st_size, stat.st_mtime_ns)
    return out


# ── where things go ──────────────────────────────────────────────────────────

def test_the_export_goes_beside_the_data_and_not_into_the_cache(dataset):
    """A DegaFile set is a deliverable, not a cache the viewer may clear.

    It has to sit outside ``deletable_roots()`` for the same reason ``plots/``
    does: a cache rebuild deletes what is inside them, and nobody expects a
    rebuild to throw away the thing they published.
    """
    from palms.utils.store_inventory import deletable_roots

    out = de.degafiles_dir(dataset)
    assert out == dataset / "degafiles"

    roots = [Path(r).resolve() for r in deletable_roots(dataset, dataset / "sdata_cached.zarr")]
    assert not any(out.resolve() == r or r in out.resolve().parents for r in roots), (
        "the published DegaFiles are inside a directory the viewer may delete"
    )


def test_a_published_export_is_reported_as_viewer_output_and_not_as_raw_data(dataset):
    """220 MB has to land in the right section of Tools -> Dataset.

    ``store_inventory`` defaults an unrecognised entry in the dataset directory
    to raw 10x output, which is the safe default and, for this one, a lie: it
    would be listed as "original 10x output, never modified by the viewer" and
    the size report would stop adding up. Listed, blocked, and told to be
    removed by hand — exactly how ``plots/`` is treated, and for the same reason.
    """
    from palms.utils import store_inventory as si

    de.degafiles_dir(dataset).mkdir()
    (de.degafiles_dir(dataset) / "landscape_parameters.json").write_text("{}")

    nodes = {n.key: n for section in si.build_inventory(dataset, None)
             for n in section.nodes}
    node = nodes["derived:degafiles"]
    assert not node.deletable and node.blocked_reason
    assert "raw" not in (node.detail or "").lower()
    assert "degafiles" not in {n.name.rstrip("/") for n in nodes.values()
                               if n.kind == si.RAW}


def test_the_staging_farm_is_reported_as_the_viewers_own(dataset):
    """Everything in ``viewer_cache/`` is viewer output, so this needs no rule.

    Asserted anyway: the farm is where a stalled export leaves its debris, and
    a user who wants the disk back has to be able to find it.
    """
    from palms.utils import store_inventory as si

    de.extract_archives(de.stage_raw_output(dataset))
    nodes = {n.key: n for section in si.build_inventory(dataset, None)
             for n in section.nodes}
    node = nodes["sidecar:dega_staging"]
    assert node.deletable and node.path == dataset / "viewer_cache" / "dega_staging"
    si.assert_deletable(node.path, si.deletable_roots(dataset, None))


def test_the_staging_farm_is_named_after_the_dataset(dataset):
    """celldega's contract is a sample *name* plus its parent directory.

    That name is what it stamps into ``landscape_parameters.json``, so a flat
    ``dega_staging`` directory would publish DegaFiles whose sample is called
    "dega_staging". Hence the extra level.
    """
    staging = de.staging_dir(dataset)
    assert staging.name == dataset.name
    assert staging.parent == dataset / "viewer_cache" / "dega_staging"


# ── the symlink farm ─────────────────────────────────────────────────────────

def test_staging_links_the_raw_entries_and_leaves_them_alone(dataset):
    before = _snapshot(dataset)
    staging = de.stage_raw_output(dataset)

    for name in ("experiment.xenium", "transcripts.parquet", "cells.csv.gz",
                 "cells.zarr.zip", "analysis.tar.gz", "morphology_focus"):
        link = staging / name
        assert link.is_symlink(), f"{name} was copied, or not staged at all"
        assert link.readlink() == dataset / name

    after = _snapshot(dataset)
    raw_before = {k: v for k, v in before.items() if not k.startswith("viewer_cache/")}
    raw_after = {k: v for k, v in after.items() if not k.startswith("viewer_cache/")}
    assert raw_before == raw_after, "staging modified the raw output"


def test_staging_does_not_link_the_viewers_own_directories(dataset):
    """``viewer_cache`` above all: the farm lives inside it.

    Linking it in would put the staging directory inside itself, and every tool
    that walks the farm would follow the loop.
    """
    staging = de.stage_raw_output(dataset)
    staged = {p.name for p in staging.iterdir()}
    assert not staged & {"viewer_cache", "plots", "transcript_cache",
                         "sdata_cached.zarr", "degafiles"}


def test_staging_is_idempotent_and_refreshes_a_link_that_moved(dataset):
    staging = de.stage_raw_output(dataset)
    (staging / "experiment.xenium").unlink()
    (staging / "experiment.xenium").symlink_to(dataset / "transcripts.parquet")

    de.stage_raw_output(dataset)
    assert (staging / "experiment.xenium").readlink() == dataset / "experiment.xenium"


def test_staging_keeps_what_a_previous_export_extracted(dataset):
    """A real file in the farm is celldega's output, and re-doing it is minutes."""
    staging = de.stage_raw_output(dataset)
    de.extract_archives(staging)
    extracted = staging / "cells.csv"
    assert extracted.is_file() and not extracted.is_symlink()
    stamp = extracted.stat().st_mtime_ns

    de.stage_raw_output(dataset)
    assert extracted.is_file() and not extracted.is_symlink()
    assert extracted.stat().st_mtime_ns == stamp


# ── extraction ───────────────────────────────────────────────────────────────

def test_extraction_reads_through_the_symlinks_that_gzip_refuses(dataset):
    """The reason :func:`extract_archives` exists at all.

    ``gzip -dk`` on a symlink is "not a regular file", exit 1 — so linking the
    archives and letting celldega's own unzipper run does not merely put the
    output in the wrong place, it does not work. Doing it in Python also drops
    their unstated dependency on the gzip, unzip and tar binaries.
    """
    staging = de.stage_raw_output(dataset)
    assert (staging / "cells.csv.gz").is_symlink()

    produced = de.extract_archives(staging)
    assert set(produced) == {"cells.csv", "cells.zarr", "analysis",
                             "cell_feature_matrix"}
    assert (staging / "cells.csv").read_text().startswith("cell_id,")
    assert (staging / "cells.zarr" / ".zgroup").exists()
    assert (staging / "analysis" / "contents.txt").read_text() == "analysis"


def test_extraction_lands_in_the_farm_and_not_beside_the_data(dataset):
    before = _snapshot(dataset)
    staging = de.stage_raw_output(dataset)
    de.extract_archives(staging)

    raw = {k: v for k, v in _snapshot(dataset).items()
           if not k.startswith("viewer_cache/")}
    assert raw == {k: v for k, v in before.items()
                   if not k.startswith("viewer_cache/")}
    for leaked in ("cells.csv", "cells.zarr", "analysis", "cell_feature_matrix"):
        assert not (dataset / leaked).exists(), (
            f"{leaked} was decompressed into the raw 10x output"
        )


def test_a_second_extraction_is_a_no_op(dataset):
    staging = de.stage_raw_output(dataset)
    de.extract_archives(staging)
    assert de.extract_archives(staging) == [], (
        "re-extracting would undo the saving that makes a second export cheap"
    )


def test_a_tar_member_cannot_escape_the_staging_directory(dataset, tmp_path):
    """``filter='data'`` is the default from 3.14 and a warning before it.

    Named here rather than left to the default, because the archives come from
    the user's dataset directory and a viewer that unpacks them is a viewer that
    can be handed a hostile one.
    """
    escape = tmp_path / "escape.txt"
    escape.write_text("should not be written")
    (dataset / "analysis.tar.gz").unlink()
    with tarfile.open(dataset / "analysis.tar.gz", "w:gz") as tf:
        tf.add(escape, arcname="../../escaped.txt")

    staging = de.stage_raw_output(dataset)
    with pytest.raises(tarfile.TarError):
        de.extract_archives(staging)
    assert not (dataset / "viewer_cache" / "escaped.txt").exists()
    assert not (dataset.parent / "escaped.txt").exists()


def test_clearing_staging_removes_the_farm_and_nothing_else(dataset):
    staging = de.stage_raw_output(dataset)
    de.extract_archives(staging)
    before = _snapshot(dataset)

    de.clear_staging(dataset)

    assert not (dataset / "viewer_cache" / "dega_staging").exists()
    assert (dataset / "viewer_cache").is_dir()
    raw = {k: v for k, v in before.items() if not k.startswith("viewer_cache/")}
    assert raw == {k: v for k, v in _snapshot(dataset).items()
                   if not k.startswith("viewer_cache/")}


# ── a dataset with nothing to export ─────────────────────────────────────────

@pytest.fixture
def crop_export(tmp_path):
    """A Crop Dataset export: experiment.xenium, a zarr cache, transcripts.

    None of celldega's raw inputs, because a crop export does not have them —
    its SpatialData zarr *is* the data. Shaped after the real ``crop_6`` demo
    dataset, where this was measured.
    """
    data_path = tmp_path / "crop_6"
    data_path.mkdir()
    (data_path / "experiment.xenium").write_text('{"run_name": "crop_6"}')
    (data_path / "transcripts.parquet").write_bytes(b"PAR1" * 32)
    (data_path / "sdata_cached.zarr").mkdir()
    (data_path / "transcript_cache").mkdir()
    return data_path


def test_a_crop_export_is_refused_before_anything_is_staged(crop_export):
    """The failure without this guard is late, misleading, and leaves debris.

    Measured on the real crop_6 demo dataset: celldega gets as far as its own
    unzipper and dies with ``CalledProcessError: Command '['gzip', '-dk',
    'cells.csv.gz']' returned non-zero exit status 1``, which reads as a broken
    gzip rather than as "this dataset has no raw output and never will". Issue
    #17's rule, applied to the one path that reads the raw output *instead of*
    the cache.
    """
    assert de.is_cache_only(crop_export)

    with pytest.raises(de.NotExportable) as excinfo:
        de.require_exportable(crop_export)
    message = str(excinfo.value)
    assert "crop_6" in message and "Crop Dataset export" in message
    # It has to say what to do instead, or it is only a better-worded dead end.
    assert "cropped from" in message

    # Nothing tried, nothing left behind.
    assert not (crop_export / "viewer_cache").exists()
    assert not de.degafiles_dir(crop_export).exists()


def test_the_export_refuses_a_crop_before_it_even_looks_for_celldega(crop_export,
                                                                    stub_celldega):
    """Order matters: the unsatisfiable refusal comes first.

    Asking "is celldega installed?" of a dataset that could never be published
    either way sends the user to install a dependency that will not help.
    """
    with pytest.raises(de.NotExportable):
        de.export_degafiles(crop_export)
    assert stub_celldega == [], "celldega was called for a dataset with no raw output"
    assert not de.degafiles_dir(crop_export).exists()


def test_a_partial_bundle_is_a_different_message_from_a_crop(dataset):
    """One missing archive is a broken copy, not a crop, and is said so.

    Collapsing the two would tell someone whose download was truncated to go and
    publish "the dataset this one was cropped from", which does not exist.
    """
    (dataset / "cells.csv.gz").unlink()
    assert not de.is_cache_only(dataset)

    with pytest.raises(de.NotExportable) as excinfo:
        de.require_exportable(dataset)
    message = str(excinfo.value)
    assert "cells.csv.gz" in message and "incomplete copy" in message
    # It names the crop case only to rule it out, and must not give the crop
    # advice — there is no dataset this one was cropped from.
    assert "cropped from" not in message
    assert "has no raw 10x output" not in message


def test_a_complete_bundle_is_allowed(dataset):
    assert de.missing_raw_inputs(dataset) == []
    de.require_exportable(dataset)  # must not raise


def test_the_publish_tab_says_so_on_arrival_for_a_crop(qapp, crop_export):
    """Stated when the tab is built, not discovered when the button fails.

    Filesystem-only, so it costs nothing at build time — which is the whole
    reason it can be answered this early.
    """
    from palms.tabs.tab_publish import build_tab
    from qtpy.QtWidgets import QLabel

    ctx = SimpleNamespace(state={}, viewer=None, data_path=crop_export)
    widget, _exports = build_tab(ctx)
    text = " ".join(lbl.text() for lbl in widget.findChildren(QLabel))
    assert "no raw 10x output" in text
    assert "Crop Dataset export" in text


def test_the_publish_tab_is_silent_about_it_for_a_normal_dataset(qapp, dataset):
    from palms.tabs.tab_publish import build_tab
    from qtpy.QtWidgets import QLabel

    ctx = SimpleNamespace(state={}, viewer=None, data_path=dataset)
    widget, _exports = build_tab(ctx)
    text = " ".join(lbl.text() for lbl in widget.findChildren(QLabel))
    assert "no raw 10x output" not in text


# ── the celldega call, against a stub ─────────────────────────────────────────

@pytest.fixture
def stub_celldega(monkeypatch):
    """Stand in for ``celldega.pre`` and ``pyvips``, recording the call.

    A stub and not a skip. What is worth pinning is the *arguments* — the
    sample/data_root_dir pair is celldega's own contract and the one thing that
    silently changes meaning when the pin moves — and running the real thing
    would cost five minutes and a 220 MB write to test an argument list.
    """
    calls = []

    def main(**kwargs):
        calls.append(kwargs)
        # celldega os.chdirs into the dataset directory. Reproduced because the
        # wrapper's cwd restore is a real promise about a process-global effect.
        os.chdir(kwargs["data_root_dir"])

    pre = SimpleNamespace(main=main)
    monkeypatch.setitem(sys.modules, "celldega", SimpleNamespace(pre=pre, __version__="0.24.2"))
    monkeypatch.setitem(sys.modules, "celldega.pre", pre)
    monkeypatch.setitem(sys.modules, "pyvips", SimpleNamespace(Image=object))
    return calls


def test_the_export_points_celldega_at_the_farm_and_not_at_the_data(dataset, stub_celldega):
    out = de.export_degafiles(dataset, tile_size=250, image_tile_layer="dapi",
                              max_workers=2)

    assert out == dataset / "degafiles" and out.is_dir()
    (kwargs,) = stub_celldega
    staging = de.staging_dir(dataset.resolve())
    assert kwargs["sample"] == staging.name == dataset.name
    assert Path(kwargs["data_root_dir"]) == staging.parent
    assert Path(kwargs["path_dega_files"]) == out
    assert (kwargs["tile_size"], kwargs["image_tile_layer"],
            kwargs["max_workers"]) == (250, "dapi", 2)
    # The sample directory it was handed is the farm, not the dataset.
    assert Path(kwargs["data_root_dir"]) != dataset.parent


def test_the_export_restores_the_working_directory(dataset, stub_celldega, tmp_path):
    """celldega's unzipper chdirs, and the wrapper cannot stop it.

    While it is chdir'd, every relative path *in this process* resolves
    somewhere else — a napari layer saved from another thread mid-export would
    land in the dataset directory. So the restore is asserted, not assumed.
    """
    here = tmp_path / "cwd"
    here.mkdir()
    before = Path.cwd()
    os.chdir(here)
    try:
        de.export_degafiles(dataset, tile_size=250, image_tile_layer="dapi",
                            max_workers=1)
        assert Path.cwd().resolve() == here.resolve()
    finally:
        os.chdir(before)


def test_the_export_extracts_before_calling_celldega(dataset, stub_celldega):
    """Their unzipper skips a target that already exists, so this makes it a no-op."""
    de.export_degafiles(dataset, tile_size=250, image_tile_layer="dapi",
                        max_workers=1)
    staging = de.staging_dir(dataset.resolve())
    assert (staging / "cells.csv").is_file()
    assert (staging / "cell_feature_matrix" / "contents.txt").exists()


def test_a_chosen_destination_is_honoured(dataset, stub_celldega, tmp_path):
    out = tmp_path / "elsewhere" / "degafiles"
    assert de.export_degafiles(dataset, out_dir=out) == out
    assert Path(stub_celldega[0]["path_dega_files"]) == out


# ── the missing-dependency path ──────────────────────────────────────────────

@pytest.mark.skipif(CELLDEGA_INSTALLED, reason="celldega is installed here")
def test_the_absence_of_celldega_is_reported_with_its_remedy():
    """The message is most of what ``dega_available`` is for.

    ``pip install celldega`` is the wrong advice: its caps are conservative and
    a plain resolve downgrades anndata, spatialdata and pandas underneath the
    viewer. Anyone who reads only the error has to be told ``--no-deps``.
    """
    ok, message = de.dega_available()
    assert not ok
    assert "--no-deps" in message
    assert de.CELLDEGA_PIN in message

    with pytest.raises(de.DegaUnavailable) as excinfo:
        de.require_dega()
    assert "--no-deps" in str(excinfo.value)


@pytest.mark.skipif(CELLDEGA_INSTALLED, reason="celldega is installed here")
def test_the_export_refuses_before_it_stages_anything(dataset):
    """Nothing is written for an export that cannot run."""
    with pytest.raises(de.DegaUnavailable):
        de.export_degafiles(dataset)
    assert not (dataset / "viewer_cache" / "dega_staging").exists()


def test_dega_available_imports_celldega_before_pyvips():
    """Source guard: the two imports in ``dega_available`` must not be sorted.

    Order-dependent collision, measured both ways on a box with Ubuntu's
    libvips. pyvips runs against the host ``libvips.so.42``, which drags in the
    system HDF5 1.10 while this env's h5py is built for 1.14:

    * pyvips then h5py -> ``ValueError: Not a datatype``; anndata is dead for
      the rest of the process.
    * h5py then pyvips -> fine; an h5ad round-trip afterwards succeeds.

    ``import celldega.pre`` pulls in anndata and therefore h5py, so putting it
    first is what keeps the *availability check* from breaking the session it
    was only supposed to ask a question about. Nothing about the code says so —
    alphabetising the imports would look like tidying and would reintroduce it
    on any host with libvips installed. Hence a test.
    """
    import ast
    import inspect

    source = inspect.getsource(de.dega_available)
    names = [alias.name for node in ast.walk(ast.parse(source.strip()))
             if isinstance(node, ast.Import) for alias in node.names]
    assert "celldega.pre" in names and "pyvips" in names
    assert names.index("celldega.pre") < names.index("pyvips"), (
        f"pyvips is imported before celldega in dega_available ({names}); that "
        f"order makes an availability check break h5py on a host with libvips"
    )


def test_the_install_hint_is_two_commands_and_does_not_use_the_pre_extra():
    """``--no-deps`` makes pip ignore extras, so ``celldega[pre]`` is a trap.

    It reads as though it fetches pyvips and does not: measured with
    ``pip install --dry-run --no-deps 'celldega[pre]==0.24.2'``, whose whole
    output is "Would install celldega-0.24.2". The remedy therefore names every
    dependency itself, on a second line without ``--no-deps``.
    """
    hint = de.install_hint()
    lines = [line.strip() for line in hint.splitlines() if line.strip()]
    assert len(lines) == 2, f"expected two commands, got {lines}"

    no_deps, rest = lines
    assert "--no-deps" in no_deps and de.CELLDEGA_PIN in no_deps
    assert "[pre]" not in hint, "the extra does nothing under --no-deps"
    assert "--no-deps" not in rest, (
        "the second command must resolve normally, or it puts nothing back"
    )
    for package in de.PYVIPS_PINS + de.CELLDEGA_RUNTIME_DEPS:
        assert package in rest, f"{package} is not in the install hint"


def test_the_remedy_a_user_sees_is_the_install_hint():
    """Whatever dega_available says, it must not disagree with install_hint."""
    ok, message = de.dega_available()
    if ok:
        pytest.skip("celldega is installed here, so there is no remedy to show")
    assert de.install_hint() in message


def test_the_pin_does_not_carry_the_pre_extra():
    """The ``[pre]`` spelling was removed because it did nothing.

    pyvips is declared only under celldega's ``pre`` extra, and their import of
    it is a try/except that leaves the module ``None`` -- so without pyvips the
    export runs for minutes, reaches "generating dapi image tiles" and dies with
    ``AttributeError: 'NoneType' object has no attribute 'Image'``. But pip
    ignores extras under ``--no-deps``, so naming the extra never prevented any
    of that. The dependency list does.
    """
    assert "[pre]" not in de.CELLDEGA_PIN, (
        "pip ignores extras under --no-deps, so [pre] promises pyvips and "
        "delivers nothing; the dependencies are named explicitly instead"
    )
    assert "==" in de.CELLDEGA_PIN, "a floating pin here was the insitucnv mistake"


# ── the template, against the wrapper it calls ───────────────────────────────

def _rendered_export_call(params: dict) -> ast.Call:
    """The ``export_degafiles(...)`` call in the rendered template."""
    from palms.utils.step_templates import builtin_assemble

    code = ast.parse(
        __import__("string").Template(
            builtin_assemble("export.degafiles", ["main"])).substitute(params))
    for node in ast.walk(code):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id == "export_degafiles"):
            return node
    raise AssertionError("the template does not call export_degafiles")


def test_the_template_calls_the_wrapper_and_not_celldega_directly():
    """A cell calling ``dega.pre.main`` would decompress into the reader's data.

    That is the whole reason this template is allowed to import ``palms``: the
    staging is not a convenience the notebook can skip.
    """
    from palms.utils.step_templates import builtin_assemble

    text = builtin_assemble("export.degafiles", ["main"])
    assert "from palms.utils.dega_export import export_degafiles" in text
    # Comments stripped: the prose above the code explains why it does not call
    # `celldega.pre.main`, and a substring search would read that as the defect.
    code = "\n".join(line for line in text.splitlines()
                     if not line.lstrip().startswith("#"))
    assert "celldega" not in code, (
        "the template reaches celldega's entry point past the wrapper"
    )


def test_every_argument_the_template_passes_is_one_the_wrapper_takes():
    """The drift this pair is exposed to: a template renders, runs and is wrong.

    The template is text and the wrapper is code, so a renamed keyword is not a
    syntax error anywhere — it is a ``TypeError`` at minute zero of a five
    minute export, or worse, a silently ignored setting if the signature ever
    grows ``**kwargs``.
    """
    import inspect

    from palms.utils.step_templates import builtin_spec

    spec = builtin_spec("export.degafiles")
    call = _rendered_export_call(spec.synth_params())

    accepted = inspect.signature(de.export_degafiles).parameters
    assert not any(p.kind is p.VAR_KEYWORD for p in accepted.values()), (
        "a **kwargs here would swallow a renamed argument silently"
    )
    for keyword in call.keywords:
        assert keyword.arg in accepted, (
            f"the template passes {keyword.arg}=, which export_degafiles does not take"
        )
    assert len(call.args) == 1, "the dataset is the one positional argument"


def test_the_publish_tab_renders_the_template_it_advertises(qapp, tmp_path):
    """The provider's params must actually render — and bind the declared output.

    The registry-wide gates check that a provider exists and answers; this is
    the one that runs its answer through the template, which is what the button
    does a moment later.
    """
    from palms.utils.step_templates import step_template
    from palms.utils.steps import Step, free_names
    from palms.utils.prov_graph import TERMINAL
    from palms.tabs.tab_publish import NODE_ID, TEMPLATE_ID, build_tab

    ctx = SimpleNamespace(state={}, viewer=None, data_path=tmp_path)
    # Held for the length of the test: the provider closes over magicgui
    # widgets, and once Qt has collected the tab every one of them raises
    # RuntimeError on read. The same trap ``live_providers`` documents.
    _held = build_tab(ctx)
    blocks, params, _note = ctx.state["template_preview"][TEMPLATE_ID]()

    assert params["out_dir"] == "degafiles", (
        "an absolute destination would pin the notebook to this machine"
    )
    step = Step(id=NODE_ID, **step_template(TEMPLATE_ID, blocks), params=params,
                deps=["preamble"], kind=TERMINAL, label="Publish DegaFiles",
                outputs=["degafiles_path"])
    code = step.render()
    assert "degafiles_path = export_degafiles(" in code
    # Only names the notebook preamble binds; `data_path` is the whole of it.
    assert free_names(code) <= {"data_path"}

def test_the_default_tiles_every_channel_celldega_names():
    """xv-11y: a published export must not silently drop three of four channels.

    ``celldega.pre.run_pre_processing.main`` defaults to ``'all'``, and for
    Xenium that is the four channels ``get_image_info('Xenium', 'all')`` lists —
    one per ``morphology_focus_000{0..3}.ome.tif`` the bundle ships. PALMS
    defaulted to ``'dapi'`` until 2026-09-05, which produced a Landscape with a
    single toggle and no indication that anything was missing.

    The narrow option stays available: this pins which one is the *default*, not
    which ones exist.
    """
    from palms.utils import dega_export as de

    assert de.IMAGE_TILE_LAYER == "all"


def test_the_tab_offers_both_layers_and_says_what_the_choice_costs():
    """The default is only half of it — the tradeoff has to be visible.

    Four pyramids is roughly four times the runtime and the output size, which
    is a real reason to pick ``dapi``; it is not a reason to make that choice
    silently.
    """
    import ast

    source = (Path(__file__).resolve().parent.parent / "src" / "palms"
              / "tabs" / "tab_publish.py").read_text()
    tree = ast.parse(source)
    call = next(
        n for n in ast.walk(tree)
        if isinstance(n, ast.Call)
        and any(getattr(kw, "arg", None) == "label"
                and getattr(kw.value, "value", None) == "Image layer"
                for kw in n.keywords)
    )
    kwargs = {kw.arg: kw.value for kw in call.keywords}
    choices = ast.literal_eval(kwargs["choices"])
    assert set(choices) == {"all", "dapi"}, "both must stay reachable"
    assert kwargs["value"].id == "IMAGE_TILE_LAYER", (
        "the widget must follow the module default, not restate it"
    )
    tooltip = ast.literal_eval(kwargs["tooltip"])
    assert "four" in tooltip.lower(), "the tooltip must say what 'all' costs"
