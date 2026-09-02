"""When the loader is allowed to throw a cache away.

The reported complaint was that it did so too easily. Three separate paths did:

* an unreadable cache was renamed aside — or ``rmtree``'d if the rename failed —
  with no repair attempt and no check for user data;
* staleness compared ``experiment.xenium``'s mtime against the cache *directory*
  mtime, so ``rsync``/``cp -p``/a re-download condemned a perfectly good cache;
* the sidecar list omitted the CNV caches, so a cache whose only user data was a
  multi-hour CopyKAT run reported "no user data" and was rebuilt with no dialog.

These are pure-function tests; none of them builds a 30 GB anything.

Run headless:  QT_QPA_PLATFORM=offscreen pytest tests/test_loader_policy.py
"""
from __future__ import annotations

import os
import shutil
import sys
import time
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

pytest.importorskip("spatialdata")
np = pytest.importorskip("numpy")

from palms import loader  # noqa: E402
from palms.loader import (  # noqa: E402
    _detect_user_data, _format_user_data_message, _has_any_user_data,
    _is_cache_stale, _restore_user_elements, write_manifest,
)


@pytest.fixture
def fake_cache(tmp_path):
    """A directory shaped like a cache, without the cost of building one.

    Includes a stand-in for the raw 10x output, because staleness and rebuild
    policy only mean anything for a dataset that *can* be rebuilt. The
    cache-only shape is ``crop_export_dir`` below.
    """
    cache = tmp_path / "sdata_cached.zarr"
    (cache / "tables" / "table" / "obs").mkdir(parents=True)
    (cache / "tables" / "table" / "uns").mkdir(parents=True)
    (cache / "tables" / "table" / "obsm").mkdir(parents=True)
    experiment = tmp_path / "experiment.xenium"
    experiment.write_text('{"pixel_size": 0.2125}')
    (tmp_path / "cells.zarr.zip").write_bytes(b"PK\x05\x06" + b"\0" * 18)
    return cache, experiment


@pytest.fixture
def crop_export_dir(tmp_path):
    """What Crop Dataset writes: no raw 10x output, so nothing to rebuild from.

    ``experiment.xenium`` *is* present — ``crop_export`` copies it and ``app.py``
    refuses to open a directory without one. What is absent is everything
    ``spatialdata_io.xenium`` actually reads.
    """
    cache = tmp_path / "sdata_cached.zarr"
    (cache / "tables" / "table" / "obs").mkdir(parents=True)
    experiment = tmp_path / "experiment.xenium"
    experiment.write_text('{"pixel_size": 0.2125}')
    (tmp_path / "transcripts.parquet").write_bytes(b"PAR1")
    (tmp_path / "transcript_cache").mkdir()
    return tmp_path


# ── the silent-rebuild case ──────────────────────────────────────────────────

def test_copykat_results_count_as_user_data(fake_cache):
    """The worst case in the report: hours of compute, rebuilt with no dialog."""
    cache, _ = fake_cache
    (cache / "adata_cnv_cache_copykat.h5ad").write_bytes(b"x")
    (cache / "cnv_copykat_result.json").write_text("{}")

    user_data = _detect_user_data(cache)
    assert _has_any_user_data(user_data)
    assert "adata_cnv_cache_copykat.h5ad" in user_data["sidecars"]

    message = _format_user_data_message(user_data)
    assert "CopyKAT" in message and "hours of compute" in message


def test_infercnv_results_count_as_user_data(fake_cache):
    cache, _ = fake_cache
    (cache / "adata_cnv_cache_infercnv.h5ad").write_bytes(b"x")
    assert _has_any_user_data(_detect_user_data(cache))
    assert "inferCNV" in _format_user_data_message(_detect_user_data(cache))


def test_an_untouched_cache_has_no_user_data(fake_cache):
    cache, _ = fake_cache
    assert not _has_any_user_data(_detect_user_data(cache))


def test_clusterings_and_rois_still_count(fake_cache):
    cache, _ = fake_cache
    (cache / "tables" / "table" / "obs" / "clustering_leiden_r1.0").mkdir()
    (cache / "shapes" / "rois").mkdir(parents=True)

    user_data = _detect_user_data(cache)
    assert _has_any_user_data(user_data)
    assert user_data["clusterings"] == ["clustering_leiden_r1.0"]
    assert "ROIs" in _format_user_data_message(user_data)


def test_external_images_and_their_landmarks_count_as_user_data(fake_cache):
    """These are named per file, so the old fixed key list could not see them.

    A dataset with a registered PhenoCycler image reported "no user data" and
    could be rebuilt over without a prompt.
    """
    cache, _ = fake_cache
    (cache / "images" / "ext_PostXenium5k_region1_ome").mkdir(parents=True)
    (cache / "shapes" / "ext_PostXenium5k_region1_ome_xenium_lm").mkdir(parents=True)

    user_data = _detect_user_data(cache)
    assert _has_any_user_data(user_data)
    assert "ext_PostXenium5k_region1_ome" in user_data["images"]
    assert "ext_PostXenium5k_region1_ome_xenium_lm" in user_data["shapes"]


def test_patch_overlays_count_as_user_data(fake_cache):
    cache, _ = fake_cache
    (cache / "shapes" / "patch_tumour_regions").mkdir(parents=True)
    assert "patch_tumour_regions" in _detect_user_data(cache)["shapes"]


def test_pipeline_elements_are_not_mistaken_for_user_data(fake_cache):
    """cell_circles and morphology_focus are rebuilt from raw, not user work."""
    cache, _ = fake_cache
    (cache / "shapes" / "cell_circles").mkdir(parents=True)
    (cache / "images" / "morphology_focus").mkdir(parents=True)
    assert not _has_any_user_data(_detect_user_data(cache))


def test_one_label_per_cnv_backend_not_per_file(fake_cache):
    """The h5ad and its result JSON both mean "CopyKAT ran"."""
    cache, _ = fake_cache
    (cache / "adata_cnv_cache_copykat.h5ad").write_bytes(b"x")
    (cache / "cnv_copykat_result.json").write_text("{}")

    message = _format_user_data_message(_detect_user_data(cache))
    assert message.count("CopyKAT CNV results") == 1


# ── staleness ────────────────────────────────────────────────────────────────

def test_without_a_manifest_staleness_is_uncertain(fake_cache):
    cache, experiment = fake_cache
    os.utime(experiment, (time.time() + 100, time.time() + 100))

    stale, certain = _is_cache_stale(cache, experiment)
    assert stale and not certain      # uncertain ⇒ the caller must ask


def test_a_manifest_makes_a_touched_source_a_non_event(fake_cache):
    """Copying a dataset bumps the mtime without changing a byte."""
    cache, experiment = fake_cache
    write_manifest(cache, experiment)
    os.utime(experiment, (time.time() + 100, time.time() + 100))

    stale, certain = _is_cache_stale(cache, experiment)
    assert not stale and certain


def test_a_manifest_still_detects_a_real_change(fake_cache):
    cache, experiment = fake_cache
    write_manifest(cache, experiment)
    experiment.write_text('{"pixel_size": 0.4250}')

    stale, certain = _is_cache_stale(cache, experiment)
    assert stale and certain


def test_a_missing_source_is_never_stale(fake_cache):
    """No experiment.xenium ⇒ nothing to compare against.

    This used to say "a Crop Dataset export has no experiment.xenium", which is
    false — ``crop_export`` copies it (``app.py`` refuses to open a directory
    without one). That wrong mental model is part of why cache-only datasets
    reached the rebuild branches at all; ``has_raw_xenium_source`` is what
    actually identifies them, and it is exercised below.
    """
    cache, experiment = fake_cache
    experiment.unlink()
    assert _is_cache_stale(cache, experiment) == (False, True)


def test_manifest_records_the_versions_it_was_built_with(fake_cache):
    from palms.utils.cache_repair import read_manifest

    cache, experiment = fake_cache
    write_manifest(cache, experiment)
    manifest = read_manifest(cache)
    assert manifest["source"] == "experiment.xenium"
    assert len(manifest["source_sha256"]) == 64
    assert "spatialdata_version" in manifest


def test_the_manifest_is_invisible_to_readers(tiny_sdata):
    """It lives in the store root, so it must not look like an element."""
    import warnings

    from spatialdata import read_zarr

    cache = Path(tiny_sdata.path)
    experiment = cache.parent / "experiment.xenium"
    experiment.write_text("{}")
    write_manifest(cache, experiment)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert sorted(read_zarr(str(cache)).tables) == ["table"]


# ── restoring user data on rebuild ───────────────────────────────────────────

def _adata_pair():
    import anndata
    import pandas as pd

    old = anndata.AnnData(np.ones((4, 2), dtype="float32"))
    old.obs_names = [f"c{i}" for i in range(4)]
    old.obs["clustering_leiden_r1.0"] = pd.Categorical(["0", "0", "1", "1"])
    old.obs["cnv_score_copykat"] = [0.1, 0.2, 0.3, 0.4]
    old.obs["copykat_leiden_res0.2"] = pd.Categorical(["a", "a", "b", "b"])
    old.obs["total_counts"] = [1, 2, 3, 4]          # rebuilt anyway
    old.uns["rank_genes_groups"] = {"names": ["g1"]}
    old.uns["cnv_runs"] = {"copykat": {"resolution": 0.2}}
    old.uns["rank_genes_groupby"] = "clustering_leiden_r1.0"
    old.obsm["X_umap"] = np.zeros((4, 2))
    old.obsm["X_cnv_umap"] = np.ones((4, 2))

    new = anndata.AnnData(np.ones((4, 2), dtype="float32"))
    new.obs_names = [f"c{i}" for i in range(4)]
    new.obs["total_counts"] = [9, 9, 9, 9]
    return old, new


class _FakeSdata:
    def __init__(self, table):
        self.tables = {"table": table}
        self.shapes: dict = {}
        self.images: dict = {}

    def __getitem__(self, key):
        return self.tables[key]

    def __setitem__(self, key, value):
        self.tables[key] = value


def test_restore_covers_cnv_and_rank_genes_keys():
    """Regression: these were dropped even on 'Rebuild and restore my data'."""
    old, new = _adata_pair()
    _restore_user_elements(_FakeSdata(old), _FakeSdata(new),
                           {"shapes": [], "images": [], "uns_keys": [],
                            "has_obsm_umap": True})

    assert "clustering_leiden_r1.0" in new.obs
    assert "cnv_score_copykat" in new.obs
    assert "copykat_leiden_res0.2" in new.obs
    assert new.uns["cnv_runs"] == {"copykat": {"resolution": 0.2}}
    assert new.uns["rank_genes_groupby"] == "clustering_leiden_r1.0"
    assert "X_umap" in new.obsm and "X_cnv_umap" in new.obsm


def test_restore_does_not_clobber_freshly_built_columns():
    """A rebuilt total_counts must win over the stale one."""
    old, new = _adata_pair()
    _restore_user_elements(_FakeSdata(old), _FakeSdata(new),
                           {"shapes": [], "images": [], "uns_keys": [],
                            "has_obsm_umap": True})
    assert list(new.obs["total_counts"]) == [9, 9, 9, 9]


def test_restore_reports_what_it_moved():
    old, new = _adata_pair()
    restored = _restore_user_elements(_FakeSdata(old), _FakeSdata(new),
                                      {"shapes": [], "images": [], "uns_keys": [],
                                       "has_obsm_umap": True})
    assert "clustering_leiden_r1.0" in restored
    assert "UMAP coordinates" in restored


# ── never discard without asking ─────────────────────────────────────────────

def test_stale_with_user_data_and_no_dialog_keeps_the_cache(monkeypatch, fake_cache):
    """The old default here was 'rebuild' — a silent 30 GB discard."""
    cache, _ = fake_cache
    (cache / "adata_cnv_cache_copykat.h5ad").write_bytes(b"x")
    monkeypatch.setattr(loader, "_qt_message_box", lambda: None)

    assert loader._ask_rebuild_preference(_detect_user_data(cache)) == "keep"


def test_an_unopenable_cache_with_no_dialog_raises_rather_than_rebuilding(
    monkeypatch, fake_cache,
):
    from palms.utils import cache_repair

    cache, _ = fake_cache
    (cache / "adata_cnv_cache_copykat.h5ad").write_bytes(b"x")
    monkeypatch.setattr(loader, "_qt_message_box", lambda: None)

    with pytest.raises(loader.CacheLoadAborted, match="Refusing to rebuild"):
        loader._ask_corrupt_cache(
            OSError("bad store"), cache_repair.verify(cache), _detect_user_data(cache),
        )


# ── the headless opt-in: palms-build-cache --on-stale ───────────────────────
#
# The rule above ("never discard without asking") left a terminal with no way to
# say yes. --on-stale is that yes, and these check it stays an opt-in.

def test_without_on_stale_a_stale_cache_with_user_data_is_still_kept(
    monkeypatch, fake_cache,
):
    """The default is unchanged: no answer given, no dialog available, keep."""
    cache, _ = fake_cache
    (cache / "adata_cnv_cache_copykat.h5ad").write_bytes(b"x")
    monkeypatch.setattr(loader, "_qt_message_box", lambda: None)

    assert loader._stale_preference(_detect_user_data(cache), certain=True) == "keep"


@pytest.mark.parametrize("answer", ["rebuild", "restore", "keep"])
def test_on_stale_answers_without_a_dialog(monkeypatch, fake_cache, answer):
    cache, _ = fake_cache
    (cache / "adata_cnv_cache_copykat.h5ad").write_bytes(b"x")
    monkeypatch.setattr(loader, "_qt_message_box", lambda: None)

    assert loader._stale_preference(
        _detect_user_data(cache), certain=True, on_stale=answer) == answer


def test_on_stale_also_covers_the_branches_that_never_prompt(fake_cache):
    """An uncertain, empty cache rebuilds silently — 'keep' must still stop it.

    This is the branch that returns before any dialog is consulted, so a preset
    threaded only through ``_ask_rebuild_preference`` would have missed it.
    """
    cache, _ = fake_cache
    empty = _detect_user_data(cache)

    assert loader._stale_preference(empty, certain=False) is None
    assert loader._stale_preference(empty, certain=False, on_stale="keep") == "keep"


def test_on_stale_keep_still_refuses_an_unopenable_cache(monkeypatch, fake_cache):
    """'keep' is not an answer here — the cache does not open, so keeping it is
    not a way to carry on."""
    from palms.utils import cache_repair

    cache, _ = fake_cache
    (cache / "adata_cnv_cache_copykat.h5ad").write_bytes(b"x")
    monkeypatch.setattr(loader, "_qt_message_box", lambda: None)

    with pytest.raises(loader.CacheLoadAborted):
        loader._ask_corrupt_cache(
            OSError("bad store"), cache_repair.verify(cache),
            _detect_user_data(cache), preset="keep",
        )

    for answer in ("rebuild", "restore"):
        assert loader._ask_corrupt_cache(
            OSError("bad store"), cache_repair.verify(cache),
            _detect_user_data(cache), preset=answer,
        ) == answer


def test_a_mistyped_on_stale_is_an_error_not_a_shrug(tmp_path):
    """It must not degrade to 'ask', which on a headless box means 'keep'."""
    with pytest.raises(ValueError, match="on_stale"):
        loader.load_sdata(tmp_path, on_stale="rebiuld")


# ── the console script ───────────────────────────────────────────────────────

def test_build_cache_parser_maps_flags_to_load_sdata_kwargs():
    parser = loader._build_parser()

    args = parser.parse_args(["/data"])
    assert args.on_stale == "ask"          # translated to on_stale=None in main()
    assert args.n_jobs == 8
    assert not args.no_pyramid and not args.no_cache and not args.check

    args = parser.parse_args(
        ["/data", "--on-stale", "restore", "--no-pyramid", "--n-jobs", "2", "--check"])
    assert (args.on_stale, args.n_jobs, args.no_pyramid, args.check) == (
        "restore", 2, True, True)


def test_build_cache_parser_rejects_an_unknown_on_stale():
    with pytest.raises(SystemExit):
        loader._build_parser().parse_args(["/data", "--on-stale", "maybe"])


def test_the_console_script_target_exists():
    """pyproject points palms-build-cache at this name."""
    entry = "palms-build-cache = \"palms.loader:main\""
    pyproject = (Path(__file__).resolve().parent.parent / "pyproject.toml").read_text()
    assert entry in pyproject
    assert callable(loader.main)


def test_check_reports_a_missing_cache_without_creating_one(tmp_path, capsys):
    assert loader._check_cache(tmp_path) == 1
    assert "No zarr cache" in capsys.readouterr().out
    assert not (tmp_path / "sdata_cached.zarr").exists()


def test_check_is_read_only_on_a_cache_it_condemns(fake_cache, capsys):
    """--check must not touch the store, least of all a broken one."""
    cache, experiment = fake_cache
    (cache / "adata_cnv_cache_copykat.h5ad").write_bytes(b"x")
    before = sorted(p.relative_to(cache) for p in cache.rglob("*"))

    # No manifest and a source newer than the cache ⇒ reported as possibly stale.
    os.utime(experiment, (time.time() + 10, time.time() + 10))
    assert loader._check_cache(cache.parent) == 1

    out = capsys.readouterr().out
    assert "stale" in out.lower()
    assert "adata_cnv_cache_copykat.h5ad" in out
    assert sorted(p.relative_to(cache) for p in cache.rglob("*")) == before


# ── a dataset with no raw source is never told to rebuild (issue #17) ────────
#
# A Crop Dataset export ships experiment.xenium + sdata_cached.zarr + derived
# transcripts, and nothing else. Without a manifest, freshness fell back to
# comparing experiment.xenium's mtime against the cache directory's — which a
# copy, an unzip or a sync reorders — and three separate branches then tried to
# rebuild from raw files that were never there. Two of them renamed the only
# copy of the data aside *first*.

def test_a_crop_export_is_recognised_as_having_no_raw_source(crop_export_dir):
    assert not loader.has_raw_xenium_source(crop_export_dir)


def test_a_normal_dataset_is_recognised_as_rebuildable(fake_cache):
    assert loader.has_raw_xenium_source(fake_cache[0].parent)


@pytest.mark.parametrize(
    "marker", ["cells.zarr.zip", "cell_feature_matrix.h5", "morphology_focus"])
def test_any_raw_marker_is_enough_to_count_as_rebuildable(crop_export_dir, marker):
    """Conservative on purpose: partial raw output is broken raw output, and
    should keep raising whatever error says so — not be reclassified."""
    (crop_export_dir / marker).write_bytes(b"x")
    assert loader.has_raw_xenium_source(crop_export_dir)


# ── the declaration beats the inference (xv-iy9) ─────────────────────────────
#
# crop_export stamps ``cache_only: True`` into the export's cache manifest, and
# the manifest writer's docstring says it is there "for readers that would
# otherwise have to infer it from absent files". Nothing read it: both loader
# sites, the dataset inventory, the rebuild guard and the recorded preamble all
# recomputed ``not has_raw_xenium_source(path)`` — the very inference the stamp
# exists to replace. Harmless only while an export writes no raw-shaped file.

def _declare_cache_only(path: Path, value) -> None:
    import json
    from palms.utils.zarr_safe import MANIFEST_FILE
    cache = path / "sdata_cached.zarr"
    cache.mkdir(parents=True, exist_ok=True)
    (cache / MANIFEST_FILE).write_text(json.dumps({"cache_only": value}))


def test_a_declared_cache_only_dataset_stays_cache_only_beside_raw_files(crop_export_dir):
    """The regression the stamp exists to prevent, one release early."""
    _declare_cache_only(crop_export_dir, True)
    (crop_export_dir / "cell_feature_matrix.h5").write_bytes(b"x")

    assert loader.has_raw_xenium_source(crop_export_dir)      # the inference flips
    assert loader.is_cache_only(crop_export_dir)              # the declaration does not


def test_without_a_stamp_the_inference_still_decides(crop_export_dir):
    assert loader.cache_only_declared(crop_export_dir) is None
    assert loader.is_cache_only(crop_export_dir)
    (crop_export_dir / "cells.zarr.zip").write_bytes(b"PK")
    assert not loader.is_cache_only(crop_export_dir)


def test_a_false_stamp_cannot_unset_cache_only(crop_export_dir):
    """The stamp may add certainty, never remove protection.

    Trusting ``false`` here would send the loader down a rebuild path on a
    dataset whose only copy of the data is the cache.
    """
    _declare_cache_only(crop_export_dir, False)
    assert loader.cache_only_declared(crop_export_dir) is False
    assert loader.is_cache_only(crop_export_dir)


@pytest.mark.parametrize("payload", ["{not json", '{"cache_only": "yes"}', "{}"])
def test_an_unreadable_or_silent_manifest_says_nothing(crop_export_dir, payload):
    """A manifest that cannot be parsed is a manifest that makes no claim —
    never an exception on the load path."""
    cache = crop_export_dir / "sdata_cached.zarr"
    cache.mkdir(parents=True, exist_ok=True)
    from palms.utils.zarr_safe import MANIFEST_FILE
    (cache / MANIFEST_FILE).write_text(payload)

    assert loader.cache_only_declared(crop_export_dir) is None
    assert loader.is_cache_only(crop_export_dir)


def test_the_rebuild_guard_describes_a_declared_export_without_claiming_files_are_absent(
    crop_export_dir,
):
    """The message enumerates missing files; at a declared export they are there."""
    from palms.tabs import tab_cache

    _declare_cache_only(crop_export_dir, True)
    (crop_export_dir / "cell_feature_matrix.h5").write_bytes(b"x")

    reason = tab_cache._rebuild_blocked_reason(
        type("Ctx", (), {"data_path": crop_export_dir})()
    )
    assert reason and "Recover from Backup" in reason
    assert "cache manifest declares it" in reason
    assert "no cell_feature_matrix" not in reason


def test_a_cache_only_dataset_loads_from_cache_and_moves_nothing(
    monkeypatch, crop_export_dir,
):
    """The reported trigger: experiment.xenium newer than the cache, no manifest.

    Asserted behaviourally with a ``shutil.move`` spy rather than by grepping the
    source, because what matters is that no rename *happens*, not that no rename
    is written down.
    """
    os.utime(crop_export_dir / "experiment.xenium",
             (time.time() + 100, time.time() + 100))

    moves = []
    monkeypatch.setattr(loader.shutil, "move",
                        lambda src, dst: moves.append((src, dst)))
    monkeypatch.setattr(loader, "_open_cache", lambda p: "SDATA")
    # Would prompt, or return a destructive default, if it were ever consulted.
    monkeypatch.setattr(loader, "_stale_preference",
                        lambda *a, **k: pytest.fail("staleness was consulted"))

    assert loader.load_sdata(crop_export_dir) == "SDATA"
    assert moves == []


def test_a_cache_only_dataset_is_stamped_so_the_next_run_need_not_re_derive(
    monkeypatch, crop_export_dir,
):
    from palms.utils.cache_repair import read_manifest

    monkeypatch.setattr(loader, "_open_cache", lambda p: "SDATA")
    loader.load_sdata(crop_export_dir)

    manifest = read_manifest(crop_export_dir / "sdata_cached.zarr")
    assert manifest["cache_only"] is True
    assert len(manifest["source_sha256"]) == 64


@pytest.mark.parametrize("answer", ["rebuild", "restore"])
def test_on_stale_cannot_authorise_an_impossible_rebuild(
    monkeypatch, crop_export_dir, answer,
):
    """--on-stale is a yes to a question that is not being asked here."""
    os.utime(crop_export_dir / "experiment.xenium",
             (time.time() + 100, time.time() + 100))
    moves = []
    monkeypatch.setattr(loader.shutil, "move",
                        lambda src, dst: moves.append((src, dst)))
    monkeypatch.setattr(loader, "_open_cache", lambda p: "SDATA")

    assert loader.load_sdata(crop_export_dir, on_stale=answer) == "SDATA"
    assert moves == []


def test_an_unopenable_cache_only_cache_refuses_instead_of_renaming(
    monkeypatch, crop_export_dir,
):
    """Every non-quit answer _ask_corrupt_cache can give moves the cache aside
    and then rebuilds — here that renames away the only copy and fails anyway."""
    def _boom(path):
        raise OSError("bad store")

    moves = []
    monkeypatch.setattr(loader.shutil, "move",
                        lambda src, dst: moves.append((src, dst)))
    monkeypatch.setattr(loader, "_open_cache", _boom)
    monkeypatch.setattr(loader, "_ask_corrupt_cache",
                        lambda *a, **k: pytest.fail("would have offered a rebuild"))

    with pytest.raises(loader.NoRawSourceError) as excinfo:
        loader.load_sdata(crop_export_dir)

    message = str(excinfo.value)
    assert "nothing to rebuild from" in message
    assert "Nothing has been moved or deleted" in message
    assert "backup" in message.lower()
    assert moves == []
    assert (crop_export_dir / "sdata_cached.zarr").exists()


def test_no_cache_flag_still_loads_a_cache_only_dataset(monkeypatch, crop_export_dir):
    """--no-cache means 'read the raw files', which do not exist here."""
    monkeypatch.setattr(loader, "_open_cache", lambda p: "SDATA")
    assert loader.load_sdata(crop_export_dir, use_cache=False) == "SDATA"


def test_a_cache_only_dataset_with_no_cache_left_fails_legibly(crop_export_dir):
    """Otherwise this surfaces as h5py failing to open a file the user never had."""
    shutil.rmtree(crop_export_dir / "sdata_cached.zarr")

    with pytest.raises(loader.NoRawSourceError, match="no raw Xenium output"):
        loader.load_sdata(crop_export_dir)


def test_check_reports_a_cache_only_dataset_as_not_applicable(
    crop_export_dir, capsys,
):
    """Reporting it stale would be advising a rebuild that cannot happen."""
    os.utime(crop_export_dir / "experiment.xenium",
             (time.time() + 100, time.time() + 100))
    loader._check_cache(crop_export_dir)

    out = capsys.readouterr().out
    assert "Freshness: n/a" in out
    assert "STALE" not in out


def test_the_loader_never_deletes_a_cache_directory():
    """Every destructive branch must be a rename. Source-level guard."""
    import re

    source = (Path(__file__).resolve().parent.parent
              / "src" / "palms" / "loader.py").read_text()
    offenders = [
        line.strip() for line in source.splitlines()
        if re.search(r"rmtree\(\s*(str\()?cache_path", line)
    ]
    assert offenders == [], (
        "loader.py must never rmtree the live cache — move it aside instead: "
        + "; ".join(offenders)
    )


def test_a_repr_that_raises_does_not_condemn_the_cache(monkeypatch, crop_export_dir):
    """An unrenderable summary is not an unopenable store.

    ``load_sdata`` printed the SpatialData summary inside the ``try`` whose
    ``except`` means "the cache could not be opened", so anything raising in
    ``__repr__`` routed a perfectly good cache to ``_ask_corrupt_cache`` and
    offered to rebuild it — while ``cache_repair.verify`` printed
    "✓ Cache is healthy" in the next breath. It was found when a non-default
    parquet reader made spatialdata's repr raise — it walks each points
    element's dask graph for its backing files — on every launch after the one
    that built the cache. That reader is gone, but the policy is the point: the
    summary is diagnostics, and diagnostics must never condemn a store.
    """
    class _UnrenderableSData:
        def __repr__(self):
            raise TypeError("argument of type 'Task' is not iterable")

    sdata = _UnrenderableSData()
    monkeypatch.setattr(loader, "_open_cache", lambda p: sdata)
    monkeypatch.setattr(loader, "_ask_corrupt_cache",
                        lambda *a, **k: pytest.fail("a good cache was condemned"))
    monkeypatch.setattr(loader.shutil, "move",
                        lambda src, dst: pytest.fail("the cache was moved aside"))

    assert loader.load_sdata(crop_export_dir) is sdata



def test_build_cache_summary_survives_a_dataset_with_no_analysis_folder(
    monkeypatch, tmp_path, capsys,
):
    """The summary print must not turn a missing UMAP into a traceback.

    ``load_umap`` returns ``None`` when there is no ``analysis/`` folder — a 10x
    output bundle whose ``analysis.tar.gz`` was never extracted, or a Crop
    Dataset export. The cache is written *before* the summary runs, so a
    dereference here reports a completed build as a failure. Found on the 10x
    public ``Xenium_V1_human_Pancreas_FFPE`` bundle, which ships secondary
    analysis only as a tarball.
    """
    class _FakeSData(dict):
        images = {"morphology_focus": None}
        labels = {"cell_labels": None, "nucleus_labels": None}
        points = {"transcripts": None}
        shapes: dict = {}

        @property
        def tables(self):
            return {"table": self["table"]}

    class _FakeAData:
        shape = (140702, 377)

    sdata = _FakeSData(table=_FakeAData())
    monkeypatch.setattr(loader, "load_sdata", lambda *a, **k: sdata)
    monkeypatch.setattr(loader, "load_umap", lambda p: None)
    monkeypatch.setattr(loader, "load_clusterings", lambda p: {})
    monkeypatch.setattr(sys, "argv", ["palms-build-cache", str(tmp_path)])

    loader.main()

    out = capsys.readouterr().out
    assert "140702 cells x 377 genes" in out
    assert "UMAP:" in out and "analysis/" in out
    assert "Clusterings: []" in out


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
