"""Inspect and repair the SpatialData zarr cache.

The loader used to treat "``read_zarr`` raised" as "the cache is worthless": it
renamed the store aside, or deleted it outright if the rename failed. But the
most common failure is not data loss at all — it is a *consolidated metadata*
entry that no longer matches what is on disk, which is fixable in a second by
re-consolidating. Discarding instead costs a full rebuild from the raw Xenium
output (30 GB on the dataset that prompted this).

:func:`verify` is read-only and deliberately does not use ``zarr.open`` — it
parses the root ``zarr.json`` with :mod:`json`, so it still works on a store too
broken to open. :func:`repair` only ever renames, unlinks debris, or
re-consolidates; it never deletes an element that has no replacement.
"""

from __future__ import annotations

import json
import logging
import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from palms.utils.zarr_safe import (
    ELEMENT_TYPES, INTERNAL_NAMES, JOURNAL_DIR, MANIFEST_FILE, STAGING_DIR,
    TRASH_DIR, consolidate, list_trash, recover_pending,
)

log = logging.getLogger(__name__)

# Sidecar files the viewer writes into the store root. They are not zarr nodes,
# so they must never be mistaken for elements — and several represent hours of
# compute, so they must never be mistaken for debris either.
SIDECAR_PATTERNS = (
    "adata_norm_cache.h5ad",
    "adata_cnv_cache_*.h5ad",
    "roi_deg_cache.parquet",
    "arms_tile_deg_cache.parquet",
    "cnv_*_result.json",
)

BACKUP_PATTERNS = (
    "sdata_cached_corrupt_*.zarr",
    "sdata_cached_prev_*.zarr",
    "sdata_cached_backup_*.zarr",
)

# Repair levels, in increasing order of intervention.
AUTO = "auto"      # journals, debris, strays, re-consolidate — all reversible
FULL = "full"      # + restore a missing element from its .xv_trash backup


@dataclass
class HealthReport:
    """What ``verify`` found. Nothing here mutates the store."""

    cache_path: Path
    exists: bool = True
    readable_metadata: bool = True
    metadata_error: Optional[str] = None
    pending_ops: list[str] = field(default_factory=list)
    on_disk: list[str] = field(default_factory=list)
    missing_on_disk: list[str] = field(default_factory=list)
    missing_in_meta: list[str] = field(default_factory=list)
    stray_elements: list[str] = field(default_factory=list)
    debris: list[str] = field(default_factory=list)
    sidecars: list[str] = field(default_factory=list)
    backups: list[Path] = field(default_factory=list)
    trash_available: dict[str, list[Path]] = field(default_factory=dict)
    manifest: Optional[dict] = None

    @property
    def ok(self) -> bool:
        """True when nothing needs doing. Sidecars and backups are not faults."""
        return (
            self.exists
            and self.readable_metadata
            and not self.pending_ops
            and not self.missing_on_disk
            and not self.missing_in_meta
            and not self.stray_elements
            and not self.debris
        )

    @property
    def repairable(self) -> bool:
        """True when every fault found has a fix that loses nothing."""
        if not self.exists or not self.readable_metadata:
            return False
        recoverable = [e for e in self.missing_on_disk if e in self.trash_available]
        return len(recoverable) == len(self.missing_on_disk)

    def summary(self) -> str:
        if not self.exists:
            return "No cache at this path."
        lines: list[str] = []
        if not self.readable_metadata:
            lines.append(f"✗ Root metadata unreadable: {self.metadata_error}")
        if self.pending_ops:
            lines.append(f"⚠ {len(self.pending_ops)} interrupted write(s) to finish:")
            lines += [f"    {op}" for op in self.pending_ops]
        if self.missing_on_disk:
            lines.append("✗ Listed in metadata but missing from disk:")
            for element in self.missing_on_disk:
                has_backup = " (backup available)" if element in self.trash_available else ""
                lines.append(f"    {element}{has_backup}")
        if self.missing_in_meta:
            lines.append("⚠ Present on disk but missing from metadata "
                         "(fixed by re-consolidating):")
            lines += [f"    {element}" for element in self.missing_in_meta]
        if self.stray_elements:
            lines.append("⚠ Not a valid element, safe to drop:")
            lines += [f"    {element}" for element in self.stray_elements]
        if self.debris:
            lines.append(f"⚠ {len(self.debris)} leftover file(s) from an interrupted run")
        if not lines:
            lines.append("✓ Cache is healthy.")
        if self.on_disk:
            lines.append("")
            lines.append(f"Elements ({len(self.on_disk)}): " + ", ".join(self.on_disk))
        if self.sidecars:
            lines.append("Sidecar files: " + ", ".join(self.sidecars))
        if self.trash_available:
            lines.append("Recoverable backups: " + ", ".join(sorted(self.trash_available)))
        if self.backups:
            lines.append("Previous caches kept aside: "
                         + ", ".join(p.name for p in self.backups))
        return "\n".join(lines)


@dataclass
class RepairResult:
    actions: list[str] = field(default_factory=list)
    failures: list[str] = field(default_factory=list)

    @property
    def changed(self) -> bool:
        return bool(self.actions)


# ── reading the metadata without zarr ────────────────────────────────────────

def _consolidated_elements(cache_path: Path) -> tuple[set[str], Optional[str]]:
    """Element paths named by the root consolidated metadata.

    zarr 3 writes a *flat* map whose keys are full paths ("tables/table"), but
    older/nested layouts exist too, so both are handled. Parsed with ``json``
    rather than ``zarr.open`` so a store that cannot be opened still reports.
    """
    root = cache_path / "zarr.json"
    if not root.exists():
        return set(), "no zarr.json at the store root"
    try:
        document = json.loads(root.read_text())
    except (OSError, ValueError) as exc:
        return set(), str(exc)

    metadata = document.get("consolidated_metadata", {}).get("metadata", {})
    found: set[str] = set()
    for key, value in metadata.items():
        if "/" in key:
            etype, _, name = key.partition("/")
            if etype in ELEMENT_TYPES and "/" not in name:
                found.add(key)
            continue
        if key in ELEMENT_TYPES:
            nested = value.get("consolidated_metadata", {}).get("metadata", {})
            found.update(f"{key}/{name}" for name in nested)
    return found, None


def _disk_elements(cache_path: Path) -> list[str]:
    found: list[str] = []
    for etype in ELEMENT_TYPES:
        type_dir = cache_path / etype
        if not type_dir.is_dir():
            continue
        for entry in sorted(type_dir.iterdir()):
            if entry.is_dir() and not entry.name.startswith("."):
                found.append(f"{etype}/{entry.name}")
    return found


def _is_stray(cache_path: Path, element: str) -> bool:
    """A directory under a type dir that is not a real SpatialData element.

    Generalises the hard-coded ``tables/adata_norm`` cleanup that used to live
    in ``app.py``: any group without spatialdata's encoding attributes was
    written by something that is no longer part of the model.
    """
    marker = cache_path / element / "zarr.json"
    if not marker.exists():
        return True
    try:
        attributes = json.loads(marker.read_text()).get("attributes", {})
    except (OSError, ValueError):
        return True
    return not any(
        key in attributes
        for key in ("spatialdata-encoding-type", "spatialdata_attrs", "ome", "encoding-type")
    )


def _find_debris(cache_path: Path) -> list[str]:
    debris: list[str] = []
    staging = cache_path / STAGING_DIR
    if staging.is_dir():
        debris += [f"{STAGING_DIR}/{p.name}" for p in staging.iterdir()]
    journal = cache_path / JOURNAL_DIR
    if journal.is_dir():
        debris += [f"{JOURNAL_DIR}/{p.name}" for p in journal.glob("*.json")]
    debris += [
        str(p.relative_to(cache_path)) for p in cache_path.rglob("*.partial")
    ]
    return sorted(debris)


def _find_sidecars(cache_path: Path) -> list[str]:
    found: list[str] = []
    for pattern in SIDECAR_PATTERNS:
        found += [p.name for p in cache_path.glob(pattern) if p.is_file()]
    return sorted(set(found))


def find_backups(cache_path: Path) -> list[Path]:
    """Previous caches this or an earlier version moved aside."""
    parent = cache_path.parent
    found: list[Path] = []
    for pattern in BACKUP_PATTERNS:
        found += [p for p in parent.glob(pattern) if p.is_dir()]
    return sorted(found, reverse=True)


def read_manifest(cache_path: Path) -> Optional[dict]:
    path = cache_path / MANIFEST_FILE
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return None


# ── the public API ───────────────────────────────────────────────────────────

def verify(cache_path: Path) -> HealthReport:
    """Inspect *cache_path* without modifying anything."""
    cache_path = Path(cache_path)
    report = HealthReport(cache_path=cache_path)
    if not cache_path.exists():
        report.exists = False
        return report

    report.manifest = read_manifest(cache_path)
    report.backups = find_backups(cache_path)
    report.trash_available = list_trash(cache_path)
    report.sidecars = _find_sidecars(cache_path)
    report.debris = _find_debris(cache_path)

    in_metadata, error = _consolidated_elements(cache_path)
    if error is not None:
        report.readable_metadata = False
        report.metadata_error = error

    on_disk = _disk_elements(cache_path)
    report.stray_elements = [e for e in on_disk if _is_stray(cache_path, e)]
    valid_on_disk = [e for e in on_disk if e not in report.stray_elements]
    report.on_disk = valid_on_disk

    if report.readable_metadata:
        report.missing_on_disk = sorted(in_metadata - set(valid_on_disk))
        report.missing_in_meta = sorted(set(valid_on_disk) - in_metadata)

    journal_dir = cache_path / JOURNAL_DIR
    if journal_dir.is_dir():
        for path in sorted(journal_dir.glob("*.json")):
            try:
                record = json.loads(path.read_text())
                report.pending_ops.append(
                    f"{record.get('op', 'write')} {record.get('name', path.stem)}")
            except (OSError, ValueError):
                report.pending_ops.append(f"unreadable journal {path.name}")
    return report


def repair(cache_path: Path, report: Optional[HealthReport] = None,
           *, level: str = AUTO) -> RepairResult:
    """Fix what *report* found. Idempotent; safe to run on a healthy cache."""
    cache_path = Path(cache_path)
    result = RepairResult()
    if report is None:
        report = verify(cache_path)
    if not report.exists:
        result.failures.append("no cache at this path")
        return result

    # 0. A root zarr.json that is not valid JSON blocks everything below —
    #    consolidate_metadata has to open the store to walk it. The root group
    #    carries no data of its own, only the spatialdata version stamp and the
    #    consolidated index, so rebuilding it loses nothing that step 4 does not
    #    then regenerate from what is actually on disk.
    if not report.readable_metadata:
        try:
            _write_minimal_root(cache_path)
            result.actions.append("rebuilt the unreadable root metadata")
        except Exception as exc:
            result.failures.append(f"could not rebuild root metadata: {exc}")
            return result

    # 1. Finish or unwind interrupted writes, and clear debris. Handles the
    #    *.partial files and abandoned staging trees too.
    try:
        result.actions += recover_pending(cache_path)
    except Exception as exc:
        result.failures.append(f"could not replay pending writes: {exc}")

    # 2. Drop groups that are not valid elements.
    for element in report.stray_elements:
        try:
            shutil.rmtree(cache_path / element, ignore_errors=True)
            result.actions.append(f"removed invalid element {element}")
        except Exception as exc:  # pragma: no cover - defensive
            result.failures.append(f"could not remove {element}: {exc}")

    # 3. Restore anything the metadata expects but disk lacks. Only at FULL,
    #    because it changes data rather than bookkeeping.
    if level == FULL:
        for element in report.missing_on_disk:
            backups = report.trash_available.get(element)
            if not backups:
                result.failures.append(f"{element} is missing and has no backup")
                continue
            try:
                restore_from_trash(cache_path, element, backups[0], consolidate_after=False)
                result.actions.append(f"restored {element} from backup")
            except Exception as exc:
                result.failures.append(f"could not restore {element}: {exc}")

    # 4. Re-consolidate. This alone fixes every missing_in_meta case, which is
    #    the most common corruption and the one that used to cost a rebuild.
    try:
        consolidate(cache_path)
        if report.missing_in_meta:
            result.actions.append(
                f"re-consolidated metadata ({len(report.missing_in_meta)} "
                f"element(s) were missing from it)")
        elif not result.actions:
            result.actions.append("re-consolidated metadata")
    except Exception as exc:
        result.failures.append(f"could not consolidate metadata: {exc}")
    return result


def _current_container_attrs() -> dict:
    """The root ``spatialdata_attrs`` stamp the installed spatialdata writes."""
    try:
        import spatialdata
        from spatialdata._io.format import CurrentSpatialDataContainerFormat
        fmt = CurrentSpatialDataContainerFormat()
        return {
            "version": str(fmt.spatialdata_format_version),
            "spatialdata_software_version": spatialdata.__version__,
        }
    except Exception:  # pragma: no cover - upstream layout change
        return {"version": "0.2"}


def _write_minimal_root(cache_path: Path) -> None:
    """Replace an unparseable root ``zarr.json`` with a bare valid group.

    Keeps whatever ``spatialdata_attrs`` version stamp can still be salvaged, so
    a reader does not fall back to defaults for a store written by a different
    spatialdata release.
    """
    from palms.utils.zarr_safe import atomic_json

    root = cache_path / "zarr.json"
    attributes: dict = {}
    try:
        salvaged = json.loads(root.read_text()).get("attributes", {})
        if isinstance(salvaged, dict) and "spatialdata_attrs" in salvaged:
            attributes = salvaged
    except (OSError, ValueError):
        pass
    if not attributes:
        # Do NOT salvage an element's spatialdata_attrs here — the container
        # and element stamps use different schemas, and mixing them makes the
        # store unreadable in a new way. Fall back to the installed format.
        attributes = {"spatialdata_attrs": _current_container_attrs()}
    atomic_json(root, {
        "attributes": attributes,
        "zarr_format": 3,
        "node_type": "group",
    })


def restore_from_trash(cache_path: Path, element: str, backup: Path,
                       *, consolidate_after: bool = True) -> None:
    """Move a ``.xv_trash`` copy of *element* back into place."""
    from palms.utils.zarr_safe import store_lock

    cache_path = Path(cache_path)
    live = cache_path / element
    with store_lock(cache_path):
        if live.exists():
            stamp = backup.name.rsplit(".", 2)[-1]
            shutil.move(str(live), str(
                cache_path / TRASH_DIR / "_replaced" / f"{Path(element).name}.{stamp}"))
        live.parent.mkdir(parents=True, exist_ok=True)
        os.rename(backup, live)
        if consolidate_after:
            consolidate(cache_path)


def salvageable_elements(backup: Path) -> list[str]:
    """Element paths in *backup* that look intact, without opening the store.

    Deliberately filesystem-level: a cache worth recovering from is one that
    failed to open, so anything that starts by reading it whole is useless here.
    Element directories are self-contained, so a broken root index or an
    unreadable table does not condemn the shapes and images beside it.
    """
    backup = Path(backup)
    return [e for e in _disk_elements(backup) if not _is_stray(backup, e)]


def read_obs_columns(cache_path: Path, prefixes: tuple[str, ...],
                     table: str = "tables/table") -> dict:
    """Read user obs columns straight out of zarr, bypassing anndata.

    This is how a clustering survives a table that ``_read_table`` refuses —
    each obs column is its own zarr array or categorical group, so losing the
    table's version attr (or anything else at table level) does not lose them.
    Returns ``{column_name: (index_array, values_array)}``.
    """
    import numpy as np
    import zarr

    cache_path = Path(cache_path)
    obs_dir = cache_path / table / "obs"
    if not obs_dir.is_dir():
        return {}

    def _read(path: Path):
        """Decode one obs column.

        anndata uses several on-disk encodings and they are not interchangeable:
        a plain array, a ``categorical`` group of categories+codes, and the
        ``nullable-*`` groups of values+mask. Assuming "array" silently returned
        nothing for the nullable index, which made every column unreadable.
        """
        marker = path / "zarr.json"
        if not marker.exists():
            return None
        try:
            document = json.loads(marker.read_text())
        except (OSError, ValueError):
            return None
        encoding = document.get("attributes", {}).get("encoding-type") or ""
        try:
            if encoding == "categorical":
                categories = np.asarray(zarr.open_array(str(path / "categories"), mode="r")[:])
                codes = np.asarray(zarr.open_array(str(path / "codes"), mode="r")[:])
                categories = np.asarray(categories.tolist(), dtype=object)
                out = np.full(len(codes), "", dtype=object)
                valid = codes >= 0
                out[valid] = categories[codes[valid]]
                return out
            if encoding.startswith("nullable"):
                values = np.asarray(zarr.open_array(str(path / "values"), mode="r")[:])
                mask_path = path / "mask"
                if mask_path.exists():
                    mask = np.asarray(zarr.open_array(str(mask_path), mode="r")[:])
                    values = values.astype(object)
                    values[mask.astype(bool)] = None
                return values
            if encoding in ("string-array", "array") or not encoding:
                values = np.asarray(zarr.open_array(str(path), mode="r")[:])
                # numpy 2 returns StringDType for string arrays, which several
                # downstream casts refuse; object dtype behaves everywhere.
                if values.dtype.kind in ("T", "U", "S", "O"):
                    values = np.asarray(values.tolist(), dtype=object)
                return values
            if document.get("node_type") == "group":
                # Unknown group encoding — try the conventional payload name.
                if (path / "values").exists():
                    return np.asarray(zarr.open_array(str(path / "values"), mode="r")[:])
                return None
            return np.asarray(zarr.open_array(str(path), mode="r")[:])
        except Exception as exc:
            log.debug("could not read obs column %s: %s", path.name, exc)
            return None

    index = _read(obs_dir / "_index")
    if index is None:
        return {}

    found: dict = {}
    for entry in sorted(obs_dir.iterdir()):
        if entry.name.startswith("_") or not entry.is_dir():
            continue
        if not entry.name.startswith(prefixes):
            continue
        values = _read(entry)
        if values is not None and len(values) == len(index):
            found[entry.name] = (index, values)
    return found


def describe_store(cache_path: Path) -> dict:
    """Size and free-space figures for the health panel."""
    cache_path = Path(cache_path)
    total = 0
    for root, dirs, files in os.walk(cache_path):
        dirs[:] = [d for d in dirs if d not in INTERNAL_NAMES]
        for name in files:
            try:
                total += os.path.getsize(os.path.join(root, name))
            except OSError:
                pass
    usage = shutil.disk_usage(cache_path if cache_path.exists() else cache_path.parent)
    return {
        "size_bytes": total,
        "free_bytes": usage.free,
        "total_bytes": usage.total,
    }


def human_bytes(value: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024 or unit == "TB":
            return f"{value:.0f} {unit}" if unit == "B" else f"{value:.1f} {unit}"
        value /= 1024
    return f"{value:.1f} TB"
