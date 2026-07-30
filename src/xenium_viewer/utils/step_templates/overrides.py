"""Writing and removing user overrides.

Saving deliberately **never refuses**. A user editing their own template file is
mid-thought; refusing to write it sends them to an external editor and out of
the loop that gives them feedback. So the text is written whatever its state,
and *activation* is the thing that is gated — an invalid file on disk is simply
rejected by :func:`~.loader.resolve`, which falls back to the shipped template
and says so. "Saved but not in effect" needs no special mechanism; it is what
the resolver already does.

Only blocks that actually differ from the shipped text are written. That is not
tidiness — it is the mechanism that makes upgrades survivable. A file recording
one changed block leaves the rest resolving against whatever the current release
ships, so a later fix to a block the user never touched still reaches them.

That also settles most of the classic conffile problem for free. ``dpkg`` has to
ask "did the user modify this file?" because it manages whole files; here an
unmodified block is simply *absent*, so it tracks upstream with nothing to
decide. What is left is the genuinely hard case — a block the user *did* change,
whose shipped version has since changed too — and telling that apart from an
ordinary customisation is exactly what :data:`MANIFEST_NAME` exists for. It
records the hash of the **shipped** text each block was forked from, so a later
release can be compared against the thing the user actually diverged from rather
than against their own edit.
"""

from __future__ import annotations

import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from xenium_viewer.utils.step_templates.loader import (
    builtin_spec,
    clear_cache,
    user_template_dir,
)
from xenium_viewer.utils.step_templates.spec import BLOCK_MARKER

_HEADER = "# xenium-viewer template"

#: Bookkeeping for the upgrade check, kept beside the ``.tmpl`` files rather than
#: inside them: a user editing their own template must not be able to corrupt
#: the record of what they forked from, which is the only way to distinguish
#: "upstream moved" from "I changed this".
MANIFEST_NAME = "overrides.json"

#: Bumped only if the manifest's own shape changes.
MANIFEST_VERSION = 1


def override_path(template_id: str, directory: Optional[Path] = None) -> Path:
    """Where *template_id*'s override lives (whether or not it exists)."""
    return (directory or user_template_dir()) / f"{template_id}.tmpl"


def manifest_path(directory: Optional[Path] = None) -> Path:
    return (directory or user_template_dir()) / MANIFEST_NAME


def read_manifest(directory: Optional[Path] = None) -> dict:
    """The fork record, or an empty one.

    Never raises. A corrupt manifest costs the upgrade *warning*, not the
    override itself — losing bookkeeping must not lose the user's work, and a
    hand-edited or truncated file is a plausible state to find it in.
    """
    try:
        data = json.loads(manifest_path(directory).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    if not isinstance(data, dict) or data.get("version") != MANIFEST_VERSION:
        return {}
    templates = data.get("templates")
    return templates if isinstance(templates, dict) else {}


def _write_manifest(templates: dict, directory: Optional[Path] = None) -> None:
    path = manifest_path(directory)
    if not templates:
        try:
            path.unlink(missing_ok=True)
        except OSError:
            pass
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"version": MANIFEST_VERSION, "templates": templates}
    _atomic_write(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def record_fork(template_id: str, blocks: dict,
                directory: Optional[Path] = None) -> None:
    """Note which shipped text each overridden block was forked from.

    The hashes are of the **builtin** block, not the user's — the question a
    later release has to answer is "has the thing they diverged from moved?",
    and hashing their own edit could not answer it.
    """
    builtin = builtin_spec(template_id)
    templates = read_manifest(directory)
    templates[template_id] = {
        "based_on_schema": builtin.schema_version,
        "forked_from": {
            name: builtin.blocks[name].hash()
            for name in blocks if name in builtin.blocks
        },
        "saved_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    _write_manifest(templates, directory)


def forget_fork(template_id: str, directory: Optional[Path] = None) -> None:
    """Drop *template_id* from the manifest — it is no longer customised."""
    templates = read_manifest(directory)
    if templates.pop(template_id, None) is not None:
        _write_manifest(templates, directory)


def changed_blocks(template_id: str, blocks: dict) -> dict:
    """The subset of *blocks* whose text differs from the shipped version."""
    builtin = builtin_spec(template_id)
    return {
        name: text for name, text in blocks.items()
        if name not in builtin.blocks or text != builtin.blocks[name].text
    }


def render_override(template_id: str, blocks: dict) -> str:
    """The ``.tmpl`` file contents for an override supplying *blocks*."""
    builtin = builtin_spec(template_id)
    lines = [
        _HEADER,
        f"# id: {template_id}",
        f"# schema-version: {builtin.schema_version}",
        "#",
        "# Customised copy. Only the blocks below override the shipped template;",
        "# every other block continues to track the version you have installed,",
        "# so fixes to them still reach you.",
        "",
    ]
    for name, text in blocks.items():
        lines.append(f"{BLOCK_MARKER}{name}")
        lines.append(text.strip("\n"))
        lines.append("")
    return "\n".join(lines).rstrip("\n") + "\n"


def save_override(template_id: str, blocks: dict,
                  directory: Optional[Path] = None) -> Optional[Path]:
    """Write an override for *template_id*. Returns the path, or None if reverted.

    Passing blocks that all match the shipped text removes the override rather
    than writing a file that says nothing — an inert file left behind would show
    up as "customised" in the UI forever.

    The write is atomic: a crash mid-save must not leave a truncated template
    that then fails to parse on next launch, which would look to the user like
    their edit corrupted something.
    """
    differing = changed_blocks(template_id, blocks)
    if not differing:
        remove_override(template_id, directory)
        return None

    path = override_path(template_id, directory)
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        _atomic_write(path, render_override(template_id, differing))
        # After the template, so a crash between the two leaves an override with
        # no fork record — which reads as "never reviewed" and prompts, rather
        # than as "reviewed and current", which would suppress a real warning.
        record_fork(template_id, differing, directory)
    finally:
        clear_cache()
        _forget_rejection(template_id)
    return path


def remove_override(template_id: str,
                    directory: Optional[Path] = None) -> bool:
    """Delete *template_id*'s override. Returns whether there was one."""
    path = override_path(template_id, directory)
    try:
        existed = path.exists()
        if existed:
            path.unlink()
        forget_fork(template_id, directory)
    except OSError:
        return False
    finally:
        clear_cache()
        _forget_rejection(template_id)
    return existed


def _atomic_write(path: Path, text: str) -> None:
    """Write *text* to *path* atomically.

    A crash mid-save must not leave a truncated template that fails to parse on
    next launch — to the user that looks like their edit corrupted something.
    """
    handle, tmp = tempfile.mkstemp(dir=str(path.parent), suffix=".tmp")
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as fh:
            fh.write(text)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def _forget_rejection(template_id: str) -> None:
    """Let a re-saved template report again.

    Rejections are remembered per session so a bad override warns once instead
    of on every step. Saving is the event that makes the old verdict stale.
    """
    try:
        from xenium_viewer.utils import reporting
        reporting.template_rejections().pop(template_id, None)
        reporting._template_rejections.pop(template_id, None)
    except Exception:  # pragma: no cover - reporting is best-effort
        pass
