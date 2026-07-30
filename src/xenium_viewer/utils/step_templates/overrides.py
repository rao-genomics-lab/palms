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
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import Optional

from xenium_viewer.utils.step_templates.loader import (
    builtin_spec,
    clear_cache,
    user_template_dir,
)
from xenium_viewer.utils.step_templates.spec import BLOCK_MARKER

_HEADER = "# xenium-viewer template"


def override_path(template_id: str, directory: Optional[Path] = None) -> Path:
    """Where *template_id*'s override lives (whether or not it exists)."""
    return (directory or user_template_dir()) / f"{template_id}.tmpl"


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
    text = render_override(template_id, differing)

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
    except OSError:
        return False
    finally:
        clear_cache()
        _forget_rejection(template_id)
    return existed


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
