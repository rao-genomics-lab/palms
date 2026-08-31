#!/usr/bin/env bash
#
# Stage a publishable copy of a dataset for scripts/capture_screenshots.py.
#
#     ./scripts/prepare_demo_dataset.sh /path/to/source/dataset [/tmp/Crop_6]
#     python scripts/capture_screenshots.py /tmp/Crop_6
#
# The screenshots are published, and two panels (Tools > Dataset, Tools > Cache)
# print the dataset's absolute path into the widget, so the path in the copy is
# the path on the wiki — hence a copy at a neutral location rather than a
# capture run against a working dataset. Three things have to happen to it, and
# the third is the one that gets forgotten:
#
#   1. copy, and drop the legacy xenium_viewer.log the Cache tab lists by name;
#   2. repoint the absolute paths recorded in the provenance graph, or the
#      Notebook tab and analysis.py show the old location — and the first launch
#      after the move marks every node stale for nothing;
#   3. replace the ARMS scan's filename, which reaches the napari layer list and
#      is therefore in *every* full-window screenshot. A real scan filename is a
#      slide identifier; this is the leak that is invisible until it is on a
#      published page.
#
# The capture run mutates the copy — it runs Leiden, ranks genes, computes ROI
# and ARMS DEG, draws annotations, and replaces the H&E's fine registration with
# a coarse alignment — so re-run this before re-capturing rather than shooting a
# used copy: the registration steps would otherwise photograph the coarse
# alignment, and the annotation steps would skip drawing because shapes already
# exist.
#
# Deliberately POSIX-ish bash, matching install.sh: macOS still ships bash 3.2.
set -euo pipefail

repo="$(cd "$(dirname "$0")/.." && pwd)"
src="${1:-${PALMS_DEMO_SOURCE:-}}"
dest="${2:-${PALMS_DEMO_DEST:-/tmp/Crop_6}}"
arms_name="${PALMS_DEMO_ARMS_NAME:-demo_he_section.svs}"

if [ -z "$src" ]; then
    echo "usage: $0 <source-dataset> [destination]" >&2
    echo "       (or set PALMS_DEMO_SOURCE; destination defaults to /tmp/Crop_6)" >&2
    exit 2
fi
if [ ! -f "$src/experiment.xenium" ]; then
    echo "error: $src is not a Xenium dataset (no experiment.xenium)" >&2
    exit 1
fi

# Resolve both sides before comparing: a destination that IS the source, or
# inside it, would have this delete the original.
abs_src="$(cd "$src" && pwd -P)"
abs_dest_parent="$(cd "$(dirname "$dest")" && pwd -P)"
abs_dest="$abs_dest_parent/$(basename "$dest")"
case "$abs_src" in
    "$abs_dest"|"$abs_dest"/*)
        echo "error: destination $abs_dest is, or contains, the source" >&2
        exit 1 ;;
esac
# And the other direction, which is the one that loses data: a destination
# *inside* the source is removed by the step below, and it is the source's own
# subdirectory being removed. cp refuses this too, but only after the rm.
case "$abs_dest" in
    "$abs_src"/*)
        echo "error: destination $abs_dest is inside the source" >&2
        exit 1 ;;
esac

# The only thing this script may delete is a previous copy of itself. Anything
# else at that path is someone's data, and refusing is the whole safety story.
if [ -e "$abs_dest" ]; then
    if [ ! -f "$abs_dest/experiment.xenium" ]; then
        echo "error: $abs_dest exists and is not a Xenium dataset — refusing to remove it" >&2
        exit 1
    fi
    echo "Removing the previous copy at $abs_dest"
    rm -rf "$abs_dest"
fi

echo "Copying $abs_src -> $abs_dest"
cp -a "$abs_src" "$abs_dest"
rm -f "$abs_dest/xenium_viewer.log"

echo "Repairing the recorded paths"
if command -v palms-rename-dataset >/dev/null 2>&1; then
    palms-rename-dataset "$abs_dest" --repair
else
    # Not installed (a checkout with no `pip install -e .`): the console script
    # is palms.scripts.rename_dataset:main, so call that module directly.
    PYTHONPATH="$repo/src${PYTHONPATH:+:$PYTHONPATH}" \
        python -m palms.scripts.rename_dataset "$abs_dest" --repair
fi

echo "Replacing the ARMS scan filename with $arms_name"
PYTHONPATH="$repo/src${PYTHONPATH:+:$PYTHONPATH}" python - "$abs_dest" "$arms_name" <<'PY'
import sys
from pathlib import Path

from palms.utils.zarr_safe import safe_group_update

dest, name = Path(sys.argv[1]), sys.argv[2]
store = dest / "sdata_cached.zarr"
if not (store / "viewer_session").exists():
    print("  no viewer_session in the cache — nothing to rename")
    raise SystemExit(0)
# safe_group_update, not a hand-written zarr.json: the store must never be
# written behind the viewer's back, even one that only lives in /tmp.
with safe_group_update(store, "viewer_session") as (group, _staging):
    had = group.attrs.get("arms_he_filename") is not None
    group.attrs["arms_he_filename"] = name
# Deliberately does not echo the previous value: it is the slide identifier this
# step exists to remove, and a pasted log is a published log.
print(f"  arms_he_filename {'replaced' if had else 'set'}: {name!r}")
PY

echo
echo "Ready: $abs_dest"
echo "  python scripts/capture_screenshots.py $abs_dest"
