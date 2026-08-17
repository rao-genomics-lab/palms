"""Process memory probes, for finding out where RAM actually goes.

Written because "RSS keeps climbing during the zarr write" is not something you
can answer by reading code. Every function here is best-effort and read-only:
nothing raises, because a probe that can take the viewer down is worse than no
probe at all.

Two things worth knowing before reading a number out of here:

* **Freeing memory in Python does not lower RSS.** glibc keeps freed blocks in
  per-thread arenas rather than returning them to the kernel, so RSS ratchets to
  the high-water mark and stays there. That alone makes a per-element write loop
  look like it leaks. :func:`malloc_trim` is the counter-measure — call it after
  dropping large objects, then re-read :func:`rss_bytes`, and the difference
  tells you how much was merely un-returned rather than still referenced.
* **napari installs a global dask cache.** ``napari.layers.base.Layer.__init__``
  calls ``configure_dask(data, cache=True)`` for any dask-backed layer, which
  sizes ``dask.cache.Cache`` at 25% of total RAM (~31 GB on a 125 GB box) and
  keeps every task result it sees. It is only registered while napari slices,
  but ``dask.callbacks`` registration is process-global, so a slice overlapping
  a long write captures that write's chunks too. :func:`dask_cache_bytes`
  reports what it is holding.
"""

from __future__ import annotations

import ctypes
import gc
from typing import Callable, Optional

_PAGE_SIZE = 4096

# Resolved once. None means "we looked and there is no usable malloc_trim", which
# is different from "not looked at yet" and stops us re-dlopening libc per call.
_libc: Optional[ctypes.CDLL] = None
_libc_resolved = False


def rss_bytes() -> Optional[int]:
    """Current resident set size, or None if it cannot be read."""
    try:
        with open("/proc/self/statm", "rb") as fh:
            return int(fh.read().split()[1]) * _PAGE_SIZE
    except Exception:
        return None


def peak_rss_bytes() -> Optional[int]:
    """High-water RSS (``VmHWM``), or None if it cannot be read.

    This is the number that matters for "did we survive": it is what the OOM
    killer saw, not what is resident once the peak has passed.
    """
    try:
        with open("/proc/self/status", "rb") as fh:
            for line in fh:
                if line.startswith(b"VmHWM:"):
                    return int(line.split()[1]) * 1024
    except Exception:
        pass
    return None


def _get_libc() -> Optional[ctypes.CDLL]:
    global _libc, _libc_resolved
    if not _libc_resolved:
        _libc_resolved = True
        try:
            lib = ctypes.CDLL("libc.so.6")
            lib.malloc_trim  # raises AttributeError on musl and other non-glibc
            _libc = lib
        except Exception:
            _libc = None
    return _libc


def malloc_trim() -> bool:
    """Ask glibc to return free heap memory to the OS. True if it released any.

    Non-glibc platforms have no ``malloc_trim``; there we return False and the
    caller carries on. This is a hint, not a guarantee — memory that is
    genuinely still referenced obviously stays put, which is exactly what makes
    the RSS delta across this call informative.
    """
    lib = _get_libc()
    if lib is None:
        return False
    try:
        return bool(lib.malloc_trim(0))
    except Exception:
        return False


def release(collect: bool = True) -> None:
    """Drop what Python can drop, then hand the freed pages back to the OS."""
    if collect:
        gc.collect()
    malloc_trim()


def dask_cache_bytes() -> Optional[int]:
    """Bytes held by napari's global opportunistic dask cache, if there is one.

    Reaches into ``napari.utils._dask_utils._DASK_CACHE``, which is private and
    may move; every failure mode returns None rather than raising. A return of
    ``1`` (not None) means the cache object exists but was never resized, i.e.
    no dask-backed layer has been added yet.
    """
    try:
        from napari.utils import _dask_utils
        return int(_dask_utils._DASK_CACHE.cache.total_bytes)
    except Exception:
        return None


def dask_cache_limit_bytes() -> Optional[int]:
    """Size napari's global dask cache has been allowed to grow to."""
    try:
        from napari.utils import _dask_utils
        return int(_dask_utils._DASK_CACHE.cache.available_bytes)
    except Exception:
        return None


def _gb(value: Optional[int]) -> str:
    return "?" if value is None else f"{value / 1e9:.2f}GB"


def format_memory(tag: str = "") -> str:
    """One-line memory summary: current RSS, peak RSS, and the dask cache."""
    parts = [f"rss={_gb(rss_bytes())}", f"peak={_gb(peak_rss_bytes())}"]
    cached = dask_cache_bytes()
    # 1 byte is dask.cache.Cache's unresized default, not a real cache.
    if cached is not None and cached > 1:
        parts.append(f"dask_cache={_gb(cached)}/{_gb(dask_cache_limit_bytes())}")
    prefix = f"[mem] {tag}: " if tag else "[mem] "
    return prefix + " ".join(parts)


def log_memory(tag: str = "", log: Optional[Callable[[str], None]] = None) -> str:
    """Emit :func:`format_memory` through *log* (default: the package logger)."""
    line = format_memory(tag)
    if log is None:
        from xenium_viewer.utils import reporting
        reporting.get_logger().info(line)
    else:
        log(line)
    return line
