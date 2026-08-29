"""Every third-party module `src/palms` imports must be declared in pyproject.toml.

`novae` was imported by the Domains tab and appeared in no dependency list at all —
it worked only because nobody had installed PALMS anywhere that did not already
have it. The same held for `platformdirs` and `pygments`, which arrived as
transitive dependencies of napari: a `pip install palms` that resolved a napari
release which had dropped either one would fail at import, on a line that names
neither napari nor PALMS.

So this walks the AST rather than trusting a list, and maps each imported module
back to the distribution that provides it (an import name is not a package name:
`cv2` comes from `opencv-python`, `sklearn` from `scikit-learn`).
"""

from __future__ import annotations

import ast
import sys
import tomllib
from importlib.metadata import packages_distributions
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
SRC = REPO / "src" / "palms"

# Import name -> declared distribution, for the cases the installed metadata
# cannot bridge: conda-forge ships OpenCV as `py-opencv`, whose dist-info is
# named `cv2`, so nothing links the import back to the `opencv-python` we declare.
ALIASES = {"cv2": "opencv-python"}

# Import names that are deliberately undeclared, each with the reason.
EXEMPT = {
    # Ships inside matplotlib, which is declared; it has no distribution of its own.
    "mpl_toolkits": "part of matplotlib",
    # Optional and guarded by try/except in adata_persistence: xarray absorbed
    # DataTree, so the standalone package must not become a requirement.
    "datatree": "optional, guarded, superseded by xarray.DataTree",
}


def _normalise(name: str) -> str:
    return name.lower().replace("_", "-").replace(".", "-")


def _declared_distributions() -> set[str]:
    """Every distribution named in [project] dependencies or any extra."""
    meta = tomllib.loads((REPO / "pyproject.toml").read_text())["project"]
    specs = list(meta.get("dependencies", []))
    for extra in meta.get("optional-dependencies", {}).values():
        specs.extend(extra)

    out = set()
    for spec in specs:
        # "insitucnv @ git+https://…", "palms[test]", "zarr>=3.0,<3.2", "pytest"
        name = spec.split("@")[0].split("[")[0]
        for sep in ("<", ">", "=", "!", "~", ";", " "):
            name = name.split(sep)[0]
        if name.strip():
            out.add(_normalise(name.strip()))
    return out


def _imported_top_level_modules() -> dict[str, set[Path]]:
    """Top-level module name -> the files importing it (absolute imports only)."""
    found: dict[str, set[Path]] = {}
    for path in sorted(SRC.rglob("*.py")):
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                names = [node.module]
            else:
                continue
            for name in names:
                found.setdefault(name.split(".")[0], set()).add(path.relative_to(REPO))
    return found


def test_every_third_party_import_is_declared():
    declared = _declared_distributions()
    provided_by = packages_distributions()

    undeclared = {}
    for module, files in _imported_top_level_modules().items():
        if module in sys.stdlib_module_names or module == "palms" or module in EXEMPT:
            continue

        dists = {_normalise(d) for d in provided_by.get(module, ())}
        # Not installed here: fall back to the import name, which is all we have.
        if not dists:
            dists = {_normalise(module)}
        if module in ALIASES:
            dists.add(_normalise(ALIASES[module]))

        if not (dists & declared):
            undeclared[module] = (sorted(dists), sorted(str(f) for f in files)[:3])

    assert not undeclared, (
        "imported but declared in no pyproject dependency list — a fresh install "
        "would fail on these:\n"
        + "\n".join(
            f"  {mod}: provided by {dists}, imported in {files}"
            for mod, (dists, files) in sorted(undeclared.items())
        )
    )


@pytest.mark.parametrize("module,reason", sorted(EXEMPT.items()))
def test_exemptions_are_still_imported(module, reason):
    """An exemption that no longer describes real code should be deleted, not kept."""
    assert module in _imported_top_level_modules(), (
        f"{module} is exempted ({reason}) but nothing in src/palms imports it any more"
    )
