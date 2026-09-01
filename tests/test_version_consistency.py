"""One version, stated in four places.

`pyproject.toml` carries a static version, `palms.__version__` a second copy,
`CITATION.cff` a third — read by Zenodo when it mints the DOI — and `CHANGELOG.md`
names the release in its newest heading. Nothing derives any of them from another,
so a release that bumps three of the four ships a package whose metadata disagrees
with itself, and a DOI whose version field is wrong. That is not fixable after the
fact: a minted DOI is permanent.

Pure stdlib and regex — no build backend, no import of the package.
"""
from __future__ import annotations

import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent


def _pyproject() -> str:
    text = (REPO / "pyproject.toml").read_text()
    return re.search(r'^version = "([^"]+)"', text, re.M).group(1)


def _dunder() -> str:
    text = (REPO / "src" / "palms" / "__init__.py").read_text()
    return re.search(r'^__version__ = "([^"]+)"', text, re.M).group(1)


def _citation() -> str:
    text = (REPO / "CITATION.cff").read_text()
    return re.search(r"^version: (\S+)", text, re.M).group(1)


def _changelog() -> str:
    """The newest release heading. `## [Unreleased]` is skipped: between releases
    it sits above the newest version and is not itself a version."""
    text = (REPO / "CHANGELOG.md").read_text()
    for match in re.finditer(r"^## \[([^\]]+)\]", text, re.M):
        if match.group(1) != "Unreleased":
            return match.group(1)
    raise AssertionError("CHANGELOG.md has no release heading")


def test_the_package_and_its_metadata_agree():
    assert _pyproject() == _dunder() == _citation()


def test_the_changelog_names_the_version_being_shipped():
    """A release whose changelog still says the previous version leaves a reader
    with no record of what changed — and `release.yml` publishes on the tag, so
    nothing else would notice."""
    assert _changelog() == _pyproject()


def test_the_version_is_a_release_number():
    assert re.fullmatch(r"\d+\.\d+\.\d+", _pyproject())
