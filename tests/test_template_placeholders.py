"""Every ``$token`` in a template must be a declared param.

``Step.render`` uses ``Template.substitute`` (not ``safe_substitute``), so a
template carrying a ``$token`` that no call site declares is a hard ``StepError``
at record time. That is the design: fail loudly when the code is recorded rather
than with a ``NameError`` when someone replays the notebook.

Two templates used to defeat that check by ``str.replace``-ing a fake token out
of the text *before* ``Step`` ever saw it — ``$n_suffix`` in the gene-correlation
tail and ``$dpi_kwarg`` in the marker-plot tail. Neither could ever become a real
param (``$n_suffix`` sits inside an f-string, where ``repr('')`` renders as ``''``
and breaks the literal), so the token was load-bearing punctuation that looked
exactly like a parameter. Both are now whole-line block variants.

This is a source guard, in the same spirit as the ones in
``test_persistence_safety.py``: it fails if the idiom is reintroduced.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path
from string import Template

import pytest

_TABS = Path(__file__).resolve().parent.parent / "src" / "palms" / "tabs"
_TAB_SOURCES = sorted(_TABS.glob("*.py"))


def _placeholders(text: str) -> set[str]:
    """Names ``Template.substitute`` would demand from *text*.

    Hand-rolled off ``Template.pattern`` rather than using
    ``Template.get_identifiers()``, which is Python 3.11+ while the package
    declares ``requires-python = ">=3.10"``. Using the stdlib method would work
    in the dev env and fail for anyone on the declared floor.
    """
    names = set()
    for match in Template.pattern.finditer(text):
        if match.group("invalid") is not None:
            raise AssertionError(f"malformed placeholder in template: {match.group(0)!r}")
        name = match.group("named") or match.group("braced")
        if name is not None:
            names.add(name)
    return names


@pytest.mark.parametrize("path", _TAB_SOURCES, ids=lambda p: p.name)
def test_no_dollar_token_is_stripped_before_rendering(path: Path):
    """No tab may ``str.replace`` a ``$token`` out of a template.

    Doing so hides the token from ``Template.substitute``, which is the only
    thing that checks a template against its declared params — so the tests that
    render every template would never exercise it.
    """
    tree = ast.parse(path.read_text())
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "replace"
                and node.args):
            continue
        first = node.args[0]
        if isinstance(first, ast.Constant) and isinstance(first.value, str):
            assert not first.value.startswith("$"), (
                f"{path.name}:{node.lineno} strips the template token "
                f"{first.value!r} before Step sees it. Make it a block variant "
                f"(two whole-line alternatives) or a real param instead."
            )


def test_the_placeholder_reader_agrees_with_substitute():
    """Guard the guard: ``_placeholders`` must match what ``substitute`` demands."""
    text = "a = $one\nb = ${two}\nc = '$$literal'\nd = {'k': $three}"
    assert _placeholders(text) == {"one", "two", "three"}
    # substitute succeeds with exactly those names and no others.
    Template(text).substitute({n: "1" for n in _placeholders(text)})


def test_a_bare_dollar_is_reported_not_ignored():
    with pytest.raises(AssertionError, match="malformed placeholder"):
        _placeholders("cost = 5 $ each")


@pytest.mark.parametrize(
    "module,fn,variants",
    [
        ("tab_gene_correlation", "_gene_corr_template",
         [(n, f) for n in ("Raw counts", "Fraction of total", "Log1p(CPM)")
          for f in (False, True)]),
        # The dpi variant is gone: the save block now loops over a ``paths``
        # list and passes dpi unconditionally, because the format is no longer
        # a per-plot choice. Two variants per plot type instead of four.
        ("tab_marker_genes", "_marker_plot_template",
         [(p, r)
          for p in ("dotplot", "heatmap", "matrixplot", "tracksplot",
                    "correlation_matrix")
          for r in (False, True)]),
    ],
)
def test_every_variant_of_the_de_hacked_templates_renders(module, fn, variants):
    """The two templates that carried fake tokens now render for every variant.

    Rendering with exactly the placeholders the template asks for is the check
    that used to be impossible: ``substitute`` raises on a token that is not
    supplied, which is what the ``str.replace`` was hiding.
    """
    mod = pytest.importorskip(f"palms.tabs.{module}")
    assemble = getattr(mod, fn)
    for variant in variants:
        template = assemble(*variant)
        names = _placeholders(template)
        assert "n_suffix" not in names and "dpi_kwarg" not in names, (
            f"{fn}{variant} still carries a fake token: {sorted(names)}"
        )
        rendered = Template(template).substitute({n: repr("x") for n in names})
        assert "$" not in re.sub(r"\$\$", "", rendered), (
            f"{fn}{variant} left a $ in the rendered source"
        )
