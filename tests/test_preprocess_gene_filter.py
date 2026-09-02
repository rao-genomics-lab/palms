"""palms-preprocess must cache the panel's genes, not everything in the parquet.

The filter used to be `if "is_gene" in df.columns: df = df[df["is_gene"]]`.
XOA 2.0.0 dropped that column from `transcripts.parquet`, so on a modern bundle
the guard turned the filter into a silent pass and every control codeword got
its own feather file — measured on Xenium_V1_human_Pancreas_FFPE, 514 files
where the panel has 377 genes. Nothing failed, and the step printed
"Filtering: is_gene=True" while doing no such thing.

The gene set now comes from `gene_panel.json`, whose per-target `descriptor`
does not depend on the format version.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

pd = pytest.importorskip("pandas")

from palms.preprocess import panel_genes, preprocess  # noqa: E402

GENES = ["AAA", "BBB"]
CONTROLS = ["NegControlCodeword_0500", "UnassignedCodeword_0001"]


def _panel(path: Path, genes=GENES, controls=CONTROLS) -> None:
    targets = [{"type": {"descriptor": "gene", "data": {"name": g, "id": f"E{g}"}}}
               for g in genes]
    targets += [{"type": {"descriptor": "negative_control", "data": {"name": c}}}
                for c in controls]
    path.write_text(json.dumps({"payload": {"targets": targets}}))


def _parquet(path: Path, *, is_gene: bool) -> None:
    """A miniature transcripts.parquet, with or without the dropped column."""
    rows = []
    for name in GENES + CONTROLS:
        for i in range(10):
            rows.append({"x_location": float(i), "y_location": float(i),
                         "qv": 40.0 if i else 5.0, "feature_name": name,
                         "cell_id": f"c{i}"})
    df = pd.DataFrame(rows)
    if is_gene:
        df["is_gene"] = ~df["feature_name"].isin(CONTROLS)
    df.to_parquet(path)


def _run(tmp_path: Path, *, is_gene: bool, panel: bool) -> set[str]:
    _parquet(tmp_path / "transcripts.parquet", is_gene=is_gene)
    if panel:
        _panel(tmp_path / "gene_panel.json")
    cache = tmp_path / "transcript_cache"
    preprocess(parquet_path=tmp_path / "transcripts.parquet", cache_dir=cache,
               genes=panel_genes(tmp_path))
    return {p.stem for p in cache.glob("*.feather")}


@pytest.mark.parametrize("is_gene", [True, False])
def test_only_panel_genes_are_cached(tmp_path, is_gene):
    """The answer must not depend on the column that XOA 2.0.0 removed."""
    assert _run(tmp_path, is_gene=is_gene, panel=True) == set(GENES)


def test_the_quality_floor_still_applies(tmp_path):
    _run(tmp_path, is_gene=False, panel=True)
    kept = pd.read_feather(tmp_path / "transcript_cache" / "AAA.feather")
    assert len(kept) == 9          # one of the ten rows is qv 5
    assert kept["qv"].min() >= 20


def test_without_a_panel_nothing_is_dropped_by_name(tmp_path):
    """No gene_panel.json means the controls cannot be told apart.

    Caching them is wasteful but harmless — the viewer's gene list comes from
    `adata.var_names`, not from the cache — so this stays a warning rather than
    a refusal. What must not happen is a *claim* to have filtered.
    """
    assert _run(tmp_path, is_gene=False, panel=False) == set(GENES + CONTROLS)


def test_panel_genes_reads_only_the_gene_descriptor(tmp_path):
    _panel(tmp_path / "gene_panel.json")
    assert panel_genes(tmp_path) == set(GENES)


@pytest.mark.parametrize("body", ["", "not json", '{"payload": {}}',
                                  '{"payload": {"targets": []}}'])
def test_an_unusable_panel_is_None_rather_than_an_exception(tmp_path, body):
    """A malformed panel must not take the whole preprocessing run down."""
    (tmp_path / "gene_panel.json").write_text(body)
    assert panel_genes(tmp_path) is None


def test_a_missing_panel_is_None(tmp_path):
    assert panel_genes(tmp_path) is None
