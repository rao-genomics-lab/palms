"""Unit tests for the pure string/dict helpers in gene_analysis.

Covers the LLM-annotation prompt builder / response parser and the cluster
label-mapping helper. No scanpy pipeline or AnnData needed.

Run with:  pytest tests/test_gene_analysis.py
"""
from __future__ import annotations

import pandas as pd
import pytest

from palms.utils.gene_analysis import (
    parse_llm_annotation_response, build_llm_annotation_prompt, _build_label_mapping,
)


def test_parse_fenced_json():
    assert parse_llm_annotation_response('```json\n{"0": "T cells"}\n```') == {"0": "T cells"}


def test_parse_fenced_without_language_tag():
    assert parse_llm_annotation_response('```\n{"1": "B cells"}\n```') == {"1": "B cells"}


def test_parse_bare_braces_amid_prose():
    stdout = 'Here is the annotation: {"2": "NK cells"} — hope that helps!'
    assert parse_llm_annotation_response(stdout) == {"2": "NK cells"}


def test_parse_raises_on_non_json():
    with pytest.raises(ValueError):
        parse_llm_annotation_response("Sorry, I could not determine cell types.")


def test_build_prompt_lists_top_genes_per_cluster():
    rank_df = pd.DataFrame({
        "group": ["0", "0", "0", "1", "1", "1"],
        "names": ["CD3D", "CD8A", "IL7R", "COL1A1", "DCN", "LUM"],
    })
    prompt = build_llm_annotation_prompt(rank_df, n_genes=2)
    assert "Cluster 0: CD3D, CD8A" in prompt
    assert "Cluster 1: COL1A1, DCN" in prompt
    assert "IL7R" not in prompt  # truncated to top 2
    assert "2 groups" in prompt


def test_label_mapping_int_keyed_labels_and_string_categories():
    mapping = _build_label_mapping(["0", "1", "2"], {0: "Tcell", 1: "Bcell"})
    assert mapping == {"0": "Tcell", "1": "Bcell", "2": "2"}


def test_label_mapping_string_keyed_labels():
    mapping = _build_label_mapping(["0", "1"], {"0": "X", "1": "Y"})
    assert mapping == {"0": "X", "1": "Y"}


def test_label_mapping_non_numeric_category():
    mapping = _build_label_mapping(["A"], {"A": "Astrocyte"})
    assert mapping == {"A": "Astrocyte"}
