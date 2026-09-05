"""The Templates tab shows what a step would run, before it runs.

The value of this view rests entirely on the preview being the *real* rendered
source rather than a reconstruction of it — a preview assembled by a second code
path would be exactly the drift the Step system exists to eliminate, just moved
somewhere less visible. So the tests here check that the preview equals
``Step.render()`` of the same template and params, and that a tab supplying live
widget values is actually consulted.

The tab is built against a stand-in context: the real ``ViewerContext`` needs a
napari viewer, and this view reads only ``ctx.state``.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from palms.utils.step_templates import Preview, builtin_ids, builtin_spec
from palms.utils.steps import Step


@pytest.fixture
def ctx():
    return SimpleNamespace(state={})


@pytest.fixture
def tab(qapp, ctx):
    from palms.tabs.tab_templates import build_tab
    widget, exports = build_tab(ctx)
    return widget, exports


def _tree(widget):
    from qtpy.QtWidgets import QTreeWidget
    return widget.findChild(QTreeWidget)


def test_every_registered_template_is_listed(tab):
    from qtpy.QtCore import Qt

    widget, _ = tab
    tree = _tree(widget)
    listed = set()
    for i in range(tree.topLevelItemCount()):
        group = tree.topLevelItem(i)
        for j in range(group.childCount()):
            listed.add(group.child(j).data(0, Qt.ItemDataRole.UserRole))
    assert listed == set(builtin_ids())


def test_templates_are_grouped_by_owner(tab):
    widget, _ = tab
    tree = _tree(widget)
    groups = {tree.topLevelItem(i).text(0) for i in range(tree.topLevelItemCount())}
    assert {"Clustering", "Genes", "Spatial", "ROI", "Annotations",
            "Transcripts", "Setup"} == groups


def test_the_preview_is_the_real_rendered_source(ctx):
    """Not a reconstruction: identical to what StepExecutor would compile."""
    from palms.tabs.tab_templates import _preview

    spec = builtin_spec("clustering.leiden")
    assembly = spec.assemblies[0]
    shown = _preview(ctx, spec, assembly)

    expected = Step(id="x", template=spec.assemble(assembly),
                    params=spec.synth_params()).render()
    assert shown.endswith(expected)
    assert shown.startswith("# preview — sample values")


def test_a_live_tab_supplies_the_real_widget_values(ctx):
    """When the owning tab registers a provider, the preview uses it."""
    from palms.tabs.tab_templates import _preview

    spec = builtin_spec("clustering.leiden")
    ctx.state["template_preview"] = {
        "clustering.leiden": lambda: Preview(
            list(spec.assemblies[0]), dict(spec.synth_params(), resolution=0.42)),
    }
    shown = _preview(ctx, spec, spec.assemblies[0])

    assert "resolution=0.42" in shown
    assert shown.startswith("# preview — current widget values")


def test_the_provider_chooses_the_blocks_not_just_the_values(ctx):
    """The half of "what will this run?" that is about shape, not numbers.

    ``_preview`` used to render ``spec.assemblies[0]`` whatever the widgets said,
    so unticking a checkbox that selects a different block moved the parameters
    while the code stayed the same. Block selection lives at the call site by
    design, which is exactly why the provider has to carry it.
    """
    from palms.tabs.tab_templates import _preview

    spec = builtin_spec("clustering.leiden")
    without_hvg = ["head", "tail"]
    assert tuple(without_hvg) in spec.assemblies

    ctx.state["template_preview"] = {
        "clustering.leiden": lambda: Preview(without_hvg, spec.synth_params()),
    }
    # Asked for an assembly that *does* include HVG selection; the provider's
    # choice has to win.
    with_hvg = next(a for a in spec.assemblies if "hvg" in a)
    shown = _preview(ctx, spec, with_hvg)

    assert "highly_variable_genes" not in shown
    assert "sc.tl.leiden(" in shown


def test_a_note_names_the_value_that_is_not_settled_yet(ctx):
    """A path the save dialog has not returned must not read as a real one."""
    from palms.tabs.tab_templates import _preview

    spec = builtin_spec("roi.export_expression")
    blocks = list(spec.blocks)
    ctx.state["template_preview"] = {
        spec.id: lambda: Preview(blocks, spec.synth_params(),
                                 note="path chosen on save"),
    }
    shown = _preview(ctx, spec, blocks)
    assert shown.startswith(
        "# preview — current widget values (path chosen on save)")


def test_a_broken_provider_does_not_break_the_view(ctx):
    """A half-built tab must degrade to sample values, not raise into the GUI."""
    from palms.tabs.tab_templates import _preview

    spec = builtin_spec("clustering.leiden")
    ctx.state["template_preview"] = {
        "clustering.leiden": lambda: (_ for _ in ()).throw(RuntimeError("not ready")),
    }
    shown = _preview(ctx, spec, spec.assemblies[0])
    assert "sc.tl.leiden(" in shown
    assert shown.startswith("# preview — sample values")


def test_a_provider_naming_an_unknown_block_reports_rather_than_raises(ctx):
    """An override may drop a block the call site still selects."""
    from palms.tabs.tab_templates import _preview

    spec = builtin_spec("clustering.leiden")
    ctx.state["template_preview"] = {
        "clustering.leiden": lambda: Preview(["head", "nonesuch"],
                                             spec.synth_params()),
    }
    shown = _preview(ctx, spec, spec.assemblies[0])
    assert shown.startswith("# could not assemble this template")
    assert "nonesuch" in shown


# ── the registry-wide gate ───────────────────────────────────────────────────

#: Templates that deliberately have no preview provider, and why. A new template
#: has to make this choice rather than defaulting silently to sample values.
_NO_PROVIDER = {
    "spatial_neighbors": (
        "n_neighs comes from whichever tab called ctx.ensure_spatial_neighbors(); "
        "the Nhood and L-R tabs each have their own slider, so a provider would "
        "have to pick one arbitrarily. Uses sample-params in the .tmpl instead."
    ),
}


def _tab_sources() -> dict[str, str]:
    """Every tab module's source, by module name."""
    import pkgutil

    from palms import tabs

    out = {}
    for info in pkgutil.iter_modules(tabs.__path__):
        path = Path(tabs.__path__[0]) / f"{info.name}.py"
        if path.exists():
            out[info.name] = path.read_text(encoding="utf-8")
    return out


def _registered_template_ids() -> set[str]:
    """Template ids some tab registers a preview provider for.

    Read out of the source rather than by building every tab: most of them need
    a real ``ViewerContext`` with a napari viewer and loaded data, which is the
    reason this file uses a stand-in context in the first place.
    """
    import ast

    found = set()
    for name, source in _tab_sources().items():
        module = ast.parse(source)
        constants = {
            target.id: node.value.value
            for node in ast.walk(module)
            if isinstance(node, ast.Assign) and isinstance(node.value, ast.Constant)
            for target in node.targets
            if isinstance(target, ast.Name) and isinstance(node.value.value, str)
        }
        for node in ast.walk(module):
            # ctx.state.setdefault("template_preview", {})[<id>] = <provider>
            if not isinstance(node, ast.Subscript):
                continue
            if "template_preview" not in ast.unparse(node.value):
                continue
            key = node.slice
            if isinstance(key, ast.Constant) and isinstance(key.value, str):
                found.add(key.value)
            elif isinstance(key, ast.Name) and key.id in constants:
                found.add(constants[key.id])
    return found


def test_every_template_has_a_provider_or_a_declared_reason():
    """The gate that makes the preview true for the whole registry.

    One tab registering a provider made the tab's headline question — what will
    this button actually run? — honest for a single template and merely
    illustrative for thirteen. This asserts the property over the registry, so
    the answer cannot quietly go back to being a sample value.
    """
    registered = _registered_template_ids()
    missing = sorted(set(builtin_ids()) - registered - set(_NO_PROVIDER))
    assert not missing, (
        f"template(s) {missing} have no preview provider. Register one in the "
        f"owning tab (see _leiden_preview in tab_clustering.py), or add the id "
        f"to _NO_PROVIDER here with the reason."
    )
    # Kept honest in the other direction too: a stale exemption for a template
    # that has since grown a provider is a comment claiming something false.
    both = sorted(registered & set(_NO_PROVIDER))
    assert not both, f"template(s) {both} are exempted but do have a provider"
    unknown = sorted(set(_NO_PROVIDER) - set(builtin_ids()))
    assert not unknown, f"_NO_PROVIDER names template(s) that do not exist: {unknown}"


def test_a_run_uses_its_tabs_provider_rather_than_rebuilding_the_params():
    """Source guard: the provider must be the *single* expression, not a copy.

    The registry gate above proves a provider exists; only this proves the run
    consults it. A callback that registered a provider and then rebuilt the same
    dict inline would pass every other test here and drift on the first edit —
    which is the failure mode the whole Step system exists to rule out.
    """
    import ast

    offenders = []
    for name, source in _tab_sources().items():
        module = ast.parse(source)
        providers = {
            key.value if isinstance(key, ast.Constant) else None
            for node in ast.walk(module) if isinstance(node, ast.Subscript)
            if "template_preview" in ast.unparse(node.value)
            for key in [node.slice]
        }
        if not providers:
            continue
        # Every function whose name ends in _preview is a provider; each must be
        # called somewhere other than at its own registration.
        for func in ast.walk(module):
            if not isinstance(func, ast.FunctionDef) or not func.name.endswith("_preview"):
                continue
            calls = [n for n in ast.walk(module)
                     if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
                     and n.func.id == func.name]
            if not calls:
                offenders.append(f"{name}.{func.name}")
    assert not offenders, (
        f"provider(s) {offenders} are registered but never called by their own "
        f"tab, so the run and the preview are two expressions of the settings"
    )


def test_every_step_resolves_user_overrides():
    """Source guard: a run must never assemble builtin text directly.

    ``genes.marker_plot`` did. It passed ``template=builtin_assemble(...)``,
    which by design cannot see an override path — so that template could be
    edited, validated and saved in this very tab with no effect, and its
    provenance nodes carried no ``template_id`` for the notebook banner or
    ``verify_notebook``'s ``stock_templates`` to notice. Every other call site
    splatted ``step_template``; nothing checked.
    """
    import ast

    #: Not run sites. The Templates tab renders a spec it has already resolved,
    #: and the ``_*_template`` helpers exist so tests can pin the *shipped* text
    #: without reading the developer's own overrides.
    allowed_builtin_callers = {"_preview"}

    offenders = []
    for name, source in _tab_sources().items():
        module = ast.parse(source)
        for node in ast.walk(module):
            if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                    and node.func.id == "Step"):
                continue
            explicit = [kw for kw in node.keywords if kw.arg == "template"]
            splatted = [kw for kw in node.keywords if kw.arg is None]
            enclosing = _enclosing_function(module, node)
            if enclosing in allowed_builtin_callers:
                continue
            if explicit:
                offenders.append(
                    f"{name}: Step(template=...) in {enclosing}() — use "
                    f"**step_template(id, blocks) so user overrides apply"
                )
            elif not any("step_template" in ast.unparse(kw.value)
                         or "_resolved" in ast.unparse(kw.value) for kw in splatted):
                offenders.append(
                    f"{name}: Step(...) in {enclosing}() gets its template from "
                    f"neither step_template nor _resolved"
                )
    assert not offenders, "\n".join(offenders)


def _enclosing_function(module, target) -> str:
    """Name of the innermost function containing *target*, or '<module>'."""
    import ast

    best = "<module>"
    for node in ast.walk(module):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if any(child is target for child in ast.walk(node)):
                best = node.name
    return best


# ── the providers, actually called ───────────────────────────────────────────

#: Tabs that register a provider, and can be built against the stub context
#: below. ``tab_cnv`` is here too — its provider covers inferCNV only, since
#: CopyKAT runs detached in another conda env and stays on ``record_node``.
_PROVIDER_TABS = (
    "tab_clustering", "tab_gene_analysis", "tab_nhood", "tab_co_occurrence",
    "tab_ligrec", "tab_gene_correlation", "tab_marker_genes", "tab_roi",
    "tab_cnv", "tab_umap", "tab_annot_nhood", "tab_annot_distance",
    "tab_transcripts", "tab_he_registration", "tab_publish", "tab_qc",
    "tab_preprocess",
)


@pytest.fixture
def stub_ctx():
    """Enough ``ViewerContext`` to build a tab and ask its provider.

    A stand-in rather than the real thing, which needs a napari viewer and a
    loaded dataset. It has to mirror the real contract closely, though —
    ``get_labels_for`` returns ``{}`` and not ``None`` because ``_helpers``
    guarantees a dict, and two tabs index it during ``build_tab``.
    """
    anndata = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")

    n = 40
    adata = anndata.AnnData(
        np.random.default_rng(0).poisson(3, (n, 8)).astype("float32"))
    adata.obs_names = [f"c{i}" for i in range(n)]
    adata.var_names = [f"Gene{i}" for i in range(8)]
    series = pd.Series(np.arange(n) % 3, index=adata.obs_names, name="leiden_r1.0")
    adata.obs["leiden_r1.0"] = pd.Categorical(series.astype(str))

    return SimpleNamespace(
        viewer=None, adata=adata, sdata=None, state={}, he_state={}, arms_state={},
        clusterings={"leiden_r1.0": series}, clustering_names=["leiden_r1.0"],
        gene_names=list(adata.var_names), data_path="/tmp/xv-not-written",
        pixel_size=0.2125, color_manager=SimpleNamespace(adata=adata),
        roi_layer=None, dataset_generation=0, no_cache=False,
        # The annotation tabs read the shapes layer while building their
        # previews. None is the honest stand-in: it is what ViewerContext holds
        # before anything has been drawn, and a provider must answer then too.
        annotation_layer=None, ensure_annotations=lambda preview: None,
        segmentation_source="xenium",
        gene_widget=SimpleNamespace(value="Gene0"),
        filter_check=SimpleNamespace(value=False),
        clustering_widget=SimpleNamespace(value="leiden_r1.0"),
        get_selected_cluster_ids=lambda: [],
        get_cluster_filter=lambda: None,
        get_labels_for=lambda key: {},
        refresh_clustering_choices=lambda *a, **k: None,
        record_clustering=lambda *a, **k: None,
        record_preamble=lambda *a, **k: None,
        record_node=lambda *a, **k: None, record_code=lambda *a, **k: None,
        set_status=lambda *a, **k: None, run_step=lambda step: {},
        ensure_normalized=lambda: adata, ensure_spatial_neighbors=lambda k: None,
        reload_dataset=None, external_images_state=[], patch_overlays_state=[],
        umap_viewer=None, plots_panel=None,
        # Plot providers name the files their step will write. Read-only in the
        # real helper too (``save_paths`` creates nothing), which is what lets a
        # provider call it — drawing a preview pane must not touch the disk.
        plot_paths=lambda stem: [f"/tmp/xv-not-written/plots/{stem}.png",
                                 f"/tmp/xv-not-written/plots/{stem}.pdf"],
        recorded_plot_paths=lambda paths: [str(p) for p in paths],
        show_plot=lambda *a, **k: [],
        apply_plot_font_size=lambda: None,
    )


@pytest.fixture
def live_providers(qapp, stub_ctx):
    """Every provider the real tabs register, with the widgets they read alive.

    The returned widgets are held for the lifetime of the fixture on purpose:
    a provider closes over magicgui widgets, and once Qt has collected the tab
    the C++ objects behind them are gone — every provider then raises, which
    ``_preview`` catches and turns into a silent fall back to sample values.
    That is precisely the failure this test exists to see.
    """
    import importlib

    held = []
    for name in _PROVIDER_TABS:
        module = importlib.import_module(f"palms.tabs.{name}")
        held.append(module.build_tab(stub_ctx))
    return stub_ctx.state.get("template_preview", {}), held


def test_every_registered_provider_answers(live_providers):
    """The gate above proves a provider is registered; this proves it works.

    ``_preview`` swallows a provider that raises and shows sample values
    instead — right for a half-built tab, and indistinguishable from a provider
    that is quietly broken. So call each one for real.
    """
    providers, _held = live_providers
    # Every template is accounted for, so this check covers the registry rather
    # than whichever tabs happened to build. Without it, a tab dropped from
    # _PROVIDER_TABS would take its provider out of the live check silently
    # while the source-level gate above still passed.
    uncovered = sorted(set(builtin_ids()) - set(providers) - set(_NO_PROVIDER))
    assert not uncovered, (
        f"template(s) {uncovered} are neither exercised here nor exempt — add "
        f"the owning tab to _PROVIDER_TABS"
    )
    failures = {}
    for template_id, provider in providers.items():
        try:
            provider()
        except Exception as exc:                        # noqa: BLE001 — reporting
            failures[template_id] = f"{type(exc).__name__}: {exc}"
    assert not failures, (
        f"provider(s) raised, so their preview silently degrades to sample "
        f"values: {failures}"
    )


def test_a_provider_selects_a_declared_assembly(live_providers):
    """Blocks a provider names must be ones the template declares.

    The same property ``test_template_registry`` asserts of the Python block
    selectors, now that the provider is what the call site actually uses.
    """
    providers, _held = live_providers
    for template_id, provider in providers.items():
        blocks = tuple(provider().blocks)
        assemblies = builtin_spec(template_id).assemblies
        assert blocks in assemblies, (
            f"{template_id}: provider selects {blocks}, which is not one of the "
            f"declared assemblies {assemblies}"
        )


def test_toggling_a_real_widget_changes_the_previews_shape(ctx, qapp, stub_ctx):
    """The end of the chain, driven through actual magicgui widgets.

    Every other test here hands ``_preview`` a stub provider. This one toggles
    the two checkboxes a user toggles and reads the rendered source, which is
    the only way to see the whole path — widget, provider, block selector,
    resolver, ``Step.render`` — agree. It is the property that was broken:
    the pane rendered ``assemblies[0]`` whatever the checkboxes said.
    """
    from qtpy.QtWidgets import QCheckBox

    from palms.tabs import tab_clustering
    from palms.tabs.tab_templates import _preview

    held, _exports = tab_clustering.build_tab(stub_ctx)  # held: Qt collects it
    spec = builtin_spec("clustering.leiden")
    ctx.state["template_preview"] = stub_ctx.state["template_preview"]

    # Found through *this* tab's widget tree, not by scanning live objects: an
    # earlier test in this file builds the same tab against another context, and
    # a global scan picks up whichever copy it happens to reach — the toggles
    # then land on a tab no provider is reading.
    boxes = {"Use HVGs only": None, "Scale (max_value=10)": None}
    for native in held.findChildren(QCheckBox):
        magic = getattr(native, "_magic_widget", None)
        if magic is not None and getattr(magic, "label", None) in boxes:
            boxes[magic.label] = magic
    assert all(boxes.values()), f"checkbox labels moved: {boxes}"

    for use_hvg in (False, True):
        for do_scale in (False, True):
            boxes["Use HVGs only"].value = use_hvg
            boxes["Scale (max_value=10)"].value = do_scale
            # assemblies[0] is passed deliberately: the provider's blocks must
            # override it, so a regression here shows up as a fixed shape.
            shown = _preview(ctx, spec, spec.assemblies[0])
            assert shown.startswith("# preview — current widget values")
            assert ("highly_variable_genes" in shown) is use_hvg
            assert ("sc.pp.scale" in shown) is do_scale
            # PCA is recomputed only when the gene set or the scaling changed.
            assert ("sc.pp.pca(adata_leiden" in shown) is (use_hvg or do_scale)
            compile(shown, "<preview>", "exec")


def test_a_provider_renders_the_template_it_selected(ctx, live_providers):
    """End to end: provider -> _preview -> a real rendered string.

    Not sample values, and not an error message — the two things the pane shows
    when something upstream has gone wrong.
    """
    from palms.tabs.tab_templates import _preview

    providers, _held = live_providers
    ctx.state["template_preview"] = providers
    for template_id, provider in providers.items():
        spec = builtin_spec(template_id)
        shown = _preview(ctx, spec, spec.assemblies[0])
        assert shown.startswith("# preview — current widget values"), (
            f"{template_id}: fell back to sample values\n{shown}"
        )
        # Parses as Python: the preview claims to be the string handed to exec.
        compile(shown, f"<preview:{template_id}>", "exec")


def test_the_contract_block_reports_the_declared_interface(tab):
    from palms.tabs.tab_templates import _contract_text

    spec = builtin_spec("roi.expression")
    text = _contract_text(spec)
    assert "roi_expr_cells" in text          # a declared output
    assert "pixel_size" in text              # a required param
    assert f"blocks    : {', '.join(spec.blocks)}" in text


def test_a_frozen_block_says_so(tab):
    """The Arrow shim is not editable; the view has to make that visible."""
    from palms.tabs.tab_templates import _contract_text

    text = _contract_text(builtin_spec("genes.cnv_infercnv"))
    assert "frozen    : arrow_shim" in text


def _select(widget, group_name="Clustering", index=0):
    tree = _tree(widget)
    group = next(tree.topLevelItem(i) for i in range(tree.topLevelItemCount())
                 if tree.topLevelItem(i).text(0) == group_name)
    tree.setCurrentItem(group.child(index))
    return group.child(index)


def test_selecting_a_template_fills_every_pane(tab):
    from qtpy.QtWidgets import QPlainTextEdit

    widget, _ = tab
    _select(widget)

    panes = widget.findChildren(QPlainTextEdit)
    assert len(panes) == 4, "contract, default, yours, preview"
    assert all(p.toPlainText().strip() for p in panes)


def test_exactly_one_pane_is_editable(tab):
    """The default must stay read-only, or there is no reference to diff against."""
    from qtpy.QtWidgets import QPlainTextEdit

    widget, _ = tab
    _select(widget)
    panes = widget.findChildren(QPlainTextEdit)
    writable = [p for p in panes if not p.isReadOnly()]
    assert len(writable) == 1, "only 'Yours' may be edited"


# ── editing ──────────────────────────────────────────────────────────────────

@pytest.fixture
def editable(qapp, tmp_path, monkeypatch):
    """A tab whose overrides land in tmp_path rather than the real config dir.

    One environment variable is enough because saving resolves its destination
    from the same search path reading uses — a write cannot land somewhere the
    reader does not look. Before that, these tests wrote into the developer's
    actual ``~/.config``.
    """
    from palms.utils.step_templates import (
        TEMPLATE_PATH_ENV, clear_cache, user_template_dir,
    )

    monkeypatch.setenv(TEMPLATE_PATH_ENV, str(tmp_path))
    clear_cache()
    assert user_template_dir() == tmp_path, "saves must be redirected too"

    from palms.tabs.tab_templates import build_tab
    widget, _ = build_tab(SimpleNamespace(state={}))
    yield widget, tmp_path
    clear_cache()


def _panes(widget):
    from qtpy.QtWidgets import QPlainTextEdit
    return widget.findChildren(QPlainTextEdit)


def _editor(widget):
    return next(p for p in _panes(widget) if not p.isReadOnly())


def _button(widget, text):
    from qtpy.QtWidgets import QPushButton
    return next(b for b in widget.findChildren(QPushButton)
                if b.text().replace("&", "") == text)


def _problems(widget):
    from qtpy.QtWidgets import QListWidget
    lst = widget.findChild(QListWidget)
    return [lst.item(i).text() for i in range(lst.count())]


@pytest.mark.parametrize("template_id", builtin_ids())
def test_the_editor_round_trips_every_template_body(template_id):
    """Parsing an edit back must agree with how the loader splits a file.

    Every template, because the bug this caught was specific to one: an
    annotation on the marker line of a frozen block came back as a *different*
    block name, so editing the CNV Arrow shim silently created a new block
    instead of modifying the protected one.
    """
    from palms.tabs.tab_templates import blocks_to_text, text_to_blocks

    spec = builtin_spec(template_id)
    assert text_to_blocks(blocks_to_text(spec)) == {
        name: block.text for name, block in spec.blocks.items()
    }


def test_validate_reports_a_broken_edit_without_saving(editable):
    widget, tmp_path = editable
    _select(widget)
    _editor(widget).setPlainText("#--- block head\nadata_leiden = (")

    _button(widget, "Validate").click()
    assert any("not valid Python" in p for p in _problems(widget))
    assert list(tmp_path.glob("*.tmpl")) == [], "Validate must not write anything"


def test_saving_a_valid_edit_activates_it(editable):
    from palms.utils.step_templates import resolve

    widget, tmp_path = editable
    _select(widget)
    spec = builtin_spec("clustering.leiden")
    from palms.tabs.tab_templates import blocks_to_text
    edited = blocks_to_text(spec).replace("max_value=10", "max_value=99")
    _editor(widget).setPlainText(edited)

    _button(widget, "Save  Activate").click()

    assert (tmp_path / "clustering.leiden.tmpl").is_file()
    resolved = resolve("clustering.leiden")
    assert resolved.is_customised
    assert resolved.changed_blocks() == ["scale"], (
        "only the block that actually differs should be written"
    )


def test_saving_a_broken_edit_keeps_the_text_but_does_not_activate(editable):
    """Never refuse the write; gate activation instead."""
    from palms.utils.step_templates import resolve

    widget, tmp_path = editable
    _select(widget)
    from palms.tabs.tab_templates import blocks_to_text
    spec = builtin_spec("clustering.leiden")
    _editor(widget).setPlainText(
        blocks_to_text(spec).replace("sc.pp.scale(adata_leiden, max_value=10)",
                                     "sc.pp.scale(adata_leiden,"))

    _button(widget, "Save  Activate").click()

    assert (tmp_path / "clustering.leiden.tmpl").is_file(), "the edit is kept"
    resolved = resolve("clustering.leiden")
    assert resolved.rejected
    assert not resolved.is_customised, "the shipped template runs instead"


def test_reverting_removes_the_override(editable):
    from palms.utils.step_templates import resolve

    widget, tmp_path = editable
    _select(widget)
    from palms.tabs.tab_templates import blocks_to_text
    spec = builtin_spec("clustering.leiden")
    _editor(widget).setPlainText(
        blocks_to_text(spec).replace("max_value=10", "max_value=99"))
    _button(widget, "Save  Activate").click()
    assert resolve("clustering.leiden").is_customised

    _button(widget, "Revert to default").click()
    assert list(tmp_path.glob("*.tmpl")) == []
    assert not resolve("clustering.leiden").is_customised


def test_saving_an_unchanged_edit_removes_rather_than_writes_an_inert_file(editable):
    widget, tmp_path = editable
    _select(widget)
    _button(widget, "Save  Activate").click()
    assert list(tmp_path.glob("*.tmpl")) == [], (
        "a file that says nothing would show as 'customised' forever"
    )


def test_a_frozen_block_edit_is_reported(editable):
    widget, _ = editable
    _select(widget, "Genes")
    from qtpy.QtCore import Qt
    tree = _tree(widget)
    cnv = next(
        tree.topLevelItem(g).child(c)
        for g in range(tree.topLevelItemCount())
        for c in range(tree.topLevelItem(g).childCount())
        if tree.topLevelItem(g).child(c).data(0, Qt.ItemDataRole.UserRole) == "genes.cnv_infercnv"
    )
    tree.setCurrentItem(cnv)

    editor = _editor(widget)
    editor.setPlainText(editor.toPlainText().replace(
        "_old_infer = pd.options.future.infer_string", "_old_infer = False"))
    _button(widget, "Validate").click()
    assert any("workaround" in p for p in _problems(widget))
