"""Per-tab UI modules for the PALMS control panel.

Import the tab module you want (``from palms.tabs import tab_cache``) and call its
``build_tab(ctx)``. This package deliberately re-exports nothing: ``app.py`` imports each
``build_tab`` *inside* the function that builds the panel, so the napari-heavy tab modules
are not loaded until a viewer is actually being constructed. Eager re-exports here would
undo that for anything that touches the package — and the ones that used to live here
covered 11 of the 26 tabs, chosen by nothing but the order they were written in.
"""
