"""Support for driving the running application from outside it.

Not part of the analysis API. It exists because two callers need to click the
same buttons — ``scripts/capture_screenshots.py`` and ``tests/e2e/`` — and a
private copy in either would drift from the other the moment one of them learned
something new about a widget.
"""
from palms.testing.rig import Rig, process_events

__all__ = ["Rig", "process_events"]
