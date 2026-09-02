"""The convention arithmetic in scripts/compare_he_registration.py.

Getting a convention backwards is the entire risk of that tool: it would report
a large disagreement between two transforms that actually agree, or — worse —
a small one between two that do not, and either would be read as a statement
about the registration rather than about the script.

So the two conventions are tested as *properties*, against transforms whose
answer is known by construction, not against remembered numbers from a dataset
that is not in the repo.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

np = pytest.importorskip("numpy")

import compare_he_registration as chr  # noqa: E402


def _similarity_xy(scale, deg, tx, ty):
    """A similarity in (x, y), the convention 10x's CSV is written in."""
    t = np.radians(deg)
    return np.array([[scale * np.cos(t), -scale * np.sin(t), tx],
                     [scale * np.sin(t), scale * np.cos(t), ty],
                     [0, 0, 1]], dtype=float)


def test_the_swap_is_the_same_map_read_in_the_other_order():
    """M_yx = P M_xy P, stated as what it has to *do* rather than as a formula.

    Transforming a point in (y, x) with the swapped matrix must equal swapping
    the point, transforming in (x, y), and swapping the answer back.
    """
    m_xy = _similarity_xy(1.3, 89.9, -847.0, 15680.6)
    m_yx = chr.to_yx(m_xy)

    pts_yx = np.array([[0.0, 0.0], [27502.0, 14896.0], [123.4, 567.8]])
    got = chr.apply(m_yx, pts_yx)

    pts_xy = pts_yx[:, ::-1]
    homo = np.hstack([pts_xy, np.ones((len(pts_xy), 1))])
    expected = (m_xy @ homo.T).T[:, :2][:, ::-1]

    np.testing.assert_allclose(got, expected, atol=1e-9)


def test_the_swap_is_its_own_inverse():
    m_xy = _similarity_xy(0.7, -33.0, 5.0, -9.0)
    np.testing.assert_allclose(chr.to_yx(chr.to_yx(m_xy)), m_xy, atol=1e-12)


@pytest.mark.parametrize("scale,deg", [(1.0, 0.0), (1.28889, 89.888), (0.5, -120.0)])
def test_decompose_recovers_what_was_put_in(scale, deg):
    """Read in (y, x), a similarity built in (x, y) keeps its scale and angle."""
    d = chr.decompose(chr.to_yx(_similarity_xy(scale, deg, 11.0, -22.0)))
    assert d["scale"] == pytest.approx(scale, rel=1e-9)
    assert d["rotation_deg"] == pytest.approx(deg, abs=1e-6)


def _write_inputs(tmp_path, m_xy, jitter_um=0.0, pixel_size=0.2125):
    """Landmarks generated *from* a transform, so the true answer is known."""
    he = np.array([[100.0, 200.0], [9000.0, 1200.0], [4000.0, 13000.0]])
    xen = chr.apply(chr.to_yx(m_xy), he)
    if jitter_um:
        xen = xen + jitter_um / pixel_size      # a uniform shift, in xenium px
    lm = tmp_path / "landmarks.json"
    lm.write_text(json.dumps({
        "xenium_landmarks_yx": xen.tolist(),
        "he_landmarks_yx": he.tolist(),
        "affine_3x3_yx": chr.to_yx(m_xy).tolist(),
    }))
    csv = tmp_path / "align.csv"
    np.savetxt(csv, m_xy, delimiter=",")
    return lm, csv


def _run(tmp_path, lm, csv):
    out = tmp_path / "report.json"
    sys.argv = ["compare_he_registration.py", str(lm), str(csv),
                "--he-shape", "27502", "14896", "--grid", "8", "--out", str(out)]
    chr.main()
    return json.loads(out.read_text())


def test_two_identical_transforms_disagree_by_nothing(tmp_path, capsys):
    m_xy = _similarity_xy(1.28889, 89.888, -847.0, 15680.6)
    report = _run(tmp_path, *_write_inputs(tmp_path, m_xy))
    assert report["disagreement_with_10x_um"]["max"] == pytest.approx(0.0, abs=1e-6)
    assert report["palms"]["scale"] == pytest.approx(report["tenx"]["scale"], rel=1e-12)


def test_a_known_offset_comes_back_as_that_offset(tmp_path, capsys):
    """A 5 um shift in the fitted transform must read as 5 um of disagreement.

    This is what makes the reported number a length rather than an arbitrary
    score — and it is where a wrong pixel_size or a missing unit conversion
    would show up.
    """
    m_xy = _similarity_xy(1.28889, 89.888, -847.0, 15680.6)
    lm, csv = _write_inputs(tmp_path, m_xy)

    shifted = json.loads(lm.read_text())
    affine = np.array(shifted["affine_3x3_yx"])
    affine[:2, 2] += 5.0 / 0.2125            # 5 um, in morphology pixels
    shifted["affine_3x3_yx"] = affine.tolist()
    lm.write_text(json.dumps(shifted))

    report = _run(tmp_path, lm, csv)
    d = report["disagreement_with_10x_um"]
    assert d["mean"] == pytest.approx(5.0 * np.sqrt(2), rel=1e-6)


def test_a_landmarks_file_without_a_registration_says_so(tmp_path):
    lm = tmp_path / "landmarks.json"
    lm.write_text(json.dumps({"xenium_landmarks_yx": [], "he_landmarks_yx": []}))
    csv = tmp_path / "align.csv"
    np.savetxt(csv, np.eye(3), delimiter=",")
    sys.argv = ["compare_he_registration.py", str(lm), str(csv),
                "--he-shape", "10", "10"]
    with pytest.raises(SystemExit, match="Compute Registration"):
        chr.main()
