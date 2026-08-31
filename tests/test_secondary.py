import numpy as np
import pytest

from regfans import secondary
from regfans.fan import Fan
from regfans.grow4d import grow4d
from regfans.vectorconfig import VectorConfiguration

# the 8 vertices of a cube: 64 fine triangulations, 46 of them regular
_VECS = np.ascontiguousarray(
    [[1,1,1], [1,1,-1], [1,-1,1], [1,-1,-1],
     [-1,1,1], [-1,1,-1], [-1,-1,1], [-1,-1,-1]], dtype=np.int32)


def _enumerate(vecs=_VECS):
    vc = VectorConfiguration(vecs)
    simps, starts, num_fans, status, _ = grow4d(vecs)
    assert status == 0
    return vc, simps, starts, num_fans


def test_walls_match_interior_facets():
    """walls() must find exactly the facets shared by two simplices."""
    vc, simps, starts, num_fans = _enumerate()
    labels = np.asarray(vc.labels)
    masks, apexes, wall_starts = secondary.walls(simps, starts, vc.ambient_dim)

    for i in range(num_fans):
        fan = Fan(vc, labels[simps[starts[i]:starts[i+1]]].tolist())
        expected = {frozenset(f) | {a for c in cs for a in c if a not in f}
                    for f, cs in fan.facets().items() if len(cs) == 2}

        got = set()
        for m in masks[wall_starts[i]:wall_starts[i+1]]:
            got.add(frozenset(labels[j] for j in range(len(labels))
                              if (int(m) >> j) & 1))
        assert got == expected

        # the apexes are always a sub-pair of their wall
        for m, a in zip(masks[wall_starts[i]:wall_starts[i+1]],
                        apexes[wall_starts[i]:wall_starts[i+1]]):
            assert int(a) & int(m) == int(a)
            assert bin(int(a)).count("1") == 2


def test_regular_mask_matches_per_fan():
    """The bulk pass must agree with Fan.is_regular, fan by fan."""
    vc, simps, starts, num_fans = _enumerate()
    labels = np.asarray(vc.labels)

    reference = np.array([
        Fan(vc, labels[simps[starts[i]:starts[i+1]]].tolist()).is_regular()
        for i in range(num_fans)])

    assert np.array_equal(secondary.regular_mask(VectorConfiguration(_VECS),
                                                 simps, starts), reference)
    assert reference.any() and not reference.all()   # a real mix, not trivial


def test_regular_mask_empty():
    vc = VectorConfiguration(_VECS)
    out = secondary.regular_mask(vc, np.zeros((0, 3), np.uint32),
                                 np.zeros(1, np.int32))
    assert out.shape == (0,)


def test_backend_agrees_with_per_fan_regularity():
    """all_triangulations must return the same fans through either path."""
    vc = VectorConfiguration(_VECS)

    bulk = vc.all_triangulations(only_fine=True, only_regular=True,
                                 backend="grow4d")
    # only_fine=False takes the per-fan path, so intersecting with the fine
    # ones reproduces the same set by a different route
    per_fan = [f for f in vc.all_triangulations(only_fine=False,
                                                only_regular=True,
                                                backend="grow4d")
               if f.is_fine()]

    assert {frozenset(f.cones()) for f in bulk} == \
           {frozenset(f.cones()) for f in per_fan}


def test_prescreened_fans_report_regular():
    """The bulk verdict is recorded, not silently dropped."""
    vc = VectorConfiguration(_VECS)
    for f in vc.all_triangulations(only_fine=True, only_regular=True,
                                   backend="grow4d"):
        assert f._is_regular is True
        assert f.is_regular()


def test_hyperplane_orientation_follows_the_apexes():
    """
    The same labels can span a wall in two fans with the apexes on opposite
    sides of the dependency, which is the same inequality reversed. So the
    hyperplane cannot be keyed on the wall alone -- doing so silently builds
    the wrong secondary cone for one of the two.

    Apex pairs drawn from the same side agree exactly; pairs from opposite
    sides differ by a sign, and nothing else is possible, since the
    dependency is unique up to scale.
    """
    vc, simps, starts, num_fans = _enumerate()
    masks, apexes, _ = secondary.walls(simps, starts, vc.ambient_dim)

    rows, row_of = secondary._hyperplanes(vc, masks, apexes)

    by_wall = {}
    for mask, apex in {(int(m), int(a)) for m, a in zip(masks, apexes)}:
        by_wall.setdefault(mask, []).append(rows[row_of[(mask, apex)]])

    reversed_somewhere = False
    for mask, got in by_wall.items():
        for other in got[1:]:
            same = np.array_equal(got[0], other)
            flipped = np.array_equal(got[0], -other)
            assert same or flipped, f"wall {mask} gave an unrelated row"
            reversed_somewhere |= flipped

    assert reversed_somewhere, (
        "no wall appeared with reversed orientation, so this test would "
        "not catch keying the hyperplane on the wall alone")
