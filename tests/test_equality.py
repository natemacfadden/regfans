# Regression tests for VectorConfiguration / Fan equality and hashing.
# These cover two equality bugs: the inverted VectorConfiguration.__ne__
# type-mismatch branch (a non-VC comparand was reported as equal), and
# Fan.__eq__ raising AttributeError on a non-Fan comparand.

import pytest
from regfans.vectorconfig import VectorConfiguration
from regfans.fan import Fan

# 7 vectors in 4D (a totally-cyclic reflexive configuration)
_PTS = [
    [-5, -2, -3, -2], [0, 1, 0, 0], [0, 1, 2, 0],
    [1, 0, 0, 0], [1, 0, 1, 2], [-2, -1, -1, 0], [0, 1, 1, 0],
]


def test_vc_eq_self_and_copy():
    vc = VectorConfiguration(_PTS)
    assert vc == vc
    assert vc == vc.copy()
    assert not (vc != vc.copy())


def test_vc_eq_non_vc():
    # a VectorConfiguration is never equal to a non-VectorConfiguration
    # (regression: the __ne__ type-mismatch branch used to be inverted)
    vc = VectorConfiguration(_PTS)
    assert (vc == 5) is False
    assert (vc != 5) is True
    assert (vc == "foo") is False
    assert (vc != "foo") is True


def test_vc_hash():
    # equal configurations must hash equally and work as set/dict keys
    vc1 = VectorConfiguration(_PTS)
    vc2 = VectorConfiguration(_PTS)
    assert vc1 == vc2
    assert hash(vc1) == hash(vc2)
    assert len({vc1, vc2}) == 1


def test_fan_eq_non_fan():
    # comparing a Fan to a non-Fan returns False rather than raising
    # (regression: Fan.__eq__ used to access o.vc and raise AttributeError)
    fan = VectorConfiguration(_PTS).triangulate()
    assert (fan == 5) is False
    assert (fan != 5) is True
    assert (fan == "foo") is False


def test_fan_eq_and_hash():
    # two triangulations of the same configuration compare and hash equal
    fan1 = VectorConfiguration(_PTS).triangulate()
    fan2 = VectorConfiguration(_PTS).triangulate()
    assert fan1 == fan2
    assert hash(fan1) == hash(fan2)
    assert len({fan1, fan2}) == 1
