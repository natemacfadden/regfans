from regfans.circuits import Circuits
from regfans.vectorconfig import VectorConfiguration


# Shared test input: 7 vectors in 4D, totally cyclic
_PTS = [
    [-5, -2, -3, -2], [0, 1, 0, 0], [0, 1, 2, 0],
    [1, 0, 0, 0], [1, 0, 1, 2], [-2, -1, -1, 0], [0, 1, 1, 0],
]


# Circuits.encode / decode
# ------------------------
def test_encode_single():
    c = Circuits()
    assert c.encode([0]) == 1     # bit 0 set
    assert c.encode([1]) == 2     # bit 1 set
    assert c.encode([3]) == 8     # bit 3 set


def test_encode_multiple():
    c = Circuits()
    assert c.encode([0, 1, 2]) == 7   # 0b0111
    assert c.encode([0, 2]) == 5      # 0b0101


def test_decode():
    c = Circuits()
    assert c.decode(7) == [0, 1, 2]
    assert c.decode(5) == [0, 2]
    assert c.decode(8) == [3]


def test_encode_decode_roundtrip():
    c = Circuits()
    for indices in [[0, 1], [1, 3, 5], [0, 2, 4, 6], [0]]:
        assert c.decode(c.encode(indices)) == indices


def test_encode_int_passthrough():
    # passing an already-encoded int is a no-op
    c = Circuits()
    assert c.encode(7) == 7


# Circuits.is_subset / is_superset
# ---------------------------------
def test_is_subset():
    c = Circuits()
    A = c.encode([0, 1])
    B = c.encode([0, 1, 2])
    assert c.is_subset(A, B)        # {0,1} subset of {0,1,2}
    assert not c.is_subset(B, A)    # {0,1,2} not subset of {0,1}


def test_is_superset():
    c = Circuits()
    A = c.encode([0, 1])
    B = c.encode([0, 1, 2])
    assert c.is_superset(B, A)      # {0,1,2} superset of {0,1}
    assert not c.is_superset(A, B)  # {0,1} not superset of {0,1,2}


def test_subset_superset_equal():
    c = Circuits()
    A = c.encode([0, 1])
    assert c.is_subset(A, A)
    assert c.is_superset(A, A)


def test_subset_superset_disjoint():
    c = Circuits()
    A = c.encode([0, 1])
    B = c.encode([2, 3])
    assert not c.is_subset(A, B)
    assert not c.is_superset(A, B)


# Circuit.reorient
# ----------------
def test_reorient_swaps_pos_neg():
    vc = VectorConfiguration(_PTS)
    for circ in vc.circuits():
        reoriented = circ.reorient()
        assert reoriented.Zpos == circ.Zneg
        assert reoriented.Zneg == circ.Zpos


def test_reorient_negates_lmbda():
    vc = VectorConfiguration(_PTS)
    for circ in vc.circuits():
        reoriented = circ.reorient()
        assert reoriented.lmbda == tuple(-x for x in circ.lmbda)


def test_reorient_swaps_signature():
    vc = VectorConfiguration(_PTS)
    for circ in vc.circuits():
        reoriented = circ.reorient()
        assert reoriented.signature == (circ.signature[1], circ.signature[0])


def test_reorient_preserves_support():
    vc = VectorConfiguration(_PTS)
    for circ in vc.circuits():
        reoriented = circ.reorient()
        assert reoriented.Z == circ.Z


def test_reorient_double_is_identity():
    vc = VectorConfiguration(_PTS)
    for circ in vc.circuits():
        assert circ.reorient().reorient().Zpos == circ.Zpos
        assert circ.reorient().reorient().lmbda == circ.lmbda
