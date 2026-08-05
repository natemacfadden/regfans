import numpy as np
from regfans import util


# util.gcd
# --------
def test_gcd_integers():
    assert util.gcd([2, 4, 6]) == 2
    assert util.gcd([1, 2, 3]) == 1
    assert util.gcd([3, 3, 3]) == 3
    assert util.gcd([6]) == 6

def test_gcd_rationals():
    assert util.gcd([2.0, 3.0, 1.5]) == 0.5
    assert util.gcd([0.5, 1.0, 1.5]) == 0.5
    assert util.gcd([1.0, 1.5]) == 0.5

# util.primitive
# --------------
def test_primitive_integers():
    assert util.primitive([2, 4, 6]) == [1, 2, 3]
    assert util.primitive([3, 6]) == [1, 2]
    assert util.primitive([1, 2, 3]) == [1, 2, 3]
    assert util.primitive([6, 4, 2]) == [3, 2, 1]
    assert util.primitive([2, 2]) == [1, 1]

def test_primitive_rationals():
    assert util.primitive([0.5, 1.0, 1.5]) == [1, 2, 3]

# util.lerp
# ---------
def test_lerp_endpoints():
    p0 = [0.0, 0.0]
    p1 = [2.0, 4.0]
    np.testing.assert_allclose(util.lerp(p0, p1, 0), [0.0, 0.0])
    np.testing.assert_allclose(util.lerp(p0, p1, 1), [2.0, 4.0])

def test_lerp_midpoint():
    p0 = [0.0, 0.0]
    p1 = [2.0, 4.0]
    np.testing.assert_allclose(util.lerp(p0, p1, 0.5), [1.0, 2.0])

def test_lerp_quarter():
    p0 = [1.0, 2.0]
    p1 = [5.0, 6.0]
    np.testing.assert_allclose(util.lerp(p0, p1, 0.25), [2.0, 3.0])

# util.contains
# -------------
def test_contains_positive_orthant():
    # positive orthant: H @ p >= 0 <=> p[0] >= 0 and p[1] >= 0
    H = [[1, 0], [0, 1]]
    assert util.contains(p=[1, 1], H=H)
    assert util.contains(p=[0, 0], H=H)    # origin (on boundary)
    assert util.contains(p=[5, 0], H=H)    # on a wall
    assert not util.contains(p=[-1, 1], H=H)
    assert not util.contains(p=[1, -1], H=H)
    assert not util.contains(p=[-1, -1], H=H)

def test_contains_general_cone():
    # cone spanned by (1,0) and (0,1) in 2D - same as positive orthant
    # but specified via rays instead
    R = [[1, 0], [0, 1]]
    assert util.contains(p=[1, 1], R=R)
    assert not util.contains(p=[-1, 1], R=R)
