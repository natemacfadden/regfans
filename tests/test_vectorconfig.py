import numpy as np
import pytest

from regfans import util
from regfans.vectorconfig import VectorConfiguration

# non-totally-cyclic VC tests
# ---------------------------
def test_acyclic():
    pts = [[0,0,1], [0,1,1], [1,0,1], [1,1,1]]
    vc  = VectorConfiguration(pts)

    assert vc.is_acyclic()
    assert not vc.is_totally_cyclic()
    assert sorted([tuple(n) for n in vc.support()]) == [(-1, 0, 1), (0, -1, 1), (0, 1, 0), (1, 0, 0)]

def test_lineality():
    pts = [[-1,0], [0,1], [1,0]]
    vc  = VectorConfiguration(pts)

    assert not vc.is_acyclic()
    assert not vc.is_totally_cyclic()
    assert sorted([tuple(n) for n in vc.support()]) == [(0,1)]

# 4D reflexive polytope tests
# ---------------------------
def test_eq():
    pts1 = [[-5, -2, -3, -2], [0, 1, 0, 0], [0, 1, 2, 0], [1, 0, 0, 0], [1, 0, 1, 2], [-2, -1, -1, 0], [0, 1, 1, 0]]
    pts2 = [[-2, -1, -2, -3], [0, 0, 1, 0], [0, 0, 1, 3], [0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 1], [0, 0, 1, 2]]

    vc1  = VectorConfiguration(pts1)
    vc2  = VectorConfiguration(pts2)

    assert not (vc1 == vc2)
    assert vc1 != vc2
    assert vc1 == vc1.copy()

def test_basic_and_triangulate():
    all_pts = [
         [[-5, -2, -3, -2], [0, 1, 0, 0], [0, 1, 2, 0], [1, 0, 0, 0], [1, 0, 1, 2], [-2, -1, -1, 0], [0, 1, 1, 0]] ,
         [[-2, -1, -2, -3], [0, 0, 1, 0], [0, 0, 1, 3], [0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 1], [0, 0, 1, 2]] ,
         [[0, 0, 0, 1], [1, -3, -2, -1], [-2, 2, 1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [-1, 1, 1, 0]] ,
         [[-2, -1, -1, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [2, -1, 0, -1], [-1, 0, 0, 0], [1, 0, 0, 0]] ,
         [[-2, -2, -1, -2], [-1, 1, -1, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [-1, -1, 0, -1]] ,
         [[0, 0, 1, 0], [0, 0, 0, 1], [0, 1, 0, 0], [1, -2, -1, -1], [-2, 0, -1, 0], [1, 0, 0, 0], [-1, 0, 0, 0]] ,
         [[-1, -1, 1, 0], [0, 0, 0, 1], [0, 0, 1, 0], [-1, 0, -1, -1], [0, -1, -1, -1], [0, 1, 0, 0], [1, 0, 0, 0]] ,
         [[0, 0, 1, 0], [1, 0, 0, 0], [-1, -1, 0, 2], [0, -1, -1, 0], [0, 1, 0, 0], [-1, 0, -1, -1], [0, 0, 0, 1]] ,
         [[-2, -1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [2, 0, -1, -1], [-1, 0, 0, 0]] ,
         [[1, 0, 0, 0], [-5, -3, -1, -2], [0, 1, 0, 0], [1, 1, 1, 2], [-2, 0, -1, 0], [0, 0, 1, 0], [-2, -1, 0, 0]] ,
         [[1, 0, 0, 0], [-3, -1, -1, -2], [-2, -1, 0, 0], [-2, 0, -1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [1, 1, 1, 2]] ,
         [[-2, -1, -2, -2], [0, 0, 1, 0], [0, 0, 1, 2], [0, 1, 0, 0], [1, 0, 0, 0], [-1, 0, -1, -1], [0, 0, 1, 1]] ,
         [[1, 0, 0, 0], [-1, 1, -1, -1], [0, 0, 0, 1], [0, 0, 1, 0], [-2, -2, -1, 0], [0, 1, 0, 0], [-1, -1, 0, 0]] ,
         [[0, -2, -2, -1], [0, 0, 0, 1], [-1, 1, 1, 0], [0, -1, -1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]] ,
         [[1, 0, 0, 0], [-2, 2, -3, -2], [0, 0, 0, 1], [0, 0, 1, 0], [-1, -1, 0, 0], [0, 1, 0, 0], [-1, 1, -1, -1]] ,
         [[-1, -1, 1, -1], [0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0], [-1, 0, -1, 0], [0, -1, -1, 0], [0, 0, 1, 0]] ,
         [[-1, -1, -1, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [2, -1, 0, -1], [-1, 0, 0, 0], [1, 0, 0, 0]] ,
         [[-2, 0, 0, -1], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, -1, -1, 0], [1, 0, 0, 0], [-1, 0, 0, 0]] ,
         [[0, 1, 0, 0], [-1, -2, -2, -1], [0, 0, 0, 1], [0, 0, 1, 0], [0, -1, -1, 0], [1, 0, 0, 0], [-1, 0, 1, 0]] ,
         [[0, 0, 1, 0], [1, 0, 0, 0], [-1, 0, -1, 2], [0, 0, 0, 1], [0, 1, 0, 0], [-2, -2, -1, -1], [-1, -1, 0, 0]] ,
         [[-2, -1, 2, 0], [-1, -1, -1, -1], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [-1, 0, 1, 0]] ,
         [[-1, -1, -1, -1], [0, 0, 0, 1], [0, 0, 1, 0], [-1, 0, 1, 0], [0, -1, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]] ,
         [[0, 0, 0, 1], [0, 0, 1, 0], [1, -2, -1, -1], [-1, 1, 0, 0], [0, 1, 0, 0], [-1, 0, 0, 0], [1, 0, 0, 0]] ,
         [[-3, -3, -3, -2], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [-2, -2, -2, -1], [-1, -1, -1, 0]] ,
         [[0, 1, 0, 0], [1, 0, 0, 0], [-2, -2, -2, -1], [-1, -1, 0, -1], [0, 0, 0, 1], [0, 0, 1, 0], [-1, -1, -1, 0]] ,
         [[0, 1, 0, 0], [0, -1, -1, -1], [0, 0, 0, 1], [0, 0, 1, 0], [-2, -1, 0, 0], [1, 0, 0, 0], [-1, 0, 0, 0]] ,
         [[-2, -1, 1, 1], [0, -1, -1, -1], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [-1, -1, 0, 0]] ,
         [[0, -1, -1, -1], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [-2, 1, 0, 0], [1, 0, 0, 0], [-1, 1, 0, 0]] ,
         [[-2, 0, -1, 0], [-1, -2, -1, -1], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [-1, 0, 0, 0]] ,
         [[-2, -2, 0, -1], [0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0], [-1, -1, 0, 0], [0, 0, 1, 0], [1, 0, -1, 0]] ,
         [[0, 1, 0, 0], [-2, -1, 2, -2], [0, 0, 0, 1], [1, 0, 0, 0], [-1, -1, -1, 0], [-1, 0, 1, -1], [0, 0, 1, 0]] ,
         [[1, 0, 0, 0], [-1, -1, -1, 0], [-1, -1, 0, -1], [-1, 0, -1, -1], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]] ,
         [[0, 0, 1, 0], [0, 1, 0, 0], [0, -1, -1, -1], [0, 0, 0, 1], [-1, -1, -1, -1], [-1, 0, 0, 1], [1, 0, 0, 0]] ,
         [[-1, -1, -1, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, -1, 0, -1], [-1, 0, 0, 0], [1, 0, 0, 0]] ,
         [[1, 0, 0, 0], [0, -1, -1, -1], [0, 0, 0, 1], [-1, -1, 0, 0], [-1, 0, -1, 0], [0, 0, 1, 0], [0, 1, 0, 0]] ,
         [[1, 0, 0, 0], [-1, -1, 0, 0], [-1, 0, -1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1], [1, 0, 0, -1]] ,
         [[0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0], [-1, -1, -1, 0], [-1, 0, 1, -1], [0, -1, 1, -1], [0, 0, 1, 0]] ,
         [[-1, 0, 1, -1], [0, -1, -1, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [-1, 1, 0, 0], [1, 0, 0, 0]] ,
         [[-1, -1, -1, 0], [0, 0, 1, 0], [-1, 2, 0, -1], [0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0], [0, 1, 0, -1]] ,
         [[0, 0, 0, 1], [0, 0, 1, 0], [-1, -2, -1, -1], [0, -1, -1, -1], [0, 1, 0, 0], [1, 0, 0, 0], [-1, 0, 0, 0]] ,
         [[0, 0, 0, 1], [0, 0, 1, 0], [1, 1, -2, -1], [-1, 0, 0, 0], [0, -1, 0, 0], [0, 1, 0, 0], [1, 0, 0, 0]] ,
         [[0, 0, 0, 1], [0, 0, 1, 0], [1, 0, -1, -1], [-1, -1, 0, 0], [-1, 0, 0, 0], [0, 1, 0, 0], [1, 0, 0, 0]] ,
         [[-1, 0, 0, -1], [0, -1, -1, 0], [0, 0, 0, 1], [0, 0, 1, 0], [1, 1, 0, 0], [0, 1, 0, 0], [1, 0, 0, 0]] ,
         [[0, 0, 0, 1], [0, 0, 1, 0], [2, 0, -1, -1], [-1, -1, 0, 0], [-1, 0, 0, 0], [0, 1, 0, 0], [1, 0, 0, 0]] ,
         [[-1, -1, 0, 0], [0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, -1, 0], [1, 0, 0, -1]] ,
         [[0, 1, 0, 0], [1, 0, 0, 0], [-1, 0, -1, 0], [0, -1, 0, -1], [0, 0, 0, 1], [0, 0, 1, 0], [-1, -1, 0, 0]] ,
         [[-1, 0, 0, -1], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [0, -1, -1, 0], [1, 1, 1, 0]] ,
         [[-2, -1, -1, -1], [0, -1, -1, 1], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [-1, -1, -1, 0]] ,
         [[0, 1, 0, 0], [0, 1, 1, 1], [1, 0, 0, 0], [-1, -1, -1, 0], [-1, -1, 0, -1], [0, 0, 0, 1], [0, 0, 1, 0]] ,
         [[0, 0, 1, 0], [0, 1, 0, 0], [0, -1, -1, -1], [0, 0, 0, 1], [1, -1, -1, 0], [1, 0, 0, 0], [-1, 0, 0, 0]] ,
         [[0, 0, 0, 1], [0, 0, 1, 0], [1, -1, -1, -1], [-1, 0, 0, 0], [0, -1, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]] ,
         [[-1, -1, -1, -1], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [-1, 0, 0, 0], [1, -1, 0, 0]] ,
         [[-1, -1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [1, 0, 0, 0], [1, 1, -1, -1], [0, -1, 1, 0], [0, 1, 0, 0]] ,
         [[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [2, -1, -1, -1], [-1, 0, 0, 0], [1, -2, -1, -1], [1, 0, 0, 0]] ,
         [[-5, -3, -3, -3], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [-3, -2, -2, -2], [-1, -1, -1, -1]] ,
         [[-9, -5, -3, -2], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [1, 1, 1, 2], [-4, -2, -1, 0], [-2, -1, 0, 0]] ,
         [[-1, -1, -1, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [-1, -1, -1, -1], [0, 0, 0, 1], [-1, -1, -1, 0]] ,
         [[-1, -1, -1, -1], [-1, -1, 1, 1], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [-1, -1, 0, 0]] ,
         [[-2, -2, -2, -1], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [-1, -1, -1, 0], [-1, 0, 0, 1], [1, 0, 0, 0]] ,
         [[-3, -2, -2, -2], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [-2, -1, -1, -1], [1, 0, 0, 0], [-1, -1, -1, -1]] ,
         [[0, -2, -2, -1], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 0], [0, -1, -1, 0], [1, 0, 0, 0]] ,
         [[-2, 0, 0, -1], [0, -1, -1, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [-1, 0, 0, 0]] ,
         [[0, 1, 0, 0], [1, 0, 0, 0], [-5, -3, -1, -2], [-4, -2, -1, 0], [0, 0, 1, 0], [1, 1, 1, 2], [-2, -1, 0, 0]] ,
         [[0, 0, 1, 0], [0, 1, 0, 0], [-1, -1, -1, -1], [0, -1, -1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [-1, 0, 0, 0]] ,
         [[0, 1, 0, 0], [0, -1, -1, 0], [0, -1, 0, -1], [0, 0, 0, 1], [0, 0, 1, 0], [-1, 0, 0, 0], [1, 0, 0, 0]] ,
         [[0, -1, -1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 0, -1], [1, 0, 0, 0]] ,
         [[0, 1, 0, 0], [0, 0, 0, 1], [1, -1, 1, -1], [0, -1, -1, 0], [0, 0, 1, 0], [1, 0, 0, 0], [-1, 0, 0, 0]] ,
         [[0, -1, -1, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 1, 1, -1], [-1, 0, 0, 0], [1, 0, 0, 0]] ,
         [[0, -1, -1, -1], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 1, 0, 0], [-1, 0, 0, 0], [1, 0, 0, 0]] ,
         [[0, -1, 0, -1], [0, 0, 0, 1], [0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 1, 0], [1, 0, -1, 0], [1, 0, 0, 0]] ,
         [[-1, 1, -1, -1], [0, 0, 0, 1], [0, 0, 1, 0], [-1, 0, 0, 0], [0, 1, 0, 0], [1, -1, 0, 0], [1, 0, 0, 0]] ,
         [[1, 0, 0, 0], [-1, 0, -1, -1], [0, 0, 0, 1], [0, 1, 1, 0], [-1, -1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0]] ,
         [[0, 1, 0, 0], [1, 0, 0, 0], [-1, -1, 0, 0], [0, -1, -1, -1], [0, 0, 0, 1], [0, 0, 1, 0], [1, 0, -1, 0]] ,
         [[1, 0, 0, 0], [0, 0, 0, 1], [1, 1, 1, -1], [-1, -1, 0, 0], [-1, 0, -1, 0], [0, 0, 1, 0], [0, 1, 0, 0]] ,
         [[0, 1, 0, 0], [1, 0, 0, 0], [-1, -1, 0, -1], [0, 0, 1, 0], [-1, -1, -1, 0], [0, 0, 0, 1], [0, 0, 1, -1]] ,
         [[1, 0, 0, 0], [-1, -1, -1, 0], [0, 0, 1, 0], [0, 1, 0, 1], [-1, 1, 0, -1], [0, 1, 0, 0], [0, 0, 0, 1]] ,
         [[-1, -1, -1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 1], [0, 0, 0, 1], [0, 1, 0, -1], [1, 0, 0, 0]] ,
         [[0, 0, -1, -1], [0, 0, 0, 1], [0, 0, 1, 0], [-1, 0, 0, 0], [0, -1, 0, 0], [0, 1, 0, 0], [1, 0, 0, 0]] ,
         [[0, -1, -1, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [2, -1, 0, -1], [-1, 0, 0, 0], [1, 0, 0, 0]] ,
         [[-1, 0, -1, 0], [0, 0, 0, 1], [0, 0, 1, 0], [1, -2, 0, -1], [1, 0, 0, 0], [0, 1, 0, 0], [1, -1, 0, 0]] ,
         [[-1, -1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [2, 0, -1, -1], [1, 0, -1, -1], [1, 0, 0, 0]] ,
         [[-1, -2, -1, -1], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [-2, 1, 0, 0], [1, 0, 0, 0], [-1, 1, 0, 0]] ,
         [[0, 1, 0, 0], [0, -1, 1, -1], [0, 0, 0, 1], [-1, -1, -1, 0], [0, 0, 1, 0], [1, 0, 0, 0], [-1, 0, 0, 0]] ,
         [[0, 0, 1, 0], [0, 0, 0, 1], [0, 1, -1, -1], [-1, -1, 0, 0], [0, 1, 0, 0], [1, 0, 0, 0], [-1, 0, 1, 0]] ,
         [[0, 0, 0, 1], [0, 0, 1, 0], [1, 0, -1, -1], [-1, 0, 0, 0], [0, -1, 0, 0], [0, 1, 0, 0], [1, 0, 0, 0]] ,
         [[0, 0, 0, 1], [0, 0, 1, 0], [1, -2, -1, -1], [0, 1, 0, 0], [1, 1, 0, 0], [-1, 0, 0, 0], [1, 0, 0, 0]] ,
         [[-1, -3, -2, -1], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [-1, 0, 0, 0], [0, -1, 0, 0]] ,
         [[-1, -1, 0, 0], [0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, -1], [0, 0, 1, 0], [1, 0, -1, 0]] ,
         [[1, 0, 0, 0], [-1, -1, 0, 0], [-1, 0, -1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 1, 0, -1]] ,
         [[1, 0, 0, 0], [-2, -1, -2, -1], [0, 0, 0, 1], [-1, -1, 0, 0], [-1, 0, -1, 0], [0, 0, 1, 0], [0, 1, 0, 0]] ,
         [[-2, -1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [1, 1, -1, -1], [-1, 0, 0, 0]] ,
         [[0, 0, 1, 0], [0, 1, 0, 0], [1, -1, -1, 0], [-1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 0, -1], [1, 0, 0, 0]] ,
         [[0, 0, 0, 1], [0, 0, 1, 0], [1, 1, -1, -1], [-1, 0, 0, 0], [0, -1, 0, 0], [0, 1, 0, 0], [1, 0, 0, 0]] ,
         [[0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1], [1, -1, -1, -1], [1, 0, 0, 0], [2, -1, -1, 0], [-1, 0, 0, 0]] ,
         [[0, 0, 1, 0], [1, 0, 0, 0], [-2, -2, -1, 0], [-2, -1, -1, -1], [0, 0, 0, 1], [0, 1, 0, 0], [-1, -1, 0, 0]] ,
         [[1, 0, 0, 0], [-2, -2, 0, -1], [0, 0, 0, 1], [0, 1, 0, 0], [-1, 0, -1, 0], [0, 0, 1, 0], [-1, -1, 0, 0]] ,
         [[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [2, -2, -2, -1], [-1, 0, 0, 0], [1, -1, -1, 0], [1, 0, 0, 0]] ,
         [[-1, -1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [2, 0, -2, -1], [1, 0, -1, 0]] ,
         [[0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [1, -1, -1, 0], [1, -1, 0, -1], [-1, 0, 0, 0], [1, 0, 0, 0]] ,
         [[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, -1, -1], [-1, 0, 0, 0], [0, 1, 0, 0], [1, -1, 0, 0], [1, 0, 0, 0]] ,
    ]

    for pts in all_pts:
        vc  = VectorConfiguration(pts)
        assert vc.ambient_dim == 4
        assert vc.dim == 4
        assert vc.is_solid()
        assert vc.is_totally_cyclic()
        assert not vc.is_acyclic()

        assert vc.labels == tuple(list(range(1,vc.size+1)))

        fan = vc.triangulate()
        assert fan.is_fine()
        assert fan.is_regular()

def test_all_triangulations():
    pts = [[-2, 2, 1, -1], [0, 0, 0, 1], [1, -2, 1, 1], [1, 1, -1, -1], [-1, 1, 1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [0, -1, -1, 0]]
    vc  = VectorConfiguration(pts)

    # construct all triangulations
    triangs = vc.all_triangulations(only_regular=False)
    assert len(triangs) == 3157

    # check the number of fine triangulations
    assert sum([t.is_fine() for t in triangs]) == 1814

    # check the number of regular triangulations
    assert sum([t.is_regular() for t in triangs]) == 1737

    # check the number of fine, regular triangulations
    assert sum([t.is_regular() and t.is_fine() for t in triangs]) == 794

    # check the number of fans that are also triangulations of the associated point configuration
    assert sum([t.respects_ptconfig() for t in triangs if t.is_regular()]) == 142

def test_all_circuits():
    pts = [[-5, -2, -3, -2], [0, 1, 0, 0], [0, 1, 2, 0], [1, 0, 0, 0], [1, 0, 1, 2], [-2, -1, -1, 0], [0, 1, 1, 0]]
    vc  = VectorConfiguration(pts)

    # construct all circuits
    circs = vc.circuits()
    supp  = {circ.Z for circ in circs}
    assert supp == {(1, 2, 3, 4, 5), (1, 4, 5, 7), (1, 5, 6), (2, 3, 4, 6), (2, 3, 7), (4, 6, 7)}

# smoke tests
def test_gale():
    pts = [[-5, -2, -3, -2], [0, 1, 0, 0], [0, 1, 2, 0], [1, 0, 0, 0], [1, 0, 1, 2], [-2, -1, -1, 0], [0, 1, 1, 0]]
    vc  = VectorConfiguration(pts)

    # smoke tests
    vc.gale()
    vc.flip_graph()
    vc.secondary_fan()
    vc.central_fan()

# custom labels
# -------------
_PTS = [[-5, -2, -3, -2], [0, 1, 0, 0], [0, 1, 2, 0], [1, 0, 0, 0], [1, 0, 1, 2], [-2, -1, -1, 0], [0, 1, 1, 0]]

def test_custom_labels_stored():
    labels = [10, 20, 30, 40, 50, 60, 70]
    vc = VectorConfiguration(_PTS, labels=labels)
    assert vc.labels == tuple(labels)

def test_custom_labels_in_fan():
    labels = [10, 20, 30, 40, 50, 60, 70]
    vc  = VectorConfiguration(_PTS, labels=labels)
    fan = vc.triangulate()
    # every label appearing in the fan cones must be one of the custom labels
    assert all(lbl in labels for c in fan.cones() for lbl in c)

def test_custom_labels_eq():
    labels = [10, 20, 30, 40, 50, 60, 70]
    vc1 = VectorConfiguration(_PTS, labels=labels)
    vc2 = VectorConfiguration(_PTS, labels=labels)
    assert vc1 == vc2
    assert not (vc1 != vc2)

def test_custom_labels_ne_default():
    vc_default = VectorConfiguration(_PTS)
    vc_custom  = VectorConfiguration(_PTS, labels=[10, 20, 30, 40, 50, 60, 70])
    assert vc_default != vc_custom

# labels_to_inds
# --------------
def test_labels_to_inds_default_single():
    vc = VectorConfiguration(_PTS)
    for i, lbl in enumerate(vc.labels):
        assert vc.labels_to_inds(lbl) == i

def test_labels_to_inds_default_batch():
    vc = VectorConfiguration(_PTS)
    assert vc.labels_to_inds([1, 3, 5]) == (0, 2, 4)

def test_labels_to_inds_custom_single():
    labels = [10, 20, 30, 40, 50, 60, 70]
    vc = VectorConfiguration(_PTS, labels=labels)
    for i, lbl in enumerate(labels):
        assert vc.labels_to_inds(lbl) == i

def test_labels_to_inds_custom_batch():
    labels = [10, 20, 30, 40, 50, 60, 70]
    vc = VectorConfiguration(_PTS, labels=labels)
    assert vc.labels_to_inds([10, 30, 50]) == (0, 2, 4)

def test_labels_to_inds_with_offset():
    vc = VectorConfiguration(_PTS)
    assert vc.labels_to_inds(1, offset=1) == 1
    assert vc.labels_to_inds([1, 2], offset=1) == (1, 2)


# grow4d backend
# --------------
# the 8 vertices of a cube, as a vector configuration: 64 fine triangulations
# and 166 in total, enough to truncate, and in dimension 3 rather than 4
GROW4D_VECS = np.ascontiguousarray(
    [[1,1,1], [1,1,-1], [1,-1,1], [1,-1,-1],
     [-1,1,1], [-1,1,-1], [-1,-1,1], [-1,-1,-1]], dtype=np.int32)


def test_all_triangulations_backends_agree():
    """The grow4d and flips backends must return the same triangulations."""
    vc = VectorConfiguration(GROW4D_VECS)

    for only_fine in (True, False):
        for only_regular in (True, False):
            g = vc.all_triangulations(only_fine=only_fine,
                                      only_regular=only_regular,
                                      backend="grow4d")
            f = vc.all_triangulations(only_fine=only_fine,
                                      only_regular=only_regular,
                                      backend="flips")
            assert {frozenset(t.cones()) for t in g} == \
                   {frozenset(t.cones()) for t in f}


def test_all_triangulations_unknown_backend():
    vc = VectorConfiguration([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],
                              [-1,-1,-1,-1]])
    with pytest.raises(ValueError):
        vc.all_triangulations(backend="nonesuch")


def test_all_triangulations_falls_back_to_flips():
    """A non-totally-cyclic VC has no complete fan, so grow4d must not run."""
    vc = VectorConfiguration([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    assert vc._grow4d_applies() is not None

    with pytest.warns(UserWarning, match="falling back"):
        got = vc.all_triangulations(backend="grow4d")
    assert {frozenset(t.cones()) for t in got} == \
           {frozenset(t.cones()) for t in vc.all_triangulations(backend="flips")}


def test_grow4d_is_seed_invariant():
    """
    The RNG only feeds the sampling stage, which is a shortcut for deciding
    that two cones overlap, never the authority on it. Undecided pairs fall
    through to an exact test, so the fans found must not depend on the seed.
    """
    from regfans.grow4d import grow4d

    vecs = GROW4D_VECS

    checksums = {grow4d(vecs, seed=sd)[4] for sd in (0, 1, 2, 12345)}
    assert len(checksums) == 1


def test_grow4d_respects_labels():
    """The kernel indexes rows; the cones handed back must be in labels."""
    vecs = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[-1,-1,-1,-1]]
    labels = [10, 20, 30, 40, 50]
    vc = VectorConfiguration(vecs, labels=labels)

    triangs = vc.all_triangulations(backend="grow4d")
    assert triangs
    for t in triangs:
        for cone in t.cones():
            assert set(cone).issubset(labels)


def test_grow4d_truncation():
    """Hitting a cap gives status -5 and a still-sliceable partial result."""
    from regfans.grow4d import grow4d

    vecs = GROW4D_VECS

    full, full_starts, num_full, status, _ = grow4d(vecs)
    assert status == 0
    assert num_full > 5

    simps, starts, num_fans, status, _ = grow4d(
        vecs, max_num_simps=10**6, max_num_fans=5)
    assert status == -5
    assert num_fans == 5

    # the closing sentinel must still be written, or the last fan's slice
    # would run off the end of the buffer
    assert starts[-1] == len(simps)

    # and every fan returned must be a genuine fan, not a partial complex
    everything = {frozenset(map(tuple, full[full_starts[i]:full_starts[i+1]]))
                  for i in range(num_full)}
    for i in range(num_fans):
        assert frozenset(map(tuple, simps[starts[i]:starts[i+1]])) in everything


def test_grow4d_count_only_matches():
    """Counting and materializing must agree, on both count and checksum."""
    from regfans.grow4d import grow4d

    vecs = GROW4D_VECS

    simps, starts, num_fans, status, checksum = grow4d(vecs)
    counted, count_status, count_checksum = grow4d(vecs, count_only=True)

    assert status == count_status == 0
    assert num_fans == counted
    assert checksum == count_checksum


def test_wall_normal_is_a_dependency_and_signed():
    """wall_normal returns the dependency, oriented by its first label."""
    vc = VectorConfiguration(GROW4D_VECS)
    labels = vc.labels[:4]                       # dim+1 = 4 labels in dim 3

    normal = vc.wall_normal(labels)
    assert normal is not None
    assert normal[0] > 0

    # it really is a linear dependency among those vectors
    combo = sum(c * v for c, v in zip(normal, vc.vectors(labels)))
    assert np.array_equal(combo, np.zeros(vc.ambient_dim))

    # primitive
    assert np.gcd.reduce(np.abs(normal)) == 1

    # reordering permutes the coefficients with it, up to the sign convention
    swapped = vc.wall_normal((labels[1], labels[0]) + tuple(labels[2:]))
    expected = [normal[1], normal[0]] + list(normal[2:])
    if expected[0] < 0:
        expected = [-c for c in expected]
    assert list(swapped) == expected


def test_wall_normal_none_when_not_spanning():
    """A degenerate set has no unique dependency, so no hyperplane."""
    vc = VectorConfiguration([[1,0,0], [0,1,0], [0,0,1], [-1,-1,-1],
                              [2,0,0], [0,2,0], [0,0,2], [-2,-2,-2]])
    # four labels whose vectors span only a 2-plane, so the dependency among
    # them is not unique
    labels = vc.vectors_to_labels([[1,0,0], [2,0,0], [0,1,0], [0,2,0]])
    assert vc.wall_normal(labels) is None


def test_spans_matches_rank():
    vc = VectorConfiguration(GROW4D_VECS)
    for cone in vc.triangulate().cones():
        assert vc.spans(cone) == util.is_full_rank(vc.vectors(cone))


def test_secondary_cone_hyperplanes_has_no_duplicate_rows():
    vc = VectorConfiguration(GROW4D_VECS)
    H = vc.triangulate().secondary_cone_hyperplanes()
    assert len(H) == len(set(map(tuple, H.tolist())))
