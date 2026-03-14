import pytest
import numpy as np
from regfans.vectorconfig import VectorConfiguration
from regfans import util

def test_neighbors():
    pts = [[-2, 2, 1, -1], [0, 0, 0, 1], [1, -2, 1, 1], [1, 1, -1, -1], [-1, 1, 1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [0, -1, -1, 0]]
    vc  = VectorConfiguration(pts)

    fan = vc.triangulate()
    sc  = fan.secondary_cone_hyperplanes()
    assert util.cone_dim(H = sc) == vc.size

    # construct neighbors
    neighbs, circs = fan.neighbors()
    for neighb in neighbs:
        if not neighb.is_regular():
            continue
        n_sc = neighb.secondary_cone_hyperplanes()
        assert util.cone_dim(H = n_sc) == vc.size
        assert util.cone_dim(H = np.vstack([sc,n_sc])) == vc.size-1

def test_flip_linear():
    pts = [[-2, 2, 1, -1], [0, 0, 0, 1], [1, -2, 1, 1], [1, 1, -1, -1], [-1, 1, 1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [0, -1, -1, 0]]
    vc  = VectorConfiguration(pts)

    # construct two different fans
    f1, f2 = list(vc.random_triangulations_fast(N=2, seed=0))
    assert f1 != f2

    # smoke test
    eps = np.random.uniform(-1e-4, 1e-4, size=vc.size)
    f1.flip_linear(h_target=f2.heights()+eps)

# Shared test input: 7 vectors in 4D, totally cyclic
_PTS = [
    [-5, -2, -3, -2], [0, 1, 0, 0], [0, 1, 2, 0],
    [1, 0, 0, 0], [1, 0, 1, 2], [-2, -1, -1, 0], [0, 1, 1, 0],
]


def test_vectors_getter():
    vc  = VectorConfiguration(_PTS)
    fan = vc.triangulate()

    # fan.vectors() returns the used vectors; fine triangulation uses all
    vecs = fan.vectors()
    np.testing.assert_array_equal(vecs, vc.vectors())

def test_cones_getter():
    vc  = VectorConfiguration(_PTS)
    fan = vc.triangulate()

    cones = fan.cones()

    # structural checks
    assert isinstance(cones, tuple)
    assert all(isinstance(c, tuple) for c in cones)

    # every label in every cone is a valid VC label
    assert all(lbl in vc.labels for c in cones for lbl in c)

    # triangulation => each cone has exactly dim vectors
    assert all(len(c) == fan.dim for c in cones)

def test_facets_getter():
    vc  = VectorConfiguration(_PTS)
    fan = vc.triangulate()

    facets = fan.facets()

    # structural checks
    assert isinstance(facets, dict)
    assert all(isinstance(f, tuple) for f in facets)
    assert all(len(f) == fan.dim - 1 for f in facets)
    assert all(isinstance(v, list) for v in facets.values())

    # every label in every facet is a valid VC label
    assert all(lbl in vc.labels for f in facets for lbl in f)

def test_is_triangulation_true():
    vc  = VectorConfiguration(_PTS)
    fan = vc.triangulate()
    assert fan.is_triangulation()

def test_is_triangulation_false():
    vc     = VectorConfiguration(_PTS)
    # a single cone containing all 7 labels in 4D is NOT a triangulation
    subdiv = vc.subdivide(cells=[list(vc.labels)])
    assert not subdiv.is_triangulation()
