import pytest
import numpy as np
from itertools import combinations
from regfans.vectorconfig import VectorConfiguration
from regfans.fan import Fan
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

def test_is_regular_false_santos_patching():
    # Santos patching example (Kaibel-Ziegler 2003): two fine *regular*
    # triangulations of [0,2]x[0,4] patched into a single *irregular*
    # triangulation of [0,4]^2... homogenize each lattice point (x,y) -> (x,y,1)
    verts  = sorted({(x, y) for x in range(5) for y in range(5)})
    points = [[x, y, 1] for (x, y) in verts]
    lbl    = {p: i + 1 for i, p in enumerate(verts)}

    triangles = [
        # left block [0,2]x[0,4] (regular on its own)
        ((0,0),(1,0),(0,1)), ((1,0),(1,1),(0,1)), ((1,1),(1,0),(2,0)), ((2,1),(2,0),(1,1)),
        ((0,1),(0,2),(1,1)), ((0,2),(1,1),(2,1)), ((0,2),(2,1),(1,2)), ((1,2),(2,2),(2,1)),
        ((0,2),(0,3),(1,3)), ((0,2),(1,2),(1,3)), ((0,3),(0,4),(1,4)), ((0,3),(1,3),(1,4)),
        ((1,2),(2,2),(2,3)), ((1,2),(2,3),(2,4)), ((1,2),(2,4),(1,3)), ((1,3),(1,4),(2,4)),
        # right block [2,4]x[0,4] (regular on its own)
        ((2,0),(3,0),(3,1)), ((2,0),(3,1),(3,2)), ((2,0),(2,1),(3,2)), ((2,1),(2,2),(3,2)),
        ((3,0),(4,0),(4,1)), ((3,0),(4,1),(3,1)), ((3,1),(4,1),(4,2)), ((3,1),(3,2),(4,2)),
        ((2,2),(2,3),(3,2)), ((2,3),(3,2),(4,2)), ((2,3),(3,3),(4,2)), ((3,3),(4,3),(4,2)),
        ((2,3),(2,4),(3,3)), ((2,4),(3,4),(3,3)), ((3,3),(3,4),(4,3)), ((3,4),(4,3),(4,4)),
    ]
    cones = [sorted(lbl[p] for p in t) for t in triangles]

    vc  = VectorConfiguration(points)
    fan = Fan(vc, cones)

    # it is a fine triangulation...
    assert fan.is_triangulation()
    assert fan.is_fine()

    # ... but irregular -- the whole point of the example
    assert not fan.is_regular()

    # patching property: each half on its own IS regular
    left  = [c for c, t in zip(cones, triangles) if all(p[0] <= 2 for p in t)]
    right = [c for c, t in zip(cones, triangles) if all(p[0] >= 2 for p in t)]
    assert Fan(vc, left).is_regular()
    assert Fan(vc, right).is_regular()

def test_cones_dim_matches_maximal_at_top_dim():
    vc  = VectorConfiguration(_PTS)
    fan = vc.triangulate()
    # simplicial => each maximal cone has exactly `dim` rays,
    # so cones(dim=fan.dim) reproduces the maximal cones
    assert set(fan.cones(dim=fan.dim)) == set(fan.cones())

def test_cones_dim_faces():
    vc  = VectorConfiguration(_PTS)
    fan = vc.triangulate()
    maximal = fan.cones()

    for d in range(1, fan.dim + 1):
        faces = fan.cones(dim=d)

        # each face has exactly d rays; results are de-duplicated
        assert all(len(f) == d for f in faces)
        assert len(set(faces)) == len(faces)

        # the faces are exactly the distinct d-subsets of the maximal cones
        expected = {c for m in maximal for c in combinations(m, d)}
        assert set(faces) == expected

def test_cones_dim_above_top_is_empty():
    vc  = VectorConfiguration(_PTS)
    fan = vc.triangulate()
    # no cone has more than `dim` rays, so there are no higher faces
    assert fan.cones(dim=fan.dim + 1) == ()

def test_cones_dim_respects_format_flags():
    vc  = VectorConfiguration(_PTS)
    fan = vc.triangulate()

    # labels: two rays per 2-cone
    assert all(len(c) == 2 for c in fan.cones(dim=2))

    # as_inds: two indices per 2-cone
    assert {len(c) for c in fan.cones(dim=2, as_inds=True)} == {2}

    # as_rays: each 2-cone is a 2-row array of generators
    assert all(len(r) == 2 for r in fan.cones(dim=2, as_rays=True))

    # ind_offset shifts every index by the offset
    base    = fan.cones(dim=2, as_inds=True)
    shifted = fan.cones(dim=2, as_inds=True, ind_offset=10)
    assert sorted(tuple(i + 10 for i in c) for c in base) == sorted(shifted)

def test_cones_dim_aliases_agree():
    vc  = VectorConfiguration(_PTS)
    fan = vc.triangulate()
    ref = set(fan.cones(dim=2))
    assert set(fan.simps(dim=2)) == ref
    assert set(fan.simplices(dim=2)) == ref
    assert set(fan.cells(dim=2)) == ref

def test_cones_dim_requires_simplicial():
    vc  = VectorConfiguration(_PTS)
    # a single cone containing all labels is non-simplicial
    sub = vc.subdivide(cells=[list(vc.labels)])
    assert not sub.is_triangulation()

    # dim=None still works on a non-simplicial fan...
    assert len(sub.cones()) == 1
    # ...but requesting faces by dim is not implemented => raises
    with pytest.raises(ValueError):
        sub.cones(dim=2)


def test_facets_getter_is_not_aliased():
    """
    facets() caches internally but must hand back the caller's own copy:
    mutating what one call returns cannot affect the next.
    """
    vc = VectorConfiguration(_PTS)
    fan = vc.triangulate()

    first = fan.facets()
    reference = {f: list(cs) for f, cs in first.items()}

    a_facet = next(iter(first))
    first[a_facet].append("scribble")
    first["not a facet"] = []

    assert fan.facets() == reference
