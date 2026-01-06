import pytest
import numpy as np
from regfans.vectorconfig import VectorConfiguration
from regfans import util

def test_neighbors():
	pts = [[-2, 2, 1, -1], [0, 0, 0, 1], [1, -2, 1, 1], [1, 1, -1, -1], [-1, 1, 1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [0, -1, -1, 0]]
	vc  = VectorConfiguration(pts)

	fan = vc.subdivide()
	sc  = fan.secondary_cone()
	assert util.cone_dim(H = sc) == vc.size

	# construct neighbors
	neighbs, circs = fan.neighbors()
	for neighb in neighbs:
		if not neighb.is_regular():
			continue
		n_sc = neighb.secondary_cone()
		assert util.cone_dim(H = n_sc) == vc.size
		assert util.cone_dim(H = np.vstack([sc,n_sc])) == vc.size-1

def test_flip_linear():
	pts = [[-2, 2, 1, -1], [0, 0, 0, 1], [1, -2, 1, 1], [1, 1, -1, -1], [-1, 1, 1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [0, -1, -1, 0]]
	vc  = VectorConfiguration(pts)

	# construct two different fans
	f1 = vc.subdivide()
	f2 = vc.subdivide(heights=[1]*vc.size)
	assert f1 != f2

	# smoke test
	eps = np.random.uniform(-1e-4, 1e-4, size=vc.size)
	f1.flip_linear(h_target=f2.heights()+eps)
