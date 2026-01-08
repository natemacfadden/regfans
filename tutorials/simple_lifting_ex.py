from regfans import VectorConfiguration

# specify the points/vectors in the configuration
# (below is from the points of a 4D reflexive polytope but that doesn't matter...)
pts = [[1, -2, -1, -1], [1, 1, -1, 2], [-2, 0, 0, -1], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1], [1, -1, -1, 0], [1, 0, -1, 1], [1, 0, 0, 0], [-1, 0, 0, 0]]

# specify the heights
heights = None            # picks a default Delaunay triangulation
#heights = [1]*len(pts)    # picks a refinemenent of the central subdivision
#heights = [-1]*len(pts)    # these heights are invalid
#heights = [...]           # enter your own heights!


# study the vector configuration
# ------------------------------
vc  = VectorConfiguration(pts)
print(f"constructed the following VC: {vc}")

print()
print('-'*60)
print()

# optional: check properties of the VC
# ------------------------------------
# the defined VC is solid (full-dimensional) and totally-cyclic (i.e., support=\mathbb{R}^4)
print(f"is the VC solid? {vc.is_solid()}")
print(f"is the VC totally cyclic? {vc.is_totally_cyclic()}")
print(f"what are the bounding hyperplanes of the VC? {vc.support()}")

print()
print('-'*60)
print()

# construct/study a triangulation via lifting
# -------------------------------------------
# via lifting
fan = vc.subdivide(heights=heights)
print(f"constructed the following fan: {fan}")

print()
print('-'*60)
print()

# read some properies of the fan
print(f"the triangulation has cones {fan.cones()}")
print(f"some heights defining the triangulation are {fan.heights().tolist()}")
print(f"is the triangulation valid? {fan.is_valid()}")
print(f"is the triangulation fine? {fan.is_fine()}")
print(f"is the triangulation regular? {fan.is_regular()}")
print(f"does the triangulation define a triangulation of the associated PC? {fan.respects_ptconfig()}")
