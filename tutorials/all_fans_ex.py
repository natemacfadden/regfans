from regfans import VectorConfiguration

# construct all fans from a simple acyclic vector configuration
# =============================================================
# (N.B.: acyclic VC ~= PC. I.e., this is analogous to a 2D square)
pts  = [[0, 0, 1], [0, 1, 1], [1, 0, 1], [1, 1, 1]]
vc   = VectorConfiguration(pts)
fans = vc.all_triangulations()

print(f"The vector configuration {vc} has the following fans:")
for f in fans:
    print(f"a fan with cones {f.cones()}")

print()
print('-'*60)
print()

# construct all fans from a simple totally-cyclic vector configuration
# ====================================================================
pts  = [[x,y] for x in range(-1,1+1) for y in range(-1,1+1)]
vc   = VectorConfiguration(pts) # the point (0,0) will automatically be removed
G, fans, labels = vc.flip_graph(compute_node_labels=True) # this time, compute the fans via the flip_graph method... this also gives a formal graph object

print(f"The vector configuration {vc} has the following fans:")
for f in fans:
    print(f"a fan with cones {f.cones()}")

# plot the flip graph
# -------------------
try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

if HAS_MATPLOTLIB:
    import networkx as nx
    import matplotlib.pyplot as plt

    # construct node labels
    def node_label(label: dict) -> str:
        return (
            ("F" if label.get("fine") else "")
            + ("R" if label.get("regular") else "")
            + "S"
            + ("T" if label.get("respects_ptconfig") else "")
        )

    node_labels = {i: node_label(labels[i]) for i in range(len(labels))}

    # plot the graph
    pos = nx.spring_layout(G)

    plt.figure()
    nx.draw(
        G,
        pos,
        with_labels=False,
        node_size=200,
        node_color="lightgray"
    )

    nx.draw_networkx_labels(
        G,
        pos,
        labels=node_labels,
        font_size=6
    )

    plt.axis("off")
    plt.show()
