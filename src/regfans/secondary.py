# regfans/secondary.py
"""
Bulk regularity testing for a whole enumeration at once.

Testing each triangulation on its own means rebuilding, per fan, a list of
hyperplanes and a fresh LP. Over an enumeration that is enormously
redundant: at 13 rays, 13579 fans draw their walls from 229 distinct
label sets and their secondary-cone hyperplanes from 191 distinct rows.
This module works over that fixed universe instead.

Three things make it fast:

1. Walls come straight out of the enumerator's arrays, found by sorting
   rather than by asking each fan for its facets.
2. Each fan is then just a set of row indices into one shared matrix,
   held as a bitset.
3. Irregularity is inherited. If some rows of the shared matrix are
   positively dependent then ANY fan containing all of them is
   irregular, by Gordan's theorem -- the dependency is a certificate of
   infeasibility, and padding it with zeros certifies the larger system
   too. Every infeasible LP hands one back (its Farkas dual ray), so
   after a while most irregular fans are settled by a bitmask test
   rather than an LP. Regularity does NOT transfer this way, and cannot:
   if one fan's rows are a subset of another's then its cone contains
   the other's, and distinct triangulations have interior-disjoint
   secondary cones, so the two fans coincide. Every regular fan is paid
   for with its own LP.
"""
from __future__ import annotations

import numpy as np
import highspy

__all__ = ["walls", "regular_mask"]


def walls(simps: np.ndarray, fan_starts: np.ndarray, dim: int) -> tuple:
    """
    The walls of every fan in an enumeration.

    A wall is an interior facet: one shared by exactly two of the fan's
    simplices. It is returned as the set of dim+1 labels spanning it, i.e.
    the shared facet plus the two opposite vertices, encoded as a bitmask
    over labels.

    Parameters
    ----------
    simps : np.ndarray of shape (num_simps, dim)
        The simplices of every fan, back to back, as returned by grow4d.
        Entries are indices, so must be below 64.
    fan_starts : np.ndarray of shape (num_fans + 1,)
        Offsets into the rows of `simps`; fan i is rows fan_starts[i] up
        to fan_starts[i+1].
    dim : int
        The ambient dimension.

    Returns
    -------
    wall_masks : np.ndarray
        One bitmask per wall, over all fans, concatenated.
    apex_masks : np.ndarray
        The two opposite vertices of each wall, as a bitmask. These fix the
        orientation of the hyperplane the wall imposes and so cannot be
        recovered from `wall_masks`: the same labels can span a wall in two
        fans with the apexes on opposite sides of the dependency, which is
        the same inequality reversed.
    wall_starts : np.ndarray of shape (num_fans + 1,)
        Offsets into `wall_masks`, in the same convention as fan_starts.
    """
    simps = np.asarray(simps, dtype=np.int64)
    fan_starts = np.asarray(fan_starts, dtype=np.int64)
    num_simps = simps.shape[0]
    num_fans = len(fan_starts) - 1

    # which fan each simplex belongs to
    fan_of = np.repeat(np.arange(num_fans), np.diff(fan_starts))

    # every facet of every simplex, as (simplex mask minus one vertex).
    # Working in bitmasks keeps a facet a single integer, so the grouping
    # below is a sort rather than a dictionary.
    simp_mask = np.bitwise_or.reduce(np.int64(1) << simps, axis=1)
    facet = simp_mask[:, None] ^ (np.int64(1) << simps)      # (num_simps, dim)
    opp = simps                                              # opposite vertex

    facet = facet.ravel()
    opp = (np.int64(1) << opp).ravel()
    owner = np.repeat(fan_of, dim)

    # sort so that the two occurrences of a shared facet land side by side.
    # One combined key beats lexsort's two passes, and both parts are small:
    # the facet mask spans the labels, of which there are at most 64
    order = np.argsort(owner * (np.int64(1) << np.int64(simps.max() + 1))
                       + facet, kind="stable")
    facet, opp, owner = facet[order], opp[order], owner[order]

    # a facet used twice within one fan is a wall; used once, it is on the
    # boundary. A complete fan has no boundary, but this does not assume so
    shared = (facet[:-1] == facet[1:]) & (owner[:-1] == owner[1:])
    idx = np.flatnonzero(shared)

    apex_masks = opp[idx] | opp[idx + 1]
    wall_masks = facet[idx] | apex_masks
    wall_owner = owner[idx]

    counts = np.bincount(wall_owner, minlength=num_fans)
    wall_starts = np.concatenate(([0], np.cumsum(counts)))

    return wall_masks, apex_masks, wall_starts.astype(np.int64)


def _wall_keys(wall_masks: np.ndarray, apex_masks: np.ndarray,
               num_labels: int) -> np.ndarray:
    """
    One integer per wall, identifying it together with its orientation.

    The apexes are a subset of the wall, so both fit side by side in an
    int64 for any configuration the enumerator accepts.
    """
    return (wall_masks << np.int64(num_labels)) | apex_masks


def _hyperplanes(vc, wall_masks: np.ndarray, apex_masks: np.ndarray) -> tuple:
    """
    The distinct secondary-cone hyperplanes imposed by a set of walls.

    Parameters
    ----------
    vc : VectorConfiguration
        The configuration the walls live in.
    wall_masks : np.ndarray
        Bitmasks of walls, as returned by `walls`. Repeats are fine.
    apex_masks : np.ndarray
        The matching apex bitmasks, which orient each hyperplane.

    Returns
    -------
    rows : np.ndarray of shape (num_rows, vc.size)
        The distinct hyperplanes, as {x: rows @ x >= 0}.
    row_of_wall : dict
        Mapping from (wall bitmask, apex bitmask) to its row in `rows`, or
        -1 for a wall that imposes nothing (its labels do not span).
    """
    labels = vc.labels
    label_to_ind = vc.labels_to_inds_dict

    rows = []
    seen = {}
    row_of_wall = {}

    n_lbl = len(labels)
    keys = np.unique(_wall_keys(wall_masks, apex_masks, n_lbl))
    for key in keys:
        mask, apex = int(key) >> n_lbl, int(key) & ((1 << n_lbl) - 1)
        spanning = [labels[i] for i in range(len(labels)) if (mask >> i) & 1]

        # orient by an apex: the two apexes always sit on the same side of
        # the dependency, so either fixes the same sign, and the lower one
        # makes the choice deterministic
        first = next(labels[i] for i in range(len(labels)) if (apex >> i) & 1)
        spanning.remove(first)
        spanning = (first,) + tuple(spanning)

        normal = vc.wall_normal(spanning)
        if normal is None:
            row_of_wall[(mask, apex)] = -1
            continue

        row = [0] * vc.size
        for lbl, coeff in zip(spanning, normal):
            row[label_to_ind[lbl]] = coeff
        row = tuple(row)

        if row not in seen:
            seen[row] = len(rows)
            rows.append(row)
        row_of_wall[(mask, apex)] = seen[row]

    return np.array(rows, dtype=float).reshape(-1, vc.size), row_of_wall


def regular_mask(vc,
    simps: np.ndarray,
    fan_starts: np.ndarray,
    verbosity: int = 0) -> np.ndarray:
    """
    Which fans of an enumeration are regular.

    Equivalent to calling `Fan.is_regular` on each fan in turn, but done
    over the whole enumeration at once. See the module docstring.

    Only valid for FINE fans. A fan that omits a vector imposes further
    conditions, from inserting that vector into each of its cones, which
    are not walls and so are not seen here.

    Parameters
    ----------
    vc : VectorConfiguration
        The configuration the fans triangulate.
    simps : np.ndarray of shape (num_simps, dim)
        The simplices of every fan, back to back, as returned by grow4d.
    fan_starts : np.ndarray of shape (num_fans + 1,)
        Offsets into the rows of `simps`.
    verbosity : int, optional
        The verbosity level. Higher is more verbose. Defaults to 0.

    Returns
    -------
    out : np.ndarray of bool, shape (num_fans,)
        True where the corresponding fan is regular.
    """
    num_fans = len(fan_starts) - 1
    if num_fans == 0:
        return np.zeros(0, dtype=bool)

    wall_masks, apex_masks, wall_starts = walls(simps, fan_starts, vc.ambient_dim)
    rows, row_of_wall = _hyperplanes(vc, wall_masks, apex_masks)
    num_rows, dim = rows.shape

    if verbosity >= 1:
        print(f"{num_fans} fans over {num_rows} distinct hyperplanes")

    # every wall's row, without a dictionary lookup per wall: the walls
    # collapse to a few hundred distinct (mask, apex) pairs, so look those up
    # and let the inverse index carry the answer back to every occurrence
    n_lbl = len(vc.labels)
    keys = _wall_keys(wall_masks, apex_masks, n_lbl)
    uniq, inv = np.unique(keys, return_inverse=True)
    uniq_row = np.array([row_of_wall[(int(k) >> n_lbl,
                                     int(k) & ((1 << n_lbl) - 1))]
                         for k in uniq])
    wall_row = uniq_row[inv.ravel()]

    owner_of_wall = np.repeat(np.arange(num_fans), np.diff(wall_starts))
    keep = wall_row >= 0

    # each fan as a set of rows, and as a bitset over them
    incidence = np.zeros((num_fans, num_rows), dtype=bool)
    incidence[owner_of_wall[keep], wall_row[keep]] = True

    words = (num_rows + 63) // 64
    fan_bits = np.zeros((num_fans, words), dtype=np.uint64)
    for j in range(num_rows):
        fan_bits[incidence[:, j], j >> 6] |= np.uint64(1) << np.uint64(j & 63)

    fan_rows = [np.flatnonzero(r) for r in incidence]

    # one LP, reused: rows are switched in and out by their bounds, so the
    # basis carries over from the previous fan rather than starting cold
    inf = highspy.kHighsInf
    lp = highspy.Highs()
    if verbosity < 2:
        lp.silent()
    # the problems are tiny and solved thousands of times over; presolving
    # each one costs more than it saves, and it discards the warm basis
    lp.setOptionValue("presolve", "off")
    lp.addVars(dim, np.full(dim, -inf), np.full(dim, inf))
    starts = (np.arange(num_rows) * dim).astype(np.int32)
    index = np.tile(np.arange(dim, dtype=np.int32), num_rows)
    lp.addRows(num_rows, np.full(num_rows, -inf), np.full(num_rows, inf),
               num_rows * dim, starts, index, rows.ravel())
    cols = np.arange(dim, dtype=np.int32)

    certs = np.zeros((0, words), dtype=np.uint64)
    out = np.zeros(num_fans, dtype=bool)
    active = set()
    n_skipped = 0

    for i in range(num_fans):
        # inherited irregularity: any cached dependency sitting inside this
        # fan's rows already proves it infeasible
        if len(certs) and np.any(np.all((certs & fan_bits[i]) == certs, axis=1)):
            n_skipped += 1
            continue

        want = set(fan_rows[i].tolist())
        for j in want - active:
            lp.changeRowBounds(j, 1.0, inf)
        for j in active - want:
            lp.changeRowBounds(j, -inf, inf)
        active = want

        idx = fan_rows[i]
        lp.changeColsCost(dim, cols,
                          np.ascontiguousarray(rows[idx].sum(axis=0) / len(idx)))
        lp.run()

        if lp.getModelStatus() == highspy.HighsModelStatus.kOptimal:
            out[i] = True
            continue

        # infeasible: keep the Farkas ray's support, which certifies every
        # fan that contains those rows
        _, exists, ray = lp.getDualRay()
        if exists:
            y = np.abs(np.asarray(ray, dtype=float))
            support = [int(j) for j in idx if y[j] > 1e-9]
            if support:
                cert = np.zeros(words, dtype=np.uint64)
                for j in support:
                    cert[j >> 6] |= np.uint64(1) << np.uint64(j & 63)
                certs = np.vstack([certs, cert])

    if verbosity >= 1:
        print(f"{int(out.sum())} regular; {n_skipped} irregular fans settled "
              f"by certificate, {num_fans - n_skipped} LPs solved")

    return out
