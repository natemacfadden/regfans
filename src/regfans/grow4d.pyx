# grow4d.pyx
# Cython wrapper for grow4d

# import C types
# --------------
from libc.stdint cimport int32_t, uint32_t, uint64_t
from libc.stdlib cimport malloc, free
import numpy as np

# declare the external C function
# -------------------------------
# grow4d.h is header-only: this translation unit defines GROW4D_IMPLEMENTATION
# (via setup.py) so the kernel is compiled straight into the extension
cdef extern from "grow4d.h":
    int _grow4d_c "grow4d" (
        int * pts,
        int num_pts,
        int dim,
        int num_samples,
        uint64_t seed,
        int only_fine,
        int max_num_simps,
        uint32_t * simps,
        int max_num_fans,
        int * fan_starts,
        int * num_simps,
        uint64_t * num_fans,
        uint64_t * hash_out
    )

    int _grow4d_alloc_c "grow4d_alloc" (
        int * pts,
        int num_pts,
        int dim,
        int num_samples,
        uint64_t seed,
        int only_fine,
        uint32_t ** simps,
        int ** fan_starts,
        int * num_simps,
        uint64_t * num_fans,
        uint64_t * hash_out
    )

# Python-exposed wrapper
# ----------------------
def grow4d(int[:, ::1] pts,
           int max_num_simps = -1,
           int max_num_fans = -1,
           int num_samples = 4096,
           uint64_t seed = 0,
           bint only_fine = True,
           bint count_only = False) -> tuple:
    """
    Enumerate every simplicial fan supported on the vector configuration
    ``pts`` by growth with backtracking, exhaustively and without flips.

    Every such fan is found, including irregular ones and any lying in a
    component of the flip graph disconnected from the regular fans -- the
    case ``VectorConfiguration.all_triangulations`` cannot rule out.

    Totally cyclic gives the complete fans; acyclic (a homogenized point
    configuration) gives that configuration's triangulations. The support
    must be full-dimensional either way.

    Parameters
    ----------
    pts : int[:, ::1] of shape (num_pts, dim)
        The vector configuration, one vector per row. At most 64 vectors,
        in dimension at most 6.
    max_num_simps : int, optional
        Size of the output buffer, in simplices summed over all fans.
        Defaults to -1, which lets the kernel allocate and grow its own
        buffer, so the enumeration runs once. Give a bound only to cap the
        memory; exceeding it then returns status -5.
    max_num_fans : int, optional
        Maximum number of fans to materialize. Defaults to -1, as above.
        Exceeding it returns status -5.
    num_samples : int, optional
        Random points used to certify overlapping cone pairs. Defaults to
        4096. Only affects speed: undecided pairs fall through to an exact
        test either way.
    seed : int, optional
        Seed for those random points. Defaults to 0.
    only_fine : bool, optional
        If True (default), keep only fans using every vector.
    count_only : bool, optional
        If True, tally the fans without materializing them (no output buffer
        is allocated). ``max_num_simps`` / ``max_num_fans`` are ignored.

    Returns
    -------
    simps : ndarray of shape (num_simps, dim), dtype uint32
        The simplices of every fan, back to back, as row indices into
        ``pts``. Omitted when ``count_only``.
    fan_starts : ndarray of shape (num_fans + 1,), dtype int32
        Offsets into the rows of ``simps``: fan ``i`` is
        ``simps[fan_starts[i]:fan_starts[i+1]]``, a view. The closing entry
        is the total simplex count. Omitted when ``count_only``.
    num_fans : int
        Number of fans found.
    status : int
        Status code:
             0 : success
            -1 : no points input
            -2 : more than 64 points, dimension above 6, or allocation failure
            -3 : no initial simplex found (the configuration is degenerate)
            -5 : exceeded max_num_simps or max_num_fans
    checksum : int
        Order-independent checksum of the fan set, for comparing two runs
        without materializing either.

    Notes
    -----
    Under truncation (status -5) the fans returned are those found before the
    cap, in enumeration order, and the checksum covers only those. It is not
    the checksum of the full set.
    """
    cdef int num_pts = pts.shape[0]
    cdef int dim     = pts.shape[1]
    cdef int num_simps = 0
    cdef uint64_t num_fans = 0
    cdef uint64_t checksum = 0
    cdef int status
    cdef int *pts_ptr

    if num_pts == 0:
        raise ValueError("pts must have at least one row")
    pts_ptr = &pts[0, 0]

    # counting pass: no buffers, the kernel just tallies
    if count_only:
        status = _grow4d_c(pts_ptr, num_pts, dim, num_samples, seed,
                           1 if only_fine else 0, 0, NULL, 0, NULL,
                           &num_simps, &num_fans, &checksum)
        return int(num_fans), status, int(checksum)

    cdef uint32_t *g_simps = NULL
    cdef int *g_starts = NULL

    # no size given: let the kernel grow its own buffers, so the enumeration
    # runs once rather than once to count and again to fill
    if max_num_simps < 0 or max_num_fans < 0:
        status = _grow4d_alloc_c(pts_ptr, num_pts, dim, num_samples, seed,
                                 1 if only_fine else 0, &g_simps, &g_starts,
                                 &num_simps, &num_fans, &checksum)
        try:
            if status != 0:
                return (np.empty((0, dim), dtype=np.uint32),
                        np.zeros(1, dtype=np.int32), 0, status, int(checksum))
            return _copy_out(g_simps, g_starts, num_simps, num_fans, dim,
                             status, checksum)
        finally:
            free(g_simps)
            free(g_starts)

    if max_num_fans <= 0 or max_num_simps <= 0:
        raise ValueError("max_num_simps and max_num_fans must be positive")

    # allocate: fan_starts holds one entry per fan plus the closing total
    cdef size_t n_elem = <size_t>max_num_simps * <size_t>dim
    if dim > 0 and n_elem // <size_t>dim != <size_t>max_num_simps:
        raise MemoryError("max_num_simps * dim overflows size_t")

    cdef uint32_t *c_simps = <uint32_t *>malloc(n_elem * sizeof(uint32_t))
    if c_simps == NULL:
        raise MemoryError("Failed to allocate c_simps")

    cdef int *c_starts = <int *>malloc((<size_t>max_num_fans + 1) * sizeof(int))
    if c_starts == NULL:
        free(c_simps)
        raise MemoryError("Failed to allocate c_starts")

    status = _grow4d_c(pts_ptr, num_pts, dim, num_samples, seed,
                       1 if only_fine else 0, max_num_simps, c_simps,
                       max_num_fans, c_starts, &num_simps, &num_fans,
                       &checksum)

    out = _copy_out(c_simps, c_starts, num_simps, num_fans, dim,
                    status, checksum)
    free(c_simps)
    free(c_starts)
    return out


cdef _copy_out(uint32_t *c_simps, int *c_starts, int num_simps,
               uint64_t num_fans, int dim, int status, uint64_t checksum):
    """Copy the kernel's buffers into numpy's own memory, trimmed to size."""
    simps = np.empty((num_simps, dim), dtype=np.uint32)
    starts = np.empty(num_fans + 1, dtype=np.int32)

    cdef uint32_t[:, ::1] simps_view
    cdef int32_t[::1] starts_view
    cdef int i, j

    if num_simps > 0:
        simps_view = simps
        for i in range(num_simps):
            for j in range(dim):
                simps_view[i, j] = c_simps[<size_t>i * dim + j]

    starts_view = starts
    for i in range(<int>num_fans + 1):
        starts_view[i] = c_starts[i]

    return simps, starts, int(num_fans), status, int(checksum)
