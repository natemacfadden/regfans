#ifndef GROW4D_H
#define GROW4D_H
#include <stdint.h>

/*
**Description:**
Exhaustively enumerates every simplicial fan supported on a vector
configuration. The C counterpart of grow4d.py, and the exhaustive counterpart
of grow2d.

Totally cyclic gives the complete fans; acyclic (a homogenized point
configuration) gives that configuration's triangulations. Only the support's
boundary is detected, so the cases between work too.

grow2d grows one triangulation: from a random unimodular triangle, repeatedly
select an exterior edge (used by only 1 simplex), scan the points and take
the first whose triangle satisfies the intersection property. It fails (-4) if some
exterior edge admits no point.

grow4d runs that same loop as a branch-and-backtrack:
    1) seed on every simplex containing a fixed point p interior to the
       support. A fan covers its support and its simplices meet only along
       faces, so p lies interior to exactly one simplex of any fan: every
       fan is reached once, with no loop over initial simplices,
    2) select an exterior face, f, used by only 1 simplex and not on the
       boundary of the support (fewest live candidates first),
    3) iterate over the points, recursing on every point whose simplex
       {f, pt} satisfies the intersection property,
    4) if some exterior face admits no point, backtrack (grow2d returns -4),
    5) record when no exterior face is left: every interior face is then used
       twice and every boundary face once, which is what covering the support
       means. Totally cyclic configurations have no boundary face, so this is
       the old "every face used exactly twice".
Fineness is checked on the completed complex: every point used by some simplex.

**Intersection property.** Cheap in 2D, the bottleneck in 4D, so precomputed
for all pairs in three stages:
    1) a separating plane from one simplex's own facet normals (bitmask test),
    2) a shared random point (an exact positive),
    3) otherwise: two simplices intersect iff 2d inequalities in d unknowns
       admit a solution; such a region if non-empty has a corner where d of the
       2d are exact, so all C(2d,d) choices are checked.

**Arguments:**
- `pts`:           Input vectors as a flat (num_pts x dim) row-major array.
- `num_pts`:       Number of input points.
- `dim`:           Ambient dimension (4 is what this is tested on).
- `num_samples`:   Random points used by stage 2.
- `seed`:          RNG seed.
- `only_fine`:     Keep only fans using every point.
// output
- `max_num_simps`: Max allowed number of simplices, summed over all fans.
- `simps`:         OUTPUT: Simplices of every fan, back to back, as a flat
                   (num_simps x dim) row-major array of point indices. Pass
                   NULL (with fan_starts NULL) to count without materializing.
- `max_num_fans`:  Max allowed number of fans. fan_starts must hold one more.
- `fan_starts`:    OUTPUT: Offsets into the rows of `simps`, so fan i is rows
                   fan_starts[i] up to fan_starts[i+1]. Length num_fans+1: the
                   closing entry is the total, hence the "+1".
- `num_simps`:     OUTPUT: Total simplices over all fans. Counted even when
                   `simps` is NULL, so a counting pass sizes the buffer for a
                   second, materializing one. May be NULL.
- `num_fans`:      OUTPUT: number of fans found.
- `hash_out`:      OUTPUT: order-independent checksum of the fan set.

**Returns:**
A status code:
     0: success
    -1: 0 points input
    -2: memory allocation problem
    -3: could not find an initial simplex containing the seed point
    -5: exceeded max_num_simps or max_num_fans

**Usage:**
Define GROW4D_IMPLEMENTATION in exactly one translation unit before including
this header; elsewhere include it plainly for the declaration alone.
*/
int grow4d(int *pts, int num_pts, int dim, int num_samples, uint64_t seed,
           int only_fine, int max_num_simps, uint32_t *simps,
           int max_num_fans, int *fan_starts, int *num_simps,
           uint64_t *num_fans, uint64_t *hash_out);

#ifdef GROW4D_IMPLEMENTATION

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TOL   1e-9      // interior tolerance
#define FTOL  1e-7      // feasibility tolerance for the corner check
#define MAXD  6         // max supported dimension


// EXTERNAL METHODS (from grow2d.h, itself from Blackman & Vigna 2019)
// ----------------
static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static uint64_t next(uint64_t rng_state[4]) {
    const uint64_t result = rotl(rng_state[0] + rng_state[3], 23) + rng_state[0];
    const uint64_t t = rng_state[1] << 17;

    rng_state[2] ^= rng_state[0];
    rng_state[3] ^= rng_state[1];
    rng_state[1] ^= rng_state[2];
    rng_state[0] ^= rng_state[3];
    rng_state[2] ^= t;
    rng_state[3] = rotl(rng_state[3], 45);

    return result;
}

static uint64_t splitmix64(uint64_t *state) {
    uint64_t z = (*state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}


// CONTEXT
// -------
/* Everything the growth recursion needs. Passed explicitly rather than held in
   file-scope state, so the enumerator is re-entrant. */
typedef struct {
    int       dim;          // ambient dimension
    int       num_pts;      // number of input vectors
    int       num_cands;    // number of candidate simplices
    int       num_faces;    // number of distinct (dim-1)-faces
    int       only_fine;    // keep only complexes using every point
    int      *simps;        // num_cands x dim   point indices
    int      *simp_faces;   // num_cands x dim   face id of each face
    int     **face_simps;   // num_faces -> simplices carrying that face
    int      *face_count;   // num_faces -> how many carry it
    uint64_t *conflict;     // num_cands x words bitset: may not coexist
    int       words;        // bitset words per simplex
    uint8_t  *used;         // num_faces -> current use count (0, 1 or 2)
    uint8_t  *is_bdry;      // num_faces -> face lies on the support boundary
    uint64_t *pt_cands;     // num_pts x words bitset: candidates using a point
    uint64_t  all_pts;      // every point, as a mask
    int      *stack;        // current complex, as simplex indices
    // output
    int       max_num_simps;  // capacity of simps_out, in rows
    int       max_num_fans;   // capacity of fan_starts, minus the sentinel
    uint32_t *simps_out;      // OUTPUT: fans back to back; NULL to only count
    int      *fan_starts;     // OUTPUT: fan i is rows [i], [i+1); NULL to skip
    int       num_simps;      // OUTPUT: rows written to simps_out
    int       overflow;       // a capacity was exceeded, so unwind and stop
    uint64_t  num_fans;       // OUTPUT
    uint64_t  hash;           // OUTPUT
} Ctx;

static inline int  get_bit(const uint64_t *bs, int i) {
    return (bs[i >> 6] >> (i & 63)) & 1ULL;
}
static inline void clear_bit(uint64_t *bs, int i) {
    bs[i >> 6] &= ~(1ULL << (i & 63));
}


// LINEAR ALGEBRA
// --------------
static double det_mat(const double *A, int n) {
    double M[MAXD * MAXD];
    memcpy(M, A, sizeof(double) * n * n);

    double det = 1.0;
    for (int c = 0; c < n; c++) {
        int piv = -1;
        double best = 0.0;
        for (int r = c; r < n; r++) {
            double v = fabs(M[r * n + c]);
            if (v > best) { best = v; piv = r; }
        }
        if (piv < 0 || best < 1e-14) return 0.0;

        if (piv != c) {
            for (int j = 0; j < n; j++) {
                double t = M[c * n + j];
                M[c * n + j] = M[piv * n + j];
                M[piv * n + j] = t;
            }
            det = -det;
        }
        det *= M[c * n + c];

        for (int r = c + 1; r < n; r++) {
            double f = M[r * n + c] / M[c * n + c];
            if (f != 0.0)
                for (int j = c; j < n; j++) M[r * n + j] -= f * M[c * n + j];
        }
    }
    return det;
}

static int inv_mat(const double *A, int n, double *out) {
    double M[MAXD * 2 * MAXD];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) M[i * 2 * n + j] = A[i * n + j];
        for (int j = 0; j < n; j++) M[i * 2 * n + n + j] = (i == j);
    }

    for (int c = 0; c < n; c++) {
        int piv = -1;
        double best = 0.0;
        for (int r = c; r < n; r++) {
            double v = fabs(M[r * 2 * n + c]);
            if (v > best) { best = v; piv = r; }
        }
        if (piv < 0 || best < 1e-12) return 0;

        if (piv != c)
            for (int j = 0; j < 2 * n; j++) {
                double t = M[c * 2 * n + j];
                M[c * 2 * n + j] = M[piv * 2 * n + j];
                M[piv * 2 * n + j] = t;
            }

        double d = M[c * 2 * n + c];
        for (int j = 0; j < 2 * n; j++) M[c * 2 * n + j] /= d;

        for (int r = 0; r < n; r++) {
            if (r == c) continue;
            double f = M[r * 2 * n + c];
            if (f != 0.0)
                for (int j = 0; j < 2 * n; j++)
                    M[r * 2 * n + j] -= f * M[c * 2 * n + j];
        }
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) out[i * n + j] = M[i * 2 * n + n + j];
    return 1;
}


// INTERSECTION PROPERTY, STAGE 1: a separating plane
// --------------------------------------------------
/* If every ray of simplex b lies on the non-positive side of one of simplex
   a's own facet normals, the two interiors cannot meet. Pure bitmask test.
   `nonpos[k*dim + i]` is the set of points on the <= 0 side of facet i of
   simplex k; `ray_mask[k]` is the set of rays of simplex k. */
static int separated(const uint64_t *ray_mask, const uint64_t *nonpos,
                     int dim, int a, int b) {
    for (int i = 0; i < dim; i++)
        if ((ray_mask[b] & ~nonpos[(size_t)a * dim + i]) == 0) return 1;
    for (int i = 0; i < dim; i++)
        if ((ray_mask[a] & ~nonpos[(size_t)b * dim + i]) == 0) return 1;
    return 0;
}


// INTERSECTION PROPERTY, STAGE 2: a shared random point
// -----------------------------------------------------
/* A sample point interior to both simplices is an exact positive. */
static int shares_sample(const uint64_t *member, int sample_words,
                         int a, int b) {
    for (int w = 0; w < sample_words; w++)
        if (member[(size_t)a * sample_words + w] &
            member[(size_t)b * sample_words + w]) return 1;
    return 0;
}


// INTERSECTION PROPERTY, STAGE 3: corner check
// --------------------------------------------
/* Two simplices meet iff {x : H_a x >= 1, H_b x >= 1} is non-empty. That
   region contains no line, so if non-empty it has a corner, found by making
   dim of the 2*dim inequalities exact. Checking all C(2*dim, dim) choices
   therefore decides it, with no solver. */
static int corner_feasible(const double *A, int dim,
                           const int (*corners)[MAXD], int num_corners) {
    for (int c = 0; c < num_corners; c++) {
        double S[MAXD * MAXD], Sinv[MAXD * MAXD], x[MAXD];

        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                S[i * dim + j] = A[corners[c][i] * dim + j];

        if (fabs(det_mat(S, dim)) <= 1e-12)  continue;
        if (!inv_mat(S, dim, Sinv))          continue;

        for (int i = 0; i < dim; i++) {
            x[i] = 0.0;
            for (int j = 0; j < dim; j++) x[i] += Sinv[i * dim + j];
        }

        int good = 1;
        for (int r = 0; r < 2 * dim && good; r++) {
            double v = 0.0;
            for (int j = 0; j < dim; j++) v += A[r * dim + j] * x[j];
            if (v < 1.0 - FTOL) good = 0;
        }
        if (good) return 1;
    }
    return 0;
}


// RECORDING
// ---------
/* A completed complex. Rejected here if not fine, i.e. if some point is used
   by no simplex -- grow2d's "doesn't cover any other points" in another
   guise. The checksum sums per-complex FNV hashes of the sorted simplex
   indices, so it does not depend on the order complexes are found in. */
static void record(Ctx *ctx, int depth, uint64_t used_pts) {
    if (ctx->only_fine && used_pts != ctx->all_pts) return;

    int *sorted = (int *)malloc(sizeof(int) * depth);
    memcpy(sorted, ctx->stack, sizeof(int) * depth);
    for (int i = 1; i < depth; i++) {           // insertion sort
        int k = sorted[i], j = i - 1;
        while (j >= 0 && sorted[j] > k) { sorted[j + 1] = sorted[j]; j--; }
        sorted[j + 1] = k;
    }

    uint64_t h = 1469598103934665603ULL;
    for (int i = 0; i < depth; i++) {
        h ^= (uint64_t)sorted[i];
        h *= 1099511628211ULL;
    }

    /* Write the fan out, sorted, as `depth` rows of dim point indices. The
       row where it starts is fan_starts[num_fans]; the matching end is
       written by the next fan, or by the sentinel once enumeration stops.
       Either capacity running out sets `overflow`, which unwinds grow(). */
    if (ctx->simps_out || ctx->fan_starts) {
        if (ctx->num_fans >= (uint64_t)ctx->max_num_fans ||
            ctx->num_simps + depth > ctx->max_num_simps) {
            ctx->overflow = 1;
            free(sorted);
            return;
        }
        if (ctx->fan_starts) ctx->fan_starts[ctx->num_fans] = ctx->num_simps;

        if (ctx->simps_out)
            for (int i = 0; i < depth; i++)
                for (int j = 0; j < ctx->dim; j++)
                    ctx->simps_out[(size_t)(ctx->num_simps + i) * ctx->dim + j]
                        = (uint32_t)ctx->simps[sorted[i] * ctx->dim + j];
    }
    ctx->num_simps += depth;    // tallied even when only counting, so a
    free(sorted);               // counting pass sizes the buffer exactly

    ctx->num_fans++;
    ctx->hash += h;
}


// GROWTH
// ------
/* grow2d selects an exterior edge, scans the points and takes the first
   whose triangle fits. Here every point that fits is recursed on, and a face with no
   candidate backtracks rather than returning -4. Faces are forward-checked, so
   a dead branch is cut at the first unsatisfiable face rather than later. */
static void grow(Ctx *ctx, int depth, uint64_t *live, uint64_t used_pts) {
    if (ctx->overflow) return;

    /* Fineness forward check. A point no cell uses yet can only be picked up
       by a candidate that is still live, so if some unused point has none, no
       extension of this complex is fine and the whole branch is dead. Without
       this the fineness test only happens at record(), i.e. after the branch
       has been walked to the end: `only_fine` would cost the full enumeration
       and merely report less of it. */
    if (ctx->only_fine) {
        uint64_t missing = ctx->all_pts & ~used_pts;
        while (missing) {
            int p = __builtin_ctzll(missing);
            missing &= missing - 1;

            const uint64_t *pc = ctx->pt_cands + (size_t)p * ctx->words;
            int reachable = 0;
            for (int w = 0; w < ctx->words; w++)
                if (live[w] & pc[w]) { reachable = 1; break; }
            if (!reachable) return;
        }
    }

    int best_face = -1, best_num = 1 << 30;

    for (int f = 0; f < ctx->num_faces; f++) {
        if (ctx->used[f] != 1) continue;        // not exterior
        if (ctx->is_bdry[f]) continue;          // on the support boundary: done

        int n = 0;
        for (int t = 0; t < ctx->face_count[f]; t++)
            if (get_bit(live, ctx->face_simps[f][t])) n++;

        if (n == 0) return;                     // grow2d would return -4
        if (n < best_num) { best_num = n; best_face = f; }
    }

    if (best_face < 0) { record(ctx, depth, used_pts); return; }  // none left

    uint64_t *next_live = (uint64_t *)malloc(sizeof(uint64_t) * ctx->words);
    for (int t = 0; t < ctx->face_count[best_face]; t++) {
        int k = ctx->face_simps[best_face][t];
        if (!get_bit(live, k)) continue;

        int ok = 1;
        for (int j = 0; j < ctx->dim; j++) {
            int fj = ctx->simp_faces[k * ctx->dim + j];
            // an interior face is shared by two cells, a boundary face by one
            if (ctx->used[fj] >= (ctx->is_bdry[fj] ? 1 : 2)) { ok = 0; break; }
        }
        if (!ok) continue;

        for (int j = 0; j < ctx->dim; j++) ctx->used[ctx->simp_faces[k * ctx->dim + j]]++;

        for (int w = 0; w < ctx->words; w++)
            next_live[w] = live[w] & ~ctx->conflict[(size_t)k * ctx->words + w];
        clear_bit(next_live, k);

        uint64_t next_used = used_pts;
        for (int j = 0; j < ctx->dim; j++)
            next_used |= 1ULL << ctx->simps[k * ctx->dim + j];

        ctx->stack[depth] = k;
        grow(ctx, depth + 1, next_live, next_used);

        for (int j = 0; j < ctx->dim; j++)
            ctx->used[ctx->simp_faces[k * ctx->dim + j]]--;
        if (ctx->overflow) break;
    }
    free(next_live);
}


// SETUP
// -----
/* All dim-subsets of the points whose vectors are linearly independent. */
static int build_simplices(Ctx *ctx, const int *pts) {
    int cap = 1;
    for (int i = 0; i < ctx->dim; i++)
        cap = cap * (ctx->num_pts - i) / (i + 1);

    ctx->simps = (int *)malloc(sizeof(int) * (size_t)cap * ctx->dim);
    if (!ctx->simps) return 0;

    int idx[MAXD];
    for (int i = 0; i < ctx->dim; i++) idx[i] = i;
    ctx->num_cands = 0;

    while (1) {
        double M[MAXD * MAXD];
        for (int i = 0; i < ctx->dim; i++)
            for (int j = 0; j < ctx->dim; j++)
                M[i * ctx->dim + j] = pts[idx[i] * ctx->dim + j];

        if (fabs(det_mat(M, ctx->dim)) > TOL) {
            for (int i = 0; i < ctx->dim; i++)
                ctx->simps[ctx->num_cands * ctx->dim + i] = idx[i];
            ctx->num_cands++;
        }

        int i = ctx->dim - 1;
        while (i >= 0 && idx[i] == ctx->num_pts - ctx->dim + i) i--;
        if (i < 0) break;
        idx[i]++;
        for (int j = i + 1; j < ctx->dim; j++) idx[j] = idx[j - 1] + 1;
    }

    ctx->words = (ctx->num_cands + 63) / 64;
    return 1;
}

/* Row i of the inverse of a simplex's ray matrix is its i-th facet normal. */
static double *build_normals(const Ctx *ctx, const int *pts) {
    int d = ctx->dim;
    double *H = (double *)malloc(sizeof(double) * (size_t)ctx->num_cands * d * d);
    if (!H) return NULL;

    for (int k = 0; k < ctx->num_cands; k++) {
        double M[MAXD * MAXD];
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                M[j * d + i] = pts[ctx->simps[k * d + i] * d + j];   // rays are columns
        inv_mat(M, d, H + (size_t)k * d * d);
    }
    return H;
}

/* The (dim-1)-subsets, and which simplices carry each.

   A face is identified by the bitmask of its labels, which is exact for every
   dim <= MAXD: num_pts <= 64, so the mask is one uint64. (Indexing a dense
   array by the labels themselves only works when a face has exactly 3 of
   them, i.e. dim 4, and would need num_pts^5 slots at dim 6.) The masks go in
   an open-addressed table; a face always has at least one label, so a zero
   mask means the slot is empty. */
static void build_faces(Ctx *ctx) {
    int d = ctx->dim;

    size_t want = (size_t)ctx->num_cands * d * 2;         // load factor <= 1/2
    size_t slots = 1;
    while (slots < want) slots <<= 1;
    uint64_t *hkey = (uint64_t *)calloc(slots, sizeof(uint64_t));
    int      *hval = (int *)malloc(sizeof(int) * slots);

    ctx->num_faces = 0;
    ctx->simp_faces = (int *)malloc(sizeof(int) * (size_t)ctx->num_cands * d);

    for (int k = 0; k < ctx->num_cands; k++)
        for (int drop = 0; drop < d; drop++) {
            uint64_t mask = 0;
            for (int i = 0; i < d; i++)
                if (i != drop) mask |= 1ULL << ctx->simps[k * d + i];

            size_t s = (size_t)((mask * 0x9E3779B97F4A7C15ULL) >> 40)
                       & (slots - 1);
            while (hkey[s] && hkey[s] != mask) s = (s + 1) & (slots - 1);
            if (!hkey[s]) { hkey[s] = mask; hval[s] = ctx->num_faces++; }
            ctx->simp_faces[k * d + drop] = hval[s];
        }
    free(hkey); free(hval);

    ctx->face_count = (int *)calloc(ctx->num_faces, sizeof(int));
    for (int k = 0; k < ctx->num_cands; k++)
        for (int j = 0; j < d; j++) ctx->face_count[ctx->simp_faces[k * d + j]]++;

    ctx->face_simps = (int **)malloc(sizeof(int *) * ctx->num_faces);
    for (int f = 0; f < ctx->num_faces; f++)
        ctx->face_simps[f] = (int *)malloc(sizeof(int) * ctx->face_count[f]);

    int *fill = (int *)calloc(ctx->num_faces, sizeof(int));
    for (int k = 0; k < ctx->num_cands; k++)
        for (int j = 0; j < d; j++) {
            int f = ctx->simp_faces[k * d + j];
            ctx->face_simps[f][fill[f]++] = k;
        }
    free(fill);
}

/* For each point, the candidates that use it. The forward check in grow()
   asks whether a point it still needs is reachable from what is left live. */
static void build_pt_cands(Ctx *ctx) {
    ctx->pt_cands = (uint64_t *)calloc((size_t)ctx->num_pts * ctx->words,
                                       sizeof(uint64_t));
    for (int k = 0; k < ctx->num_cands; k++)
        for (int i = 0; i < ctx->dim; i++) {
            int pt = ctx->simps[k * ctx->dim + i];
            ctx->pt_cands[(size_t)pt * ctx->words + (k >> 6)] |= 1ULL << (k & 63);
        }
}

/* Which faces lie on the boundary of the support.

   The face opposite ray i of simplex k spans a hyperplane whose normal is row
   i of that simplex's inverse, already computed and oriented so ray i is
   positive. The face is on the boundary of the support exactly when that
   hyperplane supports the whole configuration, i.e. no point sits strictly on
   the far side. Such a face is used once by a triangulation, not twice.

   A totally cyclic configuration has trivial dual cone, so no hyperplane
   supports it and no face comes back as boundary, leaving the growth below
   exactly the complete-fan enumeration it was. */
static void build_boundary(Ctx *ctx, const int *pts, const double *H) {
    int d = ctx->dim;
    ctx->is_bdry = (uint8_t *)calloc((size_t)ctx->num_faces, sizeof(uint8_t));
    uint8_t *done = (uint8_t *)calloc((size_t)ctx->num_faces, sizeof(uint8_t));

    for (int k = 0; k < ctx->num_cands; k++)
        for (int i = 0; i < d; i++) {
            int f = ctx->simp_faces[k * d + i];
            if (done[f]) continue;
            done[f] = 1;

            int bdry = 1;
            for (int m = 0; m < ctx->num_pts && bdry; m++) {
                double v = 0.0;
                for (int j = 0; j < d; j++)
                    v += H[((size_t)k * d + i) * d + j] * pts[m * d + j];
                if (v < -TOL) bdry = 0;
            }
            ctx->is_bdry[f] = (uint8_t)bdry;
        }
    free(done);
}

/* The pairwise intersection property, in the three stages above. */
static void build_conflicts(Ctx *ctx, const int *pts, const double *H,
                            const double *samples, int num_samples) {
    int d = ctx->dim, ns = ctx->num_cands;

    uint64_t *ray_mask = (uint64_t *)calloc(ns, sizeof(uint64_t));
    for (int k = 0; k < ns; k++)
        for (int i = 0; i < d; i++) ray_mask[k] |= 1ULL << ctx->simps[k * d + i];

    uint64_t *nonpos = (uint64_t *)calloc((size_t)ns * d, sizeof(uint64_t));
    for (int k = 0; k < ns; k++)
        for (int i = 0; i < d; i++) {
            uint64_t b = 0;
            for (int m = 0; m < ctx->num_pts; m++) {
                double v = 0.0;
                for (int j = 0; j < d; j++)
                    v += H[((size_t)k * d + i) * d + j] * pts[m * d + j];
                if (v <= TOL) b |= 1ULL << m;
            }
            nonpos[(size_t)k * d + i] = b;
        }

    int sw = (num_samples + 63) / 64;
    uint64_t *member = (uint64_t *)calloc((size_t)ns * sw, sizeof(uint64_t));
    for (int k = 0; k < ns; k++)
        for (int t = 0; t < num_samples; t++) {
            int in = 1;
            for (int i = 0; i < d && in; i++) {
                double v = 0.0;
                for (int j = 0; j < d; j++)
                    v += H[((size_t)k * d + i) * d + j] * samples[t * d + j];
                if (v <= TOL) in = 0;
            }
            if (in) member[(size_t)k * sw + (t >> 6)] |= 1ULL << (t & 63);
        }

    int num_corners = 1;
    for (int i = 0; i < d; i++) num_corners = num_corners * (2 * d - i) / (i + 1);
    int (*corners)[MAXD] = malloc(sizeof(int) * MAXD * num_corners);
    {
        int c[MAXD], m = 0;
        for (int i = 0; i < d; i++) c[i] = i;
        while (1) {
            for (int i = 0; i < d; i++) corners[m][i] = c[i];
            m++;
            int i = d - 1;
            while (i >= 0 && c[i] == d + i) i--;
            if (i < 0) break;
            c[i]++;
            for (int j = i + 1; j < d; j++) c[j] = c[j - 1] + 1;
        }
    }

    ctx->conflict = (uint64_t *)calloc((size_t)ns * ctx->words, sizeof(uint64_t));
    for (int a = 0; a < ns; a++)
        for (int b = a + 1; b < ns; b++) {
            if (separated(ray_mask, nonpos, d, a, b)) continue;

            int hit = shares_sample(member, sw, a, b);
            if (!hit) {
                double A[2 * MAXD * MAXD];
                for (int i = 0; i < d; i++)
                    for (int j = 0; j < d; j++) {
                        A[i * d + j]       = H[((size_t)a * d + i) * d + j];
                        A[(d + i) * d + j] = H[((size_t)b * d + i) * d + j];
                    }
                hit = corner_feasible(A, d, (const int (*)[MAXD])corners,
                                      num_corners);
            }

            if (hit) {
                ctx->conflict[(size_t)a * ctx->words + (b >> 6)] |= 1ULL << (b & 63);
                ctx->conflict[(size_t)b * ctx->words + (a >> 6)] |= 1ULL << (a & 63);
            }
        }

    free(ray_mask); free(nonpos); free(member); free(corners);
}


// ENTRY POINT
// -----------
int grow4d(int *pts, int num_pts, int dim, int num_samples, uint64_t seed,
           int only_fine, int max_num_simps, uint32_t *simps_out,
           int max_num_fans, int *fan_starts, int *num_simps,
           uint64_t *num_fans, uint64_t *hash_out) {
    if (num_pts <= 0)  return -1;
    if (num_pts > 64)  return -2;      // the fineness mask is 64-bit
    if (dim > MAXD)    return -2;

    Ctx ctx = {0};
    ctx.dim = dim;
    ctx.num_pts = num_pts;
    ctx.only_fine = only_fine;
    ctx.max_num_simps = max_num_simps;
    ctx.max_num_fans = max_num_fans;
    ctx.simps_out = simps_out;
    ctx.fan_starts = fan_starts;

    if (!build_simplices(&ctx, pts)) return -2;

    double *H = build_normals(&ctx, pts);
    if (!H) return -2;

    uint64_t rng_state[4], sm = seed;
    for (int i = 0; i < 4; i++) rng_state[i] = splitmix64(&sm);

    /* Sample directions inside the support. A random positive combination of
       the input vectors lies in the cone they generate, which is the
       support: all of R^dim when the configuration is totally cyclic, a
       pointed cone when it is not. A gaussian would do for the former but
       lands outside the latter almost every time, leaving the seed loop
       below with no simplex.
       Normalised because only the direction matters, these being cones, and a
       short vector would push the dot products below TOL. */
    double *samples = (double *)malloc(sizeof(double) * (size_t)num_samples * dim);
    for (int t = 0; t < num_samples; t++) {
        double *smp = samples + (size_t)t * dim;
        for (int j = 0; j < dim; j++) smp[j] = 0.0;

        for (int m = 0; m < num_pts; m++) {
            double c = (next(rng_state) >> 11) * 0x1.0p-53;
            for (int j = 0; j < dim; j++) smp[j] += c * pts[m * dim + j];
        }

        double nrm = 0.0;
        for (int j = 0; j < dim; j++) nrm += smp[j] * smp[j];
        nrm = sqrt(nrm);
        if (nrm > 0.0) for (int j = 0; j < dim; j++) smp[j] /= nrm;
    }

    build_conflicts(&ctx, pts, H, samples, num_samples);
    build_faces(&ctx);
    build_boundary(&ctx, pts, H);
    build_pt_cands(&ctx);

    ctx.all_pts = (num_pts == 64) ? ~0ULL : ((1ULL << num_pts) - 1);

    ctx.used  = (uint8_t *)calloc(ctx.num_faces, sizeof(uint8_t));
    ctx.stack = (int *)malloc(sizeof(int) * ctx.num_cands);

    /* Seed on the simplices containing the first sample point. That point is
       interior to the support, which every fan covers with simplices meeting
       only along faces, so it is interior to exactly one simplex of any fan:
       every fan is reached once, with no loop over initial simplices. */
    uint64_t *live = (uint64_t *)malloc(sizeof(uint64_t) * ctx.words);
    int num_seeds = 0;

    for (int k = 0; k < ctx.num_cands; k++) {
        int in = 1;
        for (int i = 0; i < dim && in; i++) {
            double v = 0.0;
            for (int j = 0; j < dim; j++)
                v += H[((size_t)k * dim + i) * dim + j] * samples[j];
            if (v <= TOL) in = 0;
        }
        if (!in) continue;
        num_seeds++;

        for (int w = 0; w < ctx.words; w++)
            live[w] = ~ctx.conflict[(size_t)k * ctx.words + w];
        clear_bit(live, k);

        uint64_t used_pts = 0;
        for (int j = 0; j < dim; j++) used_pts |= 1ULL << ctx.simps[k * dim + j];

        for (int j = 0; j < dim; j++) ctx.used[ctx.simp_faces[k * dim + j]]++;
        ctx.stack[0] = k;
        grow(&ctx, 1, live, used_pts);
        for (int j = 0; j < dim; j++) ctx.used[ctx.simp_faces[k * dim + j]]--;
    }

    /* Closing sentinel: fan i occupies rows fan_starts[i]..fan_starts[i+1],
       so the last fan needs an entry one past itself. Written even after an
       overflow, where num_fans counts the fans that did fit, so the fans
       returned alongside status -5 can still be sliced. */
    if (ctx.fan_starts) ctx.fan_starts[ctx.num_fans] = ctx.num_simps;

    *num_fans = ctx.num_fans;
    *hash_out = ctx.hash;
    if (num_simps) *num_simps = ctx.num_simps;

    free(live); free(H); free(samples);
    free(ctx.simps); free(ctx.simp_faces); free(ctx.face_count);
    for (int f = 0; f < ctx.num_faces; f++) free(ctx.face_simps[f]);
    free(ctx.face_simps); free(ctx.conflict); free(ctx.used); free(ctx.stack);
    free(ctx.is_bdry); free(ctx.pt_cands);

    if (ctx.overflow) return -5;
    return num_seeds ? 0 : -3;
}


#endif  /* GROW4D_IMPLEMENTATION */
#endif  /* GROW4D_H */
