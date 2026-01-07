# =============================================================================
#    Copyright (C) 2025  Nate MacFadden for the Liam McAllister Group
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
# =============================================================================
#
# -----------------------------------------------------------------------------
# Description:  This module contains utilities for vector configuration
#               computations.
# -----------------------------------------------------------------------------

# external imports
import fractions
import functools
import math
import numpy as np
from ortools.linear_solver import pywraplp
from typing import Union

# basic geometry
# --------------
def lerp(p0: "ArrayLike", p1: "ArrayLike", t: float) -> "ArrayLike":
    """
    **Description:**
    Computes the point specified by t along the line passing through p0 and p1.

    Particular values:
        -) t=0   -> p0
        -) t=0.5 -> (p0+p1)/2
        -) t=1   -> p1

    **Arguments:**
    - `p0`: One point.
    - `p1`: The other point.
    - `t`:  Parameter specifing where along the line Conv({p0, p1}) to return.

    **Returns:**
    The point p0 + t*(p1-p0).
    """
    # input sanitization
    p0 = np.array(p0)
    p1 = np.array(p1)

    # output
    return p0 + t*(p1-p0)

def hyperplane_inter(p0: "ArrayLike", p1: "ArrayLike", n: "ArrayLike") -> float:
    """
    **Description:**
    Computes the distance, t, along the line p0 + t*(p1-p0) that intersects the
    hyperplane {x: n.x=0}.

    This is simply computed as
        0 = n.[p0 + t*(p1-p0)]
        0 = n.p0 + t*n.(p1-p0)
        t = -n.p0 / n.(p1-p0)

    **Arguments:**
    - `p0`: One point.
    - `p1`: The other point.
    - `n`:  The hyperplane normal.

    **Returns:**
    The distance t such that p0 + t*(p1-p0) lands upon the hyperplane.
    """
    # input sanitization
    p0 = np.array(p0)
    p1 = np.array(p1)
    n  = np.array(n)

    # comute distance
    denom = np.dot(n,p1-p0)

    if denom==0:
        # line doesn't intersect hyperplane
        return None
    else:
        numer = np.dot(n,p0)
        return -numer/denom

def first_hit(
    p0: "ArrayLike",
    p1: "ArrayLike",
    H: "ArrayLike",
    max_dist: float=1,
    verbosity: int=0) -> (int, float):
    """
    **Description:**
    Given a point p0 in a convex cone {x: Hx>=0}, find the first hyperplane hit
    along the direction (p1-p0). I.e, the first intersection of the ray
    {p0+t*(p1-p0): t>=0} with the cones bounding hyperplanes.

    Allow violated hyperplanes (i.e., n such that n.p0 < 0) but ignore them.

    **Arguments:**
    - `p0`:         One point.
    - `p1`:         The other point.
    - `H`:          An array of hyperplane normals (as rows).
    - `max_dist`:   Only consider intersections along the line segment
                    [p0, p0+max_dist*(p1-p0)]
    - `verbosity`:  The verbosity level. Higher is more verbose.

    **Returns:**
    The index, i, of the first-hit hyperplane.
    The distance, t, such that H[i].lerp(p0,p1,t) = 0.
    """
    # input sanitization
    p0 = np.array(p0)
    p1 = np.array(p1)
    H  = np.array(H)

    # defaults
    first_hit_ind = -1
    first_hit_dist = 2*max_dist

    # iterate over hyperplanes
    for i,n in enumerate(H):
        if np.dot(n,p0)<=0:
            # skip intersections in the wrong direction
            continue

        # compute the distance along the ray
        dist = hyperplane_inter(p0, p1, n)
        if verbosity >= 1:
            print(f"i={i} has dist={dist}...")

        # if the distance is permissible (positive and <= max_dist)
        # and if the distance is shorter than previously found, save it
        if 0<dist<=max_dist and dist<=first_hit_dist:
            first_hit_ind = i
            first_hit_dist = dist
    
    # return
    if first_hit_ind == -1:
        # no valid intersections found
        return None, None
    else:
        return first_hit_ind, first_hit_dist

# cone geometry
# -------------
def dual_cone(data: "ArrayLike") -> "ArrayLike":
    """
    **Description:**
    Compute the data of the cone dual to the input 'primal' cone.

    This can be thought of in a couple of equivalent ways, summarized in the
    following table. E.g., if rays of the primal are input, then the
    hyperplanes of the primal are output (or, equivalently, the rays of the
    dual).

    INPUT       | PRIMAL OUTPUT | DUAL OUTPUT
    -----------------------------------------
    rays        | hyperplanes   | rays
    hyperplanes | rays          | hyperplanes


    For simplicitly in the following discussion, take the convention that one
    maps hyperplanes of the primal to rays of the primal.

    **Arguments:**
    - `data`: An array whose rows represent rays of the primal cone. (see table)

    **Returns:**
    An array whose rows represent hyperplanes of the primal cone. (see table)
    """
    # check the ppl install
    try:
        import ppl
    except ImportError as e:
        raise ImportError(
            "pplpy is required for computations of dual cones"
            "Please install via conda: conda install -c conda-forge pplpy"
        ) from e

    # input sanitization
    data = np.array(data)

    # define polyhedron in ppl
    cone = ppl.C_Polyhedron(data.shape[1])
    for row in data:
        ineq = ppl.Linear_Expression(row.tolist(), 0)
        cone.add_constraint(ppl.Constraint(ineq >= 0))

    # compute the rays
    rays = []
    for g in cone.generators():
        if g.is_ray():
            rays.append([int(c) for c in g.coefficients()])
        elif g.is_line():
            # lineality space... add both signs
            rays.append([int(c) for c in g.coefficients()])
            rays.append([-int(c) for c in g.coefficients()])

    # return
    return rays

def cone_dim(*, R: "ArrayLike"=None, H:"ArrayLike"=None) -> int:
    """
    **Description:**
    Return the dimension of the cone.

    The cone is either specified via rays,
        {R.T @ lambda: lambda>=0},
    or via hyperplanes,
        {x: H @ x>=0}.

    **Arguments:**
    - `R`: The rays of the cone as rows.
    - `H`: The hyperplanes defining the cone.

    **Returns:**
    The dimension of the cone
    """
    assert (R is None) ^ (H is None)

    if R is None:
        R = np.array(dual_cone(H))
    else:
        R = np.array(R)

    # return
    return np.linalg.matrix_rank(R.T)

def is_solid(*, R: "ArrayLike"=None, H:"ArrayLike"=None) -> int:
    """
    **Description:**
    Return whether the cone is full-dimensional.

    The cone is either specified via rays,
        {R.T @ lambda: lambda>=0},
    or via hyperplanes,
        {x: H @ x>=0}.

    **Arguments:**
    - `R`: The rays of the cone as rows.
    - `H`: The hyperplanes defining the cone.

    **Returns:**
    The dimension of the cone
    """
    assert (R is None) ^ (H is None)

    if R is None:
        # try to find a point in the strict interior
        H = np.array(H)
        return find_interior_point(H=H) is not None
    else:
        # know rays -> sum of rays sould be in strict interior
        R = np.array(R)
        H = np.array(dual_cone(R))
        return np.all(H@R.sum(axis=0) > 0.5)

def contains(*,
    p: "ArrayLike",
    R:"ArrayLike" = None,
    H: "ArrayLike" = None) -> bool:
    """
    **Description:**
    Return if the point p is contained in the cone.

    The cone is either specified via rays,
        {R.T @ lambda: lambda>=0},
    or via hyperplanes,
        {x: H @ x>=0}.

    **Arguments:**
    - `R`: The rays of the cone as rows.
    - `H`: The hyperplanes defining the cone.

    **Returns:**
    Whether p is contained in the cone.
    """
    assert (R is None) ^ (H is None)

    if H is None:
        H = np.array(dual_cone(R))
    else:
        H = np.array(H)

    # return
    return np.all(H@p >= 0)

def find_interior_point(*,
    R: "ArrayLike" = None,
    H: "ArrayLike" = None,
    stretching: float = 1,
    nonneg: bool = False,
    verbosity: int = 0) -> Union["ArrayLike", None]:
    """
    **Description:**
    Returns a point p in the relative interior of a cone. The cone can be
    specified either via its rays or its generators.

    If no point p exists, return `None`.

    Modified from CYTools' `Cone.find_interior_point`.

    **Arguments:**
    - `R`:          Generators defining the cone.
    - `H`:          Hyperplanes defining the cone.
    - `stretching`: How far p must be from any hyperplane.
    - `nonneg`:     Whether to restrict to non-negative vectors.
    - `verbosity`:  The verbosity level.

    **Returns:**
    A point p in the strict interior.
    """
    if not ((R is None) ^ (H is None)):
        raise ValueError("Either R or H can be set, but not both!")

    # simple method if R is set
    if H is None:
        if nonneg:
            raise NotImplementedError("nonneg=True is not allowed if R is set")
        return R.sum(axis=0)

    # preliminary/sanitization
    H   = np.array(H)
    m,n = H.shape

    # create the solver
    solver = pywraplp.Solver.CreateSolver("GLOP")
    if verbosity >= 1:
        solver.EnableOutput()

    # create the variables
    inf = solver.infinity()
    if nonneg:
        low, up = 0,    inf
    else:
        low, up = -inf, inf
    p   = [solver.NumVar(low, up, f"p_{i}") for i in range(n)]

    # impose H @ p >= stretching
    for i in range(m):
        cons = solver.Constraint(stretching, solver.infinity())
        for j in range(n):
            cons.SetCoefficient(p[j], float(H[i,j]))

    # pick a semi-arbitrary grading vector and return p with minimal degree
    grading = H.sum(axis=0) # guaranteed to be in the interior of the dual cone
    solver.Minimize(sum(gj*pj for gj,pj in zip(grading,p)))

    # solve/parse
    status = solver.Solve()
    if status in (solver.FEASIBLE, solver.OPTIMAL):
        solution = np.array([pi.solution_value() for pi in p])
        return solution
    elif status == solver.INFEASIBLE:
        if verbosity>=1:
            warnings.warn("Cone is not full-dimensional")
        return None
    else:
        warnings.warn(f"Unexpected error")
        return None


# basic math - UNUSED
# ----------
def gcd(vals: list[float], max_denom: float=10**6) -> float:
    """
    **Description:**
    Computes the 'GCD' of a collection of floating point numbers.
    This is the smallest number, g, such that g*values is integral.

    This is computed by
        1) converting `values` to be rational [n0/d0, n1/d1, ...],
        2) computing the LCM, l, of [d0, d1, ...],
        3) computing the GCD, g', of [l*n0/d0, l*n1/d1, ...], and then
        4) returning g=g'/l.

    **Arguments:**
    - `vals`:      The numbers to compute the GCD of.
    - `max_denom`: Assert |di| <= max_denom

    **Returns:**
    The minimum number g' such that g'*vals is integral.
    """
    # compute the rational representation
    rat   = [fractions.Fraction(v).limit_denominator(max_denom) for v in vals]
    numer = [r.numerator   for r in rat]
    denom = [r.denominator for r in rat]

    # get the relevant LCM, GCD
    l     = functools.reduce(math.lcm, denom)
    gprime= functools.reduce(math.gcd, [n*(l//d) for n,d in zip(numer,denom)])

    # return the GCD
    if gprime%l == 0:
        # integral
        return gprime//l
    else:
        return gprime/l

def primitive(vec: list[float], max_denom=10**10):
    """
    **Description:**
    Computes the primitive vector associated to the input ray {c*vec: c>=0}.
    Very similar to the gcd function.

    This is equivalent to
        vec/gcd(vec)
    but just uses a rational representation.

    **Arguments:**
    - `vec`:       A vector defining the ray {c*vec: c>=0}
    - `max_denom`: Assert |di| <= max_denom

    **Returns:**
    The primitive vector along the ray.
    """
    # compute the rational representation
    rat   = [fractions.Fraction(v).limit_denominator(max_denom) for v in vec]
    numer = [r.numerator   for r in rat]
    denom = [r.denominator for r in rat]

    # get the LCM of the denominators
    l     = functools.reduce(math.lcm, denom)

    # get the integral vector and scale it to be primitive
    prim  = [n*(l//d) for n,d in zip(numer,denom)]
    gprime= functools.reduce(math.gcd, prim)

    return [x//gprime for x in prim]
