import itertools
from math import comb

def bernstein_coefficients_nd(poly, box, degrees=None):
    """
    Correct multivariate Bernstein coefficient computation on a box.

    poly:   dict {(e1,...,en): coeff} in power basis w.r.t. x variables.
    box:    list/tuple of (a_j, b_j) intervals.
    degrees:
        Optional list of degrees per variable. If None, the maximal exponent
        in poly for each variable is used.

    Returns:
        (B, degrees), where B maps multi-indices (k1,...,kn) to Bernstein
        coefficients on the given box, and degrees is the degree per variable.
    """
    n = len(box)

    # 1) Determine degrees per variable
    if degrees is None:
        degrees = [0] * n
        for mon in poly:
            for j in range(n):
                degrees[j] = max(degrees[j], mon[j])

    # 2) Build a sparse tensor of power-basis coefficients in t after
    #    the change of variables x_j = a_j + (b_j - a_j) * t_j.
    #    P[k1,...,kn] is the coefficient of t_1^k1 * ... * t_n^kn.
    P = {}

    for mon, c in poly.items():
        # mon is the exponents in x
        # For each dimension j, expand x_j^{mon[j]} into t_j^{k_j}:
        #   x_j^{m_j} = sum_{k_j=0}^{m_j} comb(m_j, k_j) * a_j^{m_j-k_j} * (b_j-a_j)^{k_j} * t_j^{k_j}
        ranges = [range(mj + 1) for mj in mon]
        for ks in itertools.product(*ranges):
            coeff = c
            for j, kj in enumerate(ks):
                mj = mon[j]
                a_j, b_j = box[j]
                coeff *= comb(mj, kj) * (a_j ** (mj - kj)) * ((b_j - a_j) ** kj)
            P[ks] = P.get(ks, 0.0) + coeff

    # 3) Convert the power-basis tensor P into a Bernstein-basis tensor B
    #    along each axis using the 1D transformation:
    #
    #   For a univariate polynomial of degree d,
    #       p(t) = sum_{i=0}^d a_i t^i
    #   the Bernstein coefficients {b_k} of degree d on [0,1] satisfy:
    #       b_k = sum_{i=0}^k comb(k, i) / comb(d, i) * a_i.
    #
    B = dict(P)

    for axis in range(n):
        d = degrees[axis]
        newB = {}
        # Fix all other indices, and transform the coefficients along this axis.
        other_axes = [i for i in range(n) if i != axis]
        other_ranges = [range(degrees[i] + 1) for i in other_axes]

        for other_idx in itertools.product(*other_ranges):
            # Collect a_i for i = 0..d along this axis
            a_vec = []
            for i in range(d + 1):
                idx = []
                oi = 0
                for ax in range(n):
                    if ax == axis:
                        idx.append(i)
                    else:
                        idx.append(other_idx[oi])
                        oi += 1
                idx = tuple(idx)
                a_vec.append(B.get(idx, 0.0))

            # Power -> Bernstein conversion
            b_vec = [0.0] * (d + 1)
            for k in range(d + 1):
                s = 0.0
                for i in range(k + 1):
                    if comb(d, i) == 0:
                        continue
                    s += comb(k, i) / comb(d, i) * a_vec[i]
                b_vec[k] = s

            # Store back into newB
            for k, val in enumerate(b_vec):
                idx = []
                oi = 0
                for ax in range(n):
                    if ax == axis:
                        idx.append(k)
                    else:
                        idx.append(other_idx[oi])
                        oi += 1
                newB[tuple(idx)] = val

        B = newB

    return B, degrees

def bernstein_bounds(B):
    """Return min and max of a Bernstein coefficient tensor B (dict)."""
    vals = list(B.values())
    return min(vals), max(vals)

def bernstein_bounds_on_box(poly_num, poly_den, box):
    """
    Compute a crude lower/upper bound for f = N/D on a BoxND using Bernstein.

    This uses your existing bernstein_coefficients_nd and bernstein_bounds
    utilities. It ignores additional constraints; constraints are still
    handled by dReal in the main algorithm.

    Args:
        poly_num, poly_den: monomial dicts
        box: BoxND (lows, highs)

    Returns:
        (lb, ub): lower and upper bounds from Bernstein coefficients.
    """
    # Convert BoxND to a list of [a_i, b_i] intervals
    intervals = list(zip(box.lows, box.highs))
    n = len(intervals)

    # Determine global degrees for each variable
    degrees = [0] * n

    def update_degrees(poly):
        for mon in poly.keys():
            for j in range(n):
                degrees[j] = max(degrees[j], mon[j])

    update_degrees(poly_num)
    update_degrees(poly_den)

    # Compute Bernstein coefficients of numerator and denominator
    B_num, _ = bernstein_coefficients_nd(poly_num, intervals, degrees=degrees)
    B_den, _ = bernstein_coefficients_nd(poly_den, intervals, degrees=degrees)

    # Get min/max of Bernstein coefficients
    N_min, N_max = bernstein_bounds(B_num)
    D_min, D_max = bernstein_bounds(B_den)

    # Protect the denominator from being too small
    D_min = max(D_min, 1e-8)
    D_max = max(D_max, 1e-8)

    lb = N_min / D_max
    ub = N_max / D_min
    # print(intervals, N_min, N_max, D_min, D_max, lb)
    return lb, ub
