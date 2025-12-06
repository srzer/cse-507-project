import itertools
from math import comb

# ============================================
# =============== PARAMETERS =================
# ============================================

# Stop subdividing when every box side is <= this
MIN_BOX_SIZE = 0.01

# Small positive epsilon used to avoid division by zero
SAFE_EPS = 1e-8


# ============================================
# =========== AFFINE ARITHMETIC ==============
# ============================================

class AffineForm:
    """
    Very lightweight 'affine' form:
      x ≈ mid + sum_i coeffs[i] * ε_i,  with ε_i ∈ [-1,1].
    We are not doing full affine propagation here; we only use
    mid/coeffs to get a cheap range and rescale for sub-boxes.
    """
    def __init__(self, mid, coeffs=None):
        self.mid = mid
        self.coeffs = coeffs if coeffs is not None else []

    def range(self):
        r = sum(abs(c) for c in self.coeffs)
        return (self.mid - r, self.mid + r)

    def scale_shift(self, a_new, b_new):
        """
        Rescale this affine form from its current interval to [a_new, b_new].
        This is used when we subdivide a box and want to reuse parent AA info.
        """
        old_min, old_max = self.range()
        r_width = old_max - old_min
        r_new_width = b_new - a_new
        if r_width == 0:
            # Degenerate: collapse to midpoint of new interval
            alpha = 0.0
            beta = (a_new + b_new) * 0.5
        else:
            alpha = r_new_width / r_width
            beta = a_new - alpha * old_min

        new_coeffs = [c * alpha for c in self.coeffs]
        new_mid = self.mid * alpha + beta
        return AffineForm(new_mid, new_coeffs)


def eval_poly_aa(poly, affine_vars):
    """
    Very simple range evaluation of a polynomial using only the
    ranges of each variable (this is effectively interval arithmetic).
    poly: dict {monomial_tuple: coeff}
    affine_vars: list of AffineForm, one per variable
    Returns: (min_val, max_val) over the current box.
    """
    n = len(affine_vars)
    min_val = 0.0
    max_val = 0.0

    for mon, coeff in poly.items():
        # Start with the coefficient
        term_min = coeff
        term_max = coeff
        for j in range(n):
            a_j, b_j = affine_vars[j].range()
            e = mon[j]

            # Evaluate x_j^e on [a_j, b_j]
            if e == 0:
                # x^0 = 1
                p_min = p_max = 1.0
            elif e > 0:
                # monotone if interval doesn't cross 0;
                # otherwise, min at 0, max at max(|a_j|,|b_j|)
                if a_j >= 0:
                    p_min = a_j ** e
                    p_max = b_j ** e
                elif b_j <= 0:
                    p_min = b_j ** e
                    p_max = a_j ** e
                else:
                    # crosses zero
                    p_min = 0.0
                    p_max = max(abs(a_j), abs(b_j)) ** e
            else:
                # negative exponent not expected for polynomials
                raise ValueError("Negative exponent in polynomial monomial")

            # Multiply this term's bounds by p_min..p_max
            candidates = [
                term_min * p_min,
                term_min * p_max,
                term_max * p_min,
                term_max * p_max,
            ]
            term_min = min(candidates)
            term_max = max(candidates)

        min_val += term_min
        max_val += term_max

    return (min_val, max_val)


# ============================================
# ============ BERNSTEIN HELPERS =============
# ============================================

def bernstein_bounds(B):
    """Return min and max of a Bernstein coefficient tensor B (dict)."""
    vals = list(B.values())
    return min(vals), max(vals)


def bernstein_coefficients_nd(poly, box, degrees=None):
    """
    Compute multivariate Bernstein coefficients of a polynomial on a box.
    poly: {monomial_tuple: coeff}
    box: [(a1,b1), ..., (an,bn)]
    degrees: list of degrees per variable (same for all polynomials) or None
    Returns: (B, degrees), where B maps multi-indices -> coefficient
    """
    n = len(box)

    if degrees is None:
        degrees = [0] * n
        for mon in poly:
            for j in range(n):
                degrees[j] = max(degrees[j], mon[j])

    all_indices = list(itertools.product(*[range(d + 1) for d in degrees]))
    B = {idx: 0.0 for idx in all_indices}

    # Naive conversion from power basis to Bernstein on [a,b]
    # (this is not the most efficient formula, but OK for now)
    for mon, coeff in poly.items():
        for idx in all_indices:
            term = coeff
            for j in range(n):
                k = idx[j]
                if k <= mon[j]:
                    a_j, b_j = box[j]
                    term *= comb(mon[j], k) * ((b_j - a_j) ** k) * (a_j ** (mon[j] - k))
                else:
                    term = 0.0
                    break
            B[idx] += term

    return B, degrees


def de_casteljau_1d(coeffs):
    """
    1D de Casteljau subdivision of Bernstein coefficients at t = 1/2.
    Returns (left_coeffs, right_coeffs).
    """
    n = len(coeffs) - 1
    left = coeffs.copy()
    right = coeffs.copy()

    # Left triangle
    for r in range(1, n + 1):
        for i in range(n - r + 1):
            left[i] = 0.5 * left[i] + 0.5 * left[i + 1]

    # Right triangle
    for r in range(1, n + 1):
        for i in range(n - 1, r - 2, -1):
            right[i] = 0.5 * right[i] + 0.5 * right[i - 1]

    return left, right


def de_casteljau_nd(B, axis):
    """
    nD de Casteljau subdivision of a Bernstein tensor B along a chosen axis.
    Returns (left_B, right_B) for the two sub-boxes split at midpoint.
    """
    n_vars = len(next(iter(B.keys())))
    degrees = [max(idx[i] for idx in B.keys()) for i in range(n_vars)]

    left_B = {}
    right_B = {}

    other_axes = [i for i in range(n_vars) if i != axis]
    fixed_indices_list = list(
        itertools.product(*[range(degrees[i] + 1) for i in other_axes])
    )

    for fixed in fixed_indices_list:
        # Extract 1D coefficients along 'axis'
        coeffs_1d = []
        for k in range(degrees[axis] + 1):
            idx = []
            j = 0
            for ax in range(n_vars):
                if ax == axis:
                    idx.append(k)
                else:
                    idx.append(fixed[j])
                    j += 1
            coeffs_1d.append(B[tuple(idx)])

        left_1d, right_1d = de_casteljau_1d(coeffs_1d)

        # Insert left
        for k, val in enumerate(left_1d):
            idx = []
            j = 0
            for ax in range(n_vars):
                if ax == axis:
                    idx.append(k)
                else:
                    idx.append(fixed[j])
                    j += 1
            left_B[tuple(idx)] = val

        # Insert right
        for k, val in enumerate(right_1d):
            idx = []
            j = 0
            for ax in range(n_vars):
                if ax == axis:
                    idx.append(k)
                else:
                    idx.append(fixed[j])
                    j += 1
            right_B[tuple(idx)] = val

    return left_B, right_B


def variable_ratio_spread(B_num, B_den, safe_eps=SAFE_EPS):
    """
    For each variable, estimate how much the rational function N/D can vary
    along that variable, using Bernstein coefficients.

    Returns: list 'contribs' of length n_vars, one spread per variable.
    """
    n_vars = len(next(iter(B_num.keys())))
    contribs = []
    degrees = [max(idx[i] for idx in B_num.keys()) for i in range(n_vars)]

    for var_idx in range(n_vars):
        other_axes = [i for i in range(n_vars) if i != var_idx]
        fixed_indices_list = list(
            itertools.product(*[range(degrees[i] + 1) for i in other_axes])
        )
        max_spread = 0.0

        for fixed in fixed_indices_list:
            N_coeffs = []
            D_coeffs = []
            for k in range(degrees[var_idx] + 1):
                idx = []
                j = 0
                for ax in range(n_vars):
                    if ax == var_idx:
                        idx.append(k)
                    else:
                        idx.append(fixed[j])
                        j += 1
                N_coeffs.append(B_num[tuple(idx)])
                D_coeffs.append(B_den[tuple(idx)])

            minN, maxN = min(N_coeffs), max(N_coeffs)
            # Clip denominator coefficients away from zero
            D_clipped = [max(abs(d), safe_eps) for d in D_coeffs]
            minD, maxD = min(D_clipped), max(D_clipped)

            spread = maxN / minD - minN / maxD
            max_spread = max(max_spread, spread)

        contribs.append(max_spread)

    return contribs


def variable_spread(B):
    """
    For a single polynomial's Bernstein tensor B, compute how much it can vary
    along each variable, as a crude spread measure.
    Returns a list of length n_vars with spreads per variable.
    """
    n_vars = len(next(iter(B.keys())))
    degrees = [max(idx[i] for idx in B.keys()) for i in range(n_vars)]
    spreads = []

    for var_idx in range(n_vars):
        other_axes = [i for i in range(n_vars) if i != var_idx]
        fixed_indices_list = list(
            itertools.product(*[range(degrees[i] + 1) for i in other_axes])
        )
        max_spread = 0.0

        for fixed in fixed_indices_list:
            coeffs = []
            for k in range(degrees[var_idx] + 1):
                idx = []
                j = 0
                for ax in range(n_vars):
                    if ax == var_idx:
                        idx.append(k)
                    else:
                        idx.append(fixed[j])
                        j += 1
                coeffs.append(B[tuple(idx)])
            max_spread = max(max_spread, max(coeffs) - min(coeffs))

        spreads.append(max_spread)

    return spreads


# ============================================
# ====== BRANCH-AND-BOUND MINIMIZATION =======
# ============================================

class RationalBranchAndBoundMinimizer:
    """
    Minimize a rational function f = N/D over a box and semialgebraic constraints
    g_i(x) >= 0 using:
      - AA for cheap first pruning
      - Bernstein polynomials for tighter pruning
      - Branch-and-bound on a global best upper bound.
    """

    def __init__(self, poly_num, poly_den, constraints=None,
                 min_box_size=MIN_BOX_SIZE, safe_eps=SAFE_EPS):
        self.poly_num = poly_num
        self.poly_den = poly_den
        self.constraints = constraints or []
        self.min_box_size = min_box_size
        self.safe_eps = safe_eps

        self.dim = None
        self.degrees = None  # global degrees per variable

        # Branch-and-bound global state
        self.best_ub = float('inf')
        self.best_box = None
        self.leaf_boxes = []  # list of (box, (lb, ub))
        self.midpoint_ub = float('inf')

    def _init_degrees(self, box):
        """Compute a global degrees vector covering numerator, denominator, constraints."""
        n = len(box)
        self.dim = n
        degs = [0] * n

        def update_with_poly(poly):
            for mon in poly:
                for j in range(n):
                    degs[j] = max(degs[j], mon[j])

        update_with_poly(self.poly_num)
        update_with_poly(self.poly_den)
        for g in self.constraints:
            update_with_poly(g)

        self.degrees = degs

    def minimize(self, box):
        """
        Entry point: run branch-and-bound starting from 'box'.
        Returns (best_ub, best_box, leaf_boxes).
        """
        self._init_degrees(box)
        # Seed the search with a midpoint evaluation to get a finite UB.
        mid_point = tuple((a + b) * 0.5 for a, b in box)
        mid_val = self._eval_objective(mid_point)
        self.best_ub = mid_val
        self.midpoint_ub = mid_val
        self.best_box = box

        self._subdivide(box, affine_vars=None)
        return self.best_ub, self.best_box, self.leaf_boxes

    def _eval_objective(self, point):
        """Evaluate the rational objective N/D at a concrete point."""
        num_val = 0.0
        den_val = 0.0
        for mon, coeff in self.poly_num.items():
            term = coeff
            for j, p_j in enumerate(point):
                term *= p_j ** mon[j]
            num_val += term
        for mon, coeff in self.poly_den.items():
            term = coeff
            for j, p_j in enumerate(point):
                term *= p_j ** mon[j]
            den_val += term
        den_val = max(den_val, self.safe_eps)
        return num_val / den_val

    def _subdivide(self, box, affine_vars):
        """
        Recursive branch-and-bound subdivision routine.
        """

        # -----------------------------
        # 1. Affine Arithmetic (AA)
        # -----------------------------
        if affine_vars is None:
            affine_vars = [
                AffineForm((a + b) * 0.5, [(b - a) * 0.5]) for a, b in box
            ]
        else:
            affine_vars = [
                v.scale_shift(a, b) for v, (a, b) in zip(affine_vars, box)
            ]

        # AA bounds for numerator and denominator
        num_min, num_max = eval_poly_aa(self.poly_num, affine_vars)
        den_min, den_max = eval_poly_aa(self.poly_den, affine_vars)

        # Assume denominator > 0, but clip numerically
        den_min = max(den_min, self.safe_eps)
        den_max = max(den_max, self.safe_eps)

        aa_lb = num_min / den_max
        aa_ub = num_max / den_min

        # Keep best upper bound aggressive to help pruning.
        if aa_ub < self.best_ub:
            self.best_ub = aa_ub
            self.best_box = box

        # AA-based branch-and-bound pruning (minimization):
        # if the AA lower bound is already worse than our best UB, prune.
        if aa_lb >= self.best_ub:
            return

        # AA-based constraint pruning
        for g in self.constraints:
            g_min, g_max = eval_poly_aa(g, affine_vars)
            if g_max < 0:
                # Entire box violates g(x) >= 0
                return

        # -----------------------------
        # 2. Terminal box (leaf) check
        # -----------------------------
        if all((b - a) <= self.min_box_size for a, b in box):
            # Hook: replace this midpoint evaluation with a dReal call if desired.
            mid_point = tuple((a + b) * 0.5 for a, b in box)
            mid_val = self._eval_objective(mid_point)
            local_lb = aa_lb
            local_ub = min(aa_ub, mid_val)

            self.leaf_boxes.append((box, (local_lb, local_ub)))

            if local_ub < self.best_ub:
                self.best_ub = local_ub
                self.best_box = box
            return

        # -----------------------------
        # 3. Choose split axis (heuristic)
        # -----------------------------
        # Split along the widest interval to ensure progress.
        widths = [b - a for a, b in box]
        var_to_split = max(range(self.dim), key=lambda j: widths[j])

        # -----------------------------
        # 4. Subdivide box
        # -----------------------------
        a, b = box[var_to_split]
        mid = 0.5 * (a + b)
        box1 = list(box)
        box2 = list(box)
        box1[var_to_split] = (a, mid)
        box2[var_to_split] = (mid, b)

        # -----------------------------
        # 5. Recurse on children
        # -----------------------------
        self._subdivide(tuple(box1), affine_vars)
        self._subdivide(tuple(box2), affine_vars)


# ============================================
# ================ EXAMPLE ===================
# ============================================

if __name__ == "__main__":
    import json
    import time

    # -----------------------------------------------
    # Edit this block to set your problem by hand.
    # -----------------------------------------------
    USER_PROBLEM = {
        # Numerator: list of [coeff, [e1, e2, ...]]
        "numerator": [
            [2.0, [2, 0]],
            [1.0, [1, 1]],
            [-3.0, [0, 1]],
        ],
        # Denominator: list of [coeff, [e1, e2, ...]]
        "denominator": [
            [1.0, [1, 0]],
            [2.0, [0, 0]],
        ],
        # Box as list of [low, high] per variable
        "box": [
            [0.0, 1.0],  # x
            [1.0, 2.0],  # y
        ],
        # Polynomial inequality constraints g(x) >= 0
        "constraints": [
            [[1.0, [0, 1]], [-1.0, [1, 0]]],  # y - x >= 0
        ],
        # Stopping box width
        "min_box_size": MIN_BOX_SIZE,
    }

    def poly_from_terms(terms):
        poly = {}
        for coeff, exps in terms:
            poly[tuple(exps)] = float(coeff)
        return poly

    poly_num = poly_from_terms(USER_PROBLEM["numerator"])
    poly_den = poly_from_terms(USER_PROBLEM["denominator"])
    box = tuple((float(a), float(b)) for a, b in USER_PROBLEM["box"])
    constraints = [poly_from_terms(p) for p in USER_PROBLEM["constraints"]]
    min_box_size = float(USER_PROBLEM["min_box_size"])

    opt = RationalBranchAndBoundMinimizer(
        poly_num,
        poly_den,
        constraints=constraints,
        min_box_size=min_box_size,
        safe_eps=SAFE_EPS,
    )

    start = time.time()
    best_ub, best_box, leaf_boxes = opt.minimize(box)
    elapsed = time.time() - start

    print("Best UB (approx. minimum) =", best_ub)
    print("Best box =", best_box)
    print("Number of leaf boxes explored =", len(leaf_boxes))
    print("Elapsed time (s) =", elapsed)
