import itertools
from math import comb

# ============================================
# =============== PARAMETERS =================
# ============================================

# Stop subdividing when every box side is <= this
MIN_BOX_SIZE = 0.1

# Small positive epsilon used to avoid division by zero
SAFE_EPS = 1e-8
# Split heuristic options: "width", "objective_gradient", "constraint_gradient"
DEFAULT_SPLIT_STRATEGY = "objective_gradient"
DEFAULT_USE_BERNSTEIN = False
DEFAULT_USE_FULL_AFFINE = False
DEFAULT_EQ_TOL = 1e-6
DEFAULT_FULL_PRUNE_ANALYSIS = False


# ============================================
# =========== AFFINE ARITHMETIC ==============
# =======================================``=====

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

    def _pad_coeffs(self, other):
        m = max(len(self.coeffs), len(other.coeffs))
        a = self.coeffs + [0.0] * (m - len(self.coeffs))
        b = other.coeffs + [0.0] * (m - len(other.coeffs))
        return a, b

    def add(self, other):
        a, b = self._pad_coeffs(other)
        coeffs = [ai + bi for ai, bi in zip(a, b)]
        return AffineForm(self.mid + other.mid, coeffs)

    def sub(self, other):
        a, b = self._pad_coeffs(other)
        coeffs = [ai - bi for ai, bi in zip(a, b)]
        return AffineForm(self.mid - other.mid, coeffs)

    def mul_scalar(self, s):
        return AffineForm(self.mid * s, [c * s for c in self.coeffs])

    def mul(self, other):
        """
        Multiplication of two affine forms with one new noise symbol
        capturing the bilinear remainder.
        """
        a, b = self._pad_coeffs(other)
        c0 = self.mid * other.mid
        coeffs = [self.mid * bj + other.mid * ai for ai, bj in zip(a, b)]
        radius_self = sum(abs(ai) for ai in a)
        radius_other = sum(abs(bj) for bj in b)
        q = radius_self * radius_other
        coeffs.append(q)
        return AffineForm(c0, coeffs)

    def pow_int(self, e):
        if e == 0:
            return AffineForm(1.0, [])
        result = AffineForm(self.mid, list(self.coeffs))
        for _ in range(e - 1):
            result = result.mul(self)
        return result

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


def eval_poly_affine(poly, affine_vars):
    """
    Evaluate polynomial using affine arithmetic propagation.
    Returns range (min, max).
    """
    n = len(affine_vars)
    total = AffineForm(0.0, [])
    for mon, coeff in poly.items():
        term = AffineForm(1.0, [])
        for j in range(n):
            e = mon[j]
            if e == 0:
                continue
            term = term.mul(affine_vars[j].pow_int(e))
        term = term.mul_scalar(coeff)
        total = total.add(term)
    return total.range()


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

    def __init__(self, poly_num, poly_den, constraints=None, eq_constraints=None,
                 min_box_size=MIN_BOX_SIZE, safe_eps=SAFE_EPS,
                 split_strategy=DEFAULT_SPLIT_STRATEGY, max_depth=None,
                 use_bernstein=DEFAULT_USE_BERNSTEIN,
                 use_full_affine=DEFAULT_USE_FULL_AFFINE,
                 eq_tol=DEFAULT_EQ_TOL,
                 full_prune_analysis=DEFAULT_FULL_PRUNE_ANALYSIS):
        self.poly_num = poly_num
        self.poly_den = poly_den
        self.constraints = constraints or []
        self.eq_constraints = eq_constraints or []
        self.min_box_size = min_box_size
        self.safe_eps = safe_eps
        self.split_strategy = split_strategy
        self.max_depth = max_depth
        self.use_bernstein = use_bernstein
        self.use_full_affine = use_full_affine
        self.eq_tol = eq_tol
        self.full_prune_analysis = full_prune_analysis

        self.dim = None
        self.degrees = None  # global degrees per variable

        # Branch-and-bound global state
        self.best_ub = float('inf')
        self.best_box = None
        self.leaf_boxes = []  # list of (box, (lb, ub))
        self.midpoint_ub = float('inf')
        self.stats = {
            "aa_pruned": 0,
            "constraint_pruned": 0,
            "equality_pruned": 0,
            "bp_pruned": 0,
            "leaves": 0,
            "visited": 0,
        }

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
        for h in self.eq_constraints:
            update_with_poly(h)

        self.degrees = degs

    def minimize(self, box):
        """
        Entry point: run branch-and-bound starting from 'box'.
        Returns (best_ub, best_box, leaf_boxes).
        """
        self._init_degrees(box)
        self.stats = {k: 0 for k in self.stats}
        self.leaf_boxes = []
        # Seed the search with a midpoint evaluation to get a finite UB.
        mid_point = tuple((a + b) * 0.5 for a, b in box)
        mid_val = self._eval_objective(mid_point)
        self.best_ub = mid_val
        self.midpoint_ub = mid_val
        self.best_box = box

        self._subdivide(box, affine_vars=None, depth=0,
                        B_num=None, B_den=None, B_constraints=None, B_eq_constraints=None)
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

    def _eval_rational_grad(self, point):
        """Evaluate gradient of N/D at a concrete point."""
        n = len(point)
        num_val = 0.0
        den_val = 0.0
        grad_num = [0.0] * n
        grad_den = [0.0] * n

        for mon, coeff in self.poly_num.items():
            term = coeff
            for j, p_j in enumerate(point):
                term *= p_j ** mon[j]
            num_val += term
            for k in range(n):
                e_k = mon[k]
                if e_k == 0:
                    continue
                term_grad = coeff * e_k * (point[k] ** (e_k - 1))
                for j, p_j in enumerate(point):
                    if j == k:
                        continue
                    term_grad *= p_j ** mon[j]
                grad_num[k] += term_grad

        for mon, coeff in self.poly_den.items():
            term = coeff
            for j, p_j in enumerate(point):
                term *= p_j ** mon[j]
            den_val += term
            for k in range(n):
                e_k = mon[k]
                if e_k == 0:
                    continue
                term_grad = coeff * e_k * (point[k] ** (e_k - 1))
                for j, p_j in enumerate(point):
                    if j == k:
                        continue
                    term_grad *= p_j ** mon[j]
                grad_den[k] += term_grad

        den_val = max(den_val, self.safe_eps)
        grad_f = []
        for k in range(n):
            grad_f.append((grad_num[k] * den_val - num_val * grad_den[k]) / (den_val ** 2))
        return grad_f

    def _eval_poly_and_grad(self, poly, point):
        """Evaluate poly and its gradient at a point."""
        n = len(point)
        val = 0.0
        grad = [0.0] * n
        for mon, coeff in poly.items():
            term = coeff
            for j, p_j in enumerate(point):
                term *= p_j ** mon[j]
            val += term
            for k in range(n):
                e_k = mon[k]
                if e_k == 0:
                    continue
                term_grad = coeff * e_k * (point[k] ** (e_k - 1))
                for j, p_j in enumerate(point):
                    if j == k:
                        continue
                    term_grad *= p_j ** mon[j]
                grad[k] += term_grad
        return val, grad

    def _subdivide(self, box, affine_vars, depth,
                   B_num, B_den, B_constraints, B_eq_constraints):
        """
        Recursive branch-and-bound subdivision routine.
        """
        self.stats["visited"] += 1
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
        eval_fn = eval_poly_affine if self.use_full_affine else eval_poly_aa
        num_min, num_max = eval_fn(self.poly_num, affine_vars)
        den_min, den_max = eval_fn(self.poly_den, affine_vars)

        # Assume denominator > 0, but clip numerically
        den_min = max(den_min, self.safe_eps)
        den_max = max(den_max, self.safe_eps)

        aa_lb = num_min / den_max
        aa_ub = num_max / den_min

        # Keep best upper bound aggressive to help pruning.
        if aa_ub < self.best_ub:
            self.best_ub = aa_ub
            self.best_box = box

        # AA/equality/inequality pruning flags
        aa_prune = aa_lb >= self.best_ub
        constraint_prune = False
        equality_prune = False

        # Evaluate constraints/equalities either always (analysis) or only if AA didn't prune.
        if self.full_prune_analysis or not aa_prune:
            for g in self.constraints:
                g_min, g_max = eval_fn(g, affine_vars)
                if g_max < 0:
                    constraint_prune = True
            for h in self.eq_constraints:
                h_min, h_max = eval_fn(h, affine_vars)
                if h_min > self.eq_tol or h_max < -self.eq_tol:
                    equality_prune = True

        if self.full_prune_analysis:
            if aa_prune:
                self.stats["aa_pruned"] += 1
            if constraint_prune:
                self.stats["constraint_pruned"] += 1
            if equality_prune:
                self.stats["equality_pruned"] += 1
            if aa_prune or constraint_prune or equality_prune:
                return
        else:
            if aa_prune:
                self.stats["aa_pruned"] += 1
                return
            if constraint_prune:
                self.stats["constraint_pruned"] += 1
                return
            if equality_prune:
                self.stats["equality_pruned"] += 1
                return

        # -----------------------------
        # -----------------------------
        # 2. Optional Bernstein bounds (for tighter pruning and reuse)
        # -----------------------------
        if self.use_bernstein:
            if B_num is None:
                B_num, _ = bernstein_coefficients_nd(self.poly_num, box, degrees=self.degrees)
                B_den, _ = bernstein_coefficients_nd(self.poly_den, box, degrees=self.degrees)
                B_constraints = []
                for g in self.constraints:
                    B_g, _ = bernstein_coefficients_nd(g, box, degrees=self.degrees)
                    B_constraints.append(B_g)
                B_eq_constraints = []
                for h in self.eq_constraints:
                    B_h, _ = bernstein_coefficients_nd(h, box, degrees=self.degrees)
                    B_eq_constraints.append(B_h)

            B_num_min, B_num_max = bernstein_bounds(B_num)
            B_den_min, B_den_max = bernstein_bounds(B_den)
            B_den_min = max(B_den_min, self.safe_eps)
            B_den_max = max(B_den_max, self.safe_eps)

            bp_lb = B_num_min / B_den_max
            bp_ub = B_num_max / B_den_min

            # Update UB from BP as well
            if bp_ub < self.best_ub:
                self.best_ub = bp_ub
                self.best_box = box

            # BP-based pruning
            if bp_lb >= self.best_ub:
                self.stats["bp_pruned"] += 1
                return

            # Constraint BP pruning
            for B_g in (B_constraints or []):
                g_min, g_max = bernstein_bounds(B_g)
                if g_max < 0:
                    self.stats["constraint_pruned"] += 1
                    return
            # Equality BP pruning
            for B_h in (B_eq_constraints or []):
                h_min, h_max = bernstein_bounds(B_h)
                if h_min > self.eq_tol or h_max < -self.eq_tol:
                    self.stats["equality_pruned"] += 1
                    return
        else:
            bp_lb, bp_ub = aa_lb, aa_ub
            B_num = B_den = B_constraints = B_eq_constraints = None

        # -----------------------------
        # 2. Terminal box (leaf) check
        # -----------------------------
        widths = [b - a for a, b in box]
        max_width = max(widths)
        if max_width <= self.min_box_size or (
            self.max_depth is not None and depth >= self.max_depth
        ):
            # Hook: replace this midpoint evaluation with a dReal call if desired.
            mid_point = tuple((a + b) * 0.5 for a, b in box)
            mid_val = self._eval_objective(mid_point)
            local_lb = max(aa_lb, bp_lb)
            local_ub = min(aa_ub, bp_ub, mid_val)

            self.leaf_boxes.append((box, (local_lb, local_ub)))
            self.stats["leaves"] += 1

            if local_ub < self.best_ub:
                self.best_ub = local_ub
                self.best_box = box
            return

        # -----------------------------
        # 3. Choose split axis (heuristic)
        # -----------------------------
        splittable = [j for j, w in enumerate(widths) if w > self.min_box_size]
        if not splittable:
            splittable = list(range(self.dim))
        var_to_split = max(splittable, key=lambda j: widths[j])

        # Optional gradient-based strategies
        mid_point = tuple((a + b) * 0.5 for a, b in box)
        if self.split_strategy == "objective_gradient":
            grad_f = self._eval_rational_grad(mid_point)
            if any(abs(g) > 0 for g in grad_f):
                var_to_split = max(range(self.dim), key=lambda j: abs(grad_f[j]))
        elif self.split_strategy == "constraint_gradient" and self.constraints:
            # Focus on the most violated/closest-to-violated constraint
            best_idx = None
            best_val = float('inf')
            grad_choice = None
            for g in self.constraints:
                val, grad_g = self._eval_poly_and_grad(g, mid_point)
                if val < best_val:
                    best_val = val
                    grad_choice = grad_g
                    best_idx = 0  # placeholder to note we found one
            if grad_choice and any(abs(g) > 0 for g in grad_choice):
                var_to_split = max(range(self.dim), key=lambda j: abs(grad_choice[j]))
        # Ensure chosen axis is splittable; otherwise fall back.
        if widths[var_to_split] <= self.min_box_size and splittable:
            var_to_split = max(splittable, key=lambda j: widths[j])

        # -----------------------------
        # 4. Subdivide box (reuse Bernstein via de Casteljau)
        # -----------------------------
        a, b = box[var_to_split]
        mid = 0.5 * (a + b)
        box1 = list(box)
        box2 = list(box)
        box1[var_to_split] = (a, mid)
        box2[var_to_split] = (mid, b)

        left_B_num, right_B_num = None, None
        left_B_den, right_B_den = None, None
        left_B_constraints, right_B_constraints = None, None
        left_B_eq, right_B_eq = None, None
        if self.use_bernstein:
            left_B_num, right_B_num = de_casteljau_nd(B_num, var_to_split)
            left_B_den, right_B_den = de_casteljau_nd(B_den, var_to_split)

            if B_constraints:
                left_B_constraints = []
                right_B_constraints = []
                for B_g in B_constraints:
                    L, R = de_casteljau_nd(B_g, var_to_split)
                    left_B_constraints.append(L)
                    right_B_constraints.append(R)
            if B_eq_constraints:
                left_B_eq = []
                right_B_eq = []
                for B_h in B_eq_constraints:
                    L, R = de_casteljau_nd(B_h, var_to_split)
                    left_B_eq.append(L)
                    right_B_eq.append(R)

        # -----------------------------
        # 5. Recurse on children
        # -----------------------------
        self._subdivide(tuple(box1), affine_vars, depth + 1,
                        left_B_num, left_B_den, left_B_constraints, left_B_eq)
        self._subdivide(tuple(box2), affine_vars, depth + 1,
                        right_B_num, right_B_den, right_B_constraints, right_B_eq)


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
        # More nonlinear 10-variable problem on [0,1]^10 with higher-degree terms.
        "numerator": [
            # Cubic/quartic self terms
            [0.8, [3, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
            [1.0, [0, 3, 0, 0, 0, 0, 0, 0, 0, 0]],
            [1.2, [0, 0, 3, 0, 0, 0, 0, 0, 0, 0]],
            [1.4, [0, 0, 0, 4, 0, 0, 0, 0, 0, 0]],
            [1.6, [0, 0, 0, 0, 4, 0, 0, 0, 0, 0]],
            [1.8, [0, 0, 0, 0, 0, 3, 0, 0, 0, 0]],
            [2.0, [0, 0, 0, 0, 0, 0, 3, 0, 0, 0]],
            [2.2, [0, 0, 0, 0, 0, 0, 0, 4, 0, 0]],
            [2.4, [0, 0, 0, 0, 0, 0, 0, 0, 3, 0]],
            [2.6, [0, 0, 0, 0, 0, 0, 0, 0, 0, 4]],
            # Mixed higher-degree interactions
            [0.5, [1, 1, 1, 0, 0, 0, 0, 0, 0, 0]],
            [0.4, [0, 1, 0, 1, 1, 0, 0, 0, 0, 0]],
            [0.3, [0, 0, 1, 0, 1, 1, 0, 0, 0, 0]],
            [0.25, [0, 0, 0, 1, 0, 0, 1, 1, 0, 0]],
            [0.2, [0, 0, 0, 0, 1, 0, 0, 1, 1, 0]],
            [0.2, [0, 0, 0, 0, 0, 1, 0, 0, 1, 1]],
            [0.15, [0, 0, 0, 0, 0, 0, 1, 0, 2, 1]],
            [0.1, [0, 0, 0, 0, 0, 0, 0, 2, 1, 1]],
            # Mild linear bias
            [-0.2, [1, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
            [-0.2, [0, 1, 0, 0, 0, 0, 0, 0, 0, 0]],
        ],
        # Denominator: keep positive but nonlinear
        "denominator": [
            [1.0, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
            [0.3, [1, 1, 0, 0, 0, 0, 0, 0, 0, 0]],
            [0.25, [0, 0, 1, 1, 0, 0, 0, 0, 0, 0]],
            [0.2, [0, 0, 0, 0, 1, 1, 0, 0, 0, 0]],
            [0.2, [0, 0, 0, 0, 0, 0, 1, 1, 0, 0]],
            [0.15, [0, 0, 0, 0, 0, 0, 0, 0, 1, 1]],
            [0.1, [2, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
            [0.1, [0, 0, 0, 0, 0, 0, 0, 0, 0, 2]],
        ],
        # Box as list of [low, high] per variable (x0..x9)
        "box": [[0.0, 1.0]] * 10,
        # Polynomial inequality constraints g(x) >= 0
        "constraints": [
            # Nonlinear/multivariate constraints
            # x0*x1 + x2*x3 - 0.3 >= 0
            [[1.0, [1, 1, 0, 0, 0, 0, 0, 0, 0, 0]],
             [1.0, [0, 0, 1, 1, 0, 0, 0, 0, 0, 0]],
             [-0.3, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]],
            # x4^2 + x5^2 - 0.4 >= 0
            [[1.0, [0, 0, 0, 0, 2, 0, 0, 0, 0, 0]],
             [1.0, [0, 0, 0, 0, 0, 2, 0, 0, 0, 0]],
             [-0.4, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]],
            # 0.9 - x6*x7 - x8 >= 0  (upper bound involving product)
            [[0.9, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
             [-1.0, [0, 0, 0, 0, 0, 0, 1, 1, 0, 0]],
             [-1.0, [0, 0, 0, 0, 0, 0, 0, 0, 1, 0]]],
            # x9^3 - 0.2 >= 0
            [[1.0, [0, 0, 0, 0, 0, 0, 0, 0, 0, 3]],
             [-0.2, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]],
            # x0*x4 + x5*x6 - 0.35 >= 0
            [[1.0, [1, 0, 0, 0, 1, 0, 0, 0, 0, 0]],
             [1.0, [0, 0, 0, 0, 0, 1, 1, 0, 0, 0]],
             [-0.35, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]],
            # x2*x7 + x3*x8 - 0.25 >= 0
            [[1.0, [0, 0, 1, 0, 0, 0, 0, 1, 0, 0]],
             [1.0, [0, 0, 0, 1, 0, 0, 0, 0, 1, 0]],
             [-0.25, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]],
            # x1*x9 - x0*x8 - 0.1 >= 0
            [[1.0, [0, 1, 0, 0, 0, 0, 0, 0, 0, 1]],
             [-1.0, [1, 0, 0, 0, 0, 0, 0, 0, 1, 0]],
             [-0.1, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]],
            # x3*x4*x5 - 0.15 >= 0
            [[1.0, [0, 0, 0, 1, 1, 1, 0, 0, 0, 0]],
             [-0.15, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]],
            # (x6+x7+x8) - 1.2 >= 0
            [[1.0, [0, 0, 0, 0, 0, 0, 1, 0, 0, 0]],
             [1.0, [0, 0, 0, 0, 0, 0, 0, 1, 0, 0]],
             [1.0, [0, 0, 0, 0, 0, 0, 0, 0, 1, 0]],
             [-1.2, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]],
            # 1.5 - (x0 + x5 + x9) >= 0
            [[1.5, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
             [-1.0, [1, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
             [-1.0, [0, 0, 0, 0, 0, 1, 0, 0, 0, 0]],
             [-1.0, [0, 0, 0, 0, 0, 0, 0, 0, 0, 1]]],
            # x2^2 + x7^2 - 0.3 >= 0
            [[1.0, [0, 0, 2, 0, 0, 0, 0, 0, 0, 0]],
             [1.0, [0, 0, 0, 0, 0, 0, 0, 2, 0, 0]],
             [-0.3, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]],
            # x8*x9 + x4*x6 - 0.4 >= 0
            [[1.0, [0, 0, 0, 0, 0, 0, 0, 0, 1, 1]],
             [1.0, [0, 0, 0, 0, 1, 0, 1, 0, 0, 0]],
             [-0.4, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]],
        ],
        # Polynomial equality constraints h(x) = 0
        "eq_constraints": [
            # x0 + x1 + x2 - 1.0 = 0
            [[1.0, [1, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
             [1.0, [0, 1, 0, 0, 0, 0, 0, 0, 0, 0]],
             [1.0, [0, 0, 1, 0, 0, 0, 0, 0, 0, 0]],
             [-1.0, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]],
            # x4^2 - x5 = 0
            [[1.0, [0, 0, 0, 0, 2, 0, 0, 0, 0, 0]],
             [-1.0, [0, 0, 0, 0, 0, 1, 0, 0, 0, 0]]],
        ],
        # Stopping box width
        "min_box_size": 0.05,
        # Split strategy: "width", "objective_gradient", or "constraint_gradient"
        "split_strategy": DEFAULT_SPLIT_STRATEGY,
        # Optional max recursion depth (None for unlimited)
        "max_depth": 25,
        # Use Bernstein pruning? Set False to run AA-only for speed.
        "use_bernstein": False,
        # Use fuller affine propagation instead of interval-like AA?
        "use_full_affine": False,
        # Evaluate all prune tests each node (objective vs constraints vs equality) for stats
        "full_prune_analysis": False,
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
    eq_constraints = [poly_from_terms(p) for p in USER_PROBLEM.get("eq_constraints", [])]
    min_box_size = float(USER_PROBLEM["min_box_size"])
    split_strategy = USER_PROBLEM.get("split_strategy", DEFAULT_SPLIT_STRATEGY)
    max_depth = USER_PROBLEM.get("max_depth", None)
    use_bernstein = USER_PROBLEM.get("use_bernstein", DEFAULT_USE_BERNSTEIN)
    use_full_affine = USER_PROBLEM.get("use_full_affine", DEFAULT_USE_FULL_AFFINE)
    full_prune_analysis = USER_PROBLEM.get("full_prune_analysis", DEFAULT_FULL_PRUNE_ANALYSIS)

    opt = RationalBranchAndBoundMinimizer(
        poly_num,
        poly_den,
        constraints=constraints,
        eq_constraints=eq_constraints,
        min_box_size=min_box_size,
        safe_eps=SAFE_EPS,
        split_strategy=split_strategy,
        max_depth=max_depth,
        use_bernstein=use_bernstein,
        use_full_affine=use_full_affine,
        full_prune_analysis=full_prune_analysis,
    )

    start = time.time()
    best_ub, best_box, leaf_boxes = opt.minimize(box)
    elapsed = time.time() - start

    print("Best UB (approx. minimum) =", best_ub)
    print("Best box =", best_box)
    print("Number of leaf boxes explored =", len(leaf_boxes))
    print("Elapsed time (s) =", elapsed)
    print("Stats:", opt.stats)
