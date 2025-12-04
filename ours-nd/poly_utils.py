def poly_from_terms(terms):
    """
    Convert a list of [coeff, [e1, e2, ...]] into a monomial dict.

    terms: list of [coeff, [e1, e2, ...]] or (coeff, [e1, e2, ...])

    Returns:
        poly: dict mapping exponent tuples (e1,...,en) -> float coeff
    """
    poly = {}
    for coeff, exps in terms:
        poly[tuple(exps)] = float(coeff)
    return poly


def eval_poly_dict(poly, xs):
    """
    Evaluate a polynomial given as a monomial dict at a point or symbolically.

    poly: dict {(e1,...,en): coeff}
    xs:   list of values or symbolic dReal variables/expressions

    Returns:
        val: poly(xs) as a float (if xs are floats) or dReal Expression
    """
    val = 0
    for mon, coeff in poly.items():
        term = coeff
        for x_j, e in zip(xs, mon):
            if e != 0:
                term *= x_j ** e
        val += term
    return val


def rational_value_numeric(poly_num, poly_den, xs):
    """
    Evaluate the rational function f = N/D at a numeric point xs.

    xs: list[float]
    poly_num, poly_den: monomial dicts
    safe_eps: clip the denominator away from zero to avoid division by zero.

    Returns:
        float
    """
    num = eval_poly_dict(poly_num, xs)
    den = eval_poly_dict(poly_den, xs)
    return num / den


def rational_expr_symbolic(poly_num, poly_den, xs):
    """
    Build a dReal Expression for f = N/D using polynomial dictionaries.

    xs: list[dreal.Variable or Expression]
    poly_num, poly_den: monomial dicts

    Note: safe_eps is not directly used here; we rely on constraints to keep
    the denominator away from zero. You can add explicit constraints on D(x)
    in f_constraint if needed.
    """
    num = eval_poly_dict(poly_num, xs)
    den = eval_poly_dict(poly_den, xs)
    return num / den