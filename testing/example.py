# from box import BoxN, Point

from dreal import And, Formula, Variable
from poly.type import Polynomial, Rational, to_term


# the input rational function is represented as two lists of terms
# each term is [coefficient, [exponent_vector]]
# for example

# f(x,y) = (2x^2y + 3y^2) / (x^2 + y)
# numerator_terms = [ [2, [2,1]], [3, [0,2]] ]
# denominator_terms = [ [1, [2,0]], [1, [0,1]] ]


def rational_objective_example(dim: int) -> Rational:
    numerator_terms = []
    denominator_terms = []

    # ----- Denominator: sum_i xi^2 -----
    for i in range(dim):
        exps = [0] * dim
        exps[i] = 2
        denominator_terms.append(to_term(1.0, exps))

    # ----- Numerator: sum_i (xi^2 -2xi + 1) -----
    constant_term = 0.0
    for i in range(dim):
        # xi^2 term
        exps2 = [0] * dim
        exps2[i] = 2
        numerator_terms.append(to_term(1.0, exps2))

        # -2 xi term
        exps1 = [0] * dim
        exps1[i] = 1
        numerator_terms.append(to_term(-2.0, exps1))

        # +1 constant for each dimension, accumulated later
        constant_term += 1.0

    # add constant term n
    numerator_terms.append(to_term(constant_term, [0] * dim))

    return Rational(Polynomial(numerator_terms), Polynomial(denominator_terms))


def ball_constraint_example(vars: list[Variable]) -> Formula:
    """
    Example constraint generalized to n dimensions:

    The original 3D version had three ball constraints.
    Here we mimic a similar structure by placing several
    2-radius balls at fixed positions in n-D space.

    You can customize this if your real constraint is different.
    """

    # Build three centers similar to the 3D case:
    # (3,3,3,...), (4,4,4,...), (3,4,4,...)
    dim = len(vars)

    center1 = [3] * dim
    center2 = [4] * dim
    center3 = [3] + [4] * (dim - 1)

    def sq_dist(xs, center):
        return sum((x - c) ** 2 for x, c in zip(xs, center))

    return And(
        sq_dist(vars, center1) <= 2 * 2,
        sq_dist(vars, center2) <= 2 * 2,
        sq_dist(vars, center3) <= 2 * 2,
    )
