# from box import BoxN, Point

from dreal import And, Formula, Variable
from box import build_basic_box, BoxN
from objective import Polynomial, Rational, Term, to_term


# the input rational function is represented as two lists of terms
# each term is [coefficient, [exponent_vector]]
# for example

# f(x,y) = (2x^2y + 3y^2) / (x^2 + y)
# numerator_terms = [ [2, [2,1]], [3, [0,2]] ]
# denominator_terms = [ [1, [2,0]], [1, [0,1]] ]

# TODO: refactor these examples to use in main test suite!


def rational_objective_example(dim: int) -> Rational:
    numerator_terms: list[Term] = []
    denominator_terms: list[Term] = []

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


def initial_constraint_box(dim: int) -> BoxN:
    return build_basic_box(1, 10, dim)


def rational_objective_example_2(dim: int) -> Rational:
    """
    Construct numerator_terms and denominator_terms for the 3D rational
    function example generalized to n dimensions:

        numerator:
            sum_i x_i^4
          - 3 sum_i x_i^2
          + 0.5 sum_{i<j} x_i^2 x_j^2
          + 2

        denominator:
            1 + sum_i x_i^2

    NOTE: denominator must stay strictly > 0 on your domain.
    """

    numerator_terms: list[Term] = []
    denominator_terms: list[Term] = []

    # ----- Numerator: sum_i x_i^4 -----
    for i in range(dim):
        exps = [0] * dim
        exps[i] = 4
        numerator_terms.append(to_term(1.0, exps))

    # ----- Numerator: -3 * sum_i x_i^2 -----
    for i in range(dim):
        exps = [0] * dim
        exps[i] = 2
        numerator_terms.append(to_term(-3.0, exps))

    # ----- Numerator: 0.5 * sum_{i<j} x_i^2 x_j^2 -----
    for i in range(dim):
        for j in range(i + 1, dim):
            exps = [0] * dim
            exps[i] = 2
            exps[j] = 2
            numerator_terms.append(to_term(0.5, exps))

    # ----- Numerator: constant term = +2 -----
    numerator_terms.append(to_term(2.0, [0] * dim))

    # ----- Denominator: 1 + sum_i x_i^2 -----
    # constant term 1
    denominator_terms.append(to_term(1.0, [0] * dim))

    # sum_i x_i^2
    for i in range(dim):
        exps = [0] * dim
        exps[i] = 2
        denominator_terms.append(to_term(1.0, exps))

    return Rational(Polynomial(numerator_terms), Polynomial(denominator_terms))


def ball_constraint_example_2(vars: list[Variable]) -> Formula:
    dim = len(vars)

    center1 = [0.0] * dim
    center2 = [0.7] * dim
    center3 = [-0.7] * dim

    def sq_dist(xs, center):
        return sum((x - c) ** 2 for x, c in zip(xs, center))

    R1 = 4.0
    R2 = 4.0
    R3 = 4.0

    return And(
        sq_dist(vars, center1) <= R1 * R1,
        sq_dist(vars, center2) <= R2 * R2,
        sq_dist(vars, center3) <= R3 * R3,
    )


def initial_constraint_box_2(dim: int) -> BoxN:
    return build_basic_box(-10, 10, dim)
