from box_utils import BoxND
from dreal import *
def get_poly_terms(n):
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

    numerator_terms = []
    denominator_terms = []

    # ----- Numerator: sum_i x_i^4 -----
    for i in range(n):
        exps = [0] * n
        exps[i] = 4
        numerator_terms.append([1.0, exps])

    # ----- Numerator: -3 * sum_i x_i^2 -----
    for i in range(n):
        exps = [0] * n
        exps[i] = 2
        numerator_terms.append([-3.0, exps])

    # ----- Numerator: 0.5 * sum_{i<j} x_i^2 x_j^2 -----
    for i in range(n):
        for j in range(i + 1, n):
            exps = [0] * n
            exps[i] = 2
            exps[j] = 2
            numerator_terms.append([0.5, exps])

    # ----- Numerator: constant term = +2 -----
    numerator_terms.append([2.0, [0] * n])

    # ----- Denominator: 1 + sum_i x_i^2 -----
    # constant term 1
    denominator_terms.append([1.0, [0] * n])

    # sum_i x_i^2
    for i in range(n):
        exps = [0] * n
        exps[i] = 2
        denominator_terms.append([1.0, exps])

    return numerator_terms, denominator_terms

def f_constraint(*xs):
    n = len(xs)

    center1 = [0.0] * n
    center2 = [0.7] * n
    center3 = [-0.7] * n

    def sq_dist(xs, center):
        return sum((x - c) ** 2 for x, c in zip(xs, center))

    R1 = 4.0 
    R2 = 4.0
    R3 = 4.0 

    return And(
        sq_dist(xs, center1) <= R1 * R1,
        sq_dist(xs, center2) <= R2 * R2,
        sq_dist(xs, center3) <= R3 * R3,
    )

# def get_poly_terms(n):
#     numerator_terms = []
#     denominator_terms = []

#     # ----- Denominator: sum_i xi^2 -----
#     for i in range(n):
#         exps = [0] * n
#         exps[i] = 2
#         denominator_terms.append([1.0, exps])

#     # ----- Numerator: sum_i (xi^2 -2xi + 1) -----
#     constant_term = 0.0
#     for i in range(n):
#         # xi^2 term
#         exps2 = [0] * n
#         exps2[i] = 2
#         numerator_terms.append([1.0, exps2])

#         # -2 xi term
#         exps1 = [0] * n
#         exps1[i] = 1
#         numerator_terms.append([2.0, exps1])

#         # +1 constant for each dimension, accumulated later
#         constant_term += 1.0

#     # Add constant term n
#     numerator_terms.append([constant_term, [0] * n])

#     return numerator_terms, denominator_terms

# def f_constraint(*xs):
#     """
#     Example constraint generalized to n dimensions:

#     The original 3D version had three ball constraints.
#     Here we mimic a similar structure by placing several
#     2-radius balls at fixed positions in n-D space.

#     You can customize this if your real constraint is different.
#     """

#     # Build three centers similar to the 3D case:
#     # (3,3,3,...), (4,4,4,...), (3,4,4,...)
#     n = len(xs)

#     center1 = [3] * n
#     center2 = [4] * n
#     center3 = [3] + [4] * (n - 1)

#     def sq_dist(xs, center):
#         return sum((x - c) ** 2 for x, c in zip(xs, center))

#     return And(
#         sq_dist(xs, center1) <= 5 * 5,
#         sq_dist(xs, center2) <= 5 * 5,
#         sq_dist(xs, center3) <= 5 * 5,
#     )


def get_init_box(n: int):
    """
    Build an initial axis-aligned n-D box:

        [1,10] × [1,10] × ... × [1,10]

    Same as your original 3D version.
    """
    lows = [-10.0] * n
    highs = [10.0] * n

    return BoxND(lows, highs)