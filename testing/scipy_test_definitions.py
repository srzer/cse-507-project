from dataclasses import dataclass
from typing import Callable, List

from box import BoxN, Point
from objective import Rational, PolyBuilder, make_vars


@dataclass
class NumericProblem:
    name: str
    rational_obj: Rational
    func_lambda: Callable
    constraint_maker: Callable
    float_constraint: Callable
    init_box: BoxN
    description: str


# ==========================================
# Helpers for Intersection Tests
# ==========================================


def _get_3sphere_constraints(x, y, z):
    return [
        (x - 3) ** 2 + (y - 3) ** 2 + (z - 3) ** 2 <= 4,
        (x - 4) ** 2 + (y - 4) ** 2 + (z - 4) ** 2 <= 4,
        (x - 3) ** 2 + (y - 4) ** 2 + (z - 4) ** 2 <= 4,
    ]


def _check_3sphere_constraints(x, y, z):
    return (
        ((x - 3) ** 2 + (y - 3) ** 2 + (z - 3) ** 2 <= 4)
        and ((x - 4) ** 2 + (y - 4) ** 2 + (z - 4) ** 2 <= 4)
        and ((x - 3) ** 2 + (y - 4) ** 2 + (z - 4) ** 2 <= 4)
    )


# ==========================================
# Standardized Tests
# ==========================================


def test_sanity_poly():
    x, y, z = make_vars(3)
    r_obj = (x**2 + y**2 + z**2) / 1.0

    def constraint_maker(x, y, z, A, O):
        return A(x >= -2, x <= 2, y >= -2, y <= 2, z >= -2, z <= 2)

    def float_constraint(x, y, z):
        return -2 <= x <= 2 and -2 <= y <= 2 and -2 <= z <= 2

    return NumericProblem(
        "Sanity Poly",
        r_obj,
        lambda x, y, z: x**2 + y**2 + z**2,
        constraint_maker,
        float_constraint,
        BoxN(Point((-2, -2, -2)), Point((2, 2, 2))),
        "Simple x^2",
    )


def test_sanity_rational():
    x, y, z = make_vars(3)
    r_obj = (x + 1.0) / (y + 1.0)

    def constraint_maker(x, y, z, A, O):
        return A(x >= 0, x <= 1, y >= 0, y <= 1, z == 0)

    def float_constraint(x, y, z):
        return 0 <= x <= 1 and 0 <= y <= 1 and z == 0

    return NumericProblem(
        "Sanity Rational",
        r_obj,
        lambda x, y, z: (x + 1) / (y + 1),
        constraint_maker,
        float_constraint,
        BoxN(Point((0, 0, -1)), Point((1, 1, 1))),
        "Simple x/y",
    )


def test_rational_bowl():
    x, y, z = make_vars(3)
    r_obj = (x**2 + y**2 + z**2) / (x + y + z)

    def constraint_maker(x, y, z, A, O):
        return A(x**2 + y**2 + z**2 <= 9)

    def float_constraint(x, y, z):
        return x**2 + y**2 + z**2 <= 9

    return NumericProblem(
        "Rational Bowl",
        r_obj,
        lambda x, y, z: (x**2 + y**2 + z**2) / (x + y + z),
        constraint_maker,
        float_constraint,
        BoxN(Point((0.1, 0.1, 0.1)), Point((5, 5, 5))),
        "Simple rational",
    )


def test_himmelblau_ratio():
    x, y, z = make_vars(3)
    num = (x**2 + y - 11) ** 2 + (x + y**2 - 7) ** 2 + z**2
    den = 1.0 + (x**2 + y**2) * 0.01
    r_obj = num / den

    def constraint_maker(x, y, z, A, O):
        return A(x**2 + y**2 <= 50)

    def float_constraint(x, y, z):
        return x**2 + y**2 <= 50

    return NumericProblem(
        "Himmelblau Ratio",
        r_obj,
        lambda x, y, z: ((x**2 + y - 11) ** 2 + (x + y**2 - 7) ** 2 + z**2)
        / (1 + 0.01 * (x**2 + y**2)),
        constraint_maker,
        float_constraint,
        BoxN(Point((-5, -5, -1)), Point((5, 5, 1))),
        "Multimodal",
    )


def test_split_islands():
    x, y, z = make_vars(3)
    r_obj = (x + y + z) / 1.0

    def constraint_maker(x, y, z, A, O):
        return O((x - 2) ** 2 + y**2 + z**2 <= 0.25, (x + 2) ** 2 + y**2 + z**2 <= 0.25)

    def float_constraint(x, y, z):
        return (x - 2) ** 2 + y**2 + z**2 <= 0.25 or (x + 2) ** 2 + y**2 + z**2 <= 0.25

    return NumericProblem(
        "Split Islands",
        r_obj,
        lambda x, y, z: x + y + z,
        constraint_maker,
        float_constraint,
        BoxN(Point((-5, -2, -2)), Point((5, 2, 2))),
        "Disconnected",
    )


def test_singularity_edge():
    x, y, z = make_vars(3)
    num = y * z + x * z + x * y
    den = x * y * z
    r_obj = num / den

    def constraint_maker(x, y, z, A, O):
        return A(x**2 + y**2 + z**2 >= 1)

    def float_constraint(x, y, z):
        return x**2 + y**2 + z**2 >= 1

    return NumericProblem(
        "Singularity Edge",
        r_obj,
        lambda x, y, z: 1 / x + 1 / y + 1 / z,
        constraint_maker,
        float_constraint,
        BoxN(Point((0.1, 0.1, 0.1)), Point((3, 3, 3))),
        "Non-convex hole",
    )


def test_pole_avoidance():
    x, y, z = make_vars(3)
    r_obj = PolyBuilder.const(1.0, 3) / (x + y + z - 2.5)

    def constraint_maker(x, y, z, A, O):
        return A(x >= 1, y >= 1, z >= 1)

    def float_constraint(x, y, z):
        return x >= 1 and y >= 1 and z >= 1

    return NumericProblem(
        "Pole Avoidance",
        r_obj,
        lambda x, y, z: 1 / (x + y + z - 2.5),
        constraint_maker,
        float_constraint,
        BoxN(Point((1, 1, 1)), Point((2, 2, 2))),
        "Singularity",
    )


def test_rational_valley():
    x, y, z = make_vars(3)
    r_obj = (x**2 + y**2 + z**2 + 1) / (x * y * z + 1)

    def constraint_maker(x, y, z, A, O):
        return A(x <= 2, y <= 2, z <= 2)

    def float_constraint(x, y, z):
        return x <= 2 and y <= 2 and z <= 2

    return NumericProblem(
        "Rational Valley",
        r_obj,
        lambda x, y, z: (x**2 + y**2 + z**2 + 1) / (x * y * z + 1),
        constraint_maker,
        float_constraint,
        BoxN(Point((0.1, 0.1, 0.1)), Point((2, 2, 2))),
        "Valley",
    )


def test_positive_islands():
    x, y, z = make_vars(3)
    r_obj = (x**2 + 1) / 1.0

    def constraint_maker(x, y, z, A, O):
        return O((x - 2) ** 2 + y**2 + z**2 <= 0.25, (x + 2) ** 2 + y**2 + z**2 <= 0.25)

    def float_constraint(x, y, z):
        return (x - 2) ** 2 + y**2 + z**2 <= 0.25 or (x + 2) ** 2 + y**2 + z**2 <= 0.25

    return NumericProblem(
        "Positive Islands",
        r_obj,
        lambda x, y, z: x**2 + 1,
        constraint_maker,
        float_constraint,
        BoxN(Point((-5, -5, -5)), Point((5, 5, 5))),
        "Strict positive disconnected",
    )


def test_sparse_intersection():
    x, y, z = make_vars(3)
    num = (x - 1) ** 2 + (y - 1) ** 2 + (z - 1) ** 2
    den = x + y + z
    r_obj = num / den

    def constraint_maker(x, y, z, A, O):
        return A(*_get_3sphere_constraints(x, y, z))

    return NumericProblem(
        "Sparse Intersection",
        r_obj,
        lambda x, y, z: ((x - 1) ** 2 + (y - 1) ** 2 + (z - 1) ** 2) / (x + y + z),
        constraint_maker,
        _check_3sphere_constraints,
        BoxN(Point((1, 1, 1)), Point((10, 10, 10))),
        "Hard Intersection",
    )


def all_numeric_tests() -> List[NumericProblem]:
    """Returns a list of all legacy test problems."""
    return [
        test_sanity_poly(),
        test_sanity_rational(),
        test_rational_bowl(),
        test_himmelblau_ratio(),
        test_split_islands(),
        test_singularity_edge(),
        test_pole_avoidance(),
        test_rational_valley(),
        test_positive_islands(),
        test_sparse_intersection(),
    ]
