from typing import List

from dreal import And, Or, Formula, Variable
from box import BoxN, Point, build_basic_box
from objective import (
    Rational,
    Polynomial,
    to_term,
)

from .test_type import Problem


def _get_3sphere_constraints_formula(variables: List[Variable]) -> Formula:
    x, y, z = variables[:3]  # Assuming dim >= 3
    return And(
        (x - 3) ** 2 + (y - 3) ** 2 + (z - 3) ** 2 <= 4,
        (x - 4) ** 2 + (y - 4) ** 2 + (z - 4) ** 2 <= 4,
        (x - 3) ** 2 + (y - 4) ** 2 + (z - 4) ** 2 <= 4,
    )


def test_sanity_poly(dim: int) -> Problem:
    """
    Objective: x^2 + y^2 + z^2
    Constraints: -2 <= x,y,z <= 2
    """
    variables = [Variable(f"x{i}") for i in range(dim)]
    x, y, z = variables[0], variables[1], variables[2]

    # Objective: x^2 + y^2 + z^2
    numerator_poly = Polynomial(
        [
            to_term(1.0, [2, 0, 0]),
            to_term(1.0, [0, 2, 0]),
            to_term(1.0, [0, 0, 2]),
        ]
    )
    denominator_poly = Polynomial([to_term(1.0, [0, 0, 0])])
    objective = Rational(numerator_poly, denominator_poly)

    # Constraints: -2 <= x,y,z <= 2
    constraints = And(
        x >= -2,
        x <= 2,
        y >= -2,
        y <= 2,
        z >= -2,
        z <= 2,
    )

    initial_box = BoxN(Point([-2.0] * dim), Point([2.0] * dim))

    return Problem(
        name="Sanity Poly",
        dim=dim,
        objective=objective,
        constraints=constraints,
        initial_box=initial_box,
        variables=variables,
    )


def test_sanity_rational(dim: int) -> Problem:
    variables = [Variable(f"x{i}") for i in range(dim)]
    x, y, z = variables[0], variables[1], variables[2]

    # r_obj = (x + 1.0) / (y + 1.0)
    numerator_poly = Polynomial([to_term(1.0, [1, 0, 0]), to_term(1.0, [0, 0, 0])])
    denominator_poly = Polynomial([to_term(1.0, [0, 1, 0]), to_term(1.0, [0, 0, 0])])
    objective = Rational(numerator_poly, denominator_poly)

    # constraint_maker(x, y, z, A, O): return A(x >= 0, x <= 1, y >= 0, y <= 1, z == 0)
    constraints = And(x >= 0, x <= 1, y >= 0, y <= 1, z == 0)

    # BoxN(Point((0, 0, -1)), Point((1, 1, 1)))
    initial_box = BoxN(Point([0.0] * dim), Point([1.0] * dim))
    if dim > 2:
        initial_box = BoxN(
            Point([0.0, 0.0, -1.0] + [0.0] * (dim - 3)),
            Point([1.0, 1.0, 1.0] + [1.0] * (dim - 3)),
        )

    return Problem(
        "Sanity Rational",
        dim,
        objective,
        constraints,
        initial_box,
        variables,
    )


def test_rational_bowl(dim: int) -> Problem:
    variables = [Variable(f"x{i}") for i in range(dim)]
    x, y, z = variables[0], variables[1], variables[2]

    # r_obj = (x**2 + y**2 + z**2) / (x + y + z)
    numerator_poly = Polynomial(
        [to_term(1.0, [2, 0, 0]), to_term(1.0, [0, 2, 0]), to_term(1.0, [0, 0, 2])]
    )
    denominator_poly = Polynomial(
        [to_term(1.0, [1, 0, 0]), to_term(1.0, [0, 1, 0]), to_term(1.0, [0, 0, 1])]
    )
    objective = Rational(numerator_poly, denominator_poly)

    # constraint_maker(x, y, z, A, O): return A(x**2 + y**2 + z**2 <= 9)
    constraints = And(x**2 + y**2 + z**2 <= 9)

    # BoxN(Point((0.1, 0.1, 0.1)), Point((5, 5, 5)))
    initial_box = BoxN(Point([0.1] * dim), Point([5.0] * dim))

    return Problem(
        "Rational Bowl",
        dim,
        objective,
        constraints,
        initial_box,
        variables,
    )


def test_himmelblau_ratio(dim: int) -> Problem:
    variables = [Variable(f"x{i}") for i in range(dim)]
    x, y, z = variables[0], variables[1], variables[2]

    # num = (x**2 + y - 11) ** 2 + (x + y**2 - 7) ** 2 + z**2
    # den = 1.0 + (x**2 + y**2) * 0.01
    # r_obj = num / den
    if dim != 3:
        raise ValueError("Himmelblau Ratio test only supports dim=3")

    num_poly = Polynomial(
        [
            to_term(1.0, [4, 0, 0]),
            to_term(2.0, [2, 1, 0]),
            to_term(1.0, [0, 2, 0]),
            to_term(-22.0, [2, 0, 0]),
            to_term(-22.0, [0, 1, 0]),
            to_term(121.0, [0, 0, 0]),
            to_term(1.0, [2, 0, 0]),
            to_term(2.0, [1, 2, 0]),
            to_term(1.0, [0, 4, 0]),
            to_term(-14.0, [1, 0, 0]),
            to_term(-14.0, [0, 2, 0]),
            to_term(49.0, [0, 0, 0]),
            to_term(1.0, [0, 0, 2]),
        ]
    )
    den_poly = Polynomial(
        [to_term(1.0, [0, 0, 0]), to_term(0.01, [2, 0, 0]), to_term(0.01, [0, 2, 0])]
    )
    objective = Rational(num_poly, den_poly)

    # constraint_maker(x, y, z, A, O): return A(x**2 + y**2 <= 50)
    constraints = And(x**2 + y**2 <= 50)

    # BoxN(Point((-5, -5, -1)), Point((5, 5, 1)))
    initial_box = BoxN(Point([-5.0, -5.0, -1.0]), Point([5.0, 5.0, 1.0]))

    return Problem(
        "Himmelblau Ratio",
        dim,
        objective,
        constraints,
        initial_box,
        variables,
    )


def test_split_islands(dim: int) -> Problem:
    variables = [Variable(f"x{i}") for i in range(dim)]
    x, y, z = variables[0], variables[1], variables[2]

    # r_obj = (x + y + z) / 1.0
    numerator_poly = Polynomial(
        [to_term(1.0, [1, 0, 0]), to_term(1.0, [0, 1, 0]), to_term(1.0, [0, 0, 1])]
    )
    denominator_poly = Polynomial([to_term(1.0, [0, 0, 0])])
    objective = Rational(numerator_poly, denominator_poly)

    # O((x - 2) ** 2 + y**2 + z**2 <= 0.25, (x + 2) ** 2 + y**2 + z**2 <= 0.25)
    constraints = Or(
        (x - 2) ** 2 + y**2 + z**2 <= 0.25, (x + 2) ** 2 + y**2 + z**2 <= 0.25
    )

    # BoxN(Point((-5, -2, -2)), Point((5, 2, 2)))
    initial_box = BoxN(Point([-5.0] * dim), Point([5.0] * dim))
    if dim > 1:
        initial_box = BoxN(
            Point([-5.0, -2.0] + [-2.0] * (dim - 2)),
            Point([5.0, 2.0] + [2.0] * (dim - 2)),
        )

    return Problem(
        "Split Islands",
        dim,
        objective,
        constraints,
        initial_box,
        variables,
    )


def test_singularity_edge(dim: int) -> Problem:
    variables = [Variable(f"x{i}") for i in range(dim)]
    x, y, z = variables[0], variables[1], variables[2]

    # num = y * z + x * z + x * y
    # den = x * y * z
    # r_obj = num / den
    numerator_poly = Polynomial(
        [to_term(1.0, [0, 1, 1]), to_term(1.0, [1, 0, 1]), to_term(1.0, [1, 1, 0])]
    )
    denominator_poly = Polynomial([to_term(1.0, [1, 1, 1])])
    objective = Rational(numerator_poly, denominator_poly)

    # constraint_maker(x, y, z, A, O): return A(x**2 + y**2 + z**2 >= 1)
    constraints = And(x**2 + y**2 + z**2 >= 1)

    # BoxN(Point((0.1, 0.1, 0.1)), Point((3, 3, 3)))
    initial_box = BoxN(Point([0.1] * dim), Point([3.0] * dim))

    return Problem(
        "Singularity Edge",
        dim,
        objective,
        constraints,
        initial_box,
        variables,
    )


def test_pole_avoidance(dim: int) -> Problem:
    variables = [Variable(f"x{i}") for i in range(dim)]
    x, y, z = variables[0], variables[1], variables[2]

    # r_obj = PolyBuilder.const(1.0, 3) / (x + y + z - 2.5)
    numerator_poly = Polynomial([to_term(1.0, [0, 0, 0])])
    denominator_poly = Polynomial(
        [
            to_term(1.0, [1, 0, 0]),
            to_term(1.0, [0, 1, 0]),
            to_term(1.0, [0, 0, 1]),
            to_term(-2.5, [0, 0, 0]),
        ]
    )
    objective = Rational(numerator_poly, denominator_poly)

    # constraint_maker(x, y, z, A, O): return A(x >= 1, y >= 1, z >= 1)
    constraints = And(x >= 1, y >= 1, z >= 1)

    # BoxN(Point((1, 1, 1)), Point((2, 2, 2)))
    initial_box = BoxN(Point([1.0] * dim), Point([2.0] * dim))

    return Problem(
        "Pole Avoidance",
        dim,
        objective,
        constraints,
        initial_box,
        variables,
    )


def test_rational_valley(dim: int) -> Problem:
    variables = [Variable(f"x{i}") for i in range(dim)]
    x, y, z = variables[0], variables[1], variables[2]

    # r_obj = (x**2 + y**2 + z**2 + 1) / (x * y * z + 1)
    numerator_poly = Polynomial(
        [
            to_term(1.0, [2, 0, 0]),
            to_term(1.0, [0, 2, 0]),
            to_term(1.0, [0, 0, 2]),
            to_term(1.0, [0, 0, 0]),
        ]
    )
    denominator_poly = Polynomial([to_term(1.0, [1, 1, 1]), to_term(1.0, [0, 0, 0])])
    objective = Rational(numerator_poly, denominator_poly)

    # constraint_maker(x, y, z, A, O): return A(x <= 2, y <= 2, z <= 2)
    constraints = And(x <= 2, y <= 2, z <= 2)

    # BoxN(Point((0.1, 0.1, 0.1)), Point((2, 2, 2)))
    initial_box = BoxN(Point([0.1] * dim), Point([2.0] * dim))

    return Problem(
        "Rational Valley",
        dim,
        objective,
        constraints,
        initial_box,
        variables,
    )


def test_positive_islands(dim: int) -> Problem:
    variables = [Variable(f"x{i}") for i in range(dim)]
    x, y, z = variables[0], variables[1], variables[2]

    # r_obj = (x**2 + 1) / 1.0
    numerator_poly = Polynomial([to_term(1.0, [2, 0, 0]), to_term(1.0, [0, 0, 0])])
    denominator_poly = Polynomial([to_term(1.0, [0, 0, 0])])
    objective = Rational(numerator_poly, denominator_poly)

    # O((x - 2) ** 2 + y**2 + z**2 <= 0.25, (x + 2) ** 2 + y**2 + z**2 <= 0.25)
    constraints = Or(
        (x - 2) ** 2 + y**2 + z**2 <= 0.25, (x + 2) ** 2 + y**2 + z**2 <= 0.25
    )

    # BoxN(Point((-5, -5, -5)), Point((5, 5, 5)))
    initial_box = BoxN(Point([-5.0] * dim), Point([5.0] * dim))

    return Problem(
        "Positive Islands",
        dim,
        objective,
        constraints,
        initial_box,
        variables,
    )


def test_sparse_intersection(dim: int) -> Problem:
    variables = [Variable(f"x{i}") for i in range(dim)]
    x, y, z = variables[0], variables[1], variables[2]

    # num = (x - 1) ** 2 + (y - 1) ** 2 + (z - 1) ** 2
    # den = x + y + z
    # r_obj = num / den
    numerator_poly = Polynomial(
        [
            to_term(1.0, [2, 0, 0]),
            to_term(-2.0, [1, 0, 0]),
            to_term(1.0, [0, 0, 0]),
            to_term(1.0, [0, 2, 0]),
            to_term(-2.0, [0, 1, 0]),
            to_term(1.0, [0, 0, 0]),
            to_term(1.0, [0, 0, 2]),
            to_term(-2.0, [0, 0, 1]),
            to_term(1.0, [0, 0, 0]),
        ]
    )
    denominator_poly = Polynomial(
        [to_term(1.0, [1, 0, 0]), to_term(1.0, [0, 1, 0]), to_term(1.0, [0, 0, 1])]
    )
    objective = Rational(numerator_poly, denominator_poly)

    # constraint_maker(x, y, z, A, O): return A(*_get_3sphere_constraints(x, y, z))
    constraints = _get_3sphere_constraints_formula(variables)

    # BoxN(Point((1, 1, 1)), Point((10, 10, 10)))
    initial_box = BoxN(Point([1.0] * dim), Point([10.0] * dim))

    return Problem(
        "Sparse Intersection",
        dim,
        objective,
        constraints,
        initial_box,
        variables,
    )


def test_main_example(dim: int) -> Problem:
    variables = [Variable(f"x{i}") for i in range(dim)]
    x, y, z = variables[0], variables[1], variables[2]

    # num = (x - 1) ** 2 + (y - 1) ** 2 + (z - 1) ** 2
    # den = x**2 + y**2 + z**2
    # r_obj = num / den
    numerator_poly = Polynomial(
        [
            to_term(1.0, [2, 0, 0]),
            to_term(-2.0, [1, 0, 0]),
            to_term(1.0, [0, 0, 0]),
            to_term(1.0, [0, 2, 0]),
            to_term(-2.0, [0, 1, 0]),
            to_term(1.0, [0, 0, 0]),
            to_term(1.0, [0, 0, 2]),
            to_term(-2.0, [0, 0, 1]),
            to_term(1.0, [0, 0, 0]),
        ]
    )
    denominator_poly = Polynomial(
        [to_term(1.0, [2, 0, 0]), to_term(1.0, [0, 2, 0]), to_term(1.0, [0, 0, 2])]
    )
    objective = Rational(numerator_poly, denominator_poly)

    # constraint_maker(x, y, z, A, O): return A(*_get_3sphere_constraints(x, y, z))
    constraints = _get_3sphere_constraints_formula(variables)

    # BoxN(Point((1, 1, 1)), Point((10, 10, 10)))
    initial_box = BoxN(Point([1.0] * dim), Point([10.0] * dim))

    return Problem(
        "Main Example",
        dim,
        objective,
        constraints,
        initial_box,
        variables,
    )


# list of test suite problems
def all_standard_tests(dim: int) -> List[Problem]:
    return [
        test_sanity_poly(dim),
        test_sanity_rational(dim),
        test_rational_bowl(dim),
        # test_himmelblau_ratio(dim),  # Only supports dim=3
        test_split_islands(dim),
        test_singularity_edge(dim),
        test_pole_avoidance(dim),
        test_rational_valley(dim),
        test_positive_islands(dim),
        test_sparse_intersection(dim),
        test_main_example(dim),
    ]
