"""
Problem definitions for constrained optimization.

This module provides a unified Problem class and
example problems that demonstrate different optimization scenarios.
"""

from typing import Callable, List, Any
from dataclasses import dataclass
from dreal import *

from box import Box, Box3


@dataclass
class Problem:
    """
    Problem definition for constrained optimization.

    Attributes:
        name: Human-readable problem name
        dimension: Number of variables
        objective: Function that takes symbolic variables and returns objective expression
        constraint: Function that takes symbolic variables and returns constraint expression
        initial_box: Initial search box (Box or Box3)
        evaluate_objective: Function to numerically evaluate objective at a point
        description: Optional description of the problem
    """

    name: str
    dimension: int
    objective: Callable[..., Formula]
    constraint: Callable[..., Formula]
    initial_box: Box
    evaluate_objective: Callable[[List[float]], float]
    description: str = ""


def create_3d_rational_problem() -> Problem:
    def objective(*vars):
        x, y, z = vars[:3]
        return ((x - 1) ** 2 + (y - 1) ** 2 + (z - 1) ** 2) / (x + y + z)

    def constraint(*vars):
        x, y, z = vars[:3]
        return And(
            (x - 3) ** 2 + (y - 3) ** 2 + (z - 3) ** 2 <= 2 * 2,
            (x - 4) ** 2 + (y - 4) ** 2 + (z - 4) ** 2 <= 2 * 2,
            (x - 3) ** 2 + (y - 4) ** 2 + (z - 4) ** 2 <= 2 * 2,
        )

    def evaluate(point):
        x, y, z = point[:3]
        return ((x - 1) ** 2 + (y - 1) ** 2 + (z - 1) ** 2) / (x + y + z)

    return Problem(
        name="3D Rational Function",
        dimension=3,
        objective=objective,
        constraint=constraint,
        initial_box=Box3((1, 1, 1), (10, 10, 10)),
        evaluate_objective=evaluate,
        description="3D rational function optimization with ball constraints",
    )


def create_nd_rational_problem(n: int = 3) -> Problem:
    def objective(*vars):
        xs = vars[:n]
        num = sum((x - 1) ** 2 for x in xs)
        denom = sum(x**2 for x in xs)
        return num / denom

    def constraint(*vars):
        xs = vars[:n]

        # Build three centers similar to the 3D case
        center1 = [3] * n
        center2 = [4] * n
        center3 = [3] + [4] * (n - 1)

        def sq_dist(xs, center):
            return sum((x - c) ** 2 for x, c in zip(xs, center))

        return And(
            sq_dist(xs, center1) <= 2 * 2,
            sq_dist(xs, center2) <= 2 * 2,
            sq_dist(xs, center3) <= 2 * 2,
        )

    def evaluate(point):
        xs = point[:n]
        num = sum((x - 1) ** 2 for x in xs)
        denom = sum(x**2 for x in xs)
        return num / denom

    lows = [1.0] * n
    highs = [10.0] * n

    return Problem(
        name=f"{n}D Rational Function",
        dimension=n,
        objective=objective,
        constraint=constraint,
        initial_box=Box(tuple(lows), tuple(highs)),
        evaluate_objective=evaluate,
        description=f"{n}-dimensional rational function optimization with ball constraints",
    )


def create_bernstein_problem(n: int = 3) -> Problem:
    def objective(*vars):
        xs = vars[:n]
        # f(x) = sum_i (xi - 1)^2 / sum_i xi^2 (same as ND case for simplicity)
        num = sum((x - 1) ** 2 for x in xs)
        denom = sum(x**2 for x in xs)
        return num / denom

    def constraint(*vars):
        xs = vars[:n]

        # Same constraint as ND case
        center1 = [3] * n
        center2 = [4] * n
        center3 = [3] + [4] * (n - 1)

        def sq_dist(xs, center):
            return sum((x - c) ** 2 for x, c in zip(xs, center))

        return And(
            sq_dist(xs, center1) <= 2 * 2,
            sq_dist(xs, center2) <= 2 * 2,
            sq_dist(xs, center3) <= 2 * 2,
        )

    def evaluate(point):
        xs = point[:n]
        num = sum((x - 1) ** 2 for x in xs)
        denom = sum(x**2 for x in xs)
        return num / denom

    lows = [1.0] * n
    highs = [10.0] * n

    return Problem(
        name=f"Bernstein {n}D Function",
        dimension=n,
        objective=objective,
        constraint=constraint,
        initial_box=Box(tuple(lows), tuple(highs)),
        evaluate_objective=evaluate,
        description=f"{n}-dimensional function optimization with Bernstein polynomial bounds",
    )


def create_test_problems() -> List[Problem]:
    # Simple quadratic in 2D
    def quad_obj(*vars):
        x, y = vars[:2]
        return x**2 + y**2

    def quad_constraint(*vars):
        x, y = vars[:2]
        return And(x >= -2, x <= 2, y >= -2, y <= 2)

    def quad_eval(point):
        x, y = point[:2]
        return x**2 + y**2

    quad_problem = Problem(
        name="Simple 2D Quadratic",
        dimension=2,
        objective=quad_obj,
        constraint=quad_constraint,
        initial_box=Box((-2, -2), (2, 2)),
        evaluate_objective=quad_eval,
        description="Simple 2D quadratic for testing",
    )

    # Rational function with singularity
    def rational_obj(*vars):
        x, y = vars[:2]
        return (x**2 + y**2 + 1) / (x * y + 1)

    def rational_constraint(*vars):
        x, y = vars[:2]
        return And(x <= 2, y <= 2)

    def rational_eval(point):
        x, y = point[:2]
        return (x**2 + y**2 + 1) / (x * y + 1)

    rational_problem = Problem(
        name="2D Rational with Singularity",
        dimension=2,
        objective=rational_obj,
        constraint=rational_constraint,
        initial_box=Box((0.1, 0.1), (2, 2)),
        evaluate_objective=rational_eval,
        description="2D rational function with potential singularities",
    )

    return [
        quad_problem,
        rational_problem,
        create_3d_rational_problem(),
        create_nd_rational_problem(3),
        create_nd_rational_problem(4),
        create_bernstein_problem(3),
    ]


# get a problem by name from test suite
def get_problem_by_name(name: str) -> Problem:
    problems = {p.name: p for p in create_test_problems()}
    if name not in problems:
        raise ValueError(
            f"Problem '{name}' not found. Available: {list(problems.keys())}"
        )
    return problems[name]
