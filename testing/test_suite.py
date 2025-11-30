from dataclasses import dataclass
from typing import Callable, List
from dreal import *
import time
import math
import warnings

from scipy.optimize import shgo, dual_annealing, differential_evolution
import ours_test


@dataclass
class TestProblem:
    name: str
    f_maker: Callable
    constraint_maker: Callable
    float_constraint: Callable
    init_box: ours_test.Box3D
    description: str


# Numeric Baseline 1: SHGO
def run_scipy_shgo(p: TestProblem):
    bounds = [
        (p.init_box.xl, p.init_box.xu),
        (p.init_box.yl, p.init_box.yu),
        (p.init_box.zl, p.init_box.zu),
    ]

    def obj(arr):
        x, y, z = arr
        if not p.float_constraint(x, y, z):
            return float("inf")
        return float(p.f_maker(x, y, z))

    try:
        res = shgo(obj, bounds, n=64, iters=3, sampling_method="sobol")
        if not res.success or math.isinf(res.fun):
            return None
        return res.fun
    except Exception:
        return "Error"


# Numeric Baseline 2: Differential Evolution
def run_scipy_diff_evo(p: TestProblem):
    bounds = [
        (p.init_box.xl, p.init_box.xu),
        (p.init_box.yl, p.init_box.yu),
        (p.init_box.zl, p.init_box.zu),
    ]

    def obj(arr):
        x, y, z = arr
        if not p.float_constraint(x, y, z):
            return float("inf")
        return float(p.f_maker(x, y, z))

    try:
        res = differential_evolution(
            obj,
            bounds,
            maxiter=200,
            tol=1e-6,
            polish=True,
            updating="deferred",
        )
        return None if math.isinf(res.fun) else res.fun
    except Exception:
        return "Error"


# Numeric Baseline 3: Dual Annealing
def run_scipy_dual_annealing(p: TestProblem):
    bounds = [
        (p.init_box.xl, p.init_box.xu),
        (p.init_box.yl, p.init_box.yu),
        (p.init_box.zl, p.init_box.zu),
    ]

    def obj(arr):
        x, y, z = arr
        if not p.float_constraint(x, y, z):
            return float("inf")
        return float(p.f_maker(x, y, z))

    try:
        res = dual_annealing(obj, bounds=bounds, maxiter=200)
        return None if math.isinf(res.fun) else res.fun
    except Exception:
        return "Error"


# Test Runner
def run_test_suite(problems: List[TestProblem]):
    header = (
        f"{'TEST NAME':<20} | "
        f"{'MY SOLVER':<10} | "
        f"{'DREAL':<10} | "
        f"{'SHGO':<10} | "
        f"{'DE':<10} | "
        f"{'DA':<10} | "
        f"{'TIME (My/Dr/Sh/DE/DA)'}"
    )
    print(header)
    print("-" * 125)

    for p in problems:
        print(f"Running {p.name}...", end="\r")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # My solver
            t0 = time.time()
            try:
                dc = lambda x, y, z: p.constraint_maker(x, y, z, And, Or)
                my_res = ours_test.global_min_branch_and_bound(
                    p.f_maker, dc, p.init_box, 1e-3, 0.1, 1e-3
                )
            except Exception:
                my_res = "Error"
            t_my = time.time() - t0

            # dReal
            t1 = time.time()
            try:
                dc = lambda x, y, z: p.constraint_maker(x, y, z, And, Or)
                dreal_res = ours_test.baseline_min_dreal(
                    p.f_maker, dc, p.init_box, 1e-3
                )
            except Exception:
                dreal_res = "Error"
            t_dreal = time.time() - t1

            # SHGO builds a simplicial complex over the search domain and uses topology (homology) to systematically explore and refine promising regions. It guarantees systematic global coverage and is good for smooth problems, but can struggle with highly disconnected or very nonlinear constraints.
            t2 = time.time()
            shgo_res = run_scipy_shgo(p)
            t_sh = time.time() - t2

            # Differential Evolution is a population-based global optimizer that explores the search space by repeatedly mutating and recombining candidate solutions using vector differences. It is powerful for nonconvex, noisy, or rugged landscapes and tends to find good global minima with minimal assumptions on smoothness.
            t3 = time.time()
            de_res = run_scipy_diff_evo(p)
            t_de = time.time() - t3

            # Dual Annealing combines classical simulated annealing with a local search phase, allowing it to escape local minima and then refine good regions. It is particularly effective on multimodal landscapes where the optimum is hidden behind many local minima.
            t4 = time.time()
            da_res = run_scipy_dual_annealing(p)
            t_da = time.time() - t4

        def fmt(val):
            if val is None:
                return "None"
            if isinstance(val, str):
                return val
            return f"{val:.4f}"

        print(
            f"{p.name:<20} | "
            f"{fmt(my_res):<10} | "
            f"{fmt(dreal_res):<10} | "
            f"{fmt(shgo_res):<10} | "
            f"{fmt(de_res):<10} | "
            f"{fmt(da_res):<10} | "
            f"{t_my:.2f}/{t_dreal:.2f}/{t_sh:.2f}/{t_de:.2f}/{t_da:.2f}s"
        )


# Test Definitions (same as before)
def test_sanity_poly():
    return TestProblem(
        "Sanity Poly",
        lambda x, y, z: x**2 + y**2 + z**2,
        lambda x, y, z, A, O: A(x >= -2, x <= 2, y >= -2, y <= 2, z >= -2, z <= 2),
        lambda x, y, z: True,
        ours_test.Box3D(-2, 2, -2, 2, -2, 2),
        "Simple x^2",
    )


def test_sanity_rational():
    return TestProblem(
        "Sanity Rational",
        lambda x, y, z: (x + 1.0) / (y + 1.0),
        lambda x, y, z, A, O: A(x >= 0, x <= 1, y >= 0, y <= 1, z == 0),
        lambda x, y, z: 0 <= x <= 1 and 0 <= y <= 1,
        ours_test.Box3D(0, 1, 0, 1, -1, 1),
        "Simple x/y",
    )


def test_rational_bowl():
    return TestProblem(
        "Rational Bowl",
        lambda x, y, z: (x**2 + y**2 + z**2) / (x + y + z),
        lambda x, y, z, A, O: A(x**2 + y**2 + z**2 <= 9),
        lambda x, y, z: x**2 + y**2 + z**2 <= 9,
        ours_test.Box3D(0.1, 5, 0.1, 5, 0.1, 5),
        "Simple rational",
    )


def test_himmelblau_ratio():
    def f(x, y, z):
        num = (x**2 + y - 11) ** 2 + (x + y**2 - 7) ** 2 + z**2
        den = 1 + 0.01 * (x**2 + y**2)
        return num / den

    return TestProblem(
        "Himmelblau Ratio",
        f,
        lambda x, y, z, A, O: A(x**2 + y**2 <= 50),
        lambda x, y, z: x**2 + y**2 <= 50,
        ours_test.Box3D(-5, 5, -5, 5, -1, 1),
        "Multimodal",
    )


def test_split_islands():
    return TestProblem(
        "Split Islands",
        lambda x, y, z: x + y + z,
        lambda x, y, z, A, O: O(
            (x - 2) ** 2 + y**2 + z**2 <= 0.25,
            (x + 2) ** 2 + y**2 + z**2 <= 0.25,
        ),
        lambda x, y, z: (
            (x - 2) ** 2 + y**2 + z**2 <= 0.25
            or (x + 2) ** 2 + y**2 + z**2 <= 0.25
        ),
        ours_test.Box3D(-5, 5, -2, 2, -2, 2),
        "Disconnected",
    )


def test_singularity_edge():
    return TestProblem(
        "Singularity Edge",
        lambda x, y, z: (1 / x) + (1 / y) + (1 / z),
        lambda x, y, z, A, O: x**2 + y**2 + z**2 >= 1,
        lambda x, y, z: x**2 + y**2 + z**2 >= 1,
        ours_test.Box3D(0.1, 3, 0.1, 3, 0.1, 3),
        "Non-convex hole",
    )


def test_pole_avoidance():
    return TestProblem(
        "Pole Avoidance",
        lambda x, y, z: 1 / (x + y + z - 2.5),
        lambda x, y, z, A, O: A(x >= 1, y >= 1, z >= 1),
        lambda x, y, z: x >= 1 and y >= 1 and z >= 1,
        ours_test.Box3D(1, 2, 1, 2, 1, 2),
        "Singularity",
    )


def test_rational_valley():
    return TestProblem(
        "Rational Valley",
        lambda x, y, z: (x**2 + y**2 + z**2 + 1) / (x * y * z + 1),
        lambda x, y, z, A, O: A(x <= 2, y <= 2, z <= 2),
        lambda x, y, z: x <= 2 and y <= 2 and z <= 2,
        ours_test.Box3D(0.1, 2, 0.1, 2, 0.1, 2),
        "Valley",
    )


def test_positive_islands():
    return TestProblem(
        "Positive Islands",
        lambda x, y, z: x**2 + 1,
        lambda x, y, z, A, O: O(
            (x - 2) ** 2 + y**2 + z**2 <= 0.25,
            (x + 2) ** 2 + y**2 + z**2 <= 0.25,
        ),
        lambda x, y, z: (
            (x - 2) ** 2 + y**2 + z**2 <= 0.25
            or (x + 2) ** 2 + y**2 + z**2 <= 0.25
        ),
        ours_test.Box3D(-5, 5, -5, 5, -5, 5),
        "Strict positive disconnected",
    )


def test_sparse_intersection():
    def c(x, y, z, A, O):
        return A(
            (x - 3) ** 2 + (y - 3) ** 2 + (z - 3) ** 2 <= 4,
            (x - 4) ** 2 + (y - 4) ** 2 + (z - 4) ** 2 <= 4,
            (x - 3) ** 2 + (y - 4) ** 2 + (z - 4) ** 2 <= 4,
        )

    def cf(x, y, z):
        return (
            (x - 3) ** 2 + (y - 3) ** 2 + (z - 3) ** 2 <= 4
            and (x - 4) ** 2 + (y - 4) ** 2 + (z - 4) ** 2 <= 4
            and (x - 3) ** 2 + (y - 4) ** 2 + (z - 4) ** 2 <= 4
        )

    return TestProblem(
        "Sparse Intersection",
        lambda x, y, z: ((x - 1) ** 2 + (y - 1) ** 2 + (z - 1) ** 2) / (x + y + z),
        c,
        cf,
        ours_test.Box3D(1, 10, 1, 10, 1, 10),
        "Hard Intersection",
    )


# Run
if __name__ == "__main__":
    problems = [
        test_rational_bowl(),
        test_himmelblau_ratio(),
        test_split_islands(),
        test_singularity_edge(),
        test_pole_avoidance(),
        test_rational_valley(),
        test_positive_islands(),
        test_sparse_intersection(),
    ]
    run_test_suite(problems)
