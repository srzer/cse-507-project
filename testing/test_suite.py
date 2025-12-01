import time
import math
import warnings
import sys
import os
from dataclasses import dataclass
from typing import Callable, List
from collections import defaultdict

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

bernstein_path = os.path.join(project_root, "bernstein-nd")
if bernstein_path not in sys.path:
    sys.path.append(bernstein_path)

from dreal import Variable, And, Or

from scipy.optimize import shgo, dual_annealing, differential_evolution

from box import BoxN, Point
from poly import Rational, Polynomial, PolyBuilder, make_vars
from poly.type import to_term

from ours import global_min_branch_and_bound, baseline_min_dreal
from poly_utils import poly_from_terms
from box_utils import BoxND



@dataclass
class TestProblem:
    name: str
    rational_obj: Rational
    func_lambda: Callable
    constraint_maker: Callable
    float_constraint: Callable
    init_box: BoxN
    description: str

def run_scipy_shgo(p: TestProblem):
    bounds = list(zip(p.init_box.min, p.init_box.max))
    def obj(arr):
        if not p.float_constraint(*arr): return float("inf")
        return float(p.func_lambda(*arr))
    try:
        res = shgo(obj, bounds, n=64, iters=3, sampling_method="sobol")
        if not res.success or math.isinf(res.fun): return None
        return res.fun
    except Exception: return "Error"

def run_scipy_diff_evo(p: TestProblem):
    bounds = list(zip(p.init_box.min, p.init_box.max))
    def obj(arr):
        if not p.float_constraint(*arr): return float("inf")
        return float(p.func_lambda(*arr))
    try:
        res = differential_evolution(obj, bounds, maxiter=200, tol=1e-6, polish=True)
        return None if math.isinf(res.fun) else res.fun
    except Exception: return "Error"

def run_scipy_dual_annealing(p: TestProblem):
    bounds = list(zip(p.init_box.min, p.init_box.max))
    def obj(arr):
        if not p.float_constraint(*arr): return float("inf")
        return float(p.func_lambda(*arr))
    try:
        res = dual_annealing(obj, bounds=bounds, maxiter=200)
        return None if math.isinf(res.fun) else res.fun
    except Exception: return "Error"


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

        dim = p.init_box.dim
        dreal_vars = [Variable(f"x{i}") for i in range(dim)]
        dreal_constraint = p.constraint_maker(*dreal_vars, And, Or)
        
        # 1. Convert Polynomial (List) to Dictionary
        num_dict = poly_from_terms(p.rational_obj.num)
        den_dict = poly_from_terms(p.rational_obj.den)
        
        # 2. Convert BoxN (min/max) to BoxND (lows/highs)
        box_nd = BoxND(list(p.init_box.min), list(p.init_box.max))

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # 1. My Solver
            t0 = time.time()
            try:
                my_res = global_min_branch_and_bound(
                    initial_box=box_nd,
                    n=dim,
                    poly_num=num_dict,
                    poly_den=den_dict,
                    vars=dreal_vars,
                    constraint=dreal_constraint,
                    delta_dreal=1e-3,
                    min_box_size=1e-3,
                    eps=1e-3
                )
            except Exception as e:
                # print(e)
                my_res = "Error"
            t_my = time.time() - t0

            # 2. dReal Baseline
            t1 = time.time()
            try:
                dreal_res = baseline_min_dreal(
                    initial_box=box_nd,
                    n=dim,
                    poly_num=num_dict,
                    poly_den=den_dict,
                    vars=dreal_vars,
                    constraint=dreal_constraint,
                    delta_dreal=1e-3
                )
            except Exception:
                dreal_res = "Error"
            t_dreal = time.time() - t1

            # 3. SciPy Baselines
            t2 = time.time()
            shgo_res = run_scipy_shgo(p)
            t_sh = time.time() - t2

            t3 = time.time()
            de_res = run_scipy_diff_evo(p)
            t_de = time.time() - t3

            t4 = time.time()
            da_res = run_scipy_dual_annealing(p)
            t_da = time.time() - t4

        def fmt(val):
            if val is None: return "None"
            if isinstance(val, str): return val
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



def test_sanity_poly():
    x, y, z = make_vars(3)
    r_obj = (x**2 + y**2 + z**2) / 1.0
    return TestProblem("Sanity Poly", r_obj, lambda x,y,z: x**2+y**2+z**2, lambda x,y,z,A,O: A(x>=-2,x<=2,y>=-2,y<=2,z>=-2,z<=2), lambda x,y,z: True, BoxN(Point((-2, -2, -2)), Point((2, 2, 2))), "Simple x^2")

def test_sanity_rational():
    x, y, z = make_vars(3)
    r_obj = (x + 1.0) / (y + 1.0)
    return TestProblem("Sanity Rational", r_obj, lambda x,y,z: (x+1)/(y+1), lambda x,y,z,A,O: A(x>=0,x<=1,y>=0,y<=1,z==0), lambda x,y,z: 0<=x<=1 and 0<=y<=1, BoxN(Point((0,0,-1)), Point((1,1,1))), "Simple x/y")

def test_rational_bowl():
    x, y, z = make_vars(3)
    r_obj = (x**2 + y**2 + z**2) / (x + y + z)
    return TestProblem("Rational Bowl", r_obj, lambda x,y,z: (x**2+y**2+z**2)/(x+y+z), lambda x,y,z,A,O: A(x**2+y**2+z**2<=9), lambda x,y,z: x**2+y**2+z**2<=9, BoxN(Point((0.1,0.1,0.1)), Point((5,5,5))), "Simple rational")

def test_himmelblau_ratio():
    x, y, z = make_vars(3)
    num = (x**2 + y - 11)**2 + (x + y**2 - 7)**2 + z**2
    den = 1.0 + (x**2 + y**2) * 0.01
    r_obj = num / den
    def f_lambda(x, y, z): return ((x**2 + y - 11)**2 + (x + y**2 - 7)**2 + z**2) / (1 + 0.01*(x**2+y**2))
    return TestProblem("Himmelblau Ratio", r_obj, f_lambda, lambda x,y,z,A,O: A(x**2+y**2<=50), lambda x,y,z: x**2+y**2<=50, BoxN(Point((-5,-5,-1)), Point((5,5,1))), "Multimodal")

def test_split_islands():
    x, y, z = make_vars(3)
    r_obj = (x + y + z) / 1.0
    return TestProblem("Split Islands", r_obj, lambda x,y,z: x+y+z, lambda x,y,z,A,O: O((x-2)**2+y**2+z**2<=0.25, (x+2)**2+y**2+z**2<=0.25), lambda x,y,z: (x-2)**2+y**2+z**2<=0.25 or (x+2)**2+y**2+z**2<=0.25, BoxN(Point((-5,-2,-2)), Point((5,2,2))), "Disconnected")

def test_singularity_edge():
    x, y, z = make_vars(3)
    num = y*z + x*z + x*y
    den = x*y*z
    r_obj = num / den
    return TestProblem("Singularity Edge", r_obj, lambda x,y,z: 1/x+1/y+1/z, lambda x,y,z,A,O: x**2+y**2+z**2>=1, lambda x,y,z: x**2+y**2+z**2>=1, BoxN(Point((0.1,0.1,0.1)), Point((3,3,3))), "Non-convex hole")

def test_pole_avoidance():
    x, y, z = make_vars(3)
    r_obj = PolyBuilder.const(1.0, 3) / (x + y + z - 2.5)
    return TestProblem("Pole Avoidance", r_obj, lambda x,y,z: 1/(x+y+z-2.5), lambda x,y,z,A,O: A(x>=1,y>=1,z>=1), lambda x,y,z: x>=1 and y>=1 and z>=1, BoxN(Point((1,1,1)), Point((2,2,2))), "Singularity")

def test_rational_valley():
    x, y, z = make_vars(3)
    r_obj = (x**2 + y**2 + z**2 + 1) / (x*y*z + 1)
    return TestProblem("Rational Valley", r_obj, lambda x,y,z: (x**2+y**2+z**2+1)/(x*y*z+1), lambda x,y,z,A,O: A(x<=2,y<=2,z<=2), lambda x,y,z: x<=2 and y<=2 and z<=2, BoxN(Point((0.1,0.1,0.1)), Point((2,2,2))), "Valley")

def test_positive_islands():
    x, y, z = make_vars(3)
    r_obj = (x**2 + 1) / 1.0
    return TestProblem("Positive Islands", r_obj, lambda x,y,z: x**2+1, lambda x,y,z,A,O: O((x-2)**2+y**2+z**2<=0.25, (x+2)**2+y**2+z**2<=0.25), lambda x,y,z: (x-2)**2+y**2+z**2<=0.25 or (x+2)**2+y**2+z**2<=0.25, BoxN(Point((-5,-5,-5)), Point((5,5,5))), "Strict positive disconnected")

def test_sparse_intersection():
    x, y, z = make_vars(3)
    num = (x - 1)**2 + (y - 1)**2 + (z - 1)**2
    den = x + y + z
    r_obj = num / den

    def c_maker(x, y, z, A, O):
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

    return TestProblem("Sparse Intersection", r_obj, lambda x, y, z: ((x - 1) ** 2 + (y - 1) ** 2 + (z - 1) ** 2) / (x + y + z), c_maker, cf, BoxN(Point((1, 1, 1)), Point((10, 10, 10))), "Hard Intersection")

def test_main_example():
    dim = 3
    x, y, z = make_vars(dim)
    
    num = (x - 1)**2 + (y - 1)**2 + (z - 1)**2
    den = x**2 + y**2 + z**2
    r_obj = num / den


    def constraint_maker(x, y, z, A, O):
        return A(
            (x - 3)**2 + (y - 3)**2 + (z - 3)**2 <= 4,
            (x - 4)**2 + (y - 4)**2 + (z - 4)**2 <= 4,
            (x - 3)**2 + (y - 4)**2 + (z - 4)**2 <= 4
        )

    def float_constraint(x, y, z):
        c1 = (x - 3)**2 + (y - 3)**2 + (z - 3)**2 <= 4
        c2 = (x - 4)**2 + (y - 4)**2 + (z - 4)**2 <= 4
        c3 = (x - 3)**2 + (y - 4)**2 + (z - 4)**2 <= 4
        return c1 and c2 and c3

    func_lambda = lambda x, y, z: ((x - 1)**2 + (y - 1)**2 + (z - 1)**2) / (x**2 + y**2 + z**2)

    return TestProblem(
        "Main Example", 
        r_obj, 
        func_lambda, 
        constraint_maker, 
        float_constraint, 
        BoxN(Point((1, 1, 1)), Point((10, 10, 10))), 
        "Quadratic/Quadratic Rational"
    )
