import time
import math
import warnings
import pandas as pd
from typing import List, Dict, Any

from scipy.optimize import shgo, dual_annealing, differential_evolution

from testing.scipy_test_definitions import NumericProblem, all_numeric_tests


# def run_smt_algo(algo, problem: NumericProblem) -> float:
#     """Wrapper to run SMT-based algorithms."""
#     # Note: This runner uses the legacy `NumericProblem` definition
#     # which is different from the `Problem` type used in `main.py`
#     dim = problem.init_box.dim
#     dreal_vars = [Variable(f"x{i}") for i in range(dim)]
#     dreal_constraint = problem.constraint_maker(*dreal_vars, And, Or)

#     # The SMT algos we are using don't conform to the new interface yet,
#     # so we call them with the old, unpacked parameters.
#     # This is a temporary bridge during refactoring.
#     result = algo(
#         dim=dim,
#         init_box=problem.init_box,
#         obj=problem.rational_obj,
#         vars=dreal_vars,
#         constr=dreal_constraint,
#         min_box_size=1e-3,
#         delta=1e-3,
#         err=1e-3,
#     )
#     if isinstance(result, Success):
#         bound, _ = result.unwrap()
#         return bound
#     raise Exception(f"Algorithm failed: {result.failure()}")


def run_scipy_algo(solver_func, problem: NumericProblem) -> float:
    """Wrapper to run SciPy-based algorithms."""
    bounds = list(zip(problem.init_box.min, problem.init_box.max))

    def obj_wrapper(arr):
        # SciPy solvers don't natively handle constraints, so they must
        # be incorporated into the objective function.
        if not problem.float_constraint(*arr):
            return float("inf")
        return float(problem.func_lambda(*arr))

    # Default parameters are based on the legacy test suite
    if solver_func == shgo:
        res = solver_func(obj_wrapper, bounds, n=64, iters=3, sampling_method="sobol")
    elif solver_func == differential_evolution:
        res = solver_func(obj_wrapper, bounds, maxiter=200, tol=1e-6, polish=True)
    elif solver_func == dual_annealing:
        res = solver_func(obj_wrapper, bounds, maxiter=200)
    else:
        raise ValueError(f"Unknown SciPy solver: {solver_func.__name__}")

    if hasattr(res, "success") and not res.success:
        raise Exception("Solver failed to find a solution.")

    if math.isinf(res.fun):
        raise Exception("Solution is at infinity, likely due to constraints.")

    return res.fun


def run_comparison(
    problems: List[NumericProblem],
    solvers: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """
    Runs a list of problems against a list of solvers and returns structured results.
    """
    all_results = []
    total_runs = len(problems) * len(solvers)
    run_count = 0

    print(
        f"Starting comparison of {len(problems)} problems and {len(solvers)} solvers."
    )

    for problem in problems:
        for solver_config in solvers:
            run_count += 1
            solver_name = solver_config["name"]
            solver_func = solver_config["func"]
            skip_list = solver_config.get("skip", [])

            print(
                f"({run_count}/{total_runs}) "
                f"Running {solver_name} on {problem.name}...       ",
                end="\r",
            )

            if problem.name in skip_list:
                result_row = {
                    "problem": problem.name,
                    "algorithm": solver_name,
                    "runtime": "SKIPPED",
                    "bound": "SKIPPED",
                }
                all_results.append(result_row)
                continue

            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    t_start = time.time()
                    bound = solver_func(solver_config["instance"], problem)
                    t_end = time.time()

                result_row = {
                    "problem": problem.name,
                    "algorithm": solver_name,
                    "runtime": round(t_end - t_start, 4),
                    "bound": round(bound, 6),
                }

            except Exception as e:
                print(e)
                result_row = {
                    "problem": problem.name,
                    "algorithm": solver_name,
                    "runtime": "Error",
                    "bound": "Error",
                }

            all_results.append(result_row)

    print("\nComparison finished.")
    return all_results


def save_results_to_csv(results: List[Dict[str, Any]], filename: str):
    """Saves a list of result dictionaries to a CSV file."""
    df = pd.DataFrame(results)
    df.to_csv(filename, index=False)
    print(f"Results saved to {filename}")


def run_numeric_tests():
    # Define the list of problems from the legacy definitions
    problems_to_run = all_numeric_tests()

    # Define the solver configurations to be tested
    # The 'instance' is the object to call, 'func' is the wrapper to use.
    SOLVER_CONFIGS = [
        {
            "name": "SHGO",
            "instance": shgo,
            "func": run_scipy_algo,
            # "skip": ["Split Islands", "Sparse Intersection"],
        },
        {
            "name": "Diff. Evo.",
            "instance": differential_evolution,
            "func": run_scipy_algo,
        },
        {"name": "Dual Annealing", "instance": dual_annealing, "func": run_scipy_algo},
    ]

    # Run the comparison
    results = run_comparison(problems_to_run, SOLVER_CONFIGS)

    # Save the results to a CSV file
    save_results_to_csv(results, "report/log/numeric_comparison.csv")
