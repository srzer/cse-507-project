#!/usr/bin/env python3
"""
Usage:
    python main.py                           # Run default demos
    python main.py --problem "3D Rational Function"  # Run specific problem
    python main.py --method basic --problem "Simple 2D Quadratic"
    python main.py --list-problems           # List available problems
    python main.py --compare                 # Compare different methods

Examples:
    docker compose run --rm dreal main.py
    docker compose run --rm dreal main.py --problem "3D Rational Function"
    docker compose run --rm dreal main.py --compare
"""

import argparse
import sys
import time
from typing import List, Dict

from optimization import optimize, BranchAndBoundSolver, BasicBounds, BernsteinBounds
from optimization.problems import create_test_problems, get_problem_by_name, Problem
from optimization.solver import baseline_solve, OptimizationResult


def print_header():
    print("=" * 80)
    print("  positive rational function optimization")
    print("  branch-and-bound, interval arithmethic, bernstein polynomial form")
    print("=" * 80)
    print()


def print_result(problem: Problem, result: OptimizationResult, verbose: bool = True):
    """Print optimization result in a formatted way."""
    if verbose:
        print(f"Problem: {problem.name}")
        print(f"Description: {problem.description}")
        print(f"Dimension: {problem.dimension}")
        print(f"Initial box: {problem.initial_box}")
        print("-" * 40)

    print(f"Method: {result.method}")
    print(f"Optimal value: {result.optimal_value:.6f}")
    print(f"Solve time: {result.solve_time:.3f}s")
    print(f"Iterations: {result.iterations}")
    print()


def run_demo():
    print_header()
    print("Running demo suite...")
    print()

    # Select representative problems for demo
    demo_problems = [
        get_problem_by_name("Simple 2D Quadratic"),
        get_problem_by_name("3D Rational Function"),
        get_problem_by_name("3D Rational Function"),
    ]

    methods = ["basic", "baseline"]

    print(f"Testing {len(demo_problems)} problems with {len(methods)} methods...")
    print()

    for i, problem in enumerate(demo_problems):
        print(f"Problem {i + 1}/{len(demo_problems)}: {problem.name}")
        print(f"   {problem.description}")
        print()

        for method in methods:
            try:
                print(f"  Method: {method.upper()}")
                result = optimize(problem, method=method, min_box_size=0.1)
                print(
                    f"     Optimal value: {result.optimal_value:.6f} ({result.solve_time:.3f}s, {result.iterations} iter)"
                )
            except Exception as e:
                print(f"     Failed: {str(e)[:60]}...")
        print()


def run_single_problem(problem_name: str, method: str = "basic"):
    print_header()
    print(f"Single problem optimization: {problem_name}")
    print()

    try:
        problem = get_problem_by_name(problem_name)
        result = optimize(problem, method=method, min_box_size=0.1)
        print_result(problem, result)

    except Exception as e:
        print(f"Error: {e}")
        return


def compare_methods(problem_name: str = "3D Rational Function"):
    """Compare different optimization methods on the same problem."""
    print_header()
    print(f"Method comparison: {problem_name}")
    print()

    try:
        problem = get_problem_by_name(problem_name)
        print(f"Problem: {problem.description}")
        print(f"Dimension: {problem.dimension}, Initial box: {problem.initial_box}")
        print()

        methods = ["basic", "baseline"]
        results = []

        for method in methods:
            print(f"Running {method.upper()} method...")
            try:
                result = optimize(problem, method=method, min_box_size=0.1)
                results.append((method, result))
                print(
                    f"  {result.optimal_value:.6f} in {result.solve_time:.3f}s ({result.iterations} iter)"
                )
            except Exception as e:
                print(f"  Failed: {e}")
                results.append((method, None))

        print()
        print("Comparison summary:")
        print("-" * 60)
        print(
            f"{'Method':<12} {'Optimal Value':<15} {'Time (s)':<10} {'Iterations':<10}"
        )
        print("-" * 60)

        for method, result in results:
            if result:
                print(
                    f"{method.upper():<12} {result.optimal_value:<15.6f} {result.solve_time:<10.3f} {result.iterations:<10}"
                )
            else:
                print(f"{method.upper():<12} {'FAILED':<15} {'-':<10} {'-':<10}")

        print()

    except Exception as e:
        print(f"Error: {e}")


def list_problems():
    print_header()
    print("Available problems:")
    print()

    problems = create_test_problems()

    for i, problem in enumerate(problems, 1):
        print(f"{i:2d}. {problem.name}")
        print(f"     Dimension: {problem.dimension}")
        print(f"     {problem.description}")
        print()


def main():
    parser = argparse.ArgumentParser(
        description="Unified Constrained Optimization Framework",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Usage:")[1] if "Usage:" in __doc__ else "",
    )

    parser.add_argument("--problem", "-p", type=str, help="Specific problem to solve")
    parser.add_argument(
        "--method",
        "-m",
        type=str,
        default="basic",
        choices=["basic", "bernstein", "baseline"],
        help="Optimization method to use",
    )
    parser.add_argument(
        "--list-problems", "-l", action="store_true", help="List all available problems"
    )
    parser.add_argument(
        "--compare", "-c", action="store_true", help="Compare different methods"
    )
    parser.add_argument(
        "--compare-problem",
        type=str,
        default="3D Rational Function",
        help="Problem to use for method comparison",
    )

    args = parser.parse_args()

    if args.list_problems:
        list_problems()
    elif args.compare:
        compare_methods(args.compare_problem)
    elif args.problem:
        run_single_problem(args.problem, args.method)
    else:
        run_demo()


if __name__ == "__main__":
    main()
