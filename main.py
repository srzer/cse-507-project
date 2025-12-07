import time
from datetime import datetime
from itertools import product

from box import SplitGradient, SplitLongestSide
from objective import AffineBounds, BernsteinBounds
from returns.result import Success, Failure

from testing.example import make_example_problem_1, make_example_problem_2
from testing.standard_tests import all_standard_tests
from algorithms import (
    GlobalMinBranchAndBound,
    ImprovedGlobalMinBranchAndBound,
    BaselineMin,
)
from algorithms.log import save_logs_to_csv
from testing.logging_table import results_to_markdown, save_results_to_csv

# main config
GENERATE_DETAILED_LOGS = False
PROBLEM_DIMENSION = 3

# skip config for known problematic tests
# This list defines combinations of problem/algorithm/heuristics that
# are known to hang indefinitely due to dreal limitations.
SKIP_LIST = [
    {
        "problem": "Example Problem 2",
        "algorithm": "GlobalMinBranchAndBound",
        "splitter": "SplitLongestSide",
    },
    {
        "problem": "Sparse Intersection",
        "algorithm": "ImprovedGlobalMinBranchAndBound",
        "splitter": "SplitGradient",
    },
    {
        "problem": "Sparse Intersection",
        "algorithm": "GlobalMinBranchAndBound",
        "splitter": "SplitGradient",
    },
]

if __name__ == "__main__":
    problems = [
        # make_example_problem_1(PROBLEM_DIMENSION),
        # make_example_problem_2(PROBLEM_DIMENSION),
    ]
    problems.extend(all_standard_tests(PROBLEM_DIMENSION))

    bnb_algorithms = [GlobalMinBranchAndBound(), ImprovedGlobalMinBranchAndBound()]
    splitters = [SplitLongestSide(), SplitGradient()]
    bounders = [AffineBounds(), BernsteinBounds()]

    run_configs = [
        {"algo": algo, "splitter": splitter, "bounder": bounder}
        for algo, splitter, bounder in product(bnb_algorithms, splitters, bounders)
    ]
    run_configs.append({"algo": BaselineMin(), "splitter": None, "bounder": None})

    all_results = []
    print(f"starting test suite on {datetime.now()}")

    for problem in problems:
        print(f"\nrunning problem: {problem.name} dim {problem.dim}")
        for config in run_configs:
            algo = config["algo"]
            splitter = config["splitter"]
            bounder = config["bounder"]

            algo_name = algo.__class__.__name__
            splitter_name = splitter.__class__.__name__ if splitter else "NA"
            bounder_name = bounder.__class__.__name__ if bounder else "NA"

            print(f"  -> running: {algo_name} with {splitter_name} and {bounder_name}")

            is_skipped = False
            for skip_rule in SKIP_LIST:
                if (
                    skip_rule.get("problem") == problem.name
                    and skip_rule.get("algorithm") == algo_name
                    and skip_rule.get("splitter") == splitter_name
                ):
                    is_skipped = True
                    break

            if is_skipped:
                print("    -> SKIPPED (known to hang)")
                all_results.append(
                    {
                        "problem": problem.name,
                        "algorithm": algo_name,
                        "splitter": splitter_name,
                        "bounder": bounder_name,
                        "runtime": "SKIPPED",
                        "bound": "SKIPPED",
                    }
                )
                continue

            t_start = time.time()
            result = algo(
                problem=problem,
                splitter=splitter,
                bounder=bounder,
            )
            t_end = time.time()
            run_time = t_end - t_start

            if isinstance(result, Success):
                bound, logs = result.unwrap()
            else:
                bound, logs = -1.0, []
                print(f"    -> algorithm failed with error: {result.failure()}")

            all_results.append(
                {
                    "problem": problem.name,
                    "algorithm": algo_name,
                    "splitter": splitter_name,
                    "bounder": bounder_name,
                    "runtime": round(run_time, 4),
                    "bound": round(bound, 6),
                }
            )

            if GENERATE_DETAILED_LOGS:
                log_filename = f"logs_{problem.name}_{algo_name}_{splitter_name}_{bounder_name}.csv"
                save_logs_to_csv(logs, log_filename)
                print(
                    f"    -> finished in {run_time:.4f}s. result: {bound:.6f}. logs: {log_filename}"
                )
            else:
                print(f"    -> finished in {run_time:.4f}s. result: {bound:.6f}.")

    print("\noverall comparison summary")
    print(results_to_markdown(all_results))
    save_results_to_csv(all_results, "comparison_summary.csv")
    print("\nsummary saved to comparison_summary.csv")
