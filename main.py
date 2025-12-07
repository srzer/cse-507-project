import time
from datetime import datetime
from itertools import product

from box import SplitGradient, SplitLongestSide
from objective import AffineBounds, BernsteinBounds
from returns.result import Success
from testing.example import (
    rational_objective_example,
    ball_constraint_example,
    initial_constraint_box,
)
from algorithms import (
    GlobalMinBranchAndBound,
    ImprovedGlobalMinBranchAndBound,
    BaselineMin,
)
from algorithms.log import save_logs_to_csv
from testing.reporting import results_to_markdown, save_results_to_csv

if __name__ == "__main__":
    # configuration
    dim = 3
    min_box_size = 0.1
    delta = 1e-3
    err = 1e-4

    # problem definition
    init_box = initial_constraint_box(dim)
    fn_obj = rational_objective_example(dim)
    # assume all algos can use the same default vars for this problem
    shared_vars = GlobalMinBranchAndBound()._default_variables(dim)
    constr = ball_constraint_example(shared_vars)

    # define runs
    bnb_algorithms = [GlobalMinBranchAndBound(), ImprovedGlobalMinBranchAndBound()]
    splitters = [SplitLongestSide(), SplitGradient()]
    bounders = [AffineBounds(), BernsteinBounds()]

    run_configs = [
        {"algo": algo, "splitter": splitter, "bounder": bounder}
        for algo, splitter, bounder in product(bnb_algorithms, splitters, bounders)
    ]
    run_configs.append({"algo": BaselineMin(), "splitter": None, "bounder": None})

    # execution
    all_results = []
    print(f"=== starting algorithm comparison on {datetime.now()} ===")
    print(f"    dim          = {dim}")
    print(f"    min_box_size = {min_box_size}")
    print(f"    delta        = {delta}")
    print(f"    err          = {err}")

    for config in run_configs:
        algo = config["algo"]
        splitter = config["splitter"]
        bounder = config["bounder"]

        algo_name = algo.__class__.__name__
        splitter_name = splitter.__class__.__name__ if splitter else "N/A"
        bounder_name = bounder.__class__.__name__ if bounder else "N/A"

        print(f"running: {algo_name} with {splitter_name} and {bounder_name}")

        t_start = time.time()
        result = algo(
            dim,
            init_box,
            fn_obj,
            shared_vars,
            constr,
            splitter=splitter,
            bounder=bounder,
            min_box_size=min_box_size,
            delta=delta,
            err=err,
        )
        t_end = time.time()
        run_time = t_end - t_start

        if isinstance(result, Success):
            bound, logs = result.unwrap()
        else:
            bound, logs = -1.0, []
            print(f"  -> algorithm failed with error: {result.failure()}")

        all_results.append(
            {
                "algorithm": algo_name,
                "splitter": splitter_name,
                "bounder": bounder_name,
                "runtime (s)": round(run_time, 4),
                "final Bound": round(bound, 6),
            }
        )

        log_filename = f"logs_{algo_name}_{splitter_name}_{bounder_name}.csv"
        save_logs_to_csv(logs, log_filename)
        print(
            f"  -> finished in {run_time:.4f}s. result: {bound:.6f}. logs: {log_filename}"
        )

    print("\n=== comparison summary ===")
    print(results_to_markdown(all_results))
    save_results_to_csv(all_results, "comparison_summary.csv")
    print("\nsummary saved to comparison_summary.csv")
