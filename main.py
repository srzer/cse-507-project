import time

from box.split import SplitGradient, SplitLongestSide
from objective.bound.affine import AffineBounds, _affine_bounds
from objective.bound.bernstein import BernsteinBounds
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

# TODO: make use of asserts vs value errors consistent with whether function is exposed to client

if __name__ == "__main__":
    global_algo = GlobalMinBranchAndBound()
    improved_algo = ImprovedGlobalMinBranchAndBound()
    baseline_algo = BaselineMin()

    dim = 3
    min_box_size = 0.1
    delta = 1e-3  # δ for dReal
    # stop splitting when max side < min_box_size
    # pruning margin
    # both above values were set to 1e-3
    err = 1e-4
    init_box = initial_constraint_box(dim)
    fn_obj = rational_objective_example(dim)
    vars = global_algo._default_variables(dim)
    constr = ball_constraint_example(vars)
    splitter = SplitGradient()
    # splitter = SplitLongestSide()
    # bounder = AffineBounds()
    bounder = BernsteinBounds()

    args = (
        dim,
        init_box,
        fn_obj,
        vars,
        constr,
        splitter,
        bounder,
        min_box_size,
        delta,
        err,
    )

    print(f"running with parameters:")
    print(f"  dim          = {dim}")
    print(f"  init_box     = {init_box}")
    print(f"  obj fn       = {fn_obj}")
    print(f"  contraint    = {constr}")
    print(f"  min_box_size = {min_box_size}")
    print(f"  delta        = {delta}")
    print(f"  err          = {err}")
    print(f"  splitter     = {splitter}")
    print(f"  bounder      = {bounder}")

    t0 = time.time()
    print("\n=== Running GlobalMinBranchAndBound ===")
    result1 = global_algo(*args)
    t1 = time.time()
    print(f"Result: {result1}")
    print("Branch-and-bound time:", t1 - t0, "seconds")

    time.sleep(2)
    t2 = time.time()
    print("\n=== Running ImprovedGlobalMinBranchAndBound ===")
    result2 = improved_algo(*args)
    t3 = time.time()
    print(f"Result: {result2}")
    print("Improved Branch-and-bound time:", t3 - t2, "seconds")

    time.sleep(2)
    t4 = time.time()
    print("\n=== Running BaselineMin ===")
    result3 = baseline_algo(*args)
    t5 = time.time()
    print(f"Result: {result3}")
    print("Baseline dReal Minimize time:", t5 - t4, "seconds")
