import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  
from test_suite import *
import ours_test


def visualize_problem_3d(p, results_dict, resolution=60, save_dir="visuals"):
    """
    Creates a 3D surface plot of f(x,y,z_mid) with solver result points.
    results_dict = {"my": (x,y,z,val), "dreal": (...), ...}
    """

    os.makedirs(save_dir, exist_ok=True)

    # Determine slice (z = midpoint of box)
    box = p.init_box
    z_mid = 0.5 * (box.zl + box.zu)

    # Make grid
    xs = np.linspace(box.xl, box.xu, resolution)
    ys = np.linspace(box.yl, box.yu, resolution)
    X, Y = np.meshgrid(xs, ys)
    Z = np.full_like(X, z_mid)

    # Evaluate f and feasibility mask
    F = np.zeros_like(X, dtype=float)
    feas_mask = np.zeros_like(X, dtype=bool)

    for i in range(resolution):
        for j in range(resolution):
            x, y, z = X[i, j], Y[i, j], Z[i, j]
            feas_mask[i, j] = p.float_constraint(x, y, z)
            if feas_mask[i, j]:
                try:
                    F[i, j] = p.f_maker(x, y, z)
                except Exception:
                    F[i, j] = np.nan
            else:
                F[i, j] = np.nan

    # Plot
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    # Surface plot
    surf = ax.plot_surface(
        X, Y, F,
        rstride=1, cstride=1,
        cmap='viridis', alpha=0.85,
        linewidth=0, antialiased=True
    )

    ax.contour(
        X, Y, feas_mask,
        levels=[0.5],
        colors='white',
        linestyles='dashed'
    )

    # Plot solver result points
    colors = {
        "my": "red",
        "dreal": "blue",
        "shgo": "green",
        "de": "orange",
        "da": "purple",
    }

    for key, res in results_dict.items():
        if res is None or isinstance(res, str):
            continue

        pass 

    # Labels
    ax.set_title(f"{p.name} (z = {z_mid:.3f} slice)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("f(x, y, z_mid)")

    # Save image
    outfile = os.path.join(save_dir, f"{p.name}.png")
    plt.savefig(outfile, dpi=150)
    plt.close()
    print(f"[visual] Saved: {outfile}")



def run_with_visuals(problems):
    """
    Runs experiments exactly like run_test_suite(), but also produces 3D visualizations.
    """

    from test_suite import run_test_suite

    print("Running experiments...")
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

    results = {}  

    for p in problems:
        print(f"Running {p.name}...", end="\r")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # My solver
            import time
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

            # SHGO
            from test_suite import run_scipy_shgo, run_scipy_diff_evo, run_scipy_dual_annealing
            t2 = time.time()
            shgo_res = run_scipy_shgo(p)
            t_sh = time.time() - t2

            # Differential Evolution
            t3 = time.time()
            de_res = run_scipy_diff_evo(p)
            t_de = time.time() - t3

            # Dual Annealing
            t4 = time.time()
            da_res = run_scipy_dual_annealing(p)
            t_da = time.time() - t4

        # Print row
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

        results[p.name] = {
            "my": my_res,
            "dreal": dreal_res,
            "shgo": shgo_res,
            "de": de_res,
            "da": da_res,
        }

    print("\nGenerating 3D plots...")
    for p in problems:
        visualize_problem_3d(p, results[p.name])


#   Main with argument parser

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--visuals", action="store_true",
                        help="Generate 3D surface visualizations")
    args = parser.parse_args()

    problems = [
        test_sanity_poly(), test_sanity_rational(),
        test_rational_bowl(), test_himmelblau_ratio(),
        test_split_islands(), test_singularity_edge(),
        test_pole_avoidance(), test_rational_valley(),
        test_positive_islands(), test_sparse_intersection()
    ]

    if args.visuals:
        run_with_visuals(problems)
    else:
        run_test_suite(problems)
