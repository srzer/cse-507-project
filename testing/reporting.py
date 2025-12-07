import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Define file paths
BASE_DIR = Path(__file__).resolve().parent.parent
SMT_CSV_PATH = BASE_DIR / "report" / "log" / "comparison_summary.csv"
NUMERIC_CSV_PATH = BASE_DIR / "report" / "log" / "numeric_comparison.csv"
OUTPUT_DIR = BASE_DIR / "report" / "renders"

# Define a stable sort order for algorithms
ALGO_SORT_ORDER = [
    "Baseline(dReal)",
    "Longest-Affine",
    "Longest-Bernstein",
    "Gradient-Affine",
    "Gradient-Bernstein",
    "SHGO",
    "Diff. Evo.",
    "Dual Annealing",
]


def create_config_label(row):
    """Creates a descriptive label for each algorithm configuration."""
    algo_name = row["algorithm"]
    if algo_name in ["SHGO", "Diff. Evo.", "Dual Annealing"]:
        return algo_name
    if algo_name == "BaselineMin":
        return "Baseline(dReal)"
    if "Improved" in algo_name:
        splitter = str(row.get("splitter", "NA")).replace("SplitLongestSide", "Longest").replace("SplitGradient", "Gradient")
        bounder = str(row.get("bounder", "NA")).replace("Bounds", "")
        return f"{splitter}-{bounder}"
    return algo_name


def generate_per_problem_plots(df):
    """
    Generates a bar chart for each problem, showing runtime and annotating the bound.
    Expects df to have a 'config_label' column.
    """
    problems = df["problem"].unique()

    for problem in problems:
        problem_df = df[df["problem"] == problem].copy()
        if problem_df.empty:
            continue

        # Sort the bars for consistent plotting order
        problem_df['config_label'] = pd.Categorical(problem_df['config_label'], categories=ALGO_SORT_ORDER, ordered=True)
        problem_df = problem_df.sort_values("config_label")

        plt.style.use("ggplot")
        fig, ax = plt.subplots(figsize=(12, 8))

        plot_data = problem_df.dropna(subset=["runtime"])
        if plot_data.empty:
            print(f"Skipping plot for problem {problem}: No valid runtime data.")
            plt.close(fig)
            continue

        bars = ax.bar(plot_data["config_label"], plot_data["runtime"], color="skyblue")

        for i, bar in enumerate(bars):
            bound_val = plot_data["bound"].iloc[i]
            if pd.notna(bound_val):
                ax.text(
                    bar.get_x() + bar.get_width() / 2.0,
                    bar.get_height(),
                    f"{bound_val:.3f}",
                    ha="center",
                    va="bottom",
                    fontsize=9,
                )

        ax.set_yscale("log")
        ax.set_ylabel("Runtime (seconds, log scale)")
        ax.set_title(f"Performance on Problem: {problem}")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()

        output_path_svg = OUTPUT_DIR / f"full_comparison_{problem.replace(' ', '_')}.svg"
        plt.savefig(output_path_svg, format="svg")
        plt.close(fig)
        print(f"Saved combined plot: {output_path_svg}")


def generate_aggregate_plot(df):
    """
    Generates a box plot of normalized runtime across all problems.
    """
    plot_df = df[df["config_label"] != "Baseline(dReal)"].dropna(subset=["normalized_runtime"]).copy()
    if plot_df.empty:
        print("No valid data for aggregate plot after normalization.")
        return

    plt.style.use("ggplot")
    fig, ax = plt.subplots(figsize=(12, 8))
    
    plot_df['config_label'] = pd.Categorical(plot_df['config_label'], categories=[c for c in ALGO_SORT_ORDER if c != "Baseline(dReal)"], ordered=True)
    plot_df = plot_df.sort_values("config_label")
    configs = plot_df["config_label"].unique()

    data_to_plot = [plot_df[plot_df["config_label"] == config]["normalized_runtime"] for config in configs]

    ax.boxplot(data_to_plot, tick_labels=configs, vert=True, showfliers=False)
    ax.axhline(y=1.0, color="r", linestyle="--", linewidth=2, label="Baseline Performance")
    ax.set_yscale("log")
    ax.set_ylabel("Normalized Runtime (log scale) (runtime / baseline_runtime)")
    ax.set_title("Aggregate Performance Across All Problems")
    ax.legend()
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    output_path_svg = OUTPUT_DIR / "full_aggregate_normalized_runtime.svg"
    plt.savefig(output_path_svg, format="svg")
    plt.close(fig)
    print(f"Saved combined plot: {output_path_svg}")


def generate_aggregate_bound_plot(df):
    """
    Generates a box plot of the difference between algorithm bounds and baseline bounds.
    """
    baseline_bounds_map = df[df["config_label"] == "Baseline(dReal)"].set_index("problem")["bound"].to_dict()
    plot_df = df[df["config_label"] != "Baseline(dReal)"].copy()
    plot_df["baseline_bound"] = plot_df["problem"].map(baseline_bounds_map)
    plot_df = plot_df.dropna(subset=["bound", "baseline_bound"])
    plot_df["bound_difference"] = plot_df["bound"] - plot_df["baseline_bound"]
    
    # Per user request, filter out algorithms with high variance to improve plot readability
    removed_for_clarity = ["Gradient-Affine", "Dual Annealing"]
    plot_df = plot_df[~plot_df["config_label"].isin(removed_for_clarity)].copy()

    if plot_df.empty:
        print("No valid data for aggregate bound plot after processing.")
        return

    plt.style.use("ggplot")
    fig, ax = plt.subplots(figsize=(12, 8))

    plot_df['config_label'] = pd.Categorical(plot_df['config_label'], categories=[c for c in ALGO_SORT_ORDER if c not in ["Baseline(dReal)"] + removed_for_clarity], ordered=True)
    plot_df = plot_df.sort_values("config_label")
    configs = plot_df["config_label"].unique()

    data_to_plot = [
        plot_df[plot_df["config_label"] == config]["bound_difference"]
        for config in configs
    ]

    ax.boxplot(data_to_plot, tick_labels=configs, vert=True, showfliers=False)

    # Add a line for the baseline performance (difference of 0)
    ax.axhline(y=0.0, color="r", linestyle="--", linewidth=2, label="Baseline Bound")

    ax.set_ylabel("Bound Difference (bound - baseline_bound)")
    ax.set_title("Aggregate Bound Closeness to Baseline Across All Problems")
    ax.legend()

    plt.xticks(rotation=45, ha="right")

    # Add annotation for removed data
    fig.text(0.02, 0.02, f"Removed for clarity: {', '.join(removed_for_clarity)}",
             ha='left', va='bottom', fontsize=10, color='gray')
    
    plt.tight_layout(rect=[0, 0.05, 1, 1])  # Adjust layout to make space

    # Save the figure
    output_path_svg = OUTPUT_DIR / "full_aggregate_bound_difference.svg"
    plt.savefig(output_path_svg, format="svg")
    plt.close(fig)
    print(f"Saved combined plot: {output_path_svg}")


def main():
    """Main function to load data from two sources and generate combined plots."""
    try:
        df_smt = pd.read_csv(SMT_CSV_PATH)
        df_numeric = pd.read_csv(NUMERIC_CSV_PATH)
    except FileNotFoundError as e:
        print(f"Error: Input CSV not found. {e}")
        return

    # Filter for only the numeric solvers from the second CSV
    numeric_solvers = ["SHGO", "Diff. Evo.", "Dual Annealing"]
    df_numeric = df_numeric[df_numeric["algorithm"].isin(numeric_solvers)].copy()

    # Combine the dataframes
    df_combined = pd.concat([df_smt, df_numeric], ignore_index=True)

    # Convert columns to numeric, coercing errors
    df_combined["runtime"] = pd.to_numeric(df_combined["runtime"], errors="coerce")
    df_combined["bound"] = pd.to_numeric(df_combined["bound"], errors="coerce")

    # Filter out old SMT algo as per previous request
    df_combined = df_combined[df_combined["algorithm"] != "GlobalMinBranchAndBound"].copy()

    # Create the unified configuration label
    df_combined["config_label"] = df_combined.apply(create_config_label, axis=1)

    # --- Calculate metrics relative to the baseline ---
    baseline_df = df_combined[df_combined["config_label"] == "Baseline(dReal)"].copy()
    baseline_runtimes_map = baseline_df.set_index("problem")["runtime"].to_dict()
    df_combined["baseline_runtime_for_problem"] = df_combined["problem"].map(baseline_runtimes_map)
    
    df_combined["normalized_runtime"] = np.nan
    mask = (df_combined["config_label"] != "Baseline(dReal)") & df_combined["runtime"].notna() & df_combined["baseline_runtime_for_problem"].notna()
    df_combined.loc[mask, "normalized_runtime"] = df_combined.loc[mask, "runtime"] / df_combined.loc[mask, "baseline_runtime_for_problem"]

    # --- Generate Plots ---
    print("\nGenerating per-problem combined performance plots...")
    generate_per_problem_plots(df_combined.copy())
    
    print("\nGenerating aggregate combined performance plot...")
    generate_aggregate_plot(df_combined.copy())
    
    print("\nGenerating aggregate combined bound difference plot...")
    generate_aggregate_bound_plot(df_combined.copy())
    
    print("\nAll combined plots generated successfully.")


if __name__ == "__main__":
    main()

