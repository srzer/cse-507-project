import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Define file paths
BASE_DIR = Path(__file__).resolve().parent.parent
INPUT_CSV = BASE_DIR / "report" / "log" / "comparison_summary.csv"
OUTPUT_DIR = BASE_DIR / "report" / "renders"


def create_config_label(row):
    """Creates a descriptive label for each algorithm configuration."""
    if row["algorithm"] == "BaselineMin":
        return "Baseline"

    # Handle cases where splitter or bounder might be NaN or other non-string types
    splitter = str(row["splitter"]) if pd.notna(row["splitter"]) else "NA"
    bounder = str(row["bounder"]) if pd.notna(row["bounder"]) else "NA"

    # Per user request, the 'Improved' algorithm is the only one being compared, so its name is omitted from the label.
    splitter = splitter.replace("SplitLongestSide", "Longest").replace(
        "SplitGradient", "Gradient"
    )
    bounder = bounder.replace("Bounds", "")

    return f"{splitter}-{bounder}"


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

        problem_df = problem_df.sort_values("config_label")

        plt.style.use("ggplot")
        fig, ax = plt.subplots(figsize=(12, 8))

        # Filter out NaN runtimes for plotting, but keep for annotation if possible
        plot_data = problem_df.dropna(subset=["runtime"])
        if plot_data.empty:
            print(f"Skipping plot for problem {problem}: No valid runtime data.")
            plt.close(fig)
            continue

        bars = ax.bar(plot_data["config_label"], plot_data["runtime"], color="skyblue")

        # Annotate bars with the bound value from the original problem_df
        # We need to ensure that the annotation matches the bar it's on.
        # This requires matching config_label between plot_data and problem_df
        for i, bar in enumerate(bars):
            config = plot_data["config_label"].iloc[i]
            # Find the original bound for this configuration from the full problem_df
            original_bound_row = problem_df[problem_df["config_label"] == config]
            if not original_bound_row.empty:
                bound_val = original_bound_row["bound"].iloc[0]
                if pd.notna(bound_val):
                    ax.text(
                        bar.get_x() + bar.get_width() / 2.0,
                        bar.get_height(),  # Place annotation on top of the bar
                        f"{bound_val:.3f}",
                        ha="center",
                        va="bottom",
                        fontsize=9,
                        rotation=0,
                    )

        ax.set_yscale("log")
        ax.set_ylabel("Runtime (seconds, log scale)")
        ax.set_title(f"Performance on Problem: {problem}")

        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()

        # Save the figure
        output_path_pdf = OUTPUT_DIR / f"{problem.replace(' ', '_')}_performance.pdf"
        output_path_svg = OUTPUT_DIR / f"{problem.replace(' ', '_')}_performance.svg"
        plt.savefig(output_path_pdf, format="pdf")
        plt.savefig(output_path_svg, format="svg")
        plt.close(fig)
        print(f"Saved plots: {output_path_pdf}, {output_path_svg}")


def generate_aggregate_plot(df):
    """
    Generates a box plot of normalized runtime across all problems.
    Expects df to already have 'config_label' and 'normalized_runtime' columns.
    """
    # Filter out baseline and rows where normalization wasn't possible
    plot_df = (
        df[df["algorithm"] != "BaselineMin"]
        .dropna(subset=["normalized_runtime"])
        .copy()
    )

    if plot_df.empty:
        print("No valid data for aggregate plot after normalization.")
        return

    plt.style.use("ggplot")
    fig, ax = plt.subplots(figsize=(12, 8))

    configs = plot_df["config_label"].unique()
    configs.sort()

    data_to_plot = [
        plot_df[plot_df["config_label"] == config]["normalized_runtime"]
        for config in configs
    ]

    # Showfliers=False to hide extreme outliers, which can distort the plot for log scale
    ax.boxplot(data_to_plot, tick_labels=configs, vert=True, showfliers=False)

    # Add a line for the baseline performance (normalized runtime of 1.0)
    ax.axhline(
        y=1.0, color="r", linestyle="--", linewidth=2, label="Baseline Performance"
    )

    ax.set_yscale("log")
    ax.set_ylabel("Normalized Runtime (log scale) (runtime / baseline_runtime)")
    ax.set_title("Aggregate Performance Across All Problems")
    ax.legend()

    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    # Save the figure
    output_path_pdf = OUTPUT_DIR / "aggregate_normalized_runtime.pdf"
    output_path_svg = OUTPUT_DIR / "aggregate_normalized_runtime.svg"
    plt.savefig(output_path_pdf, format="pdf")
    plt.savefig(output_path_svg, format="svg")
    plt.close(fig)
    print(f"Saved plots: {output_path_pdf}, {output_path_svg}")


def generate_aggregate_bound_plot(df):
    """
    Generates a box plot of the difference between algorithm bounds and baseline bounds.
    """
    # Isolate baseline bounds for mapping
    baseline_bounds_map = (
        df[df["algorithm"] == "BaselineMin"].set_index("problem")["bound"].to_dict()
    )

    # Work with the improved algorithm data
    plot_df = df[df["algorithm"] != "BaselineMin"].copy()
    plot_df["baseline_bound"] = plot_df["problem"].map(baseline_bounds_map)

    # Calculate the difference; drop rows where a difference can't be computed
    plot_df = plot_df.dropna(subset=["bound", "baseline_bound"])
    plot_df["bound_difference"] = plot_df["bound"] - plot_df["baseline_bound"]

    # Per user request, filter out 'Gradient-Affine' to improve plot readability
    plot_df = plot_df[plot_df["config_label"] != "Gradient-Affine"].copy()

    if plot_df.empty:
        print("No valid data for aggregate bound plot after processing.")
        return

    plt.style.use("ggplot")
    fig, ax = plt.subplots(figsize=(12, 8))

    configs = plot_df["config_label"].unique()
    configs.sort()

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
    plt.tight_layout()

    # Save the figure
    output_path_pdf = OUTPUT_DIR / "aggregate_bound_difference.pdf"
    output_path_svg = OUTPUT_DIR / "aggregate_bound_difference.svg"
    plt.savefig(output_path_pdf, format="pdf")
    plt.savefig(output_path_svg, format="svg")
    plt.close(fig)
    print(f"Saved plots: {output_path_pdf}, {output_path_svg}")


def main():
    print("Available matplotlib styles:", plt.style.available)
    """Main function to load data, process, export, and generate plots."""
    print(f"Reading data from: {INPUT_CSV}")
    try:
        df = pd.read_csv(INPUT_CSV)
    except FileNotFoundError:
        print(f"Error: Input file not found at {INPUT_CSV}")
        return

    # Convert relevant columns to numeric, coercing errors
    df["runtime"] = pd.to_numeric(df["runtime"], errors="coerce")
    df["bound"] = pd.to_numeric(df["bound"], errors="coerce")

    # Per user request, filter out the original algorithm to focus on the improved one vs baseline
    df = df[df["algorithm"] != "GlobalMinBranchAndBound"].copy()

    # Add config_label column
    df["config_label"] = df.apply(create_config_label, axis=1)

    # Calculate normalized runtime for non-baseline entries and add to df
    baseline_df = df[df["algorithm"] == "BaselineMin"].copy()

    # Store baseline runtimes in a dict for easy lookup
    baseline_runtimes_map = baseline_df.set_index("problem")["runtime"].to_dict()

    # Add 'baseline_runtime_for_problem' column to all rows
    df["baseline_runtime_for_problem"] = df["problem"].map(baseline_runtimes_map)

    # Calculate normalized_runtime only for non-baseline rows that have valid runtime and baseline_runtime
    df["normalized_runtime"] = np.nan  # Initialize column
    mask_for_normalization = (
        (df["algorithm"] != "BaselineMin")
        & df["runtime"].notna()
        & df["baseline_runtime_for_problem"].notna()
    )
    df.loc[mask_for_normalization, "normalized_runtime"] = (
        df.loc[mask_for_normalization, "runtime"]
        / df.loc[mask_for_normalization, "baseline_runtime_for_problem"]
    )

    # Export processed CSV
    processed_csv_path = OUTPUT_DIR / "processed_comparison_summary.csv"
    df.to_csv(processed_csv_path, index=False)
    print(f"Exported processed data to: {processed_csv_path}")

    print("\nGenerating per-problem performance plots...")
    generate_per_problem_plots(
        df.copy()
    )  # Pass a copy to prevent unintended modifications

    print("\nGenerating aggregate performance plot...")
    generate_aggregate_plot(df.copy())  # Pass a copy

    print("\nGenerating aggregate bound difference plot...")
    generate_aggregate_bound_plot(df.copy())

    print("\nAll plots and processed data generated successfully.")


if __name__ == "__main__":
    main()
