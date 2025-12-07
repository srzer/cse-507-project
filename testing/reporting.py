# This file contains helper functions for reporting results.
import csv
from typing import List, Dict, Any

def results_to_markdown(results: List[Dict[str, Any]]) -> str:
    """Converts a list of result dictionaries to a markdown table."""
    if not results:
        return "No results to display."
    
    headers = results[0].keys()
    md = "| " + " | ".join(headers) + " |\n"
    md += "| " + " | ".join(["---"] * len(headers)) + " |\n"
    for res in results:
        md += "| " + " | ".join(str(v) for v in res.values()) + " |\n"
    return md

def save_results_to_csv(results: List[Dict[str, Any]], filename: str):
    """Saves a list of result dictionaries to a CSV file."""
    if not results:
        return
    with open(filename, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
