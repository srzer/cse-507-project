import csv
from typing import List, Dict, Any


# convert a list of result dictionaries to a markdown table
def results_to_markdown(results: List[Dict[str, Any]]) -> str:
    if not results:
        return "No results to display."

    headers = results[0].keys()
    md = "| " + " | ".join(headers) + " |\n"
    md += "| " + " | ".join(["---"] * len(headers)) + " |\n"
    for res in results:
        md += "| " + " | ".join(str(v) for v in res.values()) + " |\n"
    return md


# save list of result dicts to csv file
def save_results_to_csv(results: List[Dict[str, Any]], filename: str):
    if not results:
        return
    with open(filename, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
