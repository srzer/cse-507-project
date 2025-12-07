import csv
from typing import TypedDict


# log entry, attempt to add writer monad
class LogEntry(TypedDict):
    iteration: int
    action: str
    box: str
    bound: float
    volume: float
    notes: str


def save_logs_to_csv(logs: list[LogEntry], filename: str):
    if not logs:
        return
    with open(filename, "w", newline="") as csvfile:
        fieldnames = logs[0].keys()
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(logs)


def logs_to_markdown(logs: list[LogEntry]) -> str:
    if not logs:
        return "No logs were produced."
    header = logs[0].keys()
    md = "| " + " | ".join(header) + " |\n"
    md += "| " + " | ".join(["---"] * len(header)) + " |\n"
    for entry in logs:
        md += "| " + " | ".join(str(v) for v in entry.values()) + " |\n"
    return md
