#!/usr/bin/env python3

import csv
import os
from collections import defaultdict

# =========================================================
# INPUT FILES
# =========================================================

EVENT_FILES = {
    "A": "A_lineages_ARG_carrier_cpg.csv",
    "B": "B_lineages_ARG_carrier_cpg.csv",
    "E": "E_lineages_ARG_carrier_cpg.csv",
    "H": "H_lineages_ARG_carrier_cpg.csv",
    "J": "J_lineages_ARG_carrier_cpg.csv",
    "L": "L_lineages_ARG_carrier_cpg.csv",
}

POOLED_FILE = "all_lineages_ARG_carrier_cpg.csv"

TABLE1_OUT = "Table1_ARG_host_event_summary.csv"
TABLE2_OUT = "Table2_ARG_host_temporal_matrix.csv"

EVENT_ORDER = ["A", "B", "E", "H", "J", "L"]


# =========================================================
# LOAD ONE ARG-HOST FILE
# =========================================================

def load_arg_host_file(filename):

    rows = []

    with open(filename, newline="") as fh:
        reader = csv.DictReader(fh)

        required_cols = {
            "ARG",
            "carrier",
            "rank_within_ARG_carrier",
            "lineage",
            "total_cpg"
        }

        if reader.fieldnames is None:
            raise ValueError(f"{filename} has no header")

        missing = required_cols - set(reader.fieldnames)

        if missing:
            raise ValueError(
                f"{filename} is missing required columns: {sorted(missing)}"
            )

        for row in reader:

            arg = row["ARG"].strip()
            carrier = row["carrier"].strip()
            lineage = row["lineage"].strip()

            try:
                total_cpg = float(row["total_cpg"])
            except (TypeError, ValueError):
                total_cpg = 0.0

            if not arg or not lineage:
                continue

            rows.append({
                "ARG": arg,
                "carrier": carrier,
                "lineage": lineage,
                "total_cpg": total_cpg
            })

    return rows


# =========================================================
# CALCULATE TABLE 1 SUMMARY METRICS
# =========================================================

def calculate_summary(label, rows):

    unique_args = set()
    unique_hosts = set()

    # ARG + lineage
    arg_host_pairs = set()

    # plasmid ARG + lineage
    plasmid_arg_host_pairs = set()

    unique_plasmid_args = set()
    unique_plasmid_hosts = set()

    for row in rows:

        arg = row["ARG"]
        carrier = row["carrier"]
        lineage = row["lineage"]

        unique_args.add(arg)
        unique_hosts.add(lineage)

        arg_host_pairs.add(
            (arg, lineage)
        )

        if "plasmid" in carrier.lower():

            plasmid_arg_host_pairs.add(
                (arg, lineage)
            )

            unique_plasmid_args.add(arg)
            unique_plasmid_hosts.add(lineage)

    n_arg_host = len(arg_host_pairs)
    n_plasmid_arg_host = len(plasmid_arg_host_pairs)

    if n_arg_host > 0:
        plasmid_percent = (
            n_plasmid_arg_host
            / n_arg_host
            * 100
        )
    else:
        plasmid_percent = 0.0

    return {
        "sampling_event": label,
        "unique_ARGs": len(unique_args),
        "unique_hosts": len(unique_hosts),
        "unique_ARG_host_associations": n_arg_host,
        "plasmid_ARG_host_associations": n_plasmid_arg_host,
        "plasmid_ARG_host_association_percent": round(
            plasmid_percent, 2
        ),
        "unique_plasmid_ARGs": len(unique_plasmid_args),
        "unique_plasmid_hosts": len(unique_plasmid_hosts)
    }


# =========================================================
# CHECK INPUT FILES
# =========================================================

for event, filename in EVENT_FILES.items():

    if not os.path.exists(filename):
        raise FileNotFoundError(
            f"Cannot find file for event {event}: {filename}"
        )

if not os.path.exists(POOLED_FILE):
    raise FileNotFoundError(
        f"Cannot find pooled file: {POOLED_FILE}"
    )


# =========================================================
# LOAD EVENT FILES
# =========================================================

event_data = {}

for event in EVENT_ORDER:

    filename = EVENT_FILES[event]

    print(f"Loading event {event}: {filename}")

    event_data[event] = load_arg_host_file(filename)


# =========================================================
# LOAD OVERALL FILE
# =========================================================

print(f"Loading Overall: {POOLED_FILE}")

overall_data = load_arg_host_file(POOLED_FILE)


# =========================================================
# TABLE 1
# =========================================================

summary_rows = []

for event in EVENT_ORDER:

    summary_rows.append(
        calculate_summary(
            event,
            event_data[event]
        )
    )

# Add pooled dataset
summary_rows.append(
    calculate_summary(
        "Overall",
        overall_data
    )
)


table1_fields = [
    "sampling_event",
    "unique_ARGs",
    "unique_hosts",
    "unique_ARG_host_associations",
    "plasmid_ARG_host_associations",
    "plasmid_ARG_host_association_percent",
    "unique_plasmid_ARGs",
    "unique_plasmid_hosts"
]


with open(TABLE1_OUT, "w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=table1_fields
    )

    writer.writeheader()
    writer.writerows(summary_rows)


print(f"\nSaved: {TABLE1_OUT}")


# =========================================================
# TABLE 2
#
# Key:
# ARG + carrier + lineage
# =========================================================

# key -> event -> cpg
temporal_data = defaultdict(
    lambda: defaultdict(float)
)


for event in EVENT_ORDER:

    for row in event_data[event]:

        key = (
            row["ARG"],
            row["carrier"],
            row["lineage"]
        )

        temporal_data[key][event] += row["total_cpg"]


# =========================================================
# PREPARE OVERALL / POOLED VALUES
# =========================================================

overall_values = defaultdict(float)

for row in overall_data:

    key = (
        row["ARG"],
        row["carrier"],
        row["lineage"]
    )

    overall_values[key] += row["total_cpg"]


# =========================================================
# INCLUDE KEYS THAT APPEAR ONLY IN OVERALL FILE
# =========================================================

all_keys = set(temporal_data.keys()) | set(overall_values.keys())


# =========================================================
# BUILD TABLE 2
# =========================================================

table2_rows = []

for key in all_keys:

    arg, carrier, lineage = key

    row_out = {
        "ARG": arg,
        "carrier": carrier,
        "lineage": lineage
    }

    sum_of_events = 0.0

    for event in EVENT_ORDER:

        value = temporal_data[key].get(
            event,
            0.0
        )

        row_out[event] = round(value, 4)

        sum_of_events += value


    overall = overall_values.get(
        key,
        0.0
    )

    difference = (
        sum_of_events
        - overall
    )


    row_out["Overall"] = round(
        overall,
        4
    )

    row_out["sum_of_events"] = round(
        sum_of_events,
        4
    )

    row_out["difference"] = round(
        difference,
        4
    )

    table2_rows.append(row_out)


# =========================================================
# SORT TABLE 2
#
# First by ARG
# then carrier
# then Overall cpg from highest to lowest
# =========================================================

table2_rows.sort(
    key=lambda x: (
        x["ARG"],
        x["carrier"],
        -x["Overall"]
    )
)


# =========================================================
# WRITE TABLE 2
# =========================================================

table2_fields = [
    "ARG",
    "carrier",
    "lineage",
    "A",
    "B",
    "E",
    "H",
    "J",
    "L",
    "Overall",
    "sum_of_events",
    "difference"
]


with open(TABLE2_OUT, "w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=table2_fields
    )

    writer.writeheader()
    writer.writerows(table2_rows)


print(f"Saved: {TABLE2_OUT}")


# =========================================================
# QC REPORT
# =========================================================

if table2_rows:

    max_difference = max(
        abs(row["difference"])
        for row in table2_rows
    )

else:
    max_difference = 0.0


print("\nFinished.")

print(
    f"Table 1 contains {len(summary_rows)} rows "
    f"(6 sampling events + Overall)."
)

print(
    f"Table 2 contains {len(table2_rows)} unique "
    f"ARG-carrier-host combinations."
)

print(
    f"Maximum absolute difference between "
    f"sum of six event cpg and Overall cpg: "
    f"{max_difference:.4f}"
)