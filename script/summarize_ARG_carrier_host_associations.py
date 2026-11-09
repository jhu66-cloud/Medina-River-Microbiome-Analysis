#!/usr/bin/env python3

import os
import csv
from collections import defaultdict

# =========================================================
# INPUT
# =========================================================

ARGO_DIRS = [
    "./Argo_results"
]

OUT_FILE = "Table1_ARG_carrier_host_event_summary.csv"

VALID_EVENTS = {"A", "B", "E", "H", "J", "L"}


# =========================================================
# FIND ALL .sarg.tsv FILES
# =========================================================

def find_tsv_files(argo_dirs):

    files = []

    for base in argo_dirs:

        if not os.path.isdir(base):
            print(f"Warning: directory not found: {base}")
            continue

        for sub in sorted(os.listdir(base)):

            dp = os.path.join(base, sub)

            if not os.path.isdir(dp):
                continue

            for filename in sorted(os.listdir(dp)):

                if filename.endswith(".sarg.tsv"):

                    sample = (
                        filename
                        .replace(".sarg.tsv", "")
                        .replace("chopped_", "")
                    )

                    event = sample[-1].upper()

                    if event not in VALID_EVENTS:
                        print(
                            f"Warning: could not determine event "
                            f"for sample {sample}"
                        )
                        continue

                    files.append({
                        "sample": sample,
                        "event": event,
                        "path": os.path.join(dp, filename)
                    })

    return files


# =========================================================
# INITIALIZE EVENT DATA
# =========================================================

def new_event_data():

    return {
        "args": set(),
        "hosts": set(),
        "associations": set(),
        "plasmid_associations": set(),
        "plasmid_args": set(),
        "plasmid_hosts": set()
    }


# =========================================================
# MAIN
# =========================================================

def main():

    files = find_tsv_files(ARGO_DIRS)

    print(f"Found {len(files)} .sarg.tsv files.")

    if len(files) == 0:
        raise FileNotFoundError(
            "No .sarg.tsv files were found."
        )

    # One set of metrics for each event
    event_data = {
        event: new_event_data()
        for event in VALID_EVENTS
    }

    # Pooled data across all six events
    pooled = new_event_data()


    # =====================================================
    # READ FILES
    # =====================================================

    for item in files:

        event = item["event"]
        filepath = item["path"]

        with open(
            filepath,
            "r",
            newline="",
            encoding="utf-8-sig"
        ) as fh:

            reader = csv.DictReader(
                fh,
                delimiter="\t"
            )

            required = {
                "subtype",
                "carrier",
                "lineage"
            }

            if reader.fieldnames is None:
                raise ValueError(
                    f"No header found in {filepath}"
                )

            missing = required - set(reader.fieldnames)

            if missing:
                raise ValueError(
                    f"{filepath} missing columns: "
                    + ", ".join(sorted(missing))
                )

            for row in reader:

                arg = row["subtype"].strip()
                carrier = row["carrier"].strip().lower()
                lineage = row["lineage"].strip()

                if not arg:
                    continue

                if not lineage:
                    lineage = "unclassified"

                # -----------------------------------------
                # Definition of one unique association:
                #
                # ARG + carrier + lineage
                # -----------------------------------------

                association = (
                    arg,
                    carrier,
                    lineage
                )

                # =========================================
                # EVENT LEVEL
                # =========================================

                event_data[event]["args"].add(arg)

                event_data[event]["hosts"].add(lineage)

                event_data[event]["associations"].add(
                    association
                )

                if carrier == "plasmid":

                    event_data[event][
                        "plasmid_associations"
                    ].add(association)

                    event_data[event][
                        "plasmid_args"
                    ].add(arg)

                    event_data[event][
                        "plasmid_hosts"
                    ].add(lineage)


                # =========================================
                # POOLED LEVEL
                # =========================================

                pooled["args"].add(arg)

                pooled["hosts"].add(lineage)

                pooled["associations"].add(
                    association
                )

                if carrier == "plasmid":

                    pooled[
                        "plasmid_associations"
                    ].add(association)

                    pooled[
                        "plasmid_args"
                    ].add(arg)

                    pooled[
                        "plasmid_hosts"
                    ].add(lineage)


    # =====================================================
    # IDENTIFY EVENT-SPECIFIC ASSOCIATIONS
    # =====================================================

    event_specific = {}

    for event in ["A", "B", "E", "H", "J", "L"]:

        other_events = [
            e
            for e in ["A", "B", "E", "H", "J", "L"]
            if e != event
        ]

        associations_other_events = set()

        for other in other_events:

            associations_other_events.update(
                event_data[other]["associations"]
            )

        event_specific[event] = (
            event_data[event]["associations"]
            - associations_other_events
        )


    # =====================================================
    # CREATE SUMMARY
    # =====================================================

    rows = []

    pooled_total = len(
        pooled["associations"]
    )

    for event in ["A", "B", "E", "H", "J", "L"]:

        d = event_data[event]

        total_assoc = len(
            d["associations"]
        )

        plasmid_assoc = len(
            d["plasmid_associations"]
        )

        if total_assoc > 0:

            plasmid_percent = (
                plasmid_assoc
                / total_assoc
                * 100
            )

        else:

            plasmid_percent = 0.0

        specific_count = len(
            event_specific[event]
        )

        if pooled_total > 0:

            specific_percent = (
                specific_count
                / pooled_total
                * 100
            )

        else:

            specific_percent = 0.0


        rows.append({
            "Sampling_event":
                event,

            "Unique_ARGs":
                len(d["args"]),

            "Unique_hosts":
                len(d["hosts"]),

            "Unique_ARG_carrier_host_associations":
                total_assoc,

            "Plasmid_ARG_carrier_host_associations":
                plasmid_assoc,

            "Plasmid_ARG_carrier_host_association_percent":
                plasmid_percent,

            "Unique_plasmid_ARGs":
                len(d["plasmid_args"]),

            "Unique_plasmid_hosts":
                len(d["plasmid_hosts"]),

            "Event_specific_ARG_carrier_host_associations":
                specific_count,

            "Percent_of_overall_associations_unique_to_event":
                specific_percent
        })


    # =====================================================
    # ADD ALL EVENTS ROW
    # =====================================================

    pooled_plasmid_assoc = len(
        pooled["plasmid_associations"]
    )

    if pooled_total > 0:

        pooled_plasmid_percent = (
            pooled_plasmid_assoc
            / pooled_total
            * 100
        )

    else:

        pooled_plasmid_percent = 0.0


    rows.append({
        "Sampling_event":
            "All events",

        "Unique_ARGs":
            len(pooled["args"]),

        "Unique_hosts":
            len(pooled["hosts"]),

        "Unique_ARG_carrier_host_associations":
            pooled_total,

        "Plasmid_ARG_carrier_host_associations":
            pooled_plasmid_assoc,

        "Plasmid_ARG_carrier_host_association_percent":
            pooled_plasmid_percent,

        "Unique_plasmid_ARGs":
            len(pooled["plasmid_args"]),

        "Unique_plasmid_hosts":
            len(pooled["plasmid_hosts"]),

        "Event_specific_ARG_carrier_host_associations":
            "",

        "Percent_of_overall_associations_unique_to_event":
            ""
    })


    # =====================================================
    # WRITE CSV
    # =====================================================

    fieldnames = [
        "Sampling_event",
        "Unique_ARGs",
        "Unique_hosts",
        "Unique_ARG_carrier_host_associations",
        "Plasmid_ARG_carrier_host_associations",
        "Plasmid_ARG_carrier_host_association_percent",
        "Unique_plasmid_ARGs",
        "Unique_plasmid_hosts",
        "Event_specific_ARG_carrier_host_associations",
        "Percent_of_overall_associations_unique_to_event"
    ]

    with open(
        OUT_FILE,
        "w",
        newline="",
        encoding="utf-8"
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=fieldnames
        )

        writer.writeheader()

        for row in rows:

            output = row.copy()

            if (
                output[
                    "Plasmid_ARG_carrier_host_association_percent"
                ] != ""
            ):
                output[
                    "Plasmid_ARG_carrier_host_association_percent"
                ] = (
                    f"{output['Plasmid_ARG_carrier_host_association_percent']:.2f}"
                )

            if (
                output[
                    "Percent_of_overall_associations_unique_to_event"
                ] != ""
            ):
                output[
                    "Percent_of_overall_associations_unique_to_event"
                ] = (
                    f"{output['Percent_of_overall_associations_unique_to_event']:.2f}"
                )

            writer.writerow(output)


    # =====================================================
    # QC
    # =====================================================

    print()
    print("=" * 65)

    print(
        f"Pooled unique ARG-carrier-host associations: "
        f"{pooled_total}"
    )

    print(
        f"Pooled plasmid-associated combinations: "
        f"{pooled_plasmid_assoc}"
    )

    print(
        f"Unique ARGs: "
        f"{len(pooled['args'])}"
    )

    print(
        f"Unique hosts: "
        f"{len(pooled['hosts'])}"
    )

    print(
        f"Output file: {OUT_FILE}"
    )

    print("=" * 65)


if __name__ == "__main__":
    main()