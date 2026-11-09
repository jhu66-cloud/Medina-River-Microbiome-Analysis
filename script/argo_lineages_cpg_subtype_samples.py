#!/usr/bin/env python3

import os
import csv

# =========================================================
# INPUT DIRECTORIES
# =========================================================

ARGO_DIRS = [
    "./Argo_results"
]

OUT_FILE = "ARG_host_summary_54_samples.csv"

VALID_EVENTS = {"A", "B", "E", "H", "J", "L"}


# =========================================================
# FIND SARG TSV FILES
# =========================================================

def find_tsv_files(argo_dirs):

    files = []

    for base in argo_dirs:

        if not os.path.isdir(base):
            print(f"Warning: directory not found: {base}")
            continue

        # Each sample has its own subdirectory
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

                    files.append({
                        "sample": sample,
                        "path": os.path.join(dp, filename)
                    })

    return files


# =========================================================
# PARSE SITE AND SAMPLING EVENT
# =========================================================

def parse_site_event(sample):

    if len(sample) > 1 and sample[-1].upper() in VALID_EVENTS:

        event = sample[-1].upper()
        site = sample[:-1]

    else:

        event = ""
        site = sample

    return site, event


# =========================================================
# SUMMARIZE ONE SAMPLE
# =========================================================

def summarize_sample(filepath):

    unique_args = set()
    unique_hosts = set()
    arg_host_associations = set()

    plasmid_arg_host_associations = set()
    unique_plasmid_args = set()
    unique_plasmid_hosts = set()

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

        required_columns = {
            "subtype",
            "lineage",
            "carrier"
        }

        if reader.fieldnames is None:
            raise ValueError(
                f"No header found in {filepath}"
            )

        missing = (
            required_columns
            - set(reader.fieldnames)
        )

        if missing:
            raise ValueError(
                f"{filepath} is missing columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:

            arg = row["subtype"].strip()
            lineage = row["lineage"].strip()
            carrier = row["carrier"].strip().lower()

            # Skip rows without ARG name
            if not arg:
                continue

            # Keep unclassified as a lineage category
            if not lineage:
                lineage = "unclassified"

            # ---------------------------------------------
            # ALL ARG-HOST METRICS
            # ---------------------------------------------

            unique_args.add(arg)

            unique_hosts.add(lineage)

            # Same ARG + same host counted once
            # regardless of chromosome/plasmid
            arg_host_associations.add(
                (arg, lineage)
            )

            # ---------------------------------------------
            # PLASMID-ASSOCIATED METRICS
            # ---------------------------------------------

            if carrier == "plasmid":

                plasmid_arg_host_associations.add(
                    (arg, lineage)
                )

                unique_plasmid_args.add(arg)

                unique_plasmid_hosts.add(lineage)

    # =====================================================
    # CALCULATE METRICS
    # =====================================================

    n_unique_args = len(
        unique_args
    )

    n_unique_hosts = len(
        unique_hosts
    )

    n_arg_host = len(
        arg_host_associations
    )

    n_plasmid_arg_host = len(
        plasmid_arg_host_associations
    )

    n_unique_plasmid_args = len(
        unique_plasmid_args
    )

    n_unique_plasmid_hosts = len(
        unique_plasmid_hosts
    )

    if n_arg_host > 0:

        plasmid_percent = (
            n_plasmid_arg_host
            / n_arg_host
            * 100
        )

    else:

        plasmid_percent = 0.0

    return {
        "Unique_ARGs":
            n_unique_args,

        "Unique_hosts":
            n_unique_hosts,

        "Unique_ARG_host_associations":
            n_arg_host,

        "Plasmid_ARG_host_associations":
            n_plasmid_arg_host,

        "Plasmid_ARG_host_association_percent":
            plasmid_percent,

        "Unique_plasmid_ARGs":
            n_unique_plasmid_args,

        "Unique_plasmid_hosts":
            n_unique_plasmid_hosts
    }


# =========================================================
# MAIN
# =========================================================

def main():

    files = find_tsv_files(
        ARGO_DIRS
    )

    print(
        f"Found {len(files)} .sarg.tsv files."
    )

    if len(files) == 0:

        raise FileNotFoundError(
            "No .sarg.tsv files were found."
        )

    results = []

    for item in files:

        sample = item["sample"]
        filepath = item["path"]

        site, event = parse_site_event(
            sample
        )

        metrics = summarize_sample(
            filepath
        )

        row = {
            "Sample": sample,
            "Site": site,
            "Sampling_event": event,

            "Unique_ARGs":
                metrics["Unique_ARGs"],

            "Unique_hosts":
                metrics["Unique_hosts"],

            "Unique_ARG_host_associations":
                metrics[
                    "Unique_ARG_host_associations"
                ],

            "Plasmid_ARG_host_associations":
                metrics[
                    "Plasmid_ARG_host_associations"
                ],

            "Plasmid_ARG_host_association_percent":
                metrics[
                    "Plasmid_ARG_host_association_percent"
                ],

            "Unique_plasmid_ARGs":
                metrics[
                    "Unique_plasmid_ARGs"
                ],

            "Unique_plasmid_hosts":
                metrics[
                    "Unique_plasmid_hosts"
                ]
        }

        results.append(row)

        print(
            f"{sample}: "
            f"{metrics['Unique_ARGs']} ARGs, "
            f"{metrics['Unique_hosts']} hosts, "
            f"{metrics['Unique_ARG_host_associations']} "
            f"ARG-host associations, "
            f"{metrics['Plasmid_ARG_host_associations']} "
            f"plasmid ARG-host associations"
        )


    # =====================================================
    # SORT OUTPUT
    # =====================================================

    event_order = {
        "A": 1,
        "B": 2,
        "E": 3,
        "H": 4,
        "J": 5,
        "L": 6
    }

    results.sort(
        key=lambda x: (
            event_order.get(
                x["Sampling_event"],
                99
            ),
            x["Site"]
        )
    )


    # =====================================================
    # WRITE CSV
    # =====================================================

    fieldnames = [
        "Sample",
        "Site",
        "Sampling_event",
        "Unique_ARGs",
        "Unique_hosts",
        "Unique_ARG_host_associations",
        "Plasmid_ARG_host_associations",
        "Plasmid_ARG_host_association_percent",
        "Unique_plasmid_ARGs",
        "Unique_plasmid_hosts"
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

        for row in results:

            output_row = row.copy()

            output_row[
                "Plasmid_ARG_host_association_percent"
            ] = (
                f"{row['Plasmid_ARG_host_association_percent']:.2f}"
            )

            writer.writerow(
                output_row
            )


    # =====================================================
    # FINAL QC
    # =====================================================

    print()
    print("=" * 60)

    print(
        f"Samples summarized: {len(results)}"
    )

    print(
        f"Output file: {OUT_FILE}"
    )

    if len(results) == 54:

        print(
            "QC passed: 54 samples found."
        )

    else:

        print(
            f"WARNING: expected 54 samples, "
            f"but found {len(results)}."
        )

    print("=" * 60)


if __name__ == "__main__":
    main()