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

OUT_FILE = "ARG_subtype_drug_class_summary.csv"


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

                    files.append(
                        os.path.join(dp, filename)
                    )

    return files


# =========================================================
# MAIN
# =========================================================

def main():

    files = find_tsv_files(ARGO_DIRS)

    print(f"Found {len(files)} .sarg.tsv files.")

    if not files:
        raise FileNotFoundError(
            "No .sarg.tsv files were found."
        )

    # subtype -> set of drug classes
    arg_classes = defaultdict(set)

    for filepath in files:

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

            required = {"subtype", "type"}

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
                drug_class = row["type"].strip()

                if not arg:
                    continue

                if not drug_class:
                    drug_class = "unclassified"

                arg_classes[arg].add(drug_class)


    # =====================================================
    # WRITE OUTPUT
    # =====================================================

    with open(
        OUT_FILE,
        "w",
        newline="",
        encoding="utf-8"
    ) as fh:

        writer = csv.writer(fh)

        writer.writerow([
            "ARG",
            "Drug_class",
            "Number_of_drug_classes"
        ])

        for arg in sorted(arg_classes):

            classes = sorted(arg_classes[arg])

            writer.writerow([
                arg,
                "; ".join(classes),
                len(classes)
            ])


    # =====================================================
    # QC
    # =====================================================

    print()
    print(f"Unique ARG subtypes: {len(arg_classes)}")

    multi_class = {
        arg: classes
        for arg, classes in arg_classes.items()
        if len(classes) > 1
    }

    print(
        f"ARGs assigned to >1 drug class: "
        f"{len(multi_class)}"
    )

    if multi_class:

        print()
        print("ARGs with multiple drug-class assignments:")

        for arg, classes in sorted(
            multi_class.items()
        ):
            print(
                f"  {arg}: "
                + ", ".join(sorted(classes))
            )

    print()
    print(f"Output file: {OUT_FILE}")


if __name__ == "__main__":
    main()