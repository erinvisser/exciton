#!/usr/bin/env python3
"""
Compare two YANDF-0.4 internal transition rate data blocks.

Output shows two tables:
  1. Absolute difference |file1 - file2| for each lambda column.
  2. Relative error |file1 - file2| / max(|file1|,|file2|) as percentage.

Usage:
    python3 compare_lambda.py <file1> <file2> <output> [tolerance]

Default tolerance: 1e-4 (absolute) / 1% (relative)
"""

import sys
from datetime import datetime

COLUMN_NAMES = [
    "p(p)", "h(p)", "p(n)", "h(n)",
    "lambdapiplus", "lambdanuplus", "lambdapinu", "lambdanupi",
]
COLUMN_UNITS = ["", "", "", "", "sec^-1", "sec^-1", "sec^-1", "sec^-1"]

REL_COLUMN_NAMES = [
    "p(p)", "h(p)", "p(n)", "h(n)",
    "lambdapiplus", "lambdanuplus", "lambdapinu", "lambdanupi",
]
REL_COLUMN_UNITS = ["", "", "", "", "%", "%", "%", "%"]


def parse_file(path):
    with open(path) as f:
        lines = f.readlines()

    rows = []

    for line in lines:
        raw = line.strip()
        if not raw:
            continue
        if raw.startswith("#"):
            continue
        if raw.startswith("##"):
            continue
        parts = raw.split()
        if len(parts) == 8:
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                pass

    return rows


def write_header(out, file1, rows1, file2, rows2, tol):
    out.write("# lambda comparison report\n")
    out.write(f"# generated: {datetime.now().isoformat()}\n")
    out.write("#\n")
    out.write(f"# file1: {file1}\n")
    out.write(
        f"#   internal transition rate: {len(rows1)} rows, "
        f"{len(rows1[0]) if rows1 else 0} cols\n"
    )
    out.write("#\n")
    out.write(f"# file2: {file2}\n")
    out.write(
        f"#   internal transition rate: {len(rows2)} rows, "
        f"{len(rows2[0]) if rows2 else 0} cols\n"
    )
    out.write("#\n")
    out.write("# Description:\n")
    out.write("#   Table 1: Absolute difference |file1 - file2| for each lambda column.\n")
    out.write("#   Table 2: Relative error |file1 - file2| / max(|f1|,|f2|) as percentage.\n")
    out.write("#\n")
    out.write(f"# absolute tolerance: {tol:.0e}\n")
    out.write("# relative tolerance: 1%\n")
    out.write("#\n")


def write_block(out, label, col_names, col_units, rows1, rows2, tol, mode="abs"):
    if not rows1 or not rows2:
        out.write(f"# [SKIP] {label}: no data\n")
        return

    if len(rows1) != len(rows2):
        out.write(
            f"# [SKIP] {label}: row count mismatch ({len(rows1)} vs {len(rows2)})\n"
        )
        return

    ncols = min(
        min(len(r) for r in rows1),
        min(len(r) for r in rows2),
    )

    hdr_names = list(col_names)
    hdr_units = list(col_units)

    out.write(f"# --- {label} ---\n")
    out.write(f"# columns: {ncols}\n")
    out.write(f"# entries: {len(rows1)}\n")
    out.write("##" + "".join(f"{n:>14}" for n in hdr_names[:ncols]) + "\n")
    out.write("##" + "".join(f"{u:>14}" for u in hdr_units[:ncols]) + "\n")

    max_val = 0.0
    max_row = -1
    max_col = -1
    mismatch_count = 0
    total = 0

    for i, (r1, r2) in enumerate(zip(rows1, rows2)):
        out.write("".join(f"{r1[j]:>14.0f}" for j in range(4)))
        for j in range(4, ncols):
            if mode == "abs":
                val = abs(r1[j] - r2[j])
                unit_tol = tol
            else:
                denom = max(abs(r1[j]), abs(r2[j]), 1e-300)
                val = abs(r1[j] - r2[j]) / denom * 100.0
                unit_tol = tol
            total += 1
            if val > max_val:
                max_val = val
                max_row = i
                max_col = j
            if val > unit_tol:
                mismatch_count += 1
            if mode == "abs":
                out.write(f" {val:>14.6E}")
            else:
                out.write(f" {val:>13.4f}%")
        out.write("\n")

    out.write("#\n")
    if mismatch_count == 0:
        out.write(f"# PASS: {total} values match within tol={tol:.0e}\n")
    else:
        out.write(
            f"# FAIL: {mismatch_count}/{total} values differ beyond tol={tol:.0e}\n"
        )
    tag = "difference" if mode == "abs" else "relative error"
    out.write(
        f"# Max {tag}: {max_val:.6e}  "
        f"(row {max_row}, col {col_names[max_col]})\n"
    )


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    outpath = sys.argv[3]
    tol = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-4

    rows1 = parse_file(file1)
    rows2 = parse_file(file2)

    with open(outpath, "w") as out:
        write_header(out, file1, rows1, file2, rows2, tol)

        write_block(
            out,
            "internal transition rate (absolute)",
            COLUMN_NAMES,
            COLUMN_UNITS,
            rows1,
            rows2,
            tol,
            mode="abs",
        )

        out.write("\n")

        write_block(
            out,
            "internal transition rate (relative)",
            REL_COLUMN_NAMES,
            REL_COLUMN_UNITS,
            rows1,
            rows2,
            1.0,
            mode="rel",
        )

    print(f"Report written to: {outpath}")

    rel_tol = 1.0
    all_passed = True
    for label, r1, r2, tol_used, mode in [
        ("abs", rows1, rows2, tol, "abs"),
        ("rel", rows1, rows2, rel_tol, "rel"),
    ]:
        if not r1 or not r2:
            print(f"[SKIP] {label}: no data")
            continue
        if len(r1) != len(r2):
            print(f"[SKIP] {label}: row count mismatch")
            continue
        ncols = min(min(len(r) for r in r1), min(len(r) for r in r2))
        mismatches = 0
        total = 0
        max_val = 0.0
        for a, b in zip(r1, r2):
            for j in range(4, ncols):
                if mode == "abs":
                    val = abs(a[j] - b[j])
                else:
                    denom = max(abs(a[j]), abs(b[j]), 1e-300)
                    val = abs(a[j] - b[j]) / denom * 100.0
                total += 1
                max_val = max(max_val, val)
                if val > tol_used:
                    mismatches += 1
        tag = "diff" if mode == "abs" else "rel err"
        if mismatches == 0:
            print(f"PASS [{tag}]: {total} values match within tol={tol_used:.0e}")
        else:
            print(
                f"FAIL [{tag}]: {mismatches}/{total} values differ "
                f"beyond tol={tol_used:.0e}"
            )
            all_passed = False
        print(f"       max {tag} = {max_val:.6e}")

    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
