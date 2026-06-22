#!/usr/bin/env python3
"""
Compare two ph_density.out files containing two data blocks:
  Block 1: particle-hole state density  (25 rows, Ex=1..25)
  Block 2: spin distribution            (12 rows, n=1..12)

Rows are identified by YAML section headers (# quantity: ...).
Output shows two tables per block: absolute diff and relative error.

Usage:
    python3 compare_omega.py <file1> <file2> <output> [tolerance]

Default tolerance: 1e-4 (absolute) / 1% (relative)
"""

import sys
from datetime import datetime

STATE_COLUMNS = [
    "Ex", "P(3)", "g(p)", "g(n)",
    "1p1h0p0h", "0p0h1p1h", "1p1h1p0h", "1p0h1p1h",
    "2p1h0p0h", "0p0h2p1h", "2p2h0p0h", "0p0h2p2h", "1p1h1p1h",
]

SPIN_COLUMNS = ["n"] + [f"J={j}" for j in range(9)] + ["Sum"]

STATE_UNITS = ["MeV", "MeV", "MeV^-1", "MeV^-1"] + ["MeV^-1"] * 9

SPIN_UNITS = [""] * 11


def parse_file(path):
    with open(path) as f:
        lines = f.readlines()

    state_rows = []
    spin_rows = []
    current = None

    for line in lines:
        raw = line.strip()
        if not raw or raw.startswith("##"):
            continue

        if raw.startswith("#"):
            if "type:" in raw:
                if "particle-hole" in raw:
                    current = "state"
                elif "spin" in raw:
                    current = "spin"
            continue

        parts = raw.split()
        row = []
        for x in parts:
            try:
                row.append(float(x))
            except ValueError:
                row = None
                break
        if row is None:
            continue

        if current == "state":
            state_rows.append(row)
        elif current == "spin":
            spin_rows.append(row)

    return state_rows, spin_rows


def write_header(out, file1, state_rows, spin_rows, file2, state_rows2, spin_rows2, tol):
    out.write("# omega comparison report\n")
    out.write(f"# generated: {datetime.now().isoformat()}\n")
    out.write("#\n")
    out.write("# file1: {}\n".format(file1))
    out.write("#   particle-hole state density: {} rows, {} cols\n".format(
        len(state_rows), len(state_rows[0]) if state_rows else 0))
    out.write("#   spin distribution:           {} rows, {} cols\n".format(
        len(spin_rows), len(spin_rows[0]) if spin_rows else 0))
    out.write("#\n")
    out.write("# file2: {}\n".format(file2))
    out.write("#   particle-hole state density: {} rows, {} cols\n".format(
        len(state_rows2), len(state_rows2[0]) if state_rows2 else 0))
    out.write("#   spin distribution:           {} rows, {} cols\n".format(
        len(spin_rows2), len(spin_rows2[0]) if spin_rows2 else 0))
    out.write("#\n")
    out.write("# Description:\n")
    out.write("#   Table 1: particle-hole state density — absolute difference\n")
    out.write("#   Table 2: particle-hole state density — relative error (%)\n")
    out.write("#   Table 3: spin distribution           — absolute difference\n")
    out.write("#   Table 4: spin distribution           — relative error (%)\n")
    out.write("#\n")
    out.write(f"# absolute tolerance: {tol:.0e}\n")
    out.write("# relative tolerance: 1%\n")
    out.write("# values shown: abs(file1 - file2) | rel = |f1-f2|/max(|f1|,|f2|)\n")
    out.write("#\n")


def write_block(out, label, col_names, col_units, rows1, rows2, tol, mode="abs"):
    if not rows1 or not rows2:
        out.write(f"# [SKIP] {label}: no data\n")
        return

    if len(rows1) != len(rows2):
        out.write(f"# [SKIP] {label}: row count mismatch "
                  f"({len(rows1)} vs {len(rows2)})\n")
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
        if mode == "abs":
            out.write(f"  {r1[0]:>14.6E}")
        else:
            out.write(f"  {r1[0]:>14.6E}")
        for j in range(1, ncols):
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
            f"# FAIL: {mismatch_count}/{total} values differ beyond "
            f"tol={tol:.0e}\n"
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

    state_rows, spin_rows = parse_file(file1)
    state_rows2, spin_rows2 = parse_file(file2)

    with open(outpath, "w") as out:
        write_header(out, file1, state_rows, spin_rows,
                     file2, state_rows2, spin_rows2, tol)

        write_block(out, "particle-hole state density (absolute)",
                    STATE_COLUMNS, STATE_UNITS,
                    state_rows, state_rows2, tol, mode="abs")

        out.write("\n")

        write_block(out, "particle-hole state density (relative)",
                    STATE_COLUMNS, ["MeV"] + ["%"] * 12,
                    state_rows, state_rows2, 1.0, mode="rel")

        out.write("\n")

        write_block(out, "spin distribution (absolute)",
                    SPIN_COLUMNS, SPIN_UNITS,
                    spin_rows, spin_rows2, tol, mode="abs")

        out.write("\n")

        write_block(out, "spin distribution (relative)",
                    SPIN_COLUMNS, [""] + ["%"] * 10,
                    spin_rows, spin_rows2, 1.0, mode="rel")

    print(f"Report written to: {outpath}")

    # Overall summary to stdout
    rel_tol = 1.0
    all_passed = True
    for label, r1, r2, tol_used, mode in [
        ("state density abs", state_rows, state_rows2, tol, "abs"),
        ("state density rel", state_rows, state_rows2, rel_tol, "rel"),
        ("spin distribution abs", spin_rows, spin_rows2, tol, "abs"),
        ("spin distribution rel", spin_rows, spin_rows2, rel_tol, "rel"),
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
            for j in range(1, ncols):
                if mode == "abs":
                    val = abs(a[j] - b[j])
                else:
                    denom = max(abs(a[j]), abs(b[j]), 1e-300)
                    val = abs(a[j] - b[j]) / denom * 100.0
                total += 1
                max_val = max(max_val, val)
                if val > tol_used:
                    mismatches += 1
        if mismatches == 0:
            print(f"PASS [{label}]: {total} values match within tol={tol_used:.0e}")
        else:
            print(
                f"FAIL [{label}]: {mismatches}/{total} values differ "
                f"beyond tol={tol_used:.0e}"
            )
            all_passed = False
        print(f"       max {'diff' if mode == 'abs' else 'rel err'} = {max_val:.6e}")

    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
