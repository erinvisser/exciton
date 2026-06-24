#!/usr/bin/env python3
"""Run TALYS, evaluate MARLEY at TALYS's radii, and compare U(r).

Usage:
  python3 talys_vs_marley.py --energies 20.0 --talys-exe /path/to/talys \
    --marley-exe /path/to/exciton --l 0 -o comparison.dat

Output columns: E(MeV)  r(fm)  Re(MARLEY)(MeV)  Im(MARLEY)(MeV)  Re(TALYS)(MeV)  Im(TALYS)(MeV)  rel_diff_Re  rel_diff_Im

Re(TALYS) = -V - Vso_raw * ls
Im(TALYS) = -W - Wso_raw * ls

Relative difference = (MARLEY - TALYS) / |TALYS|, or 0 if TALYS == 0.
"""

import argparse
import os
import subprocess
import sys


def parse_energies(s):
    if ":" in s:
        parts = s.split(":")
        if len(parts) == 3:
            start, stop, step = float(parts[0]), float(parts[1]), float(parts[2])
            if step <= 0:
                raise argparse.ArgumentTypeError("step must be positive")
            energies = []
            e = start
            while e <= stop + 1e-10:
                energies.append(e)
                e += step
            return energies
        raise argparse.ArgumentTypeError("expected start:stop:step")
    energies = []
    for token in s.split(","):
        token = token.strip()
        if not token:
            continue
        energies.append(float(token))
    return energies


def main():
    parser = argparse.ArgumentParser(
        description="Compare MARLEY and TALYS total OMP U(r)"
    )
    parser.add_argument(
        "--energies", required=True, help="Comma-separated list or start:stop:step"
    )
    parser.add_argument(
        "--l",
        type=int,
        default=0,
        dest="l_val",
        help="l quantum number for SO eigenvalue (default: 0)",
    )
    parser.add_argument(
        "--output", "-o", default=None, help="Output file (default: stdout)"
    )
    parser.add_argument("--element", default="ar")
    parser.add_argument("--mass", type=int, default=40)
    parser.add_argument("--projectile", default="n")
    parser.add_argument(
        "--talys-dir",
        default=None,
        help="TALYS working dir (default: script directory)",
    )
    parser.add_argument("--talys-exe", default="talys")
    parser.add_argument(
        "--marley-exe",
        default="./exciton",
        help="MARLEY exciton executable (default: ./exciton)",
    )
    parser.add_argument(
        "--marley-workdir",
        default=None,
        help="MARLEY working dir (default: same as --talys-dir)",
    )

    args = parser.parse_args()
    energies = parse_energies(args.energies)

    work_dir = args.talys_dir or os.path.dirname(os.path.abspath(__file__))
    if args.marley_workdir:
        marley_workdir = args.marley_workdir
    else:
        # Infer from executable path: go up 2 dirs (executables/lambda -> project root)
        marley_workdir = os.path.dirname(
            os.path.dirname(os.path.abspath(args.marley_exe))
        )

    out_fh = open(args.output, "w") if args.output else sys.stdout

    out_fh.write("# Comparison: MARLEY vs TALYS U(r)\n")
    out_fh.write(
        f"# System: {args.projectile} + {args.element}{args.mass}  (KD03 global)\n"
    )
    out_fh.write(f"# l = {args.l_val}  (l.s eigenvalue = {args.l_val})\n")
    out_fh.write(
        "# Columns:  E(MeV)  r(fm)  "
        "Re(MARLEY)(MeV)  Im(MARLEY)(MeV)  "
        "Re(TALYS)(MeV)  Im(TALYS)(MeV)  "
        "rel_diff_Re  rel_diff_Im\n"
    )

    eval_file = os.path.join(work_dir, "_marley_eval_list.inp")
    total_points = 0

    for energy in energies:
        inp_path = os.path.join(work_dir, "talys.inp")
        omp_path = os.path.join(work_dir, "omp.out")

        # --- Step 1: Run TALYS ---
        if os.path.exists(omp_path):
            os.remove(omp_path)

        lines = [
            f"energy {energy}",
            f"element {args.element}",
            f"mass {args.mass}",
            f"projectile {args.projectile}",
            "localomp n",
            "outomp y",
            "blockomp y",
            "",
        ]
        with open(inp_path, "w") as f:
            f.write("\n".join(lines))

        try:
            with open(inp_path) as inp_f:
                result = subprocess.run(
                    [args.talys_exe],
                    cwd=work_dir,
                    stdin=inp_f,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
        except FileNotFoundError:
            print(
                f'Error: TALYS executable "{args.talys_exe}" not found.',
                file=sys.stderr,
            )
            sys.exit(1)
        if result.returncode != 0:
            print(f"Error: TALYS failed at E={energy}", file=sys.stderr)
            sys.exit(1)
        if not os.path.exists(omp_path):
            print(f"Error: omp.out not found at E={energy}", file=sys.stderr)
            sys.exit(1)

        # --- Step 2: Parse omp.out to get radii and TALYS U ---
        talys_data = []  # list of (r, k, U_re, U_im)
        with open(omp_path) as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                parts = s.split()
                if len(parts) != 5:
                    continue
                try:
                    r = float(parts[0])
                    V = float(parts[1])
                    W = float(parts[2])
                    Vso = float(parts[3])
                    Wso = float(parts[4])
                except ValueError:
                    continue
                U_re_talys = -V - Vso * args.l_val
                U_im_talys = -W - Wso * args.l_val
                talys_data.append((r, U_re_talys, U_im_talys))

        # --- Step 3: Write MARLEY eval list (E r only) ---
        with open(eval_file, "w") as f:
            for r, _, _ in talys_data:
                f.write(f"{energy} {r}\n")

        # --- Step 4: Run MARLEY ---
        try:
            result = subprocess.run(
                [
                    args.marley_exe,
                    "--mode",
                    "eval_list",
                    "--input",
                    eval_file,
                    "--l",
                    str(args.l_val),
                    "--proj",
                    args.projectile,
                ],
                cwd=marley_workdir,
                capture_output=True,
                text=True,
            )
        except FileNotFoundError:
            print(
                f'Error: MARLEY executable "{args.marley_exe}" not found.',
                file=sys.stderr,
            )
            sys.exit(1)
        if result.returncode != 0:
            print(f"Error: MARLEY failed at E={energy}", file=sys.stderr)
            print(result.stderr, file=sys.stderr)
            sys.exit(1)

        # --- Step 5: Parse MARLEY output and combine ---
        marley_points = []  # list of (r, U_re, U_im)
        for line in result.stdout.splitlines():
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) != 4:
                continue
            try:
                E_m = float(parts[0])
                r_m = float(parts[1])
                re_m = float(parts[2])
                im_m = float(parts[3])
            except ValueError:
                continue
            marley_points.append((r_m, re_m, im_m))

        if len(marley_points) != len(talys_data):
            print(
                f"Warning: E={energy}: MARLEY gave {len(marley_points)} points, "
                f"TALYS gave {len(talys_data)}",
                file=sys.stderr,
            )

        for (r_m, re_m, im_m), (r_t, re_t, im_t) in zip(marley_points, talys_data):
            rel_re = (re_m - re_t) / abs(re_t) if abs(re_t) > 1e-30 else 0.0
            rel_im = (im_m - im_t) / abs(im_t) if abs(im_t) > 1e-30 else 0.0
            out_fh.write(
                f"  {energy:9.5f}  {r_m:8.5f}  "
                f"{re_m:16.8e}  {im_m:16.8e}  "
                f"{re_t:16.8e}  {im_t:16.8e}  "
                f"{rel_re:16.8e}  {rel_im:16.8e}\n"
            )
            total_points += 1

        n_r = len(talys_data)
        print(
            f"  E = {energy:6.2f} MeV: {n_r} radial points for projectile {args.projectile}",
            file=sys.stderr,
        )

    # Clean up
    if os.path.exists(eval_file):
        os.remove(eval_file)

    if args.output:
        out_fh.close()

    print(f"Wrote {total_points} lines to {args.output or 'stdout'}", file=sys.stderr)


if __name__ == "__main__":
    main()
