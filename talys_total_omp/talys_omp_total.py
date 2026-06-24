#!/usr/bin/env python3
"""Run TALYS at specified energies and extract total U(r) from omp.out.

Usage:
  python3 talys_omp_total.py --energies 20.0 --l 0 -o talys_u.dat
  python3 talys_omp_total.py --energies 0.0,5.0,10.0,20.0 --l 1
  python3 talys_omp_total.py --energies 0.0:20.0:2.0 --l 0

Output columns: E(MeV)  r(fm)  Re(U)(MeV)  Im(U)(MeV)

Re(U) = -V - Vso_raw * l
Im(U) = -W - Wso_raw * l

where V/W/Vso/Wso are the four potential columns from TALYS omp.out,
and l is the orbital angular momentum (l.s eigenvalue = l for aligned spin j=l+1/2).
For l=0 (s-wave), the spin-orbit term is zero.
"""

import argparse
import os
import subprocess
import sys


def parse_energies(s):
    """Parse --energies argument. Accepts comma-separated list or start:stop:step."""
    if ":" in s:
        parts = s.split(":")
        if len(parts) == 3:
            try:
                start = float(parts[0])
                stop = float(parts[1])
                step = float(parts[2])
            except ValueError:
                raise argparse.ArgumentTypeError(
                    f"invalid range '{s}': expected start:stop:step"
                )
            if step <= 0:
                raise argparse.ArgumentTypeError(
                    f"invalid range '{s}': step must be positive"
                )
            energies = []
            e = start
            while e <= stop + 1e-10:
                energies.append(e)
                e += step
            return energies
        raise argparse.ArgumentTypeError(
            f"invalid range '{s}': expected start:stop:step"
        )
    energies = []
    for token in s.split(","):
        token = token.strip()
        if not token:
            continue
        try:
            energies.append(float(token))
        except ValueError:
            raise argparse.ArgumentTypeError(f"invalid energy '{token}'")
    return energies


def main():
    parser = argparse.ArgumentParser(
        description="Run TALYS and extract total OMP U(r) at multiple energies"
    )
    parser.add_argument(
        "--energies", required=True, help="Comma-separated list or start:stop:step"
    )
    parser.add_argument(
        "--l",
        type=int,
        default=0,
        help="Orbital angular momentum for SO eigenvalue l.s=l (default: 0)",
    )
    parser.add_argument(
        "--output", "-o", default=None, help="Output file (default: stdout)"
    )
    parser.add_argument("--element", default="ar", help="Target element (default: ar)")
    parser.add_argument(
        "--mass", type=int, default=40, help="Target mass number (default: 40)"
    )
    parser.add_argument("--projectile", default="n", help="Projectile (default: n)")
    parser.add_argument(
        "--talys-dir",
        default=None,
        help="TALYS working directory (default: directory of this script)",
    )
    parser.add_argument(
        "--talys-exe", default="talys", help="TALYS executable (default: talys)"
    )
    parser.add_argument(
        "--blockomp",
        default="y",
        choices=["y", "n"],
        help="Combine all energies in one omp.out file (default: y)",
    )

    args = parser.parse_args()

    energies = parse_energies(args.energies)
    if not energies:
        print("Error: no energies specified", file=sys.stderr)
        sys.exit(1)

    work_dir = (
        args.talys_dir if args.talys_dir else os.path.dirname(os.path.abspath(__file__))
    )
    talys_exe = args.talys_exe

    out_fh = open(args.output, "w") if args.output else sys.stdout

    out_fh.write(f"# TALYS total OMP: {args.projectile} + {args.element}{args.mass}\n")
    out_fh.write("# KD03 global (localomp n)\n")
    out_fh.write(f"# l = {args.l}  (l.s eigenvalue = {args.l})\n")
    out_fh.write("# Columns:  E(MeV)  r(fm)  Re(U)(MeV)  Im(U)(MeV)\n")

    total_points = 0

    for energy in energies:
        inp_path = os.path.join(work_dir, "talys.inp")
        omp_path = os.path.join(work_dir, "omp.out")

        # Remove old omp.out before each run
        if os.path.exists(omp_path):
            os.remove(omp_path)

        lines = [
            f"energy {energy}",
            f"element {args.element}",
            f"mass {args.mass}",
            f"projectile {args.projectile}",
            "localomp n",
            "outomp y",
            f"blockomp {args.blockomp}",
            "",
        ]
        with open(inp_path, "w") as f:
            f.write("\n".join(lines))

        try:
            with open(inp_path) as inp_f:
                result = subprocess.run(
                    [talys_exe],
                    cwd=work_dir,
                    stdin=inp_f,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
        except FileNotFoundError:
            print(f'Error: TALYS executable "{talys_exe}" not found.', file=sys.stderr)
            print("Use --talys-exe to specify the full path, e.g.:", file=sys.stderr)
            print(
                f"  python3 {sys.argv[0]} --energies ... --talys-exe /path/to/talys",
                file=sys.stderr,
            )
            sys.exit(1)
        if result.returncode != 0:
            print(
                f"Error: TALYS failed at energy {energy} (exit code {result.returncode})",
                file=sys.stderr,
            )
            sys.exit(1)

        if not os.path.exists(omp_path):
            print(
                f"Error: omp.out not found after TALYS run at E={energy}",
                file=sys.stderr,
            )
            sys.exit(1)

        n_parsed = 0
        with open(omp_path) as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                parts = stripped.split()
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
                U_re = -V - Vso * args.l
                U_im = -W - Wso * args.l
                out_fh.write(f"  {energy:9.5f}  {r:8.5f}  {U_re:16.8e}  {U_im:16.8e}\n")
                n_parsed += 1
        total_points += n_parsed
        print(
            f"  E = {energy:6.2f} MeV: parsed {n_parsed} radial points", file=sys.stderr
        )

    if args.output:
        out_fh.flush()
        out_fh.close()

    out_name = args.output or "stdout"
    print(f"Wrote {total_points} lines to {out_name}", file=sys.stderr)


if __name__ == "__main__":
    main()
