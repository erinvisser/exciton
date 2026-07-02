#!/usr/bin/env python3
"""
compare_diagnostics.py — Compare TALYS_DIAG vs Exciton # DIAG diagnostic lines.

Usage:
  python3 compare_diagnostics.py <talys_file> <exciton_file>
      [--states 0010 0021 1121 ...]
      [--output-dir <dir>]
      [--quiet]

Input formats:
  TALYS: lines starting with "TALYS_DIAG" (Fortran write to stdout)
  Exciton: lines starting with "# DIAG" (C++ std::cerr output)

Output:
  One formatted comparison table per selected state per term type,
  written to <output-dir>/compare_<timestamp>.log

Author: auto-generated for OMP validation
"""

import argparse
import os
import re
import sys
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Default states to compare: (p_pi, h_pi, p_nu, h_nu)
# ---------------------------------------------------------------------------
DEFAULT_STATES: List[Tuple[int, int, int, int]] = [
    (0, 0, 1, 0),  # 0010  initial
    (0, 0, 2, 1),  # 0021  pure-neutron 3-exciton
    (1, 1, 0, 0),  # 1100  pure-proton 2-exciton
    (1, 1, 2, 1),  # 1121  mixed 5-exciton
    (2, 2, 1, 0),  # 2210  mixed 5-exciton no neutron holes
    (3, 3, 2, 1),  # 3321  deep mixed 9-exciton
    (4, 4, 1, 0),  # 4410  deep mixed 9-exciton no neutron holes
    (5, 5, 1, 0),  # 5510  deepest mixed 11-exciton
]


def state_label(p_pi: int, h_pi: int, p_nu: int, h_nu: int) -> str:
    """Return a compact 4-digit state label."""
    return f"{p_pi}{h_pi}{p_nu}{h_nu}"


def state_tuple(s: str) -> Optional[Tuple[int, int, int, int]]:
    """Parse a 4-digit state string like '0021' -> (0,0,2,1)."""
    if len(s) != 4 or not s.isdigit():
        return None
    return (int(s[0]), int(s[1]), int(s[2]), int(s[3]))


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


def _extract_kv(text: str, key: str, conv=float) -> Any:
    """Extract value for *key=...* from text, converting with *conv*."""
    m = re.search(r"\b" + re.escape(key) + r"=\s*([-0-9eE.+]+)", text)
    if m:
        try:
            return conv(m.group(1))
        except Exception:
            return None
    return None


def _int(text: str, key: str) -> Optional[int]:
    return _extract_kv(text, key, int)


def _float(text: str, key: str) -> Optional[float]:
    return _extract_kv(text, key, float)


# ----- TALYS parser -------------------------------------------------------

# Matches: leading space, then TALYS_DIAG <function> state=<i4> <i4> <i4> <i4>
TALYS_RE = re.compile(
    r"^\s*TALYS_DIAG\s+(?P<func>\S+)"
    r"\s+state=\s*(?P<a>\d+)\s*(?P<b>\d+)\s*(?P<c>\d+)\s*(?P<d>\d+)"
)


def parse_talys_line(line: str) -> Optional[Dict[str, Any]]:
    m = TALYS_RE.search(line)
    if not m:
        return None
    d: Dict[str, Any] = {}
    d["_source"] = "TALYS"
    d["function"] = m.group("func")
    d["p_pi"] = int(m.group("a"))
    d["h_pi"] = int(m.group("b"))
    d["p_nu"] = int(m.group("c"))
    d["h_nu"] = int(m.group("d"))
    d["k"] = _int(line, "k")
    d["Zix"] = _int(line, "Zix")
    d["Nix"] = _int(line, "Nix")
    d["S"] = _float(line, "S")
    d["uu"] = _float(line, "uu")
    d["eopt"] = _float(line, "eopt")
    d["nen"] = _int(line, "nen")
    d["wvol"] = _float(line, "wvol")
    d["Weff"] = _float(line, "Weff")
    d["lamcol"] = _float(line, "lamcol")
    d["ratio"] = _float(line, "ratio")  # may be None for particle terms
    # j is present in lambdanuplus/lambdapiplus, not in lambdanupi/lambdapinu
    d["j"] = _int(line, "j")
    return d


# ----- Exciton parser -----------------------------------------------------

EXCITON_RE = re.compile(
    r"^#\s+DIAG\s+(?P<diag_type>pair_creation|conversion)\s+"
    r"(particle=(?P<particle>\w+)|type=(?P<conv_type>\w+))"
)


def parse_exciton_line(line: str) -> Optional[Dict[str, Any]]:
    m = EXCITON_RE.search(line)
    if not m:
        return None
    d: Dict[str, Any] = {}
    d["_source"] = "Exciton"
    d["diag_type"] = m.group("diag_type")
    if m.group("diag_type") == "pair_creation":
        d["particle"] = m.group("particle")
    else:
        d["conv_type"] = m.group("conv_type")
    d["p_pi"] = _int(line, "p_pi")
    d["h_pi"] = _int(line, "h_pi")
    d["p_nu"] = _int(line, "p_nu")
    d["h_nu"] = _int(line, "h_nu")
    d["k"] = _int(line, "k")
    d["j"] = _int(line, "j")
    d["sep"] = _float(line, "sep")
    d["uu"] = _float(line, "uu")
    d["eopt"] = _float(line, "eopt")
    d["nen"] = _int(line, "nen")
    d["wvol"] = _float(line, "wvol")
    d["Weff"] = _float(line, "Weff")
    d["lambda_col"] = _float(line, "lambda_col")
    # optional fields
    d["Z_res"] = _int(line, "Z_res")
    d["A_res"] = _int(line, "A_res")
    d["pdg"] = _int(line, "pdg")
    return d


# ---------------------------------------------------------------------------
# Matching
# ---------------------------------------------------------------------------


def _match_key_talys(entry: Dict[str, Any]) -> Tuple:
    """Build a matching key from a TALYS diagnostic entry."""
    f = entry["function"]
    st = (entry["p_pi"], entry["h_pi"], entry["p_nu"], entry["h_nu"])
    k = entry.get("k")
    j = entry.get("j")
    return (f, st, k, j)


def _match_key_exciton(entry: Dict[str, Any]) -> Tuple:
    """Build a matching key from an Exciton diagnostic entry."""
    dt = entry["diag_type"]
    particle = entry.get("particle")
    conv_type = entry.get("conv_type")
    st = (entry["p_pi"], entry["h_pi"], entry["p_nu"], entry["h_nu"])
    k = entry.get("k")
    j = entry.get("j")

    # Map Exciton particle/conv_type to the TALYS function name
    if dt == "pair_creation" and particle == "neutron":
        func = "lambdanuplus"
    elif dt == "pair_creation" and particle == "proton":
        func = "lambdapiplus"
    elif dt == "conversion" and conv_type == "NeutronToProton":
        func = "lambdanupi"
    elif dt == "conversion" and conv_type == "ProtonToNeutron":
        func = "lambdapinu"
    else:
        func = None
    return (func, st, k, j)


def match_entries(talys: List[Dict], exciton: List[Dict]) -> List[Tuple[Dict, Dict]]:
    """Return matched (talys_entry, exciton_entry) pairs."""
    # Index exciton entries by match key
    ex_index: Dict[Tuple, Dict] = {}
    for e in exciton:
        key = _match_key_exciton(e)
        if key[0] is not None:
            ex_index[key] = e

    pairs: List[Tuple[Dict, Dict]] = []
    for t in talys:
        key = _match_key_talys(t)
        e = ex_index.get(key)
        if e is not None:
            pairs.append((t, e))
    return pairs


# ---------------------------------------------------------------------------
# Table generation
# ---------------------------------------------------------------------------

_PARAM_ROWS = [
    # (label, talys_key, exciton_key, fmt)
    ("S(n) [MeV]", "S", "sep", ".6f"),
    ("uu [MeV]", "uu", "uu", ".6f"),
    ("e_kin [MeV]", "eopt", "eopt", ".6f"),
    ("nen", "nen", "nen", "d"),
    ("wvol [MeV]", "wvol", "wvol", ".6e"),
    ("Weff [MeV]", "Weff", "Weff", ".6e"),
    ("collision rate [1/s]", "lamcol", "lambda_col", ".6e"),
    ("ratio (hole term)", "ratio", None, ".6f"),
]


def fmt_val(val, fmt: str) -> str:
    if val is None:
        return "N/A"
    if fmt == "d":
        return f"{int(val):d}"
    # scientific for large numbers, fixed otherwise
    if "e" in fmt and abs(val) < 1e10:
        return f"{val:{fmt}}"
    elif "e" in fmt:
        return f"{val:{fmt}}"
    return f"{val:{fmt}}"


def generate_table(
    talys_entry: Dict, exciton_entry: Dict, state_str: str, func_label: str
) -> str:
    """Return a formatted comparison table for one matched pair."""
    lines: List[str] = []
    sep = "-" * 78
    lines.append("# " + sep)
    lines.append(
        f"# State {state_str}: p_pi={talys_entry['p_pi']}  "
        f"h_pi={talys_entry['h_pi']}  "
        f"p_nu={talys_entry['p_nu']}  "
        f"h_nu={talys_entry['h_nu']}"
    )
    lines.append(f"# {func_label}")
    lines.append("# " + sep)

    # Column widths
    col_source = 12
    col_param = 24
    col_val = 18
    header = (
        f"{'Source':>{col_source}}  {'Parameter':<{col_param}}  "
        f"{'TALYS':>{col_val}}  {'Exciton':>{col_val}}  "
        f"{'|Diff|':>{col_val}}  {'Rel%':>10}"
    )
    lines.append("#")
    lines.append(header)
    lines.append("#  " + "-" * (col_source + col_param + 3 * col_val + 12))

    max_rel = 0.0
    for label, t_key, e_key, fmt in _PARAM_ROWS:
        t_val = talys_entry.get(t_key)
        e_val = exciton_entry.get(e_key) if e_key else None
        if t_val is None and e_val is None:
            continue

        t_str = fmt_val(t_val, fmt)
        e_str = fmt_val(e_val, fmt)

        # Absolute difference
        if t_val is not None and e_val is not None:
            diff = abs(t_val - e_val)
            diff_str = f"{diff:.6e}"
            # Relative error
            denom = max(abs(t_val), abs(e_val))
            if denom > 0:
                rel = diff / denom * 100.0
            else:
                rel = 0.0
            if rel > max_rel:
                max_rel = rel
            rel_str = f"{rel:.4f}%"
        else:
            diff_str = "N/A"
            rel_str = "N/A"

        line = (
            f"{'':>{col_source}}  {label:<{col_param}}  "
            f"{t_str:>{col_val}}  {e_str:>{col_val}}  "
            f"{diff_str:>{col_val}}  {rel_str:>10}"
        )
        lines.append(line)

    lines.append("#  " + "-" * (col_source + col_param + 3 * col_val + 12))
    lines.append(f"#  State {state_str} max rel. error: {max_rel:.4f}%")
    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare TALYS vs Exciton lambda diagnostic output."
    )
    parser.add_argument("talys_file", help="TALYS diagnostic file (TALYS_DIAG lines)")
    parser.add_argument("exciton_file", help="Exciton diagnostic file (# DIAG lines)")
    parser.add_argument(
        "--states",
        nargs="+",
        default=None,
        help="State labels to compare, e.g. 0010 0021 1121",
    )
    parser.add_argument(
        "--output-dir", default=None, help="Output directory (default: macOS path)"
    )
    parser.add_argument("--quiet", action="store_true", help="Suppress progress output")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # Determine output directory
    macos_dir = "/Users/erinvisser/Downloads/opencode/exciton/lambda_printouts/OMP/TALYS_Mimic/Comparison"
    out_dir = args.output_dir or macos_dir

    # Ensure output directory exists
    os.makedirs(out_dir, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = os.path.join(out_dir, f"compare_{timestamp}.log")

    # Determine selected states
    if args.states:
        selected = []
        for s in args.states:
            t = state_tuple(s)
            if t is None:
                print(
                    f'error: invalid state label "{s}" — use 4 digits', file=sys.stderr
                )
                sys.exit(1)
            selected.append(t)
    else:
        selected = DEFAULT_STATES
    selected_set = set(selected)

    # ---- Parse files ----
    talys_entries: List[Dict] = []
    with open(args.talys_file) as f:
        for line in f:
            d = parse_talys_line(line)
            if d is not None:
                talys_entries.append(d)

    exciton_entries: List[Dict] = []
    with open(args.exciton_file) as f:
        for line in f:
            d = parse_exciton_line(line)
            if d is not None:
                exciton_entries.append(d)

    if not args.quiet:
        print(f"TALYS diagnostic entries:   {len(talys_entries)}", file=sys.stderr)
        print(f"Exciton diagnostic entries: {len(exciton_entries)}", file=sys.stderr)

    if not talys_entries:
        print("error: no TALYS_DIAG lines found", file=sys.stderr)
        sys.exit(1)
    if not exciton_entries:
        print("error: no Exciton # DIAG lines found", file=sys.stderr)
        sys.exit(1)

    # ---- Match ----
    pairs = match_entries(talys_entries, exciton_entries)

    if not args.quiet:
        print(f"Matched entries:            {len(pairs)}", file=sys.stderr)

    if not pairs:
        print("error: no matching state entries found", file=sys.stderr)
        print(
            "  Check that the diagnostic files use compatible formats.", file=sys.stderr
        )
        sys.exit(1)

    # ---- Filter to selected states ----
    filtered_pairs = []
    for t, e in pairs:
        st = (t["p_pi"], t["h_pi"], t["p_nu"], t["h_nu"])
        if st in selected_set:
            filtered_pairs.append((t, e))

    if not args.quiet:
        print(f"Entries for selected states: {len(filtered_pairs)}", file=sys.stderr)

    # ---- Build output ----
    out_lines: List[str] = []
    out_lines.append(
        "# ======================================================================"
    )
    out_lines.append("# Comparison: TALYS vs Exciton lambda diagnostics")
    out_lines.append(f"# Generated: {datetime.now().isoformat()}")
    out_lines.append(f"# TALYS file:  {args.talys_file}")
    out_lines.append(f"# Exciton file: {args.exciton_file}")
    out_lines.append(
        f"# Selected states: {' '.join(state_label(*s) for s in selected)}"
    )
    out_lines.append(f"# Matched pairs (all): {len(pairs)}")
    out_lines.append(f"# Matched pairs (selected): {len(filtered_pairs)}")
    out_lines.append(
        "# ======================================================================"
    )
    out_lines.append("")
    out_lines.append("# ---- Parameter descriptions ----")
    out_lines.append(f"# {'Parameter':<26s}  Meaning")
    out_lines.append("# " + "-" * 41)
    out_lines.append("#  S(n) [MeV]                  Neutron separation energy")
    out_lines.append("#  uu [MeV]                    Effective excitation energy available")
    out_lines.append("#  e_kin [MeV]                 Asymptotic kinetic energy fed to OMP")
    out_lines.append("#  nen                         Radial quadrature index")
    out_lines.append("#  wvol [MeV]                  Volume-integrated imaginary OMP strength")
    out_lines.append("#  Weff [MeV]                  Effective squared matrix element |M|^2")
    out_lines.append("#  collision rate [1/s]        Final transition rate lambda")
    out_lines.append("#  ratio (hole term)           Density-of-states ratio for hole terms")
    out_lines.append("# " + "-" * 41)
    out_lines.append("#  |Diff| = absolute difference between TALYS and Exciton")
    out_lines.append("#  Rel%   = |Diff| / max(|TALYS|,|Exciton|) * 100")
    out_lines.append(
        "# ======================================================================"
    )
    out_lines.append("")

    for idx, (t, e) in enumerate(filtered_pairs):
        st = (t["p_pi"], t["h_pi"], t["p_nu"], t["h_nu"])
        sl = state_label(*st)

        # Build a human-readable function label
        func = t["function"]
        j = t.get("j")
        k = t.get("k")
        label_parts = [func]
        if j is not None:
            label_parts.append(f"j={j}")
        if k is not None:
            label_parts.append(f"k={k}")

        # Map to Exciton side description
        dt = e["diag_type"]
        if dt == "pair_creation":
            label_parts.append(f"({e['particle']} pair creation)")
        else:
            label_parts.append(f"({e['conv_type']} conversion)")

        table = generate_table(t, e, sl, "  ".join(label_parts))
        out_lines.append(table)

        if idx < len(filtered_pairs) - 1:
            out_lines.append("")

    out_lines.append(
        "# ======================================================================"
    )
    out_lines.append("# End of comparison")
    out_lines.append(
        "# ======================================================================"
    )

    # ---- Write ----
    with open(out_path, "w") as f:
        f.write("\n".join(out_lines) + "\n")

    if not args.quiet:
        print(f"Output: {out_path}", file=sys.stderr)

    sys.exit(0)


if __name__ == "__main__":
    main()
