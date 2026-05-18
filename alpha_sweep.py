"""
Alpha sweep of an OpenVSP geometry through AVL.

For each angle of attack in ALPHAS, runs AVL and records CL and CDi
(Trefftz-plane induced drag). Writes one CSV per configuration:
    alpha_sweep_<config>.csv  with columns: alpha, CL, CDi, e

Usage: edit the CONFIG list, then `python alpha_sweep.py`.
"""

import csv
import os
import re
import subprocess
import sys

# ---- Config -----------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))

def _find(name, env_var, hints):
    """Locate a file by walking up from HERE and checking hint subdirs. Lets
    you drop this folder anywhere near the 'CDR - Aero stuff' tree without
    editing paths. Override with the env var if auto-discovery picks wrong."""
    if env_var and os.environ.get(env_var):
        p = os.environ[env_var]
        if os.path.isfile(p):
            return p
    cur = HERE
    for _ in range(4):
        for hint in hints:
            p = os.path.join(cur, hint, name)
            if os.path.isfile(p):
                return os.path.abspath(p)
        cur = os.path.dirname(cur)
    return None

_HINTS = [
    ".",
    "CDR - Aero stuff",
    os.path.join("CDR - Aero stuff-20260517T062817Z-3-001", "CDR - Aero stuff"),
    os.path.join("..", "CDR - Aero stuff"),
    os.path.join("..", "CDR - Aero stuff-20260517T062817Z-3-001", "CDR - Aero stuff"),
]
AVL_BIN    = _find("avl352.exe",                    "AVL_BIN",    _HINTS)
COLLAB_AVL = _find("ModWing_CorrSpan_LargeSlat.avl","COLLAB_AVL", _HINTS)

# Alpha range: -5 to 20 in 1-deg steps. Wide enough to capture both flapped and
# clean configs without wasting AVL calls.
ALPHAS = list(range(-5, 21, 1))

# Each config gets one alpha sweep that uses the collaborator's existing AVL
# geometry (ModWing_CorrSpan_LargeSlat.avl). That .avl already has the wing,
# tails, and named CONTROL surfaces (VertStabCS=d1, InnerCS=d2, Slats=d3,
# OuterCS=d4), so flap-deflected polars come from real vortex-lattice runs.

CONFIGS = [
    # Clean cruise: full Mach, no flap deflection.
    {"key": "cruise",  "avl": COLLAB_AVL, "mach": 0.85, "deflections": {}},
    # Takeoff: 30 deg of trailing-edge flap.
    {"key": "takeoff", "avl": COLLAB_AVL, "mach": 0.20, "deflections": {"d2": 30.0}},
    # Landing: 65 deg trailing-edge flap.
    {"key": "landing", "avl": COLLAB_AVL, "mach": 0.18, "deflections": {"d2": 65.0}},
]

# ---- AVL runner -------------------------------------------------------------
_MACH_HEADER_RE = re.compile(r"(#\s*Mach.*\n)\s*[-+0-9.Ee]+", re.IGNORECASE)

# In the collaborator's .avl, InnerCS (the trailing-edge flap) is only declared
# at the outer sections of its span, leaving the intermediate sections Y=7.02
# and Y=-7.02 without it. AVL spans a control between successive sections that
# both declare the same control, so the gap broke the flap into two tiny
# panels. Patch those intermediate sections by injecting the InnerCS CONTROL
# block (just under their AFILE line) before handing to AVL.
_INNERCS_HINGE = "InnerCS   1.0   0.759467   0. 0. 0.   1   | injected for continuity"

def _patch_innercs_continuity(text):
    # Section signature: "<Xle> <Yle> ..." line whose Yle is +/-7.017742, then
    # AFILE naca64206.dat, then (no InnerCS) next section. Insert CONTROL block
    # right after the AFILE line.
    out = []
    lines = text.splitlines(keepends=True)
    i = 0
    target_yle = ("7.017742", "-7.017742")
    while i < len(lines):
        out.append(lines[i])
        stripped = lines[i].strip().split()
        # Match a section data row like "17.567299   7.017742   0.000000  ..."
        if (len(stripped) >= 7
                and any(t in lines[i] for t in target_yle)
                and "Yle" not in lines[i]
                and "#" not in lines[i]):
            # Find the AFILE block, then check if next CONTROL is InnerCS
            j = i + 1
            inserted = False
            while j < len(lines) and not lines[j].lstrip().startswith("SECTION"):
                out.append(lines[j])
                if (lines[j].strip().lower().startswith("afile")
                        and j + 1 < len(lines)):
                    out.append(lines[j + 1])
                    j += 2
                    # Walk forward over any blank/comment lines until next non-blank
                    while j < len(lines) and lines[j].strip().startswith("CONTROL"):
                        # Already has a CONTROL — insert InnerCS ahead of it
                        out.append("CONTROL\n")
                        out.append(_INNERCS_HINGE + "\n")
                        inserted = True
                        out.append(lines[j])
                        out.append(lines[j + 1])
                        j += 2
                    if not inserted:
                        out.append("CONTROL\n")
                        out.append(_INNERCS_HINGE + "\n")
                        inserted = True
                    continue
                j += 1
            i = j
            continue
        i += 1
    return "".join(out)

def run_alpha_sweep(avl_path, alphas, mach, deflections, work_dir):
    """avl_path: absolute path to an existing .avl file. Returns list of forces
    file paths. `deflections` is a dict like {"d2": 30.0} for control deflection
    in degrees, where d1/d2/.. number matches the CONTROL order in the .avl."""
    with open(avl_path) as f:
        text = f.read()
    # AVL takes Mach from the .avl header. Rewriting the header is more robust
    # than driving AVL's interactive parameter editor.
    text = _MACH_HEADER_RE.sub(rf"\g<1>{mach}", text, count=1)
    if deflections:
        text = _patch_innercs_continuity(text)
    geom_file = os.path.join(work_dir, "geom.avl")
    with open(geom_file, "w") as f:
        f.write(text)

    cmd_lines = ["plop", "g", "", "load geom.avl", "oper"]
    out_files = []
    for a in alphas:
        out = f"forces_alpha_{a:+d}.txt"
        out_files.append(out)
        cmd_lines += [f"a a {a}"]
        for ctrl, val in deflections.items():
            cmd_lines += [f"{ctrl} {ctrl} {val}"]
        cmd_lines += ["x", "ft", out]
    cmd_lines += ["", "quit"]
    cmd = "\n".join(cmd_lines) + "\n"

    proc = subprocess.run(
        [AVL_BIN],
        input=cmd,
        capture_output=True,
        text=True,
        cwd=work_dir,
        timeout=600,
    )
    log_path = os.path.join(work_dir, "avl_run.log")
    with open(log_path, "w") as f:
        f.write("=== AVL COMMAND SCRIPT ===\n")
        f.write(cmd)
        f.write("\n=== AVL STDOUT ===\n")
        f.write(proc.stdout)
        f.write("\n=== AVL STDERR ===\n")
        f.write(proc.stderr)
    # AVL sometimes exits non-zero from EOF on stdin even after writing all
    # forces files. Only fail if the forces files are missing.
    missing = [f for f in out_files
               if not os.path.isfile(os.path.join(work_dir, f))]
    if missing:
        print("AVL stdout tail:", proc.stdout[-500:])
        raise RuntimeError(
            f"AVL run did not produce {len(missing)} forces files "
            f"(first missing: {missing[0]}); see {log_path}"
        )
    return [os.path.join(work_dir, f) for f in out_files]


# ---- Force-file parser ------------------------------------------------------
# AVL's "ft" output prints CLtot, CDtot, CDind, CLff, CDff (Trefftz), e on
# specific lines. We pull CLff / CDff (the inviscid Trefftz-plane values) since
# those are what go into the inviscid portion of the drag polar.
_RE_ALPHA = re.compile(r"Alpha\s*=\s*([-+0-9.Ee]+)")
_RE_CLFF = re.compile(r"CLff\s*=\s*([-+0-9.Ee]+)\s+CDff\s*=\s*([-+0-9.Ee]+)")
_RE_E = re.compile(r"e\s*=\s*([-+0-9.Ee]+)")

def parse_forces(path):
    with open(path) as f:
        text = f.read()
    a = _RE_ALPHA.search(text)
    clcd = _RE_CLFF.search(text)
    e = _RE_E.search(text)
    if not (a and clcd):
        raise ValueError(f"Could not parse alpha / CLff / CDff from {path}")
    return {
        "alpha": float(a.group(1)),
        "CL": float(clcd.group(1)),
        "CDi": float(clcd.group(2)),
        "e": float(e.group(1)) if e else float("nan"),
    }


# ---- Main -------------------------------------------------------------------
def main():
    missing_paths = []
    if not AVL_BIN:
        missing_paths.append("avl352.exe (set AVL_BIN env var or place next to 'CDR - Aero stuff/')")
    if not COLLAB_AVL:
        missing_paths.append("ModWing_CorrSpan_LargeSlat.avl (set COLLAB_AVL env var)")
    if missing_paths:
        sys.exit("Could not locate:\n  - " + "\n  - ".join(missing_paths))

    for cfg in CONFIGS:
        key = cfg["key"]
        avl_path = cfg["avl"]
        mach = cfg["mach"]
        deflections = cfg.get("deflections", {})
        if not os.path.isfile(avl_path):
            print(f"[{key}] SKIP — {avl_path} not found")
            continue

        work = os.path.join(HERE, f"avl_work_{key}")
        os.makedirs(work, exist_ok=True)
        print(f"[{key}] sweeping {os.path.basename(avl_path)} "
              f"M={mach}, deflections={deflections}, "
              f"{len(ALPHAS)} alphas -> {work}")
        out_files = run_alpha_sweep(avl_path, ALPHAS, mach, deflections, work)
        rows = []
        for f in out_files:
            try:
                rows.append(parse_forces(f))
            except (ValueError, FileNotFoundError) as e:
                print(f"  WARN: {e}")

        rows.sort(key=lambda r: r["alpha"])
        csv_path = os.path.join(HERE, f"alpha_sweep_{key}.csv")
        with open(csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["alpha", "CL", "CDi", "e"])
            w.writeheader()
            w.writerows(rows)
        print(f"[{key}] wrote {csv_path}  ({len(rows)} points)")


if __name__ == "__main__":
    main()
