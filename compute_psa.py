#!/usr/bin/env python
"""
compute_psa.py
==============

Compute 3D polar surface area (PSA).

• default  : QikProp (Schrödinger)            -- reliable, reference value
• --method open : placeholder for open-source route (to be implemented)

The script:
1) Copies the input SDF to a temporary directory.
2) Runs $SCHRODINGER/qikprop (blocking if Job Control is used).
3) Reads PSA from   *.CSV   (modern) or  *.qprop  (legacy) output.
"""

import os, shutil, subprocess, argparse, tempfile, csv
from pathlib import Path

# ─────────────────────────── QikProp branch ─────────────────────────── #
def psa_qikprop(sdf_path: str) -> float:
    sch = os.environ.get("SCHRODINGER")
    if not sch:
        raise RuntimeError("SCHRODINGER env var not set")

    sdf_path = Path(sdf_path).resolve()
    if not sdf_path.exists():
        raise FileNotFoundError(sdf_path)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_sdf = Path(tmpdir) / sdf_path.name
        shutil.copy(sdf_path, tmp_sdf)

        # Run QikProp, capture stdout to see if Job Control submits a job
        result = subprocess.run(
            [Path(sch) / "qikprop", tmp_sdf.name],
            cwd=tmpdir,
            text=True,
            capture_output=True,
            check=True,
        )

        # Wait for Job Control job if needed
        for line in result.stdout.splitlines():
            if line.startswith("JobId:"):
                job_id = line.split()[1]
                subprocess.run([Path(sch) / "util" / "jobwait", job_id],
                               cwd=tmpdir, check=True)
                break

        base       = tmp_sdf.with_suffix("")
        csv_file   = base.with_suffix(".CSV")
        qprop_file = base.with_suffix(".qprop")

        if csv_file.exists():                       # modern versions
            with csv_file.open(newline="") as fh:
                row = next(csv.DictReader(fh))
                return float(row["PSA"])

        if qprop_file.exists():                     # legacy versions
            with qprop_file.open() as fh:
                for line in fh:
                    if line.strip().startswith("PSA"):
                        return float(line.split()[1])

        raise RuntimeError("QikProp finished but produced neither CSV nor qprop")


# ────────────────────────── open-source placeholder ──────────────────── #
def psa_open(_sdf: str) -> float:
    raise NotImplementedError("Open-source PSA route not yet implemented")


# ───────────────────────────────── CLI ───────────────────────────────── #
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("sdf", help="Input SDF file with conformers")
    parser.add_argument("--method", choices=["qikprop", "open"],
                        default="qikprop",
                        help="Engine (default=qikprop)")
    args = parser.parse_args()

    try:
        psa_val = psa_qikprop(args.sdf) if args.method == "qikprop" else psa_open(args.sdf)
        print(f"3D PSA = {psa_val:.1f} Å²   ({args.method})")
    except Exception as exc:
        print(f"[ERROR] {exc}")

