#!/usr/bin/env python
"""
compute_psa.py
==============

Compute 3D polar surface area (PSA) for a ligand.

• default  : QikProp (Schrödinger) reference implementation
• --method open : placeholder for open-source route

Behavior:
1. Copy input SDF into a temp directory
2. Run $SCHRODINGER/qikprop on it
3. Detect JobId and either use jobwait (if available) or poll for output files
4. Parse PSA from modern .CSV or legacy .qprop
"""
import os
import sys
import csv
import shutil
import subprocess
import argparse
import tempfile
import time
from pathlib import Path


def psa_qikprop(sdf_path: str) -> float:
    sch = os.environ.get("SCHRODINGER")
    if not sch:
        raise RuntimeError("SCHRODINGER env var not set")

    sdf = Path(sdf_path).resolve()
    if not sdf.exists():
        raise FileNotFoundError(f"Input file not found: {sdf}")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_sdf = Path(tmpdir) / sdf.name
        shutil.copy(sdf, tmp_sdf)

        # Submit QikProp job
        result = subprocess.run(
            [Path(sch) / "qikprop", tmp_sdf.name],
            cwd=tmpdir,
            text=True,
            capture_output=True,
            check=True,
        )

        base = tmp_sdf.with_suffix("")
        csv_file = base.with_suffix(".CSV")
        qprop_file = base.with_suffix(".qprop")

        # Wait for job to finish
        # Extract job_id if Job Control is active
        job_id = None
        for line in result.stdout.splitlines():
            if line.startswith("JobId:"):
                job_id = line.split()[1]
                break

        if job_id:
            # Try jobwait locations
            jw1 = Path(sch) / "util" / "jobwait"
            jw2 = Path(sch) / "util" / "jobcontrol" / "jobwait"
            if jw1.exists():
                subprocess.run([jw1, job_id], cwd=tmpdir, check=True)
            elif jw2.exists():
                subprocess.run([jw2, job_id], cwd=tmpdir, check=True)
            else:
                print("[info] jobwait not found; polling for output…", file=sys.stderr)
                # Poll for output
                for _ in range(60):  # up to 2 minutes
                    if csv_file.exists() or qprop_file.exists():
                        break
                    time.sleep(2)

        # Parse PSA value
        if csv_file.exists():
            with csv_file.open(newline="") as fh:
                reader = csv.DictReader(fh)
                row = next(reader)
                return float(row.get("PSA"))

        if qprop_file.exists():
            with qprop_file.open() as fh:
                for l in fh:
                    if l.strip().startswith("PSA"):
                        return float(l.split()[1])

        raise RuntimeError("QikProp finished but no .CSV or .qprop output found")


def psa_open(_sdf: str) -> float:
    raise NotImplementedError("Open-source PSA route not implemented yet")


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="compute_psa.py",
        description="Compute 3D polar surface area (PSA) using QikProp or open-source route",
    )
    parser.add_argument("sdf", nargs="?", help="Input SDF file with conformers")
    parser.add_argument(
        "--method", choices=["qikprop", "open"], default="qikprop",
        help="Engine to use; default is qikprop"
    )
    args = parser.parse_args(argv)

    if not args.sdf:
        parser.print_usage(sys.stderr)
        print("Error: the <sdf> positional argument is required.", file=sys.stderr)
        sys.exit(1)

    try:
        if args.method == "qikprop":
            psa = psa_qikprop(args.sdf)
        else:
            psa = psa_open(args.sdf)
        print(f"3D PSA = {psa:.1f} Å²   ({args.method})")
    except Exception as e:
        print(f"[ERROR] {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()

