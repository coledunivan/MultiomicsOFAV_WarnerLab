#!/usr/bin/env python3
"""
run_pipeline.py — run the whole coral multiomics pipeline end to end.

    python run_pipeline.py                 # full run
    python run_pipeline.py --skip-atac     # RNA-only (fast; skips ATAC DESeq2)
    python run_pipeline.py --from 3        # resume from a given step
    python run_pipeline.py --only 5        # run a single step
    python run_pipeline.py --config my.yaml

Steps:
    1  01_differential.py   DE (RNA) + DA (ATAC)
    2  02_integration.py    RNA-ATAC links + candidates
    3  03_clustering.py     fuzzy temporal programs
    4  04_enrichment.py     GO enrichment
    5  05_figures.py        Figure1..5
"""
import argparse
import subprocess
import sys
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PIPE = os.path.join(HERE, "pipeline")

STEPS = [
    ("1", "01_differential.py"),
    ("2", "02_integration.py"),
    ("3", "03_clustering.py"),
    ("4", "04_enrichment.py"),
    ("5", "05_figures.py"),
]


def run(script, extra):
    cmd = [sys.executable, os.path.join(PIPE, script)] + extra
    print(f"\n>>> {' '.join(os.path.basename(c) for c in cmd[1:2])} {' '.join(extra)}")
    r = subprocess.run(cmd)
    if r.returncode != 0:
        print(f"\n!!! step {script} failed (exit {r.returncode}); stopping.")
        sys.exit(r.returncode)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--skip-atac", action="store_true",
                    help="skip the slow ATAC DESeq2 (step 1b)")
    ap.add_argument("--from", dest="start", default="1",
                    help="resume from this step number")
    ap.add_argument("--only", default=None, help="run only this step number")
    args = ap.parse_args()

    cfg_args = (["--config", args.config] if args.config else [])

    steps = STEPS
    if args.only:
        steps = [s for s in STEPS if s[0] == args.only]
    else:
        steps = [s for s in STEPS if s[0] >= args.start]

    for num, script in steps:
        extra = list(cfg_args)
        if script == "01_differential.py" and args.skip_atac:
            extra.append("--skip-atac")
        run(script, extra)

    print("\n" + "=" * 60)
    print("Pipeline complete. See results/tables/ and results/figures/.")
    print("=" * 60)


if __name__ == "__main__":
    main()
