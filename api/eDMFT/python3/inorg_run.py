#!/usr/bin/env python3
"""
inorg_run.py -- Bridge script to run iNORG solver from eDMFT.

coded by Jia-Ming Wang (jiaming.w@rutgers.edu, Rutgers University, U.S.) date 2025 - now

Workflow:
  1. Read mpi_prefix.dat for MPI command
  2. Read Trans.dat for sigind matrix
  3. Convert input.tsv -> solver_input.tsv
     (site-basis sigind-compressed -> K-basis nband x nband, Python complex format)
  4. Run: {mpi_prefix} {inorg_exe}
  5. Postprocess: io/seimp.txt, io/gfimp.txt -> Sig.out, Gf.out
"""

import sys
import os
import subprocess
import numpy as np

# Add this directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from inorg_edmft_io import (
    read_sigind_matrix,
    read_params_norg,
    convert_edinput_to_solver_input,
    postprocess_from_norg,
)


def main():
    print("=" * 60)
    print("inorg_run.py: iNORG solver bridge for eDMFT")
    print("=" * 60)

    # 1. Read MPI prefix
    mpi_prefix = ""
    if os.path.exists("mpi_prefix.dat"):
        with open("mpi_prefix.dat", 'r') as f:
            mpi_prefix = f.read().strip()
        print(f"MPI prefix: {mpi_prefix}")
    else:
        print("Warning: mpi_prefix.dat not found, running without MPI")

    # 2. Read Trans.dat for sigind matrix
    sigind = None
    nband = None
    if os.path.exists("Trans.dat"):
        sigind = read_sigind_matrix("Trans.dat")
        # Determine nband from sigind diagonal
        ntot = sigind.shape[0]
        half = ntot // 2
        nband = np.count_nonzero(np.diag(sigind[:half, :half]))
        print(f"Read Trans.dat: ntot={ntot}, nband={nband}")
        print(f"  sigind shape: {sigind.shape}")
    else:
        print("Warning: Trans.dat not found, using defaults")

    # 3. Check that PARAMS.norg exists and read sym parameter
    if not os.path.exists("PARAMS.norg"):
        print("ERROR: PARAMS.norg not found")
        sys.exit(1)
    params_norg = read_params_norg("PARAMS.norg")
    sym = params_norg.get('sym', 'psc')
    print(f"PARAMS.norg found, sym={sym}")

    # 4. Convert input.tsv -> solver_input.tsv
    if os.path.exists("input.tsv"):
        print("\nConverting input.tsv -> solver_input.tsv ...")
        convert_edinput_to_solver_input(
            "input.tsv", "solver_input.tsv",
            sigind=sigind, transform_to_K=True, sym=sym)
    else:
        print("ERROR: input.tsv not found")
        sys.exit(1)

    # 5. Find and run C++ solver
    # Look for inorg binary in WIEN_DMFT_ROOT or same directory as this script
    wien_root = os.environ.get("WIEN_DMFT_ROOT",
                                os.path.dirname(os.path.abspath(__file__)))
    inorg_exe = os.path.join(wien_root, "inorg")
    if not os.path.exists(inorg_exe):
        # Try current directory
        if os.path.exists("./inorg"):
            inorg_exe = "./inorg"
        else:
            print(f"ERROR: inorg binary not found at {inorg_exe} or ./inorg")
            sys.exit(1)

    cmd = f"{mpi_prefix} {inorg_exe}".strip()
    print(f"\nRunning C++ solver: {cmd}")
    print("-" * 60)
    ret = subprocess.call(cmd, shell=True)
    print("-" * 60)

    if ret != 0:
        print(f"ERROR: iNORG solver returned exit code {ret}")
        print("Check nohup_imp.out for details.")
        sys.exit(ret)

    print("\nC++ solver completed successfully")

    # 6. Postprocess: io/seimp.txt, io/gfimp.txt -> Sig.out, Gf.out
    print("\nPostprocessing output ...")
    postprocess_from_norg(".", sigind=sigind, nband=nband, sym=sym)

    print("\n inorg_run.py completed successfully")


if __name__ == "__main__":
    main()
