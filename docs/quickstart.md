# Quick start

Start with the one-shot impurity example after [building `inorg`](installation.md). These commands assume a fresh checkout and start at the repository root.

Each example directory carries the number of the example in the iNORG paper:

| Paper example | Directory | Mode | Required fresh-run inputs |
| --- | --- | --- | --- |
| Example 1: one-shot multi-orbital impurity | `tutor/1_Single-impurity_Anderson_model_with_multiple_orbitals` | `1` | `PARAMS.norg`, `Delta.inp` |
| Example 2: Bethe lattice DMFT | `tutor/2_Multi-orbital_model_on_Bethe_lattice_with_DMFT` | `0` | `PARAMS.norg` |

Both modes can read an optional `ose_hop` bath seed. Mode 1 fits a fixed hybridization; mode 0 obtains successive targets from lattice self-consistency. The [solver workflow](implementation/workflow.md) explains the inner orbital optimization shared by both.

<!-- iNORG paper: Examples / Example 1 and Example 2. src/main.cpp:39; src/main.cpp:44; src/api_edmft.cpp:26; src/api_edmft.cpp:32; src/dmft_bethe.cpp:33; src/dmft_bethe.cpp:60. -->

```bash
cd tutor/1_Single-impurity_Anderson_model_with_multiple_orbitals
mpiexec -n <N> ../../inorg > run.log 2>&1
```

Replace `<N>` with the number of physical cores on your machine (at least 2). For reference, this example takes about 2 hours on 24 cores. The executable is two directories above the tutorial, not inside `tutor/`.

The example ships without an `ose_hop` file, so the solver prints a warning and fits from its generated initial bath. After the run, the fitted `ose_hop` can serve as a seed for a later run if its bath counts still match `PARAMS.norg`. Save previous results before rerunning: mode 1 overwrites its bath and numerical output files. For separate experiments, copy the tutorial directory into a fresh run directory and adjust the binary path.

## Inspect the result

The run reads `PARAMS.norg` and `Delta.inp`, fits the bath, searches particle-number sectors, solves the impurity, and writes its Green's function and self-energy.

<!-- src/api_edmft.cpp:19; src/api_edmft.cpp:29; src/api_edmft.cpp:32; src/api_edmft.cpp:47; src/api_edmft.cpp:52; src/api_edmft.cpp:67; src/api_edmft.cpp:79 -->

```bash
tail -n 20 run.log
head -n 5 hb_read.txt hb_fit.txt
head -n 5 Gf.out Sig.out
head -n 12 nmat.txt
```

Compare `hb_fit.txt` with `hb_read.txt` to assess bath discretization. Inspect `Gf.out` and `Sig.out` for finite, smooth imaginary-frequency data. The first block of `nmat.txt` gives impurity occupations; its double-occupancy and derived correlation blocks are not reliable in the current version.

On a fresh run, a warning about missing `Sig.out` is expected: the solver tries to compare with a previous self-energy before writing the new one. This comparison does not make a one-shot run self-consistent. Read [output files](output-files.md) for exact columns and restart behavior.

<!-- src/api_edmft.cpp:73; src/green.cpp:207 -->

The [one-shot tutorial](tutorials/siam-oneshot.md) explains every input line and both input formats. The other supplied example is `tutor/2_Multi-orbital_model_on_Bethe_lattice_with_DMFT`; follow the [Bethe DMFT tutorial](tutorials/bethe-dmft.md) for that mode.

## References

- J.-M. Wang et al., *iNORG: An open-source quantum impurity solver package based on the natural orbitals renormalization group*, “Installation and usage” and “Examples.”
