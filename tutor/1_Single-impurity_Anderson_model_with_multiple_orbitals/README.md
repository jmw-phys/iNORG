# Example 1: One-shot solution of a multi-orbital impurity problem in DMFT

Three-orbital Anderson impurity model (mode 1) with U = 8, J = 0.8, Ising-type interaction and
beta = 200. The hybridization function is read from `Delta.inp` and fitted with 11 bath sites per
orbital. See Example 1 of the iNORG paper and the
[online manual](https://jmw-phys.github.io/iNORG/tutorials/siam-oneshot/).

## Run

```bash
cd tutor/1_Single-impurity_Anderson_model_with_multiple_orbitals
mpiexec -n <N> ../../inorg > run.log 2>&1
```

The inputs are `PARAMS.norg` and `Delta.inp`. Replace `<N>` with the number of physical cores (at least 2).
The example has no `ose_hop` file: the solver prints a warning and fits from its generated initial bath.

## Reference run

- Build: Intel oneAPI 2024.0 (icpx 2024.0.2) through the Open MPI 5.0.0 wrapper `mpicxx`, Intel MKL 2024.0.
- 216 MPI processes on 6 nodes (2 x Intel Xeon Gold 6140, 36 cores each): about 54 min.
  With 24 processes on one node the same example takes about 2 h.
- `important_output/` holds `Sig.out`, `Gf.out`, `g0imp.txt`, `hb_read.txt`, `hb_fit.txt`, `h0.txt`,
  `nmat.txt`, the fitted `ose_hop`, `log.norg`, `run.log`, and two figures
  (`hybridization_fit.png`, `self_energy.png`).
- The reference run started its bath fit from an older seed file, not from the generated initial bath.
  A run of the current input can therefore give slightly different numbers.

## How to compare

- The fit error printed in `run.log` (`err = ...`) should be about 2e-4 for each orbital.
- The occupations per spin-orbital in `nmat.txt` should be about 0.330, 0.342 and 0.366 for orbitals 1, 2 and 3.
- Im Sigma in `Sig.out` (columns 3, 5 and 7) should agree with the reference within about 1% for omega_n <= 2.
- The bath fit starts from random points, and the result depends on the number of processes,
  so the numbers do not agree bit for bit.
