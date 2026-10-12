# Example 2: Multi-orbital DMFT calculation on the Bethe lattice

Two-orbital Hubbard model on the Bethe lattice at half filling (mode 0), with U = 3, U' = 2.7,
hoppings t = 0.5 (wide band) and t = 0.25 (narrow band). See Example 2 of the iNORG paper and
the [online manual](https://jmw-phys.github.io/iNORG/tutorials/bethe-dmft/).

## Run

```bash
cd tutor/2_Multi-orbital_model_on_Bethe_lattice_with_DMFT
mpiexec -n <N> ../../inorg > run.log 2>&1
```

The only input is `PARAMS.norg`. Replace `<N>` with the number of physical cores (at least 2).

## Reference run

- Build: Intel oneAPI 2024.0 (icpx 2024.0.2) through the Open MPI 5.0.0 wrapper `mpicxx`, Intel MKL 2024.0.
- 36 MPI processes on one node (2 x Intel Xeon Gold 6140, 36 cores): about 50 min to iteration 30,
  about 1 h 45 min to the end of the run (62 iterations).
- `important_output/` holds the state after iteration 30 (`zic030.*`), `log.txt` and `run.log` up to
  iteration 30, `log.parameters.txt`, and two figures (`convergence.png`, `self_energy.png`).

## How to compare

- From iteration 1 on, every spin-orbital occupation in `log.txt` (columns 4 to 7) must equal 0.5 within about 1e-9. (Row 0 is the initial placeholder.)
- The error in `log.txt` (column 2) should fall from about 4 at iteration 1 to about 2e-3 at iteration 30.
- Im Sigma in `zic030.mb.seimp.txt` (columns 5 and 6) should agree with the reference within about 1%
  for omega_n <= 2. The narrow band diverges as omega_n -> 0 (Mott insulator); the wide band does not.
- The bath fit starts from random points, and the result depends on the number of processes,
  so the numbers do not agree bit for bit.
