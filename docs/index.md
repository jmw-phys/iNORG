# iNORG user manual

iNORG is an open source quantum impurity solver based on the **natural orbitals renormalization group** (NORG). It solves interacting impurity models coupled to a finite bath and computes impurity occupations, Green's functions, and self-energies. It is written in C++17 and uses MPI and Intel MKL.

The solver can solve a multi-orbital impurity with a supplied hybridization function, or repeat impurity solutions within a Bethe lattice dynamical mean-field theory (DMFT) loop. The supplied examples use density-density interactions and imaginary-frequency data.

The [Theory section](theory/impurity-model.md) develops the impurity model, natural orbitals, many-body representation, and Lanczos method. The [Implementation section](implementation/workflow.md) follows the calculation from bath fitting through occupation constraints and sparse operators to Green's functions. Equations describe the method, while input examples and source notes identify the behavior of this version.

## Choose a run mode

Put `PARAMS.norg` in the directory from which you launch `inorg`. Its `mode` entry selects the calculation.

| Setting | Calculation | Start here |
| --- | --- | --- |
| `mode 0` | Multi-orbital Bethe lattice DMFT, with repeated bath fits and impurity solutions. | [Bethe DMFT tutorial](tutorials/bethe-dmft.md) |
| `mode 1` | One-shot impurity solution using `Delta.inp`. There is no lattice self-consistency loop. | [One-shot tutorial](tutorials/siam-oneshot.md) |
| No `mode` key | The eDMFT interface path. `PARAMS.norg` must still exist. | Interface pointer below |

<!-- src/main.cpp:21; src/main.cpp:37; src/main.cpp:39; src/main.cpp:44; src/main.cpp:72; src/main.cpp:102; src/main.cpp:125 -->

Use the literal value `0` or `1` for these tutorials. Other mode values also route to the eDMFT path; they do not select further tutorial modes.

The external eDMFT workflow is in [`api/eDMFT`](https://github.com/jmw-phys/iNORG/tree/main/api/eDMFT). It contains Python drivers and interface utilities for embedding the solver in a larger calculation. This manual covers the two standalone modes; it does not yet describe setup of that interface.

## First calculation

1. [Build iNORG](installation.md) with a compatible MPI and MKL toolchain.
2. Follow the [quick start](quickstart.md). The one-shot example has no `ose_hop`; the solver warns and fits from its generated initial bath.
3. Use the [parameter reference](parameters.md) and [output-file reference](output-files.md) to inspect your run.

Source code and issue tracking are on [GitHub](https://github.com/jmw-phys/iNORG). See [Citing iNORG](citing.md) for the method paper and the [paper list](paper-list.md) for applications.

## Future developments

The iNORG paper identifies further adaptive occupancy-constraint schemes, finite-temperature extensions, two-particle responses, cluster DMFT, and Python interfaces for integration with DFT+DMFT workflows as future work. These are development directions; the current standalone manual covers ground-state, single-site impurity calculations.

<!-- iNORG paper: Future developments, inorg.tex:1050; Program summary restrictions. -->

## References

- J.-M. Wang et al., *iNORG: An open-source quantum impurity solver package based on the natural orbitals renormalization group*, “Future developments.”
