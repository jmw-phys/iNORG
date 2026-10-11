# Installation

The supplied build targets Linux and produces the standalone MPI executable. iNORG is released under the GNU Affero General Public License, version 3.

<!-- README.md:45; iNORG paper: Installation and usage / Obtaining iNORG; Overall workflow / Development Environment. -->

## Obtain the source

Clone the public repository and enter its root:

```bash
git clone https://github.com/jmw-phys/iNORG.git
cd iNORG
```

Build and run instructions below use this checkout's makefile and tutorial directories.

<!-- The iNORG paper lists QTMEC-RUC/iNORG. This manual follows the public repository identified in README.md:17: jmw-phys/iNORG. -->

## Requirements

- A C++17 compiler and an MPI C++ compiler wrapper named `mpicxx`.
- Intel MKL headers and libraries. The supplied makefile uses the compiler option `-qmkl` at compile and link time.
- GNU Make and an MPI launcher such as `mpiexec`.

The `mpicxx` wrapper must invoke a compiler that accepts `-qmkl`, with MKL configured in its environment. A wrapper around a compiler that rejects this flag will not build with the supplied settings. Use the MPI and compiler environment provided for your machine, and use the same MPI installation to build and run.

<!-- makefile:7; makefile:10; makefile:11; makefile:17; makefile:46 -->

The paper lists modern GCC, Clang, and Intel compilers as C++17 examples. C++17 support alone is not sufficient for the supplied makefile: its `-qmkl` option also has to work. Using another compiler requires compatible MKL compile and link settings.

<!-- Paper wording: GCC 7+, Clang 6+, Intel C++ Compiler 19+ as examples. Current build flags: makefile:10; makefile:11. -->

## Build

From the repository root, with the compiler, MPI, and MKL environment enabled:

```bash
make clean && make
test -x inorg
```

The executable is **`inorg` in the repository root**. The default make target compiles and links it; it does not launch a calculation. `make clean` removes previous objects, dependency files, and the executable.

<!-- makefile:7; makefile:24; makefile:46; makefile:82 -->

Run Make in the checkout root. The supplied makefile resolves `src/` relative to that directory; it does not provide a separate build-directory configuration step. It already requests parallel compilation. You can set the parallel job count explicitly with `make -j N`, replacing `N` with the available build parallelism.

<!-- Paper wording: create build/, then make, with executable in build/src or build/bin. Actual paths and parallelism: makefile:1; makefile:6; makefile:7; makefile:22; makefile:46. -->

## MPI execution

Use **at least two MPI ranks**. The program aborts with fewer ranks because bath fitting uses a master-worker scheme. Two ranks are a minimum, not a runtime or memory estimate. Larger impurity spaces can require substantial resources.

<!-- src/main.cpp:19; src/main.cpp:21; src/bath.cpp:138; src/bath.cpp:212 -->

Launch from a run directory containing `PARAMS.norg` and the inputs for its mode. Replace `<N>` with the number of physical cores on your machine (at least 2). For example, from either tutorial directory:

```bash
mpiexec -n <N> ../../inorg
```

For an InfiniBand system, build with an MPI implementation that includes UCX support. For example, Open MPI 5 can explicitly select UCX with `mpiexec --mca pml ucx -n <N> ../../inorg`. Follow your site's allocation and launch rules. See the [Open MPI InfiniBand and UCX documentation](https://docs.open-mpi.org/en/v5.0.x/tuning-apps/networking/ib-and-roce.html).

Continue with the [quick start](quickstart.md) to run the supplied one-shot example. It ships without `ose_hop`, so the solver warns and fits from its generated initial bath.

Keep the runtime MPI and MKL libraries available in the launch environment. Run from a directory containing the intended input files, redirect the log, and inspect the fit and iteration history. Restart support is specific to each mode and file; see [output files](output-files.md#logs-restarts-and-completion) before reusing a previous run directory.

<!-- iNORG paper: Installation and usage / Running iNORG. Concrete restart behavior: src/api_edmft.cpp:46; src/api_edmft.cpp:53; src/dmft_bethe.cpp:33; src/dmft_bethe.cpp:47. -->

## Preview and publish this manual

The documentation build is independent of the solver build. With Python 3 and `pip` available, run from the repository root:

```bash
python -m pip install mkdocs-material
python -m mkdocs serve
```

Open the local URL printed by MkDocs. To check a static build, use `python -m mkdocs build --strict`.

The documentation workflow builds and deploys on pushes to `main` that change `docs/**` or `mkdocs.yml`. It can also be run manually through `workflow_dispatch`. The workflow publishes the `gh-pages` branch. A repository maintainer must select **Deploy from a branch**, branch **gh-pages**, folder **/(root)** in the repository's Pages settings. See [Material for MkDocs: publishing](https://squidfunk.github.io/mkdocs-material/publishing-your-site/).

## References

- J.-M. Wang et al., *iNORG: An open-source quantum impurity solver package based on the natural orbitals renormalization group*, “Installation and usage.”
