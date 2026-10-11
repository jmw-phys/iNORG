# Natural orbitals renormalization group

NORG finds a compact representation of an interacting ground state by optimizing its single-particle orbitals. You retain the occupation fluctuations that matter most and restrict those with small weight. The orbital basis and the many-body state are refined together.

## Natural orbitals

For a normalized many-body state $|\Psi\rangle$, define the single-particle density matrix

$$
D_{ij}=\langle\Psi|c_i^\dagger c_j|\Psi\rangle.
$$

Its eigenvectors define the natural orbitals. Its eigenvalues are their occupations $n_g$. A spin-orbital can contain at most one electron, so $0\leq n_g\leq1$. A spatial-orbital occupation summed over both spins instead lies between zero and two. The NOOC rules use the spin-orbital convention.

Use a unitary matrix $U$ to change the orbital basis:

$$
d_g^\dagger=\sum_i (U^\dagger)_{gi}c_i^\dagger,
\qquad
d_g=\sum_i U_{ig}c_i.
$$

The inverse relations are

$$
c_i^\dagger=\sum_g U_{ig}d_g^\dagger,
\qquad
c_i=\sum_g(U^\dagger)_{gi}d_g.
$$

Unitarity preserves the fermion anticommutation relations. Rotating the full orbital space alone is a change of representation; the approximation enters when you restrict the allowed many-body configurations in that representation.

<!-- iNORG paper: Model and methods / Natural orbitals, inorg.tex:160, label eq2:transformation. -->

## Why occupations guide truncation

Separate the determinants in a wave function according to whether orbital $g$ is empty or occupied:

$$
|\Psi\rangle=\sum_i h_i^\circ|\phi_i^{\circ,g}\rangle
+\sum_j h_j^\bullet|\phi_j^{\bullet,g}\rangle.
$$

Orthonormality and normalization give

$$
\sum_i|h_i^\circ|^2+\sum_j|h_j^\bullet|^2=1,
\qquad
n_g=\sum_j|h_j^\bullet|^2.
$$

If you freeze a nearly empty orbital to be empty, the discarded component has squared norm $n_g$. If you freeze a nearly full orbital to be full, it has squared norm $1-n_g$. This motivates the occupation weight $\min(n_g,1-n_g)$ used by the automatic division selection. It measures a single-orbital discarded weight, not a universal bound on energy or Green's-function error.

| Occupation | Role | Treatment |
| --- | --- | --- |
| Close to one | Nearly full, inactive | Limit holes relative to full occupancy. |
| Far from zero and one | Active | Keep its occupation fluctuations explicitly. |
| Close to zero | Nearly empty, inactive | Limit electrons relative to the empty state. |

NOOC can allow a small number of excitations in an inactive group instead of freezing every orbital completely. This provides a hierarchy of spaces that you can enlarge to test accuracy. See [NOOC rules and input syntax](../implementation/nooc.md).

<!-- iNORG paper: Natural orbitals, label eq2:ketPsi; Natural orbital occupancy constraint. Occupation weight: src/api_edmft.cpp:673. -->

## Renormalization cycle

The method repeats the following steps:

1. Start from an orbital basis and a choice of occupation constraints.
2. Enumerate the allowed determinants and construct the Hamiltonian in that subspace.
3. Solve for its ground state with Lanczos.
4. Compute the single-particle density matrix of that state.
5. Diagonalize the relevant density-matrix blocks and rotate the orbitals.
6. Solve again in the updated basis. Check energy and occupations, and revise the subspace when required.

The variational energy is the lowest energy within the chosen subspace. Loosening constraints enlarges the variational freedom for a fixed orbital basis. In practice, compare observables after both the orbital optimization and the enlarged-space solve have settled.

The [implementation workflow](../implementation/workflow.md) distinguishes this orbital optimization from the preliminary particle-number search and the outer DMFT loop. Their convergence tests serve different purposes.

<!-- iNORG paper: Renormalization strategy, inorg.tex:231; NOOC variational discussion. -->

## Bath-only and full-orbital rotations

The paper describes two strategies:

| Rotation | Density matrix to diagonalize | Effect on interactions |
| --- | --- | --- |
| All orbitals | Full single-particle density matrix | Local impurity interactions spread over the transformed orbitals. |
| Bath only | Bath block, with the impurity basis fixed | Interactions stay on the impurity; the transformed bath can acquire additional one-body couplings. |

Bath-only rotation is the default in this repository. Keeping the interaction local avoids the growth of a general four-index interaction tensor under rotation. The paper describes typical bath-size scaling of $O(N_{\mathrm{bath}}^3)$ for this strategy and $O(N^4)$ for full-orbital rotation. These are method-level scaling estimates; the selected Hilbert-space dimension and constraints still determine the cost of a particular run.

In the current implementation, the bath density matrices of the two spins are averaged within each band before extracting the rotation. Separate magnetic input channels therefore do not imply independent spin-dependent orbital optimization.

<!-- iNORG paper: Renormalization strategy, inorg.tex:245. src/prmtr.cpp:69; src/densitymatrix.cpp:72. -->

## References

- P.-O. Löwdin (1955), “Quantum theory of many-particle systems. I,” *Physical Review* **97**, 1474.
- P.-O. Löwdin and H. Shull (1956), “Natural Orbitals in the Quantum Theory of Two-Electron Systems,” *Physical Review* **101**, 1730–1739.
- R.-Q. He and Z.-Y. Lu (2014), “Quantum renormalization groups based on natural orbitals,” *Physical Review B* **89**, 085108.
- J.-M. Wang et al. (2025), “Ab initio dynamical mean field theory with natural orbitals renormalization group impurity solver,” *npj Computational Materials* **11**, 86.
- J.-M. Wang et al., *iNORG: An open-source quantum impurity solver package based on the natural orbitals renormalization group*, “Natural orbitals renormalization group method.”
