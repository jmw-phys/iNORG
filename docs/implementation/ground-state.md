# Ground-state solution

Each orbital-optimization sweep solves the Hamiltonian in its current NOOC space. The [Lanczos recurrence](../theory/lanczos.md) supplies a small tridiagonal eigenproblem. iNORG then reconstructs the low-energy state in the selected many-body basis and uses it to calculate the single-particle density matrix.

## Sparse products and two-pass reconstruction

The eigensolver assembles the Hamiltonian with `find_hmlt_V2`, generates a random starting vector, and calls the MPI Lanczos routine. Dot products and normalization use reductions across ranks.

Conceptually, a ground-state solve follows these steps:

1. Normalize the starting vector and calculate its Hamiltonian product.
2. Generate the tridiagonal coefficients with the three-term recurrence.
3. Diagonalize the growing tridiagonal matrix when checking convergence.
4. Obtain its lowest eigenvalue and eigenvector.
5. Reconstruct the many-body state as the corresponding combination of Lanczos vectors.

The implementation makes two passes. The first retains the tridiagonal coefficients. The second restarts from the same initial vector, regenerates the recurrence, and accumulates the state using the small eigenvector's components. It therefore avoids retaining every full-dimensional Krylov vector for this ground-state reconstruction.

<!-- iNORG paper: Ground-state solution based on the Lanczos algorithm, inorg.tex:678, algorithm lanczos. Code: src/operator.cpp:304; src/operator.cpp:317; src/operator.cpp:326; src/gen/lanczos.h:84; src/gen/lanczos.h:94; src/gen/lanczos.h:119; src/gen/lanczos.h:123. -->

## Residual convergence

For an orthonormal Lanczos basis after $m$ steps, let $\theta$ be a Ritz value and $y$ its normalized eigenvector in that basis. With $v=Q_my$,

$$
\|Hv-\theta v\|_2=\beta_m|y_m|.
$$

The paper's convergence criterion is

$$
\beta_m|y_m|\leq\epsilon.
$$

This asks whether the projected eigenpair also approximately solves the original selected-space problem. It does not test whether the selected space is large enough for the physical observable of interest.

The current routine evaluates a last-component residual estimator using its stored tridiagonal coefficients. Its standard threshold is $10^{-9}$; its fast path uses $10^{-8}$. Checks begin after the internal Krylov index reaches 60. The caller passes a nominal Krylov limit of 1000; the routine warns at that limit and has an additional guard twenty index steps later. These are compiled implementation settings, not `PARAMS.norg` keywords. A limit warning needs inspection even when the routine returns a state.

<!-- iNORG paper: labels eq3:lanczosResidual and eq3:ResidualCriterion. Thresholds and iteration guards are code values, not paper benchmark claims: src/gen/lanczos.h:18; src/gen/lanczos.h:25; src/gen/lanczos.h:98; src/gen/lanczos.h:100; src/operator.cpp:327. -->

## Orthogonality and degeneracy

Residual convergence and numerical orthogonality are separate issues. The paper discusses careful summation, optional reorthogonalization, spectral shifts, and block Lanczos as ways to improve difficult calculations. These describe numerical strategies; they are not all user-selectable options in the current program.

The current routine obtains multiple low-energy states sequentially. It orthogonalizes subsequent searches against states already found and can expand the requested eigenpair count when it detects a degenerate level. This is not a block-Lanczos implementation. The ground-state routine also does not retain a full set of Krylov vectors for general reorthogonalization within its first search.

<!-- Paper wording: optional partial reorthogonalization, Kahan summation, spectral shifts, and block methods. Current sequential deflation: src/gen/lanczos.h:34; src/gen/lanczos.h:76; src/gen/lanczos.h:86; src/gen/lanczos.h:108; src/gen/lanczos.h:144. -->

The input name `pred_gs_deg` needs care:

| Input behavior | Mode 0 | Mode 1 |
| --- | --- | --- |
| Sets the density-matrix branch selector | Yes | Yes |
| Initializes the stored degeneracy from this value | Yes | No; the setup resets the stored value. |
| Selects the explicit multi-state density-matrix branch | Only when the value is exactly `1` | Only when the value is exactly `1` |

Thus `pred_gs_deg 2` does not impose a twofold ground state. Later eigensolves can increase the stored degeneracy. Preserve the distinction between the selector, the number of states sought internally, and the detected degeneracy. The [parameter reference](../parameters.md#degeneracy-and-magnetic-input) gives the exact input behavior.

<!-- Paper wording: pred_gs_deg is the expected ground-state degeneracy. Code: src/main.cpp:248; src/main.cpp:249; src/api_edmft.cpp:291; src/api_edmft.cpp:498; src/densitymatrix.cpp:102; src/densitymatrix.cpp:111; src/operator.cpp:331. -->

## Returning to orbital optimization

The state supplies $D_{ij}=\langle\psi|c_i^\dagger c_j|\psi\rangle$. Its bath blocks determine the next natural-orbital rotation. The subsequent NORG sweep changes Hamiltonian coefficients and solves another selected-space eigenproblem.

Monitor the reported energy and occupation changes across sweeps. The NORG stopping function combines energy thresholds, stability history, and an occupation test; it differs from the Lanczos residual criterion. It can also reach the iteration cap. The [workflow page](workflow.md) explains the distinction from outer DMFT convergence.

<!-- src/densitymatrix.cpp:100; src/densitymatrix.cpp:114; src/norg.cpp:68; src/norg.cpp:92; src/norg.cpp:529. -->

## References

- G. H. Golub and C. F. Van Loan (2013), *Matrix Computations*, JHU Press.
- K. Wu, A. Canning, H. D. Simon, and L.-W. Wang (1999), “Thick-restart Lanczos method for electronic structure calculations,” *Journal of Computational Physics* **154**, 156–173.
- J.-M. Wang et al., *iNORG: An open-source quantum impurity solver package based on the natural orbitals renormalization group*, “Ground-state solution based on the Lanczos algorithm.”
