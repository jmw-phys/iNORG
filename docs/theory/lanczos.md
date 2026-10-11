# Lanczos algorithm

Lanczos finds low-energy eigenstates of a large sparse Hamiltonian through repeated Hamiltonian-vector products. It replaces a full diagonalization with a much smaller projected eigenproblem. In NORG, the large Hamiltonian already acts within the selected occupation-constrained space.

## Krylov subspace

Starting with a normalized vector $|\varphi_0\rangle$, form

$$
\mathcal K_m(H,|\varphi_0\rangle)=
\operatorname{span}\{|\varphi_0\rangle,H|\varphi_0\rangle,
\ldots,H^{m-1}|\varphi_0\rangle\}.
$$

Lanczos generates an orthonormal basis for this subspace without explicitly calculating and storing the powers of $H$. With $\beta_0=0$ and $|\varphi_{-1}\rangle=0$, the recurrence is

$$
\begin{aligned}
\alpha_i&=\langle\varphi_i|H|\varphi_i\rangle,\\
|r_i\rangle&=(H-\alpha_i)|\varphi_i\rangle
-\beta_i|\varphi_{i-1}\rangle,\\
\beta_{i+1}&=\|r_i\|,\qquad
|\varphi_{i+1}\rangle=|r_i\rangle/\beta_{i+1}.
\end{aligned}
$$

If the residual vector vanishes, the generated subspace is invariant and there is no next vector to normalize. In exact arithmetic, the recurrence makes the projected Hamiltonian tridiagonal:

$$
T_m=\begin{pmatrix}
\alpha_0&\beta_1&0&\cdots\\
\beta_1&\alpha_1&\beta_2&\cdots\\
0&\beta_2&\alpha_2&\ddots\\
\vdots&\vdots&\ddots&\ddots
\end{pmatrix}.
$$

<!-- iNORG paper: Lanczos algorithm, inorg.tex:283, label eq:lancKryl; continued-fraction section for T. Residual vector notation expands the paper's three-term recurrence. -->

## Ritz values and ground-state reconstruction

Diagonalizing $T_m$ gives approximate eigenvalues, called Ritz values. If $Q_m$ contains the Lanczos vectors as columns and $y$ is the lowest eigenvector of $T_m$, then $Q_my$ approximates the ground state in the original selected space.

The paper gives roughly 200 iterations as a typical illustration of rapid extremal-eigenvalue convergence. This is not a fixed input setting or a guarantee for a model. Near-degenerate states, the starting vector, and the selected space affect the work required.

A small change in energy is useful evidence but does not fully characterize an approximate eigenvector. The Ritz residual provides a more direct check:

$$
\|H v-\theta v\|_2=\beta_m|y_m|,
\qquad v=Q_my.
$$

Here $y_m$ is the last component when the components of $y$ are numbered from one. The identity assumes an orthonormal Lanczos basis and the corresponding three-term recurrence. It lets you estimate the residual without another large matrix-vector multiplication.

<!-- iNORG paper: Lanczos algorithm; Ground-state solution based on the Lanczos algorithm, labels eq3:lanczosResidual and eq3:ResidualCriterion. -->

## Numerical accuracy and scope

Finite-precision vectors can lose orthogonality. A residual test monitors convergence of an eigenpair; it does not itself restore orthogonality. Reorthogonalization and careful summation address a different source of error. The paper also discusses spectral shifts and block methods as possible aids for difficult or degenerate spectra.

For the concrete solver, see [ground-state implementation](../implementation/ground-state.md). That page explains the current two-pass reconstruction, stopping controls, and `pred_gs_deg` behavior. Convergence of this inner eigenproblem must be assessed separately from convergence with respect to the NOOC space and the outer DMFT loop.

Lanczos also gives projected resolvents for Green's functions. Their starting vectors are particle-addition or particle-removal states, rather than the random starting vector used for the ground-state search. See [Green's functions](../implementation/greens-functions.md).

## References

- A. Weiße et al. (2006), “The kernel polynomial method,” *Reviews of Modern Physics* **78**, 275–306.
- K. Wu and H. Simon (2000), “Thick-restart Lanczos method for large symmetric eigenvalue problems,” *SIAM Journal on Matrix Analysis and Applications* **22**, 602–616.
- G. H. Golub and C. F. Van Loan (2013), *Matrix Computations*, JHU Press.
- J.-M. Wang et al., *iNORG: An open-source quantum impurity solver package based on the natural orbitals renormalization group*, “Lanczos algorithm.”
