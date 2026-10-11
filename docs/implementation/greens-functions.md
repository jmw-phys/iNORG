# Green's functions

After obtaining the ground state, iNORG evaluates the response to adding or removing one electron. Lanczos projects these neighboring particle-number problems into Krylov spaces. Diagonal Green's functions can then be evaluated as continued fractions; matrix-valued functions also require projections onto distinct source and sink states.

## Zero-temperature resolvent

For ground state $|\psi_0\rangle$ and energy $E_0$, the retarded Green's function is

$$
\begin{aligned}
G^R_{ij\sigma}(z)={}&
\langle\psi_0|c_{i\sigma}
(z+E_0-H)^{-1}c_{j\sigma}^\dagger|\psi_0\rangle\\
&+\langle\psi_0|c_{j\sigma}^\dagger
(z-E_0+H)^{-1}c_{i\sigma}|\psi_0\rangle,
\qquad z=\omega+i\eta,\quad\eta>0.
\end{aligned}
$$

The first term is the particle-addition contribution $G^>$; the second is the particle-removal contribution $G^<$. Here these symbols label the two resolvent contributions, following the paper. Addition propagates in a sector with one extra particle, and removal in a sector with one fewer. The ground-state sector alone is insufficient.

The same resolvent expressions can be evaluated at $z=i\omega_n$. The supplied standalone paths write imaginary-frequency results. They do not request real-axis output through a parameter-file switch. In this context, `beta` sets the Matsubara sampling mesh for the ground-state Green's function.

<!-- iNORG paper: Continued-fraction evaluation of diagonal Green's functions, inorg.tex:721. src/norg.cpp:276; src/api_edmft.cpp:58; src/api_edmft.cpp:64; src/dmft_bethe.cpp:72. -->

## Continued fraction for a diagonal element

For a selected spin-orbital, define the source states and their weights:

$$
|f^-\rangle=c|\psi_0\rangle,\qquad
|f^+\rangle=c^\dagger|\psi_0\rangle,
\qquad W_\pm=\langle f^\pm|f^\pm\rangle.
$$

Normalize a nonzero source to start its Lanczos run. The tridiagonal matrices $T_-$ and $T_+$ are built separately in the removal and addition sectors. Then

$$
G^<(z)=W_-[(z+T_--E_0)^{-1}]_{00},
\qquad
G^>(z)=W_+[(z-T_++E_0)^{-1}]_{00}.
$$

The weights are essential. Normalizing a source changes its norm; it must not change the physical spectral weight. A vanishing source contributes zero and should not be normalized in the mathematical construction.

<!-- Paper wording: the normalized-source derivation writes a unit numerator, omitting the original source norm. The active code restores it through upper_fraction: src/crrltfun.cpp:69; src/crrltfun.cpp:77; src/crrltfun.cpp:115. W notation makes this factor explicit. -->

With diagonal coefficients $\alpha_i^\pm$ and off-diagonal coefficients $\beta_i^\pm$, the diagonal inverse becomes

$$
G^<(z)=\cfrac{W_-}{z-E_0+\alpha_0^-
-\cfrac{(\beta_1^-)^2}{z-E_0+\alpha_1^-
-\cfrac{(\beta_2^-)^2}{z-E_0+\alpha_2^- -\cdots}}},
$$

$$
G^>(z)=\cfrac{W_+}{z+E_0-\alpha_0^+
-\cfrac{(\beta_1^+)^2}{z+E_0-\alpha_1^+
-\cfrac{(\beta_2^+)^2}{z+E_0-\alpha_2^+ -\cdots}}}.
$$

Terminate at a finite Krylov dimension and evaluate from the deepest denominator upward. Once the coefficients are available, this scalar evaluation avoids a large matrix inverse at every frequency. Add the two contributions to obtain the diagonal Green's function.

<!-- iNORG paper: Continued-fraction evaluation of diagonal Green's functions, continued-fraction expressions. src/crrltfun.cpp:110; src/crrltfun.cpp:113; src/crrltfun.cpp:114. -->

### Current diagonal implementation

`CrrltFun::find_gf` uses the same evaluation routine with opposite energy signs for addition and removal. It starts with a compiled number of iterations depending on whether the mesh is imaginary or real, then extends the recurrence and checks coefficient or Green's-function changes. These controls are not exposed as `PARAMS.norg` keys.

The current recurrence retains Krylov vectors for its orthogonalization step. Do not interpret the small storage needed for the scalar continued-fraction coefficients as the memory use of the whole Green's-function solve. Its Hamiltonian, source state, and retained vectors still live in the selected many-body space.

<!-- Paper wording: O(N) memory for continued-fraction evaluation. Full current routine: src/crrltfun.cpp:79; src/crrltfun.cpp:87; src/crrltfun.cpp:101; src/crrltfun.cpp:120; src/crrltfun.cpp:476. -->

## Matrix-valued Green's functions

A diagonal continued fraction evaluates a resolvent between the same source and sink. A generic element needs distinct vectors,

$$
G_{\alpha\beta}^>(\omega)=
\langle\phi'_\alpha|(z-H)^{-1}|\phi_\beta\rangle,
\qquad z=E_0+\omega+i\eta.
$$

In this subsection, $z$ includes the ground-state energy, as in the paper's matrix formulation. For each source column $\beta$, build a Lanczos basis $V$ and projected Hamiltonian $T=V^\dagger HV$. The approximation is

$$
G_{\alpha\beta}^>\approx
\langle\phi'_\alpha|V(z-T)^{-1}V^\dagger|\phi_\beta\rangle.
$$

Diagonalize the small matrix with $S^\dagger TS=D$. Define

$$
|m_\beta\rangle=S^\dagger V^\dagger|\phi_\beta\rangle,
\qquad
\langle n_\alpha|=\langle\phi'_\alpha|VS.
$$

Each frequency then requires the small sum

$$
G_{\alpha\beta}^>(\omega)\approx
\sum_k\frac{(n_\alpha)_k(m_\beta)_k}{z-D_k}.
$$

Here $(n_\alpha)_k$ is already a component of a bra; its conjugation is included in the definition. Keep the original source and sink amplitudes in these projections. The removal part follows with the appropriate removal-sector Hamiltonian and energy signs.

<!-- iNORG paper: Krylov-subspace evaluation of matrix-valued Green's functions, inorg.tex:784. -->

One source-resolved Lanczos run supplies one projected source and all projected sinks, hence one matrix column. An $N_{\mathrm{orb}}\times N_{\mathrm{orb}}$ matrix requires $N_{\mathrm{orb}}$ such runs for each particle/removal part. The frequency grid reuses those projections. For the pre-rotated form above, the paper gives assembly cost $O(N_KN_{\mathrm{orb}}^2N_\omega)$ after projection and diagonalization. That cost excludes constructing the Krylov spaces.

### Comparison with correction vectors

The alternative correction-vector method solves

$$
(\omega+i\eta+E_0-H)|\chi(\omega)\rangle
=A|\psi_0\rangle
$$

for every source and frequency, then contracts with the desired sink. Each solve acts in the full selected Hilbert space. The source-resolved Krylov approach replaces repeated large linear solves with reusable projected data. The total advantage depends on projection cost, the number of sources, the number of frequencies, and the accuracy required.

<!-- iNORG paper: Krylov-subspace evaluation of matrix-valued Green's functions, correction-vector equation and assembly-cost discussion. -->

### Scope of the current code

The matrix routine is `Crrvec::krylov_space_for_green_matrix`. It builds a separate recurrence for each source and gathers sink overlaps. Its implementation currently retains the full source vectors, repeats small diagonalizations for source/sink pairs, and applies small eigenvector matrices inside the frequency loop. Thus the fully pre-rotated assembly formula above describes the paper's algorithmic organization, not an exact operation or memory count for this routine.

The standalone tutorials call the diagonal continued-fraction path. Describing the matrix algorithm does not imply that an off-diagonal `Delta.inp` layout or a new mode selector is accepted by those inputs.

<!-- Paper wording: retain only projected source/sink quantities and pre-rotate before frequency assembly. Current code: src/crrvec.cpp:257; src/crrvec.cpp:267; src/crrvec.cpp:294; src/crrvec.cpp:298; src/crrvec.cpp:306; src/crrvec.cpp:311. Standalone dispatch: src/api_edmft.cpp:65; src/api_edmft.cpp:66; src/norg.cpp:300; src/dmft_bethe.cpp:72. -->

## Output and checks

For each calculation, compare the target and fitted hybridization before assessing the resulting $G$ and $\Sigma$. Check smoothness on the imaginary axis and stability as you enlarge the bath or relax constraints. The code forms the self-energy through $\Sigma=G_0^{-1}-G^{-1}$, so both noninteracting and interacting Green's functions enter the result.

See [output files](../output-files.md) for mode-dependent names and column ordering. The tutorial [reference outputs](../tutorials/siam-oneshot.md#reference-output) give results to compare with.

<!-- src/api_edmft.cpp:70; src/dmft_bethe.cpp:73. -->

## References

- E. R. Gagliano and C. A. Balseiro (1987), “Dynamical Properties of Quantum Many-Body Systems at Zero Temperature,” *Physical Review Letters* **59**, 2999–3002.
- R. Haydock, V. Heine, and M. J. Kelly (1975), “Electronic structure based on the local atomic environment for tight-binding bands. II,” *Journal of Physics C: Solid State Physics* **8**, 2591.
- T. D. Kühner and S. R. White (1999), “Dynamical correlation functions using the density matrix renormalization group,” *Physical Review B* **60**, 335–343.
- R. Peters (2011), “Spectral functions for single- and multi-impurity models using density matrix renormalization group,” *Physical Review B* **84**, 075139.
- J.-M. Wang et al., *iNORG: An open-source quantum impurity solver package based on the natural orbitals renormalization group*, “Continued-fraction evaluation of diagonal Green's functions” and “Krylov-subspace evaluation of matrix-valued Green's functions.”
