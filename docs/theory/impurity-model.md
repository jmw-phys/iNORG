# Anderson impurity model

An Anderson impurity model separates an interacting local system from a noninteracting environment. The local system may contain several correlated orbitals. Electrons move between these orbitals and the bath through hybridization. In dynamical mean-field theory (DMFT), you solve this impurity problem repeatedly while adjusting the bath to represent the lattice environment.

## Hamiltonian and conventions

Write the impurity Hamiltonian as

$$
H_{\mathrm{imp}}=H_{\mathrm{loc}}+H_{\mathrm{bath}}+H_{\mathrm{hyb}},
$$

with

$$
\begin{aligned}
H_{\mathrm{loc}} & = \sum_{\alpha\beta}E_{\alpha\beta}
d_\alpha^\dagger d_\beta
+\sum_{\alpha\beta\gamma\delta}U_{\alpha\beta\gamma\delta}
d_\alpha^\dagger d_\beta^\dagger d_\gamma d_\delta,\\
H_{\mathrm{bath}} & = \sum_{\alpha k}\epsilon_{\alpha k}
c_{\alpha k}^\dagger c_{\alpha k},\\
H_{\mathrm{hyb}} & = \sum_{k\alpha\beta}V_k^{\alpha\beta}
d_\alpha^\dagger c_{\beta k}+\mathrm{h.c.}
\end{aligned}
$$

Here $d_\alpha^\dagger$ creates an impurity electron, $c_{\alpha k}^\dagger$ creates a bath electron, and $\alpha,\beta$ label spin-orbitals. The matrix $E$ contains local one-electron energies and mixing. The tensor $U$ contains the interaction coefficients in the operator ordering shown above. The bath energies $\epsilon$ and couplings $V$ specify the environment.

The tensor expression describes the general model. The current standalone inputs construct **density-density interactions**. In mode 1, use `CoulombF 'Ising'`; changing this token does not enable spin-flip or pair-hopping terms. The conversion from `U` and `J` to the actual interaction coefficients is explained in the [one-shot tutorial](../tutorials/siam-oneshot.md).

<!-- iNORG paper: Model and methods / Anderson impurity model, inorg.tex:136, label eq1:aim. Current interaction support: src/api_edmft.cpp:138; src/impurity.cpp:184; src/impurity.cpp:394. -->

## Replacing the continuous bath

For the diagonal hybridization used by the tutorials, a finite bath gives

$$
\Delta_\alpha(z)=\sum_{k=1}^{N_b}
\frac{|V_{\alpha k}|^2}{z-\epsilon_{\alpha k}}.
$$

The paper denotes this function by $\Gamma_{\mathrm{imp},\alpha}$. This manual uses $\Delta$ to match `Delta.inp`. Fitting adjusts the bath energies and couplings so that this rational function approximates the requested hybridization. Increasing the bath size improves the available frequency resolution, but also enlarges the many-body problem. See [hybridization fitting](../implementation/hybridization-fitting.md) for the objective and its diagnostics.

<!-- iNORG paper: Hybridization function fitting algorithm, label eq:hyb_imp. src/impurity.cpp:76; src/hyberr.cpp:249. -->

## What NORG approximates

Even a noninteracting bath contributes many fermionic degrees of freedom. Full exact diagonalization therefore becomes expensive as the bath grows. NORG exploits the fact that the interactions are confined to the impurity: after optimizing the bath orbitals, many of them are almost full or almost empty. Their occupation fluctuations can be restricted while retaining more freedom in the active orbitals.

There are two separate approximations to inspect. Bath fitting replaces a continuous hybridization by a finite bath. The [natural-orbital constraints](../implementation/nooc.md) then select a subspace of that finite model. A good bath fit alone does not establish convergence of the many-body solution.

The current solver computes ground-state properties. Matsubara frequencies provide a useful sampling grid for Green's functions; setting `beta` does not turn the calculation into a finite-temperature ensemble.

<!-- iNORG paper: Program summary; Model and methods / Natural orbitals renormalization group method. src/densitymatrix.cpp:100; src/prmtr.cpp:249. -->

## References

- P. W. Anderson (1961), “Localized Magnetic States in Metals,” *Physical Review* **124**, 41–53.
- A. Georges et al. (1996), “Dynamical mean-field theory of strongly correlated fermion systems and the limit of infinite dimensions,” *Reviews of Modern Physics* **68**, 13–125.
- J.-M. Wang et al., *iNORG: An open-source quantum impurity solver package based on the natural orbitals renormalization group*, “Anderson impurity model.”
