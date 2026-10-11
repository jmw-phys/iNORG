# Hybridization function fitting

Bath fitting turns the hybridization supplied by a DMFT calculation into bath energies and impurity-bath couplings. You choose the bath size and, in mode 1, the frequency rows used by the fit. The fitter adjusts the parameters of the finite bath.

## Discrete bath and least-squares objective

For a diagonal hybridization channel,

$$
\Delta_{\mathrm{fit}}(z)=\sum_{k=1}^{N_b}
\frac{|V_k|^2}{z-\epsilon_k}.
$$

In matrix notation this is $V(zI-\epsilon_{\mathrm{bath}})^{-1}V^\dagger$. The fitted parameters $\boldsymbol\theta$ contain the bath energies and couplings. The residual vector contains real and imaginary parts of the difference between the target and fitted functions, together with any additional penalties.

<!-- iNORG paper: Hybridization function fitting algorithm / Fitting details, inorg.tex:359. The paper's residual displays V(epsilon_bath-z)^(-1)V^dagger, opposite to its scalar hybridization formula. Use the scalar formula and code sign: src/hyberr.cpp:249; src/impurity.cpp:88. -->

The paper describes an objective with three contributions:

$$
\begin{aligned}
\chi^2={}&\sum_n w_n
\left|\Delta_{\mathrm{target}}(z_n)-\Delta_{\mathrm{fit}}(z_n)\right|^2\\
&+w_{\mathrm{mom}}
\left|VV^\dagger-\left(K_{0,\mathrm{loc}}-H_{0,\mathrm{loc}}^2\right)\right|^2\\
&+w_{\mathrm{reg}}\left(\sum_k\epsilon_k^4\right)^2.
\end{aligned}
$$

For a matrix target, the absolute-value terms represent sums over the fitted matrix components. The weights set the relative importance of fitting the frequency curve, reproducing its high-frequency moment, and keeping bath energies bounded. The local moments are

$$
H_{0,\mathrm{loc}}=\int_{\mathrm{BZ}}\frac{d^dk}{(2\pi)^d}H_0(k),
\qquad
K_{0,\mathrm{loc}}=\int_{\mathrm{BZ}}\frac{d^dk}{(2\pi)^d}H_0(k)^2.
$$

The leading high-frequency coefficient of the hybridization is $VV^\dagger$. Matching it to the local second central moment imposes the bath sum rule discussed by Koch, Sangiovanni, and Gunnarsson.

**The active fitter uses a different penalty implementation.** Its explicit moment penalty is commented out. It instead adds the real and imaginary hybridization at row `nmesh - 2` as a high-frequency residual. Each bath energy itself is a regularization residual with target zero, so the resulting energy penalty is quadratic after squaring the residuals. The formula above explains the paper's design; the table below describes the objective you actually run.

| Active contribution | Construction | Output diagnostic |
| --- | --- | --- |
| Frequency curve | Selected real and imaginary values; components are normalized using the target magnitudes and weights. | `err_crv` |
| High-frequency point | A separately weighted point near the end of the full input mesh. | `err_hight_ev` |
| Bath-energy regularization | One weighted residual per bath energy. | `err_regE` |
| Combined objective | Euclidean norm of all weighted residuals, $\sqrt{\chi^2}$. | `err` |

<!-- Paper wording: explicit bath sum rule and (sum epsilon^4)^2 regularization. Code: src/hyberr.cpp:70 (disabled sum rule); src/hyberr.cpp:9; src/hyberr.cpp:86; src/hyberr.cpp:97; src/hyberr.cpp:130; src/hyberr.cpp:312; src/hyberr.h:65; src/impurity.cpp:210. -->

## PH-symmetric fit in mode 0

Mode 0 defaults to `bath_ph_symm 1`. This requires `bethe_mu 0` and odd entries in `fit_nbaths`; violations are input errors. For one impurity channel with $N_b=2h+1$ bath levels, the fit uses $N_b$ real parameters instead of $2N_b$:

$$
\mathbf r=(v_0,\ldots,v_h,e_0,\ldots,e_{h-1}).
$$

The $h+1$ hoppings and $h$ energy-pair parameters expand as

$$
\begin{aligned}
V_j=V_{2h-j}&=v_j,&
\epsilon_j&=-e_j,& \epsilon_{2h-j}&=e_j
&& (0\leq j<h),\\
V_h&=v_h,& \epsilon_h&=0.
\end{aligned}
$$

Thus the middle energy stays at zero while its hopping remains a fitted parameter. To initialize the reduced vector, the code sorts the full seed by energy, keeping each hopping with its energy. For these sorted seed values $\widetilde V_j,\widetilde\epsilon_j$, it takes

$$
v_j=\sqrt{\frac{\widetilde V_j^2+\widetilde V_{2h-j}^2}{2}},
\qquad
e_j=\frac{\widetilde\epsilon_{2h-j}-\widetilde\epsilon_j}{2}
\quad (0\leq j<h),
\qquad v_h=\widetilde V_h.
$$

`ph_reduce` and `ph_expand` implement these maps. `HybErrPH` evaluates the same objective on the expanded bath and combines its analytic derivatives by the chain rule. MPI trials optimize the reduced vector, then expand the best result for the impurity solve.

<!-- src/main.cpp:154; src/main.cpp:176; src/main.cpp:283; src/main.cpp:286; src/hyberr.h:123; src/hyberr.h:136; src/hyberr.h:149; src/bath.cpp:91; src/bath.cpp:149; src/bath.cpp:168; src/bath.cpp:215; src/bath.cpp:236 -->

At half filling the Bethe problem is particle-hole (PH) symmetric. In an insulating solution the zero-energy bath level decouples. An unconstrained fit can move it away from zero, break PH symmetry, change the impurity ground state, and stop DMFT convergence. Fixing its energy prevents this drift. The driver also sets $\operatorname{Re}\Sigma_{aa}(i\omega_n)=0$ for the initial self-energy before storing it and for each computed impurity self-energy before writing it, evaluating the error, or updating the next input. The stored iterates and next self-energy remain projected.

Mode 0 defaults to plain substitution (`se = seimp`) with `pulay_mixing 0`. Pulay mixing is used only with `pulay_mixing 1`, starting once an iteration after the first has self-energy error below `8e-3` and remaining active thereafter. The residuals and mixed self-energy also stay projected when PH symmetry is enabled. Mode 1 and the eDMFT path do not use this loop setting.

<!-- src/main.cpp:155; src/main.cpp:281; src/dmft_bethe.cpp:53; src/dmft_bethe.cpp:74; src/dmft_bethe.cpp:88; src/dmft_bethe.cpp:91; src/dmft_bethe.cpp:96 -->

Setting `bath_ph_symm 0` restores the unconstrained fit and omits the self-energy projection. Mode 1 and the eDMFT path do not set this flag and retain the unconstrained fit. An existing `ose_hop` bypasses the first mode-0 fit, so it should already be PH symmetric when the constraint is enabled.

<!-- src/main.cpp:39; src/main.cpp:44; src/prmtr.h:45; src/bath.cpp:221; src/dmft_bethe.cpp:47; src/dmft_bethe.cpp:60; src/dmft_bethe.cpp:53; src/dmft_bethe.cpp:74; src/dmft_bethe.cpp:96 -->

## Levenberg–Marquardt updates

The paper presents the damped least-squares step

$$
\boldsymbol\theta_{\mathrm{new}}=
\boldsymbol\theta_{\mathrm{old}}-
(J^\mathsf TJ+\lambda I)^{-1}J^\mathsf T\mathbf r,
$$

where $J$ is the Jacobian of the residual vector. Large damping makes the step more conservative. Reducing the damping near a minimum allows faster local convergence. The code uses analytic derivatives of the hybridization with respect to bath energies and couplings.

The implemented variant multiplies the diagonal of the normal-equation matrix by $1+\lambda$, rather than adding $\lambda I$. It lowers damping after an improvement and raises it after a rejected step. It checks changes in the objective, an iteration limit, and singularity of the linear solve. The paper's discussion of Jacobian condition-number monitoring should not be read as an exposed diagnostic in the current fitter.

<!-- Paper wording: J^T J + lambda I and Jacobian condition-number monitoring. Code: src/gen/fitmrq.h:74; src/gen/fitmrq.h:134; src/gen/fitmrq.h:138; src/gen/fitmrq.h:154; src/gen/fitmrq.h:155; src/hyberr.cpp:275. -->

The objective can have multiple local minima. The bath module tries perturbed starting parameters, distributes trials among MPI workers, sorts the candidates, and refines promising fits. The best parameters are broadcast for the impurity solve. This master-worker stage is why the executable requires at least two MPI ranks.

<!-- src/bath.cpp:138; src/bath.cpp:184; src/bath.cpp:204; src/bath.cpp:212; src/bath.cpp:232; src/bath.cpp:240; src/main.cpp:21. -->

## Inputs and fit assessment

Set `fit_nbaths` to an odd bath count for every spin-orbital set, with matching counts for equivalent channels. Mode 1 accepts `fit_points` as **zero-based row indices**, not frequency values. If you supply them, the fit samples those rows while the Green's-function reader still consumes the full mesh. Without them, `fit_range` selects a contiguous low-frequency region. See [parameters](../parameters.md) for the exact bounds and defaults.

The one-shot example's sparse row selection is motivated by the discrete Lehmann representation. iNORG reads the supplied indices; this input path does not generate a discrete Lehmann grid. An existing `ose_hop` seeds the fit, but its row lengths must match the requested bath sizes. The example ships without `ose_hop`, so the solver warns and fits from its generated initial bath; see the [quick start](../quickstart.md).

<!-- iNORG paper: Example 1, fit_points description. src/api_edmft.cpp:521; src/prmtr.cpp:241; src/hyberr.cpp:42; src/green.cpp:236; src/bath.cpp:361. -->

Compare target and fitted hybridizations over the full mesh, especially the low-energy region relevant to your calculation. Also inspect the separate objective contributions in `h0.txt` or its iteration-prefixed equivalent. A small combined objective depends on the selected rows and weights; it is not an independent bound on the final self-energy error.

The paper's fitting-performance example compares wide and narrow bands with $U=3$, $\Delta=0.3$, and $t_2=0.5t_1$. Here $\Delta$ is the model parameter in the figure caption, not this page's hybridization symbol. As reported in the paper, the fit errors decrease exponentially with bath count in that example and fall below $10^{-4}$ for $n_b\geq7$. This describes that fit comparison, not a bath-size guarantee for arbitrary input.

<!-- iNORG paper: Fitting performance, inorg.tex:393, label fig:bathfit. No numerical results have been generated for this manual. -->

## References

- E. Koch, G. Sangiovanni, and O. Gunnarsson (2008), “Sum rules and bath parametrization for quantum cluster theories,” *Physical Review B* **78**, 115102.
- X. Wang et al. (2023), “Efficient and robust Levenberg–Marquardt Algorithm based on damping parameters for parameter inversion in underground metal target detection,” *Computers & Geosciences* **176**, 105354.
- J. Kaye, K. Chen, and O. Parcollet (2022), “Discrete Lehmann representation of imaginary time Green's functions,” *Physical Review B* **105**, 235115.
- J.-M. Wang et al., *iNORG: An open-source quantum impurity solver package based on the natural orbitals renormalization group*, “Hybridization function fitting algorithm.”
