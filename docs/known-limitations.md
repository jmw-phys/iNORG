# Known limitations

These notes describe the current source and supplied examples. Check the input and iteration history before treating a run as a physical result.

## Bath seed row length is not checked

When you supply an `ose_hop` file, every `the` and `epi` row must contain exactly as many values as the matching `fit_nbaths` entry. The reader uses the configured array sizes, reads values without checking stream failure afterward, and does not reject a wrong count. A run can continue with a partly misread seed. A missing `ose_hop` produces only a warning; the generated starting bath remains available and mode 1 still fits it.

<!-- src/bath.cpp:25; src/bath.cpp:29; src/bath.cpp:361; src/bath.cpp:364; src/bath.cpp:374; src/api_edmft.cpp:32 -->

## Mode-0 settings that do not take effect

Mode 0 does not parse `beta` or `NOOC`. Its frequency mesh retains $\beta=50\pi$, and its constraint retains `phss_v2`. Adding those lines to `PARAMS.norg` does not change the calculation. Other unsupported keywords are also silently ignored.

With the mode-0 default `bath_ph_symm 1`, a nonzero `bethe_mu` is an error. Setting `bath_ph_symm 0` allows it through that check, but subsequent bath setup resets the chemical potential to zero before the DMFT driver reads it. A nonzero `bethe_mu` therefore does not provide a working doping control in the current standalone path.

<!-- src/main.cpp:168; src/main.cpp:240; src/prmtr.cpp:33; src/prmtr.cpp:82; src/main.cpp:245; src/main.cpp:283; src/main.cpp:289; src/prmtr.cpp:129; src/dmft_bethe.cpp:17 -->

## Double-occupancy output

The $\langle n_i n_j\rangle$ block in `nmat` output is zero in the current bath-only rotation path. The operator table contains one-body off-diagonal entries, while the double-occupancy routine only accumulates entries with four-fermion indices. Those entries are absent in this path, so its initially zero accumulator stays zero.

Do not interpret this as vanishing physical double occupancy. The printed spin correlations, $\langle S_z^2\rangle$, $\langle S^2\rangle$, and charge correlations derived from that matrix are also unreliable. The first block of one-particle occupations is calculated separately.

<!-- src/prmtr.cpp:69; src/operator.cpp:11; src/operator.cpp:199; src/operator.cpp:203; src/operator.cpp:215; src/densitymatrix.cpp:204; src/densitymatrix.cpp:213; src/densitymatrix.cpp:231; src/norg.cpp:752; src/norg.cpp:791; src/norg.cpp:803; src/norg.cpp:811; src/norg.cpp:816; src/norg.cpp:831 -->

## Bethe DMFT convergence

At half filling, the Bethe problem is particle-hole (PH) symmetric. An odd symmetric bath has a level at zero energy. In an insulating solution that level decouples. With `bath_ph_symm 0`, the unconstrained fit can move it away from zero, break PH symmetry, change the impurity ground state, and prevent the DMFT loop from converging.

Mode 0 defaults to `bath_ph_symm 1`, including the supplied tutorial, which omits this key. The fit uses mirrored hoppings and opposite-energy pairs with the middle energy fixed at zero. The driver sets $\operatorname{Re}\Sigma_{aa}(i\omega_n)=0$ for the initial self-energy before storing it and for each computed impurity self-energy before writing it, evaluating the error, or updating the next input. This keeps the stored iterates and next input projected, including Pulay residuals when mixing is enabled. It prevents symmetry drift in the fit and update. The input checks require `bethe_mu 0` and odd `fit_nbaths` entries. Mode 1 and the eDMFT path retain the unconstrained fit. See the [fit parameterization](implementation/hybridization-fitting.md#ph-symmetric-fit-in-mode-0).

<!-- src/main.cpp:39; src/main.cpp:44; src/main.cpp:154; src/main.cpp:280; src/main.cpp:283; src/main.cpp:286; src/prmtr.h:45; src/hyberr.h:123; src/hyberr.h:136; src/hyberr.h:149; src/bath.cpp:91; src/bath.cpp:215; src/dmft_bethe.cpp:53; src/dmft_bethe.cpp:74; tutor/2_Multi-orbital_model_on_Bethe_lattice_with_DMFT/PARAMS.norg:1 -->

An existing `ose_hop` bypasses the first fit, so supply a PH-symmetric seed when this constraint is enabled. Inspect the full iteration history to assess convergence. The tutorial page shows a [reference convergence curve](tutorials/bethe-dmft.md#convergence) and an approximate run time.

<!-- src/dmft_bethe.cpp:47; src/dmft_bethe.cpp:60; src/dmft_bethe.cpp:63 -->

The mode-0 driver defaults to `pulay_mixing 0`, which uses plain substitution (`se = seimp`) at every iteration. Pulay mixing is used only with `pulay_mixing 1`: it starts once an iteration after the first reaches a self-energy error below $8\times10^{-3}$ and remains active thereafter. Neither choice guarantees convergence. The stopping function also includes tests on the recent error history, so termination alone does not establish physical self-consistency. The inner message `THE NORG CONVERGED` is printed even if that inner loop reaches its iteration limit, and does not certify convergence of the outer DMFT loop.

<!-- src/main.cpp:42; src/main.cpp:155; src/main.cpp:281; src/dmft_bethe.cpp:84; src/dmft_bethe.cpp:88; src/dmft_bethe.cpp:91; src/dmft_bethe.cpp:394; src/norg.cpp:68; src/norg.cpp:103 -->

Start from a better matched bath or consistent restart data from a nearby converged model. Compare the complete self-energy and occupation histories, not only the last error. The `pulay_mixing` key selects plain substitution or permits Pulay mixing; it does not set a damping strength. Changing the damping requires a source change and rebuild. The Pulay implementation currently uses unit residual weights.

<!-- src/dmft_bethe.cpp:33; src/dmft_bethe.cpp:47; src/dmft_bethe.cpp:48; src/dmft_bethe.cpp:308; src/dmft_bethe.cpp:329; src/main.cpp:183; src/main.cpp:187 -->

## Frequency-mesh cost in mode 0

The hybridization update calls `g_loc.inverse()` inside the frequency loop. That method inverts the matrices at every frequency each time it is called. For $N$ frequencies and a fixed band count, this part therefore performs $O(N^2)$ frequency work. It adds overhead to iterations that construct a new bath target; it is not a measured total-runtime scaling for the solver.

<!-- src/dmft_bethe.cpp:60; src/dmft_bethe.cpp:61; src/dmft_bethe.cpp:207; src/dmft_bethe.cpp:211; src/green.h:152; src/green.h:155 -->

## Mode-1 interaction and input restrictions

Use `CoulombF 'Ising'`. Although the supplied input comment mentions `'Full'`, changing the token only triggers a warning and does not enable rotationally invariant interactions. The interaction builder accepts 2, 3, or 5 spatial orbitals and constructs density-density terms.

<!-- src/api_edmft.cpp:138; src/impurity.cpp:184; src/impurity.cpp:394; src/impurity.cpp:505 -->

Several mode-1 inputs have no usable default, and list lengths are not comprehensively checked. In particular, supply `Ed`, `deg_idx` (or `Deg`), `U`, `J`, `beta`, `fit_nbaths`, `restrain`, `weight_nooc`, and `weight_freze`. The frequency values in `Delta.inp` are consumed but not checked against the internal mesh. Use the dimensions and units in the [parameter reference](parameters.md) and [one-shot tutorial](tutorials/siam-oneshot.md).

<!-- src/api_edmft.cpp:123; src/api_edmft.cpp:149; src/api_edmft.cpp:164; src/api_edmft.cpp:452; src/green.cpp:222; src/green.cpp:240 -->
