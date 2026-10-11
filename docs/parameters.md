# PARAMS.norg reference

The solver reads `PARAMS.norg` from its working directory. Keywords are case-sensitive. Write one keyword and its value per line, with whitespace between them; do not use `=`. Lists use square brackets, with spaces or commas between entries. Keep each list on one line and close it with `]` before an optional `#` comment.

```text
mode 1
U 8.0
deg_idx [1, 2, 3]
```

Mode selection strips `#` comments. The two parameter parsers read the first token and then the expected scalar or list; neither reports unknown keys. A misspelled or unsupported key can therefore have no effect. Avoid duplicate scalar keys and duplicate lists. Lists are appended during parsing, including when both aliases of a key are present.

<!-- src/main.cpp:83; src/main.cpp:163; src/api_edmft.cpp:452; src/api_edmft.cpp:463; src/api_edmft.cpp:474 -->

Defaults below describe the effective input path, not just constructor values. **Required** means there is no usable fallback: some missing entries are rejected, while others leave empty arrays or uninitialized scalars. The mode-1 parser does not provide a complete input validation step. Supply every required field.

## Keys shared by both modes

| Keyword | Type | Default | Meaning | Used in mode |
| --- | --- | --- | --- | --- |
| `mode` | Integer selector | required for a standalone tutorial | `0`: Bethe DMFT; `1`: one-shot impurity. A missing key or any other value selects the eDMFT path, not a default tutorial. | 0, 1 <!-- src/main.cpp:102; src/main.cpp:106; src/main.cpp:125; src/main.cpp:168 --> |
| `fit_nbaths` | List of integers | required; parser starts with an empty list | Bath count for each spin-orbital set, in interleaved order `1up, 1dn, 2up, 2dn, ...`. Supply `2 * bethe_band` entries in mode 0 or twice the length of `deg_idx` in mode 1. See bath-count constraints below. | 0, 1 <!-- src/main.cpp:159; src/main.cpp:204; src/main.cpp:277; src/api_edmft.cpp:132; src/api_edmft.cpp:543; src/prmtr.cpp:133 --> |
| `restrain` | List of integers | required; parser starts with an empty list | Limits used to constrain occupations across natural-orbital divisions. Use the eight-entry form shown in the tutorials for the current bath-only rotation setup. | 0, 1 <!-- src/main.cpp:160; src/main.cpp:214; src/main.cpp:269; src/api_edmft.cpp:131; src/api_edmft.cpp:554; src/prmtr.cpp:131 --> |
| `pred_gs_deg` | Integer | `0` | Controls the density-matrix eigensolve branch. Exactly `1` enables its explicit multi-state branch. Mode 0 also initializes the stored degeneracy to this value; mode 1 does not. See the degeneracy note below. | 0, 1 <!-- src/main.cpp:153; src/main.cpp:174; src/main.cpp:248; src/api_edmft.cpp:498; src/prmtr.cpp:68; src/api_edmft.cpp:291; src/densitymatrix.cpp:102 --> |

## Mode 0: Bethe lattice DMFT

Use these keys together with the shared keys.

| Keyword | Type | Default | Meaning | Used in mode |
| --- | --- | --- | --- | --- |
| `bath_ph_symm` | Exact token `0` or `1` | `1` in mode 0 | `1`: fit a particle-hole (PH) symmetric bath with mirrored hoppings, paired energies, and a middle level fixed at zero; project the initial self-energy before storing it and each computed impurity self-energy before output, error evaluation, and mixing by setting their real diagonal parts to zero. Requires `bethe_mu 0` and odd `fit_nbaths` entries; either violation is an error. `0`: use the unconstrained fit and omit this self-energy projection. | 0 <!-- src/main.cpp:154; src/main.cpp:176; src/main.cpp:181; src/main.cpp:280; src/main.cpp:283; src/main.cpp:286; src/hyberr.h:123; src/bath.cpp:215; src/dmft_bethe.cpp:53; src/dmft_bethe.cpp:74 --> |
| `pulay_mixing` | Exact token `0` or `1` | `0` in mode 0 | `0`: use plain substitution (`se = seimp`) at every iteration. `1`: start with plain substitution, then retain Pulay mixing once an iteration after the first has self-energy error below `8e-3`. Only trailing whitespace or a `#` comment may follow the token; any other value or extra token is an error. | 0 <!-- src/main.cpp:155; src/main.cpp:183; src/main.cpp:185; src/main.cpp:187; src/main.cpp:188; src/main.cpp:281; src/dmft_bethe.cpp:55; src/dmft_bethe.cpp:88; src/dmft_bethe.cpp:91 --> |
| `bethe_mu` | Real | `0.0` | Chemical-potential offset. A nonzero value is rejected when `bath_ph_symm 1`. With `bath_ph_symm 0`, it is parsed but reset to zero by bath setup before DMFT uses it. It cannot currently control doping. | 0 <!-- src/main.cpp:151; src/main.cpp:170; src/main.cpp:245; src/main.cpp:283; src/main.cpp:289; src/prmtr.cpp:129; src/dmft_bethe.cpp:17 --> |
| `bethe_band` | Integer | `0` (invalid; supply a positive value) | Number of spatial orbitals per lattice site. There are twice as many spin-orbital sets. | 0 <!-- src/main.cpp:152; src/main.cpp:172; src/main.cpp:246; src/main.cpp:275 --> |
| `bethe_t` | List of reals | required; parser starts with an empty list | One hopping scale $t_a$ per band. The semicircular noninteracting band has half-bandwidth $2\lvert t_a\rvert$. Use nonzero hopping scales. | 0 <!-- src/main.cpp:156; src/main.cpp:190; src/main.cpp:276; src/dmft_bethe.cpp:256 --> |
| `bethe_u` | Real | `0.0` | Intra-orbital density-density interaction $U$. | 0 <!-- src/main.cpp:157; src/main.cpp:200; src/main.cpp:252; src/impurity.cpp:270 --> |
| `bethe_uprim` | Real | `0.0` | Inter-orbital density-density interaction $U'$. Mode 0 has zero Hund coupling in its current setup. | 0 <!-- src/main.cpp:158; src/main.cpp:202; src/main.cpp:253; src/dmft_bethe.cpp:21; src/prmtr.cpp:50; src/impurity.cpp:272 --> |
| `distribute` | List of integers | required with a nonempty `restrain`; parser starts empty | Number of orbitals in each division, repeated for every spin-orbital set. Its length must match `restrain`; its sum must equal the number of sites in every set, including the impurity. | 0 <!-- src/main.cpp:161; src/main.cpp:227; src/main.cpp:278; src/main.cpp:292; src/main.cpp:306 --> |

The `bath_ph_symm` value must be exactly `0` or `1`, followed only by end of line, whitespace, or a `#` comment (which may immediately follow the value). Missing values, numeric prefixes such as `0.5` or `1junk`, `2`, and extra non-comment tokens report `bath_ph_symm must be 0 or 1`. When enabled, the input checks report `bath_ph_symm = 1 requires bethe_mu = 0` for nonzero chemical potential and `bath_ph_symm requires odd fit_nbaths entries` for any even bath count. See the [PH-symmetric fit](implementation/hybridization-fitting.md#ph-symmetric-fit-in-mode-0) for the parameterization. Mode 1 and the eDMFT path do not read this key and retain the unconstrained fit.

<!-- src/main.cpp:39; src/main.cpp:44; src/main.cpp:181; src/main.cpp:283; src/main.cpp:286; src/prmtr.h:45; src/bath.cpp:91; src/bath.cpp:221 -->

`pulay_mixing` controls only the mode-0 self-energy loop. Mode 1 and the eDMFT path do not use this loop setting. The parameter member is initialized to `1`, but the mode-0 reader replaces it with the input value or its default of `0`.

<!-- src/prmtr.h:46; src/main.cpp:39; src/main.cpp:42; src/main.cpp:44; src/main.cpp:155; src/main.cpp:281 -->

`NOOC`, `beta`, `nmesh`, `fit_points`, and mode-1 interaction keys are **not read in mode 0**. The `NOOC phss_v2` line in the supplied mode-0 tutorial agrees with the compiled default but does not change it. The imaginary-frequency mesh remains `nmesh = 8193`, with $\beta = 50\pi$ and $\omega_n=(2n+1)\pi/\beta$. There is no input keyword for DMFT mixing strength, convergence tolerance, or maximum iteration count. The compiled DMFT iteration limit is 200.

<!-- src/main.cpp:168; src/main.cpp:240; src/prmtr.cpp:26; src/prmtr.cpp:33; src/prmtr.cpp:35; src/prmtr.cpp:82; src/prmtr.cpp:249 -->

## Mode 1: one-shot impurity solution

Use these keys together with the shared keys. Aliases are listed separately for clarity; supply only one spelling of each alias pair.

| Keyword | Type | Default | Meaning | Used in mode |
| --- | --- | --- | --- | --- |
| `Ed` | List of reals | required | Impurity levels. For `magnetic 0`, each orbital uses `Ed[deg_idx - 1]`. With `magnetic 1`, supply one value per spatial orbital; the two spins share it. | 1 <!-- src/api_edmft.cpp:123; src/api_edmft.cpp:456; src/api_edmft.cpp:158; src/api_edmft.cpp:161 --> |
| `deg_idx` | List of integers | required | One 1-based degeneracy-group label per spatial orbital. Its length sets the number of bands. Use consecutive groups in nondecreasing order, starting at 1, so input and output column mapping agree. | 1 <!-- src/api_edmft.cpp:124; src/api_edmft.cpp:467; src/api_edmft.cpp:477; src/green.cpp:197; src/green.cpp:228 --> |
| `Deg` | List of integers | required unless `deg_idx` is supplied | Alias of `deg_idx`; do not supply both, since their entries would be appended to the same list. | 1 <!-- src/api_edmft.cpp:467; src/api_edmft.cpp:474 --> |
| `U` | Real | required | Coulomb input parameter, identified as $F^0$ in the tutorial. The interaction builder converts `U` and `J` to orbital-dependent density-density coefficients; `U` alone is not the final intra-orbital coefficient. | 1 <!-- src/api_edmft.cpp:128; src/api_edmft.cpp:495; src/api_edmft.cpp:141; src/impurity.cpp:399; src/impurity.cpp:468; tutor/1_Single-impurity_Anderson_model_with_multiple_orbitals/PARAMS.norg:4 --> |
| `J` | Real | required | Hund/Coulomb exchange input used with `U` to construct the interaction coefficients. | 1 <!-- src/api_edmft.cpp:125; src/api_edmft.cpp:480; src/impurity.cpp:469; src/impurity.cpp:486 --> |
| `CoulombF` | String token | required; use `'Ising'` | Only density-density (Ising) interactions are implemented on this path. Include the single quotes: the comparison is against the literal token `'Ising'`. Other tokens warn and still use the same interaction builder. | 1 <!-- src/api_edmft.cpp:126; src/api_edmft.cpp:483; src/api_edmft.cpp:138; src/impurity.cpp:184 --> |
| `beta` | Real | required | Sets the imaginary-frequency mesh through $\omega_n=(2n+1)\pi/\beta$. Use a positive value consistent with `Delta.inp`. This is a ground-state solver; this mesh parameter does not request a thermal ensemble. | 1 <!-- src/api_edmft.cpp:127; src/api_edmft.cpp:492; src/api_edmft.cpp:149; src/prmtr.cpp:249; src/densitymatrix.cpp:100 --> |
| `nmesh` | Integer | `8193` | Number of positive Matsubara rows read and written. Supply at least that many rows in `Delta.inp`. It must cover all `fit_points`; keep at least two rows because fitting also accesses the high-frequency tail. | 1 <!-- src/api_edmft.cpp:504; src/prmtr.cpp:35; src/prmtr.cpp:252; src/green.cpp:236; src/hyberr.cpp:97 --> |
| `fit_range` | Real | `5.0` | Sets `bandw = 4 * fit_range`. Without `fit_points`, the fit uses the first `round(fit_range * beta / (2*pi))` mesh points. The fitter also uses a high-frequency tail point. | 1 <!-- src/api_edmft.cpp:16; src/api_edmft.cpp:486; src/api_edmft.cpp:285; src/prmtr.cpp:241; src/hyberr.cpp:97 --> |
| `fit_points` | List of integer indices (read as reals, stored as integers) | `[]` | Zero-based Matsubara row indices used for fitting. They are not frequencies. Nonempty input replaces the contiguous fit-point selection. Use integers in `0 ... nmesh-1`; fractional values would be converted to integers. | 1 <!-- src/api_edmft.cpp:521; src/api_edmft.cpp:525; src/api_edmft.cpp:528; src/prmtr.h:47; src/prmtr.cpp:242; src/hyberr.cpp:42 --> |
| `NOOC` | String token | `phss_v2` | Selects the natural-orbital occupation constraint. Recognized names are `nooc`, `cpnooc`, `cnooc`, `phess`, and `phss_v2`; see below. | 1 <!-- src/api_edmft.cpp:489; src/prmtr.cpp:82; src/nocspace.cpp:297 --> |
| `Minsulator` | Integer | `0` | Alias of `pred_gs_deg` in mode 1. It sets the same eigensolve branch selector, not a separate insulating-state condition. | 1 <!-- src/api_edmft.cpp:498; src/prmtr.cpp:68; src/densitymatrix.cpp:102 --> |
| `magnetic` | Integer flag | `0` | `0`: use spin-degenerate hybridization groups. Nonzero: read and fit separate spin channels, with all up-spin bands followed by all down-spin bands in frequency files. See the qualification below. | 1 <!-- src/api_edmft.cpp:17; src/api_edmft.cpp:501; src/api_edmft.cpp:151; src/api_edmft.cpp:26; src/green.cpp:275 --> |
| `ful_pcl_sch` | Integer flag | `1` | Nonzero enables the particle-number-sector search before the final solve. Zero uses the initial/restarted sector and a preliminary solve to obtain natural-orbital weights. | 1 <!-- src/api_edmft.cpp:16; src/api_edmft.cpp:507; src/api_edmft.cpp:374; src/api_edmft.cpp:387 --> |
| `weight_nooc` | List of reals | required | Per-band threshold for assigning bath natural orbitals to constrained divisions. Supply at least one value per spatial orbital. The constructor's five `1e-4` entries are overwritten by the parsed list, even if empty. | 1 <!-- src/api_edmft.cpp:17; src/api_edmft.cpp:129; src/api_edmft.cpp:164; src/api_edmft.cpp:510; src/api_edmft.cpp:416 --> |
| `weight_freze` | List of reals | required | Per-band threshold for assigning bath natural orbitals to the outer divisions. Supply at least one value per spatial orbital. The constructor's five `1e-13` entries are overwritten by the parsed list, even if empty. Preserve this spelling. | 1 <!-- src/api_edmft.cpp:17; src/api_edmft.cpp:130; src/api_edmft.cpp:165; src/api_edmft.cpp:532; src/api_edmft.cpp:415 --> |

`distribute` is not read in mode 1. The solver derives division sizes from `fit_nbaths`, then updates them using natural-orbital weights. `bethe_*` keys have no effect here. The current interaction builder accepts 2, 3, or 5 spatial orbitals.

<!-- src/api_edmft.cpp:568; src/api_edmft.cpp:287; src/api_edmft.cpp:423; src/impurity.cpp:394 -->

## Bath counts and array sizes

Use positive **odd** bath counts in the supplied setup. For a requested count $B$, bath setup builds the eight-entry row

```text
[1, floor(B/2), 0, 0, 1, 0, 0, floor((B-1)/2)]
```

and derives the bath count from the sum minus the impurity entry. It equals $B$ for odd $B$; an even input becomes one extra bath site if it reaches this setup. Mode 0 with `bath_ph_symm 1` rejects even entries before bath setup. Keep the two spin entries for a band equal. In the spin-degenerate mode-1 fit, bands in the same degeneracy group also need matching bath sizes. The fitter takes the count from the up-spin entry and copies fitted parameters between equivalent bands.

<!-- src/main.cpp:284; src/main.cpp:289; src/prmtr.cpp:134; src/prmtr.cpp:135; src/prmtr.cpp:140; src/bath.cpp:21; src/bath.cpp:25; src/bath.cpp:45 -->

For mode 0, use equal counts across all sets when supplying a single `distribute` template. Each template row must sum to $B+1$, including its impurity orbital. The input checks enforce matching lengths and set sizes, but do not validate every physical or array constraint. Keep the eight-entry `restrain` and `distribute` structure of the example.

<!-- src/main.cpp:275; src/main.cpp:278; src/main.cpp:297; src/prmtr.cpp:135 -->

## NOOC, restrain, and distribute

Natural-orbital occupation constraints restrict which many-body configurations enter the solver. The first row of `control_divs` holds `restrain`. Each later row holds the number of orbitals in each division for one spin-orbital set. In mode 0, `distribute` supplies those later rows. Division sizes are counts of orbitals, not occupation probabilities.

<!-- src/main.cpp:317; src/main.cpp:320; src/main.cpp:324; src/nocspace.cpp:284; src/nocspace.cpp:290 -->

In the current bath-only rotation setup, the first division contains the impurity and the middle division is treated as active. Occupied-side bath divisions constrain holes; empty-side divisions constrain particles. Negative restraint entries express hole limits, and positive entries express particle limits. Zero is position-dependent, so do not treat every zero as a free parameter. The tutorial templates retain the required first and middle entries.

<!-- src/prmtr.cpp:69; src/prmtr.cpp:135; src/prmtr.cpp:155; src/nocspace.cpp:373; src/nocspace.cpp:386; src/nocspace.cpp:770 -->

The code recognizes these constraint names:

- `nooc`: applies the individual column checks, without the extra correlated-column checks.
- `cpnooc`: adds paired hole/particle checks for a subset of the divisions.
- `cnooc`: adds those paired checks across the occupied-side divisions.
- `phess`: selects the `check_if_PHSs` constraint.
- `phss_v2`: selects cumulative hole-plus-particle bounds across paired bath divisions. It is the compiled default and the tutorial choice.

<!-- src/nocspace.cpp:79; src/nocspace.cpp:95; src/nocspace.cpp:122; src/nocspace.cpp:296; src/nocspace.cpp:792 -->

Mode 1 first estimates natural-orbital occupations $n$ and uses weights $\min(n,1-n)$, setting very small weights below $10^{-14}$ to zero. It compares the weights against `weight_freze` first and then `weight_nooc`. The resulting eight-entry row is `[1, frozen_occupied, constrained_occupied, kept_occupied, 1, kept_empty, constrained_empty, frozen_empty]`. The actual allowed excitations still depend on `restrain`; a threshold alone does not define the whole many-body space.

<!-- src/api_edmft.cpp:413; src/api_edmft.cpp:415; src/api_edmft.cpp:423; src/api_edmft.cpp:673 -->

## Degeneracy and magnetic input

The name `pred_gs_deg` suggests a count, but its behavior differs between modes. Mode 0 sets both the initial stored degeneracy and the branch selector to this integer. Mode 1 sets only the branch selector and then resets the stored degeneracy to zero. `DensityMat::update` takes the explicit multi-state branch only when the selector is exactly `1`; `0`, `2`, and other values take its default eigensolve branch. Thus the supplied mode-1 value `2` does **not** impose a twofold ground state. Later eigensolves can increase the stored degeneracy when they detect more states.

<!-- src/main.cpp:248; src/main.cpp:249; src/api_edmft.cpp:498; src/api_edmft.cpp:291; src/densitymatrix.cpp:102; src/densitymatrix.cpp:111; src/operator.cpp:331 -->

`magnetic` provides separate input/output spin channels. For this path, use `deg_idx [1, 2, ..., nband]` and one `Ed` per band. The current natural-orbital rotation still averages the two bath density matrices of each band, so separate input channels do not imply a fully independent spin-dependent rotation. The supplied tutorials use `magnetic 0`.

<!-- src/api_edmft.cpp:151; src/api_edmft.cpp:155; src/densitymatrix.cpp:72 -->

## Units

The solver performs no unit conversion. Use one energy unit consistently for impurity levels, hoppings, interactions, hybridization, and frequency; `beta` has the inverse unit. The one-shot tutorial labels its inputs in eV. Bethe tutorial energies are model units. Green's functions have inverse-energy units, while self-energy and hybridization have energy units. Particle counts, occupation weights, and restraint entries are dimensionless.

<!-- tutor/1_Single-impurity_Anderson_model_with_multiple_orbitals/PARAMS.norg:7; src/prmtr.cpp:249; src/impurity.cpp:41; src/impurity.cpp:88; src/api_edmft.cpp:70 -->
