# Output files

Both standalone modes write into the directory from which `inorg` is launched. Names are case-sensitive. Mode 0 normally includes an iteration prefix; mode 1 normally overwrites fixed filenames. Save results before another run.

<!-- src/main.cpp:37; src/green.cpp:103; src/api_edmft.cpp:67 -->

## Frequency-table layouts and units

All frequency files produced by the active tutorial paths are on the positive Matsubara mesh, $\omega_n=(2n+1)\pi/\beta$. Their columns contain diagonal matrix elements only.

**Mode 0 layout (A):** a text header followed by one row per frequency. For $B$ spatial bands there are $2(B+1)$ columns:

```text
omega  ReF_11 ... ReF_BB  omega  ImF_11 ... ImF_BB
```

The frequency appears twice. For the two-band tutorial, the imaginary parts are columns 5 and 6, not alternating with the real parts. Headers contain the quantity name, iteration, and matrix index.

<!-- src/green.cpp:99; src/green.cpp:108; src/green.cpp:116; src/green.cpp:124 -->

**Mode 1 layout (B):** no header; each row is frequency followed by real/imaginary pairs:

```text
omega  ReF_1 ImF_1  ReF_2 ImF_2 ...
```

For `magnetic 0`, there is one pair per degeneracy group, in the ordered groups described by `deg_idx`. For `magnetic 1`, there is one pair per spin-orbital channel: all up-spin bands, then all down-spin bands. The supplied three-band example has seven columns; a magnetic three-band input would have thirteen.

<!-- src/green.cpp:188; src/green.cpp:197; src/green.cpp:305; src/api_edmft.cpp:20 -->

Frequencies, hybridizations $\Delta$, and self-energies $\Sigma$ have energy units. Green's functions $G$ and $G_0$ have inverse-energy units. The code applies no unit conversion: mode 0 uses the energy unit chosen for its model parameters; the supplied mode-1 input uses eV. Occupations and particle numbers are dimensionless. Fit-error measures are normalized diagnostics, not energies or uncertainty estimates.

<!-- src/prmtr.cpp:249; src/impurity.cpp:41; src/impurity.cpp:88; src/api_edmft.cpp:70; src/hyberr.h:73; tutor/1_Single-impurity_Anderson_model_with_multiple_orbitals/PARAMS.norg:7 -->

## Mode 0 files

`NNN` is the outer DMFT iteration number, padded to three digits. Initial files use `000`; impurity solutions start at `001`. These are the actual names, including `.mb.` where applicable.

| File | Written when | Contents and layout |
| --- | --- | --- |
| `log.txt` | Start and each self-energy update; appends | Iteration, latest relative error, hyphen-separated particle counts for the spin-orbital sets, one occupation per impurity spin orbital, timestamp, event label. See log notes below. <!-- src/dmft_bethe.h:60; src/dmft_bethe.h:80; src/dmft_bethe.cpp:46; src/dmft_bethe.cpp:85 --> |
| `log.parameters.txt` | Once per invocation; appends | Labeled derived settings, including process count, bath count, constraint matrix, interactions, frequency scales, and fit settings. Mixed scalar/vector/matrix text, not a fixed-column table. <!-- src/dmft_bethe.cpp:24; src/prmtr.cpp:260 --> |
| `zic000.mb.g0loc.txt` | Fresh start without `mb.gfloc` | Noninteracting local lattice Green's function; layout A, inverse energy. <!-- src/dmft_bethe.cpp:50; src/green.cpp:103 --> |
| `zic000.mb.gfloc.txt` | Start with `mb.gfloc` | The read-in local Green's function; layout A, inverse energy. <!-- src/dmft_bethe.cpp:49; src/green.cpp:103 --> |
| `zic000.mb.seimp.txt` | Start with `mb.seimp` | The read-in self-energy; layout A, energy. <!-- src/dmft_bethe.cpp:48; src/green.cpp:103 --> |
| `zicNNN.mb.hb.txt` | Before each bath fit | Target hybridization; layout A, energy. Absent at iteration 1 when an existing `ose_hop` supplies the bath directly. <!-- src/dmft_bethe.cpp:60; src/dmft_bethe.cpp:61; src/green.cpp:103 --> |
| `zicNNN.mb.hb-fit.txt` | After impurity construction | Hybridization reconstructed from the discrete bath; layout A, energy. <!-- src/dmft_bethe.cpp:67; src/impurity.cpp:76; src/green.cpp:103 --> |
| `zicNNN.mb.g0imp.txt` | Each impurity solution | Noninteracting impurity Green's function from the stored one-body matrix; layout A, inverse energy. <!-- src/dmft_bethe.cpp:71; src/green.cpp:103 --> |
| `zicNNN.mb.gfimp.txt` | Each impurity solution | Interacting impurity Green's function; layout A, inverse energy. <!-- src/dmft_bethe.cpp:72; src/green.cpp:103 --> |
| `zicNNN.mb.seimp.txt` | Each impurity solution | Impurity self-energy $G_0^{-1}-G^{-1}$ before the outer self-energy update, with real diagonal parts set to zero when `bath_ph_symm 1`; layout A, energy. <!-- src/dmft_bethe.cpp:73; src/dmft_bethe.cpp:74; src/green.cpp:103 --> |
| `zicNNN.mb.gfloc.txt` | After the self-energy update | Local lattice Green's function built using the updated self-energy; plain substitution by default, with Pulay mixing used only with `pulay_mixing 1` once its activation condition is met; layout A, inverse energy. <!-- src/dmft_bethe.cpp:88; src/dmft_bethe.cpp:91; src/dmft_bethe.cpp:96; src/green.cpp:103 --> |
| `zicNNN.ose_hop.txt` | After each bath fit | Bath parameter blocks described below; energies. It is not written for a first iteration that skips fitting. <!-- src/dmft_bethe.cpp:60; src/dmft_bethe.cpp:63; src/bath.cpp:336 --> |
| `zicNNN.h0.txt` | After impurity construction | Fit diagnostics and stored one-body matrix blocks; see below. <!-- src/dmft_bethe.cpp:65; src/impurity.cpp:208 --> |
| `zicNNN.nmat.txt` | After each impurity solution | Occupation and correlation blocks described below. <!-- src/dmft_bethe.cpp:70; src/norg.cpp:746; src/norg.cpp:762 --> |

The files often described informally as `hb.txt`, `hb-fit.txt`, `gfloc`, `gfimp`, `g0imp`, `nmat.txt`, `ose_hop.txt`, and `h0.txt` carry the prefixes shown above in this mode. There is no automatically updated unprefixed `mb.gfloc`, `mb.seimp`, or `ose_hop` restart file. The bath constructor calls its writer with iteration zero, but that writer opens no file for zero, so there is no `zic000.ose_hop.txt` from that call.

<!-- src/dmft_bethe.cpp:33; src/bath.cpp:33; src/bath.cpp:332; src/bath.cpp:336 -->

## Mode 1 files

| File | Written when | Contents and layout |
| --- | --- | --- |
| `hb_read.txt` | After reading `Delta.inp` | Hybridization as represented on the internal mesh; layout B, energy. It is a read-back diagnostic, not a byte copy of the input. <!-- src/api_edmft.cpp:26; src/api_edmft.cpp:29; src/green.cpp:188; src/green.cpp:305 --> |
| `hb_fit.txt` | After fitting | Hybridization reconstructed from the fitted bath; layout B, energy. Note the underscore, unlike mode 0's `hb-fit`. <!-- src/api_edmft.cpp:43; src/api_edmft.cpp:44; src/impurity.cpp:76 --> |
| `ose_hop` | After fitting | Fitted bath blocks; creates the file or overwrites an existing seed. <!-- src/api_edmft.cpp:32; src/api_edmft.cpp:35; src/bath.cpp:333 --> |
| `h0.txt` | After impurity construction | Fit diagnostics and one-body matrix blocks; energies in the matrices. <!-- src/api_edmft.cpp:39; src/impurity.cpp:207 --> |
| `nmat.txt` | After the final impurity solve | Occupation and correlation blocks. <!-- src/api_edmft.cpp:52; src/norg.cpp:759 --> |
| `g0imp.txt` | Before interacting Green's-function output | Noninteracting impurity Green's function; layout B, inverse energy. <!-- src/api_edmft.cpp:58; src/api_edmft.cpp:61 --> |
| `Gf.out` | After Green's-function calculation | Interacting impurity Green's function; layout B, inverse energy. <!-- src/api_edmft.cpp:64; src/api_edmft.cpp:67 --> |
| `Sig.out` | At the end | Self-energy $G_0^{-1}-G^{-1}$; layout B, energy. An existing file is read first to calculate the reported change. <!-- src/api_edmft.cpp:70; src/api_edmft.cpp:73; src/api_edmft.cpp:80 --> |
| `log.norg` | At the self-energy comparison; appends | Counter, relative self-energy change, particle-number sector label, per-band `fit_error~occupation` pairs, timestamp, event label. Magnetic output groups all up-spin bands before down-spin bands. <!-- src/api_edmft.h:98; src/api_edmft.cpp:314; src/api_edmft.cpp:338; src/api_edmft.cpp:342 --> |
| `edmft_back_up` | After the impurity solve | Two labeled rows: `iter_count` followed by the incremented invocation counter; `npartical` followed by one integer per interleaved spin-orbital set. Spelling is literal. Read again on the next invocation. <!-- src/api_edmft.cpp:46; src/api_edmft.cpp:55; src/api_edmft.cpp:620; src/api_edmft.cpp:643; src/api_edmft.cpp:646 --> |

## Bath blocks: ose_hop

Both bath writers use the same text structure:

```text
set[0]
the
V1 V2 ... VB
epi
epsilon1 epsilon2 ... epsilonB
```

Each numeric row has exactly the bath count for that part. `the` is the row of real hybridization amplitudes; `epi` is the row of bath energies. Values in the same column describe the same bath site. There is one block per spatial band for spin-degenerate fitting. Magnetic mode uses separate blocks for all up-spin bands followed by all down-spin bands.

<!-- src/bath.cpp:8; src/bath.cpp:18; src/bath.cpp:329; src/bath.cpp:423; src/bath.cpp:448 -->

In mode 0, an iteration's `zicNNN.ose_hop.txt` can be used as `ose_hop` for a compatible restart. In mode 1, `ose_hop` is already both input and output. Always check the bath count; the reader does not validate it.

## One-body matrices: h0

The `h0` files start with fit diagnostics for each band or degeneracy group. The labels are `nmin`, `err`, `err_crv`, `err_hight_ev`, `err_regE`, and `norm`; mode 1 also prints `nbath`. `nmin` counts fits near the best minimum. The error terms describe the total objective, hybridization curve, high-frequency term, and bath-energy regularization. `norm` is the norm of the combined hopping/energy parameter vector.

<!-- src/impurity.cpp:210; src/bath.cpp:70; src/bath.cpp:184; src/hyberr.cpp:86; src/hyberr.h:73; src/hyberr.h:85; src/hyberr.h:98 -->

Each following `band` label introduces a square matrix. For $B$ baths per spin there are $2(B+1)$ rows and columns: the up-spin impurity and baths, then the down-spin impurity and baths. Within each spin block, the impurity is first and bath sites are sorted by energy. Matrix values have energy units. Each row begins with `row[i] =`, with a zero-based index, followed by its values. The matrix printer emits only a size message if there are more than 1024 rows or 64 columns.

In mode 0, the printed matrix is the stored `h0`; the separate half-filling shift is applied to the solver's Hamiltonian copy. Account for that convention when comparing its diagonal with the interacting model.

<!-- src/impurity.cpp:189; src/impurity.cpp:193; src/impurity.cpp:222; src/impurity.cpp:229; src/impurity.cpp:306; src/bath.cpp:425; src/bath.cpp:440; src/bath.cpp:449; src/gen/mat.h:142 -->

## Occupations and correlations: nmat

The first block has two columns: a spin-orbital label such as `1up` or `1dn`, followed by $\langle n_i\rangle$. Rows `sup`, `sdn`, and `sum` give spin-up, spin-down, and total occupations. These quantities are dimensionless.

<!-- src/norg.cpp:787; src/norg.cpp:790; src/norg.cpp:795 -->

The subsequent labeled blocks are matrices for $\langle s_{z,a}s_{z,b}\rangle$, $\langle N_aN_b\rangle$, connected orbital-number correlations, $\langle n_i n_j\rangle$, and connected spin-orbital-number correlations. Each matrix row starts with `row[i] =` and then its column values. Orbital matrices follow band order; spin-orbital matrices follow `1up, 1dn, 2up, 2dn, ...`. There are also scalar values labeled $\langle S_z^2\rangle$ and $\langle S^2\rangle$. Here $N_a=n_{a\uparrow}+n_{a\downarrow}$, and the spin convention uses $s_z=(n_\uparrow-n_\downarrow)/2$.

<!-- src/norg.cpp:801; src/norg.cpp:814; src/norg.cpp:827; src/gen/mat.h:155 -->

!!! warning "Double-occupancy output"

    The current operator-table path leaves the printed $\langle n_i n_j\rangle$ matrix zero. It is not a physical result. Do not use this block or the spin/charge correlation quantities derived from it. The one-particle occupation block is computed separately. See [known limitations](known-limitations.md#double-occupancy-output).

<!-- src/norg.cpp:752; src/norg.cpp:799; src/norg.cpp:827; src/prmtr.cpp:69; src/operator.cpp:11; src/operator.cpp:199; src/operator.cpp:203; src/operator.cpp:215; src/densitymatrix.cpp:204; src/densitymatrix.cpp:213 -->

## Logs, restarts, and completion

`run.log` in the tutorials is created by shell redirection; it is not an additional solver-managed file. It contains stdout and stderr, including bath-fit progress, natural occupations, ground-state energies, and warnings. Its format is diagnostic text, not a stable numeric table.

<!-- src/bath.cpp:79; src/norg.cpp:95; src/norg.cpp:111; src/main.cpp:25 -->

Mode 0's `log.txt` initial error is a sentinel, not a measured residual. Later `sigerr_update` rows compare the current input self-energy with the newly computed impurity self-energy, before the update. With the default `pulay_mixing 0`, plain substitution makes these successive impurity self-energies. Pulay mixing is used only with `pulay_mixing 1`; it starts once an iteration after the first has error below `8e-3` and remains active thereafter. The next input may then be a mixed self-energy. Compare several iterations of the curves and occupations as well as this number. A stable error or a stopped loop alone is insufficient to establish convergence.

<!-- src/main.cpp:155; src/main.cpp:281; src/dmft_bethe.cpp:29; src/dmft_bethe.cpp:46; src/dmft_bethe.cpp:84; src/dmft_bethe.cpp:88; src/dmft_bethe.cpp:91; src/dmft_bethe.cpp:394 -->

Mode 1's first comparison with `Sig.out` is not meaningful when there was no earlier file. On later invocations, keep `Sig.out`, `edmft_back_up`, and the bath only when you intend to continue the same calculation. Optional `ru<sector>.bi` rotation files can be read if supplied, but the active tutorial paths do not save new ones. Real-axis spectra and the additional development-only diagnostic writers are not called by these tutorial paths.

<!-- src/api_edmft.cpp:50; src/api_edmft.cpp:53; src/api_edmft.cpp:73; src/api_edmft.cpp:84; src/norg.cpp:604; src/norg.cpp:614; src/norg.cpp:63; src/occler.cpp:56; src/dmft_bethe.cpp:100 -->
