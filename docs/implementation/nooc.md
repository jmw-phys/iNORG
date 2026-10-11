# Natural orbital occupancy constraints

Natural orbital occupancy constraints (NOOC) select the determinants retained in a NORG calculation. You group orbitals according to their natural occupations, then limit departures from nearly full or nearly empty groups. The active impurity and bath divisions retain unrestricted occupations, subject to the total particle counts.

## Group notation

Each orbital in this notation is a **spin-orbital**, with capacity one electron.

| Notation | Meaning for a group of $m$ spin-orbitals |
| --- | --- |
| $m^*$ | Active: from zero to $m$ electrons are allowed. |
| $m^{l+}$ | At most $l$ holes relative to full occupancy; at least $m-l$ electrons. |
| $m^{l-}$ | At most $l$ electrons relative to the empty state. |
| $m^{l\pm}$ | At most $l$ changes relative to a preset full/empty pattern; count holes and electrons together. |

The signs in the paper's superscripts are **not the signs in `restrain`**. The input uses a negative number for a hole limit and a positive number for an electron limit. Division position determines whether a zero means an active division or a frozen bath division.

<!-- iNORG paper: Natural orbital occupancy constraint, inorg.tex:407. Input signs: src/nocspace.h:75; src/nocspace.cpp:296. -->

## Direct-product rule

The paper considers one impurity spatial orbital coupled to fifteen bath spatial orbitals. These give two impurity spin-orbitals and thirty bath spin-orbitals. Group the bath as

$$
8^{1+}\;6^{2+}\;2^*\;6^{2-}\;8^{1-}.
$$

The groups contain nearly full, moderately full, active, moderately empty, and nearly empty orbitals. Each group has its own constraint. Their direct product, together with the active impurity space, gives the retained space before imposing the fixed particle counts. The split into five bath groups is an example, not a requirement of the direct-product rule.

The corresponding algorithm name in the input is `nooc`. It performs the individual division checks without adding a joint excitation limit.

<!-- iNORG paper: Direct product constraint rule, inorg.tex:423. src/nocspace.cpp:296; src/nocspace.cpp:297. -->

## Joint rule

Define a change quantity as

$$
x=n_h+n_e,
$$

where $n_h$ counts holes on the preset full side and $n_e$ counts electrons on the preset empty side. Pairing the occupied and empty groups gives the paper's example

$$
2^*\;12^{2\pm}\;16^{1\pm}.
$$

The inner constrained group contains six preset full and six preset empty orbitals. It allows at most two changes in total. The outer group contains eight of each and allows at most one change. A hole and an electron each consume one unit of this budget; a particle-hole pair consumes two.

The input name is `cnooc`. The code first checks each division separately, then checks the combined holes and electrons in each mirrored pair of bath divisions.

<!-- iNORG paper: Joint constraint rule, inorg.tex:444. src/nocspace.cpp:79; src/nocspace.cpp:299. -->

## PHES rule

The particle-hole excitation state (PHES) rule also limits the cumulative excitations across groups. In the paper's example, the outer group allows at most one change, and the inner and outer groups together allow at most two:

$$
2^*\left(12+16^{1\pm}\right)^{2\pm}.
$$

This rejects the last configuration in the following table, even though it passes the joint rule.

| Inner-group changes | Outer-group changes | Total changes | Joint | PHES |
| --- | --- | --- | --- | --- |
| 0 | 0 | 0 | Keep | Keep |
| 1 | 0 | 1 | Keep | Keep |
| 2 | 0 | 2 | Keep | Keep |
| 0 | 1 | 1 | Keep | Keep |
| 1 | 1 | 2 | Keep | Keep |
| 2 | 1 | 3 | Keep | Reject |

Active-group occupations are unrestricted in every row. The nested rule removes combinations of individually allowed excitations that together exceed the desired budget. The paper identifies these high-excitation configurations as having small ground-state weight in its illustration.

Use `phss_v2` for this cumulative rule. The older token `phess` selects a different, specialized test; it is not an alias for `phss_v2`. Another token, `cpnooc`, pairs only a subset of the divisions. For the three rules described here, use `nooc`, `cnooc`, or `phss_v2`.

<!-- iNORG paper: PHES constraint rule, inorg.tex:464, table change_classification. src/nocspace.cpp:95; src/nocspace.cpp:129; src/nocspace.cpp:298; src/nocspace.cpp:300; src/nocspace.cpp:301. -->

## Eight-division layout in PARAMS.norg

The standalone bath-only setup uses this layout for **each spin-orbital set**:

| Zero-based division | Contents | Tutorial `restrain` entry | Mirrored pair |
| --- | --- | --- | --- |
| 0 | Impurity, active | `0` | — |
| 1 | Outermost nearly full bath | `-2` | 7 |
| 2 | Next nearly full bath | `-3` | 6 |
| 3 | Inner occupied-side bath | `-4` | 5 |
| 4 | Active bath | `0` | — |
| 5 | Inner empty-side bath | `4` | 3 |
| 6 | Next nearly empty bath | `3` | 2 |
| 7 | Outermost nearly empty bath | `2` | 1 |

`distribute` contains the **number of orbitals** in each division. `restrain` contains the **allowed changes**, summed over all spin-orbital sets in a division. The excitation budget is not repeated independently for each band or spin. The solver also enforces the particle number of each set.

<!-- src/main.cpp:317; src/main.cpp:320; src/main.cpp:324; src/nocspace.cpp:261; src/nocspace.cpp:292; src/nocspace.cpp:939. -->

Let $x_1$ count holes plus electrons in divisions $(1,7)$, $x_2$ in $(2,6)$, and $x_3$ in $(3,5)$. With the tutorial limits:

- `nooc` limits holes and electrons separately, with outer-to-inner bounds 2, 3, and 4 on each side.
- `cnooc` additionally requires $x_1\leq2$, $x_2\leq3$, and $x_3\leq4$.
- `phss_v2` requires $x_1\leq2$, $x_1+x_2\leq3$, and $x_1+x_2+x_3\leq4$, in addition to the individual division checks.

The positive-side entries supply the paired or cumulative bounds. Keep the symmetric sign pattern shown here when translating the paper's rules.

<!-- src/nocspace.h:81; src/nocspace.cpp:87; src/nocspace.cpp:112; src/nocspace.cpp:122; src/nocspace.cpp:296. These inequalities spell out the current eight-division implementation. -->

## Selecting each rule in mode 1

Start with the complete [one-shot input](../tutorials/siam-oneshot.md#input-parameters-line-by-line). Replace its `NOOC` and `restrain` lines with **one** of these alternatives; the snippets are not complete parameter files.

Direct product:

```text
NOOC      nooc
restrain  [0, -2, -3, -4, 0, 4, 3, 2]
```

Joint:

```text
NOOC      cnooc
restrain  [0, -2, -3, -4, 0, 4, 3, 2]
```

Cumulative PHES:

```text
NOOC      phss_v2
restrain  [0, -2, -3, -4, 0, 4, 3, 2]
```

Mode 1 reads these tokens directly. It does **not** read `distribute`. Instead, supply the bath counts and the occupation thresholds, as in the three-orbital example:

```text
fit_nbaths    [11, 11, 11, 11, 11, 11]
weight_nooc   [1e-05, 1e-05, 1e-05, 1e-05, 1e-05]
weight_freze  [1e-07, 1e-07, 1e-07, 1e-07, 1e-07]
```

Only the first three threshold entries are used for this example. After a preliminary solve, weights $\min(n,1-n)$ below `weight_freze` go into the outer divisions; weights below `weight_nooc` go into the next divisions. The remaining occupied-side and empty-side orbitals go into the inner constrained divisions. One impurity and one active bath orbital occupy divisions 0 and 4 of each set.

The generated row is

```text
[1, frozen_occupied, constrained_occupied, kept_occupied,
 1, kept_empty, constrained_empty, frozen_empty]
```

The word `frozen` names the threshold category. Those orbitals are strictly frozen only when their allowed excitations are zero. The tutorial permits excitations in its outer divisions. The preliminary sector search, or the preliminary fixed-sector solve when `ful_pcl_sch 0`, uses `cnooc` internally; your requested rule applies to the final space.

<!-- src/api_edmft.cpp:489; src/api_edmft.cpp:554; src/api_edmft.cpp:568 (distribute parser disabled); src/api_edmft.cpp:374; src/api_edmft.cpp:390; src/api_edmft.cpp:415; src/api_edmft.cpp:423; src/occler.cpp:23. -->

These inputs select the three algorithms with the current eight-division layout. They do not impose the paper's illustrative group sizes: those sizes depend on the calculated occupations in mode 1.

## Setting division sizes in mode 0

Mode 0 reads `restrain` and `distribute`. The supplied five-bath example uses

```text
fit_nbaths  [5, 5, 5, 5]
restrain    [0, -2, -3, -4, 0, 4, 3, 2]
distribute  [1, 0, 0, 2, 1, 2, 0, 0]
```

The row sums to six sites: one impurity and five baths. It is repeated for all four spin-orbital sets. Keep the same number of entries in both lists and make each row sum to the bath count plus one. Empty divisions are allowed.

**Mode 0 ignores the `NOOC` keyword.** Its compiled rule is `phss_v2`, even if you write `NOOC nooc` or `NOOC cnooc`. There is no parameter-file-only way to select the other two rules in this mode. Use mode 1 to select them through `PARAMS.norg`. The `NOOC phss_v2` line in the supplied Bethe input documents the default but does not set it.

<!-- Paper wording: Example 2 lists NOOC as a selectable parameter and describes distribute as orbital positions. The code reads counts and no NOOC token: src/main.cpp:168; src/main.cpp:227; src/main.cpp:240; src/main.cpp:278; src/main.cpp:297; src/main.cpp:324; src/prmtr.cpp:82. -->

## Choosing a space

Use the [parameter reference](../parameters.md#nooc-restrain-and-distribute) to check list sizes and defaults. Inspect the generated divisions and natural occupations in the log, then compare results with looser constraints. Ground-state energy, occupations, and Green's functions answer different convergence questions. A stable orbital rotation does not establish that a tight constraint is adequate.

## References

- J.-M. Wang et al. (2025), “Ab initio dynamical mean field theory with natural orbitals renormalization group impurity solver,” *npj Computational Materials* **11**, 86.
- J.-M. Wang et al., *iNORG: An open-source quantum impurity solver package based on the natural orbitals renormalization group*, “Natural orbital occupancy constraint.”
