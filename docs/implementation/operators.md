# Rapid generation of many-body operators

The many-body Hamiltonian is large but sparse. Most determinants cannot be connected by a single one-body or interaction term. iNORG enumerates allowed transitions and caches their relation to the single-particle coefficients, then assembles a sparse matrix for Lanczos.

## Hamiltonian construction

For a Hamiltonian written as

$$
H=\sum_{ij}t_{ij}c_i^\dagger c_j
+\sum_{ijkl}U_{ijkl}c_i^\dagger c_j^\dagger c_kc_l,
$$

the matrix element between determinants is $H_{mn}=\langle m|H|n\rangle$. Starting from each retained determinant, the construction finds the allowed operator actions, calculates their fermion signs, and locates the destination determinants using the [configuration encoding](encoding.md). Terms whose destination lies outside the selected NOOC space do not enter the projected Hamiltonian.

Several operator terms can contribute to the same matrix element. Their coefficients must be accumulated with the appropriate signs. The current bath-only path handles one-body hopping connections through its table and calculates diagonal one-body and density-density interaction contributions separately.

<!-- iNORG paper: Rapid generation of many-body operators / Hamiltonian construction, inorg.tex:639. src/state.cpp:17; src/operator.cpp:173; src/operator.cpp:245; src/operator.cpp:280. -->

## Compressed sparse row storage

Compressed sparse row (CSR) format stores three arrays: the start of each row, the column index of each stored element, and its value. The paper gives

$$
A=\begin{pmatrix}1&0&0\\0&2&0\\3&0&4\end{pmatrix}
$$

as a small example:

| CSR array | Values |
| --- | --- |
| Row pointers | `[0, 1, 2, 4]` |
| Column indices | `[0, 1, 0, 2]` |
| Nonzero values | `[1, 2, 3, 4]` |

The final row pointer marks the total number of stored elements. A Hamiltonian-vector product visits these stored entries, rather than all pairs of determinants. MPI partitions the Hamiltonian rows among ranks; the resulting sparse operator supplies the products required by the eigensolver.

<!-- iNORG paper: CSR footnote in Hamiltonian construction. src/operator.cpp:237; src/operator.cpp:298; src/operator.cpp:317. -->

## NORG-Table

The NORG-Table separates operator connectivity from numerical coefficients. A record contains a many-body row, a many-body column, and an operator tag. Reusing these records avoids finding the same destination determinant after every coefficient update.

For $N_{\mathrm{sp}}$ spin-orbitals, the paper defines tags with one-based orbital labels:

$$
t_{\mathrm{tag}}=
\begin{cases}
(i-1)N_{\mathrm{sp}}+(j-1),
&c_i^\dagger c_j,\\[2mm]
N_{\mathrm{sp}}^2+
\bigl([(i-1)N_{\mathrm{sp}}+(j-1)]N_{\mathrm{sp}}+(k-1)\bigr)
N_{\mathrm{sp}}+(l-1),
&c_i^\dagger c_j^\dagger c_kc_l.
\end{cases}
$$

The one-body and four-fermion tags occupy distinct ranges below $N_{\mathrm{sp}}^2+N_{\mathrm{sp}}^4$. The tag identifies an operator coefficient; it is not a wave-function index.

The code uses zero-based orbital labels and stores a **signed, one-offset tag**. For a hopping term, its magnitude is `i * Nsp + j + 1`; its sign stores the fermion sign. The offset lets even the first coefficient carry a sign. Four-fermion coefficients follow the one-body range with the same offset convention. Each rank's table rows are local row indices, while destination columns refer to the full selected space.

<!-- iNORG paper: NORG-Table data structure, inorg.tex:661. Paper tags start at zero and omit an explicit sign. Code: src/operator.cpp:179; src/operator.cpp:185; src/operator.cpp:199; src/operator.cpp:221; src/operator.cpp:285. -->

During an orbital-optimization sweep, transformed coefficients are inserted into the existing table structure. A change in allowed configurations or index ordering requires a new table. A numerical orbital rotation can reuse the table when the selected occupation basis and possible transitions retain the same labels. The current NORG loop constructs its operator object once and updates coefficients during the sweep.

<!-- Paper wording: rebuild when the subspace or orbital basis changes. The code reuses connectivity across numerical rotations with fixed divisions: src/norg.cpp:56; src/norg.cpp:72; src/norg.cpp:77; src/norg.cpp:80; src/operator.cpp:11. -->

## Bath-only specialization and observables

Leaving the impurity basis fixed keeps density-density interactions diagonal in the occupation basis. The default table therefore records off-diagonal hopping connections. `find_hmlt_V2` calculates the diagonal energies and combines them with those connections to construct the Hamiltonian.

This specialization matters when reusing a table for another observable. The one-particle density matrix is calculated with a separate diagonal contribution. The current double-occupancy routine instead expects four-fermion table entries that are absent in this path, producing an unreliable zero block. Do not use it to interpret correlations.

<!-- src/prmtr.cpp:69; src/operator.cpp:11; src/operator.cpp:203; src/operator.cpp:215; src/operator.cpp:253; src/operator.cpp:261; src/densitymatrix.cpp:138; src/densitymatrix.cpp:204; src/densitymatrix.cpp:213. -->

## References

- Intel Corporation (2024), *OneMKL Developer Reference: Sparse Storage Formats*.
- J.-M. Wang et al., *iNORG: An open-source quantum impurity solver package based on the natural orbitals renormalization group*, “Rapid generation of many-body operators.”
