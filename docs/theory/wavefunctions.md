# Representation of many-body wave functions

A single-particle basis specifies the orbitals. A many-body basis specifies which of those orbitals are occupied. NORG changes the first basis to make the second basis smaller after applying occupation constraints.

## Fock states and particle-number sectors

For $M$ spatial orbitals there are $2M$ spin-orbitals, each with occupation $n_i\in\{0,1\}$. A determinant is labeled by its occupation vector $|n_0,n_1,\ldots,n_{2M-1}\rangle$. The full Fock space is

$$
\mathcal H=\mathcal H_0\oplus\mathcal H_1\oplus\cdots\oplus\mathcal H_{2M},
\qquad \dim\mathcal H=2^{2M}.
$$

Particle-number conservation allows the solver to work in a sector $\mathcal H_N$. Within a fixed sector, there are $\binom{2M}{N}$ unconstrained determinants. NOOC then retains only the determinants satisfying the chosen group limits. A wave function is stored as amplitudes on those retained determinants; the orbital occupations alone do not specify it.

In the density-density models with diagonal hybridization used here, the implementation further organizes the space by the particle count in each spin-orbital set. The preliminary mode-1 search compares sectors before the final impurity solve. Green's functions require neighboring sectors because creation or annihilation changes a particle count.

<!-- iNORG paper: Representation of many-body wave functions, inorg.tex:253; Many-body wave function search, inorg.tex:611. src/nocspace.cpp:258; src/occler.cpp:17; src/norg.cpp:276. -->

## Binary representation

Choose a fixed ordering of orbitals and let bit $i$ represent $n_i$. With orbital zero assigned to the least significant bit, the integer code is

$$
b=\sum_{i=0}^{2M-1}n_i2^i.
$$

The paper's example has six spin-orbitals and electrons on orbitals 2 and 3. Its occupation vector in increasing orbital order is $|0,0,1,1,0,0\rangle$, and its integer code is $2^2+2^3=12$.

An occupation vector printed from orbital zero upward is not generally a conventional binary numeral, which is printed from the most significant bit downward. The code's `Onb` printer uses increasing orbital order. Use the bit positions when interpreting a state rather than assuming that the leftmost printed digit is the highest bit.

<!-- iNORG paper: Representation of many-body wave functions, six-orbital example. The later C(6,2) table reads printed strings as conventional binary numbers, which reverses the earlier LSB convention. Follow src/gen/onb.h:61 and src/gen/onb.h:130 here. -->

The paper introduces a global ordering with all up-spin orbitals followed by all down-spin orbitals. The standalone solver instead groups sites into interleaved band/spin sets: `1up, 1dn, 2up, 2dn, ...`. Each set contains its impurity and bath divisions. Keep this distinction in mind when translating a schematic state into input arrays or internal indices. Frequency-file column order is a separate convention described in [output files](../output-files.md).

<!-- Paper wording: all up-spin orbitals first, then all down-spin orbitals. Code: src/api_edmft.cpp:160; src/state.cpp:25; src/state.cpp:30; src/prmtr.cpp:135. -->

## Fermion signs

Creation requires an empty orbital, and annihilation requires an occupied orbital. The implementation checks which transitions are allowed, flips the relevant bits, and attaches the sign from fermion ordering. For a single creation or annihilation, the sign depends on the parity of the occupied orbitals preceding that position. For a hop, it depends on the occupied orbitals strictly between the two positions.

The bit operation and the fermion sign are distinct parts of the operator action. Ignoring the sign would produce an incorrect many-body Hamiltonian even if the destination occupation vector were correct.

<!-- src/gen/onb.h:65; src/gen/onb.h:85; src/state.cpp:27; src/state.cpp:37; src/state.cpp:43. -->

## From bit strings to compact storage

Integer bit codes are sparse within a fixed-particle sector. They are useful for applying operators but waste space if used directly as amplitude-array indices. iNORG uses combinatorial ranks within each orbital division and mixed-radix indices across divisions. The [encoding page](../implementation/encoding.md) works through the paper's examples and explains how the solver locates the destination of an operator action.

## References

- J.-M. Wang et al., *iNORG: An open-source quantum impurity solver package based on the natural orbitals renormalization group*, “Representation of many-body wave functions” and “Encoding of many-body wave functions.”
