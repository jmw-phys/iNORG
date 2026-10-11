# Encoding many-body wave functions

iNORG needs to convert quickly between a determinant and its location in a wave-function array. Binary occupations make operator actions inexpensive. Combinatorial ranks and mixed-radix indices make storage compact. Particle-number sectors and NOOC determine which combinations exist in the selected space.

## Combinatorial number encoding

For $k$ electrons in $n$ spin-orbitals, there are $\binom nk$ configurations. Their bit codes are not consecutive, but their ranks can be consecutive integers from zero to $\binom nk-1$. Rank the sorted occupied-orbital lists in lexicographic order.

For the paper's $\binom62=15$ example:

| Rank | Occupied orbitals | Occupations in increasing orbital order |
| --- | --- | --- |
| 0 | $\{0,1\}$ | $(1,1,0,0,0,0)$ |
| 1 | $\{0,2\}$ | $(1,0,1,0,0,0)$ |
| 9 | $\{2,3\}$ | $(0,0,1,1,0,0)$ |
| 13 | $\{3,5\}$ | $(0,0,0,1,0,1)$ |
| 14 | $\{4,5\}$ | $(0,0,0,0,1,1)$ |

The state on orbitals 2 and 3 has bit code 12 but rank 9. Only the rank is a compact array index. The distinction also matters for NOOC: $6^{2-}$ allows **up to** two electrons, so it includes the zero- and one-electron sectors as well as these fifteen two-electron configurations.

<!-- iNORG paper: Encoding of many-body wave functions / Combinatorial number encoding, inorg.tex:509, table c62_configurations. The table above retains the paper's ranks and occupations but omits its inconsistent decimal column. Bit convention: src/gen/onb.h:61; src/gen/onb.h:130. -->

### Ranking and unranking

At a candidate occupied position, count how many configurations would occur if that position were chosen and the remaining particles occupied later positions. Binomial coefficients provide this count. Ranking adds the blocks skipped before reaching the actual occupied positions. Unranking subtracts whole blocks until it finds the one containing the requested rank.

The following zero-based pseudocode states the common counting operation. `occupied` is sorted, and `C(a,b)` denotes a binomial coefficient.

```text
rank = 0
previous = -1
for i = 0, ..., k-1:
    for j = previous+1, ..., occupied[i]-1:
        rank += C(n-j-1, k-i-1)
    previous = occupied[i]
```

For unranking, start with the requested rank and test candidate positions in increasing order. If the remaining rank is at least `C(n-j-1, k-i-1)`, subtract that block and advance the candidate. Otherwise occupy that position and continue with one fewer particle. `Onb` implements the equivalent include-or-skip construction, with special cases for empty, full, one-electron, and one-hole divisions.

<!-- Paper wording: ranking/unranking pseudocode uses C(n-j,k-i) with one-based i, inconsistent with its zero-based orbital labels. The count here excludes the candidate orbital, as in src/gen/onb.h:146; src/gen/onb.h:157; src/gen/onb.h:180; src/gen/onb.h:192. -->

Neither direction requires a stored list of every bit string. `BnmlFast` supplies the binomial counts. `Onb` stores an individual division's occupation pattern; it does not place the whole impurity problem in a single integer bit field.

<!-- src/gen/onb.h:35; src/gen/onb.h:38; src/gen/comdivs.cpp:47. -->

## Mixed-radix encoding across groups

A mixed-radix representation assigns a different base to each digit. Here a digit is a group's configuration rank, and its base is the number of configurations allowed in that group.

The paper illustrates this with $6^{2+}6^{2-}$. Each group permits zero, one, or two changes relative to full or empty occupancy, giving

$$
R=\binom60+\binom61+\binom62=1+6+15=22.
$$

If the first and second group ranks are $A$ and $B$, and the second group is the less significant digit, the combined code is

$$
q=A R_B+B.
$$

The paper's example $A=3$, $B=14$, $R_B=22$ gives $q=80$. Recover the digits by integer quotient and remainder. If the second group instead contains exactly two electrons, its base is fifteen, not twenty-two.

<!-- iNORG paper: Mixed-radix numeral system, inorg.tex:583. -->

The implementation first separates the space into fixed occupation-count patterns. Within one pattern, division $i$ contains exactly $N_i$ particles on $m_i$ sites, so its base is $\binom{m_i}{N_i}$. `ComDivs` makes the **first stored division** the least significant digit. It constructs cumulative base products and sums each division rank times its product. Decoding uses division and remainder. This is the same mixed-radix idea with a different digit order from the two-group illustration.

<!-- Paper illustration: second group is least significant; current code: src/gen/comdivs.cpp:47; src/gen/comdivs.cpp:93; src/gen/comdivs.cpp:107; src/gen/comdivs.cpp:120. -->

## Searching the allowed space

Particle-number conservation gives $[\hat N,H]=0$. The paper's search example uses sixteen spin-orbitals at half filling, with

$$
6^{2+}\;4^*\;6^{2-},\qquad N=8.
$$

The active group contains the impurity and active bath orbitals. Enumerate the first and third group counts, then obtain the active count from conservation:

$$
N_1\in\{4,5,6\},\qquad N_3\in\{0,1,2\},
\qquad N_2=8-N_1-N_3.
$$

Keep only patterns with $0\leq N_2\leq4$. For the pattern $(N_1,N_2,N_3)=(5,2,1)$, the number of determinants is

$$
\binom65\binom42\binom61=6\times6\times6=216.
$$

<!-- iNORG paper: Many-body wave function search, inorg.tex:611. -->

The implementation applies this idea across all sets and divisions. It stores allowed occupation-count patterns with their starting offsets. To identify a determinant, locate its pattern and add its within-pattern mixed-radix index. To apply an operator, update occupations, reject destinations outside the selected space, and rank the retained destination. The fermion sign is calculated separately.

This lets [Hamiltonian construction](operators.md) traverse allowed transitions without allocating a dense matrix or scanning all unconstrained Fock states.

<!-- src/nocspace.cpp:759; src/nocspace.cpp:860; src/state.cpp:8; src/state.cpp:24; src/state.cpp:40; src/state.cpp:43. -->

## References

- J.-M. Wang et al., *iNORG: An open-source quantum impurity solver package based on the natural orbitals renormalization group*, “Encoding of many-body wave functions.”
