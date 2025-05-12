# iNORG - First version

Now this package is currently under heavy development. 

## Version

v1.10.01

## Installation

To install iNORG, follow these steps:

1. Clone the repository:
```bash
$ git clone https://github.com/jmw-phys/iNORG.git
$ cd iNORG
```

2. Compile the code:
```bash
$ make clean
$ make
```

3. Run the program:
```bash
$ mpirun -n [number_of_processes] ./testing/inorg
```

Note: 
- The code requires MPI and MKL libraries
- Make sure you have the necessary compilers and libraries installed
- The number of processes should be adjusted based on your system configuration


## Tutorial — Two Typical Examples

Example 1: One-shot solution of multi-orbital impurity problem in DMFT

Example 2: Multi-orbital DMFT calculation on Bethe lattice


## License

GNU AFFERO General Public License Version 3

## Documentation

A tutorial will be provided in the `docs` folder. Additionally, a new tutorial article detailing the specific implementation algorithms and user manual for iNORG will soon be available on arXiv—stay tuned! [See author publications on arXiv](https://arxiv.org/search/?query=Wang%2C+Jia-Ming&searchtype=author&abstracts=show&order=-announced_date_first&size=50)


## Research Publications with iNORG

1. Jia-Ming Wang, *et al*., *Low-energy inter-band kondo bound states in orbital-selective mott phases*, [Phys. Rev. B, **111**, 155107 (2025)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.111.155107).

2. Jia-Ming Wang, *et al*., *Ab initio dynamical mean-field theory with natural orbitals renormalization group impurity solver: Formalism and applications*, [npj Comput Mater 11, 86 (2025)](https://www.nature.com/articles/s41524-025-01586-6).

3. Zhenfeng Ouyang, *et al*., *DFT+DMFT study of correlated electronic structure in the monolayer-trilayer phase of La₃Ni₂O₇*, [Phys. Rev. B, **111**(12), 125111 (2025)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.111.125111).

4. Yin Chen, *et al*., *Non-fermi liquid and antiferromagnetic correlations with hole doping in the bilayer two-orbital hubbard model of La₃Ni₂O₇ at zero temperature*, [Phys. Rev. B, **110**(23), 235119 (2024)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.110.235119).

5. Yi-Heng Tian, *et al*., *Correlation effects and concomitant two-orbital s±-wave superconductivity in La₃Ni₂O₇ under high pressure*, [Phys. Rev. B, **109**(16), 165154 (2024)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.109.165154).



-------------------------------------------------------------------------------------------
## Recent Updates:

Version: v1.10.00 @ 2025.04.27
    realized: Implemented parameter reading from the PARAMS.norg file, with the option to choose between a single NORG iteration or standard Bethe lattice multi-orbital DMFT calculation.

Version: v1.10.01 @ 2025.04.30
    realized: Successfully tested controlling the impurity solver via the PARAMS.norg file, allowing the user to choose between a single NORG iteration or a standard Bethe lattice multi-orbital DMFT calculation.
    note: To enhance computational stability, the program implements two types of degeneracy determination: one is the actual degeneracy, and the other is a user-specified preset degeneracy. During calculations, the program always uses the larger of these two values to determine the degenerate ground state basis, but afterwards, the preset degeneracy (p.degel) is updated according to the actual degeneracy. In DMFT mode, once a degenerate ground state is detected during the iteration process, the degeneracy count can only increase and will not decrease.
