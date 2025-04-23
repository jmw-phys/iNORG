# iNORG - First version

Now this package is currently under heavy development. 

## Version

v1.9.16

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


## License

GNU AFFERO General Public License Version 3

## Documentation

A tutorial is planned to be included in the `docs` folder.

## Research Publications with iNORG

1. Jia-Ming Wang, *et al*., *Low-energy inter-band kondo bound states in orbital-selective mott phases*, [Phys. Rev. B, **111**, 155107 (2025)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.111.155107).

2. Jia-Ming Wang, *et al*., *Ab initio dynamical mean-field theory with natural orbitals renormalization group impurity solver: Formalism and applications*, [npj Comput Mater 11, 86 (2025)](https://www.nature.com/articles/s41524-025-01586-6).

3. Zhenfeng Ouyang, *et al*., *DFT+DMFT study of correlated electronic structure in the monolayer-trilayer phase of La₃Ni₂O₇*, [Phys. Rev. B, **111**(12), 125111 (2025)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.111.125111).

4. Yin Chen, *et al*., *Non-fermi liquid and antiferromagnetic correlations with hole doping in the bilayer two-orbital hubbard model of La₃Ni₂O₇ at zero temperature*, [Phys. Rev. B, **110**(23), 235119 (2024)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.109.165154).

5. Yi-Heng Tian, *et al*., *Correlation effects and concomitant two-orbital s±-wave superconductivity in La₃Ni₂O₇ under high pressure*, [Phys. Rev. B, **109**(16), 165154 (2024)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.109.165154).



-------------------------------------------------------------------------------------------
## Recent Updates:

Version: v1.9.00.p3 @ 2024.10.29
    realized: Fix the bug for the 2 or 3 orbital cases

Version: v1.9.11.p3 @ 2024.11.12
    realized: The position that needs to be fitted is controlled by passing parameters through eDMFT.

Version: v1.9.13.p3 @ 2024.11.28
    realized: Different impurity orbitals can now use different numbers of bath sites during NORG process.

Version: v1.9.14.p3 @ 2024.12.16
    realized: Implemented orthogonalization of all orbitals during Green's function calculation.

Version: v1.9.15.p3 @ 2025.03.26
    realized: Improved Lanczos stability and added the functionality to specify the number of degeneracies.

Version: v1.9.16 @ 2025.04.20
    realized: Beging begin the example code for the CPC paper.