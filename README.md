# iNORG - First fock version

Now this package is currently under developement. 

## Version

v1.9.16

## Installation

To install iNORG, follow these steps:

1. Clone the repository:
```bash
git clone https://github.com/jmw-phys/iNORG.git
cd iNORG
```

2. Compile the code:
```bash
make clean
make
```

3. Run the program:
```bash
mpirun -n [number_of_processes] ./testing/inorg
```

Note: 
- The code requires MPI and MKL libraries
- Make sure you have the necessary compilers and libraries installed
- The number of processes should be adjusted based on your system configuration


## License

GNU AFFERO General Public License Version 3

## Documentation

The reference guide is placed in the `docs` folder.

-------------------------------------------------------------------------------------------
recent big up data:

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