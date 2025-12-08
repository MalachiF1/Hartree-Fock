# Hartree-Fock Program

![Language](https://img.shields.io/badge/language-C%2B%2B-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)

This is a Hartree-Fock program implemented in C++. While written primarily for pedagogical purposes, the code is optimized for performance using integral screening methods, an efficient integral-engine algorthim, modern C++ techniques and OpenMP parallelization (However, there is still much room for improvement).

## Features

- **Methods:**
  - Restricted Hartree-Fock (RHF) for closed-shell molecules.
  - Unrestricted Hartree-Fock (UHF) for open-shell molecules.
- **Integral Engine:**
  - McMurchie-Davidson scheme for one- and two-electron integrals.
  - Schwartz screening and density screening for performance.
  - Evaluation of Boys function using an efficient interpolation-recurrence algorithm.
  - Support for Cartesian Gaussian-type orbitals.
- **SCF Convergence:**
  - Direct Inversion in the Iterative Subspace (DIIS).
  - Damping and level shifting.
  - Direct SCF method (recomputing integrals on the fly) to save memory.
- **Symmetry:**
  - Automatic point group detection.
  - Automatic reorientation to canonical axes.

## Building the Program

### Dependencies:

- Make
- GCC (version > 10 recommended)
- OpenBLAS (Optional, but recommended)

### Installation (Ubuntu/Debian):

```bash
sudo apt install make gcc libopenblas-dev
```

### Build Steps

```bash
git clone --recurse-submodules https://github.com/MalachiF1/Hartree-Fock.git
cd Hartree-Fock
make -j $nproc
```

## How to Use

Run the program by providing an input file and an output destination.

```Bash
./bin/Hartree-Fock <input file> <output file>
```

### Input File Example

The input format is inspired by Q-Chem. It is case-insensitive, and comments begin with #.

```Plaintext
$molecule    # H2O molecule
    0    1
    H    1.4327076887   -0.9610871670    0.0000000000
    H   -1.4327076887   -0.9610871670    0.0000000000
    O   -0.0000000000    0.2403536480   -0.0000000000
$end

$scf
    UNRESTRICTED       TRUE
    GUESS_MIX          2      # break symmetry
    MAX_ITER           50
    ENERGY_TOL         1.0e-8
    DIRECT             TRUE
    DENSITY_THRESH     1.0e-8
    PRINT_FULL_MOS     true
$end

$basis
    6-31G
$end
```

## Input Documentation

Exact whitespace and linebreaks are not enforced, and the input may work well without the recommended structure. However, in order to avoid issues it is best to work by the following guidelines.

The input is seperated into three blocks: the `$molecule` block, `$scf` block and `$basis` block. The order of the blocks does not matter.

### Molecule Block (`$molecule`)

Contains the charge, spin multiplicity and geometry.

- **Line 1:** `Charge` `Multiplicity`
- **Subsequent Lines:** `Atom` `X` `Y` `Z` (Coordinates in **Bohr**, the atom may be specified by atom number or atomic symbol).

### Basis Block (`$basis`)

Specifies the basis set. Currently supports standard basis sets avialable in JSON format from the [Basis Set Exchange](www.basissetexchange.org).

- Example: `STO-3G`, `6-31G**`, `cc-pVDZ`

In the future, you will be able to add custom basis sets, but this is not yet implemented.
If a basis set is not available and it exists in the Basis Set Exchange, you may download (in JSON format) and add it to the `basis_sets` directory. It will then be available for use.

### SCF Block (`$SCF`)

Controls the parameters of the calculation
|Keyword|Type|Default|Description|
|:---|:---:|:---:|:---|
|**MAX_ITER**|`int`|`50`|Maximum number of SCF iterations.|
|**ENERGY_TOL**|`float`|`1e-6`|Convergence threshold for total energy.|
|**DENSITY_TOL**|`float`|`1e-6`|Convergence threshold for the density matrix (RMSD).|
|**DIIS**|`bool`|`true`|Enable DIIS acceleration.|
|**DIIS_MAX_SIZE**|`int`|`6`|Max error vectors (and previous Fock matrices) stored for DIIS subspace.|
|**DIIS_ERROR_TOL**|`float`|`1e-6`|DIIS error norm convergence threshold.|
|**DIRECT**|`bool`|`false`|Enable Direct SCF (intergals computed on-the-fly)|
|**SCHWARTZ_THRESH**|`float`|`1e-10`|Threshold for Schwartz screening of integrals|
|**DENSITY_THRESH**|`float`|`1e-10`|Threshold for density screning (Direct SCF only).|
|**SYMMETRY**|`bool`|`true`|Detect point group and reorient molecule.|
|**SYMMETRY_TOL**|`float`|`1e-5`|Tolerance for symmetry detection. It is recommended to not decrease this value, and lowering of the tolerance may be needed based on the user input|
|**PRINT_FULL_MOS**|`bool`|`false`|Print all MO coefficients to output.|
|**UNRESTRICTED**|`bool`|`false`.| Perform UHF calculation (default `true` for open shell systems).|
|**GUESS_MIX**|`int`|`0`|Mix the $\alpha$ HOMO/LUMO to break $\alpha$-$\beta$ symmetry (0-10). The new $\alpha$ HOMO is given by: $\frac{1}{\sqrt{1+k^2}}\left(\mathrm{HOMO}+k\cdot\mathrm{LUMO}\right)$, and the $\alpha$ LUMO is given by: $\frac{1}{\sqrt{1+k^2}}\left(\mathrm{LUMO}-k\cdot\mathrm{HOMO}\right)$, where $k =\frac{\mathrm{GUESS\\_MIX}}{10}$.|
|**DAMP**|`int`|`0`| Static density damping percentage (0-100). The new density matrix is given by $D\_{damped} = \left(1 - \frac{\mathrm{DAMP}}{100}\right) \cdot D_{n} + \frac{\mathrm{DAMP}}{100} \cdot D_{n-1}$.|
|**MAX_DAMP_CYCLES**|`int`|`30`|Iteration limit for damping.|
|**STOP_DAMP_THRESH**|`float`|`0`|Disables damping when $\Delta E$ is less than **STOP_DAMP_THRESH**. If set to `0`, damping will not be disabled by energy convergence.|
|**LEVEL_SHIFT**|`float`|`0.0`|Shift virtual orbital energies (in Hartree).|
|**LSHIFT_GAP_TOL**|`float`|`infinity`|Only applies level-shifting if the HOMO-LUMO gap is less than **LSHIFT_GAP_TOL**.|
|**MAX_LSHIFT_CYCLES**|`int`|`30`|Maximum iterations of level-shifting.|
|**STOP_LSHIFT_THRESH**|`double`|`0`|Disables level-shifting when $\Delta E$ is less than **STOP_LSHIFT_THRESH**. If set to `0`, level-shifting will not be disabled by energy convergence.|

## Implementation Details

### Integral Evaluation

The integral engine is based on the **McMurchie-Davidson scheme** [[1]](#1).
Two electron integrals are calculated a shell-quartet at a time.
A "shell" is defined as a group of basis functions sharing the same exponent, center, and total angular momentum.
This structure allows many intermediate quantities to be reused independent of the specific angular momentum distribution (e.g., $d_{xy}$ vs $d_{z^2}$).

**Note on Future Work:** I am currently working on migrating to Gill's and Pople's **PRISM algorithm** [[2]](#2)-[[4]](#4) on the PRISM branch.
PRISM dynamically chooses the optimal contraction/transformation path for each batch of shell-quartets to maximize efficiency. Shell-quartets should be calculated by batches of similar quartets at a time (i.e. have the same angular momenta and degrees of contraction) in order to maximize utilization of vectorization (SIMD).

The integral workload is divided into muliple parallel threads using OpemMP.
For the initial $\left[0\right]_{m}$ integrals, the Boys function is evaluated using an efficient interpolation-recurrence algorithm based on Ref [[5]](#5).

### Integral Screening

There are $N^{4}$ possible shell-quartets, where $N$ is the number of basis functions.
We exploit the 8-fold permuational symmetry of the two-electron integrals:

$$\langle ab | cd \rangle = \langle ba | cd \rangle = \langle ab | dc \rangle = \langle ba | dc \rangle = \langle cd | ab \rangle = \langle dc | ab \rangle = \langle cd | ba \rangle = \langle dc | ba \rangle$$

This reduces the unique integrals to $\frac{n^2(n+1)^2}{8}$ (roughly $\frac{1}{8}N^{4}$).
To further reduce cost, we apply **Schwartz Screening** [[6]](#6). An upper bound is established for each quartet:

$$|\langle ab | cd \rangle| < \sqrt{\langle ab | ab \rangle \langle cd | cd \rangle}$$

If the bound is below `SCWARTZ_THRESH`, the quartet is skipped.
The $\langle ij|ij \rangle$ terms can be calculated in comparatively negligible $O(n^2)$ time. This is a significant speedup for large systems, where the majority of integrals are small and can be skipped.

In Direct SCF, we also use **Density Screening** [[7]](#7).
Integrals are skipped if their contribution to the Fock matrix would be negligible based on the magnitude of the corresponing density elements.
The Fock matrix is built incrementally from the previous iteration using the difference density matrix $\Delta D = D_{n} - D_{n-1}$.
This allows for more efficient density screening, especially in later iterations where the density matrix does not change significantly.

### DIIS Convergence Acceleration

Direct Inversion in the Iterative Subspace (DIIS) extrapolates the Fock matrix from previous iterations to minimize an error residual [[8]](#8).
This is achieved by obtaining a set of optimal coefficients in each iteration, given by the solution of a set of linear equations.
The error is defined as the commutator of the Fock and Density matrices (which is zero for the SCF solution):

$$e = \mathbf{SDF} - \mathbf{FDS}$$

For unrestricted calculations, the alpha and beta errors are combined: $e = e_{\alpha} + e_{\beta}$.

### Damping and Level-Shifting

To handle oscillating or non-converging systems:

1. **Level-Shifting:** Adds a constant shift to the virtual block of the Fock matrix. This increases the HOMO-LUMO gap and prevents oscillations due to switching of the HOMO and LUMO between iterations.
2. **Damping:** Mixes the current density matrix with the previous one: $D_{new} = (1-\alpha)D_{n} + \alpha D_{n-1}$. This reduces the effective step size [[10]](#10).

### Symmetry Detection

The program detects the molecular point group by classifying atoms into sets of Symmetrically Equivalent Atoms (SEAs) by computing the interatomic distance matrix, and with taking into acount the atomic numbers, finding atoms with equivalent distance profiles.
The SEAs are used to generate a small set of possible symmetry elements, which are tested by applying the corresponding symmetry operations on the molecular geometry. If the transformed geometry is a permutation of the original geometry (up to a tolerance), the symmetry operation is valid. A decision tree algorithm utilizes the information from the SEAs and principal moments of inertia to determine the point group and principal axes.[[11]](#11)

Once identified, the molecule is rotated into a canonical orientation (Primary Axis $\rightarrow z$, Secondary axis $\rightarrow x$) following the conventions of the MolSym package. [[12]](#12)

_Note: Currently, symmetry is used for point group detection and geometry oreintation only. Use of Symmetry adapted Linear Combinations (SALCs) for block-diagonalizing the Fock matrix and integral screening (two electron-integrals are only non-zero if the direct product of the irreducible representations of the four orbitals contains the totally symmetric irreducible representation of the molecule's point group) are planned for future releases._

## References

<a id="1">[1]</a>
McMurchie, L.; Davidson, E. R. One- and Two-Electron Integrals over Cartesian Gaussian Functions. Journal of Computational Physics 1978, 26 (2), 218–231. https://doi.org/10.1016/0021-9991(78)90092-x.

<a id="2">[2]</a>
Peter; Pople, J. A. The Prism Algorithm for Two‐Electron Integrals. International Journal of Quantum Chemistry 1991, 40 (6), 753–772. https://doi.org/10.1002/qua.560400605.

<a id="2">[3]</a>
Gill, P. M. W. Molecular Integrals over Gaussian Basis Functions; Academic Press, 2008; Vol. 25, pp. 141–205. https://doi.org/10.1016/S0065-3276(08)60019-2.

<a id="4">[4]</a>
Gill, P.; Head-Gordon, M.; Pople, J. A. Efficient Computation of Two-Electron - Repulsion Integrals and Their Nth-Order Derivatives Using Contracted Gaussian Basis Sets. 1990, 94 (14), 5564–5572. https://doi.org/10.1021/j100377a031.

<a id="5">[5]</a>
Alexander; Ochsenfeld, C. A Rigorous and Optimized Strategy for the Evaluation of the Boys Function Kernel in Molecular Electronic Structure Theory. Journal of Computational Chemistry 2015, 36 (18), 1390–1398. https://doi.org/10.1002/jcc.23935.

<a id="6">[6]</a>
Horn, H.; Weiß, H.; Háser, M.; Ehrig, M.; Ahlrichs, R. Prescreening of Two-Electron Integral Derivatives in SCF Gradient and Hessian Calculations. Journal of Computational Chemistry 1991, 12 (9), 1058–1064. https://doi.org/10.1002/jcc.540120903.

<a id="7">[7]</a>
Häser, M.; Ahlrichs, R. Improvements on the Direct SCF Method. Journal of Computational Chemistry 1989, 10 (1), 104–111. https://doi.org/10.1002/jcc.540100111.

<a id="8">[8]</a>
Pulay, P. Convergence Acceleration of Iterative Sequences. The Case of Scf Iteration. Chemical Physics Letters 1980, 73 (2), 393–398. https://doi.org/10.1016/0009-2614(80)80396-4.

<a id="9">[9]</a>
Saunders, V. R.; Hillier, I. H. A “Level–Shifting” Method for Converging Closed Shell Hartree–Fock Wave Functions. International Journal of Quantum Chemistry 1973, 7 (4), 699–705. https://doi.org/10.1002/qua.560070407.

<a id="10">[10]</a>
Zerner, M. C.; Hehenberger, M. A Dynamical Damping Scheme for Converging Molecular Scf Calculations. Chemical Physics Letters 2002, 62 (3), 550–554. https://doi.org/10.1016/0009-2614(79)80761-7.

<a id="11">[11]</a>
Beruski, O.; Vidal, L. N. Algorithms for Computer Detection of Symmetry Elements in Molecular Systems. Journal of Computational Chemistry 2013, 35 (4), 290–299. https://doi.org/10.1002/jcc.23493.

<a id="12">[12]</a>
Goodlett, S. M.; Kitzmiller, N. L.; Turney, J. M.; Schaefer, H. F. MolSym: A Python Package for Handling Symmetry in Molecular Quantum Chemistry. The Journal of Chemical Physics 2024, 161 (2). https://doi.org/10.1063/5.0216738.
