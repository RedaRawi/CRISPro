# CRISPro
An Automated Pipeline for Protein Conformation Stabilization by Proline

## Motivation
Recent studies have shown that the yield, antigenicity, and immunogenicity of an immunogen can be enhanced by stabilizing it into a specific conformation. Such stabilization often involves the engineering of proline mutations at residue positions where a proline is structurally compatible with the target conformation, but not with the alternative conformation. However, there is no publicly available tool that can design proline mutations for this purpose automatically. Here, we implemented an automated tool, CRISP, that inputs structural coordinates of a desirable conformation and/or an alternative conformation, and outputs a list of residue positions where proline mutations are predicted to stabilize the target conformation based on compatibility of phi-psi angles, secondary structure, and steric clashes. Thus, CRISP can be used to engineer immunogens into specific conformation and to design serologic probes, capable of isolating antibodies that recognize a target shape.

## Installation

### Requirements
- R (https://www.r-project.org)
  - R libraries
    - bio3d
    - foreach
    - doParallel
- PyMOL (version 1.7.2.1) (https://sourceforge.net/projects/pymol/files/pymol/1.7/)
- DSSP command-line version (http://swift.cmbi.ru.nl/gv/dssp/)

## Run 
CRISPro can be run in the command line

Several input arguments are necessary to run the CRISPro.
  1.  Overlap cutoff, 1.5 is recommended
  2.  Distance cutoff, 10 is recommended
  3.  Phi-Psi cutoff, 0.0005 is recommended
  4.  Number of cores
  5.  Path to file containing trans-proline phi-psi angles
  6.  Path to PyMOL script that mutates positions
  7.  Output prefix, for instance "OUTPUT"
  8.  Option (either 1, 2, or 3)
  9.  PDB file of target conformation (https://www.rcsb.org/pdb/static.do?p=file_formats/pdb/index.html) (Option 1 and 2)
  10. Chain within PDB file of desired conformation (Option 1 and 2)
  11. PDB file of alternative conformation (https://www.rcsb.org/pdb/static.do?p=file_formats/pdb/index.html) (Option 1 and 3)
  12. Chain within PDB file of alternative conformation (Option 1 and 3)


### Execute in the command line
R --vanilla < CRISPro_v.1.1-1.R 1.5 10 0.0005 8 rama8000-transpro.data mutate.py OUTPUT 1 target.pdb A alternative.pdb A


### Result
- Results will be saved in OUTPUT_DF.csv
- Explanation of Results
  - Option 1 columns include:
    - Residue number
    - Chain of target conformation
    - Amino acid of target conformation
    - Phi angle of target conformation
    - Psi angle of target conformation
    - Trans proline angle of target conformation
    - Helix of target conformation
    - DSSP secondary structure prediction of target conformation
    - Number of clashes of target conformation
    - Chain of alternative conformation
    - Amino acid of alternative conformation
    - Phi angle of alternative conformation
    - Psi angle of alternative conformation
    - Trans proline angle of alternative conformation
    - Helix of alternative conformation
    - DSSP secondary structure prediction of alternative conformation
    - Number of clashes of alternative conformation
    - Compatible with target conformation
    - Compatible with alternative conformation
    - Compatible with target and not with alternative conformation
  - Option 2 columns include:
    - Chain of target conformation
    - Residue number of target conformation
    - Amino acid of target conformation
    - Phi angle of target conformation
    - Psi angle of target conformation
    - Trans proline angle of target conformation
    - Helix of target conformation
    - DSSP secondary structure prediction of target conformation
    - Number of clashes of target conformation
    - Compatible with target conformation
  - Option 3 columns include: 
    - Chain of alternative conformation
    - Residue number of alternative conformation
    - Amino acid of alternative conformation
    - Phi angle of alternative conformation
    - Psi angle of alternative conformation
    - Trans proline angle of alternative conformation
    - Helix of alternative conformation
    - DSSP secondary structure prediction of alternative conformation
    - Number of clashes of alternative conformation
    - Compatible with alternative conformation


### Example/Test

#### RSV F example
Copy all required files into a directory (Files: CRISPro_v.1.1-1.R, rama8000-transpro.data, mutate.py, 4jhw_trimer.pdb, and 3rrr_trimer.pdb)

- Option 1: R --vanilla < CRISPro_v.1.1-1.R 1.5 10 0.0005 8 rama8000-transpro.data mutate.py RSV_F_option-1 1 4jhw_trimer.pdb A 3rrr_trimer.pdb A
- Option 2: R --vanilla < CRISPro_v.1.1-1.R 1.5 10 0.0005 8 rama8000-transpro.data mutate.py RSV_F_option-2 2 4jhw_trimer.pdb A 
- Option 3: R --vanilla < CRISPro_v.1.1-1.R 1.5 10 0.0005 8 rama8000-transpro.data mutate.py RSV_F_option-3 3 3rrr_trimer.pdb A
