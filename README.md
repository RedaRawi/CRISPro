### CRISPro
An Automated Pipeline for Protein Conformation Stabilization by Proline

## Motivation
Recent studies have shown that the yield, antigenicity, and immunogenicity of an immunogen can be enhanced by stabilizing it into a specific conformation. Such stabilization often involves the engineering of proline mutations at residue positions where a proline is structurally compatible with the target conformation, but not with the alternative conformation. However, there is no publicly available tool that can design proline mutations for this purpose automatically. Here, we implemented an automated tool, CRISP, that inputs structural coordinates of a desirable conformation and/or an alternative conformation, and outputs a list of residue positions where proline mutations are predicted to stabilize the target conformation based on compatibility of phi-psi angles, secondary structure, and steric clashes. Thus, CRISP can be used to engineer immunogens into specific conformation and to design serologic probes, capable of isolating antibodies that recognize a target shape.

## Installation

# Requirements
- R (https://www.r-project.org)
  - R libraries
    - bio3d
    - foreach
    - doParallel
- PyMOL (version 1.7.2.1) (HTML LINK
- DSSP command-line version (http://swift.cmbi.ru.nl/gv/dssp/)

# Run 
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
  9. PDB file of desired conformation (https://www.rcsb.org/pdb/static.do?p=file_formats/pdb/index.html)
  10. Chain within PDB file of desired conformation
  11. PDB file of alternative conformation (https://www.rcsb.org/pdb/static.do?p=file_formats/pdb/index.html)
  12. Chain within PDB file of alternative conformation
          
  
# Execute in the command line

