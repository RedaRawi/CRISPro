### CRISPro
An Automated Pipeline for Protein Conformation Stabilization by Proline

## Motivation
Recent studies have shown that the yield, antigenicity, and immunogenicity of an immunogen can be enhanced by stabilizing it into a specific conformation. Such stabilization often involves the engineering of proline mutations at residue positions where a proline is structurally compatible with the target conformation, but not with the alternative conformation. However, there is no pub-licly available tool that can design proline mutations for this purpose automatically. Here we im-plemented an automated tool, CRISP, that inputs structural coordinates of a desirable confor-mation and/or an alternative conformation, and outputs a list of residue positions where proline mutations are predicted to stabilize the target conformation based on compatibility of phi-psi an-gles, secondary structure, and steric clashes. Thus, CRISP can be used to engineer immuno-gens into specific conformation and to design serologic probes, capable of isolating antibodies that recognize a target shape.

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
  1.  PDB file of desired conformation (https://www.rcsb.org/pdb/static.do?p=file_formats/pdb/index.html)
  2.  Chain within PDB file of desired conformation
  3.  PDB file of alternative conformation (https://www.rcsb.org/pdb/static.do?p=file_formats/pdb/index.html)
  4.  Chain within PDB file of alternative conformation
  5.  Output prefix, for instance "OUTPUT"
  6.  Path to file containing trans-proline phi-psi angles
  7.  Path to PyMOL script that mutates positions
  8.  Overlap cutoff, if not provided default of 1.5 is used
  9.  Distance cutoff, if not provided default of 10 is used
  10. Number of cores, if not provided default of 1 is used
  
# Execute in the command line

