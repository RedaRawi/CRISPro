#from pymol import cmd
import sys

pdb = sys.argv[1]
pdb_output = sys.argv[2]
selection = sys.argv[3]
mutant = sys.argv[4]
angstroem_distance = sys.argv[5]

selection_distance = "br. all within " + angstroem_distance + " of " + selection

print ("pdb: " + pdb)
print ("name: " + pdb_output)
print ("selection: " + selection)
print ("mutant: " + mutant)
print ("selection_distance: " + selection_distance)

# Load PDB
print ("Load PDB")
cmd.load(pdb)

# Mutate PDB
print ("Mutate PDB")
cmd.wizard("mutagenesis")
cmd.refresh_wizard()
cmd.get_wizard().do_select(selection)
cmd.get_wizard().set_mode(mutant)
cmd.get_wizard().apply()
cmd.set_wizard()

## Save mutated  structure
#print("{0}".format(pdb_output))
#cmd.save("{0}".format(pdb_output))


# Select all residues with a maximum distance of 10A from mutant residue
print ("Select all residues with a maximum distance of 10A from mutant residue")
cmd.select("selection_distance",selection_distance)


# Save mutated and subselected structure
print ("Save mutated and subselected structure")
print("{0}".format(pdb_output))
cmd.save("{0}".format(pdb_output),selection_distance)
