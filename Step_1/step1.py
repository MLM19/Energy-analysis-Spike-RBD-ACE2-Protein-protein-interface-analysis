from Bio.PDB import *

def get_interface_residues(structure, chain1_id, chain2_id, distance_threshold):
	inter_res = set()
# Iterate through the atoms of each chain
	for model in structure:
		for chain1 in model[chain1_id]:
			for chain2 in model[chain2_id]:
				for atom1 in chain1:
					for atom2 in chain2:
						if atom1.get_parent().id == atom1.get_parent().id and atom2.get_parent().id == atom2.get_parent().id:
							distance = atom1 - atom2
							if distance <= distance_threshold:
								inter_res.add(atom1.get_parent().id[1])
								inter_res.add(atom2.get_parent().id[1])

	return sorted(list(inter_res))



# Load the structure (assuming 6m0j.pdb is in the same directory)
parser = PDBParser(QUIET=True)
structure = parser.get_structure("6m0j", "6m0j_fixed.pdb")

# Step 1: Visual inspection in PyMOL to choose a distance threshold (let's say 5 Å + 1-2 Å = 6-7 Å)
distance_threshold = 4.3  # Change this value based on your visual inspection

# Get interface residues
chain1_id = "A"
chain2_id = "E"
all_inter_res = get_interface_residues(structure, chain1_id, chain2_id, distance_threshold)

# Print the list of interface residues for both chains
print(f"Interface residues : {all_inter_res}")
