import pyrosetta
pyrosetta.init()

# Load a protein structure from a local PDB file
pose = pyrosetta.pose_from_pdb("your_protein.pdb")

# Get a standard Rosetta energy function
scorefxn = pyrosetta.get_fa_scorefxn()

# Calculate and print total energy
total_energy = scorefxn(pose)
print(f"Total Rosetta Energy: {total_energy} REU")
