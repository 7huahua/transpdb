from Bio.PDB import PDBParser, PDBIO, Superimposer
from Bio.Seq import Seq
import numpy as np
import os

def read_pdb(file_path):
    """Read PDB file and return structure object"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', file_path)
    return structure

def calculate_centroid(atoms):
    """Calculate the centroid of a set of atoms"""
    coords = [atom.get_coord() for atom in atoms]
    centroid = np.mean(coords, axis=0)
    return centroid

def rotate_y_180(atoms, centroid):
    """Rotate 180 degrees around the y-axis about the centroid"""
    rotation_matrix = np.array([
        [-1, 0, 0],  # x = -x
        [0, 1, 0],   # y = y
        [0, 0, -1]   # z = -z
    ])
    for atom in atoms:
        coord = atom.get_coord() - centroid
        new_coord = np.dot(rotation_matrix, coord) + centroid
        atom.set_coord(new_coord)

def translate(atoms, distance_x, distance_y, distance_z):
    """Translate atoms by a given distance along the x-axis"""
    for atom in atoms:
        coord = atom.get_coord()
        coord[0] += distance_x
        coord[1] += distance_y
        coord[2] += distance_z
        atom.set_coord(coord)
    
def translate_x(atoms, distance):
    """Translate atoms by a given distance along the x-axis"""
    for atom in atoms:
        coord = atom.get_coord()
        coord[0] += distance
        atom.set_coord(coord)

def rotate_x(atoms, angle, centroid):
    """Rotate atoms by a given angle around the x-axis about the centroid"""
    rad = np.deg2rad(angle)  # Convert angle to radians
    rotation_matrix = np.array([
        [1, 0, 0],
        [0, np.cos(rad), -np.sin(rad)],
        [0, np.sin(rad), np.cos(rad)]
    ])
    for atom in atoms:
        coord = atom.get_coord() - centroid  # Translate to origin based on centroid
        new_coord = np.dot(rotation_matrix, coord) + centroid  # Rotate and translate back
        atom.set_coord(new_coord)

def write_pdb(structure, output_file):
    """Save the structure to a PDB file"""
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file)

def get_residue_atoms(structure, chain_id, start, end):
    """Get all atoms of residues from start to end in the specified chain"""
    atoms = []
    for residue_id in range(start, end + 1):
        residue = structure[0][chain_id][residue_id]
        atoms.extend(residue.get_atoms())
    return atoms

def calculate_rotation_angle(coord1, coord2):
    """Calculate the rotation angle around the x-axis between two coordinates"""
    # Calculate vector
    vector = coord2 - coord1
    
    # Calculate rotation angle in the y-z plane
    angle_rad = np.arctan2(vector[2], vector[1])
    angle_deg = np.degrees(angle_rad)
    return angle_deg


if __name__ == "__main__":
    # folder path
    folder_path = "C:/Users/Han/Downloads/New/New/"  # Replace with your folder path
    # Load PDB file
    long_pdb_path = os.path.join(folder_path, "M13mp18-start48-end48.pdb")  # Replace with your PDB file path
    short_5_3 = os.path.join(folder_path, "TAATAGTAGTAGCATTCCCCCCCC.pdb")  # Replace with your PDB file path
    short_3_5 = os.path.join(folder_path, "ACTAAATAACCTACAAGGGGGGGG.pdb")  # Replace with your PDB file path
    # pdb_path_A = os.path.join(folder_path, "M13mp18-start48-end48.pdb")  # Replace with your PDB file path
    # pdb_path_B = os.path.join(folder_path, "TAATAGTAGTAGCATTCCCCCCCC.pdb")  # Replace with your PDB file path
    output_long_path = os.path.join(folder_path, "output_long.pdb")  # Replace with the output file path
    output_short_5_3_path = os.path.join(folder_path, "output_short_5_3.pdb")  # Replace with the output file path
    output_short_3_5_path = os.path.join(folder_path, "output_short_3_5.pdb")  # Replace with the output file path

    # Load PDB files
    structure_long = read_pdb(long_pdb_path)
    structure_short_5_3 = read_pdb(short_5_3)
    structure_short_3_5 = read_pdb(short_3_5)

    # Select the chain to be rotated, assume it's chain A
    chain_long = structure_long[0]['A']  # Ensure the chain ID is correct
    chain_short_5_3 = structure_short_5_3[0]['A']  # Ensure the chain ID is correct
    chain_short_3_5 = structure_short_3_5[0]['A']  # Ensure the chain ID is correct

    # The starting positions of A and B are the same. We need to calculate the matching sequence position for A and B, then calculate the centroid based on the position
    chain_long_sequence = ''.join([residue.get_resname()[1] for residue in chain_long])
    chain_short_5_3_sequence = ''.join([residue.get_resname()[1] for residue in chain_short_5_3])
    chain_short_3_5_sequence = ''.join([residue.get_resname()[1] for residue in chain_short_3_5])

    # Use Biopython's Seq to get the matching positions
    chain_long_seq = Seq(chain_long_sequence)
    chain_short_5_3_seq = Seq(chain_short_5_3_sequence)
    chain_short_3_5_seq = Seq(chain_short_3_5_sequence)

    # Get the reverse complement of chain A
    chain_long_seq_complement = chain_long_seq.complement()
    chain_short_5_3_seq_reverse = chain_short_5_3_seq[::-1]

    # # 5->3 should shift the short_5_3 sequence to the left
    # # distance = the first x coords of long - the first x coords of short_5_3
    # distance_x = chain_long[1]['P'].get_coord()[0] - chain_short_5_3[1]['P'].get_coord()[0]
    # distance_y = chain_long[1]['P'].get_coord()[1] - chain_short_5_3[1]['P'].get_coord()[1]
    # distance_z = chain_long[1]['P'].get_coord()[2] - chain_short_5_3[1]['P'].get_coord()[2]
    # translate(list(chain_short_5_3.get_atoms()), distance_x, distance_y, distance_z)

    # # then, rotate 180 degrees around the y-axis about the centroid, which is the 8th, 9th residue
    # centroid = calculate_centroid(chain_short_5_3[8].get_list() + chain_short_5_3[9].get_list())
    # rotate_y_180(list(chain_short_5_3.get_atoms()), centroid)

    # # 3->5 should shift the short_3_5 sequence to the right
    # # distance = the last x coords of long - the last x coords of short_3_5
    # distance_x = chain_long[len(chain_long)]['P'].get_coord()[0] - chain_short_3_5[16]['P'].get_coord()[0]
    # distance_y = chain_long[len(chain_long)]['P'].get_coord()[1] - chain_short_3_5[16]['P'].get_coord()[1]
    # distance_z = chain_long[len(chain_long)]['P'].get_coord()[2] - chain_short_3_5[16]['P'].get_coord()[2]
    # translate(list(chain_short_3_5.get_atoms()), distance_x, distance_y, distance_z)
    


    # 5->3 should shift the short_5_3 sequence to the left
    sup = Superimposer()
    fixed_atoms = structure_long[0]['A'][1].get_list()
    moving_atoms = structure_short_5_3[0]['A'][1].get_list()
    print('this is moving_atoms', moving_atoms[0].get_coord())
    sup.set_atoms(fixed_atoms, moving_atoms)
    print('this is rotran',sup.rotran)
    print(sup.rms)
    all_moving_atoms = structure_short_5_3[0]['A'].get_atoms()
    sup.apply(all_moving_atoms)
    print(fixed_atoms[0].get_coord())
    print(structure_short_5_3[0]['A'][1].get_list()[0].get_coord())

    max_match = ""  
    max_match_position = -1  

    for i in range(len(chain_short_5_3_seq), 0, -1):
        subseq = chain_short_5_3_seq[:i]  
        subseq_reverse = subseq[::-1]

        match_position = chain_long_seq_complement.find(subseq_reverse)
        if match_position != -1:
            max_match = subseq_reverse
            max_match_position = match_position
            break  

    # Get the position of the match, which is the offset of the first appearance of chain B
    shift = i

    # Calculate the centroid of the short_5_3 chain
    if shift % 2 == 0:
        centroid = calculate_centroid(chain_short_5_3[shift//2].get_list() + chain_short_5_3[shift//2+1].get_list())
    else:
        centroid = calculate_centroid(chain_short_5_3[shift//2].get_list())

    # Rotate chain B around the y-axis by 180 degrees about the centroid
    rotate_y_180(list(chain_short_5_3.get_atoms()), centroid)


    # 3->5 should shift the short_3_5 sequence to the right
    max_match = ""  
    max_match_position = -1  
    for i in range(len(chain_short_3_5_seq), 0, -1):
        subseq = chain_short_3_5_seq[:i]

        match_position = chain_long_seq_complement.find(subseq)
        if match_position != -1:
            max_match = subseq
            max_match_position = match_position
            break  

    shift = i

    # centroid = calculate_centroid(chain_short_3_5.get_atoms())
    # rotate_y_180(list(chain_short_3_5.get_atoms()), centroid)


    sup = Superimposer()
    fixed_atoms = structure_long[0]['A'][len(chain_long_seq)].get_list()
    moving_atoms = structure_short_3_5[0]['A'][shift].get_list()
    min_len = min(len(fixed_atoms), len(moving_atoms))
    fixed_atoms = fixed_atoms[:min_len]
    moving_atoms = moving_atoms[:min_len]
    print(structure_long[0]['A'][len(chain_long_seq)].get_resname())
    print(structure_short_3_5[0]['A'][shift].get_resname())
    sup.set_atoms(fixed_atoms, moving_atoms)

    sup.apply(structure_short_3_5[0]['A'].get_atoms())

    rotate_x(list(chain_short_3_5.get_atoms()), 180, calculate_centroid(chain_short_3_5.get_atoms()))

    # chain_short_5_3.id = 'B'
    # structure_long[0].add(chain_short_5_3)
    # chain_short_3_5.id = 'C'
    # structure_long[0].add(chain_short_3_5)

    write_pdb(structure_long, output_long_path)
    write_pdb(structure_short_5_3, output_short_5_3_path)
    write_pdb(structure_short_3_5, output_short_3_5_path)

