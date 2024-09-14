from Bio.PDB import PDBParser, PDBIO, Superimposer
from Bio.Seq import Seq
from Bio.PDB.Residue import Residue
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

class RoTran:
    def __init__(self, rotation_matrix = None, translation_vector = None, compute=False,chain=None):
        if rotation_matrix:
            self.rotation_matrix = rotation_matrix
        if translation_vector:
            self.translation_vector = translation_vector
        elif compute == False:
            self.rotation_matrix = np.array(
                [[ 9.99999993e-01,  2.15569388e-06, -2.86627892e-06],
                [-1.66738053e-07,  8.26238644e-01,  5.63320239e-01],
                [ 3.58331382e-06, -5.63320239e-01,  8.26238644e-01]]
            )
            # [ 3.38000009e+00 -5.65847633e-04  1.09018337e-03]
            self.translation_vector = np.array([ 3.38000009e+00, -5.65847633e-04,  1.09018337e-03])
        else:
            self.get_chain_rotran(chain)

    def get_chain_rotran(self, chain_long):
        group_atom_names = ["P", "O1P", "O2P", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'"]

        # for every residue in the long chain, get the atoms of the group_atom_names
        groups = []
        for residue in chain_long:
            group = [atom for atom in residue.get_list() if atom.get_name() in group_atom_names]
            groups.append(group)

        assert len(groups) == len(chain_long)  

        # get the average rotran
        rots, trans = [], []
        sup = Superimposer()
        for i in range(len(groups) - 1):
            moving_atoms = groups[i]
            fixed_atoms = groups[i + 1]
            sup.set_atoms(fixed_atoms, moving_atoms)
            rot,tran = sup.rotran
            rots.append(rot)
            trans.append(tran)

        # get the average rotran
        self.rotation_matrix = np.mean(rots, axis=0)
        self.translation_vector = np.mean(trans, axis=0)
        print(self.rotation_matrix)
        print(self.translation_vector)

def rotran_coords(atoms, rotran):
    moving_coord = np.zeros((len(atoms), 3))
    for j in range(len(atoms)):
        moving_coord[j] = atoms[j].get_coord()
    moving_coord = np.dot(moving_coord, rotran.rotation_matrix) + rotran.translation_vector
    for j in range(len(atoms)):
        atoms[j].set_coord(moving_coord[j])

def create_strcture(dna_seq, rotran):

    folder_path = "D:/Han/projects/transpdb"
    A_pdb_path = os.path.join(folder_path, "A.pdb")  # Replace with your PDB file path
    T_pdb_path = os.path.join(folder_path, "T.pdb")  # Replace with your PDB file path
    C_pdb_path = os.path.join(folder_path, "C.pdb")  # Replace with your PDB file path
    G_pdb_path = os.path.join(folder_path, "G.pdb")  # Replace with your PDB file path

    # Load PDB files
    structure_A = read_pdb(A_pdb_path)
    structure_T = read_pdb(T_pdb_path)
    structure_C = read_pdb(C_pdb_path)
    structure_G = read_pdb(G_pdb_path)

    res_A = structure_A[0]['A'][1]
    res_T = structure_T[0]['A'][1]
    res_C = structure_C[0]['A'][1]
    res_G = structure_G[0]['A'][1]
    # create a new structure
    if dna_seq[0] == 'A':
        structure = structure_A
    elif dna_seq[0] == 'T':
        structure = structure_T
    elif dna_seq[0] == 'C':
        structure = structure_C
    elif dna_seq[0] == 'G':
        structure = structure_G

    last_A, last_T, last_C, last_G = res_A.copy(), res_T.copy(), res_C.copy(), res_G.copy()

    # compute the pdb file
    for i,res_name in enumerate(dna_seq):
        if i == 0:
            continue

        if i!=0:
            if res_name == 'A':
                new_res = last_A.copy()
            elif res_name == 'T':
                new_res = last_T.copy()
            elif res_name == 'C':
                new_res = last_C.copy()
            elif res_name == 'G':
                new_res = last_G.copy()
                
            rotran_coords(new_res.get_list(), rotran)

            # update the last residue
            for r in [last_A, last_T, last_C, last_G]:
                rotran_coords(r.get_list(), rotran)

        # create a new residue
        res_id = i+1
        res = Residue((' ', res_id, ' '), 'D'+res_name, res_id)

        # add the atoms to the new residue
        if i!= len(dna_seq)-1:
            for atom in new_res.get_list():
                if atom.get_name() not in ['HTER','OXT','HCAP']:
                    # new_atom = Atom(atom.get_name(), atom.get_coord(), atom.get_bfactor(), atom.get_occupancy(), atom.get_altloc(), atom.get_fullname(), atom.get_serial_number())
                    res.add(atom)
        else:
            for atom in new_res.get_list():
                if atom.get_name() not in ['HTER','OXT']:
                    res.add(atom)

        # add the residue to the structure
        structure[0]['A'].add(res)

    return structure

if __name__ == "__main__":
    # folder path
    folder_path = "C:/Users/Han/Desktop"  # Replace with your folder path
    # Load PDB file
    long_pdb_path = os.path.join(folder_path, "M13mp18-start48-end48.pdb")  # Replace with your PDB file path
    structure_long = read_pdb(long_pdb_path)

    # "D:\Han\projects\transpdb\A.pdb"


    # Select the chain to be rotated, assume it's chain A
    chain_long = structure_long[0]['A']  # Ensure the chain ID is correct


    # The starting positions of A and B are the same. We need to calculate the matching sequence position for A and B, then calculate the centroid based on the position
    chain_long_sequence = ''.join([residue.get_resname()[1] for residue in chain_long])

    # Use Biopython's Seq to get the matching positions
    chain_long_seq = Seq(chain_long_sequence)

    # Get the reverse complement of chain A
    chain_long_seq_complement = chain_long_seq.complement()

    # compute a pdb file from the sequence, using the rotation matrix and translation vector
    rotran = RoTran(compute=True, chain=chain_long)
    dna_seq = chain_long_sequence
    dna_seq = dna_seq.replace('U', 'T')

    pdb_seq = chain_long_seq

    straight_chain = create_strcture(dna_seq, rotran)

    # write_pdb(structure, os.path.join(folder_path, "new_generate.pdb"))  # Replace with the output file path

    n = len(dna_seq)  # 核苷酸数量
    theta_deg = 360 / n  # 每步旋转的角度，单位为度
    theta_rad = np.radians(theta_deg)  # 转换为弧度

    # 旋转矩阵
    R = np.array([
        [np.cos(theta_rad), -np.sin(theta_rad), 0],
        [np.sin(theta_rad), np.cos(theta_rad), 0],
        [0, 0, 1]
    ])

    # 计算平移距离
    # structure 的最大x坐标-最小x坐标即为一周的周长
    atoms = get_residue_atoms(straight_chain, 'A', 1, n)
    coords = np.array([atom.get_coord() for atom in atoms])
    max_x = np.max(coords[:, 0])
    min_x = np.min(coords[:, 0])

    R_circle = max_x - min_x
    d = (2 * np.pi * R_circle) / n
    T = np.array([
        d * np.cos(theta_rad / 2),
        d * np.sin(theta_rad / 2),
        0
    ])

    print("旋转矩阵 R:")
    print(R)
    print("平移向量 T:")
    print(T)









