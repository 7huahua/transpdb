from Bio.PDB import PDBParser, PDBIO, Superimposer
from Bio.Seq import Seq
from Bio.PDB.Residue import Residue
import numpy as np
from datetime import datetime
import os
import argparse


def read_pdb(file_path):
    """Read PDB file and return structure object"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', file_path)
    return structure

def write_pdb(structure, output_file):
    """Save the structure to a PDB file"""
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file)


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

def generate_chain(dna_seq, structure, chain_name, rotran, res_A, res_T, res_C, res_G):

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
        structure[0][chain_name].add(res)

    return structure

def init(folder_path,dna_seq, chain_name):

    rotran = RoTran()
    # dna_seq = chain_long_sequence
    # convert atcg to ATCG
    dna_seq = dna_seq.upper()
    dna_seq = dna_seq.replace('U', 'T')
    # folder_path = "D:/Han/projects/transpdb/init_pdb"
    # A_pdb_path = os.path.join(folder_path, "A.pdb")  # Replace with your PDB file path
    # T_pdb_path = os.path.join(folder_path, "T.pdb")  # Replace with your PDB file path
    # C_pdb_path = os.path.join(folder_path, "C.pdb")  # Replace with your PDB file path
    # G_pdb_path = os.path.join(folder_path, "G.pdb")  # Replace with your PDB file path

    if chain_name == 'A':
        # Load PDB files
        structure_A_A = read_pdb(os.path.join(folder_path, "A.pdb"))
        structure_A_T = read_pdb(os.path.join(folder_path, "T.pdb"))
        structure_A_C = read_pdb(os.path.join(folder_path, "C.pdb"))
        structure_A_G = read_pdb(os.path.join(folder_path, "G.pdb"))

        chain_A_res_A = structure_A_A[0]['A'][1]
        chain_A_res_T = structure_A_T[0]['A'][1]
        chain_A_res_C = structure_A_C[0]['A'][1]
        chain_A_res_G = structure_A_G[0]['A'][1]

        if dna_seq[0] == 'A':
            structure = structure_A_A
        elif dna_seq[0] == 'T':
            structure = structure_A_T
        elif dna_seq[0] == 'C':
            structure = structure_A_C
        elif dna_seq[0] == 'G':
            structure = structure_A_G

        new_structure = generate_chain(dna_seq, structure, 'A', rotran, chain_A_res_A, chain_A_res_T, chain_A_res_C, chain_A_res_G)


    elif chain_name == 'B':

        structure_B_A = read_pdb(os.path.join(folder_path, "TA.pdb"))
        structure_B_T = read_pdb(os.path.join(folder_path, "AT.pdb"))
        structure_B_C = read_pdb(os.path.join(folder_path, "GC.pdb"))
        structure_B_G = read_pdb(os.path.join(folder_path, "CG.pdb"))

        # Select the chain to be rotated, assume it's chain A
        # chain_long = structure_long[0]['A']  # Ensure the chain ID is correct

        chain_B_res_A = structure_B_A[0]['B'][1]
        chain_B_res_T = structure_B_T[0]['B'][1]
        chain_B_res_C = structure_B_C[0]['B'][1]
        chain_B_res_G = structure_B_G[0]['B'][1]
        
        if reverse_dna_seq[0] == 'A':
            structure = structure_B_A
        elif reverse_dna_seq[0] == 'T':
            structure = structure_B_T
        elif reverse_dna_seq[0] == 'C':
            structure = structure_B_C
        elif reverse_dna_seq[0] == 'G':
            structure = structure_B_G
    
        new_structure = generate_chain(reverse_dna_seq, structure, 'A', rotran, chain_B_res_A, chain_B_res_T, chain_B_res_C, chain_B_res_G)

    return new_structure

def write_config(config_path, chain_A_path, chain_B_path):
    header = '''
# Auto Generated Config
#
# VR Model Specification file <Check SimulationInput.exe for further explanation>
#
# Property     | Example            | Explination
#--------------------------------------------------------------------------------------------
# P (path)       "P path"           Currently a max of 32 files, ask if you need more
#
# E (path)       "E path"           Path to the EEM text file
#
# T (type)       "T pdb", "T obj"   PDB: Molecules, OBJ: 3D graphics, PDBC: Custom PDB with charge.
#
# C (count)      "C n"              Number of L/R should match or exceed the count specified
#
# F (flag)       "F 0"              Bit flag: 1 = sticky grid, 2 = sticky phosphorus, 
#                                               4 = don't center, 8 = ball & stick, 16 = ?
#
# B(bind)        "B 1"              Bind the chains together within the PDB file. 0 to not bind.
#
# S (Scale)      "S 1.1"            Some PDB files have tight atom placing and the structure is sensitive
#
# L (location)   "L x y z"          Location is an offset from the center of the simulation space.
#                                      -x is left, -y is down, -z is closer to the camera 
#
# R (rotation )  "R 90 0 330"       x-axis: left to right, y-axis: vertical, z-axis: depth
#                                     x: top rotates back, y: front rotates right, z: top rotates right
#
# *Anything line that does not start with P, T, C, F, S, L or R will be ignored*
#
# There can be extra L and R that are not used. Only 'C' (count) number will be used.
    '''
    config_chain_A = f'''
#-------------------------------------
# B type DNA 10.5 bpt, generated by Avogadro
# Red and green
P {chain_A_path}
# E 
T pdb
C 1
F 4p
S 1.0
L 0 0 0
R 0 0 0
'''
    config_chain_B = f'''
#-------------------------------------
# B type DNA 10.5 bpt, generated by Avogadro
# Red and green
P {chain_B_path}
# E 
T pdb
C 1
F 4p
S 1.0
L 0 0 0
R 0 0 0'''
    
    config = header + config_chain_A + config_chain_B

    with open(config_path, "w") as f:
        f.write(header + config)




if __name__ == "__main__":

    # long_pdb_path = os.path.join(folder_path, "M13mp18-start48-end48.pdb")  # Replace with your PDB file path

    folder_path = "D:/Han/projects/transpdb/init_pdb"

    # parse the DNA sequence from arguements
    parser = argparse.ArgumentParser()
    parser.add_argument("--dna_seq", help="DNA sequence to generate PDB file")

    dna_seq = parser.parse_args().dna_seq


    reverse_dna_seq = Seq(dna_seq).reverse_complement()[::-1]
    print(reverse_dna_seq)

    structure_A = init(folder_path, dna_seq, 'A')
    structure_B = init(folder_path, reverse_dna_seq, 'B')

    # # generate the incremented pdb file path name
    # output_file_A = os.path.join(folder_path, "chain_A.pdb")
    # if os.path.exists(output_file_A):
    #     i = 1
    #     while os.path.exists(output_file_A):
    #         output_file_A = os.path.join(folder_path, f"chain_A_{i}.pdb")
    #         i += 1

    # write_pdb(structure_A, output_file_A)  # Replace with the output file path

    # output_file_B = os.path.join(folder_path, "chain_B.pdb")
    # if os.path.exists(output_file_B):
    #     i = 1
    #     while os.path.exists(output_file_B):
    #         output_file_B = os.path.join(folder_path, f"chain_B_{i}.pdb")
    #         i += 1
    
    # write_pdb(structure_B, output_file_B)  # Replace with the output file path


    # # generate the config file to D:\2024.02.22_VRSim-Gutmann\SimulationFolder\PDB_Files\Configuration_Files
    # config_folder = "D:/2024.02.22_VRSim-Gutmann/SimulationFolder/PDB_Files/Configuration_Files"
    # datetime = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    # config_file = os.path.join(config_folder, f"{datetime}_config.txt")

    # write_config(config_file, output_file_A, output_file_B)








