from Bio.PDB import PDBParser, PDBIO
import os

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

def compute_config(folder_path, long_pdb_path, short_5_3, short_3_5):
    long_pdb_path = os.path.join(folder_path, long_pdb_path)
    short_5_3 = os.path.join(folder_path, short_5_3)
    short_3_5 = os.path.join(folder_path, short_3_5)

    # Load PDB files
    structure_long = read_pdb(long_pdb_path)
    structure_short_5_3 = read_pdb(short_5_3)
    structure_short_3_5 = read_pdb(short_3_5)

    # Select the chain to be rotated, assume it's chain A
    chain_long = structure_long[0]['A']  # Ensure the chain ID is correct
    chain_short_5_3 = structure_short_5_3[0]['A']  # Ensure the chain ID is correct
    chain_short_3_5 = structure_short_3_5[0]['A']  # Ensure the chain ID is correct

    # # The starting positions of A and B are the same. We need to calculate the matching sequence position for A and B, then calculate the centroid based on the position
    # chain_long_sequence = ''.join([residue.get_resname()[1] for residue in chain_long])
    # chain_short_5_3_sequence = ''.join([residue.get_resname()[1] for residue in chain_short_5_3])
    # chain_short_3_5_sequence = ''.join([residue.get_resname()[1] for residue in chain_short_3_5])

    # print(chain_long_sequence)
    # print(chain_short_5_3_sequence)
    # print(chain_short_3_5_sequence)

    # longは5'->3'なので、short_5_3を回転させて一部を揃え、short_3_5は直接揃えた後、
    # 具体的には、これらの三つの鎖は中央で自動的に揃います。
    # したがって、short_5_3を左にlen(long)/2 - (len(short_5_3) - 超過部分)/2 * 3.4だけシフトし、その値をintに丸める必要があります。
    # そして、short_3_5を右にlen(long)/2 - (len(short_3_5) - 超過部分)/2 * 3.4だけシフトし、その値をintに丸める必要があります。

    # 超過部分は手動で指定されており、ここでは8です。
    outer = 8

    align_5_3 = (len(chain_long) - len(chain_short_5_3)) // 2 + 1 + outer
    align_3_5 = len(chain_long) - (len(chain_long) - len(chain_short_5_3)) // 2 - outer

    shift_5_3 = chain_long[1]["P"].get_coord()[0] - chain_long[align_5_3]["P"].get_coord()[0]
    shift_3_5 = chain_long[len(chain_long)]["P"].get_coord()[0] - chain_long[align_3_5]["P"].get_coord()[0]

    # 設定ファイルの内容は以下の通りです：

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
    config_short_5_3 = f'''
#-------------------------------------
# B type DNA 10.5 bpt, generated by Avogadro
# Red and green
P {short_5_3}
T pdb
C 1
F 8p
S 1.0
L {shift_5_3} 0  0
R 0 180 0'''
    config_long = f'''
#-------------------------------------
# B type DNA 10.5 bpt, generated by Avogadro
# Red and green
P {long_pdb_path}
T pdb
C 1
F 8p
S 1.0
L 0 0 0
R 0 0 0'''
    
    config_short_3_5 = f'''
#-------------------------------------
# B type DNA 10.5 bpt, generated by Avogadro
# Red and green
P {short_3_5}
T pdb
C 1
F 8p
S 1.0
L {shift_3_5} 0  0
R 0 0 0'''
    
    config = config_short_5_3 + "\n" + config_long + "\n" + config_short_3_5

    with open(os.path.join(folder_path, "config.txt"), "w") as f:
        f.write(header + config)

    return os.path.join(folder_path, "config.txt")


if __name__ == "__main__":

    # folder path
    folder_path = "C:/Users/Han/Downloads/10.5/10.5/"
    # Load PDB file
    long_pdb_path = os.path.join(folder_path, "M13mp18-start48-end48.pdb")  # Replace with your PDB file path
    short_5_3 = os.path.join(folder_path, "TAATAGTAGTAGCATTAAAAAAAA-10.5.pdb")  # Replace with your PDB file path
    short_3_5 = os.path.join(folder_path, "ACTAAATAACCTACAATTTTTTTT-10.5.pdb")  # Replace with your PDB file path
   
    short_5_3_charge_path = os.path.join(folder_path, "TAATAGTAGTAGCATTCCCCCCCC-EEM-Geidl 2015 (Cheminf_b3lyp_aim).txt")  # Replace with your charge file path
    short_3_5_charge_path = os.path.join(folder_path, "ACTAAATAACCTACAAGGGGGGGG-EEM-Geidl 2015 (Cheminf_b3lyp_aim).txt")  # Replace with your charge file path
    long_charge_path = os.path.join(folder_path, "M13mp18-start48-end48-EEM-Geidl 2015 (Cheminf_b3lyp_aim).txt")  # Replace with your charge file path  

