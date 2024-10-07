import sys
import os

def calculate_shift(space_count):
    return space_count * 3.4

def generate_config_space(space_counts, pdb_paths, config_folder):
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
    
    config_content = header + "\n"
    
    for i, (space_count, pdb_path) in enumerate(zip(space_counts, pdb_paths)):
        shift_left = calculate_shift(space_count)
        shift_down = i * 20
        config_content += f'''
#-------------------------------------
# B type DNA 10.5 bpt, DNA sequence {i+1}
P {pdb_path}
T pdb
C 1
F 8
S 1.0
L {shift_left:.2f} {shift_down:.2f} 0
R 0 0 0
'''
    
    config_path = os.path.join(config_folder, "dna_alignment_config.txt")
    with open(config_path, "w") as f:
        f.write(config_content)
    
    return config_path


def generate_config(pdb_paths, config_folder):
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
    
    config_content = header + "\n"
    
    for i, pdb_path in enumerate(pdb_paths):
        config_content += f'''
#-------------------------------------
# B type DNA 10.5 bpt, DNA sequence {i+1}
P {pdb_path}
T pdb
C 1
F 4
S 1.0
L 0 0 0
R 0 0 0
'''
    
    config_path = os.path.join(config_folder, "dna_alignment_config.txt")
    with open(config_path, "w") as f:
        f.write(config_content)
    return config_path

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("用法: python generate_config.py <空格数1> <空格数2> ... <PDB路径1> <PDB路径2> ...")
        sys.exit(1)
    
    # search for all the pdb files
    pdb_paths = [arg for arg in sys.argv if arg.endswith('.pdb')]
    pdb_start_index = sys.argv.index(pdb_paths[0])
    pdb_end_index = sys.argv.index(pdb_paths[-1])

    # 分割参数列表
    space_counts = [int(arg) for arg in sys.argv[1:pdb_start_index]]

    config_folder = sys.argv[-1]

    
    if len(space_counts) != len(pdb_paths):
        print("错误：空格数量与PDB文件数量不匹配")
        sys.exit(1)
    
    config_path = generate_config(pdb_paths, config_folder)
    print(config_path)

