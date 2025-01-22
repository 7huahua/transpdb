from Bio.PDB import PDBParser, PDBIO, Superimposer
from Bio.Seq import Seq
from Bio.PDB.Residue import Residue
import numpy as np
from datetime import datetime
import os
import argparse
import math

# %% read and write pdb file
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

def generate_chain(dna_seq, structure, chain_name, rotran, res_A, res_T, res_C, res_G, left_padding, right_padding):
    last_A, last_T, last_C, last_G = res_A.copy(), res_T.copy(), res_C.copy(), res_G.copy()

    # compute the pdb file
    for i,res_name in enumerate(dna_seq):
        print(res_name)
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
        if chain_name == 'A':
            # 处理A链的特殊情况
            if i == int(left_padding) + 1:  # left padding之后的第一个残基
                for atom in new_res.get_list():
                    try:
                        if atom.get_name() not in ['HCAP']:  # 只删除HCAP
                            res.add(atom)
                    except:
                        print(atom)
            elif i == len(dna_seq) - int(right_padding):  # right padding之前的最后一个残基
                for atom in new_res.get_list():
                    try:
                        if atom.get_name() not in ['HTER', 'OXT']:  # 只删除HTER和OXT
                            res.add(atom)
                    except:
                        print(atom)
            else:  # 中间残基
                for atom in new_res.get_list():
                    try:
                        if atom.get_name() not in ['HTER','OXT','HCAP']:
                            res.add(atom)
                    except:
                        print(atom)
        else:
            # B chain
            if i == len(dna_seq)-1:
                for atom in new_res.get_list():
                    try:
                        if atom.get_name() not in ['HCAP']:
                            res.add(atom)
                    except:
                        print(atom)
            else:
                for atom in new_res.get_list():
                    try:
                        if atom.get_name() not in ['HTER','OXT','HCAP']:
                            res.add(atom)
                    except:
                        print(atom)

        # add the residue to the structure
        structure[0][chain_name].add(res)

    return structure

def calculate_center(chain):
    """计算第一条链的质心，作为旋转轴的中点"""
    coords = [atom.coord for residue in chain for atom in residue.get_list()]
    center = np.mean(coords, axis=0)
    return center

def rotation_matrix_around_axis(axis, theta):
    """
    计算围绕任意轴旋转的旋转矩阵
    :param axis: 旋转轴 (x, y, z)
    :param theta: 旋转角度 (单位：弧度)
    :return: 3x3 旋转矩阵
    """
    axis = axis / np.linalg.norm(axis)  # 单位化旋转轴
    ux, uy, uz = axis
    cos_theta = math.cos(theta)
    sin_theta = math.sin(theta)
    return np.array([
        [cos_theta + ux**2 * (1 - cos_theta), ux * uy * (1 - cos_theta) - uz * sin_theta, ux * uz * (1 - cos_theta) + uy * sin_theta],
        [uy * ux * (1 - cos_theta) + uz * sin_theta, cos_theta + uy**2 * (1 - cos_theta), uy * uz * (1 - cos_theta) - ux * sin_theta],
        [uz * ux * (1 - cos_theta) - uy * sin_theta, uz * uy * (1 - cos_theta) + ux * sin_theta, cos_theta + uz**2 * (1 - cos_theta)]
    ])

def rotate_chain(chain, center, axis, theta):
    """围绕给定轴旋转链中的所有原子"""
    rotation_mat = rotation_matrix_around_axis(axis, theta)
    for residue in chain:
        for atom in residue.get_list():
            relative_coord = atom.coord - center  # 平移到中心
            rotated_coord = np.dot(rotation_mat, relative_coord)  # 应用旋转矩阵
            atom.coord = rotated_coord + center  # 平移回原始位置
    return chain

def init(folder_path,dna_seq, chain_name, left_padding, right_padding):
    rotran = RoTran()
    dna_seq = dna_seq.upper()
    dna_seq = dna_seq.replace('U', 'T')

    if chain_name == 'A':
        # Load PDB files for chain A
        structure_A_A = read_pdb(os.path.join(folder_path, "A.pdb"))
        structure_A_T = read_pdb(os.path.join(folder_path, "T.pdb"))
        structure_A_C = read_pdb(os.path.join(folder_path, "C.pdb"))
        structure_A_G = read_pdb(os.path.join(folder_path, "G.pdb"))

        chain_A_res_A = structure_A_A[0]['A'][1]
        chain_A_res_T = structure_A_T[0]['A'][1]
        chain_A_res_C = structure_A_C[0]['A'][1]
        chain_A_res_G = structure_A_G[0]['A'][1]
    else:
        # Load PDB files for chain B
        structure_A_A = read_pdb(os.path.join(folder_path, "TA.pdb"))
        structure_A_T = read_pdb(os.path.join(folder_path, "AT.pdb"))
        structure_A_C = read_pdb(os.path.join(folder_path, "GC.pdb"))
        structure_A_G = read_pdb(os.path.join(folder_path, "CG.pdb"))

        chain_A_res_A = structure_A_A[0]['B'][1]
        chain_A_res_T = structure_A_T[0]['B'][1]
        chain_A_res_C = structure_A_C[0]['B'][1]
        chain_A_res_G = structure_A_G[0]['B'][1]

    if dna_seq[0] == 'A':
        structure = structure_A_A
    elif dna_seq[0] == 'T':
        structure = structure_A_T
    elif dna_seq[0] == 'C':
        structure = structure_A_C
    elif dna_seq[0] == 'G':
        structure = structure_A_G

    # filter the HCAP atom in the first residue if len(dna_seq) > 1
    if len(dna_seq) > 1:
        if chain_name == 'A':
            for atom in structure[0]['A'][1].get_list():
                if atom.get_name() == 'HCAP':
                    structure[0]['A'][1].detach_child(atom.id)
        else:
            for atom in structure[0]['B'][1].get_list():
                if atom.get_name() in ['HTER', 'OXT']:
                    structure[0]['B'][1].detach_child(atom.id)

    new_structure = generate_chain(dna_seq, structure, 
                                   chain_name, rotran, chain_A_res_A, chain_A_res_T, chain_A_res_C, chain_A_res_G, left_padding, right_padding)

    if chain_name == 'A':
        # Load PDB files
        pass
    elif chain_name == 'B':
        # 重新排列残基编号，使其符合3'-5'方向
        # 删除A链
        # 删除A链
        if 'A' in new_structure[0]:
            new_structure[0].detach_child('A')
        chain = new_structure[0]['B']
        residue_list = list(chain.child_list)
        total_residues = len(residue_list)
        
        # 先从链中移除所有残基
        for residue in residue_list:
            chain.detach_child(residue.id)
        
        # 以相反的顺序重新添加残基，并赋予新的编号
        for i, residue in enumerate(reversed(residue_list)):
            new_number = i + 1
            residue.id = (' ', new_number, ' ')
            chain.add(residue)
        
        chain.id = 'B'
        
        
        # 对每个残基中的原子重新排序和编号
        for residue in chain:
            # 对原子按名称排序
            # residue.child_list.sort(key=lambda x: x.get_name())
            # 重新编号
            serial_number = 1
            for atom in residue:
                atom.serial_number = serial_number
                serial_number += 1


        # # 沿着整条链的中心坐标，x轴旋转180度
        # center = calculate_center(new_structure[0]['A'])  # 计算第一条链的质心
        # rotation_axis = np.array([0, 1, 0])  # 平行于 y 轴的旋转轴
        # rotate_chain(new_structure[0]['A'], center, rotation_axis, math.pi)  # 绕旋转轴旋转 180°


        # # rename chain B to A
        # new_structure[0]['A'].id = 'B'
        # renumber the serial number of atoms in chain A
        # for residue in new_structure[0]['B']:
        #     for i, atom in enumerate(residue.get_list(), start=1):
        #         atom.serial_number = i
        
    if left_padding:
        # delete the first spaces_num residues
        if chain_name == 'A':
            chain = new_structure[0][chain_name]
            residues_to_remove = list(chain.child_list)[:int(left_padding)] + list(chain.child_list)[-int(right_padding):]
            for residue in residues_to_remove:
                chain.detach_child(residue.id)
            
            # 重新编号剩余的残基
            for i, residue in enumerate(chain, start=1):
                residue.id = (' ', i, ' ')
        else:
            # 目前B链似乎不存在需要删除的情况
            pass
            # chain = new_structure[0][chain_name]
            # residues_to_remove = list(chain.child_list)[-int(spaces_num):]
            # for residue in residues_to_remove:
            #     chain.detach_child(residue.id)

    return new_structure

def old_init(folder_path,dna_seq, chain_name, spaces_num=0):

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

        # filter the HCAP atom in the first residue if len(dna_seq) > 1
        if len(dna_seq) > 1:
            for atom in structure[0]['A'][1].get_list():
                if atom.get_name() == 'HCAP':
                    structure[0]['A'][1].detach_child(atom.id)
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
        
        if dna_seq[0] == 'A':
            structure = structure_B_A
        elif dna_seq[0] == 'T':
            structure = structure_B_T
        elif dna_seq[0] == 'C':
            structure = structure_B_C
        elif dna_seq[0] == 'G':
            structure = structure_B_G
        
        if len(dna_seq) > 1:
            for atom in structure[0]['B'][1].get_list():
                if atom.get_name() == 'HCAP':
                    structure[0]['B'][1].detach_child(atom.id)
        new_structure = generate_chain(dna_seq, structure, 'B', rotran, chain_B_res_A, chain_B_res_T, chain_B_res_C, chain_B_res_G)
        # delete chain A
        delete_chain(new_structure, 'A')
        # rename chain B to A
        # new_structure[0]['B'].id = 'A'
        # renumber the serial number of atoms in chain A
        for residue in new_structure[0]['B']:
            for i, atom in enumerate(residue.get_list(), start=1):
                atom.serial_number = i
        
    if spaces_num:
        # delete the first spaces_num residues
        chain = new_structure[0][chain_name]
        residues_to_remove = list(chain.child_list)[:int(spaces_num)]
        for residue in residues_to_remove:
            chain.detach_child(residue.id)
        
        # 重新编号剩余的残基
        for i, residue in enumerate(chain, start=1):
            residue.id = (' ', i, ' ')
    return new_structure

def delete_chain(structure, chain_id):
    # 获取第一个模型
    model = structure[0]
    
    # 检查链是否存在
    if chain_id in model:
        # 删除指定的链
        del model[chain_id]
        print(f"已删除链 {chain_id}")
    else:
        print(f"链 {chain_id} 不存在")

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


# %% add connect information
def extract_connect_info(pdb_path):
    """从 PDB 文件提取 CONNECT 信息"""
    connect_dict = {}
    with open(pdb_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith("CONECT"):
                # 提取中心原子和连接原子
                parts = line.split()
                center_atom = int(parts[1])
                connected_atoms = list(map(int, parts[2:]))
                connect_dict[center_atom] = connected_atoms
    return connect_dict

def generate_connect_info(chain, connect_info_template):
    """
    根据 DNA 链和模板 CONNECT 信息生成完整的 CONNECT 数据
    :param chain: 当前链对象，包含所有残基和原子信息
    :param connect_info_template: 模板 CONNECT 信息，按残基类型组织
    :return: CONNECT 信息列表，每行一个字符串
    """
    connect_lines = []
    previous_o3_serial = None  # 记录上一个碱基的 O3' 原子编号
    global_atom_serial = 1  # 全局原子编号，用于生成 PDB 文件的原子编号

    for residue in chain:
        residue_name = residue.get_resname()[-1]  # 获取残基名称 (A/T/C/G)
        if residue_name not in connect_info_template:
            continue  # 如果模板中没有该残基类型，跳过

        residue_template = connect_info_template[residue_name]  # 获取模板
        shift = global_atom_serial - residue.get_list()[0].serial_number  # 计算原子编号偏移量

        for atom in residue.get_list():
            connected_atoms = []
            for connected_atom in residue_template[atom.serial_number]:
                if connected_atom == 2 and previous_o3_serial is not None:
                    # 如果模板中的编号是 2，将其替换为上一个 O3' 原子编号
                    connected_atoms.append(previous_o3_serial)
                else:
                    connected_atoms.append(connected_atom + shift)
            
            # line = f"CONECT{residue_first_atom_serial + original_atom - 1:>5}"
            # line += "".join(f"{atom:>5}" for atom in adjusted_line)
            line = (f"CONECT{atom.serial_number + shift:>5}")
            line += "".join(f"{atom:>5}" for atom in connected_atoms)
            connect_lines.append(line)

            # 更新 previous_o3_serial
            if atom.get_name() == "O3'":
                previous_o3_serial = global_atom_serial

            global_atom_serial += 1

    return connect_lines


def write_pdb_with_connect(output_path, connect_lines):
    """将结构和 CONNECT 信息写入 PDB 文件"""
    with open(output_path, 'a+') as pdb_file:
        # # 写入 ATOM 信息
        # for atom in structure.get_atoms():
        #     pdb_file.write(atom.get_pdb_string())
        
        # 写入 CONNECT 信息
        for line in connect_lines:
            pdb_file.write(line + '\n')


if __name__ == "__main__":

    # long_pdb_path = os.path.join(folder_path, "M13mp18-start48-end48.pdb")  # Replace with your PDB file path

    folder_path = "./init_pdb"

    # parse the DNA sequence from arguements
    parser = argparse.ArgumentParser()
    parser.add_argument("--dna_seq", help="DNA sequence to generate PDB file")
    # add A or B（indicate the rotation direction)
    parser.add_argument("--chain_type", help="select from A or B")
    parser.add_argument("--output_path", help="output path")
    parser.add_argument("--left_padding",help="spaces on the left")
    parser.add_argument("--right_padding",help="spaces on the right")

    dna_seq = parser.parse_args().dna_seq
    chain_type = parser.parse_args().chain_type
    output_path = parser.parse_args().output_path
    left_padding = parser.parse_args().left_padding
    right_padding = parser.parse_args().right_padding

    structure_A = init(folder_path, dna_seq, chain_type, left_padding, right_padding)
    # structure_B = init(folder_path, reverse_dna_seq, 'B')

    # # generate the incremented pdb file path name
    # output_file_A = os.path.join(folder_path, "chain_A.pdb")
    # if os.path.exists(output_file_A):
    #     i = 1
    #     while os.path.exists(output_file_A):
    #         output_file_A = os.path.join(folder_path, f"chain_A_{i}.pdb")
    #         i += 1

    write_pdb(structure_A, output_path)  # Replace with the output file path

# %% add connect information
    connect_info = {
        'A': extract_connect_info(os.path.join(folder_path, "A.pdb")),
        'T': extract_connect_info(os.path.join(folder_path, "T.pdb")),
        'C': extract_connect_info(os.path.join(folder_path, "C.pdb")),
        'G': extract_connect_info(os.path.join(folder_path, "G.pdb")),
    }
    connect_lines = generate_connect_info(structure_A[0]['A'], connect_info)

    # 写入输出文件
    write_pdb_with_connect(output_path, connect_lines)

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








