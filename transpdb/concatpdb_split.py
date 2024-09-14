from Bio.PDB import PDBParser, PDBIO, Superimposer
import numpy as np
from Bio.Seq import Seq

def read_pdb(file_path):
    """读取PDB文件并返回结构对象"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', file_path)
    return structure

def calculate_centroid(atoms):
    """计算原子集合的质心"""
    coords = [atom.get_coord() for atom in atoms]
    centroid = np.mean(coords, axis=0)
    return centroid

def rotate_y_180(atoms, centroid):
    """围绕质心沿y轴旋转180度"""
    rotation_matrix = np.array([
        [-1, 0, 0],  # x = -x
        [0, 1, 0],   # y = y
        [0, 0, -1]   # z = -z
    ])
    for atom in atoms:
        coord = atom.get_coord() - centroid
        new_coord = np.dot(rotation_matrix, coord) + centroid
        atom.set_coord(new_coord)

def translate_x(atoms, distance):
    """在x轴方向平移给定距离"""
    for atom in atoms:
        coord = atom.get_coord()
        coord[0] += distance  # 只在x轴方向加上平移距离
        atom.set_coord(coord)

def rotate_x(atoms, angle, centroid):
    """围绕质心沿x轴旋转给定角度"""
    rad = np.deg2rad(angle)  # 将角度转换为弧度
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
    """将结构保存到PDB文件"""
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file)


def get_residue_atoms(structure, chain_id, start, end):
    """从指定链中获取从start到end的残基的所有原子"""
    atoms = []
    for residue_id in range(start, end + 1):
        residue = structure[0][chain_id][residue_id]
        atoms.extend(residue.get_atoms())
    return atoms

def calculate_rotation_angle(coord1, coord2):
    """计算两个坐标之间围绕x轴的旋转角度"""
    # 计算向量
    vector = coord2 - coord1
    
    # 计算围绕x轴的旋转角度
    angle_rad = np.arctan2(vector[2], vector[1])  # 在y-z平面上的角度
    angle_deg = np.degrees(angle_rad)
    return angle_deg

# 加载PDB文件
pdb_path = 'D:/Han/projects/transpdb/random.pdb'  # 替换为你的PDB文件路径
structure_A = read_pdb(pdb_path)

# 选择要旋转的链，假设是链B
chain_A = structure_A[0]['A']  # 确保链ID正确

# 加载PDB文件
pdb_path = 'D:/Han/projects/transpdb/random_double.pdb'  # 替换为你的PDB文件路径
structure_B = read_pdb(pdb_path)

# 选择要旋转的链，假设是链B
chain_A = structure_A[0]['A']  # 确保链ID正确
chain_B = structure_B[0]['B']  # 确保链ID正确


# A B的起始位置是相同的。我们需要计算出AB匹配的sequence的位置，然后根据位置计算质心
chain_A_sequence = ''.join([residue.get_resname()[1] for residue in chain_A])
chain_B_sequence = ''.join([residue.get_resname()[1] for residue in chain_B])

# 利用bipython的Seq获取匹配的位置
chain_A_seq = Seq(chain_A_sequence)
chain_B_seq = Seq(chain_B_sequence)

# 获取chain A的reverse_complement()
chain_A_seq_reverse = chain_A_seq.reverse_complement()

target_seq = ""

# 获取匹配的位置，也就是chain_b第一次出现的位置的偏移
shift = chain_A_seq_reverse.find(chain_B_seq)
print(shift)

# 假设每个残基之间的距离是3.4A
# diatance = shift * 3.4
# 计算chain A的第一个残基的X坐标
chain_A_first_residue_atoms = chain_A[1].get_atoms()
first_coords = calculate_centroid(chain_A_first_residue_atoms)

chain_A_shift_residue_atoms = chain_A[shift+1].get_atoms()
shift_coords = calculate_centroid(chain_A_shift_residue_atoms)

distance = shift_coords[0] - first_coords[0]

# shift chain b
translate_x(list(chain_B.get_atoms()), distance)  # distance是你计算的或估计的值

chain_id = 'A'  # 请根据实际情况替换链ID
# 计算质心
# tg_atoms = get_residue_atoms(chain_B, 'A', start_residue, end_residue)  # 需要正确指定起始和结束的残基索引
centroid = calculate_centroid(chain_B.get_atoms())
print(centroid)


# 计算旋转角度
rotation_angle = calculate_rotation_angle(shift_coords, first_coords)
print(f"Rotation angle around the x-axis is {rotation_angle} degrees")


angle = -rotation_angle  # 这里假设已经计算出需要旋转30度
rotate_x(list(chain_B.get_atoms()), angle, centroid)
# 应用旋转
rotate_y_180(list(chain_B.get_atoms()), centroid)

# 合并AB链

# 合并两个结构,并且B链的链名改为B
chain_B.id = 'B'  # 修改链ID
structure_A[0].add(chain_B)


# 保存合并后的PDB文件
output_path = 'D:/Han/projects/transpdb/random_merged.pdb'  # 输出文件名
write_pdb(structure_A, output_path)



