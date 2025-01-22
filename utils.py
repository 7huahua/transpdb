import math
import numpy as np

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

def reverse_chain_x(chain):
    # 沿着整条链的中心坐标，x轴旋转180度
    center = calculate_center(chain)  # 计算第一条链的质心
    rotation_axis = np.array([0, 1, 0])  # 平行于 y 轴的旋转轴
    rotate_chain(chain, center, rotation_axis, math.pi)  # 绕旋转轴旋转 180°