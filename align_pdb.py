import argparse
from generate_pdb import init, extract_connect_info, generate_connect_info, write_pdb_with_connect
from Bio.PDB import PDBParser, PDBIO, Superimposer
import os
from utils import *

def write_pdb(structure, output_file):
    """Save the structure to a PDB file"""
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file)

class AlignPDB:

    def __init__(self, structure, dna_seq, chain_type, output_path, align_position):
        self.structure = structure
        self.dna_seq = dna_seq
        self.chain_type = chain_type
        self.output_path = output_path
        self.align_position = align_position

    def get_residue_rotran(self, structure_long, structure_short):
        group_atom_names = ["P", "O1P", "O2P", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'"]

        # for every residue in the long chain, get the atoms of the group_atom_names
        long_groups = []
        chain_long=structure_long[0]['B']
        for residue in chain_long:
            group = [atom for atom in residue.get_list() if atom.get_name() in group_atom_names]
            long_groups.append(group)

        target=long_groups[self.align_position]

        chain_short=structure_short[0]['B']
        for residue in chain_short:
            source = [atom for atom in residue.get_list() if atom.get_name() in group_atom_names]
            sup = Superimposer()
            sup.set_atoms(target, source)
            rot,tran = sup.rotran
            break

        return rot, tran

    def align_pdb(self):
        init_folder_path = 'init_pdb'
        structure_short=init(init_folder_path,self.dna_seq, self.chain_type,0,0)
        structure_long=self.structure

        if self.chain_type == "B":
            # # 沿着整条链的中心坐标，x轴旋转180度
            # reverse_chain_x(structure_short[0]['B'])
            # reverse_chain_x(structure_long[0]['B'])

            rot, tran = self.get_residue_rotran(structure_long, structure_short)

            # apply the rotran to the short structure
            for residue in structure_short:
                for atom in residue.get_list():
                    atom.transform(rot, tran)

            # # reverse back
            # reverse_chain_x(structure_short[0]['B'])
            # reverse_chain_x(structure_long[0]['B'])

            # write the short structure to the output path
            structure_short[0]['B'].id = 'A'
            write_pdb(structure_short, self.output_path)

            connect_info = {
                'A': extract_connect_info(os.path.join(init_folder_path, "A.pdb")),
                'T': extract_connect_info(os.path.join(init_folder_path, "T.pdb")),
                'C': extract_connect_info(os.path.join(init_folder_path, "C.pdb")),
                'G': extract_connect_info(os.path.join(init_folder_path, "G.pdb")),
            }
            connect_lines = generate_connect_info(structure_short[0]['A'], connect_info)

            # 写入输出文件
            write_pdb_with_connect(self.output_path, connect_lines)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--structure", help="structure file")
    parser.add_argument("--dna_seq", help="DNA sequence")
    parser.add_argument("--chain_type", help="chain type")

    # 'python', 'align_pdb.py', '--structure',complement_seq_structure,'--dna_seq', sequence_data['sequence'], '--chain_type', chain_type, '--output_path', output_path, '--align_position', align_position
    structure = parser.parse_args().structure
    dna_seq = parser.parse_args().dna_seq
    chain_type = parser.parse_args().chain_type
    output_path = parser.parse_args().output_path
    align_position = parser.parse_args().align_position

    align_pdb = AlignPDB(structure, dna_seq, chain_type, output_path, align_position)
    align_pdb.align_pdb()

