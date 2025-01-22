import tkinter as tk
from tkinter import messagebox
import subprocess
import os
import Bio.Align
from generate_pdb import *
from align_pdb import *

import logging
logging.basicConfig(filename="error.log", level=logging.DEBUG)


DEBUG_MODE = False
# ---------------------------------------------------
# D:\Han\projects\transpdb\han\generated_pdb
# D:\2024.02.22_VRSim-Gutmann\SimulationFolder\PDB_Files\Configuration_Files
# -----------------------------------
#       gctagctagcta
# atcgatcgatcgatcgatcgatcgatcgatcg
#                 tagctagctagc
# --------------------------------------
# gctagctagcta
#       atcgatcgatcgatcgatcgatcgatcgatcg
#                                   tagctagctagc
# ----------------------------------------------

def parse_dna_sequence(input_text):
    """
    Parse user input DNA sequence, record the number of spaces and sequence information for each line
    :param input_text: Multi-line sequence text input by the user
    :return: List containing space count and sequence information
    """
    sequence_data = []
    lines = input_text.splitlines()

    for line in lines:
        leading_spaces = len(line) - len(line.lstrip(' '))
        dna_sequence = line.strip()

        if dna_sequence:
            sequence_data.append({
                'spaces': leading_spaces,
                'sequence': dna_sequence
            })

    return sequence_data

def clean_filename(filename):
    """Clean filename, remove illegal characters"""
    illegal_chars = ['<', '>', ':', '"', '/', '\\', '|', '?', '*', '\n', '\r', '\t']
    for char in illegal_chars:
        filename = filename.replace(char, '')
    filename = filename.strip('. ')
    return filename if filename else 'unnamed_sequence'

#  %% current code
def run_generate_pdb_config_align():
    global pdb_folder, config_folder
    if DEBUG_MODE:

        pdb_folder = r"D:\Han\projects\transpdb\han\generated_pdb"
        config_folder = r"D:\2024.10.04_VR-Sim_Gutmann_v1.2.0\2024.10.04_VR-Sim_Gutmann_v1.2.0\SimulationFolder\PDB_Files\Configuration_Files"
        # config_folder = r"D:\2024.02.22_VRSim-Gutmann\SimulationFolder\PDB_Files\Configuration_Files"
        debug_input = "      gcta\natcgatcgatcgatcgatcgatcgatcgatcg\n                tagc"
        debug_input="gctagctagcta\n      atcgatcgatcgatcgatcgatcgatcgatcg\n                                  tagctagctagc"
        
        # debug_input = "tagctagctagctagc\natcgatcgatcgatcgatcgatcgatcgatcg"

        entry_pdb_folder.delete(0, tk.END)
        entry_pdb_folder.insert(0, pdb_folder)
        entry_config_folder.delete(0, tk.END)
        entry_config_folder.insert(0, config_folder)
        text_input.delete("1.0", tk.END)
        text_input.insert("1.0", debug_input)

    if not get_folder_paths():
        return

    user_input = text_input.get("1.0", tk.END)
    if not user_input.strip():
        messagebox.showwarning("Input Error", "Please enter DNA sequence")
        return
    
    parsed_sequences = parse_dna_sequence(user_input)
    pdb_paths = []

    output_text.delete("1.0", tk.END)  # Clear output text box

    # 找到最长序列并产生对应的B链的pdb文件并标记为complement seq pdb
    # get the longest sequence
    longest_sequence = max(parsed_sequences, key=lambda x: len(x['sequence']))
    longest_sequence['sequence'] = longest_sequence['sequence'].upper()
    
    if longest_sequence['spaces'] > 0:
        # 找到左边没有空格的短链
        left_aligned_sequence = next((seq for seq in parsed_sequences if seq['spaces'] == 0), None)
        
        if left_aligned_sequence:
            # 获取需要填充的碱基数量
            padding_length = min(longest_sequence['spaces'], len(left_aligned_sequence['sequence']))
            
            # 获取短链中前k个元素的互补序列
            padding = left_aligned_sequence['sequence'][:padding_length].upper()
            padding = padding.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c')
            
            # 将填充添加到最长链的左边
            padded_sequence = padding + longest_sequence['sequence']
        else:
            # 如果没有找到左对齐的序列，使用'N'填充
            padded_sequence = 'N' * longest_sequence['spaces'] + longest_sequence['sequence']
    else:
        padded_sequence = longest_sequence['sequence']

    # 处理右侧填充
    longest_end_position = longest_sequence['spaces'] + len(longest_sequence['sequence'])
    right_padding = ''
    
    # 查找最右端的序列
    for seq in parsed_sequences:
        seq_end_position = seq['spaces'] + len(seq['sequence'])
        if seq_end_position > longest_end_position:
            # 计算需要填充的长度
            padding_needed = seq_end_position - longest_end_position
            # 获取该序列末尾部分作为填充
            padding_sequence = seq['sequence'][-padding_needed:].upper()
            # 获取互补序列
            right_padding = padding_sequence.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()
            break
    
    # 添加右侧填充
    padded_sequence = padded_sequence + right_padding
    
    # the dna base pair: A-T, C-G, get the complement of the longest sequence
    complement_longest_sequence = padded_sequence.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c')
    # upper case
    complement_longest_sequence = complement_longest_sequence.upper()
    folder_path = 'init_pdb'
    complement_seq_structure = init(folder_path, complement_longest_sequence, 'B', 0, 0)

    # 计算整体序列的总长度（包括空格）
    total_length = max(seq['spaces'] + len(seq['sequence']) for seq in parsed_sequences)
    
    for index, sequence_data in enumerate(parsed_sequences):
        try:
            # 如果是最长序列，生成chain A的pdb文件
            if sequence_data['sequence'] == longest_sequence['sequence']:
                chain_type = 'A'
                sequence_name = clean_filename(sequence_data['sequence'])
                
                # 计算左右padding的长度
                left_padding = longest_sequence['spaces'] if longest_sequence['spaces'] > 0 else 0
                right_padding = len(right_padding) if right_padding else 0
                
                output_path = os.path.join(pdb_folder, f'{sequence_name}_{chain_type}_left{left_padding}_right{right_padding}.pdb')
                
                # 使用统一的命令调用方式
                result = subprocess.run(
                    ['python', 'generate_pdb.py',
                     '--dna_seq', padded_sequence if left_padding or right_padding else sequence_data['sequence'],
                     '--chain_type', chain_type,
                     '--output_path', output_path,
                     '--left_padding', str(left_padding),
                     '--right_padding', str(right_padding)],
                    capture_output=True,
                    text=True
                )
                
                if result.returncode == 0:
                    pdb_paths.append(output_path)
                    output_text.insert(tk.END, f"PDB file {index + 1} generated successfully: {output_path}\n")
                else:
                    output_text.insert(tk.END, f"Failed to generate PDB file {index + 1}, please check!\n")
                    output_text.insert(tk.END, result.stderr)
                    return  # Stop processing if any PDB generation fails
            else:
                # 计算该序列的右端位置
                seq_end_position = sequence_data['spaces'] + len(sequence_data['sequence'])
                # 计算该序列距离最右端的距离
                right_padding_length = total_length - seq_end_position
                
                chain_type = 'B'
                sequence_name = clean_filename(sequence_data['sequence'][::-1].upper())
                output_path = os.path.join(pdb_folder, 
                    f'{sequence_name}_{chain_type}_space{sequence_data["spaces"]}_right{right_padding_length}.pdb')
                
                # 更新AlignPDB调用，添加right_padding_length参数
                align_pdb = AlignPDB(
                    complement_seq_structure, 
                    sequence_data['sequence'], 
                    chain_type, 
                    output_path,
                    right_padding_length
                )
                
                try:
                    align_pdb.align_pdb()
                    pdb_paths.append(output_path)
                    output_text.insert(tk.END, f"PDB file {index + 1} generated successfully: {output_path}\n")
                except Exception as e:
                    output_text.insert(tk.END, f"Failed to generate PDB file {index + 1}, please check!\n")
                    output_text.insert(tk.END, str(e) + "\n")
                    return
        
        except Exception as e:
            messagebox.showerror("Error", f"Problem running generate_pdb.py: {e}")
            return

    if len(pdb_paths) == len(parsed_sequences):
        try:
            command = [
                'python', 
                'generate_config.py'
            ] + [str(seq['spaces']) for seq in parsed_sequences] + pdb_paths + [
                '--config_folder', 
                config_folder
            ]
            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=True
            )

            
            if result.returncode == 0:
                config_path = result.stdout.strip()
                output_text.insert(tk.END, f"\nConfig file generated successfully: {config_path}\n")
            else:
                output_text.insert(tk.END, "Failed to generate Config file, please check!\n")
                output_text.insert(tk.END, result.stderr)
        
        except Exception as e:
            messagebox.showerror("Error", f"Problem running generate_config.py: {e}")
    else:
        output_text.insert(tk.END, "Config file not generated due to incomplete PDB file generation.\n")

try:
    root = tk.Tk()
    root.title("DNA Alignment Tool")

    pdb_folder_label = tk.Label(root, text="PDB Folder Path:")
    pdb_folder_label.pack()

    entry_pdb_folder = tk.Entry(root, width=50)
    entry_pdb_folder.pack()

    config_folder_label = tk.Label(root, text="Config Folder Path:")
    config_folder_label.pack()

    entry_config_folder = tk.Entry(root, width=50)
    entry_config_folder.pack()

    def get_folder_paths():
        global pdb_folder, config_folder
        pdb_folder = entry_pdb_folder.get()
        config_folder = entry_config_folder.get()
        pdb_folder = pdb_folder.replace('\n', '')
        config_folder = config_folder.replace('\n', '')
        if not pdb_folder or not config_folder:
            # messagebox.showwarning("Warning", "Please enter PDB folder and config folder paths")
            # return False
            pdb_folder = r"D:\Han\projects\transpdb\han\generated_pdb"
            config_folder = r"D:\2024.10.04_VR-Sim_Gutmann_v1.2.0\2024.10.04_VR-Sim_Gutmann_v1.2.0\SimulationFolder\PDB_Files\Configuration_Files"

        return True

    label_input = tk.Label(root, text="Please enter DNA sequence:")
    label_input.pack()

    text_input = tk.Text(root, height=10, width=50)
    text_input.pack()

    parse_button = tk.Button(root, text="Generate PDB and Config", command=run_generate_pdb_config_align)
    parse_button.pack()

    label_output = tk.Label(root, text="Parsing Result:")
    label_output.pack()

    output_text = tk.Text(root, height=10, width=50)
    output_text.pack()

    root.mainloop()

except Exception as e:
    logging.exception("An error occurred")
    input("Press Enter to exit...")  # 暂停程序，便于查看错误