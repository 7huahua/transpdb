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

def run_generate_pdb_and_config():
    if DEBUG_MODE:
        pdb_folder = r"D:\Han\projects\transpdb\han\generated_pdb"
        config_folder = r"D:\2024.02.22_VRSim-Gutmann\SimulationFolder\PDB_Files\Configuration_Files"
        debug_input = "      gcta\natcgatcgatcgatcgatcgatcgatcgatcg\n                tagc"
        
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

    for index, sequence_data in enumerate(parsed_sequences):
        try:
            sequence_name = clean_filename(sequence_data['sequence'])
            output_path = os.path.join(pdb_folder, f'{sequence_name}.pdb')
            if os.path.exists(output_path):
                output_text.insert(tk.END, f"PDB file {index + 1} already exists: {output_path}\n")
                pdb_paths.append(output_path)
                continue
            result = subprocess.run(
                ['python', 'generate_pdb.py', '--dna_seq', sequence_data['sequence'], '--chain_type', 'A', '--output_path', output_path],
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

def test_align():
    if DEBUG_MODE:
        pdb_folder = r"D:\Han\projects\transpdb\han\generated_pdb"
        config_folder = r"D:\2024.02.22_VRSim-Gutmann\SimulationFolder\PDB_Files\Configuration_Files"
        debug_input = "      gcta\natcgatcgatcgatcgatcgatcgatcgatcg\n                tagc"
        
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
    # the dna base pair: A-T, C-G, get the complement of the longest sequence
    complement_longest_sequence = longest_sequence['sequence'].replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c')
    # upper case
    complement_longest_sequence = complement_longest_sequence.upper()
    folder_path = 'init_pdb'
    complement_seq_structure = init(folder_path, complement_longest_sequence, 'B')

    # 对于每条序列：
    #     如果是最长序列，则生成chain A 的pdb文件
    #     否则：
    #         获取对齐位置和短链的序列
    #         对齐序列第一个核苷酸骨架和对应位置的最长链B链的核苷酸骨架
    #         根据骨架的偏移获取第一个核苷酸位置
    #         根据rotran生成chain B 的pdb文件
    #         返回文件path

    # 对于所有path，生成config文件
    sequence_data=longest_sequence

    try:
        # if it is the longest seq, generate the pdb file 
        if sequence_data['sequence'] == longest_sequence['sequence']:
            chain_type = 'A'
            sequence_name = clean_filename(sequence_data['sequence'])
            output_path = os.path.join(pdb_folder, f'{sequence_name}_{chain_type}.pdb')
            if os.path.exists(output_path):
                output_text.insert(tk.END, f"PDB file {sequence_name} already exists: {output_path}\n")
                pdb_paths.append(output_path)
            result = subprocess.run(
                ['python', 'generate_pdb.py', '--dna_seq', sequence_data['sequence'], '--chain_type', chain_type, '--output_path', output_path],
                capture_output=True,
                text=True
            )
            
            if result.returncode == 0:
                pdb_paths.append(output_path)
                output_text.insert(tk.END, f"PDB file {sequence_name} generated successfully: {output_path}\n")
            else:
                output_text.insert(tk.END, f"Failed to generate PDB file {sequence_name}, please check!\n")
                output_text.insert(tk.END, result.stderr)
                return  # Stop processing if any PDB generation fails
        else:
            # get the align position and the short seq
            chain_type = 'B'
            sequence_name = clean_filename(complement_longest_sequence)
            output_path = os.path.join(pdb_folder, f'{sequence_name}_{chain_type}.pdb')
            if os.path.exists(output_path):
                output_text.insert(tk.END, f"PDB file {sequence_name} already exists: {output_path}\n")
                pdb_paths.append(output_path)
            result = subprocess.run(
                ['python', 'generate_pdb.py', '--dna_seq', complement_longest_sequence, '--chain_type', chain_type, '--output_path', output_path],
                capture_output=True,
                text=True
            )

            if result.returncode == 0:
                pdb_paths.append(output_path)
                output_text.insert(tk.END, f"PDB file {sequence_name} generated successfully: {output_path}\n")
            else:
                output_text.insert(tk.END, f"Failed to generate PDB file {sequence_name}, please check!\n")
                output_text.insert(tk.END, result.stderr)
                return  # Stop processing if any PDB generation fails
    
    except Exception as e:
        messagebox.showerror("Error", f"Problem running generate_pdb.py: {e}")
        return

    if len(pdb_paths) == 2:
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

def run_generate_pdb_config_align():
    global pdb_folder, config_folder
    if DEBUG_MODE:

        pdb_folder = r"D:\Han\projects\transpdb\han\generated_pdb"
        config_folder = r"D:\2024.02.22_VRSim-Gutmann\SimulationFolder\PDB_Files\Configuration_Files"
        # debug_input = "      gcta\natcgatcgatcgatcgatcgatcgatcgatcg\n                tagc"
        debug_input="gctagctagcta\n      atcgatcgatcgatcgatcgatcgatcgatcg\n                                  tagctagctagc"
        
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
    
        # 检查最长序列左边是否有空格
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

    
    # the dna base pair: A-T, C-G, get the complement of the longest sequence
    complement_longest_sequence = padded_sequence.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c')
    # upper case
    complement_longest_sequence = complement_longest_sequence.upper()
    folder_path = 'init_pdb'
    complement_seq_structure = init(folder_path, complement_longest_sequence, 'B')

    # 对于每条序列：
    #     如果是最长序列，
    #         如果最长序列前面空格为空：
    #           则生成chain A 的pdb文件
    #         如果不为空：
    #           需要根据最左端补齐chain a的pdb文件
    #     否则：
    #         获取对齐位置和短链的序列
    #         对齐序列第一个核苷酸骨架和对应位置的最长链B链的核苷酸骨架
    #         根据骨架的偏移获取第一个核苷酸位置
    #         根据rotran生成chain B 的pdb文件
    #         返回文件path

    # 对于所有path，生成config文件

    for index, sequence_data in enumerate(parsed_sequences):
        try:

            # if it is the longest seq, generate the pdb file 
            if sequence_data['sequence'] == longest_sequence['sequence']:
                chain_type = 'A'
                spaces_num=sequence_data['spaces']
                sequence_name = clean_filename(sequence_data['sequence'])
                output_path = os.path.join(pdb_folder, f'{sequence_name}_{chain_type}_space{spaces_num}.pdb')
                if os.path.exists(output_path):
                    output_text.insert(tk.END, f"PDB file {index + 1} already exists: {output_path}\n")
                    pdb_paths.append(output_path)
                    continue

                if not spaces_num:
                    result = subprocess.run(
                        ['python', 'generate_pdb.py', '--dna_seq', sequence_data['sequence'], 
                        '--chain_type', chain_type, '--output_path', output_path, 
                        "--spaces_num", str(0)],
                        capture_output=True,
                        text=True
                    )
                else:
                    result = subprocess.run(
                        ['python', 'generate_pdb.py', '--dna_seq', padded_sequence, 
                        '--chain_type', chain_type, '--output_path', output_path, 
                        "--spaces_num", str(spaces_num)],
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
                # get the align position and the short seq
                align_position = sequence_data['spaces']
                chain_type = 'B'
                sequence_name = clean_filename(sequence_data['sequence'].upper())
                output_path = os.path.join(pdb_folder, f'{sequence_name}_{chain_type}_space{align_position}.pdb')
                if os.path.exists(output_path):
                    output_text.insert(tk.END, f"PDB file {index + 1} already exists: {output_path}\n")
                    pdb_paths.append(output_path)
                    continue
                # result = subprocess.run(
                #     ['python', 'align_pdb.py', '--structure',complement_seq_structure,'--dna_seq', sequence_data['sequence'], '--chain_type', chain_type, '--output_path', output_path, '--align_position', align_position],
                #     capture_output=True,
                #     text=True
                # )

                # structure, dna_seq, chain_type, output_path, align_position
                align_pdb = AlignPDB(complement_seq_structure, sequence_data['sequence'], chain_type, output_path, align_position)
                
                try:
                    align_pdb.align_pdb()

                # if result.returncode == 0:
                    pdb_paths.append(output_path)
                    output_text.insert(tk.END, f"PDB file {index + 1} generated successfully: {output_path}\n")
                except Exception as e:
                    output_text.insert(tk.END, f"Failed to generate PDB file {index + 1}, please check!\n")
                    output_text.insert(tk.END, str(e) + "\n")
                    return  # Stop processing if any PDB generation fails
        
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
            messagebox.showwarning("Warning", "Please enter PDB folder and config folder paths")
            return False
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