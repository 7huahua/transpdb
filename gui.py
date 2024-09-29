import tkinter as tk
from tkinter import messagebox
import subprocess
import os

DEBUG_MODE = False

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

parse_button = tk.Button(root, text="Generate PDB and Config", command=run_generate_pdb_and_config)
parse_button.pack()

label_output = tk.Label(root, text="Parsing Result:")
label_output.pack()

output_text = tk.Text(root, height=10, width=50)
output_text.pack()

root.mainloop()