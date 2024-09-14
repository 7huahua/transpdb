import tkinter as tk
from tkinter import messagebox
import subprocess
from compute_config import *
from generate_pdb import *

# 假设这是你的DNA对齐脚本函数
def align_dna(folder, sequence1, sequence2, seuqence3):
    # 简单示例对齐逻辑
    # 实际逻辑需要根据你的需求实现
    config_path = compute_config(folder, sequence1, sequence2, seuqence3)
    return config_path

# 点击按钮时调用的函数
def run_alignment():
    folder = entry_folder.get()
    seq1 = entry_seq1.get()
    seq2 = entry_seq2.get()
    seq3 = entry_seq3.get()
    
    if not seq1 or not seq2 or not seq3:
        messagebox.showwarning("Fault", "Please input the DNA PDB file path")
        return
    
    result = align_dna(folder, seq1, seq2, seq3)
    result_label.config(text=f"Result Path: {result}")

# 解析DNA序列的函数
def parse_dna_sequence(input_text):
    """
    解析用户输入的DNA序列，记录每行的空格数量和序列信息
    :param input_text: 用户输入的多行序列文本
    :return: 包含空格数量和序列信息的列表
    """
    sequence_data = []
    lines = input_text.splitlines()

    for line in lines:
        # 计算每行开头的空格数量
        leading_spaces = len(line) - len(line.lstrip(' '))
        # 提取有效的DNA序列（A, T, C, G）
        dna_sequence = line.strip()

        # 只处理包含有效DNA字符的行
        if dna_sequence:
            sequence_data.append({
                'spaces': leading_spaces,
                'sequence': dna_sequence
            })

    return sequence_data

# # 点击解析按钮时调用的函数
# def run_parsing():
#     user_input = text_input.get("1.0", tk.END)
#     if not user_input.strip():
#         messagebox.showwarning("输入错误", "请输入DNA序列")
#         return
    
#     # 调用解析函数
#     parsed_sequences = parse_dna_sequence(user_input)

#     # 清空输出框
#     output_text.delete("1.0", tk.END)
    
#     # 显示解析结果
#     for item in parsed_sequences:
#         output_text.insert(tk.END, f"空格数量: {item['spaces']}, 序列: {item['sequence']}\n")


# 解析DNA序列的函数
def parse_dna_sequence(input_text):
    """
    解析用户输入的DNA序列，记录每行的空格数量和序列信息
    :param input_text: 用户输入的多行序列文本
    :return: 包含空格数量和序列信息的字符串
    """
    sequence_data = []
    lines = input_text.splitlines()

    for line in lines:
        # 计算每行开头的空格数量
        leading_spaces = len(line) - len(line.lstrip(' '))
        # 提取有效的DNA序列（A, T, C, G）
        dna_sequence = line.strip()

        # 只处理包含有效DNA字符的行
        if dna_sequence:
            # 保留空格并记录序列
            sequence_data.append(' ' * leading_spaces + dna_sequence)

    # 返回包含空格信息的完整序列
    return '\n'.join(sequence_data)

# 点击解析并生成PDB按钮时调用的函数
def run_generate_pdb():
    user_input = text_input.get("1.0", tk.END)
    if not user_input.strip():
        messagebox.showwarning("输入错误", "请输入DNA序列")
        return
    
    # 调用解析函数，获取包含空格的DNA序列
    parsed_sequences = parse_dna_sequence(user_input)

    for index, sequence in enumerate(parsed_sequences):

        # 将解析出的序列作为参数传递给generate_pdb.py
        try:
            result = subprocess.run(
                ['python', 'generate_pdb.py', '--dna_seq', sequence],
                capture_output=True,
                text=True
            )
            
            if result.returncode == 0:
                output_text.delete("1.0", tk.END)
                output_text.insert(tk.END, "PDB文件生成成功！\n")
                output_text.insert(tk.END, result.stdout)
            else:
                output_text.delete("1.0", tk.END)
                output_text.insert(tk.END, "生成失败，请检查！\n")
                output_text.insert(tk.END, result.stderr)
        
        except Exception as e:
            messagebox.showerror("错误", f"运行generate_pdb.py时出现问题: {e}")



# 创建主窗口
root = tk.Tk()
root.title("DNA Alignment Tool")

# 创建并放置标签和文本框
folder_path_label = tk.Label(root, text="Folder Path:")
folder_path_label.pack()

entry_folder = tk.Entry(root, width=50)
entry_folder.pack()

label_seq1 = tk.Label(root, text="DNA PDB file path 1:")
label_seq1.pack()

entry_seq1 = tk.Entry(root, width=50)
entry_seq1.pack()

label_seq2 = tk.Label(root, text="DNA PDB file path 2:")
label_seq2.pack()

entry_seq2 = tk.Entry(root, width=50)
entry_seq2.pack()

label_seq3 = tk.Label(root, text="DNA PDB file path 3:")
label_seq3.pack()

entry_seq3 = tk.Entry(root, width=50)
entry_seq3.pack()

# 创建并放置按钮
run_button = tk.Button(root, text="RUN", command=run_alignment)
run_button.pack()

# 创建并放置结果标签
result_label = tk.Label(root, text="Result Path:")
result_label.pack()

# -----------------------------------------------------

# 创建文本框用于输入DNA序列
label_input = tk.Label(root, text="请输入DNA序列:")
label_input.pack()

text_input = tk.Text(root, height=10, width=50)
text_input.pack()

# 创建并放置解析按钮
parse_button = tk.Button(root, text="解析序列", command=run_generate_pdb)
parse_button.pack()

# 创建文本框用于显示解析结果
label_output = tk.Label(root, text="解析结果:")
label_output.pack()

output_text = tk.Text(root, height=10, width=50)
output_text.pack()



# 启动主事件循环
root.mainloop()
