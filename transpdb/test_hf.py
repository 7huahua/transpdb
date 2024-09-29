import huggingface_hub


# read the pdb file as a string
with open("D:/2024.02.22_VRSim-Gutmann/SimulationFolder/PDB_Files/nucleotide_A.pdb", "r") as f:
    pdb_file = f.read()

# import transformers and llama tokenizers
from transformers import LlamaTokenizer

# load the llama tokenizer
tokenizer = LlamaTokenizer.from_pretrained("meta-llama/Llama-2-7b-hf")

# tokenize the pdb file
inputs = tokenizer(pdb_file, return_tensors="pt")
print(len(inputs["input_ids"]))
