This is a tool to align the PDB files.

## Debug

if you want to debug the code, you can set the DEBUG_MODE to True in the gui.py file

once DEBUG_MODE is set to True, click the "Generate PDB and Config" button will directly generate the PDB files and Config files

## Usage

to run the gui, use the following command:

```
python gui.py
```

1. input the PDB folder path, the generated Configuration folder path
2. input the DNA sequence and the position to align like this: 

      gcta
atcgatcgatcgatcgatcgatcgatcgatcg
                tagc


3. click the "Generate PDB and Config" button

the generated PDB files and Config files will be saved in the corresponding folders

please load the generated Config file in the Configuration folder to the simulation software
