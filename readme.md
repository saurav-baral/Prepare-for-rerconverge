This script need Biopython, MAFFT and Clipkit to be installed.

If conda is already present:
```bash
conda env create -f environment.yml
```
This creates an environment called bioenv. You can actviate this environment and run the python script there.
```bash
conda activate bioenv
python Prepare_for_rerconverge.py
```
Prepare_for_rerconverge creates MSA using MAFFT and trims it using ClipKit
If you do not have a species tree run Prepare_species_tree.py, which selects 200 random genes with all species present and creates a species tree.
```
python Prepare_for_rerconverge.py
```
