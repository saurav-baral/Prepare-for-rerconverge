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