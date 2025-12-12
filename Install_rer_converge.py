import os
from Bio import SeqIO
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed


print("Step 0: Installing RERconverge")
script_dir = os.path.dirname(os.path.abspath(__file__))
if " " in script_dir:
    print(f"Space in location {script_dir}")
    assert False
result = subprocess.run(f"Rscript {script_dir}/Install_rerconverge.r", 
                        shell=True)            # return output as string

# Check if the command succeeded
if result.returncode == 0:
    print("Step 0: Done successfully!")
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)
else:
    print("Step 0: Failed!")
    print("Return code:", result.returncode)
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)
    assert False
