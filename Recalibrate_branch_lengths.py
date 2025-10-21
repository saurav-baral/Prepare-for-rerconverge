import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

location = "/mnt/griffin/saubar/Species_sequences"


print(f"Step 7: Recalibrating Branch lengths")

list_of_folders =  os.listdir(f"{location}")
if "9.For_recalibration" not in list_of_folders:
    os.mkdir(f"{location}/9.For_recalibration")
else:
    subprocess.run(f"rm -r {location}/9.For_recalibration/*", shell = True)


list_of_sequence = os.listdir(f"{location}/7.Clipkit_trimmed/")

counter = 1
seq_counter = 0

tree_folders = []

if len(list_of_sequence) > 100:
    num_gene = len(list_of_sequence)/10
else:
    num_gene = 100
for i,seq_file in enumerate(list_of_sequence):
    seq_counter += 1
    if seq_counter >= num_gene:
        counter += 1
        seq_counter = 0
    target_location = f"tree_{counter}"
    if target_location not in tree_folders:
        tree_folders.append(target_location)
    list_of_tree_folders =  os.listdir(f"{location}/9.For_recalibration/")
    if target_location not in list_of_tree_folders:
        os.mkdir(f"{location}/9.For_recalibration/{target_location}")
    

    subprocess.run(f"cp '{location}/7.Clipkit_trimmed/{seq_file}' '{location}/9.For_recalibration/{target_location}/{seq_file}'", shell = True)

subprocess.run (f"cp '{location}/species_tree.tre' '{location}/9.For_recalibration/'", shell = True)



def run_rscript(folder_name):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    result = subprocess.run(f"Rscript {script_dir}/Tree_recalibrate.r {folder_name} {location}/9.For_recalibration/", 
                        shell=True, 
                        capture_output=True,  # capture stdout/stderr
                        text=True)            # return output as string
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)
    return folder_name

with ThreadPoolExecutor(max_workers=10) as executor:
    futures = [executor.submit(run_rscript, f) for f in tree_folders]

    for i, future in enumerate(as_completed(futures), 1):
        print(f"{i}/{len(futures)} completed: {future.result()}")

list_of_folders =  os.listdir(f"{location}")
if "10.Run_rerconverge" not in list_of_folders:
    os.mkdir(f"{location}/10.Run_rerconverge")
else:
    subprocess.run(f"rm -r {location}/10.Run_rerconverge/*", shell = True)
subprocess.run(f"cd {location}/9.For_recalibration/\ncat *.recal.tre > ../10.Run_rerconverge/recalibrated_trees.tree", shell = True)

print(f"Step 7: Done")

print(f"Step 8: Running Rerconverge")

list_of_species = [species_file_name.split(".")[0] for species_file_name in os.listdir(f"{location}/1.Species_sequences")]
species_phenotype = {}

with open(f"{location}/phenotype.csv", 'r') as pheno_file:
    for lines in pheno_file.readlines()[1:]:
        lines_split = lines.strip().split(",")
        species_phenotype[lines_split[0]]= lines_split[1]
foreground_species = []
for species_name in list_of_species:
    if species_name not in species_phenotype:
        print(f"{species_name} missing from phenotype data")
        assert False
    else:
        if species_phenotype[species_name] == "1":
            foreground_species.append(species_name)
foreground_species = ",".join(foreground_species)


script_dir = os.path.dirname(os.path.abspath(__file__))
result = subprocess.run(f"Rscript {script_dir}/Run_rerconverge_binary.r {foreground_species} {location}/10.Run_rerconverge/", 
                        shell=True, 
                        capture_output=True,  # capture stdout/stderr
                        text=True)            # return output as string
print("STDOUT:", result.stdout)
print("STDERR:", result.stderr)