print(f"Step 6: Make Tree with random 200 genes")

import os
from Bio import SeqIO
import subprocess
import random

location = "/mnt/griffin/saubar/Species_sequences" # location of the root folder where all data and analysis will be done 

list_of_folders =  os.listdir(f"{location}")
if "8.For_tree" not in list_of_folders:
    os.mkdir(f"{location}/8.For_tree")
else:
    subprocess.run(f"rm -r {location}/8.For_tree/*", shell = True)


species_dict = {}

with open(f'{location}/counts.csv', 'r') as count_file:
    count_lines = count_file.readlines()[1:]

max_species = 0
for lines in count_lines:
    lines_split = lines.strip().split(",")
    if int(lines_split[1]) > max_species:
            max_species = int(lines_split[1])
gene_to_be_used = []
for lines in count_lines:
    lines_split = lines.strip().split(",")
    if int(lines_split[1]) == max_species:
           gene_to_be_used.append(lines_split[0]) 

count_gene = 0  
for gene_name in   gene_to_be_used:
                  
    if count_gene > 200:
        break
    count_gene += 1
    
    # print(lines_split[0])
    fasta_file = SeqIO.parse(f"{location}/7.Clipkit_trimmed/{random.choice(gene_to_be_used)}", "fasta")
    for records in fasta_file:
        species_dict.setdefault(records.id, '')
        species_dict[records.id] = species_dict[records.id] + records.seq
output = ''
for key, value in species_dict.items():
    output += f">{key}\n{value}\n\n"

with open(f"{location}/8.For_tree/merged_{max_species}sp.fas", 'w') as out_file:
    out_file.write(output)


subprocess.run(f"cd '{location}/8.For_tree'\niqtree -s merged_{max_species}sp.fas -m MFP -nt AUTO", shell = True)

subprocess.run (f"cd '{location}/8.For_tree'\ncp merged_{max_species}sp.fas.treefile ../species_tree.tre", shell = True)

print(f"Step 6: Done")