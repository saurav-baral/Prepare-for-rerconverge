
import os
from Bio import SeqIO
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed


location = "/mnt/griffin/saubar/Species_sequences" # location of the root folder where all data and analysis will be done 


list_of_fasta_files = os.listdir(f"{location}/1.Species_sequences") 


'''
"1.Species_sequences" is the folder that contains species specific fasta files with sequences.
This folder is present in the location i.e., inside Species_sequences folder.
This contains files for all the species.
Each fasta file has the following name Species.fas
Each fasta file has the following structure:
>BUSCOid1
Sequence
>BUSCOid2
Sequences

The BUSCO ids are shared between the species.

'''
print("Step 1: Preparing Sequences")
gene_dictionary = {}
gene_name_list = []
species_name_list = []
for file_name in list_of_fasta_files:
    species_name = file_name.split(".")[0]
    if species_name not in species_name_list:
        species_name_list.append(species_name)
    gene_dictionary[species_name] = SeqIO.to_dict(SeqIO.parse(f"{location}/1.Species_sequences/{file_name}", "fasta"))
    for gene_name in gene_dictionary[species_name]:
        if gene_name not in gene_name_list:
            gene_name_list.append(gene_name)

list_of_folders =  os.listdir(f"{location}")
if "2.Gene_sequences" not in list_of_folders:
    os.mkdir(f"{location}/2.Gene_sequences")
else:
    subprocess.run(f"rm -r {location}/2.Gene_sequences/*", shell = True)

for gene_name in gene_name_list:
    output = ''
    for species_name in species_name_list:
        if gene_name in gene_dictionary[species_name]:
            sequence = gene_dictionary[species_name][gene_name].seq
            output += f">{species_name}\n{sequence}\n"
    with open(f"{location}/2.Gene_sequences/{gene_name}.fas", 'w') as out_file:
        out_file.write(output)
print("Step 1: Done")

print("Step 2: Preparing Sequences for Mafft alignment")

list_of_fasta_files = os.listdir(f"{location}/2.Gene_sequences/")

list_of_folders =  os.listdir(f"{location}")
if "3.Sequences_no_gap" not in list_of_folders:
    os.mkdir(f"{location}/3.Sequences_no_gap")
else:
    subprocess.run(f"rm -r {location}/3.Sequences_no_gap/*", shell = True)


if "4.Translated_sequences" not in list_of_folders:
    os.mkdir(f"{location}/4.Translated_sequences")
else:
    subprocess.run(f"rm -r {location}/4.Translated_sequences/*", shell = True)


for i,fasta_file_name in enumerate(list_of_fasta_files):
    if i % 100 == 0:
        print(i,"/",len(list_of_fasta_files))
    output = ''
    output_nogap = ''
    fasta_file = SeqIO.parse(f"{location}/2.Gene_sequences/{fasta_file_name}", 'fasta')
    for records in fasta_file:
        no_gap =  records.seq.replace("-",'')
        translated_sequence = no_gap.translate()
        if translated_sequence[-1] == "*":
            if "*" not in translated_sequence[:-1]:
                output += f">{records.id}\n{translated_sequence[:-1]}\n"
                output_nogap += f">{records.id}\n{no_gap}\n"
        else:
           if "*" not in translated_sequence:
               output += f">{records.id}\n{translated_sequence}\n"
               output_nogap += f">{records.id}\n{no_gap}\n"
               
    with open(f'{location}/4.Translated_sequences/{fasta_file_name}', 'w') as out_file:
        out_file.write(output)
    with open(f'{location}/3.Sequences_no_gap/{fasta_file_name}', 'w') as out_file:
        out_file.write(output_nogap)


print("Step 2: Done")


print("Step 3: Aligning")
list_of_folders =  os.listdir(f"{location}")
if "5.Aligned_MAFFT_protein" not in list_of_folders:
    os.mkdir(f"{location}/5.Aligned_MAFFT_protein")
else:
    subprocess.run(f"rm -r {location}/5.Aligned_MAFFT_protein/*", shell = True)

list_of_fasta_files = os.listdir(f"{location}/4.Translated_sequences")
list_of_done_files = os.listdir(f"{location}/5.Aligned_MAFFT_protein")



def run_mafft(file_name):
    cmd = (
        f"mafft --auto --thread 1 "
        f"'{location}/4.Translated_sequences/{file_name}' > "
        f"'{location}/5.Aligned_MAFFT_protein/{file_name}'"
    )
    subprocess.run(cmd, shell=True)
    return file_name

with ThreadPoolExecutor(max_workers=10) as executor:
    futures = [
        executor.submit(run_mafft, f)
        for f in list_of_fasta_files
        if f not in list_of_done_files
    ]

    for i, future in enumerate(as_completed(futures), 1):
        print(f"{i}/{len(futures)} completed: {future.result()}")


def map_gaps_to_untranslated(translated_with_gaps, untranslated_no_gaps):
    """
    Map the gaps in a translated sequence onto its corresponding untranslated sequence.

    Parameters
    ----------
    translated_with_gaps : str
        Translated sequence containing '-' for gaps (1 char per amino acid/codon).
    untranslated_no_gaps : str
        Untranslated nucleotide sequence without any gaps.

    Returns
    -------
    str
        Untranslated nucleotide sequence with gaps inserted to match the translated alignment.
    """
    untranslated_with_gaps = []
    nt_index = 0  # position in the nucleotide sequence

    for aa in translated_with_gaps:
        if aa == "-":
            untranslated_with_gaps.append("---")  # gap at codon level
        else:
            codon = untranslated_no_gaps[nt_index:nt_index+3]
            untranslated_with_gaps.append(codon)
            nt_index += 3

    return "".join(untranslated_with_gaps)



list_of_untranslated_files = os.listdir(f'{location}/3.Sequences_no_gap')

list_of_aligned_files = os.listdir(f'{location}/5.Aligned_MAFFT_protein')

list_of_folders =  os.listdir(f"{location}")
if "6.Aligned_MAFFT_untranslated" not in list_of_folders:
    os.mkdir(f"{location}/6.Aligned_MAFFT_untranslated")
else:
    subprocess.run(f"rm -r {location}/6.Aligned_MAFFT_untranslated/*", shell = True)

for files in list_of_aligned_files:
    aligned_file_dictionary = SeqIO.to_dict(SeqIO.parse(f'{location}/5.Aligned_MAFFT_protein/{files}', "fasta"))
    untrans_file_dictionbary = SeqIO.to_dict(SeqIO.parse(f'{location}/3.Sequences_no_gap/{files}', "fasta"))
    output = ''
    for sequence_name in aligned_file_dictionary:
        aligned_sequence = aligned_file_dictionary[sequence_name].seq
        untrans_sequence =  untrans_file_dictionbary[sequence_name].seq
        # print(aligned_sequence)
        output += f">{sequence_name}\n{(map_gaps_to_untranslated(str(aligned_sequence), str(untrans_sequence)))}\n"
    with open(f'{location}/6.Aligned_MAFFT_untranslated/{files}', 'w') as out_file:
        out_file.write(output)

print("Step 3: Done")

print(f"Step 4: Counting species per gene and removing those with < {int(len(species_name_list)/2)} missing")
list_of_files = os.listdir(f'{location}/6.Aligned_MAFFT_untranslated')
output = "Sequence,Count"
for i,files in enumerate(list_of_files):
    if i % 100 == 0:
        print(i, "Done")
    aligned_file_dictionary = SeqIO.to_dict(SeqIO.parse(f'{location}/6.Aligned_MAFFT_untranslated/{files}', "fasta"))
    output += f"\n{files},{len(aligned_file_dictionary)}"
    if len(aligned_file_dictionary) < len(species_name_list)/2:
        subprocess.run(f"rm {location}/6.Aligned_MAFFT_untranslated/{files}", shell = True)

with open(f'{location}/counts.csv', 'w') as out_file:
    out_file.write(output)

print("Step 4: Done")

print(f"Step 5: Trimming with Clipkit")

list_of_folders =  os.listdir(f"{location}")
if "7.Clipkit_trimmed" not in list_of_folders:
    os.mkdir(f"{location}/7.Clipkit_trimmed")
else:
    subprocess.run(f"rm -r {location}/7.Clipkit_trimmed/*", shell = True)


list_of_files = os.listdir(f"{location}/6.Aligned_MAFFT_untranslated")

def run_clipkit(file_name):
    cmd = (
        f"clipkit '{location}/6.Aligned_MAFFT_untranslated/{file_name}' "
        f"--codon -o '{location}/7.Clipkit_trimmed/{file_name}'"
    )
    subprocess.run(cmd, shell=True)
    return file_name

with ThreadPoolExecutor(max_workers=15) as executor:
    futures = [executor.submit(run_clipkit, f) for f in list_of_files]

    for i, future in enumerate(as_completed(futures), 1):
        print(f"{i}/{len(futures)} completed: {future.result()}")

print(f"Step 5: Done")