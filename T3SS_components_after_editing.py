import os
import subprocess
import csv
from Bio import SeqIO
from sys import argv

RUN_WITH_CONDA = True
working_directory = argv[1]
bacterial_proteome = os.path.join(working_directory, argv[2])
T3SS_data = os.path.join(working_directory, "T3SS_data")
tmp_mmseqs = os.path.join(working_directory, "tmp_mmseqs")
output_mmseqs = os.path.join(working_directory, "output_mmseqs")
output_file = os.path.join(working_directory, "T3SS.csv")

os.makedirs(tmp_mmseqs, exist_ok=True)
os.makedirs(output_mmseqs, exist_ok=True)
os.makedirs(T3SS_data, exist_ok=True)


def run_mmseqs(query_files_directory):
    for query_file in os.listdir(query_files_directory):
        query_file_path = os.path.join(query_files_directory, query_file)
        output_path = os.path.join(
            output_mmseqs, os.path.basename(query_file) + "_output")
        conda_activate_command = ". ~/miniconda3/etc/profile.d/conda.sh; conda activate test;"
        mmseqs_command = f"mmseqs easy-search {query_file_path} {bacterial_proteome} {output_path} {tmp_mmseqs}"
        run_mmseqs =  conda_activate_command +  mmseqs_command
        if RUN_WITH_CONDA == True:
            subprocess.run(run_mmseqs, shell=True)
        else:
            subprocess.run(mmseqs_command, shell=True)

run_mmseqs(T3SS_data)


def get_T3SS_data_protein_lengths(T3SS_data_file):

    T3SS_data_protein_lengths_dict = {}

    for record in SeqIO.parse(T3SS_data_file, "fasta"):
        protein_name = record.id
        protein_length = len(record.seq)
        T3SS_data_protein_lengths_dict[protein_name] = protein_length

    return T3SS_data_protein_lengths_dict


def get_mmseqs_results_dictionary(mmseqs_results_file, T3SS_data_protein_lengths_dict):
    mmseq_results_dict = {}

    with open(mmseqs_results_file, 'r') as mmseqs_results:
        for line in mmseqs_results:

            columns = line.split()
            T3SS_protein = columns[0]
            bacterial_protein = columns[1]
            alignment_length = float(columns[3])
            e_value = float(columns[-2])
            bit_score = float(columns[-1])

            T3SS_protein_length = T3SS_data_protein_lengths_dict[T3SS_protein]
            alignment_coverage_percentage = alignment_length / T3SS_protein_length

            e_value_cut_off = 1e-10
            alignment_coverage_percentage_cut_off = 0.3

            if e_value < e_value_cut_off and alignment_coverage_percentage > alignment_coverage_percentage_cut_off:
                if T3SS_protein not in mmseq_results_dict:
                    mmseq_results_dict[T3SS_protein] = (
                        bacterial_protein, bit_score)
                elif bit_score > mmseq_results_dict[T3SS_protein][1]:
                    mmseq_results_dict[T3SS_protein] = (
                        bacterial_protein, bit_score)

    return mmseq_results_dict


def get_all_subsystems_dict(T3SS_data, output_mmseqs):

    all_subsystems_dict = {}

    for output_mmseqs_file in os.listdir(output_mmseqs):
        subsystem_name = output_mmseqs_file.split(".")[0]
        mmseqs_results_file = os.path.join(output_mmseqs, output_mmseqs_file)
        T3SS_data_file_name = output_mmseqs_file.replace(
            ".fasta_output", ".fasta")
        T3SS_data_file = os.path.join(
            T3SS_data, T3SS_data_file_name)

        T3SS_data_protein_lengths = get_T3SS_data_protein_lengths(
            T3SS_data_file)

        mmseqs_results_dictionary = get_mmseqs_results_dictionary(
            mmseqs_results_file, T3SS_data_protein_lengths)

        all_subsystems_dict[subsystem_name] = mmseqs_results_dictionary

    return all_subsystems_dict


all_subsystems_dict = get_all_subsystems_dict(T3SS_data, output_mmseqs)


def get_T3SS_homologous_bacterial_genes_list(all_subsystems_dict):
    T3SS_homologous_bacterial_genes = []
    for subsystem_dict in all_subsystems_dict.values():
        for bacterial_T3SS_tuple in subsystem_dict.values():
            bacterial_gene = bacterial_T3SS_tuple[0]
            if bacterial_gene not in T3SS_homologous_bacterial_genes:
                T3SS_homologous_bacterial_genes.append(bacterial_gene)
    return T3SS_homologous_bacterial_genes


def get_best_bacterial_T3SS_match_dict(all_subsystems_dict, T3SS_homologous_bacterial_genes):

    best_bacterial_T3SS_match = {}

    for bacterial_gene in T3SS_homologous_bacterial_genes:
        max_bit_score = 0
        for subsystem_name, subsystem_dict in all_subsystems_dict.items():
            for subsystem_gene, bacterial_gene_bit_score_tup in subsystem_dict.items():
                if bacterial_gene in bacterial_gene_bit_score_tup:
                    bit_score = bacterial_gene_bit_score_tup[1]
                    if bit_score > max_bit_score:
                        max_bit_score = bit_score
                        system_gene = subsystem_gene
                        subsystem = subsystem_name
        best_bacterial_T3SS_match[bacterial_gene] = (system_gene, subsystem)

    return best_bacterial_T3SS_match

T3SS_homologous_bacterial_genes = get_T3SS_homologous_bacterial_genes_list(all_subsystems_dict)
best_bacterial_T3SS_match_dict = get_best_bacterial_T3SS_match_dict(all_subsystems_dict, T3SS_homologous_bacterial_genes)


def get_full_bacterial_T3SS_dict(T3SS_data, best_bacterial_T3SS_match_dict):
    full_bacterial_T3SS_dict = {}
    i = 1
    for T3SS_data_file in os.listdir(T3SS_data):
        subsystem_name = T3SS_data_file.split(".")[0]
        if any(subsystem_name in value for value in best_bacterial_T3SS_match_dict.values()):
            path_to_T3SS_data_file = os.path.join(T3SS_data, T3SS_data_file)
            T3SS_proteins_names = set(
                [rec.id for rec in SeqIO.parse(path_to_T3SS_data_file, "fasta")])
            T3SS_proteins_tup_list = []
            for T3SS_protein in T3SS_proteins_names:
                T3SS_protein_tup = (T3SS_protein, subsystem_name)
                T3SS_proteins_tup_list.append(T3SS_protein_tup)
            for T3SS_protein_tup in T3SS_proteins_tup_list:
                if T3SS_protein_tup in best_bacterial_T3SS_match_dict.values():
                    for key, value in best_bacterial_T3SS_match_dict.items():
                        if T3SS_protein_tup == value:
                            full_bacterial_T3SS_dict[key] = value
                else:
                    full_bacterial_T3SS_dict[int(i)] = T3SS_protein_tup
                    i = i + 1
    return full_bacterial_T3SS_dict

full_bacterial_T3SS_dict = get_full_bacterial_T3SS_dict(T3SS_data, best_bacterial_T3SS_match_dict) 


def write_dict_to_output_file(full_bacterial_T3SS_dict):
    flagella_system = "Flagellar"
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        for bacterial_gene, T3SS_protein_subsystem_tup in full_bacterial_T3SS_dict.items():
            T3SS_protein = T3SS_protein_subsystem_tup[0]
            subsystem = T3SS_protein_subsystem_tup[1]
            if subsystem != flagella_system:
                if isinstance(bacterial_gene, int):
                    writer.writerow([f"{subsystem}_{T3SS_protein}", None])
                else:
                    writer.writerow([f"{subsystem}_{T3SS_protein}", bacterial_gene])
        for bacterial_gene, T3SS_protein_subsystem_tup in full_bacterial_T3SS_dict.items():
            T3SS_protein = T3SS_protein_subsystem_tup[0]
            subsystem = T3SS_protein_subsystem_tup[1]
            if subsystem == flagella_system:
                if isinstance(bacterial_gene, int):
                    writer.writerow([f"{subsystem}_{T3SS_protein}", None])
                else:
                    writer.writerow([f"{subsystem}_{T3SS_protein}", bacterial_gene])

write_dict_to_output_file(full_bacterial_T3SS_dict)
