import argparse
import os
import subprocess
import pandas as pd
from Bio import SeqIO

RUN_WITH_CONDA = True
E_VALUE_CUT_OFF = 1e-10
QUERY_COVERAGE_PERCENTAGE_CUT_OFF = 0.3
FLAGELLA_SYSTEM = "Flagellar"
MMSEQS_OUTPUT_FORMAT = 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,qcov,tstart,tend,evalue,bits'

# DONE: Don't put code in the global scope, put it in a main function

# DONE: I would suggest to use argparse instead of sys.argv

def run_mmseqs(query_files_directory, output_mmseqs, bacterial_proteome, tmp_mmseqs):
    for query_file in os.listdir(query_files_directory):
        query_file_path = os.path.join(query_files_directory, query_file)
        # DONE: Change mmseqs output file name to be with .m8 extension.
        output_path = os.path.join(
            output_mmseqs, os.path.basename(query_file).split(".")[0] + ".m8")
        # DONE: Move the lines relevant only for conda to the if statement
        mmseqs_command = f"mmseqs easy-search {query_file_path} {bacterial_proteome} {output_path} {tmp_mmseqs} --format-output {MMSEQS_OUTPUT_FORMAT}"
        if RUN_WITH_CONDA == True:
            conda_activate_command = ". ~/miniconda3/etc/profile.d/conda.sh; conda activate test;"
            run_mmseqs =  conda_activate_command +  mmseqs_command
            subprocess.run(run_mmseqs, shell=True)
        else:
            subprocess.run(mmseqs_command, shell=True)


def get_mmseqs_results_dictionary(mmseqs_results_file):
    # DONE: Add comment with example of the dict content
    mmseq_results_dict = {} # {T3SS_protein: (bacterial_protein, bit_score)}

    # DONE: Use pandas to read and filter the mmseqs result
    df = pd.read_csv(mmseqs_results_file, sep='\t', header=None)
    
    df.columns = ['T3SS_protein', 'bacterial_protein', "identity_percent", "alignment length", "mismatch", "gapopen", "query_start", "query_end", "query_ length", 'query_coverage_percentage', "target_ start", "target_end", 'e_value', 'bit_score']
    
    # DONE: rename to query coverage percentage
    # DONE: You can ask mmseqs to output the query coverage percentage, and then you don't need
    #     to calculate it yourself (and you don't need the merhod get_T3SS_data_protein_lengths)
    # DONE: Use constants for the cutoff values

    filtered_df = df[(df['e_value'] < E_VALUE_CUT_OFF) & (df['query_coverage_percentage'] > QUERY_COVERAGE_PERCENTAGE_CUT_OFF)]
    
    for row in filtered_df.iterrows():
        T3SS_protein = row[1]['T3SS_protein']
        bacterial_protein = row[1]['bacterial_protein']
        bit_score = row[1]['bit_score']
        
        if T3SS_protein not in mmseq_results_dict:
            mmseq_results_dict[T3SS_protein] = (bacterial_protein, bit_score)
        elif bit_score > mmseq_results_dict[T3SS_protein][1]:
            mmseq_results_dict[T3SS_protein] = (bacterial_protein, bit_score)
    
    return mmseq_results_dict

def get_all_subsystems_dict(output_mmseqs):

    all_subsystems_dict = {} # {subsystem_name: {T3SS_protein: (bacterial_protein, bit_score)}}

    for output_mmseqs_file in os.listdir(output_mmseqs):
        subsystem_name = output_mmseqs_file.split(".")[0]
        mmseqs_results_file = os.path.join(output_mmseqs, output_mmseqs_file)
        T3SS_data_file_name = output_mmseqs_file.replace(
            ".m8", ".fasta")

        mmseqs_results_dictionary = get_mmseqs_results_dictionary(
            mmseqs_results_file)

        all_subsystems_dict[subsystem_name] = mmseqs_results_dictionary

    return all_subsystems_dict


def get_T3SS_homologous_bacterial_genes_list(all_subsystems_dict):
    T3SS_homologous_bacterial_genes = []
    for subsystem_dict in all_subsystems_dict.values():
        for bacterial_T3SS_tuple in subsystem_dict.values():
            bacterial_gene = bacterial_T3SS_tuple[0]
            if bacterial_gene not in T3SS_homologous_bacterial_genes:
                T3SS_homologous_bacterial_genes.append(bacterial_gene)
    return T3SS_homologous_bacterial_genes


def get_best_bacterial_T3SS_match_dict(all_subsystems_dict):

    T3SS_homologous_bacterial_genes = get_T3SS_homologous_bacterial_genes_list(all_subsystems_dict)

    best_bacterial_T3SS_match = {} # {bacterial_gene: (system_gene, subsystem)}

    for bacterial_gene in T3SS_homologous_bacterial_genes:
        max_bit_score = 0
        for subsystem_name, subsystem_dict in all_subsystems_dict.items():
            for subsystem_gene, bacterial_gene_bit_score_tup in subsystem_dict.items():
                if bacterial_gene in bacterial_gene_bit_score_tup:
                    bit_score = bacterial_gene_bit_score_tup[1]
                    if bit_score > max_bit_score:
                        max_bit_score = bit_score
                        # DONE: Use more informative names ("best_system_gene")
                        best_system_gene = subsystem_gene
                        best_subsystem = subsystem_name
        best_bacterial_T3SS_match[bacterial_gene] = (best_system_gene, best_subsystem)

    return best_bacterial_T3SS_match


def get_full_bacterial_T3SS_dict(T3SS_data, best_bacterial_T3SS_match_dict):
    full_bacterial_T3SS_dict = {} # {bacterial_gene: (system_gene, subsystem), number: (system_gene, subsystem)}
    i = 1
    for T3SS_data_file in os.listdir(T3SS_data):
        subsystem_name = T3SS_data_file.split(".")[0]
        if any(subsystem_name in value for value in best_bacterial_T3SS_match_dict.values()):
            path_to_T3SS_data_file = os.path.join(T3SS_data, T3SS_data_file)
            T3SS_proteins_names = set(
                [rec.id for rec in SeqIO.parse(path_to_T3SS_data_file, "fasta")])
            # DONE: Delete the list T3SS_proteins_tup_list and the for-loop after it.
            for T3SS_protein in T3SS_proteins_names:
                if (T3SS_protein, subsystem_name) in best_bacterial_T3SS_match_dict.values():
                    for key, value in best_bacterial_T3SS_match_dict.items():
                        if (T3SS_protein, subsystem_name) == value:
                            full_bacterial_T3SS_dict[key] = value
                else:
                    full_bacterial_T3SS_dict[int(i)] = (T3SS_protein,subsystem_name)
                    i = i + 1
    return full_bacterial_T3SS_dict


# DONE: Consider adding a header line to the output csv that describes the column.
# DONE: Conside using pandas to write the csv instead of csv.writer
def write_dict_to_output_file(full_bacterial_T3SS_dict, output_file):
    headers = ["Subsystem_T3SS Protein", "Bacterial Protein ID"]
    data = []

    for bacterial_gene, (T3SS_protein, subsystem) in full_bacterial_T3SS_dict.items():
        if subsystem != FLAGELLA_SYSTEM:
            if isinstance(bacterial_gene, int):
                    data.append([f"{subsystem}_{T3SS_protein}", None])
            else:
                    data.append([f"{subsystem}_{T3SS_protein}", bacterial_gene])
        if subsystem == FLAGELLA_SYSTEM:
            if isinstance(bacterial_gene, int):
                    data.append([f"{subsystem}_{T3SS_protein}", None])
            else:
                    data.append([f"{subsystem}_{T3SS_protein}", bacterial_gene])
    
    df = pd.DataFrame(data, columns=headers)
    
    df.to_csv(output_file, index=False)

def main():
    parser = argparse.ArgumentParser(description='Run mmseqs and process mmseqs results files')
    parser.add_argument('working_directory', type=str, help='Path to the working directory')
    parser.add_argument('bacterial_proteome', type=str, help='Name of the bacterial proteome file')

    args = parser.parse_args()

    working_directory = args.working_directory
    bacterial_proteome = os.path.join(working_directory, args.bacterial_proteome)
    
    T3SS_data = os.path.join(working_directory, "T3SS_data")
    tmp_mmseqs = os.path.join(working_directory, "tmp_mmseqs")
    output_mmseqs = os.path.join(working_directory, "output_mmseqs")
    output_file = os.path.join(working_directory, "T3SS.csv")

    os.makedirs(tmp_mmseqs, exist_ok=True)
    os.makedirs(output_mmseqs, exist_ok=True)
    os.makedirs(T3SS_data, exist_ok=True)

    run_mmseqs(T3SS_data, output_mmseqs, bacterial_proteome, tmp_mmseqs)
    all_subsystems_dict = get_all_subsystems_dict(output_mmseqs)
    # DONE: Put the next line inside get_best_bacterial_T3SS_match_dict
    best_bacterial_T3SS_match_dict = get_best_bacterial_T3SS_match_dict(all_subsystems_dict)
    full_bacterial_T3SS_dict = get_full_bacterial_T3SS_dict(T3SS_data, best_bacterial_T3SS_match_dict) 
    write_dict_to_output_file(full_bacterial_T3SS_dict, output_file)

if __name__ == "__main__":
    main()
