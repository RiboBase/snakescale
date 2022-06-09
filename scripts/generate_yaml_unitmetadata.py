import argparse
import os
import pandas as pd
import numpy as np
import sqlite3
import yaml 
from pprint import pprint

###############################################################################
#                          GENERATE YAML FRO RIBOFLOW
###############################################################################
# This script generates yaml files for the RiboFlow pipeline
# It takes a study (or a list of studies) and obtains study unitmetadata
# from the sqlite file of RiboBase
# It outputs a yaml file to be used for this study
#
# Sample Run 1:
#     python generate_yaml.py --study GSE101760_dedup --template ../project.yaml \
#                             --output hebele --db ../db/db.sqlite3

# Sample Run 2:
#    python generate_yaml.py --text stuy_list.txt --template ../project.yaml \
#                            --output hebele --db ../db/db.sqlite3
# This run gets the list of studies from the text file stuy_list.txt
# where each study is given in a separate line.
#
# Dedup is a special keyword indicating RiboFlow to perform deduplication 
# See Sample Run 1 as an example.
##################################################################################

# helper method to assemble the clip sequence
def generate_clip_sequence(clipping_param_base, experiment_dict, experiment_type):
    # extract adapters across experiments
    cur_threep_set = set(experiment_dict[exp]['threep_adapter'] for exp in experiment_dict)
    cur_fivep_set  = set(experiment_dict[exp]['fivep_adapter'] for exp in experiment_dict)
    
    # check for the case of multiple adapters
    if len(cur_threep_set) > 1:
        print("ERROR - Multiple Different Adapters within the + " + experiment_type + " Experiments")
        print("The following adapters were found: ")
        for val in cur_threep_set:
            print(val)
    clip_sequence = clipping_param_base
    
    # compute the overlap required in the case of 'N' within 5' end of 3' adapter
    if len(cur_threep_set) >= 1:
        candidate_adapter = list(cur_threep_set)[0]
        overlap           = 4
        index             = 0
        reached_end       = False

        while index < len(candidate_adapter) and not reached_end:
            if not reached_end and candidate_adapter[index] == 'N':
                overlap += 1
            elif candidate_adapter[index] != 'N':
                reached_end = True

            index += 1
        # having an adapter requires additional parameters 
        
        if len(candidate_adapter) > 0:
            clip_sequence = clip_sequence + " -a " + candidate_adapter + " --overlap=" + str(overlap)
            
            if experiment_type == 'Ribo-Seq':
                clip_sequence = clip_sequence + " --trimmed-only"

    return clip_sequence


def generate_yaml(study, template, output, db, download_path, 
                  reference_file    = "scripts/references.yaml", 
                  reference_folder = "reference"):

    conn      = sqlite3.connect(db)
    
    dedup_val = False
    gse_only  = study
    
    """
    # Previous version
    # We made this part a bit more verbose when it fails
    if '_dedup' in study:
        gse_only  = study.split("_")[0]
        dedup_val = True
    """

    study_contents = study.split("_")

    if len(study_contents) > 1:
        gse_only = study_contents[0]

        if study_contents[1] == "dedup":
            dedup_val = True
        elif study_contents[1] != "test":
            raise RuntimeError("Invalid run type {} for the study {}".format(study_contents[1], study_contents[0]) )

    study_query = 'SELECT "unitmetadata_study"."id", "unitmetadata_study"."creation_date",'\
                  '"unitmetadata_study"."notes", "unitmetadata_study"."unitmetadata_checked", "unitmetadata_study"."modifier",'\
                  '"unitmetadata_study"."geo_accession", "unitmetadata_study"."study_accession", "unitmetadata_study"."study_title",'\
                  '"unitmetadata_study"."study_type", "unitmetadata_study"."study_abstract", "unitmetadata_study"."study_description",'\
                  '"unitmetadata_study"."xref_link", "unitmetadata_study"."submission_accession", "unitmetadata_study"."sradb_updated",'\
                  '"unitmetadata_study"."soft_deleted" FROM "unitmetadata_study" WHERE "unitmetadata_study"."geo_accession" = "{study_gse}"'.\
                      format(study_gse = gse_only)

    study_df = pd.read_sql_query(study_query, conn)

    if study_df.empty:
        raise Exception("Study was not found in the database.")

    study_id  = study_df['id'].values[0]
    ribo_yaml = {}
    yaml_name = study

    exp_query = '''SELECT "unitmetadata_experiment"."id", "unitmetadata_experiment"."creation_date", "unitmetadata_experiment"."notes", "unitmetadata_experiment"."unitmetadata_checked",\
                "unitmetadata_experiment"."modifier", "unitmetadata_experiment"."study_id", "unitmetadata_experiment"."matched_experiment_id", "unitmetadata_experiment"."experiment_alias",\
                "unitmetadata_experiment"."experiment_accession", "unitmetadata_experiment"."type", "unitmetadata_experiment"."title", "unitmetadata_experiment"."study_name",\
                "unitmetadata_experiment"."design_description", "unitmetadata_experiment"."sample_accession", "unitmetadata_experiment"."sample_attribute", "unitmetadata_experiment"."library_strategy",\
                "unitmetadata_experiment"."library_layout", "unitmetadata_experiment"."library_construction_protocol", "unitmetadata_experiment"."platform", "unitmetadata_experiment"."platform_parameters",\
                "unitmetadata_experiment"."xref_link", "unitmetadata_experiment"."experiment_attribute", "unitmetadata_experiment"."submission_accession", "unitmetadata_experiment"."sradb_updated",\
                "unitmetadata_experiment"."organism", "unitmetadata_experiment"."cell_line", "unitmetadata_experiment"."group", "unitmetadata_experiment"."threep_adapter", "unitmetadata_experiment"."fivep_adapter",\
                "unitmetadata_experiment"."threep_umi_length", "unitmetadata_experiment"."fivep_umi_length", "unitmetadata_experiment"."read_length", "unitmetadata_experiment"."is_paired_end", "unitmetadata_experiment"."experiment_file"
                FROM "unitmetadata_experiment" WHERE "unitmetadata_experiment"."study_id" = {id}'''\
                    .format(id = study_id)

    experiment_df = pd.read_sql_query(exp_query, conn)
    rnaseq_df     = experiment_df[experiment_df['type'] == 'RNA-Seq']
    riboseq_df    = experiment_df[experiment_df['type'] == 'Ribo-Seq']

    with open(template) as file:
        ribo_yaml = yaml.load(file, Loader=yaml.FullLoader)

    print("Currently Generating YAML File for: " + yaml_name)
    
    # Grab the study entry and its experiments
    ribo_dict     = {} # ribo id only
    ribo_rna_dict = {} # maps ribo experiment to the entire matched_experiment_id
    rna_dict      = {} # rna id only

    # For each accession number, link the GSM to a ribosome profiling experiment in ribo_id
    # Then, link the GSM to the matching RNA-Seq experiment in ribo_rna_id
    # Note that many experiments have additional RNA-Seq experiments that are not matched to a Ribo-Seq experiment,
    # but this is not of major importance 
    # first iterate over all RNA-Seq experiments, a subset of these will be matched to Ribo-Seq

    for index, row in riboseq_df.iterrows():
        cur_ribo_gsm            = row['experiment_alias']
        ribo_dict[cur_ribo_gsm] = row.to_dict()

        matched_rnaseq_id = row['matched_experiment_id']
        matched_df        = rnaseq_df[rnaseq_df['id'] == matched_rnaseq_id]
        
        for rna_index, rna_row in matched_df.iterrows():
            cur_rna_gsm                 = rna_row['experiment_alias']
            rna_row_dict                = rna_row.to_dict()
            ribo_rna_dict[cur_ribo_gsm] = rna_row_dict
            rna_dict[cur_rna_gsm]       = rna_row_dict

    # Generate the clip arugments
    ribo_clip_argument_base               = '-u 1 --maximum-length=40 --minimum-length=15 --quality-cutoff=28'
    ribo_yaml['clip_arguments']           = generate_clip_sequence(ribo_clip_argument_base, ribo_dict, 'Ribo-Seq')
    rna_clip_argument_base                = '-u 5 -l 40 --quality-cutoff=28'
    ribo_yaml['rnaseq']['clip_arguments'] = generate_clip_sequence(rna_clip_argument_base, ribo_rna_dict, 'RNA-Seq')

    organism_set = set(exp['organism'] for exp in ribo_dict.values())

    if len(organism_set) == 0:
        raise Exception("No organism detected")

    cur_organism = list(organism_set)[0].lower()

    # We need to insert a new clause for each additional organism
    # to indicate its reference files
    if len(organism_set) != 1:
        print("More than one organism detected.")
        raise Exception("More than one organism detected")

    with open( reference_file, 'rt' ) as yaml_stream:
        yaml_contents = yaml.load(yaml_stream, Loader=yaml.FullLoader)

    if cur_organism not in yaml_contents.keys():
        raise Exception("The organism {} is not in the set of supported organisms: ", list( yaml_contents.keys() ))

    # Read the reference paths from the dictionary coming from the reference yaml file
    ref_contents = yaml_contents[cur_organism]

    for ref_key in ["filter", "regions", "transcript_lengths", "transcriptome"]:
        ribo_yaml['input']['reference'][ref_key] = os.path.join(reference_folder , ref_contents[ref_key])

    #   OLD CODE:  TO BE DELETED
    # elif cur_organism.lower() == 'homo sapiens':
    #     # Homo sapiens is the default organism
    #     pass
    # elif cur_organism == 'mus musculus':
    #     ribo_yaml['input']['reference']['filter']             = './reference/filter/mouse/mouse_rtRNA*'
    #     ribo_yaml['input']['reference']['regions']            = './reference/transcriptome/mouse/appris_mouse_v2_actual_regions.bed'
    #     ribo_yaml['input']['reference']['transcript_lengths'] = './reference/transcriptome/mouse/appris_mouse_v2_transcript_lengths.tsv'
    #     ribo_yaml['input']['reference']['transcriptome']      = './reference/transcriptome/mouse/appris_mouse_v2_selected*'
    # elif cur_organism == 'arabidopsis thaliana':
    #     ribo_yaml['input']['reference']['filter']             = './reference/filter/arabidopsis/trna_rrna_seqs.fa*'
    #     ribo_yaml['input']['reference']['regions']            = './reference/transcriptome/arabidopsis/mrna_regions_unique_genes.bed'
    #     ribo_yaml['input']['reference']['transcript_lengths'] = './reference/transcriptome/arabidopsis/mrna_lengths_unique_genes.tsv'
    #     ribo_yaml['input']['reference']['transcriptome']      = './reference/transcriptome/arabidopsis/mrna_seqs_unique_genes.fa*'
    # elif cur_organism == 'caenorhabditis elegans':
    #     ribo_yaml['input']['reference']['filter']             = './reference/filter/celegans/celegans_rRNA_new*'
    #     ribo_yaml['input']['reference']['regions']            = './reference/transcriptome/celegans/appris_celegans_v1_actual_regions_new.bed'
    #     ribo_yaml['input']['reference']['transcript_lengths'] = './reference/transcriptome/celegans/appris_celegans_v1_transcript_lengths_new.tsv'
    #     ribo_yaml['input']['reference']['transcriptome']      = './reference/transcriptome/celegans/appris_celegans_v1_new*'
    # else:
    #     error_msg = "The organism {} is not supported".format( cur_organism )
    #     print( error_msg )
    #     raise Exception(error_msg)
    
    # map SRR back to the GSM to build file paths 
    srr_gsm_ribo_dict = {}
    srr_gsm_rna_dict  = {}

    # SRR_dict will map Ribo GSM to either Ribo-Seq or matched RNA-Seq, similar to the 
    # final .yaml file
    # It nests the files in the structure of the ['input']['fastq'] files
    SRR_dict = {}

    # Note that the keys of ribo_id and ribo_rna_id are the same.
    # Multiple GSM within One Study
    for cur_gsm in ribo_dict:
        cur_ribo_accessions = []
        cur_rna_accessions  = []
        # do the ribo-seq cases
        cur_ribo_id         = ribo_dict[cur_gsm]['id']
        SRR_dict[cur_gsm]   = {}

        srr_query = 'SELECT "unitmetadata_srr"."id", "unitmetadata_srr"."experiment_id",\
                    "unitmetadata_srr"."sra_accession", "unitmetadata_srr"."creation_date" FROM "unitmetadata_srr"\
                    WHERE "unitmetadata_srr"."experiment_id" = {id}'.\
                        format(id = cur_ribo_id)

        srr_df = pd.read_sql_query(srr_query, conn)

        for index, row in srr_df.iterrows():
            cur_ribo_accessions.append(row['sra_accession'])
            srr_gsm_ribo_dict[row['sra_accession']] = cur_gsm

        SRR_dict[cur_gsm]['riboseq'] = list(cur_ribo_accessions)
        
        if cur_gsm in ribo_rna_dict:
            cur_rna_experiment = ribo_rna_dict[cur_gsm]
            # do the rna-seq cases
            rna_gsm            = cur_rna_experiment['experiment_alias']
            cur_rna_id         = ribo_rna_dict[cur_gsm]['id']
            srr_query          = 'SELECT "unitmetadata_srr"."id", "unitmetadata_srr"."experiment_id",\
                                    "unitmetadata_srr"."sra_accession", "unitmetadata_srr"."creation_date" FROM "unitmetadata_srr"\
                                    WHERE "unitmetadata_srr"."experiment_id" = {id}'.\
                                        format(id = cur_rna_id)
            srr_df             = pd.read_sql_query(srr_query, conn)

            for index, row in srr_df.iterrows():
                cur_rna_accessions.append(row['sra_accession'])
                srr_gsm_rna_dict[row['sra_accession']] = rna_gsm

            SRR_dict[cur_gsm]['rnaseq'] = list(cur_rna_accessions)
            

    base_path            = os.path.join(download_path, gse_only)
    ribo_fastq_path_dict = {}
    rna_fastq_path_dict  = {}
    
    # generate the desired fastq paths
    for val in SRR_dict:
        rna_fastq_path_dict[val] = []

        if 'rnaseq' in SRR_dict[val]:
            for srr_accession in SRR_dict[val]['rnaseq']:
                cur_rna_gsm = srr_gsm_rna_dict[srr_accession]
                rna_fastq_path_dict[val].append(os.path.join(base_path, cur_rna_gsm, srr_accession + '_1.fastq.gz'))

        ribo_fastq_path_dict[val] = []

        for srr_accession in SRR_dict[val]['riboseq']:
            cur_ribo_gsm = srr_gsm_ribo_dict[srr_accession]
            ribo_fastq_path_dict[val].append(os.path.join(base_path, cur_ribo_gsm, srr_accession + '_1.fastq.gz'))
    
    if not len(ribo_fastq_path_dict.values()):
        raise Exception("No Ribo-Seq Experiments Found")
    
    empty_experiments = []

    for experiment in ribo_fastq_path_dict:
        if not ribo_fastq_path_dict[experiment]:
            empty_experiments.append(experiment)
    
    if empty_experiments:
        print("The following experiments have no sequencing files.")

        for exp in empty_experiments:
            print(exp)
        raise Exception("Ribo-Seq GSMs are missing sequencing files.")

    ribo_yaml['input']['fastq'] = ribo_fastq_path_dict

    # delete the entries corresponding to RNA-Seq if they do not exist
    ribo_yaml['rnaseq']['fastq'] = rna_fastq_path_dict
    delete_rnaseq_list           = []

    for experiment in ribo_yaml['rnaseq']['fastq']:
        if len(ribo_yaml['rnaseq']['fastq'][experiment]) == 0:
            delete_rnaseq_list.append(experiment)

    ribo_yaml['rnaseq']['deduplicate'] = dedup_val

    for cur_exp in delete_rnaseq_list:
        del ribo_yaml['rnaseq']['fastq'][cur_exp]
    if len(ribo_yaml['rnaseq']['fastq']) == 0:
        del ribo_yaml['rnaseq']

    ribo_yaml['do_rnaseq']                       = 'rnaseq' in ribo_yaml
    ribo_yaml['output']['output']['base']        = 'output_unitmetadata' + "/" + yaml_name
    ribo_yaml['output']['intermediates']['base'] = "intermediates_unitmetadata"+ "/" + yaml_name
    ribo_yaml['deduplicate']                     = dedup_val

    dir_base   = os.path.join(output, yaml_name.split("_")[0])
    final_name = os.path.join(dir_base, yaml_name + ".yaml")
    os.makedirs(os.path.dirname(final_name), exist_ok=True)

    with open(final_name, mode = 'w') as f:
        yaml.dump(ribo_yaml, f)



def get_parameters():
    parser = argparse.ArgumentParser(
                  formatter_class = argparse.RawTextHelpFormatter,
                  description = \
                      
    """
    This script generates the necessary .yaml input for RiboFlow based on RiboBase databaase entries.
    The user should either provide a single study (for example GSE101760)
    or a text file containing a list of studies where each line contains a separate study.
    Sample Use:\n\n

    python generate_yaml.py --study GSE101760 --template ../project.yaml --output hebele --db ../db/db.sqlite3

    A note on deduplication:
    "deduplicate" is a special keyword for the script ( and for this particular Snakemake workflow ).
    It tells the script to generate the yaml file with deduplication turned on.
    For example, the following command sets the flag "deduplicate: true" in the RiboFlow parameters yaml file.
    As a result, the reads will be deduplicated based on their mapped position and length.

    python generate_yaml.py --study GSE101760_dedup --template ../project.yaml --output hebele --db ../db/db.sqlite3
    """)

    parser.add_argument("--template",
                        type     = str,
                        required = True,
                        help     = "File Path of the Template Yaml File")

    parser.add_argument("--output",
                        type     = str,
                        required = True,
                        help     = "Output Directory of the Generated Yaml")

    parser.add_argument("--study",
                        type     = str,
                        required = False,
                        help     = "Study GSE of Interest")
    
    parser.add_argument("--text",
                        type     = str,
                        required = False,
                        help     = "Text File of containing the list of studies")

    parser.add_argument("--db",
                type     = str,
                required = True,
                help     = "SQLite3 File corresponding to the RiboBase Server")

    parser.add_argument("--download_path",
                type     = str,
                required = False,
                default  = 'input/fastq',
                help     = "File path describing the base fastq download directory")

    parser.add_argument("--reference_folder",
                type     = str,
                required = False,
                default  = 'reference',
                help     = "Folder containing the reference folder")
    
    parser.add_argument("--reference_file",
                type     = str,
                required = False,
                default  = 'scripts/references.yaml',
                help     = "Yaml file containing refrence paths")

    args = parser.parse_args()
    return args


def main():
    params  = get_parameters()

    if bool(params.study) == bool(params.text):
        print("Please specify one of either a study or a text file of studies.")
        exit()
    
    if params.study:
        # Generate yaml file for a single study
        generate_yaml(db               = params.db, 
                      output           = params.output, 
                      template         = params.template, 
                      study            = params.study, 
                      download_path    = params.download_path,
                      reference_file   = params.reference_file, 
                      reference_folder = params.reference_folder)
    else:
        # Generate Ribo File for a list of studies
        # where list is read from the text file given
        with open(params.text) as f:
            for line in f:
                study = line.strip()
                generate_yaml(db               = params.db, 
                              output           = params.output, 
                              template         = params.template, 
                              study            = study, 
                              download_path    = params.download_path,
                              reference_file   = params.reference_file, 
                              reference_folder = params.reference_folder)

if __name__ == "__main__":
    main()
