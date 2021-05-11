# Snakefile Information to Automate RiboFlow Runs of Many Studies
The Snakefile in the repository defines a workflow that automatically parses a project.yaml file (necessary for a RiboFlow run). This workflow will download a study's associated fastq files, gzip them, perform pre-run adapter checking, and run RiboFlow.

## Parameters
The main parameters supported by Snakemake are  `-j` (aliases include `--cores` and `--jobs`), `-n` (aliases include `--dry-run` and `--dryrun`), and `--config`. You can find out more about them [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html).

In our use case on [Lonestar5](https://portal.tacc.utexas.edu/user-guides/lonestar5) of TACC, we typically use the `-j 48` to specify 48 cores. For general testing, dry runs can be used by specifying `-n` in lieu of `-j <int>`. To get started, the user must specify a list of studies in the `studies` field in `config/config.yaml`. 

A usage template of this workflow is as follows.

    snakemake -j <number_of_cores> 

To have a dry run of what jobs this workflow will run, use the following command.
    
    snakemake -n

Lastly, in cases where the workflow needs to be split across different computing environments, the `--until` flag is useful. For example, in the case where the batch download of files needs to occur separate of the computing environment, adding `--until download_fastq_files` will stop the workflow after downloading the fastq files. Additionally, in the case where rules need to rerun, the `-R` flag is useful. Consider the case where changes were made in the `config/config.yaml` file or input files edited, adding `-R classify_studies` will allow for a rerun of all pre-run checks on the current implementation.

## Input Requirements

### RiboFlow Yaml File
Aside from the Snakefile, the settings file, which is the `project.yaml` file required as input into (RiboFlow)[https://github.com/ribosomeprofiling/riboflow], needs to be specified. Also note that this workflow will use the path `input/project/<GSE>/<GSE>.yaml` to find the structure. This means that the `project.yaml` file corresponding to a given study with some accession number `<GSE>` will be named as `<GSE>.yaml`. 

From this yaml file, the workflow will download the necessary sequencing files at the base path of `input/fastq/<GSE>/`. However, note that within the `<GSE>` directory, the different experiments (usually in the form of GSM accessions) will be listed, based on what is specified in the `<GSE>.yaml` file. Within each `<GSM>` directory, there will be the fastq files of the SRR accessions associated with the study. See the `Directory Structure` below for a schematic. 

### Reference Files
A list of available reference files can be found [here](https://github.com/RiboBase/riboflow_references). It is important to make sure that the paths specified in the `<GSE>.yaml` file lead to the correct files.


## Directory Structure

The following is an overview of paths assumed through the various steps of the workflow.

```
Workflow
│   README.md
│   Snakefile
│   Snakefile.md
│
└───config
│   │   config.yaml
│
└───envs
│   │   environment.yaml
│
└───scripts
│   │   ...
│
└───riboflow
│   │   RiboFlow.groovy
│   |   ...
│ 
└───db
│   │   db.sqlite3
│   
└───reference
│   └───filter
│   │   └───human
│   │   │   │   ...
│   │   └───mouse
│   │       │   ...
│   └───transcriptome
│       └───human
│       │   │   ...
│       └───mouse
│           │   ...
│
└───input
│   │
│   └───fastq
│   │   └───study_GSE_1
│   │   │   └───Ribo_Seq_GSM_1
│   │   │   │   │   Ribo_Seq_SRR_1.fastq
│   │   │   │   │   ...
│   │   │   │
│   │   │   └───Ribo_Seq_GSM_2
│   │   │   │    │   ...
│   │   │   │
│   │   │   └───RNA_Seq_GSM_1
│   │   │   │   │   Matched_RNA_Seq_SRR_1.fastq
│   │   │   │   │   ...
│   │   │   │
│   │   │   └───Ribo_Seq_GSM_2
│   │   │       │   Matched_RNA_Seq_SRR_2.fastq
│   │   │       |   ...
│   │   │  
│   │   └───study_GSE_2
│   │   │   │
│   │   │   └───Ribo_Seq_GSM_1
│   │   │   │   │   Ribo_Seq_SRR_1.fastq
│   │   │   │   │   ...
│   │   │   │
│   │   │   └───Ribo_Seq_GSM_2
│   │   │   │    │   ...
│   │   │   │
│   │   │   └───RNA_Seq_GSM_1
│   │   │   │   │   Matched_RNA_Seq_SRR_1.fastq
│   │   │   │   │   ...
│   │   │   │
│   │   │   └───Ribo_Seq_GSM_2
│   │   │       │   Matched_RNA_Seq_SRR_2.fastq
│   │   │       |   ...
│   │   │   
│   │   └─── ...
│   │
│   │
│   └───project
│       └───study_GSE_1
│       │   │   study_GSE_1.yaml
│       │   │   study_GSE_1_dedup.yaml
│       │   │   ...
│       │
│       └───study_GSE_2
│       │   │   study_GSE_2.yaml
│       │   │   study_GSE_2_dedup.yaml
│       │   │   ...
│       │
│       └───study_GSE_3
│       │   │   study_GSE_3.yaml
│       │   │   study_GSE_3_dedup.yaml
│       │   │   ...
│       │
│       └───study_GSE_4
│       │   │   study_GSE_4.yaml
│       │   │   study_GSE_4_dedup.yaml
│       │   │   ...
│       │  
│       └─── ...
│   │
│   │
│   └───modified_project
│       └───study_GSE_1
│       │   │   study_GSE_1_modified.yaml
│       │   │   study_GSE_1_dedup_modified.yaml
│       │   │   ...
│       │
│       └───study_GSE_2
│       │   │   study_GSE_2_modified.yaml
│       │   │   study_GSE_2_dedup_modified.yaml
│       │   │   ...
│       │
│       └───study_GSE_3
│       │   │   study_GSE_3_modified.yaml
│       │   │   study_GSE_3_dedup_modified.yaml
│       │   │   ...
│       │
│       └───study_GSE_4
│       │   │   study_GSE_4_modified.yaml
│       │   │   study_GSE_4_dedup_modified.yaml
│       │   │   ...
│       │  
│       └─── ...
│   │
│   │
│   └───modifications
│       └───study_GSE_1
│       │   │   modification.yaml
│       │
│       └───study_GSE_2
│       │   │   modification.yaml
│       │  
│       └─── ...
│   
└───intermediates
│   └───study_GSE_1
│   │   │   ...
│   └───study_GSE_2
│   │   │   ...   
│   └─── ...
│
└───output
    └───study_GSE_1
    │   |   ...
    └───study_GSE_2
    │   │   ...
    └─── ...
```


