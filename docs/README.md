# SnakeScale Documentation


--------------------------------

## Running Snakemake

The users are expected to be familiar with the snakemake workflow management system. [Official snakemake documentation](https://snakemake.readthedocs.io/en/stable/) is a good place to start 

### README.md and Snakefile.md

The README.md will specify the required steps to set up a conda environment designed for running the snakemake workflow described below. Snakefile.md will also provide more of the details necessary to get started.

### Usage

The main parameters supported by Snakemake are -j (aliases include --cores and --jobs), -n (aliases include --dry-run and --dryrun), and --config. You can find out more about them here. In our use case on our servers, we typically use ‘-j 48’ to specify 48 cores. For general testing, dry runs can be used by specifying -n in lieu of -j. To get started, the user must specify a list of studies in the studies field in `config/config.yaml`.

A usage template of this workflow is as follows.

```
snakemake -j <number_of_cores>
```

To have a dry run of what jobs this workflow will run, use the following command.

```
snakemake -n
```

When running snakemake, make sure that the file `Snakefile` is your working directory as Snakemake looks for this file to execute the steps. 

In cases where the workflow needs to be split across different computing environments, the --until flag is useful. For example, in the case where the batch download of files needs to occur separate of the computing environment, adding --until download_fastq_files will stop the workflow after downloading the fastq files.

In the case where rules need to rerun, the -R flag is useful. Consider the case where changes were made in the config/config.yaml file or input files edited, adding -R classify_studies will allow for a rerun of all pre-run checks on the current implementation.

### A Note on Pipes in Snakemake Shell Commands

In the case of using pipes in chained shell commands, the workflow will break without the use of "set +o pipefail;". This is a [frequently asked question](https://snakemake.readthedocs.io/en/stable/project_info/faq.html#my-shell-command-fails-with-exit-code-0-from-within-a-pipe-what-s-wrong), and the previously specified command is the provided solution from the developers of Snakemake.

---------------------------------

## Config File

The Snakemake workflow requires an updated config/config.yaml file. The main field to typically change between runs is the studies field, which accepts a list of studies. Note that there is a special ‘_dedup” keyword that enables deduplication to occur during the RiboFlow run. The following fields are present in the config file.

```
Studies:
  # Sample list of studies
  - GSE12345
  - GSE66929
  - GSE66929_dedup
  - GSE56148
  - GSE56148_dedup
  - GSE59815
  - GSE59815_dedup
override: False
riboflow_config: stampede_local
adapter_threshold: 50
threads: # specify number of threads to give each rule
    download_fastq_files: 6
    gzip_fastq: 1
    check_adapter: 4
    check_adapter_stats: 1
    check_lengths: 1
    run_riboflow: 48
check_adapter:
    skipped_reads: 1000000 # default is 1,000,000 reads
    sample_size: 25000 # default is 25,000 reads
guess_adapter:
    skipped_reads: 1000000 # default is 1,000,000 reads
    sample_size: 10000 # default is 10,000 reads 
    min_length: 10 # specifies the minimum length for a valid adapter guess, default is 10 
    max_length: 15 # specifies when to stop extending adapter guess, default is 15 
    skipped_nucleotides: 10 # specifies number of nucleotides to skip within a read, since we are guessing the 3' adapter, default is 10
    seed_length: 6 # specifies the initial seed length to begin extending from, default is 6
    match_ratio: 0.5 # extend the candidate adapter sequence by 1 nt if there is a single base that appears at or more than this fraction, default is 0.5
```


## Pipeline Steps (Snakemake Rules)

### Common Rule - Generating the Input into the Snakemake Workflow

The common rule is the snippet of code that exists independent of the rules, and it is run at the beginning of the workflow. Typically, any import statements and initial starter code for the other rules to use is defined here. In the case of this workflow, the common rule defines all of the helper functions and [input functions](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#functions-as-input-files) used in downstream rules.

The main purpose of the common rule is generating parameters files, for each study, to be used by the RiboFlow pipeline. The list of studies are read from the congfig file (`config/config.yaml`) and their metadata is obtained from the sqlite file (`db/db.sqlite3`). 

#### **Study Information/Setting YAML File (Project.Yaml file Used in RiboFlow)**

This workflow limits the additional overhead required by taking in the same input as RiboFlow. The project.yaml file used in RiboFlow allows the user to specify experiments, file paths to sequencing files, reference files, matched RNA-Seq files, and more.

Before any Snakemake rules are run, an initial block of code (referred to as the common rule) will detect the presence of these project.yaml files. For a study, GSE12345, the workflow will check for the following path (relative to the Snakefile), ‘input/project/GSE12345/GSE12345.yaml’. If it exists, then nothing is changed, but if the file does not exist, a custom script is run to generate the .yaml file as input for the remainder of the workflow.


### Rule download_fastq_files

The first step of the Snakemake workflow is to download all of the files necessary according to the lists of SRR accessions found in the riboflow .yaml files for each study. These are found in the folder input/project/\<GSE\>/ relative to the base path of the Snakefile.

In our experience, we observed that on-the-fly (download of sra combined with fastq conversion) generation of fastq files using fasterq-dump may result in failures occasionally. In order to achieve a more stable download, the prefetch command from [SRA Toolkit](https://www.ncbi.nlm.nih.gov/books/NBK242621/) is used before the [fasterq-dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump) command.


#### **A Note On the “_1.fastq” Suffix:**
It is quite common for RNA-Seq experiments to undergo paired-end sequencing. In order to eliminate any systematic biases between the Ribo-Seq data and matching RNA-Seq data of a given experiment, we only take the first mate pair for paired end sequencing data. When downloading paired-end sequencing files from SRA, fasterq-dump automatically separates each mate pair into separate files (ending in “_1.fastq” and “_2.fastq”). To keep consistency and equal behavior in downstream steps, any files that are single-end are renamed to have a “_1.fastq” suffix. 

### Rule gzip_fastq
After downloading the files, the Snakemake workflow will gzip each of the files for storage efficiency, and these are the final files used as input into RiboFlow. 

### **Preface: Performing Pre-Run Checks on Sequencing Files Before Processing**

During manual curation of all publicly available Ribo-Seq data from published studies, it was increasingly difficult to find the 3’ adapter used during Ribo-Seq library preparation. As a result, many of the early runs of RiboFlow failed due to an incorrect adapter. 

Alongside problems with finding the adapter, there were multiple studies with pre-trimmed fastq files. In these cases, downloaded fastq files had read lengths with widely variable ranges, usually due to adapter trimming done before uploading the data to SRA.

To resolve the difficulties in finding the 3’ adapter, we added additional checks on the downloaded fastq files to check for the low presence of an adapter and guess the adapter, if applicable. In the case of studies with uneven lengths, we added in checks on the read length of the fastq files and removed the adapter (-a) parameter in cutadapt (as part of the assumption that these sequencing files contained adapter-trimmed reads), if applicable.


### Rule check_adapter


During RiboFlow run, cutadapt is used to remove the 3’ adapter and filter these trimmed reads from the raw sequencing file inputs. To prevent wasted compute, we sample a small portion of both our RNA-Seq and Ribo-Seq fastq files and test our cutadapt parameter on them. More specifically, we typically index 1,000,000 reads into the file and take 25,000 continuous reads as input to the cutadapt command of the RiboFlow .yaml file of a given study.

In doing so, we check to make sure that a sufficient percentage of these sampled ribosome profiling reads had the candidate 3’ adapter. In the case of RNA-Seq, we found that roughly one third of the studies with matching RNA-Seq experiments still had a 3’ adapter on the RNA-Seq data, but our pre-run checks resolve other cases of RNA-Seq experiments with no adapter or a different adapter from the one specified in the cutadapt command.

The minimal cutadapt report (specified through the --report=minimal flag) of .tsv format is generated for every experiment of a given study and then fed into the check_adapter_stats rule.

A separate check_adapter rule runs on each SRR accession of a study. Note that a single experiment can be associated with multiple SRR accessions. All of these files are aggregated 


### Rule check_adapter_stats

Once the check_adapter rule runs on all sequencing files of a given study, the generated cutadapt reports are parsed to extract the percentage of the sampled reads with a trimmed adapter.
The output of this rule is a temporary .yaml file used in the modify_yaml rule that is deleted upon completion of the workflow. However, the output of this rule is dumped to a study-specific modifications.yaml file located in the modifications/ directory. The file stores which SRR accessions have a low presence of the adapter. The output of this file is fed into modify_yaml but also fed into guess_adapter. The adapter guessing script will run on any sequencing files below the adapter threshold (specified in the ‘adapter_threshold’ field in config/config.yaml) during the guess_adapter rule.


```
adapter_check:
  has_low_adapter: false
  initial_clip_arguments: -u 1 --maximum-length=40 --minimum-length=15 --quality-cutoff=28
    -a TGGAATTCTCGGGTGCCAAGG --overlap=4 --trimmed-only
  low_adapter_files: []
  rnaseq:
    clip_arguments: -u 5 -l 40 --quality-cutoff=28 --overlap=4 --trimmed-only
    has_rnaseq_adapter: false
    initial_clip_arguments: -u 5 -l 40 --quality-cutoff=28 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
      --overlap=4 --trimmed-only
    low_adapter_files:
    - SRR578266
    - SRR578268
```

#### **Behavior Differences between Ribo-Seq and RNA-Seq Experiments:**
The above example is a snippet of the modifications.yaml corresponding to the output of the check_adapter_stats rule. Note the presence of the ‘initial_clip_arguments field’ in the default fields and the nested ‘rnaseq’ field. A unique element in the ‘rnaseq’ field is the presence of the ‘clip_arguments” field, which potentially differs from the initial arguments due to the [-a parameter](https://cutadapt.readthedocs.io/en/stable/guide.html#regular-3-adapters). If all RNA-Seq experiments have a low adapter presence, then it is empirically likely that there is no adapter at all (or that there is an adapter but the initially provided one is incorrect). This is a crucial step during the modify_yaml rule, where all of the pre-run steps are pooled to generate a modified .yaml file for RiboFlow.


### Rule guess_adapter

The guess_adapter rule iterates over the sequencing files that were flagged with a low adapter presence and makes a guess on the adapter. In the case of there being multiple and different adapter guesses, the study is rendered invalid. Similar to check_adapter_stats, the rule creates a temporary output file that is dumped in the modifications.yaml file output from the modify_yaml rule.

```
adapter_guess:
  GSM1634443:
    SRR1916542: TGGAATTCTCGGGTGCCAAGG
  GSM1634444:
    SRR1916543: TGGAATTCTCGGGTGCCAAGG
  GSM1634445:
    SRR1916544: TGGAATTCTCGGGTGCCAAGG
  consensus_adapter: TGGAATTCTCGGGTGCCAAGG
  detected_adapters:
  - TGGAATTCTCGGGTGCCAAGG
  has_multiple_adapters: false
```

The example above represents a case in which the guessed adapter can be updated as the valid adapter. Since all of the adapter guesses are the same, there is a consensus adapter that will be used in lieu of the initial candidate adapter specified in the study’s project.yaml file.

### Rule check_lengths

We found several cases in which the downloaded fastq files had uneven lengths, hinting at the lack of a need to trim the adapter and provide the `-a` parameter for cutadapt. The check_lengths rule stores the minimum and maximum length of the first 1,000 reads for each sequencing file in a study.

The output of check_lengths follows similar to the other pre-run rules in generating a temporary file that is parsed during the modify_yaml rule. The following is an example of what is stored.

```
length_check:
  clip_arguments: -u 1 --maximum-length=40 --minimum-length=15 --quality-cutoff=28
    -a CTGTAGGCACCATCAAT --overlap=4 --trimmed-only
  has_all_uneven_lengths: false
  initial_clip_arguments: -u 1 --maximum-length=40 --minimum-length=15 --quality-cutoff=28
  lengths:
    SRR1158798:
      max_length: 51
      min_length: 51
    SRR1158799:
      max_length: 51
      min_length: 51
    SRR1514697:
      max_length: 51
      min_length: 18
    SRR837789:
      max_length: 36
      min_length: 36
  uneven_files:
    - SRR1514697
```

The minimum and maximum read lengths for each sequencing file are stored, and the list of uneven files is stored as well. Note the presence of the ‘has_all_uneven_lengths’ field which, if true, allows the workflow to make the assumption that there is no adapter. That is, if all of the sequencing files of a given type (Ribo-Seq or RNA-Seq) have uneven lengths, then the cutadapt parameter can be considered for removal. 

### Rule modify_yaml

The modify_yaml aggregates the outputs of check_adapter_stats, check_lengths, and guess_adapter and resolves any inconsistencies between these pre-run checks and the initial provided parameters. This rule outputs a modified project.yaml file that is used for the RiboFlow run. 


### Commonly Encountered Cases

#### **A Study with No Sequencing Files of Uneven Length, No Sequencing Files with a Low Presence of the Adapter**

This is the ideal straightforward case in which all of the adapters for all of the sequencing files, both RNA-Seq and Ribo-Seq. In this case, the modified yaml file that is the output of this rule will have no difference from the original file.

#### **A Study with No Sequencing Files of Uneven Length and Incorrect Adapters in the RNA-Seq Experiments**

In this commonly encountered case, the guess_adapter script would be run on all of the RNA-Seq experiments. 

The most common result is that the adapter guessing script is unable to find a commonly occurring adapter for all of the RNA-Seq fastq files (because there is not one). In this case, the -a parameter is removed in the modified yaml file.

However, there have been cases where there was an incorrect adapter for all of the RNA-Seq files, and there is a 3’ adapter used during the library prep of these experiments. Given sufficiently long read lengths (typically >50 nt), the output of the adapter guessing script is further examined.

If there are multiple detected adapters among a given sequencing type(note that failing to find an adapter counts as a guess), then the study is considered to be invalid. RiboFlow will not be run (unless the override parameter is specified in the config file).

If there is a single consensus adapter for a sequencing type, then the guessed adapter will replace the original adapter present in the study’s settings YAML file.

#### **A Study with No Sequencing Files of Uneven Length and Incorrect Adapters in the Ribo-Seq Experiments**

In this case, the guess adapter script would be run on all of the Ribo-Seq fastq files. If there are multiple detected adapters, then the study is considered to be invalid. RiboFlow will not be run. If there is a single consensus adapter for a sequencing type, then the guessed adapter will replace the original adapter present in the study’s settings YAML file.

### **Rule classify_studies**

This rule takes in all of the outputs from the pre-run checks (check_adapter_stats, check_lengths, guess_adapter) of all of the studies. This is the last step before the RiboFlow runs. During this rule, the studies are filtered as successfully passing or failing the pre-run checks. The log files are generated depending on the different pre-run checks (listing fastq files with low adapter percentages, uneven lengths, etc), and placed in either log/success/ or log/failed/ directory. From there, only studies placed in the success directory will run RiboFlow (unless the ‘override’ parameter in the config file is set to True).


### **Rule run_riboflow**

This rule is essentially a wrapper around RiboFlow. For any valid studies, RiboFlow will be run, and the outputs will be stored in output/\<GSE\>/ or output/\<GSE\>_\<run_label\>. For example, initially we generated runs of the form output/GSE12345 and output/GSE12345_dedup where GSE12345 is the run with the default parameters withOUT deduplication and GSE12345_dedup is the run with deduplication. 

---------------------------------

## Technical Details of the Snakemake Workflow

### A Note on Snakemake Wildcards

To allow for a rule to be more generalizable, Snakemake features wildcards, which are essentially variables/fields specified in the workflow that are automatically resolved during runtime. The workflow described here uses these wildcards to be more generalizable across multiple studies 

However, providing sufficient information within the rules is critical to allow Snakemake to resolve ambiguities. In this workflow, we use wildcards typically on rules that only apply to a single fastq file followed by a rule that aggregates each of the fastq files in a study together. An example of this is the check_adapter rule, which is run on a single SRR accession (a single fastq file), coming before the check_adapter_stats rule, which aggregates all of the outputs (cutadapt summaries) of check_adapter corresponding to a single study together. In our experience, interweaving this aggregating steps with individual steps makes the code more maintainable, and having these steps also trains a better intuition for the programmer in developing the workflow. 

### A Note on the classify_studies Rule

The outputs of this rule will be reset every time that the Snakemake workflow is run. That is, the log/success and log/failed will only contain information on the pre-run checks of the studies specified in the most recent Snakemake workflow invocation.

This is done so that the directories do not get overfilled with studies that are not immediately relevant or being currently run. Instead, the user gets a quick glimpse of the number of studies that passed the pre-run checks and the number of studies that failed. When curious about the studies that failed, the modifications.log file that is the output of this rule will provide additional explanation on why the study failed.
In subsequent reruns where new studies of interest need to undergo pre-run checks, reconfiguring the config.yaml to include any studies of interest and running the pre-run checks again is convenient using the following command. The -R flag allows the user to specify what rules would like to be rerun.

```
snakemake -j <cores> -R classify_studies
```

### A Note on Studies with Multiple Organsims

Currently, this system doesn’t support studies with multiple organisms.

If a study contains multiple organisms ( one experiment from Human samples and another experiment from mouse ), then the study is discarded.

```
     if len(organism_set) != 1:
        print("More than one organism detected.")
        raise Exception("More than one organism detected")
```

For now, we can assemble the project.yaml files of these studies manually and run them. For example, for the above study, there would be the following runs

  - GSE12345_homosapiens-dedup
  - GSE12345_homosapiens
  - GSE12345_musmusculus-dedup
  - GSE12345_musmusculus

---------------------------------

## FAQ

**Is there a way to skip previously processed studies? Since we don’t keep track of riboflow runs in the database, this might be a bit non-trivial.** 

Yes, we keep track of our processed studies in our Google doc. We can also use the contents of the outputs to determine what has been run, and what has not. *Soon, we may integrate tracking into RiboBase-web altogether.*

**In the workflow schematic, we see branching at several rules. For example: “gzip fastq” can go to check lengths or check adapters. How is this determined? Or do we have both paths being taken regardless?**

The branching is conditional. So it would be helpful to update the figure to show conditional branching. The rules are all entered and run, but what happens in the rules depends on the current state of the data. For example, during guess_adapter, every study will run this rule, but we will not try to guess the adapter unless the earlier check_adapter_stats rule states that there is low presence of the adapter. 

**How special is the “_dedup” keyword? Can’t we generalize it to any type of run? For example, a year from now, we might want to process everything using a different reference and keep the previous reference. To distinguish the cases, we could simply do “GSE12345_newhumanref”. Can we do that? How hard (or easy) is this?**

In the pipeline “dedup” is a special keyword and it is hardcoded. We need to update the rule for another keyword. Answer: At this point, dedup is a hardcoded special keyword. If we want to add another keyword, we need to add clauses to our implementation, 