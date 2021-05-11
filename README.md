# SnakeScale 
A Snakemake pipeline to automate RiboFlow runs at scale.

## Introduction
This repository provides many of necessary tools to automate RiboFlow runs coming from our curated set of GEO studies. It depends on RiboBase for the metadata of many published studies that generated ribosome profiling data. To find out more about the Snakemake workflow, please read the `Snakefile.md` file.

---------------------------------

## Setup


### Clone this Repository

Start by cloning the repository and going into the repository folder. 

```
git clone git@github.com:RiboBase/snakescale.git
cd snakescale
```

### Obtain the Database (Sqlite) File

Metadata is provided in an sqlite file, included in [**another** repository: db-sqlite](https://github.com/RiboBase/db-sqlite).
Clone [db-sqlite](https://github.com/RiboBase/db-sqlite) and move the sqlite file to the `db` folder.

```
git clone git@github.com:RiboBase/db-sqlite.git
mv db-sqlite/db.sqlite snakescale/db/
rm -rf db-sqlite
```

### Conda Environment
Running the Snakefile in this repository requires the `snakemake-ribo` conda environment. The `environment.yaml` file in this repository lists the dependencies used.

Set up the conda environment by typing in the following command.

```
conda env create -f envs/environment.yaml
```

Next, activate the environment.

```
conda activate snakemake-ribo
```



### Modify The Config File

This pipeline essentially requires a list of studies where each study is identified by its GEO accession (for example *GSE102375*). This list is provided in the file `config/config.yaml` under the section *studies*. 

In the list of studies, the keyword **dedup** has a special meaning. If it is appended to the GEO accession ID, seperated by an underscore (_), then the study will be processed in RiboFlow with the *dedup* flag on. By default, in other words, if *dedup* is not provided, then RiboFlow works with *dedup* parameter set to False. 

### Download the References
There is a curated set of RiboFlow references in [RiboBase Github Organization](https://github.com/RiboBase/riboflow_references/blob/main/README.md). The following command will download and organize these references under the folder `reference`

```
python scripts/download_reference.py --target reference --yaml scripts/references.yaml
```

After a succesful download, the following command lists the bowtie2 references (with the *bt2* extension), for the human transcriptome. 

```
ls reference/transcriptome/human/
```

If no files are seen, something must have gone wrong in downloading the references.


### A note on *libtbb*

Make sure that bowtie2 is working properly. The following command lists the help of bowtie2. 

```
bowtie2 -h
```

If you see an error message, you may need to install libtbb package. If you see the help page, then you can skip the next paragraph and proceed to the next section. 

We observed that, on some systems (for example Ubuntu 20.04), *libtbb* is required to run bowtie2. The error is observed first at the filtering step of RiboFlow. To solve this problem, *libtbb* library needs to be installed ( [see also](http://www.metagenomics.wiki/tools/bowtie2/install/libtbb-so-2) ). The following installs the required library on Ubuntu 20.04:

```
sudo apt-get install libtbb-dev
```


### Run Snakemake

Before running the pipeline, it is a good idea to check some of the required files.

1) The name of the working folder must be `SnakeScale`:
   ```
   pwd
   ```

2) Make sure the sqlite file is present:
   ```
   ls db/db.sqlite3
   ```

3) Make sure reference files are in place. For example, the following should list the bowtie2 refrences for the human transcriptome. 
   ```
   ls reference/transcriptome/human/*bt2
   ```

Finally, you can run Snakemake.

```
snakemake -p --cores 48
```

Adjust the parameter *--cores* according to the number of available cores in your system.


---------------------------------

## Supported Organisms

A list of supported organisms can be found [here](https://github.com/RiboBase/riboflow_references/blob/main/README.md).

New organisms can be introduced by adding them to the file `scripts/references.yaml` in this repository.

---------------------------------

## Documentation

More documentation is available [here](https://github.com/ribobase/snakescale/tree/main/docs).

---------------------------------

## Link to Google Sheet
We keep a track of the status of studies at googledocs.

[Here is the link to Google Sheets for run logs.](https://docs.google.com/spreadsheets/d/1kVdRq7d5-IM2tz3EG0POHNPtBjI_wxeLjw7kZr2N2mg/edit#gid=346131063)
