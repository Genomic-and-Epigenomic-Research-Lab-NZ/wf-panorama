<div align="center">  

# Panorama

</div>

## 👇 Contents

- [Panorama](#panorama)
  - [👇 Contents](#-contents)
  - [🧬 Introduction](#-introduction)
    - [Stage One: In a Pre-Clinical Trial for a defined disease and classification development 📋](#stage-one-in-a-pre-clinical-trial-for-a-defined-disease-and-classification-development-)
    - [Stage Two: Patient sample reporting following a successful clinical trial 💊](#stage-two-patient-sample-reporting-following-a-successful-clinical-trial-)
  - [🖥️ Compute requirements](#️-compute-requirements)
  - [🥳 Third-party requirements](#-third-party-requirements)
  - [🏁 Install and run](#-install-and-run)
    - [Additional install requirements](#additional-install-requirements)
      - [Containers](#containers)
      - [Human Genome Reference](#human-genome-reference)
      - [Chromosome size file](#chromosome-size-file)
    - [Run command](#run-command)
      - [Nextflow Profiles](#nextflow-profiles)
  - [🎈 Usage](#-usage)
    - [General Input Parameters](#general-input-parameters)
      - [Tissue type options](#tissue-type-options)
    - [General Output Parameters](#general-output-parameters)
    - [Mode `make_target_bed`](#mode-make_target_bed)
      - [Input](#input)
      - [Output](#output)
    - [Mode `clin_trial_mode`](#mode-clin_trial_mode)
      - [Input](#input-1)
    - [Mode `clinical_mode`](#mode-clinical_mode)
      - [Input](#input-2)
    - [Outputs for clin\_trial and clinical modes](#outputs-for-clin_trial-and-clinical-modes)
      - [Output directory `<sample name>/wf-humvar-run`](#output-directory-sample-namewf-humvar-run)
      - [Output directory `<sample name>/mod_calling`](#output-directory-sample-namemod_calling)
      - [Output directory `<sample name>/snv_annotation`](#output-directory-sample-namesnv_annotation)
      - [Output directory `<sample name>/sv_annotation`](#output-directory-sample-namesv_annotation)
      - [Output directory `<sample name>/immune_infiltrate`](#output-directory-sample-nameimmune_infiltrate)
  - [🎯 Biomarker panel input](#-biomarker-panel-input)
  - [✍️ Authors](#️-authors)
  - [🗝️ Licence](#️-licence)
  - [🤝 Contributing](#-contributing)
  - [🪢 Related protocols](#-related-protocols)
  - [🗺️ Roadmap](#️-roadmap)
  - [📜 Pipeline History](#-pipeline-history)
  - [🎉 Acknowledgements](#-acknowledgements)
  - [📚 References](#-references)
  - [🤖 AI Assistance](#-ai-assistance)
  - [🔧 Troubleshooting](#-troubleshooting)
    - [Error: method optional() on null object](#error-method-optional-on-null-object)
    - [Error: Building .sif images for wf-human-variation fails](#error-building-sif-images-for-wf-human-variation-fails)
    - [Error: scripts from wf-human-variation "Permission denied"](#error-scripts-from-wf-human-variation-permission-denied)
    - [Error: wf-human-variation has stopped partway through and you want to resume it](#error-wf-human-variation-has-stopped-partway-through-and-you-want-to-resume-it)

## 🧬 Introduction

This is a bioinformatic pipeline utilising the capabilities of nanopore sequencing ([ONT](https://nanoporetech.com/)) to combine multiple biomarkers of different sources; specifically, mutations, methylation and tumour immune infiltrate. This can enable prediction of drug compatibility in cancer tumours. This pipeline was developed as a clinical bioinformatic workflow with the goal of reducing the bioinformatic expertise required to process data. The results of such a tool, with the appropriate pre-clinical trial, can be integrated as part of a protocol aimed to assist molecular pathologists and clinicians to swiftly develop personalised treatment plans.  

Panorama is built using the Nextflow workflow language and is intended to be used in the [EPI2ME](https://github.com/epi2me-labs) framework from [Oxford Nanopore Technologies](https://community.nanoporetech.com).

Panorama is intended for use as follows:  

### Stage One: In a Pre-Clinical Trial for a defined disease and classification development 📋

   1. Biomarkers are input via a csv file, as detailed in the section [Biomarker panel input](#-biomarker-panel-input) below. This will form the basis of the bed file used for adaptive sampling, and is also used for downstream pipeline processes. Biomarkers can include SNVs, methylation markers, and certain immune infiltrate markers as determined by immune deconvolution using methylation markers.
   2. Tumour samples are nanopore sequenced using adaptive sampling and a bed file containing the target regions required for biomarker analysis. This bed file is produced by this pipeline with the flag `--make_target_bed` and the biomarker csv file using the input parameter `--panel_metadata`. Resulting bam files are used as input to Panorama. Bam files must be basecalled using modified base calling (5mC+5hmC contexts), and aligned to the hg38 genome.
   3. Using the flag `--clin_trial_mode` the pipeline will output sample data as a csv file corresponding to the input biomarkers.
   4. Following patient sequence data collection, the researchers will carry out their own classifier development using the biomarker results identified by this pipeline. This is outside this tool's scope and presumably will involve some kind of machine learning. Resulting models can be input to the `--clinical_mode` as is. Alternatively, results can be added to the original biomarker input csv file if in hazard ratio format or similar and used as input for clinical reporting.  
   
Please get in contact with the developers if you want to use this tool in a clinical trial and we will assist/collaborate with trial design, clinical reporting development, and ensuring `--clin_trial_mode` is suitable.

### Stage Two: Patient sample reporting following a successful clinical trial 💊

   1. Using the results of the pre-clinical trial described above, and the final version of the biomarker input with classifier results, a single patient sample can be submitted to this pipeline as before, and will output a report based on the clinical trial results.  
   2. The pipeline can be run for a clinical sample using the flag `--clinical_mode`.

Please note that Stage Two has not been fully developed as there has been no such clinical trial carried out yet. This is expected to be finalised during/after an actual clinical trial. 
<!-- eg read depth requirements -->

## 🖥️ Compute requirements

Recommended requirements:

- CPUs = 32
- Memory = 128GB

Minimum requirements:

- CPUs = 16
- Memory = 32GB

*Based on [wf-human-variation](https://github.com/epi2me-labs/wf-human-variation) as this is the most computationally heavy component of the workflow.*

## 🥳 Third-party requirements

This program uses MethylCIBERSORT and CIBERSORTx for immune infiltrate deconvolution. CIBERSORTx is provided as a Docker container by the developers, the Alizadeh and Newman labs at <https://cibersortx.stanford.edu/>. You will need to create an account and obtain the Docker token through the "Downloads" page. You do not need to install the Docker container, as it has been provided as an Apptainer build called `cibersort_fractions.sif` (see section [Containers](#containers)). You only need to provide your CIBERSORTx username (email) and token to Panorama for it to run. Please note it will take a few days for the CIBERSORTx developers to process your access request.  
Please follow all requirements required by the CIBERSORTx developers as stated when you register. Use of Panorama does not override the CIBERSORTx rules and requirements.  

Also note that internet access is required for the CIBERSORT container to authenticate during the run.  
Please get in contact with us if you are interested in using an alternative immune deconvolution package.  

## 🏁 Install and run

<!-- TODO: Integrate with EPI2ME -->
<!-- You can also access the workflow via the
[EPI2ME Desktop application](https://labs.epi2me.io/downloads/). -->

The workflow uses [Nextflow](https://www.nextflow.io/) to manage compute and software resources, therefore Nextflow will need to be installed before attempting to run the workflow.

> [!WARNING]  
> This workflow has been tested successfully using Nextflow v25.04.6 - v25.10.4. Panorama currently breaks on Nextflow v26.

The workflow can be run using [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html), [Apptainer](https://apptainer.org/) (the open-source fork of Singularity, common on newer HPC systems).

> [!TIP]  
> We recommend using a clean environment with nextflow and apptainer installed.
> ```
> conda create -n pano_env
> conda activate pano_env
> conda install bioconda::nextflow==25.10.4
> conda install conda-forge::apptainer

Once nextflow and apptainer is installed, the following command can be used to obtain the workflow. This will pull the repository in to the assets folder of Nextflow and provide a list of all parameters available for the workflow as well as an example command:

```sh
nextflow run Genomic-and-Epigenomic-Research-Lab-NZ/wf-panorama --help
```

or

```sh
nextflow run http://github.com/Genomic-and-Epigenomic-Research-Lab-NZ/wf-panorama
```

To update a workflow to the latest version on the command line use
the following command:

```sh
nextflow pull Genomic-and-Epigenomic-Research-Lab-NZ/wf-panorama
```

<!-- A demo dataset is provided for testing of the workflow.
It can be downloaded and unpacked using the following commands:
```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-template/wf-template-demo.tar.gz
tar -xzvf wf-template-demo.tar.gz
```
The workflow can then be run with the downloaded demo data using: -->

### Additional install requirements
<!-- TODO: update figshare.com link to actual url once it's confirmed -->

#### Containers

The container components of this tool are very large. These are hosted on [figshare.com](figshare.com):  
<https://doi.org/10.6084/m9.figshare.32165103>  
You will need to download these containers separately in order to use this tool.  
On the plus side, you shouldn't need to download or set up any other packages!

Containers:  

- general.sif (0.5 GB)
- methylcibersort.sif (1.8 GB)
- cibersortx_fractions.sif (0.3 GB)  

Place these in the directory: `wf-panorama/containers/`  

#### Human Genome Reference

You will need to obtain your preferred human genome and associated index file and put it in the `wf-panorama/resources/` directory.  
During development the genome build `GCA_000001405.15_GRCh38_no_alt_analysis_set.fna` was used. Any genome build of hg38 should work, though other builds have not been tested.  

> [!TIP]  
> You should be able to use a symlink (aka symbolic link, alias, shortcut) to avoid having multiple copies of the human genome scattered around your system. To add a symlink:  
>
> ```sh
> ln -s /path/to/<genome.fna> /path/to/wf-panorama/resources/.
> ln -s /path/to/<genome.fna.fai> /path/to/wf-panorama/resources/.
> ```
>
> Please note symlinks haven't been tested in this workflow yet. If there are mysterious errors try copying the genome to this location instead. Please let us know if you use symlinking wth the workflow and it works!

> [!WARNING]  
> Methylation sites have been annotated by their Illumina array probe names, and mapped to hg38 genome locations. Unless you rebuild all the files in this workflow that use genome locations, using the T2T genome build will NOT work for the immune infiltrate section and/or will give wrong results!

To generate the genome index file using samtools (recommended):  

```sh
cd wf-panorama/resources
ref="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"  # or your reference file
samtools faidx ${ref}
# this creates a file called "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai"
```

The index file is required to be pre-built. Ensure that this file is in the `resources` directory.

#### Chromosome size file

A chrom_sizes file is also required. While there is one included in the resources directory, it should be replaced by a fresh file created from the above index file like so:  

```sh
fai="/path/to/<genome.fna.fai>"
cd wf-panorama/resources
cut -f1,2 ${fai} > hg38_no_alt.chrom_sizes
```

This is then passed to the pipeline using the following parameter:  

```sh
--chrom_sizes_file hg38_no_alt.chrom_sizes
```

Ensure that this file is in the `resources` directory. If it is not specified, the pipeline will default to looking for `hg38_no_alt.chrom_sizes` in the `resources` directory.

### Run command

The workflow can be run using:
<!-- TODO: check profile options (epi2me uses "standard") -->
```sh
nextflow run Genomic-and-Epigenomic-Research-Lab-NZ/wf-panorama \
    --clin_trial_mode \
    --project_name Project_Name \
    --sample Test1 \
    --bam_directory /path/to/passed_bams \
    --panel_metadata /path/to/demo_input/panel_metadata.csv \
    --target_bedfile /path/to/demo_input/targets.bed \
    --mCS_cancer_type bladder \
    --cibersortx_username <email> \
    --cibersortx_token <token> \
    -profile slurm,singularity \
    -resume
```

For further information about running a workflow on
the command line see <https://labs.epi2me.io/wfquickstart/>

> [!TIP]  
> If you are new to command line:  
> Install Nextflow (<https://docs.seqera.io/nextflow/#get-started>). Once you have done that, come back here. Note that if you are on a cluster compute you may already have nextflow installed and just need to load it with something like `module load nextflow`.
> "Change directory" (`cd`) to a directory (another name for "folder") where you want to do your analysis.  
>
> ```sh
> cd /home/<your user name>/projects/run_panorama
> ```
>
> This will conduct the analysis inside this directory.
>
> ```sh
> nextflow run Genomic-and-Epigenomic-Research-Lab-NZ/wf-panorama --help
> ```
>
> Note that nextflow will keep things in `/home/<your user name>/.nextflow`

#### Nextflow Profiles

Profiles are selected with the `-profile` flag on the command line. Combine an **executor** profile with a **container** profile using a comma (e.g. `-profile slurm,singularity`).

| Executor profile | Description |
| ---------------- | ----------- |
| `local` | Explicit local execution |
| `hpc_slurm` | SLURM HPC cluster |

| Container profile | Description |
| ----------------- | ----------- |
| `singularity` | Singularity (traditional HPC) |
| `apptainer` | Apptainer (newer HPC systems, open-source Singularity fork) |

> [!TIP]  
> The `local` profile uses apptainer/singularity which should be included in your conda environment.

> [!WARNING]  
> At this stage, the workflow has been tested on a Slurm HPC cluster using Apptainer. While Nextflow should allow the pipeline to be transferable, there may be some issues with other compute setups, and is not currently set up for Docker.

Resource limits can be adjusted with `--max_memory`, `--max_cpus`, and `--max_time` to match your system's available resources.

## 🎈 Usage

There are three modes of operation, selected by the following flags:  

1. `--make_target_bed` - This uses the input biomarker panel to create a bed file with buffer regions, suitable for MinKNOW adaptive sampling. This is different to ONT's "Bed Bugs" tool as it adds the necessary regions for immune infiltrate calculation specific to this tool. You are welcome to double check the output with ONT's tool, accessible [here](https://epi2me.nanoporetech.com/bed-bugs/)
2. `--clin_trial_mode` - This runs sequence data processing, using input bam files, for a single sample. It generates a csv file of biomarker results in relation to the input biomarker panel, and can be used for machine learning classifier training with the goal of predicting response to treatment.
3. `--clinical_mode` - This runs sequence data processing using input bam files and a classifier, presumably generated during the above clinical trial. The input biomarker panel must contain all targets the classifier needs, and these must be captured by adaptive sampling.  

> [!IMPORTANT]
> A biomarker metadata csv file is required for all three modes of operation. See section [Biomarker panel input](#-biomarker-panel-input) for details.

### General Input Parameters

<!-- TODO: maybe add a watch path option once we get up to clinical implementation -->
<!-- TODO: add multiple sample processing -->
| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| make_target_bed | boolean | Generate a target bed file for MinKNOW adaptive sampling. | Uses the biomarker panel metadata to produce buffered target regions. | False |
| clin_trial_mode | boolean | Run sample processing for clinical trial mode. | Processes a single patient sample and outputs data suitable for classifier training. | False |
| clinical_mode | boolean | Generate an individual patient sample report using a trained classifier. | Requires a completed clinical trial and classifier results in the panel metadata. | False |
| base_working_dir | string | Base working directory for the workflow. | This is where all output will be stored. | `~/wf-panorama_run` |
| project_name | string | A project name that will be used for containing all the samples processed during the clinical trial. If using in clinical mode, this will be used for containing all samples processed using the same classifier. | This structure is necessary in order for the biomarker metadata to be processed appropriately. | (Required input for all modes) |
| panel_metadata | string | The path to `<panel_metadata>.csv` | See section [Biomarker panel input](#-biomarker-panel-input) for details | (Required input for all modes) |
| mCS_cancer_type | string | The cancer/tissue type of the project. | Required. See [Tissue type options](#tissue-type-options). | `bladder` |  

#### Tissue type options

The MethylCIBERSORT process has a specific set of genomic locations it uses to generate reference data for deconvolution. Pick the most appropriate one for your project, and input exactly as below into `--mCS_cancer_type`.  

    "acute_myeloid_leukaemia"
    "B_cell_leukemia"
    "B_cell_lymphoma"
    "biliary_tract"
    "bladder"
    "breast"
    "Burkitt_lymphoma"
    "chronic_myeloid_leukaemia"
    "endometrium"
    "Glioma"
    "haematopoietic_neoplasm.other"
    "head_and_neck"
    "Hodgkin_lymphoma"
    "kidney"
    "large_intestine"
    "liver"
    "lung_NSCLC_adenocarcinoma"
    "lung_NSCLC_large_cell"
    "lung_NSCLC_not_specified"
    "lung_NSCLC_squamous_cell_carcinoma"
    "lung_small_cell_carcinoma"
    "lymphoblastic_leukemia"
    "lymphoblastic_T_cell_leukaemia"
    "lymphoid_neoplasm_other"
    "melanoma"
    "mesothelioma"
    "MethylCIBERSORT_KoestlerRuns"
    "Myeloma"
    "neuroblastoma"
    "not_specified"
    "oesophagus"
    "ovary"
    "pancreas"
    "prostate"
    "sarcoma"
    "soft_tissue_other"
    "stomach"
    "T_cell_leukemia"
    "thyroid"

### General Output Parameters

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| project_outdir | string | Project-level output directory (primarily for `make_target_bed`). | Can be overridden manually, but using the default layout is recommended. | `<project_name>` |
| sample_outdir | string | Sample-level output directory for `clin_trial_mode` and `clinical_mode`. | Can be overridden manually, but using the default layout is recommended. | `<project_outdir>/<sample>` |

### Mode `make_target_bed`

This mode has two outputs, for use in MinKNOW adaptive sampling. You must use this mode to generate your adaptive sampling bed file, as it combines your specific targets with regions identified for immune deconvolution by methylation. If your bed file does not contain these regions then the tool will not be able to perform immune deconvolution.  
This mode may need to be rerun multiple times in order to create an optimal bed file that covers all regions adequately while also covering a suitable percentage of the genome. It will check if the resulting bed file meets all the requirements by using the `min_genome_coverage`, `max_genome_coverage` and `buffersize_bp` parameters. If the initial check fails, (it will tell you on the terminal) and/or you want different thresholds for these parameters, adjust them as desired and re-run until you get a successful message.  
For single site targets such as SNPs, buffer regions for adaptive sampling will be added appropriately. For targets that are larger regions, such as an entire gene, this mode adds 2000bp upstream and 1000bp downstream as the adaptive sampling "target region" and THEN adds standard buffer regions on top of these surrounding regions.  
Please note that the MethylCIBERSORT reference data it uses to build the bed file is based on hg38 genome coordinates.  

#### Input

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| chrom_sizes_file | string | The path to a chrom_sizes file, generated from the index of the reference genome. | Optional. If not provided, the bundled `resources/hg38_no_alt.chrom_sizes` file is used. It is strongly recommended to generate this from the same reference genome you are using for sequencing — see [Chromosome size file](#chromosome-size-file). | `resources/hg38_no_alt.chrom_sizes` (bundled) |  

*Example run command*  

```sh
nextflow run Genomic-and-Epigenomic-Research-Lab-NZ/wf-panorama \
    --make_target_bed \
    --base_working_dir ~/wf-panorama_run \
    --project_name Project_Name \
    --panel_metadata /path/to/demo_input/panel_metadata.csv \
    --mCS_cancer_type bladder \
    --chrom_sizes_file hg38_no_alt.chrom_sizes \
    -profile slurm,singularity \
    -resume
```

Expected runtime: < 5 min  

#### Output

| Title | File path | Description |
|-------|-----------|-------------|
| Targets with buffered regions | `minknow_input/<project_name>.targets_buffed.bed` | bed file for adaptive sampling |
| Targets for alignment stats | `minknow_input/<project_name>.targets_for_align.bed` | A bed file provided for optional alignment of your target regions during sequencing. This file does NOT have buffered regions, do not use it in the adaptive sampling input. You may use it in the alignment only section, and it can help monitor read depth in your desired regions. If you are not confident with this do NOT use it. This file is also used for analysis during sample processing. |

### Mode `clin_trial_mode`

#### Input

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample | string | A single sample name or identifier. Must start with an alphabet letter (i.e. not a digit or symbol). | Required for `--clin_trial_mode` and `--clinical_mode`. | null |
| bam_directory | string | The path to a single directory containing bams to process. | Usually the `bam_pass` directory in MinKNOW- or Dorado-processed data. Required for `--clin_trial_mode` and `--clinical_mode`. | null |
| cibersortx_username | string | CIBERSORTx account username (email). | Required to authenticate with the CIBERSORTx service. Required for `--clin_trial_mode` and `--clinical_mode`. | null |
| cibersortx_token | string | CIBERSORTx authentication token. | Obtain from the CIBERSORTx Downloads page. Required for `--clin_trial_mode` and `--clinical_mode`. | null |
| mCS_cancer_type | string | Cancer type signature to use for MethylCIBERSORT deconvolution. | Must match one of the available signature sets bundled with MethylCIBERSORT. See `nextflow_schema.json` for valid options. Required for `--clin_trial_mode` and `--clinical_mode`. | bladder |
| target_bedfile | string | Path to a bed file of target regions to use in the analysis. | Defaults to the bed file generated by `--make_target_bed`. Can be overridden to point to an existing file. | `<project_outdir>/minknow_input/<project_name>.targets_for_align.bed` |
| meth_coverage_threshold | integer | Minimum coverage (depth) threshold for a CpG site to be included in methylation analysis and immune infiltrate deconvolution. | The default of 30 is conservatively high for confident results. For development and testing, a lower value (e.g. 20) can be used to include more sites. | 30 |
<!-- | watch_path | boolean | Enable to continuously watch the input directory for new input files. | This option enables the use of Nextflow's directory watching feature to constantly monitor input directories for new files. | False | -->
<!-- | sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a folder containing sub-folders with FASTQ files. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. A `type` column is required for certain workflows and should have the following values; `test_sample`, `positive_control`, `negative_control`, `no_template_control`. An optional `analysis_group` column is used by some workflows to combine the results of multiple samples. If the `analysis_group` column is present, it needs to contain a value for each sample. |  | -->

*Example run command*  

```sh
nextflow run Genomic-and-Epigenomic-Research-Lab-NZ/wf-panorama \
    --clin_trial_mode \
    --base_working_dir ~/wf-panorama_run \
    --project_name Project_Name \
    --sample Test1 \
    --bam_directory /path/to/bam_pass/ \
    --panel_metadata /path/to/demo_input/panel_metadata.csv \
    --meth_coverage_threshold 5 \
    --cibersortx_username person@place.com \
    --cibersortx_token 12a3bc  \
    -profile slurm,singularity \
    -resume
```

Expected runtime: a couple of hours for `wf-human-variation`, a couple of hours for `modification_calling`. Test data has taken ~4-8 hours total.

### Mode `clinical_mode`

#### Input

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample | string | A single sample name or identifier. Must start with an alphabet letter (i.e. not a digit or symbol). | Required for `--clin_trial_mode` and `--clinical_mode`. | null |
| bam_directory | string | The path to a single directory containing bams to process. | Usually the `bam_pass` directory in MinKNOW- or Dorado-processed data. Required for `--clin_trial_mode` and `--clinical_mode`. | null |
| cibersortx_username | string | CIBERSORTx account username (email). | Required to authenticate with the CIBERSORTx service. Required for `--clin_trial_mode` and `--clinical_mode`. | null |
| cibersortx_token | string | CIBERSORTx authentication token. | Obtain from the CIBERSORTx Downloads page. Required for `--clin_trial_mode` and `--clinical_mode`. | null |
| mCS_cancer_type | string | Cancer type signature to use for MethylCIBERSORT deconvolution. | Must match one of the available signature sets bundled with MethylCIBERSORT. See `nextflow_schema.json` for valid options. Required for `--clin_trial_mode` and `--clinical_mode`. | bladder |
| target_bedfile | string | Path to a bed file of target regions to use in the analysis. | Defaults to the bed file generated by `--make_target_bed`. Can be overridden to point to an existing file. | `<project_outdir>/minknow_input/<project_name>.targets_for_align.bed` |
| meth_coverage_threshold | integer | Minimum coverage (depth) threshold for a CpG site to be included in methylation analysis and immune infiltrate deconvolution. | The default of 30 is conservatively high for confident results. For development and testing, a lower value (e.g. 20) can be used to include more sites. | 30 |
| report_title | string | Title to display on the generated clinical report. | Sets the title heading used in the per-sample clinical report. Defaults to `--project_name` if not provided. | null (uses project_name) |
<!-- | watch_path | boolean | Enable to continuously watch the input directory for new input files. | This option enables the use of Nextflow's directory watching feature to constantly monitor input directories for new files. | False | -->

*Example run command*  

```sh
nextflow run Genomic-and-Epigenomic-Research-Lab-NZ/wf-panorama \
    --clinical_mode \
    --base_working_dir ~/wf-panorama_run \
    --project_name Project_Name \
    --report_title "BCG on NMIBC – Patient Report" \
    --sample Test1 \
    --bam_directory /path/to/bam_pass/ \
    --panel_metadata Genomic-and-Epigenomic-Research-Lab-NZ/demo_input/panel_metadata.csv \
    --meth_coverage_threshold 5 \
    --cibersortx_username person@place.com \
    --cibersortx_token 12a3bc  \
    -profile slurm,singularity \
    -resume
```
Expected runtime: a couple of hours for `wf-human-variation`, a couple of hours for `modification_calling`. Test data has taken ~4-8 hours total.

### Outputs for clin_trial and clinical modes

Raw results for epi2melabs wf-human-variation, SNV, methylation and immune infiltrate sections can be found in their respective directories.  
In `clin_trial_mode` the outputs from these sections (excluding wf-human-variation) are collated into one csv that can be used as input to a machine learning classifier.  
In `clinical_mode` they are processed into a pdf report that uses `resources/template.md` to make.  

#### Output directory `<sample name>/wf-humvar-run`

Along with nextflow-generated directories in `work/` (that can be removed when finished to save space), there are the following files:

- sample.flagstat.tsv
- sample.gene_summary.tsv
- sample.haplotagged.bam
- sample.haplotagged.bam.bai
- sample.mosdepth.global.dist.txt
- sample.mosdepth.summary.txt
- sample.readstats.tsv.gz
- sample.regions.filt.bed.gz
- sample.stats.json
- sample.thresholds.bed.gz
- sample.wf-human-alignment-report.html
- sample.wf-human-snp-report.html
- sample.wf-human-str-report.html
- sample.wf-human-sv-report.html
- sample.wf_mods.1.bedmethyl.gz
- sample.wf_mods.2.bedmethyl.gz
- sample.wf_mods.ungrouped.bedmethyl.gz
- sample.wf_snp.haploblocks.gtf
- sample.wf_snp.vcf.gz
- sample.wf_snp.vcf.gz.tbi
- sample.wf_snp_clinvar.vcf
- sample.wf_str.straglr.tsv
- sample.wf_str.vcf.gz
- sample.wf_str.vcf.gz.tbi
- sample.wf_sv.snf
- sample.wf_sv.vcf.gz
- sample.wf_sv.vcf.gz.tbi

> [!WARNING]  
> Please note that on `-resume`, Nextflow checks if `<project_name>/<sample>/.nextflow_cache/wf-humvar` exists. If it does, the whole process is skipped, regardless of the task hash. This is to prevent unnecessary rerunning of this particular computationally heavy task. If you change `--bed`, `--bam`, etc., the cached wf-humvar directory will still be used. If you need to rerun the wf-human-variation part of the workflow, we recommend you manually delete the entire `<sample_name>` directory to enable a clean restart. If you are familiar with Nextflow and feel confident, you can `cd` into the correct directory and set this workflow to `-resume` independent of Panorama if you wish.

#### Output directory `<sample name>/mod_calling`

- sample.methatlas.csv
- sample.mod_results.csv
- sample.post_beta.csv
- sample.pre_beta.csv
- sample.rawmod_results.csv
- sample.wf_mods.all.bedmethyl.bed
- sample.wf_mods.all.dss_format.tsv

#### Output directory `<sample name>/snv_annotation`

- sample.raw_snv_results.csv
- sample.snv_results.csv

#### Output directory `<sample name>/sv_annotation`
- sample.raw_sv_results.csv
- sample.sv_results.csv

<!-- #### Output directory `<sample name>/methylCS` -->

#### Output directory `<sample name>/immune_infiltrate`
<!-- TODO: Check this output goes to immune infiltrate -->

- CIBERSORTx_sample_Results.csv
- sample.CS_mix_matrix.txt
- sample.bladder.mCS_ref.txt
- sample.immune_panel_results.csv

## 🎯 Biomarker panel input

This is a csv file containing metadata for each biomarker.  
This can be generated using the template excel file provided in [demo_input](demo_input/panel_metadata_template.xlsx), then "Save As" a csv. Ensure the ID column remains 3 digits long (it can be string/character format, not number).  
An example of the csv can be found here: [demo_input/panel_metadata.csv](demo_input/panel_metadata.csv)  

> [!IMPORTANT]
> - ID numbers in the 300's are reserved for molecular characterists comprised of multiple individual biomarkers/features, such as TNBC, HR+, subtypes, etc. These are not directly measured by sequencing but are calculated from multiple biomarker results.
> - ID numbers in the 400's are reserved for immune parameters
> - ID numbers in the 500's are reserved for additional non-molecular factors such as demographic or clinicopathologic indicators that you wish to include in the classifiers but cannot be measured by nanopore sequencing
> - The biomarker types "immune_inf" must not be changed  

> [!TIP]   
> A great gene list to bed file converter is called BED-Craft and can be found here: 
> https://keio-cmg.jp/BED-Craft/index.php
> This can be used to get coordinates to input into the excel metadata file. Please ensure the buffer size is set to 0, as `make_target_bed` will add the appropriate upstream and downstream regions BEFORE the buffer is added.

<!-- TODO: Add references to the clinical report -->

This workflow expects a CSV metadata table describing every biomarker in the panel. A validated example is provided at `demo_input/panel_metadata.csv`.

**Important rules**  
- ID must be a 3-character string (e.g. "001"). Keep it as text (not numeric) so leading zeros are preserved.
- Do not change rows with Biomarker Type `immune_inf` — these are required for immune deconvolution.
- ID ranges:
  - 001–299: individual molecular biomarkers (snv, sv, mod, area_mutations, expression, microsatellite, etc.)
  - 301–399: composite markers (derived from multiple panel entries)
  - 401–499: immune-related outputs / immune_ratio / immune_inf
  - 501–599: demographic / clinicopathology
- Coordinate columns in the CSV may contain commas in the numbers (e.g. "1,234,567"). The pipeline will parse these — but prefer plain integers if possible.
- If "Is this record the whole gene?" is "Yes" for mod/area_mutations entries, the pipeline will extend coordinates by promoter (2000 bp upstream) and downstream (1000 bp) before adding buffering for adaptive sampling. This is why the "strand" is required.

Required CSV columns (as used in the provided CSV)   
| CSV column header | Type | Required? | Notes / valid values |
|-------------------|------|-----------|----------------------|
| ID | string | Yes | 3-digit unique identifier, keep as text. |
| Biomarker name | string | Yes | Human-readable name (e.g. BRCA1, ESR1). |
| Biomarker Type (...) | string | Yes | Must match one of the supported types: snv, sv, mod, area_mutations, expression, exp_ratio, immune_ratio, immune_inf, microsatellite, composite, demographic, clinicopathology. |
| Panel or Area of Interest? | string | Yes | "Panel" or "AOI". "AOI" will not be parsed during clin_trial_mode or clinical_mode. |
| Molecular Group | string | Optional | Text grouping used by downstream reporting/classification. |
| Is this record the whole gene? | string | Conditional | "Yes" or "No". If "Yes", strand must be provided. |
| Is variant in coding region? | string | Conditional | "Yes" / "No" (applies to snv/mod/area_mutations as relevant). |
| SNP ID | string | Optional | e.g. rsIDs for SNVs. |
| Result Options | string | Required for scoring | For genotypic: "Allele1\|Allele2/...". For continuous: "min-max" (e.g. 0.0-1.0). For categorical: options separated by "/". |
| Illumina EPIC ID | string | Optional | EPIC probe ID for mod entries when available. |
| Expression Ratio Components | string | Conditional | For exp_ratio entries, format "GeneA/GeneB" (components referenced by Biomarker name). |
| Combined Marker Components | string | Optional | For composite markers, list component IDs or names. |
| Fusion Components | string | Optional | For fusion (sv) entries, format "GeneA/GeneB" (components referenced by Biomarker name). |
| DNA methylation region | string | Conditional | For mod entries: one of "position", "promoter", "intragenic", "downstream". |
| Scoring Type (categorical, continuous, genotypic) | string | Conditional | Defines how the result will be calculated. |
| Notes | string | Optional | Free text notes for the biomarker. For user reference only. |
| WARNINGS | string | Optional | Any warnings or special handling notes. For user reference only. |
| References | string | Optional | Semicolon-separated references for the biomarker. For user reference only. Future releases may incorporate these into the clinical report. |
| chrom | string | Conditional | Chromosome name, e.g. "chr1". |
| start pos | integer | Conditional | Start coordinate (may contain commas). |
| end pos | integer | Conditional | End coordinate (may contain commas). |
| length | integer | Optional | Can be blank; pipeline will calculate if missing. |
| strand | string | Conditional | "+" or "-". Required if "Is this record the whole gene?" is "Yes". |
| Result | string | Optional | For clinicopathologic / demographic values. Using a copy of the main panel sheet with these values filled in per sample will pass the results through the clin_trial_mode to be output with final csv output file. |

Practical tips  
- Keep the CSV column headers exactly as in the example file to avoid parsing errors.
- For methylation (mod) entries you can supply either a single position (for single CpG) or a broader region (promoter/intragenic/downstream). The mode --make_target_bed will expand the target appropriately for adaptove sampling.
- Use the example CSV to see how composite markers (IDs in the 300s) reference component biomarkers (e.g. HOXB13/IL17BR for an exp_ratio).
- If you plan adaptive sampling, first run --make_target_bed with your panel CSV to generate the MinKNOW target files; those files include extra regions required for immune deconvolution.

Example  
- See `demo_input/panel_metadata.csv` for a full, real-world example that includes snv, sv, mod, expression, exp_ratio, microsatellite and composite entries.
- See `demo_input/panel_metadata_template.xlsx` for a modifiable excel sheet to fill out. "Save As" csv to output to appropriate format, and place in your project directory to be called during the Panorama run.


<!-- ## Troubleshooting -->

<!---Any additional tips.--->
<!-- + If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug. -->
<!-- + See how to interpret some common nextflow exit codes [here](https://labs.epi2me.io/trouble-shooting/). -->

<!-- ## FAQ's -->

<!---Frequently asked questions, pose any known limitations as FAQ's.--->

<!-- If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-template/issues) page or start a discussion on the [community](https://community.nanoporetech.com/). -->

<!-- ## Related blog posts

+ [Importing third-party workflows into EPI2ME Labs](https://labs.epi2me.io/nexflow-for-epi2melabs/)

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts. -->

## ✍️ Authors

- Lucy Picard
  - Postdoctoral Fellow, University of Otago, Wellington, NZ 🇳🇿
  - <https://github.com/lucy924>
  - <lucy.picard@otago.ac.nz>
- Aaron Stevens
  - PI, Genomic and Epigenomic Research Lab (GERL), University of Otago, Wellington, NZ 🇳🇿
  - <aaron.stevens@otago.ac.nz>

## 🗝️ Licence

This project is licensed under the PolyForm Noncommercial License.
Commercial use requires a separate licence from the authors.

## 🤝 Contributing

Contributions are welcome!
Please reach out if you have ideas for this project. We are open to collaboration to improve this tool!  
By contributing, you agree to the [Contributor License Agreement (CLA)](CLA.md), which allows the maintainer to relicense the project in the future.  
Please note the intent of this project is to actively encourage clinical research and development. We want to make a positive impact on clinical care and improve patient outcomes, as freely as possible. If this project gets to a point where it makes sense to commercialise it, we will discuss with all contributors before making significant decisions.  

## 🪢 Related protocols

<!---Hyperlinks to any related protocols that are directly related to this workflow, check the community for any such protocols.--->

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices.  
This protocol currently uses [epi2me-labs/wf-human-variation v2.6.0](https://github.com/epi2me-labs/wf-human-variation/releases/tag/v2.6.0) for initial analysis of bam files.  

## 🗺️ Roadmap

- [ ]  Allow running of multiple samples concurrently
- [ ]  Integrate with Epi2ME Labs
- [ ]  Add option for custom immune infiltrate references
- [ ]  Improve fusion processing
  - [ ]  Output merged bam produced during wf-human-variation to enable manual visualisation
- [ ] Improve "Molecular Group" handling (requires clinical input)
- [ ] Pre-clinical trial enabling practical usage optimisation
<!-- [ ] Add immune infiltrate barplot to outputs and clinical report -->
<!-- [ ] update to more recent version of wf-human-variation -->

<!-- Note:  
Reasons not to use wf-somatic-variation
- not a comparison to patient normal 
- In tumour-only mode, SNV analysis is reliant on ClairS and needs to be tuned for each basecaller model, wf-human-variation is more up-to-date with available models 
- NEW INFO: a matching basecaller model isn't critical for ClairS to run, therefore in the future using wf-somatic-variation in tumour-only momde could be considered -->


## 📜 Pipeline History

Panorama was initially designed as part of the Doctoral thesis entitled:  
"**Epigenetic Consequences of BCG Immunotherapy In Bladder Cancer**"  
For full details on the development of this pipeline the relevant chapter is "Chapter Four: Multi-biomarker discovery using nanopore sequencing technology: proof-of-concept" and can be found in the [Otago archives](https://hdl.handle.net/10523/48489). The original development of the pipeline used [Snakemake](https://snakemake.readthedocs.io/en/stable/) and can be found on [Lucy's github repo](https://github.com/lucy924/nanopore_multiBM_pipeline) with associated material [here](https://github.com/lucy924/Multi-biomarker-ONT-project).

## 🎉 Acknowledgements

Many thanks go to the funders of this project, the Barbara Basham Medical Charitable Trust. Read about the origin of the trust [here](https://wellington.govt.nz/arts-and-culture/heritage/historic-public-memorials/aunt-daisy).

## 📚 References

- CIBERSORT
  - Newman, A.M., Liu, C.L., Green, M.R., Gentles, A.J., Feng, W., Xu, Y., Hoang, C.D., Diehn, M., Alizadeh, A.A., 2015. Robust enumeration of cell subsets from tissue expression profiles. Nat Methods 12, 453–457. <https://doi.org/10.1038/nmeth.3337>
- MethylCIBERSORT
  - Chakravarthy, A., Furness, A., Joshi, K., Ghorani, E., Ford, K., Ward, M.J., King, E.V., Lechner, M., Marafioti, T., Quezada, S.A., Thomas, G.J., Feber, A., Fenton, T.R., 2018. Pan-cancer deconvolution of tumour composition using DNA methylation. Nat Commun 9, 3220. <https://doi.org/10.1038/s41467-018-05570-1>
- BED-Craft
  - Miya, F., Kosaki, K. BED-Craft for nanopore adaptive sampling: a tool for generating bed files with gene names as input data for enrichment sequencing. BMC Res Notes 18, 79 (2025). https://doi.org/10.1186/s13104-025-07152-z
<!-- Citation style: Elsevier (author-date/Harvard, with titles) -->
- Adaptive sampling reference file modified from:
  - Stephane Plaisance (VIB-NC) 2021
  - <https://github.com/Nucleomics-VIB>

## 🤖 AI Assistance

Development of this project used AI-assisted coding tools:

- GitHub Copilot (VS Code extension)
- Claude Sonnet v4.x
- GPT-5 mini
- Prompt layering tool: @Sequera

AI was used for language refactoring from Snakemake to Nextflow, and implementing industry-standard conventions. It was also used for wiring of Nextflow components and debugging.  
AI tools were used interactively. Outputs were evaluated, edited, and validated by the author/s. No code was directly incorporated without review.  

## 🔧 Troubleshooting

### Error: method optional() on null object
**Error:**   
```sh
ERROR ~ Cannot invoke method optional() on null object
```
when wf-human-variation tries to run  

**Problem:** incompatibile nextflow version  

**Fix:** Ensure a compatible nextflow version is being used (examples are: v25.04.6, v25.04.8, v25.10.4)

### Error: Building .sif images for wf-human-variation fails
**Error:**  
On first run of clin_trial_mode or clinical_mode, the workflow will download wf-human-variation and all the associated containers into the `wf-panorama/containers/singularity` directory. You may come across errors such as this:
```
FATAL:   While making image from oci registry: error fetching image to cache: while building SIF from layers: while creating squashfs: /usr/libexec/apptainer/bin/mksquashfs command failed: exit status 137
```
**Problem:** This may be due to an unstable internet connection or similar.   

**Fix:** You can manually pull the images that show as failed in the .nextflow.log *from the work directory stated in the error*:  
```sh
cd /path/to/wf-panorama/containers/singularity
singularity pull --name ontresearch-wf-human-variation-snp-sha8cc7e88.img docker://ontresearch/wf-human-variation-snp:sha8cc7e88ff71bf593d7852309a31d3adb29a7caeb
# ...and any other images from the log
```

### Error: scripts from wf-human-variation "Permission denied"

**Error:**  
```
    .command.sh: line 2: /projects/uow/GERL/dejlu879/Panorama/wf-panorama/modules/wf_humvar/wf-human-variation/bin/resolve_clair3_model.py: Permission denied
```
or   
```
    .command.sh: line 2: /projects/uow/GERL/dejlu879/Panorama/wf-panorama/modules/wf_humvar/wf-human-variation/bin/workflow-glue: Permission denied
```  

**Problem:**  Unsure, potentially something to do with manual downloading the .sif containers above.

**Fix:**  
```
chmod a+rx /path/to/wf-panorama/modules/wf_humvar/wf-human-variation/bin/*
```
Note - you shouldn't need sudo permissions for this as the user who ran the workflow

### Error: wf-human-variation has stopped partway through and you want to resume it

**Problem:** It is difficult to ensure a nested workflow runs with resume working properly, therefore the following steps are recommended if you don't want to begin from the top.  

**Fix:**  
1. Navigate to the work directory that ran the workflow. This should be detailed on the terminal and in the `.nextflow.log` file.
2. Edit the `.command.sh` file to ensure `-resume` is at the end of the nextflow command
3. From that directory, run `bash .command.run`
4. Navigate back to your original working directory and continue as normal, using the `-resume` flag on the Panorama run.

