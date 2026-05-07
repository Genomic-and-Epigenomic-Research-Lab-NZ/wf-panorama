<div align="center">  

# Panorama

</div>

## 👇 Contents

- [Panorama](#panorama)
  - [👇 Contents](#-contents)
  - [🧬 Introduction](#-introduction)
    - [Stage One: In a Pre-Clinical Trial for a defined disease and single\* treatment option 📋](#stage-one-in-a-pre-clinical-trial-for-a-defined-disease-and-single-treatment-option-)
    - [Stage Two: Patient sample reporting following a successful clinical trial 💊](#stage-two-patient-sample-reporting-following-a-successful-clinical-trial-)
  - [🖥️ Compute requirements](#️-compute-requirements)
  - [🥳 Third-party requirements](#-third-party-requirements)
  - [🏁 Install and run](#-install-and-run)
    - [Additional install requirements](#additional-install-requirements)
    - [Containers](#containers)
    - [Human Genome Reference](#human-genome-reference)
    - [Chromosome size file](#chromosome-size-file)
    - [Run command](#run-command)
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
      - [Output directory `<sample name>/methylCS`](#output-directory-sample-namemethylcs)
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

## 🧬 Introduction

This is a bioinformatic pipeline utilising the capabilities of nanopore sequencing ([ONT](https://nanoporetech.com/)) to combine multiple biomarkers of different sources; specifically, mutations, methylation and tumour immune infiltrate. This can enable prediction of drug compatibility in cancer tumours. This pipeline was developed as a clinical bioinformatic workflow that requires little bioinformatic expertise to use. The results of such a tool, with the appropriate pre-clinical trial, can be integrated as part of a protocol aimed to assist molecular pathologists and clinicians to swiftly develop personalised treatment plans.  

Panorama is built using the Nextflow workflow language and is intended to be used in the [Epi2ME](https://github.com/epi2me-labs) framework from [Oxford Nanopore Technologies](https://community.nanoporetech.com).

Panorama is intended for use as follows:  

### Stage One: In a Pre-Clinical Trial for a defined disease and single* treatment option 📋
   1. Biomarkers are input via a csv file, as detailed in the section [Biomarker panel input](#biomarker-panel-input) below. This will form the basis of the bed file used for adaptive sampling, and is also used for downstream pipeline processes. Biomarkers can include SNVs, methylation markers, and certain immune infiltrate markers as determined by immune deconvolution using methylation markers.
   2. Tumour samples are nanopore sequenced using adaptive sampling and a bed file containing the target regions required for biomarker analysis. This bed file is produced by this pipeline with the flag `--make_target_bed` and the biomarker csv file using the input parameter `--panel_metadata`. Resulting bam files are used as input to Panorama. Bam files must be basecalled using modified base calling (5mC+5hmC contexts), and aligned to the hg38 genome.
   3. Using the flag `--clin_trial_mode` the pipeline will output sample data as a csv file corresponding to the input biomarkers.
   4. Following patient sequence data collection, the researchers will carry out their own classifier development using the biomarker results identified by this pipeline. This is outside this tool's scope and presumably will involve some kind of machine learning. Results are added to the original biomarker input csv file and used as input for clinical reporting.

**Multiple treatment options will be available in a future update*

### Stage Two: Patient sample reporting following a successful clinical trial 💊
   1. Using the results of the pre-clinical trial described above, and the final version of the biomarker input with classifier results, a single patient sample can be submitted to this pipeline as before, and will output a report based on the clinical trial results.  
   2. The pipeline can be run for a clinical sample using the flag `--clinical_mode`.
   
Please note that Stage Two has not been fully tested as there has been no such clinical trial carried out yet. This is expected to be finalised during/after an actual clinical trial. Please get in contact with the developers if you decide to use this tool in a clinical trial and we will assist/collaborate with trial design and clinical reporting development.
<!-- eg read depth requirements -->

## 🖥️ Compute requirements

Recommended requirements:

+ CPUs = 32
+ Memory = 128GB

Minimum requirements:

+ CPUs = 16
+ Memory = 32GB

*Based off [wf-human-variation](https://github.com/epi2me-labs/wf-human-variation) as this is the most computationally heavy component of the workflow.*

## 🥳 Third-party requirements

This program uses MethylCIBERSORT and CIBERSORTx for immune infiltrate deconvolution. CIBERSORTx is provided as a Docker container by the developers, the Alizadeh and Newman labs at https://cibersortx.stanford.edu/. You will need to create an account and obtain the Docker token through the "Downloads" page. You do not need to install the container yourself, Panorama will handle that for you. You only need to provide your CIBERSORTx username (email) and token to Panorama for it to run. Please note it will take a few days for the CIBERSORTx developers to process your access request.  
Please follow all requirements required by the CIBERSORTx developers as stated when you register. Use of Panorama does not override the CIBERSORTx rules and requirements.  

Please get in contact with us if you are interested in using an alternative immune deconvolution package.  

## 🏁 Install and run

<!-- TODO: check if this is doable -->
<!-- You can also access the workflow via the
[EPI2ME Desktop application](https://labs.epi2me.io/downloads/). -->

The workflow uses [Nextflow](https://www.nextflow.io/) to manage
compute and software resources,
therefore Nextflow will need to be
installed before attempting to run the workflow.

The workflow can currently be run using [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html).
<!-- This is automated provided Singularity is installed.
This is controlled by the
[`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles)
parameter as exemplified below. -->

It is not required to clone or download the git repository in order to run the workflow.
<!-- More information on running EPI2ME workflows can be found on the [EPI2ME website](https://labs.epi2me.io/wfindex). -->

Once nextflow is installed, the following command can be used to obtain the workflow. This will pull the repository in to the assets folder of Nextflow and provide a list of all parameters available for the workflow as well as an example command:

```
nextflow run lucy924/wf-panorama --help
```
To update a workflow to the latest version on the command line use
the following command:
```
nextflow pull lucy924/wf-panorama
```

<!-- A demo dataset is provided for testing of the workflow.
It can be downloaded and unpacked using the following commands:
```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-template/wf-template-demo.tar.gz
tar -xzvf wf-template-demo.tar.gz
```
The workflow can then be run with the downloaded demo data using: -->

### Additional install requirements
<!-- TODO: Set up Github Actions to get this to work properly -->
<!-- TODO: update figshare.com link to actual url -->

<!-- BEGIN:external_downloads -->
### Containers
The container components of this tool are very large. These are hosted on [figshare.com](figshare.com):  
https://figshare.com/account/articles/32165103  
You will need to download these containers separately in order to use this tool.  
On the plus side, you shouldn't need to download or set up any other packages!

Containers:  
- general.sif (0.4 GB)
- methylcibersort.sif (1.8 GB)
Place these in the directory: `wf-panorama/containers/`  

### Human Genome Reference
You will need to obtain your preferred human genome and associated index file and put it in the `wf-panorama/resources/` directory.  
During development the genome build `GCA_000001405.15_GRCh38_no_alt_analysis_set.fna` was used. Any genome build of hg38 should work, though other builds have not been tested.  

> [!TIP]
> You should be able to use a symlink (aka symbolic link, alias, shortcut) to avoid having multiple copies of the human genome scattered around your system. To add a symlink:  
> ```
> ln -s /path/to/<genome.fna> /path/to/wf-panorama/resources/.
> ln -s /path/to/<genome.fna.fai> /path/to/wf-panorama/resources/.
> ```
> Please note symlinks haven't been tested in this workflow yet. If there are mysterious errors try copying the genome to this location instead. Please let us know if you use symlinking wth the workflow and it works!

> [!WARNING] 
> Methylation sites have been identified by their Illumina array probe names, and mapped to hg38 genome locations. Unless you rebuild all the files in this workflow that use genome locations, using the T2T genome build will NOT work and/or will give wrong results!

To generate the genome index file using samtools (recommended):  
```sh
cd wf-panorama/resources
ref="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"  # or your reference file
samtools faidx ${ref}
# this creates a file called "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai"
```
The index file is required to be pre-built. Ensure that this file is in the `resources` directory.   

### Chromosome size file
A chrom_sizes file is also required. While there is one included in the resources directory, it should be replaced by a fresh file created from the above index file like so:  
```sh
fai="/path/to/<genome.fna.fai>"
cd wf-panorama/resources
cut -f1,2 ${fai} > hg38_no_alt.chrom_sizes
```
This is then passed to the pipeline using the following parameter:  
```
--chrom_sizes_file hg38_no_alt.chrom_sizes
```
Ensure that this file is in the `resources` directory. If it is not specified, the pipeline will default to looking for `hg38_no_alt.chrom_sizes` in the `resources` directory.
<!-- END:external_downloads -->

### Run command

The workflow can be run using:
<!-- TODO: fix profile options (epi2me uses "standard") -->
```
nextflow run lucy924/wf-panorama \
    --project_name Project-BCG_on_NMIBC \
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
the command line see https://labs.epi2me.io/wfquickstart/

> [!Tip] 
> If you are new to command line:  
> Install Nextflow (https://docs.seqera.io/nextflow/#get-started). Once you have done that, come back here. Note that if you are on a cluster compute you may already have nextflow installed and just need to load it with something like `module load nextflow`.
> "Change directory" (`cd`) to a directory (another name for "folder") where you want to do your analysis.  
> ```
> cd /home/<your user name>/projects/run_panorama
> ```
> This will conduct the analysis inside this directory. 
> ```
> nextflow run lucy924/wf-panorama --help
> ```
> Note that nextflow will keep things in `/home/<your user name>/.nextflow`


## 🎈 Usage
There are three modes of operation, selected by the following flags:  
1. `--make_target_bed` - This uses the input biomarker panel to create a bed file with buffer regions, suitable for MinKNOW adaptive sampling. This is different to ONT's "Bed Bugs" tool as it adds the necessary regions for immune infiltrate calculation specific to this tool. You are welcome to double check the output with ONT's tool, accessible [here](https://epi2me.nanoporetech.com/bed-bugs/)
2. `--clin_trial_mode` - This runs sequence data processing, using input bam files, for a single sample. It generates output suitable for classifer training.
3. `--clinical_mode` - This runs sequence data processing using input bam files and a classifier, presumably generated during the above clinical trial. The input biomarker panel must ensure all targets the classifier needs. and these must be captured by adaptive sampling.  

> [!IMPORTANT]   
> A biomarker metadata csv file is required for all three modes of operation. See section [Biomarker panel input](#biomarker-panel-input) for details.

### General Input Parameters

<!-- TODO: maybe add a watch path option once we get up to clinical implementation -->
<!-- TODO: add multiple sample processing -->
| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| make_target_bed | boolean | Generate a target bed file for MinKNOW adaptive sampling. | Uses the biomarker panel metadata to produce buffered target regions. | False |
| clin_trial_mode | boolean | Run sample processing for clinical trial mode. | Processes a single patient sample and outputs data suitable for classifier training. | False |
| clinical_mode | boolean | Generate an individual patient sample report using a trained classifier. | Requires a completed clinical trial and classifier results in the panel metadata. | False |
| project_name | string | A project name that will be used for containing all the samples processed during the clinical trial. If using in clinical mode, this will be used for containing all samples processed using the same classifier. | This structure is necessary in order for the biomarker metadata to be processed appropriately. | (Required input for all modes) |
| panel_metadata | string | The path to `<panel_metadata>.csv` | See section [Biomarker panel input](#biomarker-panel-input) for details | (Required input for all modes) |
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
| out_dir | string | Directory for output of all workflow results. |  | <project_name> |

### Mode `make_target_bed`
This mode has two outputs, for use in MinKNOW adaptive sampling. You must use this mode to generate your adaptive sampling bed file, as it combines your specific targets with regions identified for immune deconvolution by methylation. If your bed file does not contain these regions then the tool will not be able to perform immune deconvolution.  
This mode may need to be rerun multiple times in order to create an optimal bed file that covers all regions adequately while also covering a suitable percentage of the genome. It will check if the resulting bed file meets all the requirements by using the `min_genome_coverage`, `max_genome_coverage` and `buffersize_bp` parameters. If the initial check fails, (it will tell you on the terminal) and/or you want different thresholds for these parameters, adjust them as desired and re-run until you get a successful message.  
Please note that MethylCIBERSORT reference data it uses to build the bed file is based on hg38 genome coordinates.  

#### Input 
| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| chrom_sizes_file | string | The path to a chrom_sizes file, generated from the index of the reference genome. | Optional. If not provided, the bundled `resources/hg38_no_alt.chrom_sizes` file is used. It is strongly recommended to generate this from the same reference genome you are using for sequencing — see [Chromosome size file](#chromosome-size-file). | `resources/hg38_no_alt.chrom_sizes` (bundled) |  

*Example run command*  
<!-- TODO: update!!! -->
```sh
nextflow run ../wf-panorama \
    --make_target_bed \
    --project_name 20260418-BCG_on_NMIBC \
    --panel_metadata ../wf-panorama/demo_input/panel_metadata.csv \
    --mCS_cancer_type bladder \
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
| target_bedfile | string | Path to a bed file of target regions to use in the analysis. | Defaults to the bed file generated by `--make_target_bed`. Can be overridden to point to an existing file. | `<out_dir>/minknow_input/<project_name>.targets_for_align.bed` |
| meth_coverage_threshold | integer | Minimum coverage (depth) threshold for a CpG site to be included in methylation analysis and immune infiltrate deconvolution. | The default of 30 is conservatively high for confident results. For development and testing, a lower value (e.g. 20) can be used to include more sites. | 30 |
<!-- | watch_path | boolean | Enable to continuously watch the input directory for new input files. | This option enables the use of Nextflow's directory watching feature to constantly monitor input directories for new files. | False | -->
<!-- | sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a folder containing sub-folders with FASTQ files. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. A `type` column is required for certain workflows and should have the following values; `test_sample`, `positive_control`, `negative_control`, `no_template_control`. An optional `analysis_group` column is used by some workflows to combine the results of multiple samples. If the `analysis_group` column is present, it needs to contain a value for each sample. |  | -->

*Example run command*  
<!-- TODO: update!!! -->
```sh
nextflow run ../wf-panorama \
    --clin_trial_mode \
    --project_name 20260418-BCG_on_NMIBC \
    --sample Test1 \
    --bam_directory /path/to/bam_pass/ \
    --panel_metadata ../wf-panorama/demo_input/panel_metadata.csv \
    --meth_coverage_threshold 5 \
    --cibersortx_username person@place.com \
    --cibersortx_token 12a3bc  \
    -profile slurm,singularity \
    -resume
```

### Mode `clinical_mode`

#### Input
| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample | string | A single sample name or identifier. Must start with an alphabet letter (i.e. not a digit or symbol). | Required for `--clin_trial_mode` and `--clinical_mode`. | null |
| bam_directory | string | The path to a single directory containing bams to process. | Usually the `bam_pass` directory in MinKNOW- or Dorado-processed data. Required for `--clin_trial_mode` and `--clinical_mode`. | null |
| cibersortx_username | string | CIBERSORTx account username (email). | Required to authenticate with the CIBERSORTx service. Required for `--clin_trial_mode` and `--clinical_mode`. | null |
| cibersortx_token | string | CIBERSORTx authentication token. | Obtain from the CIBERSORTx Downloads page. Required for `--clin_trial_mode` and `--clinical_mode`. | null |
| mCS_cancer_type | string | Cancer type signature to use for MethylCIBERSORT deconvolution. | Must match one of the available signature sets bundled with MethylCIBERSORT. See `nextflow_schema.json` for valid options. Required for `--clin_trial_mode` and `--clinical_mode`. | bladder |
| target_bedfile | string | Path to a bed file of target regions to use in the analysis. | Defaults to the bed file generated by `--make_target_bed`. Can be overridden to point to an existing file. | `<out_dir>/minknow_input/<project_name>.targets_for_align.bed` |
| meth_coverage_threshold | integer | Minimum coverage (depth) threshold for a CpG site to be included in methylation analysis and immune infiltrate deconvolution. | The default of 30 is conservatively high for confident results. For development and testing, a lower value (e.g. 20) can be used to include more sites. | 30 |
| report_title | string | Title to display on the generated clinical report. | Sets the title heading used in the per-sample clinical report. Defaults to `--project_name` if not provided. | null (uses project_name) |
<!-- | watch_path | boolean | Enable to continuously watch the input directory for new input files. | This option enables the use of Nextflow's directory watching feature to constantly monitor input directories for new files. | False | -->

*Example run command*  
<!-- TODO: update!!! -->
```sh
nextflow run ../wf-panorama \
    --clinical_mode \
    --project_name 20260418-BCG_on_NMIBC \
    --report_title "BCG on NMIBC – Patient Report" \
    --sample Test1 \
    --bam_directory /path/to/bam_pass/ \
    --panel_metadata ../wf-panorama/demo_input/panel_metadata.csv \
    --meth_coverage_threshold 5 \
    --cibersortx_username person@place.com \
    --cibersortx_token 12a3bc  \
    -profile slurm,singularity \
    -resume
```

### Outputs for clin_trial and clinical modes
Raw results for epi2melabs wf-human-variation, SNV, methylation and immune infiltrate sections can be found in their respective directories.  
In `clin_trial_mode` the outputs from these sections (excluding wf-human-variation) are collated into one csv that can be used as input to a machine learning classifier.  
In `clinical_mode` they are processed into a pdf report that uses `resources/template.md` to make.  

#### Output directory `<sample name>/wf-humvar-run`
Along with nextflow-generated directories (that can be removed when finished to save space), there are the following files:                     
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

<!-- #### Output directory `<sample name>/sv_annotation`
This module is not yet implemetned, there is currently an empty placeholder output -->
#### Output directory `<sample name>/methylCS`
<!-- TODO: put this output in immune infiltrate -->

- CIBERSORTx_sample_Results.csv
- sample.CS_mix_matrix.txt
- sample.bladder.mCS_ref.txt

#### Output directory `<sample name>/immune_infiltrate`
- sample.immune_panel_results.csv


## 🎯 Biomarker panel input
This is a csv file containing metadata for each biomarker.  
This can be generated using the template excel file provided in [demo_input](demo_input/panel_metadata_template.xlsx), then "Save As" a csv. Ensure the ID column remains 3 digits long (it can be string/character format, not number).  
An example of the csv can be found here: [demo_input/panel_metadata.csv](demo_input/panel_metadata.csv)  

> [!IMPORTANT]   
> * ID numbers in the 400's are reserved for immune parameters
> * ID numbers in the 500's are reserved for additional non-molecular factors such as demographic or clinicopathologic indicators that you wish to include in the classifiers but cannot be measured by nanopore sequencing
> * The biomarker types "immune_inf" must not be changed  

<!-- TODO: Add references to the clinical report -->

| Parameter name                      | Type    | Required? | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
|-------------------------------------|---------|-----------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| ID                                  | string  | Required  | A unique ID number 3 characters long (e.g. 001).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| Biomarker Name                      | string  | Required  | A name suitable for the biomarker. Initially used for gene names. Can be used for multiple biomarkers.                                                                                                                                                                                                                                                                                                                                                                                                                          |
| Biomarker Type                      | string  | Required  | One of: snv, sv, mod, area_mutations, expression, exp_ratio, immune_ratio, immune_inf, microsatellite, demographic, clinicopathology.                                                                                                                                                                                                                                                                                                                                                                                           |
| Panel or Area of Interest?          | string  | Required  | Is this part of the Panel or is it an extra region (Area of Interest) that's been individually requested? Options are "Panel" or "AOI". AOI is irrelevant during clinical trial, everything is part of the "Panel".                                                                                                                                                                                                                                                                                                                                                                                          |
| chrom, start pos, end pos           | strings | Required  | Genome coordinates of the desired area. Entries in "chrom" column must be in the format "chr1". Columns "start pos" and "end pos" handle both comma separated numbers and normal numbers.                                                                                                                                                                                                                                                                                                                                                                                         |
| length                              | number  | Optional  | Length of the genomic region. Will auto calculate.                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| strand                              | string  | Optional  | Which strand is the feature of interest on? "+" or "-". Required if the answer for "Is this record the whole gene" is "Yes".                                                                                                                                                                                                                                                                                                                                                            |
| Scoring Type                        | string  | Required  | Options are: "genotypic", "continuous" or "categorical". Others may be added in the future. If this is not included there will be no automatic analysis as part of the panel.                                                                                                                                                                                                                                                                    |
| Result Options                      | string  | Required  | All possible results. For genotypic data, single genotype in the format “Allele 1\|Allele 2”, or multiple genotypes separated by a “/“ character “Allele 1\|Allele 2/Allele 1\|Allele 2/Allele 1\|Allele 2”. For continuous data, a range in the format: 0.0-1.0. For categorical data, the options separated by "/". If this is not included there will be no automatic analysis as part of the panel.                                             |
| Result                              | string  | Optional  | Available for Clinicopathologic and Demographic results.                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| Notes                               | string  | Optional  | Any notes the user wants to add.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| References                          | string  | Optional  | The references where the biomarker was found. Will be added to the final report in the future.                                                                                                                                                                                                                                                                                                                                                                                          |
| Is variant in coding region? (snv)   | string  | Required  | "Yes" or "No".                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| SNP ID (snv)                        | string  | Optional  | rs ID number.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| Is variant in coding region? (mod)   | string  | Required  | "Yes" or "No".                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| Is this record the whole gene? (mod) | string  | Required  | "Yes" or "No". If "Yes", the target region will be extended to include promoter (2000 bp) and downstream (1000 bp). Also the "strand" is required.                                                                                                                                                                                                                                                                                                                                     |
| Illumina EPIC ID (mod)               | string  | Optional  | The EPIC id associated with the modification site.                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| DNA methylation region (mod)         | string  | Required  | Options are "position" (for a single site), "promoter", "intragenic", or "downstream".*                                                                                                                                                                                                                                                                                                                                                                                                |
| Is variant in coding region? (area_mutations)   | string  | Optional  | "Yes" or "No".                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| Is this record the whole gene? (area_mutations) | string  | Required  | "Yes" or "No". If "Yes", the target region will be extended to include promoter (2000 bp) and downstream (1000 bp). Also the "strand" is required.                                                                                                                                                                                                                                                       |
| Expression Ratio Components (exp_ratio)         | string  | Required  | Each component of the ratio has a separate biomarker panel ID number, and the “Biomarker Name” must match the entry in “Expression Ratio Components”. Format is "Biomarker Name 1"/"Biomarker Name 2".                                                                                                                     |

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
  - https://github.com/lucy924
  - lucy.picard@otago.ac.nz
- Aaron Stevens
  - PI, Genomic and Epigenomic Research Lab (GERL), University of Otago, Wellington, NZ 🇳🇿
  - aaron.stevens@otago.ac.nz

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
<!-- TODO: update to more recent version, investigate using wf-somatic-variation instead -->


## 🗺️ Roadmap
- [ ]  Add multiple treatment options
- [ ]  Integrate with Epi2ME Labs
- [ ]  Add option for custom immune infiltrate references

## 📜 Pipeline History

Panorama was initially designed as part of the Doctoral thesis entitled:  
**Epigenetic Consequences of BCG Immunotherapy In Bladder Cancer**  
For full details on the development of this pipeline the relevant chapter is "Chapter Four: Multi-biomarker discovery using nanopore sequencing technology: proof-of-concept" and can be found in the [Otago archives](https://hdl.handle.net/10523/48489). The original development of the pipeline used [Snakemake](https://snakemake.readthedocs.io/en/stable/) and can be found on [Lucy's github repo](https://github.com/lucy924/nanopore_multiBM_pipeline) with associated material [here](https://github.com/lucy924/Multi-biomarker-ONT-project).

## 🎉 Acknowledgements
Many thanks go to the funders of this project, The Barbara Basham Medical Charitable Trust. Read about the origin of the trust [here](https://wellington.govt.nz/arts-and-culture/heritage/historic-public-memorials/aunt-daisy).

## 📚 References
- CIBERSORT
  - Newman, A.M., Liu, C.L., Green, M.R., Gentles, A.J., Feng, W., Xu, Y., Hoang, C.D., Diehn, M., Alizadeh, A.A., 2015. Robust enumeration of cell subsets from tissue expression profiles. Nat Methods 12, 453–457. https://doi.org/10.1038/nmeth.3337
- MethylCIBERSORT
  - Chakravarthy, A., Furness, A., Joshi, K., Ghorani, E., Ford, K., Ward, M.J., King, E.V., Lechner, M., Marafioti, T., Quezada, S.A., Thomas, G.J., Feber, A., Fenton, T.R., 2018. Pan-cancer deconvolution of tumour composition using DNA methylation. Nat Commun 9, 3220. https://doi.org/10.1038/s41467-018-05570-1
<!-- Citation style: Elsevier (author-date/Harvard, with titles) -->
- Adaptive sampling reference file modified from:
  - Stephane Plaisance (VIB-NC) 2021
  - https://github.com/Nucleomics-VIB
