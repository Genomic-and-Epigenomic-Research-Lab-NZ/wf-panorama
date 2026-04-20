# Panorama

## Introduction

This is a bioinformatic pipeline utilising the capabilities of nanopore sequencing ([ONT](https://nanoporetech.com/)) to combine multiple biomarkers of different sources; specifically, mutations, methylation and tumour immune infiltrate. This can enable prediction of drug compatibility in cancer tumours. This pipeline was developed as a clinical bioinformatic workflow that requires little bioinformatic expertise to use. The results of such a tool, with the appropriate pre-clinical trial, can be integrated as part of a protocol aimed to assist molecular pathologists and clinicians to swiftly develop personalised treatment plans.  

Panorama is built using the Nextflow workflow language and is intended to be used in the [Epi2ME](https://github.com/epi2me-labs) framework [Oxford Nanopore Technologies](https://community.nanoporetech.com).

### Panorama is intended for use as follows:
#### Stage One: In a Pre-Clinical Trial for a defined disease and single* treatment option
   1. Biomarkers are input via a csv file, as detailed in the section [Biomarker panel input](#biomarker-panel-input) below. This will form the basis of the bed file used for adaptive sampling, and is also used for downstream pipeline processes. Biomarkers can include SNVs, methylation markers, and certain immune infiltrate markers as determined by immune deconvolution using methylation markers.
   2. Tumour samples are nanopore sequenced using adaptive sampling and a bed file containing the target regions required for biomarker analysis. This bed file is produced by this pipeline with the flag `--make_target_bed` and the biomarker csv file using the input parameter `--panel_metadata`. Resulting bam files are used as input to Panorama. Bam files must be basecalled using modified base calling (5mC+5hmC contexts), and aligned to the hg38 genome.
   3. In the default mode, the pipeline will output sample data as a csv file corresponding to the input biomarkers.
   4. Following patient sequence data collection, the researchers will carry out their own classifier development using the biomarker results identified by this pipeline. This is outside this tool's scope and presumably will involve some kind of machine learning. Results are added to the original biomarker input csv file and used as input for clinical reporting.

**Multiple treatment options will be available in a future update*

#### Stage Two: Patient sample reporting following a successful clinical trial
   1. Using the results of the pre-clinical trial described above, and the final version of the biomarker input with classifier results, a single patient sample can be submitted to this pipeline as before, and will output a report based on the clinical trial results.  
   
Please note that Stage Two has not been fully tested as there has been no such clinical trial carried out yet. This is expected to be finalised during/after an actual clinical trial. Please get in contact with the developers if you decide to use this tool in a clinical trial and we will assist/collaborate with trial design and clinical reporting development.
<!-- eg read depth requirements -->


## Compute requirements

Recommended requirements:

+ CPUs = 32
+ Memory = 128GB

Minimum requirements:

+ CPUs = 16
+ Memory = 32GB

*Based off [wf-human-variation](https://github.com/epi2me-labs/wf-human-variation) as this is the most computationally heavy component of the workflow.*

## Install and run

These are instructions to install and run the workflow on command line.
<!-- TODO: check if this is doable -->
<!-- You can also access the workflow via the
[EPI2ME Desktop application](https://labs.epi2me.io/downloads/). -->

The workflow uses [Nextflow](https://www.nextflow.io/) to manage
compute and software resources,
therefore Nextflow will need to be
installed before attempting to run the workflow.

The workflow can currently be run using [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html)
to provide isolation of the required software.
This is automated provided Singularity is installed.
This is controlled by the
[`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles)
parameter as exemplified below.

It is not required to clone or download the git repository
in order to run the workflow.
More information on running EPI2ME workflows can
be found on the [EPI2ME website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow.
This will pull the repository in to the assets folder of
Nextflow and provide a list of all parameters
available for the workflow as well as an example command:

```
nextflow run lucy924/wf-panorama --help
```
To update a workflow to the latest version on the command line use
the following command:
```
nextflow pull lucy924/wf-panorama --help
```

<!-- A demo dataset is provided for testing of the workflow.
It can be downloaded and unpacked using the following commands:
```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-template/wf-template-demo.tar.gz
tar -xzvf wf-template-demo.tar.gz
```
The workflow can then be run with the downloaded demo data using: -->

The workflow can be run using:
<!-- TODO: fix profile options (epi2me uses "standard") -->
```
nextflow run lucy924/wf-panorama \
    --project_name Project-BCG_on_NMIBC \
    --sample Test1 \
    --bam_directory /path/to/passed_bams \
    --panel_metadata /path/to/demo_input/panel_metadata.csv \
    --target_bedfile /path/to/demo_input/targets.bed \
    -profile slurm,singularity
```

For further information about running a workflow on
the command line see https://labs.epi2me.io/wfquickstart/


## Related protocols

<!---Hyperlinks to any related protocols that are directly related to this workflow, check the community for any such protocols.--->

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices.
This protocol currently uses epi2me-labs/wf-human-variation v2.6.0 for initial analysis of bam files.  
<!-- TODO: update to more recent version, investigate using wf-somatic-variation instead -->

Find related protocols in the [Nanopore community](https://community.nanoporetech.com/docs/).


## Input example
There are three modes of operation, selected by the following flags:  
1. `--make_target_bed` - This uses the input biomarker panel to create a bed file with buffer regions, suitable for MinKNOW adaptive sampling. This is different to ONT's "Bed Bugs" tool as it adds the necessary regions for immune infiltrate calculation specific to this tool. You are welcome to double check the output with ONT's tool, accessible [here](https://epi2me.nanoporetech.com/bed-bugs/)
2. `--clin_trial_mode` - This runs sequence data processing, using input bam files, for a single sample. It generates output suitable for classifer training.
3. `--clinical_mode` - This runs sequence data processing using input bam files and a classifier, presumably generated during the above clinical trial. The input biomarker panel must ensure all targets the classifier needs. and these must be captured by adaptive sampling.
A biomarker metadata csv file is required for all three modes of operation. This can be generated using the template excel file provided in [demo_input](demo_input/panel_metadata_template.xlsx), then "Save As" a csv. Ensure the ID column remains 3 digits long.
Notes:  
* ID numbers in the 400's are reserved for immune parameters
* ID numbers in the 500's are reserved for additional non-molecular factors such as demographic or clinicopathologic indicators that you wish to include in the classifiers but cannot be measured by nanopore sequencing
* The biomarker types "immune_inf" must not be changed
The  accepts a single folder containing BAM files as input.  

## Input parameters

### Input Options

<!-- TODO: maybe add a watch path option once we get up to clinical implementation -->
<!-- TODO: add multiple sample processing -->
| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| make_target_bed | boolean |  |  | False |
| clin_trial_mode | boolean |  |  | False |
| clinical_mode | boolean |  |  | False |
| project_name | string | A project name that will be used for containing all the samples processed during the clinical trial. If using in clinical mode, this will be used for containing all samples processed using the same classifier. | This structure is necessary in order for the biomarker metadata to be processed appropriately. | (Required input for all modes) |
| panel_metadata | string | The path to `<panel_metadata>.csv` | See section [Biomarker panel input](#biomarker-panel-input) for details | (Required input for all modes) |
| sample | string | A single sample name or identifier. Must start with an alphabet letter (i.e. not a digit or symbol). |  | False (Required input for `--clin_trial_mode` and `--clinical_mode`) |
| bam_directory | string | The path to a directory containing bams to process. | Usually the `bam_pass` directory in MinKNOW- or Dorado-processed data. | False (Required input for `--clin_trial_mode` and `--clinical_mode`) |
<!-- | watch_path | boolean | Enable to continuously watch the input directory for new input files. | This option enables the use of Nextflow’s directory watching feature to constantly monitor input directories for new files. | False | -->
<!-- | sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a folder containing sub-folders with FASTQ files. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. A `type` column is required for certain workflows and should have the following values; `test_sample`, `positive_control`, `negative_control`, `no_template_control`. An optional `analysis_group` column is used by some workflows to combine the results of multiple samples. If the `analysis_group` column is present, it needs to contain a value for each sample. |  | -->


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all workflow results. |  | <project_name> |



## Outputs

### Mode `make_target_bed`
This mode has two outputs, for use in MinKNOW adaptive sampling. You must use this mode to generate your adaptive sampling bed file, as it combines your specific targets with regions identified for immune deconvolution by methylation. If your bed file does not contain these regions then the tool will not be able to perform immune deconvolution.  
This mode may need to be rerun in order to create an optimal bed file that covers all regions adequately while also covering a suitable percentage of the genome. It will check if the resulting bed file meets all the requirements by using the `min_genome_coverage`, `max_genome_coverage` and `buffersize_bp` parameters. If the initial check fails, (it will tell you on the terminal) and/or you want different thresholds for these parameters, adjust them as desired and re-run until you get a successful message.

| Title | File path | Description | 
|-------|-----------|-------------| 
| Targets with buffered regions | ./minknow_input/<project_name>.targets_buffed.bed | bed file for adaptive sampling |
| Targets for alignment stats | ./minknow_input/<project_name>.targets_for_align.bed | A bed file provided for optional alignment of your target regions only. This file does NOT have buffered regions, do not use it in the adaptive sampling input. You may use it in the alignment only section, and it can help monitor read depth in your desired regions. If you are not confident with this do NOT use it. |

### Mode `clin_trial_mode`
This mode generates a lot of files. Raw results for each SNV, Methylation and immune infiltrate files can be found in their respective sections. These are collated into one csv that can be used as input to a machine learning classifier.

| Title | File path | Description | 
|-------|-----------|-------------| 
| Targets with buffered regions | ./minknow_input/targets_for_.bed | bed file for adaptive sampling |

## Biomarker panel input
This is a csv file containing metadata for each biomarker.  
Note that an example can be found here: [demo_input/panel_metadata.csv](demo_input/panel_metadata.csv)  

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

## Target bed file
This file is generated using the `--make_target_bed` flag and requires the biomarker panel csv `--panel_metadata` and a project name `--project_name` only.


## Pipeline overview

<!---High level numbered list of main steps of the workflow and hyperlink to any tools used. If multiple workflows/different modes perhaps have subheadings and numbered steps. Use nested numbering or bullets where required.--->
### 1. Concatenates input files and generate per read stats.

The [fastcat/bamstats](https://github.com/epi2me-labs/fastcat) tool is used to concatenate multifile samples to be processed by the workflow. It will also output per read stats including average read lengths and qualities.



## Troubleshooting

<!---Any additional tips.--->
+ If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ See how to interpret some common nextflow exit codes [here](https://labs.epi2me.io/trouble-shooting/).



## FAQ's

<!---Frequently asked questions, pose any known limitations as FAQ's.--->

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-template/issues) page or start a discussion on the [community](https://community.nanoporetech.com/).



## Related blog posts

+ [Importing third-party workflows into EPI2ME Labs](https://labs.epi2me.io/nexflow-for-epi2melabs/)

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.

README
======

wf-panorama - Nextflow conversion of nanopore_multiBM_pipeline

Quick start

- Provide config/panel_metadata.csv (kept as CSV) and a samplesheet at config/samplesheet.csv or pass -params.sample and -params.project_name.
- Ensure bin/ scripts are executable: chmod +x bin/*
- Provide Apptainer/Singularity images for the containers referenced in nextflow.config and modules (place in accessible registry or local .sif paths).

Run locally:

nextflow run main.nf -profile local -params-file config/config.yaml

Run on SLURM:

nextflow run main.nf -profile slurm -params-file config/config.yaml

Notes
- The pipeline uses DSL2 modules in modules/ and CLI scripts in bin/ adapted from the original Snakemake pipeline.
- panel_metadata.csv remains a CSV and is passed to modules that require it.
- Containers are placeholders; build Apptainer images later and update nextflow.config or process directives. TODO

## Acknowledgements
Many thanks go to the funders of this project, The Barbara Basham Medical Charitable Trust. Read about the origin of the trust [here](https://wellington.govt.nz/arts-and-culture/heritage/historic-public-memorials/aunt-daisy).

## Pipeline History

Panorama was initially conceived as part of the Doctoral thesis entitled:  
**Epigenetic Consequences of BCG Immunotherapy In Bladder Cancer**  
For full details on the development of this pipeline the relevant chapter is "Chapter Four: Multi-biomarker discovery using nanopore sequencing technology: proof-of-concept" and can be found in the [Otago archives](https://hdl.handle.net/10523/48489). The original development of the pipeline used [Snakemake](https://snakemake.readthedocs.io/en/stable/) and can be found on [Lucy's github repo](https://github.com/lucy924/nanopore_multiBM_pipeline) with associated material [here](https://github.com/lucy924/Multi-biomarker-ONT-project).
