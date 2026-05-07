This directory comes with the following resources:  
- tandem repeat reference file for alignment `hg38.trf.bed`
- chrom_sizes file made by make_adaptive_ref from the genome index file `hg38_no_alt.chrom_sizes`
- EPIC probe coordinates used in immune infiltrate section: `IlluminaEPIC_genomic_locations_hg38.csv`
- template file for the clinical report `template.md`
- tab-separated file of the available methylation signatures from MethylCIBERSORT: `methylcibersort_signatures.tsv`

This directory is where the following external resources should live:
- human genome reference `GCA_000001405.15_GRCh38_no_alt_analysis_set.fna`
- human genome reference index` GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai`

If using a different human genome build, these resources will need to be updated:  
- tandem repeat reference file for alignment `hg38.trf.bed`
- chrom_sizes file made by make_adaptive_ref from the genome index file `hg38_no_alt.chrom_sizes`
- EPIC probe coordinates used in immune infiltrate section: `IlluminaEPIC_genomic_locations_hg38.csv`

Resources to update:
- Flow charts showing how this pipeline works, and suggested workflow for pre-clinical trial (`SectionFlowCharts`)
- `MinKNOW_settings` directory containing images of Recommended MinKNOW settings to use (older version of MinKNOW)
  
Resources to remove:
- CIBERSORT setup records (`CIBERSORT`) 
- `MethylCIBERSORT_Release` directory containing raw MethylCIBERSORT tar file
- bladder reference from methatlas: `ref_atlas_bladder.csv`

