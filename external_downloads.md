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
