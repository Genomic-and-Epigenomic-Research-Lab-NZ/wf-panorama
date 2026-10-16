This directory is where the following external resources should live:
- general container `general.sif`
- methylcibersort container `methylcibersort.sif`
- cibersortx_fractions container `cibersortx_fractions.sif`

These pre-built containers can be obtained from figshare.com:  
https://figshare.com/account/articles/32165103
<!-- TODO: update once containers are published -->

To get the containers:  
1. `wget --content-disposition https://ndownloader.figshare.com/articles/32165103/versions/[x]`  
2. `unzip 32165103.zip`

Note the "ndownloader" is BEFORE "figshare", unlike the link that is copied from the website   

Also note that you will want the latest version in the place of `[x]` in the address above

Thanks to https://github.com/Shaoyi-Zhang96/Figshare_wget_download for instructions


The apptainer def files for `general.sif` and `methylcibersort.sif` are included in this directory in case you want to or need to rebuild from scratch.  
Terminal command:  
```sh
cd /path/to/wf-panorama/containers
apptainer build general.sif general.def
```

We do not recommend rebuilding the methylcibersort.sif container, this could be tricky to get working.

The CIBERSORTx container was built using the following command (2026-05-27):

```sh
apptainer pull cibersortx_fractions.sif docker://cibersortx/fractions
```
