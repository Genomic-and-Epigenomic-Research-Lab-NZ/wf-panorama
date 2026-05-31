This directory is where the following external resources should live:
- general container `general.sif`
- methylcibersort container `methylcibersort.sif`
- cibersortx_fractions container `cibersortx_fractions.sif`

These pre-built containers can be obtained from figshare.com:  
https://figshare.com/account/articles/32165103
<!-- TODO: update once containers are published -->

The apptainer def files for `general.sif` and `methylcibersort.sif` are included in this directory in case you want to or need to rebuild from scratch.  
We do not recommend rebuilding the methylcibersort.sif container, this could be tricky to get working.

The CIBERSORTx container was built using the following command (2026-05-27):

```sh
apptainer pull cibersortx_fractions.sif docker://cibersortx/fractions
```
