# Final module map

Implemented modules in wf-panorama:
- panel_prep (modules/panel_prep.nf)
- sample_processing (modules/sample_processing.nf)

Implemented CLI scripts are in wf-panorama/bin/ and are standalone adaptations of the original Snakemake scripts.

Notes and differences:
- Nested Nextflow runs (epi2me-labs/wf-human-variation and wf-basecalling) are invoked via shell calls; future work could integrate them as true subworkflows.
- Scoring and modification calling were simplified to produce placeholder outputs; these can be enhanced to exactly replicate Snakemake behavior if needed.
- Containers are placeholders and must be built; no container images were created in this phase.

