# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).  
<!-- Version format of X.Y.Z (Major.Minor.Patch). Bug fixes not affecting the API increment the patch version, backward compatible API additions/changes increment the minor version, and backward incompatible API changes increment the major version. -->

## [v0.2.0] - 2026-10-16

Main goal of this change: To enable biomarkers composed of multiple measurements, e.g. the Hormone Receptor Positive (HR+) category in Breast Cancer, and adding a "Molecular Groups" category to the input biomarker panel file which enables the assigning of a sequenced tumour into an established group.  
Note that the Molecular Groups category may not be suitable for input into machine learning algorithms, and is primarily intended for enabling comparison to pre-existing clinical categorisations. 

### Added
- Added support for genome-coordinate generation from a list of genes.
- Added the option to include the full EPIC array genomic locations in panel BED generation.
- Added an msk-impact panel option for panel BED generation.
- Added SV module support, including fusion-related `NOTES` output fields.
- Added plotting support for immune infiltrate outputs.

### Changed
- Updated result processing and presentation workflows, including `collate_results` and composite results handling.
- Restructured result output directories to better separate project and sample outputs, enabling concurrent sample runs without file conflicts.
- Updated `snv_annotation` and `sv_annotation` to use revised shared variables and improved result handling.
- Improved panel preparation workflows, including exact gene-name filtering for genome coordinates and updated coverage checks.
- Updated local execution and profile configuration to support Singularity/Apptainer and newer workflow/tool versions.
- Improved handling of `wf-human-variation` acquisition to avoid path issues on HPC systems.
- Renamed the final output parameter to `--out_final`.
- Updated wf-human-variation version used to v2.7.2

### Fixed
- Fixed crashes in `snv_annotation.py` when no panel SNVs are found.
- Fixed missing index handling in `sample_processing.nf`.
- Fixed a hang in `snv_prep` when input VCFs are already compressed.
- Fixed hardcoded testing paths in `modification_calling`.
- Fixed `sv_annotation` to report `absent` when no variation is detected.
- Fixed the `epic_locs` parameter and related panel template issues.
- Fixed `modification_calling` and other workflow integration issues.
- Fixed colour consistency in immune infiltrate plotting.

### Deprecated
- None.

### Removed
- Removed resource constraints from `RUN_WF_HUMVAR` due to crashes.


## [v0.1.1] - 2026-07-17

### Fixed

- Fix to modification calling section, large portion of script had been deleted at some point so it only output the header. Script is now restored.
- README improvements.
- snv_prep resources now dynamically adjusted
- Fixed basecaller override passing through to wf-human-variation

## [v0.1.0] - 2026-06-02

### Added

- Initial public release of wf-panorama.

### Notes

- First tagged release for the project.


## [v0.0.1] - 2026-05-11

- Pre-release.
