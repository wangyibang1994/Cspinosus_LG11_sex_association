# Figure workflow status

All numerical figure inputs are deposited in `source_data/`, and their central counts are tested by `workflow/reproduce_key_numbers.py`. The exact historical R renderers are retained under `workflow/legacy_lsf/historical_figures/` because they still contain workstation-specific paths and a local panel-alignment helper.

A formal `v1.0.0` release requires portable renderers that read only repository-relative source data. Until those are written and verified against the submitted PDF/SVG/TIFF files, figure regeneration is a documented release blocker rather than a falsely claimed capability.
