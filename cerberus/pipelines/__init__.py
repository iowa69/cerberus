"""Pipeline modules compose stages into mode-specific flows.

A pipeline takes the CerberusConfig + tuned params + the QC output, and
returns a PipelineResult describing the files it produced and per-stage
read counts.
"""
