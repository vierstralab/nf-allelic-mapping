# nf-allelic-mapping

## Run
```nextflow run /home/sabramov/projects/ENCODE4/nf-allelic-mapping/main.nf -profile Altius -resume```
## Output format
Upon completion pipeline creates 3 folders
- `target_variants` - consists of per-individual files with all SNVs passing thresholds specified in params.config. Each file has the following columns:
`#chr, start, end, ID, ref, alt, AAF, RAF, GT, GQ, DP, ref_counts, alt_counts`

- `remapped_files` - per-sample files with all WASP analyzed variants. Each file has the following columns:
`#chr, start, end, ID, ref, alt, AAF, RAF, GT, ref_counts_after_wasp, alt_counts_after_wasp, initial_total_count, n_failed_mapping, n_failed_genotyping, n_failed_bias` 

- `babachi_files` - per-sample `count_reads` files in format required for BABACHI. Includes variants passing thresholds specified in params.config. Each file has the following columns:
`#chr, start, end, ID, ref, alt, ref_counts_after_wasp, alt_counts_after_wasp, sample_id, AAF, RAF, FMR`



## Notes
`filter_reads.py` file was forked from https://github.com/StamLab/stampipes/blob/main/scripts/bwa/filter_reads.py
