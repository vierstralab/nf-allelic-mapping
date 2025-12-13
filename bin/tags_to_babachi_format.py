import sys
import pandas as pd
import os



def main(in_file):
    name = os.path.basename(in_file).split('.')[0]
    variants = pd.read_table(in_file)
    variants['name'] = name
    variants['fmr'] = variants['n_failed_mapping'] / variants['n_original_reads'].replace(0, pd.NA)

    return variants[['chrom', 'start', 'end', 'ID', 'ref', 'alt', 'n_ref', 'n_alt', 'name', 'AAF', 'RAF',]]
    for line in vars_file.fetch():
        chrom, start, end, ID, ref, alt, AAF, RAF, GT, n_ref, n_alt, n_original_reads, n_failed_mapping, n_failed_genotyping, n_failed_bias = line.strip('\n').split('\t')
        variant = SNV([chrom, start, end, ID, ref, alt, AAF, RAF, GT])
        if not variant.is_het:
            continue
        n_ref, n_alt, n_original_reads, n_failed_mapping, n_failed_genotyping, n_failed_bias = map(int, [n_ref, n_alt, n_original_reads, n_failed_mapping, n_failed_genotyping, n_failed_bias])
        assert n_original_reads == n_alt + n_ref + n_failed_bias + n_failed_genotyping + n_failed_mapping

        print('\t'.join(map(str, [chrom, start, end, ID, ref, alt, n_ref, n_alt, name, variant.aaf, variant.raf, n_failed_mapping / n_original_reads if n_original_reads > 0 else 0])))


if __name__ == '__main__':
    main(sys.argv[1]).to_csv(sys.stdout, sep='\t', index=False, header=False)
