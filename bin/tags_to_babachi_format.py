import pysam
import sys
import os

from count_tags_pileup import SNV


def main(in_file):
    name = os.path.basename(in_file).split('.')[0]
    with pysam.TabixFile(in_file) as vars_file:
        for line in vars_file.fetch():
            chrom, start, end, ID, ref, alt, AAF, RAF, GT, n_ref, n_alt, n_original_reads, n_failed_mapping, n_failed_genotyping = line.strip('\n').split('\t')
            variant = SNV([chrom, start, end, ID, ref, alt, AAF, RAF, GT])
            if not variant.is_het:
                continue
            n_ref, n_alt, n_original_reads, n_failed_mapping, n_failed_genotyping, n_failed_bias = map(int, [n_ref, n_alt, n_original_reads, n_failed_mapping, n_failed_genotyping])
            assert n_original_reads == n_alt + n_ref + n_failed_bias + n_failed_genotyping + n_failed_mapping

            print('\t'.join(map(str, [chrom, start, end, ID, ref, alt, n_ref, n_alt, name, variant.aaf, variant.raf, n_failed_mapping / n_original_reads if n_original_reads > 0 else 0])))


if __name__ == '__main__':
    main(sys.argv[1])
