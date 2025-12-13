#!/usr/bin/env python3
# Jeff Vierstra 2018
# TODO:
# --add filters/etc. as option
import sys
import logging
import pandas as pd
from argparse import ArgumentParser
from tqdm import tqdm

from genome_tools.data.extractors import AllelicReadsExtractor
from genome_tools import df_to_intervals, VariantInterval
from genome_tools.helpers import df_to_tabix

logging.basicConfig(stream = sys.stderr, level='WARNING')


def transform_af(value):
	if value == '.' or value == 'None' or value is None or pd.isna(value) or value == '':
		return None
	else:
		return float(value)


def get_total_counts(all_variants, bam_file):
    all_variants = df_to_intervals(all_variants)
    result = []
    for variant in tqdm(all_variants):
        with AllelicReadsExtractor(bam_file) as read_extractor:
            read_pairs = read_extractor.extract_allelic_reads(variant)
            n_remapped_reads = len(read_pairs)
            result.append(n_remapped_reads)
			
    return pd.Series(result)


def get_allelic_counts(all_variants, bam_file):
    all_variants = df_to_intervals(all_variants)
    vals = []
    for variant in tqdm(all_variants):
        variant: VariantInterval
        counter = {
            'n_ref': 0,
            'n_alt': 0,
            'deletion_or_refskip': 0,
            'bias': 0,
            'genotype_error': 0,
            'other': 0,
        }
        with AllelicReadsExtractor(bam_file) as read_extractor:
            read_pairs = read_extractor.extract_allelic_reads(variant)

            for _, read_info in read_pairs.items():
                status = read_info['status']

                allele = read_info['read_allele']
                if allele == variant.ref:
                    counter['n_ref'] += 1
                elif allele == variant.alt:
                    counter['n_alt'] += 1

                if status == 'valid':
                    continue
                if status not in counter:
                    status = 'other'
                counter[status] += 1

            vals.append(counter)
    return pd.DataFrame(vals)


def main(all_variants, cram_file, dedup_vals) -> pd.DataFrame:
    allelic_counts_df = get_allelic_counts(all_variants, cram_file)
    all_variants['n_original_reads'] = dedup_vals
    all_variants[allelic_counts_df.columns] = allelic_counts_df.values
    return all_variants

   
if __name__ == "__main__":
    parser = ArgumentParser(description = "Count tags by allele")

    parser.add_argument("var_file", metavar = "var_file", type = str,
                        help = "Path to variant file (must have corresponding index)")

    parser.add_argument("cram_file", metavar = "cram_file", type = str, 
                        help = "Path to BAM-format tag sequence file")
    
    parser.add_argument("output_file", metavar = "output_file", type = str, 
                        help = "Path to output file")
    
    parser.add_argument("--original_dedup_cover", metavar = "counts from original deduped bam file",
                        type = str, default=None,
                        help = "Path to bed format file with total read counts")
    
    parser.add_argument('--only_coverage', default=False, action="store_true",
                        help="Specify to emit only coverage for the SNPs")

    args = parser.parse_args()
    all_variants = pd.read_table(args.var_file)
    cram_file = args.cram_file
    if args.only_coverage:
        result = get_total_counts(all_variants, cram_file).to_frame(name='coverage')
        result.to_csv(
            sys.stdout,
            sep='\t',
            index=False
        )
    else:
        dedup_df = pd.read_table(args.original_dedup_cover)['coverage'].values
        result = main(
            all_variants,
            cram_file,
            dedup_df
        )
        df_to_tabix(
            result,
            sys.stdout
        )

