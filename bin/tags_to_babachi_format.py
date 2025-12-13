import sys
import pandas as pd
import os
import numpy as np


def main(in_file):
    name = os.path.basename(in_file).split('.')[0]
    variants = pd.read_table(in_file)
    variants['name'] = name
    variants['n_failed_mapping'] = variants.eval(
        'n_original_reads - n_ref - n_alt - deletion_or_refskip - bias - genotype_error - other'
    )
    variants['fmr'] = np.where(
        variants['n_original_reads'] > 0,
        variants.eval('n_failed_mapping / n_original_reads'),
        0
    )

    return variants[['#chr', 'start', 'end', 'ID', 'ref', 'alt', 'n_ref', 'n_alt', 'name', 'RAF', 'AAF', 'fmr']]


if __name__ == '__main__':
    data = main(sys.argv[1])
    data.to_csv(sys.argv[2], sep='\t', index=False)
