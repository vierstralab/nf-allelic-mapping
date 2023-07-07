import sys
import pandas as pd
import os


def get_snps_file_name(value, prefix):
    if pd.isna(value):
        return value
    return os.path.join(prefix, str(value) + '.snps.bed')

def main(old_meta, output, file_prefix):
    df = pd.read_table(old_meta)
    df = df[~pd.isna(df['indiv_id'])]
    df['snps_file'] = df['ag_id'].apply(lambda x: get_snps_file_name(x, file_prefix))
    df.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    main(*sys.argv[1:])

