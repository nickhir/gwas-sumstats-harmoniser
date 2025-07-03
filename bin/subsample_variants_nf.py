#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import shutil
from common_constants import PVAL_DSET


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', required=True, help='Input summary statistics file')
    parser.add_argument('-o', required=True, help='Output subsampled file')
    parser.add_argument('-l', type=int, default=10000000, help='Maximum number of variants to keep')
    return parser.parse_args()


def main():
    args = parse_args()
    limit = args.l

    # read only p-value column to get counts
    pvals = pd.read_table(args.f, usecols=[PVAL_DSET], sep=None, engine='python')
    total = len(pvals)
    if total <= limit:
        shutil.copyfile(args.f, args.o)
        return

    sig_mask = pvals[PVAL_DSET] < 1e-5
    n_sig = int(sig_mask.sum())

    df = pd.read_table(args.f, sep=None, engine='python', dtype=str)
    sig_df = df[sig_mask]
    nonsig_df = df[~sig_mask]

    n_non_sig_sample = max(0, limit - len(sig_df))
    if len(nonsig_df) > n_non_sig_sample:
        nonsig_df = nonsig_df.sample(n=n_non_sig_sample, random_state=42)

    out_df = pd.concat([sig_df, nonsig_df])
    out_df.to_csv(args.o, sep='\t', index=False, na_rep='NA')


if __name__ == '__main__':
    main()
