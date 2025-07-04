#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import numpy as np
import shutil
from common_constants import PVAL_DSET


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", required=True, help="Input summary statistics file")
    parser.add_argument("-o", required=True, help="Output subsampled file")
    parser.add_argument(
        "-l", type=int, default=10000000, help="Maximum number of variants to keep"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    limit = args.l

    # read only p-value column to get counts
    pvals = pd.read_table(args.f, usecols=[PVAL_DSET], sep=None, engine="python")
    total = len(pvals)
    # if the total number of variants is less than or equal to the limit dont do anything.
    if total <= limit:
        shutil.copyfile(args.f, args.o)
        return

    sig_mask = pvals[PVAL_DSET] < 1e-5

    # turns the bool list into a list of indices (row positions)
    sig_idx = sig_mask[sig_mask].index.tolist()
    nonsig_idx = sig_mask[~sig_mask].index.tolist()

    # how many random rows do we have to sample to get to the limit?
    n_non_sig_sample = max(0, limit - len(sig_idx))
    rng = np.random.default_rng(42)
    if len(nonsig_idx) > n_non_sig_sample:
        sample_idx = rng.choice(nonsig_idx, n_non_sig_sample, replace=False)
    else:
        sample_idx = nonsig_idx

    keep_rows = set(sig_idx) | set(sample_idx)

    with open(args.f) as fin, open(args.o, "w") as fout:
        fout.write(next(fin))  # header
        for i, line in enumerate(fin):
            if i in keep_rows:
                fout.write(line)


if __name__ == "__main__":
    main()
