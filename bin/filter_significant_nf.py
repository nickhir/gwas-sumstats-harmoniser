#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
from common_constants import PVAL_DSET, NEG_LOG_PVAL_DSET


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", required=True, help="Input merged file")
    parser.add_argument(
        "-o", required=True, help="Output file with significant variants"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Determine which p-value column is present
    header = pd.read_table(args.f, nrows=0, sep=None, engine="python").columns
    if PVAL_DSET in header:
        pval_col = PVAL_DSET
        is_neglog = False
    elif NEG_LOG_PVAL_DSET in header:
        pval_col = NEG_LOG_PVAL_DSET
        is_neglog = True
    else:
        raise ValueError(
            f"Input must contain either {PVAL_DSET} or {NEG_LOG_PVAL_DSET}"
        )

    pvals = pd.read_table(args.f, usecols=[pval_col], sep=None, engine="python")
    if is_neglog:
        sig_idx = pvals[pvals[pval_col] > 5].index.tolist()
    else:
        sig_idx = pvals[pvals[pval_col] < 1e-5].index.tolist()

    with open(args.f) as fin, open(args.o, "w") as fout:
        fout.write(next(fin))  # header
        for i, line in enumerate(fin):
            if i in sig_idx:
                fout.write(line)


if __name__ == "__main__":
    main()
