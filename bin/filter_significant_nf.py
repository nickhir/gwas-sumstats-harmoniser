#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
from common_constants import PVAL_DSET


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", required=True, help="Input merged file")
    parser.add_argument(
        "-o", required=True, help="Output file with significant variants"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    pvals = pd.read_table(args.f, usecols=[PVAL_DSET], sep=None, engine="python")
    sig_idx = pvals[pvals[PVAL_DSET] < 1e-5].index.tolist()

    with open(args.f) as fin, open(args.o, "w") as fout:
        fout.write(next(fin))  # header
        for i, line in enumerate(fin):
            if i in sig_idx:
                fout.write(line)


if __name__ == "__main__":
    main()
