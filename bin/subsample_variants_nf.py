#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import numpy as np
import shutil
import gzip
from common_constants import PVAL_DSET


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", required=True, help="Input summary statistics file")
    parser.add_argument("-o", required=True, help="Output subsampled file")
    parser.add_argument(
        "-l", type=int, default=200000, help="Maximum number of variants to keep"
    )
    return parser.parse_args()


def detect_delim(header_line: str) -> str:
    """Return the delimiter used in the file."""
    if "\t" in header_line:
        return "\t"
    if "," in header_line:
        return ","
    return None  # whitespace


def read_pval(line: str, idx: int, delim: str) -> float:
    parts = line.split(delim) if delim else line.split()
    return float(parts[idx])


def open_file(file_path: str, mode: str = "rt"):
    """Open a file, handling gzipped files automatically."""
    if file_path.endswith(".gz"):
        return gzip.open(file_path, mode)
    return open(file_path, mode)


def main():
    args = parse_args()
    limit = args.l

    # read in file line by line to not overload memory...
    with open_file(args.f) as f:
        header_line = next(f)

    delim = detect_delim(header_line)
    header_cols = header_line.strip().split(delim) if delim else header_line.split()
    # get the index of the p-value column
    pval_idx = header_cols.index(PVAL_DSET)
    print(pval_idx)
    # first pass: count total rows and collect significant indices
    sig_idx = []
    total = 0
    with open_file(args.f) as f:
        next(f)  # skip header
        for i, line in enumerate(f):
            total += 1
            try:
                pval = read_pval(line, pval_idx, delim)
            except (IndexError, ValueError):
                continue
            if pval < 1e-5:
                sig_idx.append(i)
    print(total)
    print(len(sig_idx))
    # if the total number of variants is less than or equal to the limit, copy the file directly
    if total <= limit:
        with open_file(args.f, "rb") as fin, open(args.o, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        return

    if len(sig_idx) >= limit:
        keep_rows = set(sig_idx)
        with open_file(args.f) as fin, open(args.o, "w") as fout:
            fout.write(next(fin))  # header
            for i, line in enumerate(fin):
                if i in keep_rows:
                    fout.write(line)
        return

    # calculate the number of random rows to sample
    n_random_sample = limit - len(sig_idx)

    # randomly sample indices
    rng = np.random.default_rng(42)
    random_idx = rng.choice(total, n_random_sample, replace=False)

    # combine significant indices and random indices
    keep_rows = set(sig_idx) | set(random_idx)

    # write the output file
    with open_file(args.f) as fin, open(args.o, "w") as fout:
        fout.write(next(fin))  # header
        for i, line in enumerate(fin):
            if i in keep_rows:
                fout.write(line)


if __name__ == "__main__":
    main()
