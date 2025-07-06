#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import gzip
import sys
from common_constants import (
    PVAL_DSET,
    COLUMN_ALIASES,
    ALIAS_LOOKUP,
    CHR_DSET,
    BP_DSET,
    EFFECT_DSET,
    OTHER_DSET,
    OR_DSET,
    ZSCORE_DSET,
    BETA_DSET,
)

# Columns that must be present for downstream steps
# We only need the p-value column for subsampling but downstream
# modules expect these canonical columns to exist in the output.
# Therefore ensure we can map them from any header aliases before
# proceeding.
REQUIRED_COLUMNS = [
    CHR_DSET,
    BP_DSET,
    EFFECT_DSET,
    OTHER_DSET,
    PVAL_DSET,
]

# At least one of these effect size columns must be present
EFFECT_SIZE_COLUMNS = [OR_DSET, ZSCORE_DSET, BETA_DSET]


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


# go through each column name and look it up in the alias_lookup.
# if its found, then we replace it with the "first" choice name (for example p_value).
# if not, we just keep the original one.
def canonicalise_header(header):
    """Replace known aliases with canonical column names."""
    return [ALIAS_LOOKUP.get(col.lower(), col) for col in header]


# all it does is go through potential names of the "pval" column
# and if the potential name matches any name in the header, return that index which is used for the read_pval function
def detect_pval_idx(header):
    """Return the index of the p-value column using aliases."""
    for candidate in [PVAL_DSET] + COLUMN_ALIASES.get(PVAL_DSET, []):
        if candidate in header:
            return header.index(candidate)
    raise ValueError("P-value column not found")


def check_required_columns(header):
    """Ensure all required columns are present after header normalisation."""
    missing = [c for c in REQUIRED_COLUMNS if c not in header]
    has_effect_size = any([c in header for c in EFFECT_SIZE_COLUMNS])
    if missing or not has_effect_size:
        messages = []
        if missing:
            messages.append(f"Missing required columns: {', '.join(missing)}.")
        if not has_effect_size:
            messages.append(
                "Missing effect size column: one of "
                f"{', '.join(EFFECT_SIZE_COLUMNS)} must be present."
            )
        messages.append(f"Available columns: {', '.join(header)}")
        sys.exit(" ".join(messages))


def main():
    args = parse_args()
    limit = args.l

    # read in file line by line to not overload memory...
    with open_file(args.f) as f:
        header_line = next(f)

    delim = detect_delim(header_line)
    header_cols = header_line.strip().split(delim) if delim else header_line.split()
    # get the index of the p-value column
    canon_header = canonicalise_header(header_cols)
    # create the output header
    header_out = (delim or " ").join(canon_header) + "\n"
    # check if required columns were deteceted, otherwise exit
    check_required_columns(canon_header)
    pval_idx = canon_header.index(PVAL_DSET)
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
    # if the total number of variants is less than or equal to the limit, copy the file directly
    if total <= limit:
        with open_file(args.f) as fin, open(args.o, "w") as fout:
            next(fin)  # drop header
            fout.write(header_out)
            for line in fin:
                fout.write(line)
        return

    if len(sig_idx) >= limit:
        keep_rows = set(sig_idx)
        with open_file(args.f) as fin, open(args.o, "w") as fout:
            next(fin)  # discard original header
            fout.write(header_out)
            for i, line in enumerate(fin):
                if i in keep_rows:
                    fout.write(line)
        return

    # calculate the number of random rows to sample
    n_random_sample = limit - len(sig_idx)

    # randomly sample indices
    rng = np.random.default_rng(42)
    random_idx = rng.choice(total, n_random_sample, replace=False)

    # combine significant indices and random indices.
    # now it can happen that we randomly sample significant rows.
    # but this shouldnt happen to often and doesnt really matter.
    # i.e. its no essential to have exactly "limit" number of rows in the output
    keep_rows = set(sig_idx) | set(random_idx)

    # write the output file
    with open_file(args.f) as fin, open(args.o, "w") as fout:
        next(fin)  # discard header
        fout.write(header_out)
        for i, line in enumerate(fin):
            if i in keep_rows:
                fout.write(line)


if __name__ == "__main__":
    main()
