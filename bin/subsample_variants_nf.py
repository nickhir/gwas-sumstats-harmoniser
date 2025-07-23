#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import gzip
import sys
from common_constants import (
    PVAL_DSET,
    NEG_LOG_PVAL_DSET,
    COLUMN_ALIASES,
    ALIAS_LOOKUP,
    CHR_DSET,
    BP_DSET,
    EFFECT_DSET,
    OTHER_DSET,
    OR_DSET,
    RSID,
    ZSCORE_DSET,
    BETA_DSET,
    SE_DSET,
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
]

# At least one of these effect size columns must be present
EFFECT_SIZE_COLUMNS = [OR_DSET, ZSCORE_DSET, BETA_DSET]

DISPLAY_NAMES = {
    CHR_DSET: "Chromosome",
    BP_DSET: "BP",
    EFFECT_DSET: "Effect allele",
    OTHER_DSET: "Other allele",
    OR_DSET: "Odds ratio",
    ZSCORE_DSET: "Zscore",
    BETA_DSET: "Beta",
    SE_DSET: "Standard error",
    PVAL_DSET: "P value",
    NEG_LOG_PVAL_DSET: "Neg log10 P value",
}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", required=True, help="Input summary statistics file")
    parser.add_argument("-o", required=True, help="Output subsampled file")
    parser.add_argument(
        "-l", type=int, default=200000, help="Maximum number of variants to keep"
    )
    parser.add_argument(
        "--alias_log",
        help="Optional file to write detected header aliases",
    )
    return parser.parse_args()


def detect_delim(header_line):
    """Return the delimiter used in the file."""
    if "\t" in header_line:
        return "\t"
    if "," in header_line:
        return ","
    if ";" in header_line:
        return ";"
    if "|" in header_line:
        return "|"
    return None  # whitespace


def read_pval(line: str, idx: int, delim: str) -> float:
    parts = line.split(delim) if delim else line.split()
    return float(parts[idx])


def open_file(file_path: str, mode: str = "rt"):
    """Open a file, handling gzipped files automatically."""
    if file_path.endswith(".gz"):
        return gzip.open(file_path, mode)
    return open(file_path, mode)


def format_line(line: str, delim: str, indices):
    """Return tab-separated line with selected columns."""
    parts = line.rstrip("\n").split(delim) if delim else line.split()
    selected = [parts[i] if i < len(parts) else "" for i in indices]
    return "\t".join(selected) + "\n"


# go through each column name and look it up in the alias_lookup.
# if its found, then we replace it with the "first" choice name (for example p_value).
# IF NO MATCH IS FOUND THE COLUMN WILL BE REMOVED!!
def canonicalise_header(header):
    """Return canonicalised header, mapping and indices of recognised columns."""
    canonical_header = []
    mapping = {}
    keep_idx = []
    for i, col in enumerate(header):
        canonical = ALIAS_LOOKUP.get(col.lower())
        if canonical is None:
            # skip unrecognised columns
            continue
        canonical_header.append(canonical)
        keep_idx.append(i)
        # record the first appearance of this canonical name
        mapping.setdefault(canonical, col)
    return canonical_header, mapping, keep_idx


def resolve_rsid_duplicate(file_path, delim, indices):
    """Return the index that contains rsIDs from a list of indices."""
    with open_file(file_path) as f:
        next(f)  # skip header
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split(delim) if delim else line.split()
            values = [parts[i] if i < len(parts) else "" for i in indices]
            matches = [i for i, v in enumerate(values) if v.lower().startswith("rs")]
            if len(matches) == 1:
                return indices[matches[0]]
    # default to the first if no distinguishing row found
    return indices[0]


# all it does is go through potential names of the "pval" column
# and if the potential name matches any name in the header, return that index which is used for the read_pval function
def detect_pval_col(header):
    """Return (column_name, index, is_neglog) for p-value information."""
    if PVAL_DSET in header:
        return PVAL_DSET, header.index(PVAL_DSET), False
    if NEG_LOG_PVAL_DSET in header:
        return NEG_LOG_PVAL_DSET, header.index(NEG_LOG_PVAL_DSET), True
    for candidate in COLUMN_ALIASES.get(PVAL_DSET, []):
        if candidate in header:
            return candidate, header.index(candidate), False
    raise ValueError("P-value column not found")


def check_required_columns(header):
    """Ensure all required columns are present after header normalisation."""
    missing = [c for c in REQUIRED_COLUMNS if c not in header]
    has_p = PVAL_DSET in header
    has_neglog = NEG_LOG_PVAL_DSET in header
    has_effect_size = any([c in header for c in EFFECT_SIZE_COLUMNS])
    if missing or not has_effect_size or not (has_p or has_neglog):
        messages = []
        if missing:
            messages.append(f"Missing required columns: {', '.join(missing)}.")
        if not (has_p or has_neglog):
            messages.append(
                f"Either {PVAL_DSET} or {NEG_LOG_PVAL_DSET} must be present."
            )
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
    canon_header, orig_map, keep_idx = canonicalise_header(header_cols)
    idx_map = {}
    for i, canon in enumerate(canon_header):
        # in the dictonary idx_map try and get the key for canon. if it doesn't exist, create an empty list 
        # and to the list append the index i
        # IF a canon header name occures twice, then the list will contain two indices!
        idx_map.setdefault(canon, []).append(i)
    for canon, idx_list in list(idx_map.items()):
        if len(idx_list) > 1:
            # only if the duplicated column is the RSID column, we resolve it, otherwise ERROR
            if canon == RSID:
                candidate_indices = [keep_idx[i] for i in idx_list]
                chosen = resolve_rsid_duplicate(args.f, delim, candidate_indices)
                keep_canon_idx = idx_list[candidate_indices.index(chosen)]
                for j in sorted(idx_list, reverse=True):
                    if j != keep_canon_idx:
                        del canon_header[j]
                        del keep_idx[j]
                orig_map[canon] = header_cols[chosen]
            else:
                cols = [header_cols[keep_idx[i]] for i in idx_list]
                sys.exit(f"Multiple columns mapped to {canon}: {', '.join(cols)}")

    if args.alias_log:
        with open(args.alias_log, "w") as fh:
            # determine which columns to report
            cols_to_report = list(REQUIRED_COLUMNS)
            cols_to_report += [c for c in EFFECT_SIZE_COLUMNS if c in canon_header]
            if NEG_LOG_PVAL_DSET in canon_header:
                cols_to_report.append(NEG_LOG_PVAL_DSET)
            elif PVAL_DSET in canon_header:
                cols_to_report.append(PVAL_DSET)

            for canon in cols_to_report:
                if canon not in orig_map:
                    continue
                orig = orig_map[canon]
                name = DISPLAY_NAMES.get(canon, canon)
                fh.write(
                    f"{name} column was identified as {orig} and renamed to {canon}\n"
                )

    # create the output header
    header_out = "\t".join(str(x) for x in canon_header) + "\n"
    # check if required columns were deteceted, otherwise exit
    check_required_columns(canon_header)
    pval_col, pval_idx, is_neglog = detect_pval_col(canon_header)
    pval_idx_orig = keep_idx[pval_idx]
    # first pass: count total rows and collect significant indices
    sig_idx = []
    total = 0
    with open_file(args.f) as f:
        next(f)  # skip header
        for i, line in enumerate(f):
            total += 1
            try:
                val = read_pval(line, pval_idx_orig, delim)
            except (IndexError, ValueError):
                continue
            if is_neglog:
                if val > 5:
                    sig_idx.append(i)
            else:
                if val < 1e-5:
                    sig_idx.append(i)

    # if the total number of variants is less than or equal to the limit, copy the file directly
    if total <= limit:
        with open_file(args.f) as fin, open(args.o, "w") as fout:
            next(fin)  # drop header
            fout.write(header_out)
            for line in fin:
                fout.write(format_line(line, delim, keep_idx))
        return

    # if more significant variants than the limit, write all significant variants
    if len(sig_idx) >= limit:
        keep_rows = set(sig_idx)
        with open_file(args.f) as fin, open(args.o, "w") as fout:
            next(fin)  # discard original header
            fout.write(header_out)
            for i, line in enumerate(fin):
                if i in keep_rows:
                    fout.write(format_line(line, delim, keep_idx))
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
                fout.write(format_line(line, delim, keep_idx))


if __name__ == "__main__":
    main()
