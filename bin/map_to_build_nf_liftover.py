#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Liftover-only variant location mapping.

This script is an alternative to ``map_to_build_nf.py`` but skips the
RSID-based matching step entirely. All input variants are mapped via
liftover regardless of whether an RSID is present.
"""

import pandas as pd
import liftover as lft

# in common_constants there are a lot of variables defined...
import sys
from common_constants import *

import os
import argparse
from ast import literal_eval

# Allow very large fields in input file-------------
import csv

maxInt = sys.maxsize
while True:
    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt / 10)


def liftover_ss(ss, from_build, to_build, chroms, coordinate):
    """Liftover variants in ``ss`` to ``to_build``."""
    ssdf = pd.read_table(ss, sep=None, engine="python", dtype=str)
    add_fields_if_missing(df=ssdf)
    ssdf[HM_CC_DSET] = "lo"

    build_map = (
        lft.LiftOver(lft.ucsc_release.get(from_build), lft.ucsc_release.get(to_build))
        if from_build != to_build
        else None
    )
    if build_map:
        ssdf[BP_DSET] = [
            lft.map_bp_to_build_via_liftover(
                chromosome=x, bp=y, build_map=build_map, coordinate=coordinate[0]
            )
            for x, y in zip(ssdf[CHR_DSET], ssdf[BP_DSET])
        ]

    for chrom in chroms:
        df = ssdf.loc[ssdf[CHR_DSET].astype("str") == chrom]
        df = df.dropna(subset=[BP_DSET])
        df[BP_DSET] = df[BP_DSET].astype("str").str.replace("\..*$", "")
        outfile = os.path.join(f"{chrom}.merged")
        if os.path.isfile(outfile):
            df.to_csv(
                outfile,
                sep="\t",
                mode="a",
                header=False,
                index=False,
                na_rep="NA",
            )
        else:
            df.to_csv(
                outfile,
                sep="\t",
                mode="w",
                index=False,
                na_rep="NA",
            )

    print("liftover complete")

    no_chr_df = ssdf[ssdf[CHR_DSET].isnull()]
    no_pos_df = ssdf[ssdf[BP_DSET].isnull()]
    no_loc_df = pd.concat([no_chr_df, no_pos_df])
    no_loc_df[HM_CC_DSET] = "NA"
    outfile = os.path.join("unmapped")
    no_loc_df.to_csv(outfile, sep="\t", index=False, na_rep="NA")


def listify_string(string):
    if type(string) is str:
        if "[" and "]" in string:
            new = (
                string.replace(" ", "")
                .replace("[", '["')
                .replace("]", '"]')
                .replace(",", '","')
            )
            listified = literal_eval(new)
        else:
            listified = list(string)
    elif type(string) is list:
        listified = string
    else:
        listified = list(str(string))
    return listified


def add_fields_if_missing(df):
    add_column_to_df(df=df, column=RSID)
    add_column_to_df(df=df, column=CHR_DSET)
    add_column_to_df(df=df, column=BP_DSET)


def add_column_to_df(df, column, value="NA"):
    if column not in df.columns:
        df[column] = value


def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-f", help="The name of the file to be processed", required=True)
    argparser.add_argument("-vcf", help="(unused) The name of the vcf file")
    argparser.add_argument("--log", help="The name of the log file")
    argparser.add_argument(
        "-from_build",
        help='The original build e.g. "36" for NCBI36 or hg18',
        required=True,
    )
    argparser.add_argument(
        "-to_build", help='The latest (desired) build e.g. "38"', required=True
    )
    argparser.add_argument(
        "-chroms", help="A list of chromosomes to process", default=DEFAULT_CHROMS
    )
    argparser.add_argument(
        "-coordinate", help="index", nargs="?", const="1-based", required=True
    )
    args = argparser.parse_args()

    ss = args.f
    from_build = args.from_build
    to_build = args.to_build
    chroms = listify_string(args.chroms)
    coordinate = args.coordinate

    liftover_ss(ss, from_build, to_build, chroms, coordinate)


if __name__ == "__main__":
    main()