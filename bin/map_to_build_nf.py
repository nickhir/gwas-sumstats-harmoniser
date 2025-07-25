#!/usr/bin/env python
# -*- coding: utf-8 -*-

# merge on hm_variant_id with vcf of desired build
# if variant_id is rsid and ID != variant_id or not synonym --> drop and create discrep df


import pandas as pd
import liftover as lft


# in common_constants there are a lot of variables defined...
import sys
from common_constants import *

import os
import glob
import argparse
from ast import literal_eval

import pyarrow.parquet as pq
import re

# Allow very large fields in input file-------------
import csv

maxInt = sys.maxsize

while True:
    # decrease the maxInt value by factor 10
    # as long as the OverflowError occurs.

    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt / 10)


# map_to_build----------------------------------------------------
def merge_ss_vcf(ss, vcf, from_build, to_build, chroms, coordinate):
    vcfs = glob.glob(vcf)

    ssdf = pd.read_table(ss, sep=None, engine="python", dtype=str)
    add_fields_if_missing(df=ssdf)

    # If the input is already on the desired build, skip any mapping/liftover
    if from_build == to_build:
        print(
            "source and target builds are identical - skipping RSID mapping and liftover"
        )
        ssdf[HM_CC_DSET] = "NA"
        for chrom in chroms:
            df = ssdf.loc[ssdf[CHR_DSET].astype(str) == chrom]
            df = df.dropna(subset=[BP_DSET])
            df[BP_DSET] = df[BP_DSET].astype(str).str.replace(r"\..*$", "", regex=True)
            outfile = os.path.join(f"{chrom}.merged")
            df.to_csv(outfile, sep="\t", mode="w", index=False, na_rep="NA")

        no_chr_df = ssdf[ssdf[CHR_DSET].isnull()]
        no_pos_df = ssdf[ssdf[BP_DSET].isnull()]
        no_loc_df = pd.concat([no_chr_df, no_pos_df])
        no_loc_df[HM_CC_DSET] = "NA"
        outfile = os.path.join("unmapped")
        no_loc_df.to_csv(outfile, sep="\t", index=False, na_rep="NA")
        return

    # creates a list of True False depending on if there is an rsid
    rsid_mask = ssdf[RSID].str.startswith("rs").fillna(False)
    ssdf_with_rsid = ssdf[rsid_mask]
    ssdf_without_rsid = ssdf[~rsid_mask]
    header = list(ssdf.columns.values)
    print("starting rsid mapping")
    print("ssdf with rsid empty?: {}".format(ssdf_with_rsid.empty))
    # if there are records with rsids

    # keep track of rsids for which mapping worked so we can remove them later.
    mapped_rsids = set()
    if not ssdf_with_rsid.empty:
        for vcf in vcfs:
            # from the vcf file name, extract the chromosome we are currently processing.
            pattern = re.compile(r"chr([0-9]+|X|Y|MT)")
            match = pattern.search(vcf)
            if not match:
                sys.exit(f"Error: could not extract chromosome from path: {vcf}")
            file_chr = match.group(1)
            if not file_chr:
                sys.exit(f"Error: could not extract chromosome from path: {vcf}")

            # if the chromosome is not in the list of chromosomes to process, skip it.
            if file_chr not in chroms:
                print(
                    "skipping chromosome {} as it is not in the list of chromosomes to process".format(
                        file_chr
                    )
                )
                continue

            # subset the ssdf_with_rsid to the current chromosome
            # technically, in rare cases where not just the position, but also the chromosome changes.
            # this kind of mapping will fail. but they will then be liftovered with the "normal" tool anyway.
            ssdf_with_rsid_subset = ssdf_with_rsid[ssdf_with_rsid[CHR_DSET] == file_chr]

            # connect to the parquet file
            vcf_pf = pq.ParquetFile(vcf)
            # now read the parquet file into memory in chunks of 1M rows -> doesnt overwhelm memory
            for batch in vcf_pf.iter_batches(batch_size=1_000_000):
                vcf_df = batch.to_pandas()
                chrom = vcf_df.CHR.iloc[0]
                # merge on rsid, update chr and position from vcf ref
                mergedf = ssdf_with_rsid_subset.merge(
                    vcf_df, left_on=RSID, right_on="ID", how="left"
                )
                mapped = mergedf.dropna(subset=["ID"]).drop([CHR_DSET, BP_DSET], axis=1)

                # Track the successfully mapped rsIDs
                if not mapped.empty:
                    mapped_rsids.update(mapped[RSID].tolist())

                mapped[CHR_DSET] = (
                    mapped["CHR"].astype("str").str.replace(r"\..*$", "", regex=True)
                )
                mapped[BP_DSET] = (
                    mapped["POS"].astype("str").str.replace(r"\..*$", "", regex=True)
                )
                mapped = mapped[header]
                mapped[HM_CC_DSET] = "rs"
                outfile = os.path.join("{}.merged".format(chrom))
                # if the file does not exist, we will write the header, otherwise just the content
                write_header = not os.path.exists(outfile)
                mapped.to_csv(
                    outfile,
                    sep="\t",
                    index=False,
                    na_rep="NA",
                    mode="a",
                    header=write_header,
                )
                # remove them successfully mapped rsids from the ssdf_with_rsid_subset to speed up merge step
                ssdf_with_rsid_subset = mergedf.loc[mergedf["ID"].isnull(), header]

    # Remove all successfully mapped rsIDs from the original dataframe
    num_mapped = len(mapped_rsids)
    print(f"finished rsid mapping, successfully mapped {num_mapped} variants")
    ssdf_with_rsid = ssdf_with_rsid[~ssdf_with_rsid[RSID].isin(mapped_rsids)]

    # liftover the snps without rsids and those with unrecognised rsids
    print("liftover remaining variants")
    ssdf = pd.concat([ssdf_with_rsid, ssdf_without_rsid])
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
    # the liftover resulst are appended if rsid mapping worked.
    for chrom in chroms:
        df = ssdf.loc[ssdf[CHR_DSET].astype("str") == chrom]
        df = df.dropna(subset=[BP_DSET])
        df[BP_DSET] = df[BP_DSET].astype("str").str.replace(r"\..*$", "")
        outfile = os.path.join("{}.merged".format(chrom))
        if os.path.isfile(outfile):
            print("df to {}".format(outfile))
            df.to_csv(
                outfile, sep="\t", mode="a", header=False, index=False, na_rep="NA"
            )
        else:
            print("df to {}".format(outfile))
            df.to_csv(outfile, sep="\t", mode="w", index=False, na_rep="NA")
    print("liftover complete")
    # files for which liftover failed completely are also recorded.
    no_chr_df = ssdf[ssdf[CHR_DSET].isnull()]
    no_pos_df = ssdf[ssdf[BP_DSET].isnull()]
    no_loc_df = pd.concat([no_chr_df, no_pos_df])
    no_loc_df[HM_CC_DSET] = "NA"
    outfile = os.path.join("unmapped")
    no_loc_df.to_csv(outfile, sep="\t", index=False, na_rep="NA")


def listify_string(string):
    """
    listify the input. If it's a list leave it.
    It it looks like at list, make it a list.
    Otherwise convert the input into a list - which is probably not what's wanted.
    :param string:
    :return: a list
    """
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
    argparser.add_argument(
        "-f", help="The name of the file to be processed", required=True
    )
    argparser.add_argument("-vcf", help="The name of the vcf file", required=True)
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
    vcf = args.vcf
    from_build = args.from_build
    to_build = args.to_build
    chroms = listify_string(args.chroms)
    coordinate = args.coordinate

    merge_ss_vcf(ss, vcf, from_build, to_build, chroms, coordinate)


if __name__ == "__main__":
    main()
