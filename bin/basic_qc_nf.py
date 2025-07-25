#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import sys
import argparse
import logging

from common_constants import *

csv.field_size_limit(sys.maxsize)


logger = logging.getLogger("basic_qc")
logger.setLevel(logging.INFO)
formatter = logging.Formatter("%(message)s")

# 1) must have SNP, PVAL, CHR, BP
# 2) *coerce headers
# 3) *if variant_id is None: try to set it to hm_variant_id or set variant_id to hm_variant_id?
# ---- The consequence of that is that we are both filling in and also removing if the rsids don't map any more
# ---- Do we want to keep rsids even if they don't map - and certainly don't map to the coordinates
# ---- OR if hm_rsid is not None, set variant_id = hm_rsid, thus updating all rsids where possible and keeping
# ---- ones that where not harmonised.
# 3) remove rows with blank values for (1)
# 4) conform to data types:
#   - if pval not floats: remove row
#   - if chr and bp not ints: remove row
# 5) set chr 'x' and 'y' to 23 and 24

REQUIRED_COLUMNS = [
    CHR_DSET,
    BP_DSET,
    EFFECT_DSET,
    OTHER_DSET,
    PVAL_DSET,
]

# At least one of these effect size columns must be present
EFFECT_SIZE_COLUMNS = [OR_DSET, ZSCORE_DSET, BETA_DSET]

BLANK_SET = {"", " ", "-", ".", "na", None, "none", "nan", "nil"}

# hm codes to drop
HM_CODE_FILTER = {9, 14, 15, 16, 17, 18}


def remove_row_if_required_is_blank(row, header):
    # 1) Any REQUIRED_COLUMN blank → drop
    for col in REQUIRED_COLUMNS:
        val = row[header.index(col)]
        if val == "NA" or (isinstance(val, str) and val.lower() in BLANK_SET):
            return True

    # 2) At least one EFFECT_SIZE_COLUMN must be non-blank
    effect_found = False
    for col in EFFECT_SIZE_COLUMNS:
        if col not in header:
            continue
        val = row[header.index(col)]
        if val != "NA" and not (isinstance(val, str) and val.lower() in BLANK_SET):
            effect_found = True
            break

    # if we found no valid effect size → drop
    if not effect_found:
        return True

    # otherwise keep
    return False


def remove_row_if_unharmonisable(row, header):
    try:
        hm_val = int(row[header.index(HM_CODE)])
    except (ValueError, TypeError):
        return False
    return hm_val in HM_CODE_FILTER


def blanks_to_NA(row):
    for n, i in enumerate(row):
        if i.lower() in BLANK_SET:
            row[n] = "NA"
    return row


def remove_row_if_wrong_data_type(row, header, col, data_type):
    try:
        data_type(row[header.index(col)])
        return False
    except ValueError:
        return True


def map_chr_values_to_numbers(row, header):
    index_chr = header.index(CHR_DSET)
    chromosome = row[index_chr].lower()
    if "x" in chromosome:
        chromosome = "23"
    if "y" in chromosome:
        chromosome = "24"
    if "mt" in chromosome:
        chromosome = "25"
    row[index_chr] = chromosome
    return row


def drop_last_element_from_filename(filename):
    filename_parts = filename.split("-")
    return "-".join(filename_parts[:-1])


def get_csv_reader(csv_file):
    dialect = csv.Sniffer().sniff(csv_file.readline())
    csv_file.seek(0)
    return csv.reader(csv_file, dialect)


def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "-f", help="The name of the file to be processed", required=True
    )
    argparser.add_argument("-o", help="The name of the output file", required=True)
    argparser.add_argument(
        "--print_only",
        help="only print the lines removed and do not write a new file",
        action="store_true",
    )
    argparser.add_argument("--log", help="The name of the log file")
    args = argparser.parse_args()

    file = args.f
    filename = args.o
    log_file = args.log

    new_filename = filename

    header = None
    is_header = True

    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    with open(file) as csv_file:
        writer = None
        if not args.print_only:
            result_file = open(new_filename, "w")
            writer = csv.writer(result_file, delimiter="\t")
        csv_reader = get_csv_reader(csv_file)

        for index, row in enumerate(csv_reader):
            if is_header:
                header = row
                # Sanity check: ensure all required columns and at least one effect size column are present
                missing_required = [
                    col for col in REQUIRED_COLUMNS if col not in header
                ]
                if missing_required:
                    raise ValueError(f"Missing required columns: {missing_required}")
                print(header)
                present_effect_size_columns = [
                    col for col in EFFECT_SIZE_COLUMNS if col in header
                ]
                if len(present_effect_size_columns) == 0:
                    raise ValueError(
                        f"At least one effect size column must be present: {EFFECT_SIZE_COLUMNS}"
                    )

                is_header = False
                if not args.print_only:
                    writer.writerows([header])
            else:
                # First try to replace an invalid variant_id with the hm_rsid
                # Checks for blanks, integers and floats:
                row = blanks_to_NA(row)
                row = map_chr_values_to_numbers(row, header)
                unharmonisable = remove_row_if_unharmonisable(row, header)
                blank = remove_row_if_required_is_blank(row, header)
                wrong_type_chr = remove_row_if_wrong_data_type(
                    row, header, CHR_DSET, int
                )
                wrong_type_bp = remove_row_if_wrong_data_type(row, header, BP_DSET, int)
                wrong_type_pval = remove_row_if_wrong_data_type(
                    row, header, PVAL_DSET, float
                )
                remove_row_tests = [
                    unharmonisable == False,
                    blank == False,
                    wrong_type_chr == False,
                    wrong_type_bp == False,
                    wrong_type_pval == False,
                ]
                if not args.print_only:
                    writer.writerows([row])

                if not all(remove_row_tests):
                    hm_code = int(row[header.index(HM_CODE)])
                    if hm_code not in HM_CODE_FILTER:
                        hm_code = 19
                    logger.info(
                        f"Removing record number {index}, with hm_code {hm_code}"
                    )


if __name__ == "__main__":
    main()
