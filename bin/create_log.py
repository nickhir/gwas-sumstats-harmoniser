#!/usr/bin/env python
# coding: utf-8

import argparse
import os
import subprocess
import sys
from datetime import datetime

CODE_TABLE = {
    1: "Palindromic; Infer strand; Forward strand; Correct orientation; Already harmonised",
    2: "Palindromic; Infer strand; Forward strand; Flipped orientation; Requires harmonisation",
    3: "Palindromic; Infer strand; Reverse strand; Correct orientation; Already harmonised",
    4: "Palindromic; Infer strand; Reverse strand; Flipped orientation; Requires harmonisation",
    5: "Palindromic; Assume forward strand; Correct orientation; Already harmonised",
    6: "Palindromic; Assume forward strand; Flipped orientation; Requires harmonisation",
    7: "Palindromic; Assume reverse strand; Correct orientation; Already harmonised",
    8: "Palindromic; Assume reverse strand; Flipped orientation; Requires harmonisation",
    9: "Palindromic; Drop palindromic; Will not harmonise",
    10: "Forward strand; Correct orientation; Already harmonised",
    11: "Forward strand; Flipped orientation; Requires harmonisation",
    12: "Reverse strand; Correct orientation; Already harmonised",
    13: "Reverse strand; Flipped orientation; Requires harmonisation",
    14: "Required fields are not known; Cannot harmonise",
    15: "No matching variants in reference VCF; Cannot harmonise",
    16: "Multiple matching variants in reference VCF (ambiguous); Cannot harmonise",
    17: "Palindromic; Infer strand; EAF or reference VCF AF not known; Cannot harmonise",
    18: "Palindromic; Infer strand; EAF < --maf_palin_threshold; Will not harmonise",
    19: "QC failure; Any of rsID mismatch with reference, data type mismatch, missing data",
}

HM_CODE_FILTER = {9, 14, 15, 16, 17, 18, 19}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate harmonisation running log", add_help=False
    )
    parser.add_argument("-r", "--reference", required=True, help="Reference VCF file")
    parser.add_argument("-i", "--input", required=True, help="Input raw data file")
    parser.add_argument("-c", "--count", required=True, help="Total strand count file")
    parser.add_argument(
        "-d", "--deleted", required=True, help="Deleted sites file from QC"
    )
    parser.add_argument(
        "-h", "--harmonized", required=True, help="Harmonised result file"
    )
    parser.add_argument("-u", "--unmapped", required=True, help="Unmapped sites file")
    parser.add_argument("-o", "--output", required=True, help="Output log file")
    parser.add_argument(
        "-p", "--pipeline_version", required=True, help="Pipeline version"
    )
    parser.add_argument("-a", "--alias", help="Header alias log file")
    parser.add_argument("--help", action="help", help="Show this help message and exit")
    return parser.parse_args()


def count_records(path):
    with open(path) as f:
        return max(sum(1 for _ in f) - 1, 0)


def parse_count_file(path):
    palin_mode = None
    ratio_line = None
    ratio_number = None
    with open(path) as f:
        for line in f:
            if "palin_mode" in line:
                parts = line.strip().split("\t")
                if len(parts) > 1:
                    palin_mode = parts[1]
            if "ratio" in line:
                ratio_line = line.strip()
                parts = line.strip().split()
                if len(parts) > 1:
                    try:
                        ratio_number = float(parts[1])
                    except ValueError:
                        pass
    return palin_mode, ratio_number, ratio_line


def read_hm_codes(stdin):
    codes = []
    for line in stdin:
        line = line.strip()
        if not line:
            continue
        try:
            codes.append(int(line))
        except ValueError:
            pass
    return codes


def summarise_hm_codes(codes):
    counts = {}
    for code in codes:
        counts[code] = counts.get(code, 0) + 1

    with open("report.txt") as f:
        fails = f.readlines()

    total = sum(counts.values()) + len(fails)
    lines = []
    lines.append("\n6. Successfully harmonised variants\n")
    success_codes = [c for c in counts if c < 14 and c != 9]
    success_total = sum(counts[c] for c in success_codes)
    success_ratio = success_total / total if total else 0
    lines.append(
        f"{success_ratio:.2%} ({success_total} of {total}) sites successfully harmonised.\n"
    )
    lines.append("hm_code\tNumber\tPercentage\tExplanation\n")
    for c in sorted(success_codes):
        count = counts[c]
        per = count / total if total else 0
        lines.append(f"{c}\t{count}\t{per:.2%}\t{CODE_TABLE[c]}\n")
    lines.append("\n################################################################\n")

    lines.append("\n7. Failed harmonisation\n")
    hm_code_fail_dict = {}
    for line in fails:
        if "hm_code" in line:
            try:
                hm_code = int(line.split("hm_code")[-1].strip())
            except ValueError:
                continue
            if hm_code in HM_CODE_FILTER:
                hm_code_fail_dict[hm_code] = hm_code_fail_dict.get(hm_code, 0) + 1
    fail_total = sum(hm_code_fail_dict.values())
    fail_ratio = fail_total / total if total else 0
    lines.append(
        f"{fail_ratio:.2%} ({fail_total} of {total}) sites failed to harmonise.\n"
    )
    lines.append("hm_code\tNumber\tPercentage\tExplanation\n")
    for c in sorted(hm_code_fail_dict):
        count = hm_code_fail_dict[c]
        per = count / total if total else 0
        lines.append(f"{c}\t{count}\t{per:.2%}\t{CODE_TABLE[c]}\n")
    lines.append("\n################################################################\n")

    lines.append("\n7. Overview\n")
    lines.append(
        "Result\tFAILED_HARMONIZATION\n"
        if fail_ratio == 1.0
        else "Result\tSUCCESS_HARMONIZATION\n"
    )
    return "".join(lines)


def main():
    args = parse_args()
    hm_codes = read_hm_codes(sys.stdin)
    palin_mode, ratio_num, ratio_line = parse_count_file(args.count)

    unmapped_sites = count_records(args.unmapped)
    mapped_sites = count_records(args.harmonized)
    total_sites = unmapped_sites + mapped_sites
    unmapped_rate = unmapped_sites / total_sites * 100 if total_sites else 0
    mapped_rate = mapped_sites / total_sites * 100 if total_sites else 0

    ref_source = ref_reference = ref_dbsnp = ""
    try:
        header = subprocess.check_output(["tabix", "-H", args.reference], text=True)
        for line in header.splitlines():
            if "source" in line and not ref_source:
                ref_source = line
            if "reference" in line and not ref_reference:
                ref_reference = line
            if "dbSNP" in line and not ref_dbsnp:
                ref_dbsnp = line.replace("INFO=<", "").replace(">", "")
    except Exception:
        pass

    out_lines = []
    out_lines.append(
        "################################################################\n"
    )
    out_lines.append("HARMONISATION RUNNING REPORT\n")
    out_lines.append(
        "################################################################\n"
    )
    out_lines.append("1. Pipeline details\n")
    out_lines.append(f"    A. Pipeline Version: {args.pipeline_version}\n")
    out_lines.append(f"    B. Running date: {datetime.now().strftime('%b %d %Y')}\n")
    out_lines.append(f"    C. Input file: {os.path.basename(args.input)}\n")
    out_lines.append(
        "################################################################\n\n"
    )

    out_lines.append("2. Reference data\n")
    if ref_source:
        out_lines.append(f"{ref_source}\n")
    if ref_reference:
        out_lines.append(f"{ref_reference}\n")
    if ref_dbsnp:
        out_lines.append(f"{ref_dbsnp}\n")
    out_lines.append(
        "################################################################\n\n"
    )

    if args.alias and os.path.isfile(args.alias):
        with open(args.alias) as a:
            out_lines.append("3. Header column mapping\n")
            out_lines.append(a.read())
            out_lines.append("\n")
    out_lines.append(
        "################################################################\n\n"
    )

    out_lines.append("4. Mapping result\n")
    out_lines.append(
        f"{unmapped_rate:.2f}% ({unmapped_sites} sites out of {total_sites}) were dropped because they could not be mapped.\n{mapped_rate:.2f}% ({mapped_sites} sites) were carried forward.\n"
    )
    out_lines.append(
        "\n################################################################\n"
    )
    out_lines.append("5. Palindromic SNPs\n")
    out_lines.append(f"palin_mode: {palin_mode}\n")
    if palin_mode == "drop":
        if ratio_num is None:
            out_lines.append(
                "Palindromic SNPs could not be harmonized because the direction of palindromic SNPs cannot be inferred from consensus direction.\n"
            )
        else:
            out_lines.append(
                f"Palindromic SNPs could not be harmonized because the direction of palindromic SNPs cannot be inferred from consensus direction (forward sites ratio ={ratio_num}).\n"
            )
    else:
        if ratio_line and "Full" in ratio_line:
            out_lines.append(
                f"Direction of palindromic SNPs inferred as {palin_mode} by establishing consensus direction of all sites (forward sites ratio = {ratio_num:.4f}).\n"
            )
        elif ratio_line and "50_percent" in ratio_line:
            out_lines.append(
                f"Direction of palindromic SNPs inferred as {palin_mode} by establishing consensus direction of 10% of all sites (forward sites ratio = {ratio_num:.4f}).\n"
            )
        else:
            out_lines.append(
                f"Direction of palindromic SNPs inferred as {palin_mode} by establishing consensus direction (forward sites ratio = {ratio_num:.4f}).\n"
            )
    out_lines.append(
        "\n################################################################\n"
    )

    out_lines.append(summarise_hm_codes(hm_codes))

    with open(args.output, "w") as out:
        out.writelines(out_lines)


if __name__ == "__main__":
    main()
