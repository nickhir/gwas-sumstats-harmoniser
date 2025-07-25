#!/usr/bin/env python3
import argparse
import re
import yaml


def parse_log(path):
    with open(path) as f:
        lines = f.readlines()
    text = "".join(lines)
    data = {}
    m = re.search(r"\((\d+) sites out of (\d+)\) were dropped", text)
    if m:
        data["dropped_variants"] = int(m.group(1))
        data["total_variants"] = int(m.group(2))
    m = re.search(r"\((\d+) of (\d+)\) sites successfully harmonised", text)
    if m:
        data["successful_variants"] = int(m.group(1))
        data.setdefault("total_variants", int(m.group(2)))
    m = re.search(r"\((\d+) of (\d+)\) sites failed to harmonise", text)
    if m:
        data["failed_variants"] = int(m.group(1))
        data.setdefault("total_variants", int(m.group(2)))
    fail_codes = {}
    in_fail = False
    for line in lines:
        line = line.strip()
        if line.startswith("7. Failed harmonisation"):
            in_fail = True
            continue
        if in_fail:
            if line.startswith("hm_code") or not line:
                continue
            if line.startswith("###"):
                break
            m = re.match(r"(\d+)\t(\d+)\t", line)
            if m:
                fail_codes[m.group(1)] = int(m.group(2))
    if fail_codes:
        data["failed_hm_codes"] = fail_codes
    return data


def main():
    ap = argparse.ArgumentParser(description="Update metadata YAML with log summary")
    ap.add_argument("yaml_in")
    ap.add_argument("log")
    ap.add_argument("yaml_out")
    args = ap.parse_args()
    with open(args.yaml_in) as fh:
        meta = yaml.safe_load(fh) or {}
    meta.update(parse_log(args.log))
    with open(args.yaml_out, "w") as fh:
        yaml.dump(meta, fh, allow_unicode=True)


if __name__ == "__main__":
    main()
