#!/usr/bin/env python

import argparse
import yaml
import re


class MetadataClient:
    def __init__(self, meta_dict=None, in_file=None, out_file=None) -> None:
        self._meta_dict = meta_dict or {}
        self._in_file = in_file
        self._out_file = out_file

    def from_file(self) -> None:
        with open(self._in_file, "r") as fh:
            self._meta_dict = yaml.safe_load(fh) or {}

    def to_file(self) -> None:
        with open(self._out_file, "w") as fh:
            yaml.dump(self._meta_dict, fh, encoding="utf-8", allow_unicode=True)

    def update_metadata(self, data_dict) -> None:
        """
        Create a copy of the model and update (no validation occurs)
        """
        for k, v in data_dict.items():
            if v is not None:
                self._meta_dict[k] = v

    def __repr__(self) -> str:
        """
        YAML str representation of metadata.
        """
        return yaml.dump(self._meta_dict)


def parse_log(path):
    """Extract summary statistics from a harmonisation log file."""
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
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-i", help="Input metadata yaml file")
    argparser.add_argument("-o", help="output metadata yaml file")
    argparser.add_argument(
        "-e",
        help="Edit mode, provide params to edit e.g. `-GWASID GCST123456` to edit/add that value",
        action="store_true",
    )
    argparser.add_argument("--log", help="Harmonisation log file to parse")
    _, unknown = argparser.parse_known_args()
    for arg in unknown:
        if arg.startswith(("-", "--")):
            arg_pair = arg.split("=")
            argparser.add_argument(arg_pair[0])
    args = argparser.parse_args()
    known_args = {"i", "o", "e", "log"}
    in_file = args.i
    out_file = args.o
    edit_mode = args.e
    log_file = args.log

    m = MetadataClient(in_file=in_file, out_file=out_file)
    if in_file:
        m.from_file()
        print(f"===========\nMetadata in\n===========\n{m}")
    if edit_mode:
        data_dict = {k: v for k, v in vars(args).items() if k not in known_args}
        m.update_metadata(data_dict)
    if log_file:
        m.update_metadata(parse_log(log_file))
    if out_file:
        print(f"============\nMetadata out\n============\n{m}")
        m.to_file()


if __name__ == "__main__":
    main()
