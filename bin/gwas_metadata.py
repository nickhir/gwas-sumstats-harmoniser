#!/usr/bin/env python

import argparse
import yaml


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


def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-i", help="Input metadata yaml file")
    argparser.add_argument("-o", help="output metadata yaml file")
    argparser.add_argument(
        "-e",
        help="Edit mode, provide params to edit e.g. `-GWASID GCST123456` to edit/add that value",
        action="store_true",
    )
    _, unknown = argparser.parse_known_args()
    for arg in unknown:
        if arg.startswith(("-", "--")):
            arg_pair = arg.split("=")
            argparser.add_argument(arg_pair[0])
    args = argparser.parse_args()
    known_args = {"i", "o", "e"}
    in_file = args.i
    out_file = args.o
    edit_mode = args.e
    m = MetadataClient(in_file=in_file, out_file=out_file)
    if in_file:
        m.from_file()
        print(f"===========\nMetadata in\n===========\n{m}")
    if edit_mode:
        data_dict = {k: v for k, v in vars(args).items() if k not in known_args}
        m.update_metadata(data_dict)
    if out_file:
        print(f"============\nMetadata out\n============\n{m}")
        m.to_file()


if __name__ == "__main__":
    main()
