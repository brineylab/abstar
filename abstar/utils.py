# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

import click


def parse_dict_from_string(ctx, param, value):
    if not value:
        return {}
    try:
        # Split the string into key-value pairs
        kv_pairs = value.split(",")
        # Convert the pairs into a dictionary
        return dict(kv_pair.split("=") for kv_pair in kv_pairs)
    except ValueError:
        raise click.BadParameter("Format must be 'key1=val1,key2=val2'")
