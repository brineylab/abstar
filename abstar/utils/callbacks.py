# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

from typing import Optional

import click


def parse_dict_from_string(
    ctx: click.Context, param: click.Parameter, value: Optional[str]
) -> dict:
    """
    Parse a string of key-value pairs into a dictionary. Used as a ``click`` callback function
    for parsing command-line arguments.

    Parameters
    ----------
    ctx : click.Context
        The click context.

    param : click.Parameter
        The click parameter.

    value : Optional[str]
        The value to parse. If ``None``, return an empty dictionary.

    Returns
    -------
    dict
        The parsed dictionary.

    Raises
    ------
    click.BadParameter
        If the format is not similar to ``key1=val1,key2=val2``.

    """
    if not value:
        return {}
    try:
        # split the string into key-value pairs
        kv_pairs = value.split(",")
        # convert the pairs into a dictionary
        return dict(kv_pair.split("=") for kv_pair in kv_pairs)
    except ValueError:
        raise click.BadParameter("Format must be 'key1=val1,key2=val2'")
