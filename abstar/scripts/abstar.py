# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

import click

from ..core.abstar import run
from ..core.germline import build_germline_database


@click.group()
def cli():
    pass


cli.add_command(run)
cli.add_command(build_germline_database)
