#!usr/env/python
# filename: test_isotype.py

import csv
import os

from ..utils.isotype import get_isotype


class Antibody:
    def __init__(
        self,
        isotype: str,
        vdj_nt: str,
        oriented_input: str,
        species: str = "human",
        germ_db: str = "human",
    ):
        self.isotype = isotype
        self.vdj_nt = vdj_nt
        self.oriented_input = oriented_input
        self.species = species
        self.germ_db = germ_db


def test_get_isotype():
    with open(os.path.abspath("abstar/test_data/test_get_isotype.csv"), "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            antibody = Antibody(**row)
            isotype = get_isotype(antibody)
            assert isotype.isotype == antibody.isotype
