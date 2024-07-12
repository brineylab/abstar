# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

from .annotation.germline import *


#  TODO:
#
#  should think about the pros/cons of building a single germline database per species that
#  contains both BCR and TCR germline segments
#
#  the biggest con is the fact that OGRDB only has BCR data, so all databases would need to be
#  hybrid of IMGT and OGRDB datasets -- maybe this just needs some more time??
#
#  functions:
#    - build_database() -- has `schema` kwarg, which can be either "ogrdb" or "imgt"
#    - build_database_ogrdb()
#       - should use their JSON files by default, but constant regions may need to be FASTA (and come from IMGT)
#    - build_database_imgt()
#
#    - show_manifest()
#       - prints the manifest for a given germline database
