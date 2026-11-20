"""Spawn-pickleable record failure injection for recovery contract tests."""
import importlib

from abstar.annotation.annotator import annotate, annotate_single_sequence

# A real assignable light-chain sequence marks records selected for injection.
# The failure is deliberate and independent of its biological annotation.
INJECTED_FAILURE_SEQUENCE = 'GCTGGGGTCTCAGGAGGCAGCGCTCTCAGGACATCTCCACCATGGCCTGGGCTCTGCTGCTCCTCACCCTCCTCACTCAGGGCACAGGGTCCTGGGCCCAGTCTGCCCTGACTCAGCCTGCCTCCGTGTCTGGGTCTCCTGGACAGTCGATCACCATCTCCTGCGCTGGAACCAGCAGTGACGTTGGTGGTTATAACTCTGTCTCCTGGTACCAACAATACCCAGGCAAAGCCCCCAAACTCGTGATTTATGAGGTCAGTAATCGGCCCTCAGGGGTTACTTATCGCTTCTCTGGCTCCAAGTCTGGCAACACGGCCTCCCTGACCATCTCTGGGCTCCAGGCTGAGGACGAGGCTGATTATTATTGCGCATATGCAACTGACGGCACTCTCGACTTCGGCGGAGGGACGAAGCTGACCGTCCTTGGTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCACCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGTGTGTCTCATAAGTGACTTCTACCCGGGAGCCGTGACAGTGGCCTGGAAGGCAGATAGCAGCCCCGTCAAGGCGGGAGTGGAGACCACCACACCCTCCAAACAAAGCAACAACAAGTACGCGGCCAGCAGCTA'
FAILURE_MESSAGE = "deliberate annotation failure for recovery contract"


def _fail_marked_record(**kwargs):
    if kwargs["ab"].sequence_input == INJECTED_FAILURE_SEQUENCE:
        raise RuntimeError(FAILURE_MESSAGE)
    return annotate_single_sequence(**kwargs)


def annotate_with_record_failure(input_file, **kwargs):
    # Patch inside the worker so spawn multiprocessing receives the injection.
    # Chunk loading, record recovery, diagnostics and output writing remain real.
    module = importlib.import_module("abstar.annotation.annotator")
    original = module.annotate_single_sequence
    module.annotate_single_sequence = _fail_marked_record
    try:
        return annotate(input_file, **kwargs)
    finally:
        module.annotate_single_sequence = original
