import warnings
from Bio import BiopythonWarning

warnings.simplefilter("ignore", BiopythonWarning)

from .core.abstar import run, run_standalone, main, create_parser, validate_args
from .preprocess import fastqc, adapter_trim, quality_trim

from .version import __version__
