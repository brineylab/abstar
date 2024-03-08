import warnings

from Bio import BiopythonWarning

warnings.simplefilter("ignore", BiopythonWarning)
warnings.simplefilter(action="ignore", category=FutureWarning)

from .core.abstar import create_parser, main, run, run_standalone, validate_args
from .utils.preprocess import adapter_trim, fastqc, quality_trim
from .version import __version__
