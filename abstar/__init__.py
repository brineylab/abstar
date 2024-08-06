import warnings

from Bio import BiopythonWarning

warnings.simplefilter("ignore", BiopythonWarning)
warnings.simplefilter(action="ignore", category=FutureWarning)

from .core.abstar import run
from .version import __version__
