from abstar import run, run_standalone, main, parse_arguments, validate_args
from preprocess import fastqc, adapter_trim, quality_trim

from pkg_resources import get_distribution, DistributionNotFound
import os.path

try:
    _dist = get_distribution('abstar')
    # Normalize case for Windows systems
    dist_loc = os.path.normcase(_dist.location)
    here = os.path.normcase(__file__)
    if not here.startswith(os.path.join(dist_loc, 'abstar')):
        # not installed, but there is another version that *is*
        raise DistributionNotFound
except DistributionNotFound:
    __version__ = 'Please install AbStar before checking the version'
else:
    __version__ = _dist.version
