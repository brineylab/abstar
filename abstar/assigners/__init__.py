from os.path import dirname, basename, isfile
import glob

modules = glob.glob(dirname(__file__) + "/*.py")
modules = [basename(f)[:-3] for f in modules if all([isfile(f), f not in ['__init__', 'assigner', 'registry']])]

__all__ = modules
