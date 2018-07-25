import os
import sys

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reference')

if sys.version_info[0] > 2:
    import pickle
    bnab_hc_antibodies_file = 'bnab_hc_antibodies.pkl'
    bnab_lc_antibodies_file = 'bnab_lc_antibodies.pkl'
else:
    import cPickle as pickle
    bnab_hc_antibodies_file = 'bnab_hc_antibodies.pkl'
    bnab_lc_antibodies_file = 'bnab_lc_antibodies.pkl'
    # bnab_hc_antibodies_file = 'bnab_hc_antibodies_py2.pkl'
    # bnab_lc_antibodies_file = 'bnab_lc_antibodies_py2.pkl'

with open(os.path.join(data_dir, bnab_hc_antibodies_file), 'rb') as f:
    REFERENCE_BNAB_HC_ANTIBODIES = pickle.load(f)

with open(os.path.join(data_dir, bnab_lc_antibodies_file), 'rb') as f:
    REFERENCE_BNAB_LC_ANTIBODIES = pickle.load(f)
