import os
import pickle

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reference')

with open(os.path.join(data_dir, 'bnab_hc_antibodies.pkl'), 'rb') as f:
    REFERENCE_BNAB_HC_ANTIBODIES = pickle.load(f)

with open(os.path.join(data_dir, 'bnab_lc_antibodies.pkl'), 'rb') as f:
    REFERENCE_BNAB_LC_ANTIBODIES = pickle.load(f)

