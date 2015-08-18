import pytest
from hashlib import md5
import os
from io import StringIO
import re
from copy import deepcopy
from MyFuncs import TempFile

try:
    import workshop.PhyloBuddy as Pb
except ImportError:
    import PhyloBuddy as Pb

def phylo_to_hash(_phylobuddy, mode='hash'):
    if mode != "hash":
        return str(_phylobuddy)
    _hash = md5("{0}\n".format(str(_phylobuddy).rstrip()).encode()).hexdigest()
    return _hash

root_dir = os.getcwd()

def resource(file_name):
    return "{0}/unit_test_resources/{1}".format(root_dir, file_name)

pb_files = ['sample_newick_multi_tree.newick', 'sample']
