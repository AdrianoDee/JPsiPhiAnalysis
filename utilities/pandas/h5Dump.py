import os

import root_pandas
from root_pandas import read_root, to_root

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=Warning)


import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--path', action='store_true')
args = parser.parse_args()

for path, subdirs, files in os.walk(args.path):
    for name in files:
        print os.path.join(path, name)
        tt = os.path.join(path, name)
        theTest = read_root(os.path.join(path, name),ignore=["*p4"],key="rootuple/JPsiPhiTree")
        ##print os.path.join(path, name)
        theTest.to_hdf(os.path.join(path, name)[:-4] + "h5","data")
