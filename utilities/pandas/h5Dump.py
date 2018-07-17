import os

import root_pandas
from root_pandas import read_root, to_root

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=Warning)


import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--path', type=str,default=".")
args = parser.parse_args()

for path, subdirs, files in os.walk(args.path):
    for name in files:
        print os.path.join(path, name)
        tt = os.path.join(path, name)
        if os.path.isfile(tt[:-4]+"h5") or tt.endswith("h5"):
            print "done"
            continue
        try:
            theTest = read_root(os.path.join(path, name),ignore=["*p4"],key="rootuple/JPsiPhiTree") as theTest:
            theTest.to_hdf(tt[:-4] + "h5","data")
            ##print os.path.join(path, name)
        except:
            print "Not good file"
