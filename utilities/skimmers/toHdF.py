import uproot
import pandas as pd
import os

import sys

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--path', type=str, default="./")
parser.add_argument('--split', type=int, default=-1)
args = parser.parse_args()


data_path = args.path
data_files = [data_path + f for f in os.listdir(data_path) if f.endswith(".root")]

if args.split >= len(data_files):
    print "No. data_files: " +  str(len(data_files))
    print "Splitting : " + str(args.split)
    print "Returning"
    sys.exit()
if args.split > -1 and args.split < len(data_files):
    data_files = [data_files[args.split]]

print "HDFing"

for ff in data_files:
    if not os.path.isfile(ff[:-4]+"h5"):
            print ff[:-4]+"h5"
            tree = uproot.open(ff)["2mu2kSkimmedTree"]
            k = [k for k in tree.keys() if "p4"not in k and "five" not in k]
            a = tree.pandas.df(k)
            s = a.columns.to_series()
            a.columns = s + s.groupby(s).cumcount().astype(str).replace({'0':''})
            a.to_hdf(ff[:-4]+"h5","data",append=False)
            a = 0
    else:
            print "Already Exists"
~
