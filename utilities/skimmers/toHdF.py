import uproot
import pandas as pd
import os

import sys

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--path',  type=str, default="./")
parser.add_argument('--split', type=int, default=-1)
parser.add_argument('--dir',   type=str, default=None)
parser.add_argument('--tree',  type=str, default="2mu2kSkimmedTree")
args = parser.parse_args()


data_path = args.path
data_files = [data_path + f for f in os.listdir(data_path) if f.endswith(".root")]

if args.split >= len(data_files):
    print ("No. data_files: " +  str(len(data_files)))
    print ("Splitting : " + str(args.split))
    print ("Returning")
    sys.exit()
if args.split > -1 and args.split < len(data_files):
    data_files = [data_files[args.split]]

print("HDFing")

for ff in data_files:
    if not os.path.isfile(ff[:-4]+"h5"):
            name = ff[:-5]
            print(name)
            tree = uproot.open(ff)

            if args.dir is not None:
                tree = tree[args.dir][args.tree]
                name = name + "_" + args.dir + "_" + args.tree
            else:
                tree = tree[args.tree]
                name = name + "_" + args.tree
            print (name + ".h5")
            #k = [k for k in tree.keys() if "p4"not in k and "five" not in k]
            tree = tree.pandas.df()
            #s = a.columns.to_series()
            #a.columns = s + s.groupby(s).cumcount().astype(str).replace({'0':''})

            tree.to_hdf(name + ".h5","data",append=False,complevel=0)
            tree = 0
    else:
            print ("Already Exists")
