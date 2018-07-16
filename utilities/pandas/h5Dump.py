import ROOT

from ROOT import TFile,TCanvas
from ROOT import TH1,RooDataSet,gROOT,gDirectory

import rootpy
import root_pandas
from root_pandas import read_root, to_root

import matplotlib
from matplotlib import pyplot as plt

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=Warning)

from rootpy.plotting import Hist, HistStack, Legend, Canvas, Graph
from rootpy.plotting.shapes import Line
from rootpy.plotting.style import get_style, set_style
from rootpy.plotting.utils import draw

style = get_style('ATLAS')
style.SetEndErrorSize(3)
set_style(style)

import pandas as pd
import numpy as np

tt = "/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018A-PromptReco-v3_MINIAOD___20180716_135155/2mu2k-Run2018A-PromptReco-v3.root"
theTest = read_root(tt,ignore=["*p4"],key="rootuple/JPsiPhiTree")

theTest.to_hdf("test.h5","data",append=False)
