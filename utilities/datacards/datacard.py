import ROOT

import pandas as pd
import numpy as np

import argparse

from ROOT import TFile,TCanvas
from ROOT import TH1,RooDataSet,gROOT,gDirectory

import rootpy
import root_pandas
from root_pandas import read_root, to_root

import matplotlib
from matplotlib import pyplot as plt
%matplotlib inline

from ROOT import TMinuit, RooWorkspace

from ROOT import RooGaussian,RooPolynomial,RooDataSet,RooDataHist,TH1,TH1F,RooRealVar,RooArgSet,RooArgList,RooExponential
from ROOT import RooAddPdf,kGreen,kRed,kBlue,RooFit,gROOT,TList,RooGenericPdf,RooProdPdf,RooVoigtian, RooChebychev
from ROOT.RooFit import Components,LineColor,LineStyle,Name,Normalization,Layout,Format,Label,Parameters,Range,Title,Rename, Extended


parser = argparse.ArgumentParser()
parser.add_argument('--path', type=str, default="test.h5")
args = parser.parse_args()

x_min = b_min
x_max = b_max
x_mass = b_mass

data = pd.read_hdf(args.path)

mass = RooRealVar("mass","M(KK)[GeV]",4.5,5.75)

mean = RooRealVar("mean","mean of gaussian",x_mass,x_mass-0.005,x_mass+0.005);
sigma = RooRealVar("sigma","width of gaussian",0.0013);
gamma = RooRealVar("gamma","gamma of bw",0.004253)#,0.001,0.01);

sigma_2 = RooRealVar("sigma","width of gaussian",0.0013);
gFrac = RooRealVar("gFrac","gFrac",5E5,0.,5.0E6)

mean_3 = RooRealVar("mean","mean of gaussian",x_mass,x_mass-0.005,x_mass+0.005);
sigma_3 = RooRealVar("sigma","width of gaussian",0.0013);

sFrac = RooRealVar("sFrac","sFrac",5E5,0.,5.0E6)

bumpFrac = RooRealVar("bumpFracbumpFrac","bumpFrac",5E5,0.,5.0E6)

alpha = RooRealVar("alpha","alpha",-0.1,-1.0,1.0)

#cheb = RooChebychev("cheb","Background",mass,aset)
bkg = RooExponential("bkg","bkg",mass,alpha)
#gauss = RooGaussian("gauss","gaussian PDF ",mass,mean,sigma)
sig_1 = RooGaussian("sig_1","sig_1",mass,mean,sigma)
sig_2 = RooGaussian("sig_1","sig_1",mass,mean,sigma_2)

sig_3 = RooGaussian("bump","bump",mass,mean_3,sigma_3)

sig = RooAddPdf("sig","g+g",sig_1,sig_2,gFrac)

#sig_1 = RooGaussian("sig_1","sig_1",mass,mean,sigma)


tot = RooAddPdf("tot","g+cheb",RooArgList(sig,sig_3,bkg),RooArgList(sFrac,bumpFrac))

h1 = TH1F("hist","hist",00, 4.5,5.75)
map(h1.Fill, spectrum)

masslist = RooArgList(mass)
dh = RooDataHist("dh","dh",masslist,h1)
