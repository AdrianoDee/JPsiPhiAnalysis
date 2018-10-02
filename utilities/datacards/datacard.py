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
from ROOT import RooAddPdf,kGreen,kRed,kBlue,RooFit,gROOT,TList,RooGenericPdf,RooProdPdf,RooVoigtian, RooChebychev, RooBernstein
from ROOT.RooFit import Components,LineColor,LineStyle,Name,Normalization,Layout,Format,Label,Parameters,Range,Title,Rename, Extended


parser = argparse.ArgumentParser()
parser.add_argument('--path', type=str, default="test.h5")
args = parser.parse_args()

x_min = 4.0
x_max = 5.5
x_mass = 5.35

data = pd.read_hdf(args.path)

mass = RooRealVar("mass","M(KK)[GeV]",4.05,5.75)

mean = RooRealVar("mean","mean of gaussian",x_mass,x_mass-0.05,x_mass+0.05);
sigma = RooRealVar("sigma","width of gaussian",0.0013);
gamma = RooRealVar("gamma","gamma of bw",0.004253)#,0.001,0.01);

sigma_2 = RooRealVar("sigma","width of gaussian",0.0013);
gFrac = RooRealVar("gFrac","gFrac",5E5,0.,5.0E6)

mean_3 = RooRealVar("mean","mean of gaussian",x_mass,x_mass-0.005,x_mass+0.005);
sigma_3 = RooRealVar("sigma","width of gaussian",0.0013);

sFrac = RooRealVar("sFrac","sFrac",5E5,0.,5.0E6)

bumpFrac = RooRealVar("bumpFracbumpFrac","bumpFrac",5E5,0.,5.0E6)

alpha = RooRealVar("alpha","alpha",-0.1,-1.0,1.0)

a0 = RooRealVar("a0","a0",0.1,-5.0,5.0)
a1 = RooRealVar("a1","a1",0.1,-5.0,5.0)
a2 = RooRealVar("a2","a2",0.1,-5.0,5.0)
a3 = RooRealVar("a3","a3",0.01,-5.0,5.0)
a4 = RooRealVar("a4","a4",0.01,-5.0,5.0)
a5 = RooRealVar("a5","a5",0.01,-5.0,5.0)
a6 = RooRealVar("a6","a6",0.01,-5.0,5.0)
a7 = RooRealVar("a7","a7",0.001,-5.0,5.0)
a8 = RooRealVar("a8","a8",0.001,-5.0,5.0)
aset = RooArgList(a0,a1,a2,a3,a4,a5,a6,a7,a8)
bkg = RooBernstein("cheb","Background",mass,aset)
#bkg = RooExponential("bkg","bkg",mass,alpha)

#gauss = RooGaussian("gauss","gaussian PDF ",mass,mean,sigma)
sig_1 = RooGaussian("sig_1","sig_1",mass,mean,sigma)
sig_2 = RooGaussian("sig_1","sig_1",mass,mean,sigma_2)

sig_3 = RooGaussian("bump","bump",mass,mean_3,sigma_3)

sig = RooAddPdf("sig","g+g",sig_1,sig_2,gFrac)

#sig_1 = RooGaussian("sig_1","sig_1",mass,mean,sigma)

nSig = RooRealVar("nSig","nSig",100,100,len(data["mass"].values))
nBkg = RooRealVar("nBkg","nBkg",1000,100,len(data["mass"].values))

#tot = RooAddPdf("tot","g+cheb",RooArgList(sig,sig_3,bkg),RooArgList(sFrac,bumpFrac))
tot = RooAddPdf("tot","g+cheb",RooArgList(sig_1,bkg),RooArgList(nSig,nBkg))
h1 = TH1F("hist","hist",200, 4.05,5.75)
map(h1.Fill, data["mass"].values)

masslist = RooArgList(mass)
dh = RooDataHist("dh","dh",masslist,h1)

tot.fitTo(dh)

canvas = TCanvas("c","c",1200,1000)
kkFrame = mass.frame()
dh.plotOn(kkFrame)

tot.plotOn(kkFrame)#,RooFit.Normalization(1.0/float(nfit)))
dh.plotOn(kkFrame)
#tot.paramOn(kkFrame,RooFit.Layout(0.57,0.99,0.65))

kkFrame.Draw()
canvas.Draw()

canvas.SaveAs("test.png")
