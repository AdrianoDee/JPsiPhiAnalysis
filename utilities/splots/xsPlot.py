import ROOT
from ROOT import TFile,TH1,TH1F,TCanvas,TNtuple,TTreeReader,TTreeReaderValue
from ROOT import RooFit
from ROOT.RooFit import Layout
from ROOT import RooStats
from ROOT import RooAbsData
RooAbsData.setDefaultStorageType ( RooAbsData.Tree )
from array import array
import sys

import argparse
import numpy as np

from string import split

from ROOT import RooRealVar,RooAbsPdf,RooChebychev,RooExponential,RooGaussian,TLine
from ROOT import RooBernstein,RooAbsPdf,RooPlot,RooAddPdf,RooDataHist,RooArgSet,RooArgList
from ROOT import kGreen,kRed,kBlack,kBlue,kDashed,kDotted,kMagenta,RooVoigtian
from ROOT.RooFit import Components,LineColor,LineStyle,Name,Normalization,Range,SelectVars, Title
from ROOT import RooDataSet,RooFormulaVar,RooLinkedList

parser = argparse.ArgumentParser()
parser.add_argument('--path', type=str, default="./sPlot_psi_2016.root")
parser.add_argument('--ncpu',type=int,default=20)
args = parser.parse_args()

phimean = 1.020
phimin = 1.020-0.015
phimax = 1.020+0.015

rootfile = args.path
inputfile = TFile(rootfile,"READ")
xTuple = (inputfile.Get("tree"))

outname = split(split(rootfile,"/")[-1],".")[0]

# In[6]:



nentries = xTuple.GetEntries()
print nentries

# In[8]:
binning = 400
massbins = (phimax - phimean)/0.005
dimuonditrk_m_rf_c = RooRealVar("dimuonditrk_m_rf_c","M(#mu#muKK)[GeV]",4.0,6.0)
dimuonditrk_m_rf_c.setBins(binning)
masskk = RooRealVar("ditrak_m","M(KK) [GeV]",1.00,1.045);
masskk.setBins(int(200))

dimuonditrk_ctauPV    = RooRealVar("dimuonditrk_ctauPV","dimuonditrk_ctauPV",-1000.0,1000.0)
dimuon_pt             = RooRealVar("dimuon_pt_x","dimuon_pt_x",0.0,1000.0)
dimuonditrk_ctauErrPV = RooRealVar("dimuonditrk_ctauErrPV","dimuonditrk_ctauErrPV",-1000.0,1000.0)
ditrak_pt             = RooRealVar("ditrak_pt","ditrak_pt",0.0,1000.0)

theSet = RooArgSet(masskk,dimuonditrk_m_rf_c,dimuonditrk_ctauPV,dimuonditrk_ctauErrPV,dimuon_pt,ditrak_pt)
splotData = RooDataSet("alldata","alldata",xTuple,theSet)
#
print "Tree entries %d"%(splotData.numEntries())



a0 = RooRealVar("a_{0}","a0",0.001,-10.,10.)
a1 = RooRealVar("a_{1}","a1",0.001,-5.0,5.0)
a2 = RooRealVar("a_{2}","a2",-0.00001,-2.,2.)
a3 = RooRealVar("a_{3}","a3",0.0,-0.5,0.5)
a4 = RooRealVar("a_{4}","a4",0.0,-0.2,0.2)
a5 = RooRealVar("a5","a5",0.0,-0.025,0.05)
a6 = RooRealVar("a6","a6",0.0,-0.001,0.001)

aset = RooArgList(a0,a1,a2,a3,a4,a5)

sigma = RooRealVar("#sigma","width of gaussian",0.01,0.001,0.05)
gamma = RooRealVar("#Gamma","gamma of bw",0.0042,0.001,0.01)
mean = RooRealVar("#mu","mean of gaussian",phimean,phimean-0.1,phimean+0.1);

nSig = RooRealVar("n_{Sig}","nSig",float(nentries)*0.1,0.,float(nentries))
nBkg = RooRealVar("n_{Bkg}","nBkg",float(nentries)*0.9,0.,float(nentries))
cheb = RooChebychev("cheb","Background",masskk,aset)
gauss = RooGaussian("gauss","gaussian PDF ",masskk,mean,sigma)
#signal = RooVoigtian("signal","signal",psiPrimeMass,mean,gamma,sigma)
signal = gauss

B_1     = RooRealVar ( "B_{1}"    , "B_1 "   , 0.3  , -20   , 100   )
B_2     = RooRealVar ( "B_{2}"    , "B_2"    , 0.3  , -20   , 100   )
B_3     = RooRealVar ( "B_{3}"    , "B_3"    , 0.3  , -20   , 100   )
B_4     = RooRealVar ( "B_{4}"    , "B_4"    , 0.3  , -20   , 100   )

bkg    = RooChebychev("pdfB" , "pdfB"    , masskk   , RooArgList(aset))

tot = RooAddPdf("tot","g+cheb",RooArgList(signal,bkg),RooArgList(nSig,nBkg))

nfits = 0

mean.setConstant(True)
gamma.setConstant(True)
rPhifit = tot.fitTo(splotData,Range(phimin,phimax),RooFit.NumCPU(args.ncpu),RooFit.Verbose(False))
nfits = nfits + 1

mean.setConstant(True)
gamma.setConstant(False)
rPhifit = tot.fitTo(splotData,Range(phimin,phimax),RooFit.NumCPU(args.ncpu),RooFit.Verbose(False))
nfits = nfits + 1

mean.setConstant(False)
gamma.setConstant(False)
rPhifit = tot.fitTo(splotData,Range(phimin,phimax),RooFit.NumCPU(args.ncpu),RooFit.Verbose(False))
nfits = nfits + 1

c = TCanvas("canvas","canvas",1200,800)
phiFrame = masskk.frame(Range(phimin,phimax),Normalization((nSig.getValV() + nBkg.getValV())), Title("#Phi Mass"))
splotData.plotOn(phiFrame)
ratio = 1.0/float(nfits)

tot.plotOn(phiFrame,Normalization(ratio))
bFrac = (nBkg.getValV())/(nSig.getValV() + nBkg.getValV())
bkg.plotOn(phiFrame,LineColor(kRed),Normalization(bFrac),LineStyle(kDashed))
signal.plotOn(phiFrame,LineColor(kGreen),Normalization(1.0-bFrac))
tot.paramOn(phiFrame,RooFit.Layout(0.57,0.99,0.65))

phiFrame.Draw()

sidesigma = np.sqrt(gamma.getValV()**2 + sigma.getValV()**2)

plotmax = 1.5 * float(nentries/binning)
lowside = -3.*sidesigma + mean.getValV()
upside = +3.*sidesigma + mean.getValV()

linelow = TLine(lowside,0.0,lowside,plotmax)
lineup = TLine(upside,0.0,upside,plotmax)
linelow.SetLineColor(kGreen)
lineup.SetLineColor(kGreen)
linelow.SetLineWidth(2)
lineup.SetLineWidth(2)

#linelow.Draw()
#lineup.Draw()

#tot.paramOn(phiFrame,RooFit.Layout(0.57,0.99,0.65))


c.SaveAs('phimassSPlot_'+ outname + '.png')
c.SaveAs('phimassSPlot_' + outname + '.root')
c.Clear()


cD=TCanvas("cD","cD",750,600)
cD.cd()
splot   = RooStats.SPlot ("sPlot","sPlot",splotData, tot, RooArgList(nSig,nBkg))
dstree  = splotData.store().tree()


#S plot hist for signal
shistSig   = TH1F('shistSig','shistSig', 50, 4.0, 6.0)
shistSig.Sumw2()
shistSig.SetLineColor(2)
shistSig.SetMarkerColor(2); shistSig.SetMinimum(0.)
dstree.Project('shistSig','dimuonditrk_m_rf_c','nSig_sw');

shistSig.Draw('e0');
cD.SaveAs('SigSPlotPhi_' + outname + '_' + outname + '.gif  ')
cD.SaveAs('SigSPlotPhi_' + outname + '.root')


sweightSig   = TH1F('sweightSig','sweightSig', 50, 4.0, 6.0)
sweightSig.Sumw2()
sweightSig.SetLineColor(2)
sweightSig.SetMarkerColor(2); sweightSig.SetMinimum(0.)
dstree.Project('sweightSig','nSig_sw');

sweightSig.Draw('e0');
cD.SaveAs('SigWeightPlotPhi_' + outname + '.gif ')
cD.SaveAs('SigWeightPlotPhi_' + outname + '.root')


#S plot hist for bkg
shistBkg   = TH1F('shistBkg','shistBkg', 50, 4.0, 6.0)
shistBkg.Sumw2()
shistBkg.SetLineColor(2)
shistBkg.SetMarkerColor(2); shistBkg.SetMinimum(0.)
dstree.Project('shistBkg','dimuonditrk_m_rf_c','nBkg_sw');

shistBkg.Draw('e0');
cD.SaveAs('BkgSPlotPhi_' + outname + '.gif ')
cD.SaveAs('BkgSPlotPhi_' + outname + '.root')


sweightBkg   = TH1F('sweightBkg','sweightBkg', 50, 4.0, 6.0)
sweightBkg.Sumw2()
sweightBkg.SetLineColor(2)
sweightBkg.SetMarkerColor(2); sweightBkg.SetMinimum(0.)
dstree.Project('sweightBkg','nBkg_sw');

sweightBkg.Draw('e0');
cD.SaveAs('BkgWeightPlotPhi_' + outname + '.gif ')
cD.SaveAs('BkgWeightPlotPhi_' + outname + '.root')


outputFile = TFile('sPlot2016_outtree_psi_' + outname + '.root',"RECREATE")
outputFile.cd()

dstree.Write()
outputFile.Close()
