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

phimass = 1.019461

phimean = 1.020
phimin = 0.99
phimax = 1.048

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
masskk = RooRealVar("ditrak_m","M(KK) [GeV]",phimin,phimax);
masskk.setBins(int(200))

dimuonditrk_ctauPV    = RooRealVar("dimuonditrk_ctauPV","dimuonditrk_ctauPV",-1000.0,1000.0)
dimuon_pt             = RooRealVar("dimuon_pt","dimuon_pt",0.0,1000.0)
dimuonditrk_ctauErrPV = RooRealVar("dimuonditrk_ctauErrPV","dimuonditrk_ctauErrPV",-1000.0,1000.0)
ditrak_pt             = RooRealVar("ditrak_pt","ditrak_pt",0.0,1000.0)

theSet = RooArgSet(masskk,dimuonditrk_m_rf_c,dimuonditrk_ctauPV,dimuonditrk_ctauErrPV,dimuon_pt,ditrak_pt)
splotData = RooDataSet("alldata","alldata",xTuple,theSet)
#
print "Tree entries %d"%(splotData.numEntries())



a0 = RooRealVar("a_{0}","a0",3.92267e-01,-10.,10.)
a1 = RooRealVar("a_{1}","a1",-1.14362e-01,-5.0,5.0)
a2 = RooRealVar("a_{2}","a2",2.43201e-02,-2.,2.)
a3 = RooRealVar("a_{3}","a3",-2.42999e-02,-0.5,0.5)
a4 = RooRealVar("a_{4}","a4",0.01,-0.2,0.2)
a5 = RooRealVar("a_{5}","a5",0.0,-0.025,0.05)
a6 = RooRealVar("a6","a6",0.0,-0.001,0.001)

aset = RooArgList(a0,a1,a2,a3,a4)#,a5)

''' 2018 low cuts fit
1  #Gamma       6.41369e-03   4.22079e-05   2.09351e-03   6.48765e+00
   2  #mu          1.01959e+00   7.91850e-06   4.05304e-04  -4.13901e-02
   3  #sigma       1.00002e-03   9.55853e-08   5.00000e-01  -1.56957e+00
   4  a_{0}        3.92267e-01   4.74749e-04   2.09470e-05   3.92368e-02
   5  a_{1}       -1.14362e-01   7.12943e-04   4.41214e-05  -2.28745e-02
   6  a_{2}        2.43201e-02   4.14706e-04   9.36506e-05   1.21604e-02
   7  a_{3}       -2.42999e-02   6.27522e-04   3.86717e-04  -4.86190e-02
   8  n_{Bkg}      1.87466e+07   1.36843e+04   3.48814e-04   8.52813e-01
   9  n_{Sig}      1.92959e+06   1.30579e+04   3.90748e-04   3.83492e+00
'''
sigma_1 = RooRealVar("#sigma_{1}","width of gaussian",1.6002e-03,0.001,0.004)
sigma_2 = RooRealVar("#sigma_{2}","width of gaussian",1.6002e-03,0.001,0.004)

rFrac = RooRealVar("f_{1}","f1",0.5,0.0,1.0)
gamma = RooRealVar("#Gamma","gamma of bw",4.41369e-03,0.004,0.005)
mean = RooRealVar("m","mean of gaussian",1.01959e+00,phimean-0.005,phimean+0.005);

nSig = RooRealVar("n_{Sig}","nSig",float(nentries)*0.01,0.,float(nentries)*0.5)
nBkg = RooRealVar("n_{Bkg}","nBkg",float(nentries)*0.9,0.,float(nentries))
cheb = RooChebychev("cheb","Background",masskk,aset)
res = RooGaussian("gauss","gaussian PDF ",masskk,mean,sigma_1)
peak = RooVoigtian("peak","peak",masskk,mean,gamma,sigma_2)
signal = RooAddPdf("signal","signal",res,peak,rFrac)

B_1     = RooRealVar ( "B_{1}"    , "B_1 "   , 0.3  , -20   , 100   )
B_2     = RooRealVar ( "B_{2}"    , "B_2"    , 0.3  , -20   , 100   )
B_3     = RooRealVar ( "B_{3}"    , "B_3"    , 0.3  , -20   , 100   )
B_4     = RooRealVar ( "B_{4}"    , "B_4"    , 0.3  , -20   , 100   )

bkg    = RooChebychev("pdfB" , "pdfB"    , masskk   , RooArgList(aset))

tot = RooAddPdf("tot","g+cheb",RooArgList(signal,bkg),RooArgList(nSig,nBkg))

nfits = 0

mean.setConstant(True)
gamma.setConstant(True)
sigma_1.setConstant(True)
sigma_2.setConstant(True)
rPhifit = tot.fitTo(splotData,Range(phimin,phimax),RooFit.NumCPU(args.ncpu),RooFit.Verbose(False))
nfits = nfits + 1

#gamma.setConstant(False)
#sigma.setConstant(False)
#rPhifit = tot.fitTo(splotData,Range(phimin,phimax),RooFit.NumCPU(args.ncpu),RooFit.Verbose(False))
#nfits = nfits + 1

mean.setConstant(False)
gamma.setConstant(False)
sigma_1.setConstant(False)
sigma_2.setConstant(False)
rPhifit = tot.fitTo(splotData,Range(phimin,phimax),RooFit.NumCPU(args.ncpu),RooFit.Verbose(False))
nfits = nfits + 1

c = TCanvas("canvas","canvas",1200,900)
phiFrame = masskk.frame(Range(phimin,phimax),Normalization((nSig.getValV() + nBkg.getValV())), Title("#phi Mass"))
splotData.plotOn(phiFrame)
ratio = 1.0/float(nfits)

tot.plotOn(phiFrame,Normalization(ratio))
bFrac = (nBkg.getValV())/(nSig.getValV() + nBkg.getValV())
bkg.plotOn(phiFrame,LineColor(kRed),Normalization(bFrac),LineStyle(kDashed))
signal.plotOn(phiFrame,LineColor(kGreen),Normalization(1.0-bFrac))

a0.setConstant(True)
a1.setConstant(True)
a2.setConstant(True)
a3.setConstant(True)
a4.setConstant(True)
nBkg.setConstant(True)

tot.paramOn(phiFrame,RooFit.Layout(0.57,0.99,0.65))

phiFrame.Draw()

sidesigma = np.sqrt((rFrac.getValV())*sigma_1.getValV()**2 + (1.0 - (rFrac.getValV()))*sigma_2.getValV()**2)

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

sigevents = ()


masskk.setRange("phiregion_five",phimass-0.005,phimass+0.005)
masskk.setRange("phiregion_ten",phimass-0.01,phimass+0.01)
masskk.setRange("totrange",phimin,phimax)

totIntegralSig = signal.createIntegral(RooArgSet(masskk),"totrange").getVal()
totIntegralBkg = bkg.analyticalIntegral(bkg.getAnalyticalIntegral(RooArgSet(masskk),RooArgSet(masskk)),"totrange")
totIntegralTot = tot.createIntegral(RooArgSet(masskk),"totrange").getVal()

fiveIntegralSig = signal.createIntegral(RooArgSet(masskk),"phiregion_five").getVal()
tenIntegralSig  = signal.createIntegral(RooArgSet(masskk),"phiregion_ten").getVal()

fiveIntegralBkg = bkg.analyticalIntegral(bkg.getAnalyticalIntegral(RooArgSet(masskk),RooArgSet(masskk)),"phiregion_five")
tenIntegralBkg  = bkg.analyticalIntegral(bkg.getAnalyticalIntegral(RooArgSet(masskk),RooArgSet(masskk)),"phiregion_ten")

fiveIntegralTot = tot.createIntegral(RooArgSet(masskk),"phiregion_five").getVal()
tenIntegralTot  = tot.createIntegral(RooArgSet(masskk),"phiregion_ten").getVal()

fiveIntegralSig = fiveIntegralSig * nSig.getValV() / totIntegralSig
tenIntegralSig  = tenIntegralSig * nSig.getValV() / totIntegralSig

fiveIntegralBkg = fiveIntegralBkg * nBkg.getValV() / totIntegralBkg
tenIntegralBkg  = tenIntegralBkg * nBkg.getValV() / totIntegralBkg

fiveIntegralTot = fiveIntegralTot * (nBkg.getValV()+nSig.getValV()) / totIntegralTot
tenIntegralTot  = tenIntegralTot * (nBkg.getValV()+nSig.getValV()) / totIntegralTot

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
cD.SaveAs('SigSPlotPhi_' + outname + '_' + outname + '.gif')
cD.SaveAs('SigSPlotPhi_' + outname + '.root')


sweightSig   = TH1F('sweightSig','sweightSig', 50, 4.0, 6.0)
sweightSig.Sumw2()
sweightSig.SetLineColor(2)
sweightSig.SetMarkerColor(2); sweightSig.SetMinimum(0.)
dstree.Project('sweightSig','nSig_sw');

sweightSig.Draw('e0');
cD.SaveAs('SigWeightPlotPhi_' + outname + '.gif')
cD.SaveAs('SigWeightPlotPhi_' + outname + '.root')


#S plot hist for bkg
shistBkg   = TH1F('shistBkg','shistBkg', 50, 4.0, 6.0)
shistBkg.Sumw2()
shistBkg.SetLineColor(2)
shistBkg.SetMarkerColor(2); shistBkg.SetMinimum(0.)
dstree.Project('shistBkg','dimuonditrk_m_rf_c','nBkg_sw');

shistBkg.Draw('e0');
cD.SaveAs('BkgSPlotPhi_' + outname + '.gif')
cD.SaveAs('BkgSPlotPhi_' + outname + '.root')


sweightBkg   = TH1F('sweightBkg','sweightBkg', 50, 4.0, 6.0)
sweightBkg.Sumw2()
sweightBkg.SetLineColor(2)
sweightBkg.SetMarkerColor(2); sweightBkg.SetMinimum(0.)
dstree.Project('sweightBkg','nBkg_sw');

sweightBkg.Draw('e0');
cD.SaveAs('BkgWeightPlotPhi_' + outname + '.gif')
cD.SaveAs('BkgWeightPlotPhi_' + outname + '.root')


outputFile = TFile('sPlot_' + outname + '_outtree.root',"RECREATE")
outputFile.cd()

dstree.Write()
outputFile.Close()

print "Five MeV sig : %.2f"%(fiveIntegralSig)
print "Ten  MeV sig : %.2f"%(tenIntegralSig)

print "Five MeV bkg : %.2f"%(fiveIntegralBkg)
print "Ten  MeV bkg : %.2f"%(tenIntegralBkg)

print "Five MeV tot : %.2f"%(fiveIntegralTot)
print "Ten  MeV tot : %.2f"%(tenIntegralTot)

