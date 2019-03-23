import ROOT
from ROOT import TFile,TH1,TH1F,TCanvas,TNtuple,TTreeReader,TTreeReaderValue
from ROOT import RooFit
from ROOT.RooFit import Layout
from ROOT import RooStats
from ROOT import RooAbsData
RooAbsData.setDefaultStorageType ( RooAbsData.Tree )
from array import array
import sys


import numpy as np

# In[3]:


from ROOT import RooRealVar,RooAbsPdf,RooChebychev,RooExponential,RooGaussian,TLine
from ROOT import RooBernstein,RooAbsPdf,RooPlot,RooAddPdf,RooDataHist,RooArgSet,RooArgList
from ROOT import kGreen,kRed,kBlack,kBlue,kDashed,kDotted,kMagenta,RooVoigtian,kOrange
from ROOT.RooFit import Components,LineColor,LineStyle,Name,Normalization,Range,SelectVars,Parameters,Title
from ROOT import RooDataSet,RooFormulaVar,RooLinkedList,RooGenericPdf

import argparse


# In[4]:


parser = argparse.ArgumentParser()
parser.add_argument('--path', type=str, default="./sPlot_psi_2016.root")
parser.add_argument('--ncpu',type=int,default=4)
parser.add_argument('--phsps',action="store_true")
args = parser.parse_args()


rootfile = args.path
inputfile = TFile(rootfile,"READ") 
xTuple = (inputfile.Get("tree")) 


# In[6]:



nentries = xTuple.GetEntries()
print nentries

# In[8]:
binning = 400
bumean = 5.27
bumin = 5.12
bumax = 5.279 + 0.180
massbins = (6.0 - 4.0)/0.02
dimuonditrk_m_rf_c = RooRealVar("dimuonditrk_m_rf_c","M(#mu#muKK)[GeV]",4.0,6.0)
dimuonditrk_m_rf_c.setBins(100)


dimuonditrk_charge    = RooRealVar("dimuonditrk_charge","dimuonditrk_charge",-5.0,5.0)
dimuonditrk_ctauPV    = RooRealVar("dimuonditrk_ctauPV","dimuonditrk_ctauPV",-1000.0,1000.0)
dimuon_pt             = RooRealVar("dimuon_pt","dimuon_pt",0.0,1000.0)
dimuon_m              = RooRealVar("dimuon_m","dimuon_m",2.5,3.5)
dimuonditrk_ctauErrPV = RooRealVar("dimuonditrk_ctauErrPV","dimuonditrk_ctauErrPV",-1000.0,1000.0)
ditrak_pt             = RooRealVar("ditrak_pt","ditrak_pt",0.0,1000.0)
ditrak_m              = RooRealVar("ditrak_m","ditrak_m",0.9,1.5)
mass_ref_c_kkk        = RooRealVar("mass_ref_c_kkk","M(#mu#muKKK)",bumin,bumax)
dimuonditrk_ctauPV    = RooRealVar("dimuonditrk_ctauPV","dimuonditrk_ctauPV",-1000.0,1000.0)
kkk		      = RooRealVar("kkk_m","kkk_m",0.0,1000.0)
jpsik         = RooRealVar("jpsik_m_ref","jpsik_m",0.0,1000.0)


theSet = RooArgSet(mass_ref_c_kkk,dimuonditrk_m_rf_c,dimuonditrk_charge,kkk,jpsik)
   #RooArgSet(RooArgSet(mass_ref_c_kkk,dimuonditrk_m_rf_c,dimuonditrk_ctauPV,dimuonditrk_ctauErrPV,dimuon_pt),RooArgSet(ditrak_pt,mass_ref_c_kkk,dimuonditrk_charge,ditrak_m,dimuon_m))

DATA = RooDataSet("alldata","alldata",xTuple,theSet)

splotData = RooDataSet("alldata_sig","alldata_sig",DATA,theSet,"dimuonditrk_charge==0")
#splotBkgData = RooDataSet("alldata_bkg","alldata_bkg",DATA,theSet,"dimuonditrk_charge!=0 && dimuonditrk_m_rf_c<5.0 && dimuonditrk_m_rf_c>4.0")

print "Tree entries %d"%(splotData.numEntries())

print "PHSP fit"

BkgTotalMPdf = RooGenericPdf("BkgPdf","BkgPdf","sqrt( pow(dimuonditrk_m_rf_c,4) + pow(3.0967,4) + pow(1.01946,4) - 2*pow(dimuonditrk_m_rf_c,2)*pow(3.0967,2) - 2*pow(3.0967,2)*pow(1.01946,2) - 2*pow(dimuonditrk_m_rf_c,2)*pow(1.01946,2) ) * sqrt( pow(5.279,4) + pow(dimuonditrk_m_rf_c,4) + pow(0.493677,4) - 2*pow(5.279,2)*pow(dimuonditrk_m_rf_c,2) - 2*pow(5.279,2)*pow(0.493677,2) - 2*pow(dimuonditrk_m_rf_c,2)*pow(0.493677,2) ) / (dimuonditrk_m_rf_c)", RooArgList(dimuonditrk_m_rf_c));

dimuonditrk_m_rf_c.setBins(80)
dimuonditrk_m_rf_c.setRange("baserange",4.0,5.0)
s = BkgTotalMPdf.createIntegral(RooArgSet(dimuonditrk_m_rf_c),"baserange").getVal()

#bkgFit = BkgTotalMPdf.fitTo(splotBkgData,Range(4.0,5.0),RooFit.NumCPU(args.ncpu),RooFit.Verbose(False))

cb = TCanvas("canvas_b","canvas_b",1200,800) 
print s
mumukkFrame = dimuonditrk_m_rf_c.frame(Title("Phase Space Fit"),Range(4.0,5.0),Normalization(1.0))
splotData.plotOn(mumukkFrame)

BkgTotalMPdf.plotOn(mumukkFrame,Normalization(1.65))

mumukkFrame.Draw()

if args.phsps:
    cb.SaveAs(args.path[:-5] + '_bu_phsp_plot.root')
    cb.SaveAs(args.path[:-5] + '_bu_phsp_plot.png')

    sys.exit()

print "SPLOT FIT"

a0 = RooRealVar("a0","a0",0.001,-10.,10.)
a1 = RooRealVar("a1","a1",0.001,-5.0,5.0)
a2 = RooRealVar("a2","a2",-0.00001,-2.,2.)
a3 = RooRealVar("a3","a3",0.0,-0.5,0.5)
a4 = RooRealVar("a4","a4",0.0,-0.2,0.2)
a5 = RooRealVar("a5","a5",0.0,-0.025,0.05)
a6 = RooRealVar("a6","a6",0.0,-0.001,0.001)

aset = RooArgList(a0,a1,a2,a3,a4)

sigma = RooRealVar("#sigma","width of gaussian",0.011,0.005,0.05)
gamma = RooRealVar("#Gamma","gamma of bw",0.0042,0.001,0.01)
mean = RooRealVar("m","mean of gaussian",bumean,bumean-0.1,bumean+0.1);

nSig = RooRealVar("n_{Sig}","nSig",2E4,0.,6.0E6)
nBkg = RooRealVar("nBkg","nBkg",1E5,0.,6.0E6)
cheb = RooChebychev("cheb","Background",mass_ref_c_kkk,aset)
gauss = RooGaussian("gauss","gaussian PDF ",mass_ref_c_kkk,mean,sigma)

sigma_2 = RooRealVar("#sigma_{2}","width of gaussian",0.01,0.002,0.04);
gFrac_1 = RooRealVar("f_{1}","gFrac",0.1,0.,1.0)
gauss_2 = RooGaussian("gauss_2","gauss_2",mass_ref_c_kkk,mean,sigma_2)

signal = RooAddPdf("sig","g+g",RooArgList(gauss,gauss_2),RooArgList(gFrac_1))

signal = gauss
#signal = RooVoigtian("signal","signal",mass_ref_c_kkk,mean,gamma,sigma)
#signal = gauss

B_1     = RooRealVar ( "B_{1}"    , "B_1 "   , 0.3  , -20   , 100   )
B_2     = RooRealVar ( "B_{2}"    , "B_2"    , 0.3  , -20   , 100   )
B_3     = RooRealVar ( "B_{3}"    , "B_3"    , 0.3  , -20   , 100   )
B_4     = RooRealVar ( "B_{4}"    , "B_4"    , 0.3  , -20   , 100   )

bkg    = RooChebychev("pdfB" , "pdfB"    , mass_ref_c_kkk   , RooArgList(aset))

tot = RooAddPdf("tot","g+cheb",RooArgList(signal,bkg),RooArgList(nSig,nBkg))

nfits = 0

mean.setConstant(True)
sigma.setConstant(True)
rPhifit = tot.fitTo(splotData,Range(bumin,bumax),RooFit.NumCPU(args.ncpu),RooFit.Verbose(False))
nfits = nfits + 1

mean.setConstant(False)
sigma.setConstant(False)
rPhifit = tot.fitTo(splotData,Range(bumin,bumax),RooFit.NumCPU(args.ncpu),RooFit.Verbose(False))
nfits = nfits + 1

a0.setConstant(True)
a1.setConstant(True)
a2.setConstant(True)
a3.setConstant(True)
a4.setConstant(True)
#nBkg.setConstant(True)

c = TCanvas("canvas","canvas",1200,800) 
phiFrame = mass_ref_c_kkk.frame(Title("#mu#muKKK spectrum - B^{+} candidates"),Range(bumin,bumax),Normalization((nSig.getValV() + nBkg.getValV())))
splotData.plotOn(phiFrame)
ratio = 1.0/float(nfits)

tot.plotOn(phiFrame,Normalization(ratio))
bFrac = (nBkg.getValV())/(nSig.getValV() + nBkg.getValV())
bkg.plotOn(phiFrame,LineColor(kRed),Normalization(bFrac),LineStyle(kDashed))
signal.plotOn(phiFrame,LineColor(kGreen),Normalization(1.0-bFrac))
tot.paramOn(phiFrame,RooFit.Layout(0.72,0.99,0.4))#,Parameters(RooArgSet(nSig)))
#signal.paramOn(phiFrame,RooFit.Layout(0.57,0.99,0.65),Parameters(RooArgSet(mean,sigma)))

lowranges = []
uppranges = []

sides = [3.0,5.0,8.0]
for i in sides:
    mass_ref_c_kkk.setRange("sigma_"+ str(i),mean.getVal()-float(i)*sigma.getVal(),mean.getVal()+float(i)*sigma.getVal())
    lowranges.append(mean.getVal()-float(i)*sigma.getVal())
    uppranges.append(mean.getVal()+float(i)*sigma.getVal())
    
    
totIntegralSig = signal.createIntegral(RooArgSet(mass_ref_c_kkk)).getVal()
totIntegralBkg = bkg.analyticalIntegral(bkg.getAnalyticalIntegral(RooArgSet(mass_ref_c_kkk),RooArgSet(mass_ref_c_kkk)))
totIntegralTot = tot.createIntegral(RooArgSet(mass_ref_c_kkk)).getVal()

print "Tot sig: " + str(totIntegralSig)
print "Tot bkg: " + str(totIntegralBkg)
print "Tot tot: " + str(totIntegralTot)
 
sigIntegrals = []
totIntegrals = []
bkgIntegrals = []



phiFrame.Draw()

#sidesigma = np.sqrt(gamma.getValV()**2 + sigma.getValV()**2)
#
plotmax = 1.5 * float(nentries/binning)
#lowside = -3.*sidesigma + mean.getValV()
#upside = +3.*sidesigma + mean.getValV()
#
#linelow = TLine(lowside,0.0,lowside,plotmax*0.95)
##lineup = TLine(upside,0.0,upside,plotmax*0.95)
##linelow.SetLineColor(kGreen)
#lineup.SetLineColor(kGreen)
#linelow.SetLineWidth(2)
#lineup.SetLineWidth(2)

#linelow.Draw()
#lineup.Draw()

#tot.paramOn(phiFrame,RooFit.Layout(0.57,0.99,0.65))


c.SaveAs(args.path[:-5] + '_bumassSPlot.png')
c.SaveAs(args.path[:-5] + '_bumassSPlot.root')

for i in [3.0,5.0,8.0]:
    s = signal.createIntegral(RooArgSet(mass_ref_c_kkk),"sigma_"+ str(i)).getVal()
    b = bkg.createIntegral(RooArgSet(mass_ref_c_kkk),"sigma_"+ str(i)).getVal()
    t = tot.createIntegral(RooArgSet(mass_ref_c_kkk),"sigma_"+ str(i)).getVal()
    
    sigIntegrals.append(s * nSig.getValV() / totIntegralSig)
    totIntegrals.append(b * nBkg.getValV() / totIntegralBkg)
    bkgIntegrals.append(t * (nBkg.getValV()+nSig.getValV()) / totIntegralTot)
    
    print "Sigma " + str(i)
    print "sig : " + str(sigIntegrals[-1])
    print "bkg : " + str(bkgIntegrals[-1])
    print "tot : " + str(totIntegrals[-1])

    
cols = [kGreen,kRed,kRed]#,kMagenta,kOrange,kBlack,kGreen+1]
l_dw = []
l_up = []

for m,M,cc in zip(lowranges,uppranges,cols):
    print m
    print M

    l_dw.append(TLine(m,0.0,m,plotmax*0.95))
    l_up.append(TLine(M,0.0,M,plotmax*0.95))
    
    l_dw[-1].SetLineColor(cc)
    l_up[-1].SetLineColor(cc)
    l_dw[-1].SetLineWidth(2)
    l_up[-1].SetLineWidth(2)

    l_dw[-1].Draw()
    l_up[-1].Draw()

tot.paramOn(phiFrame,RooFit.Layout(0.65,0.99,0.4))
c.SaveAs(args.path[:-5] + '_bumassSPlot_sb.png')
c.SaveAs(args.path[:-5] + '_bumassSPlot_sb.root')

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
cD.SaveAs(args.path[:-5] + '_SigSPlotBu.gif')
cD.SaveAs(args.path[:-5] + '_SigSPlotBu.root')


sweightSig   = TH1F('sweightSig','sweightSig', 50, 4.0, 6.0)
sweightSig.Sumw2()
sweightSig.SetLineColor(2)    
sweightSig.SetMarkerColor(2); sweightSig.SetMinimum(0.)
dstree.Project('sweightSig','nSig_sw');  

sweightSig.Draw('e0');
cD.SaveAs(args.path[:-5] + '_SigWeightPlotBu.gif')
cD.SaveAs(args.path[:-5] + '_SigWeightPlotBu.root')


#S plot hist for bkg
shistBkg   = TH1F('shistBkg','shistBkg', 50, 4.0, 6.0)
shistBkg.Sumw2()
shistBkg.SetLineColor(2)    
shistBkg.SetMarkerColor(2); shistBkg.SetMinimum(0.)
dstree.Project('shistBkg','dimuonditrk_m_rf_c','nBkg_sw');  

shistBkg.Draw('e0');
cD.SaveAs(args.path[:-5] + '_BkgSPlotBu.gif')
cD.SaveAs(args.path[:-5] + '_BkgSPlotBu.root')


sweightBkg   = TH1F('sweightBkg','sweightBkg', 50, 4.0, 6.0)
sweightBkg.Sumw2()
sweightBkg.SetLineColor(2)    
sweightBkg.SetMarkerColor(2); sweightBkg.SetMinimum(0.)
dstree.Project('sweightBkg','nBkg_sw');  

sweightBkg.Draw('e0');
cD.SaveAs(args.path[:-5] + '_BkgWeightPlotBu.gif')
cD.SaveAs(args.path[:-5] + '_BkgWeightPlotBu.root')


outputFile = TFile( args.path + "_outtree_bu.root","RECREATE")
outputFile.cd()

dstree.Write()
outputFile.Close()


for i in range(1,8):
    s = signal.createIntegral(RooArgSet(mass_ref_c_kkk),"sigma_"+ str(i)).getVal()
    b = bkg.createIntegral(RooArgSet(mass_ref_c_kkk),"sigma_"+ str(i)).getVal()
    t = tot.createIntegral(RooArgSet(mass_ref_c_kkk),"sigma_"+ str(i)).getVal()
    
    sigIntegrals.append(s * nSig.getValV() / totIntegralSig)
    totIntegrals.append(b * nBkg.getValV() / totIntegralBkg)
    bkgIntegrals.append(t * (nBkg.getValV()+nSig.getValV()) / totIntegralTot)
    
    print "Sigma " + str(i)
    print "sig : " + str(sigIntegrals[-1])
    print "bkg : " + str(bkgIntegrals[-1])
    print "tot : " + str(totIntegrals[-1])
    
