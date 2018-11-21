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
from ROOT import kGreen,kRed,kBlack,kBlue,kDashed,kDotted,kMagenta,RooVoigtian
from ROOT.RooFit import Components,LineColor,LineStyle,Name,Normalization,Range,SelectVars
from ROOT import RooDataSet,RooFormulaVar,RooLinkedList


# In[4]:


rootfile = "./sPlot_2018_bu.root"
inputfile = TFile(rootfile,"READ") 
xTuple = (inputfile.Get("tree")) 


# In[6]:



nentries = xTuple.GetEntries()
print nentries

# In[8]:
binning = 400
bumean = 5.27
bumin = 5.27 - 0.2
bumax = 5.27 + 0.2
massbins = (6.0 - 4.0)/0.005
dimuonditrk_m_rf_c = RooRealVar("dimuonditrk_m_rf_c","M(#mu#muKK)[GeV]",4.0,6.0)
dimuonditrk_m_rf_c.setBins(binning)
masskk = RooRealVar("ditrak_m_x","M(KK)",0.99,1.05);
masskk.setBins(int(200))

dimuonditrk_ctauPV    = RooRealVar("dimuonditrk_ctauPV","dimuonditrk_ctauPV",-1000.0,1000.0)
dimuon_pt             = RooRealVar("dimuon_pt_x","dimuon_pt_x",0.0,1000.0)
dimuonditrk_ctauErrPV = RooRealVar("dimuonditrk_ctauErrPV","dimuonditrk_ctauErrPV",-1000.0,1000.0)
ditrak_pt             = RooRealVar("ditrak_pt","ditrak_pt",0.0,1000.0)
mass_ref_c_kkk        = RooRealVar("mass_ref_c_kkk","mass_ref_c_kkk",bumin,bumax)

theSet = RooArgSet(masskk,dimuonditrk_m_rf_c,dimuonditrk_ctauPV,dimuonditrk_ctauErrPV,dimuon_pt,ditrak_pt,mass_ref_c_kkk)
splotData = RooDataSet("alldata","alldata",xTuple,theSet)
#
print "Tree entries %d"%(splotData.numEntries())



a0 = RooRealVar("a0","a0",0.001,-10.,10.)
a1 = RooRealVar("a1","a1",0.001,-5.0,5.0)
a2 = RooRealVar("a2","a2",-0.00001,-2.,2.)
a3 = RooRealVar("a3","a3",0.0,-0.5,0.5)
a4 = RooRealVar("a4","a4",0.0,-0.2,0.2)
a5 = RooRealVar("a5","a5",0.0,-0.025,0.05)
a6 = RooRealVar("a6","a6",0.0,-0.001,0.001)

aset = RooArgList(a0,a1,a2,a3,a4,a5)

sigma = RooRealVar("sigma","width of gaussian",0.01,0.005,0.05)
gamma = RooRealVar("#Gamma","gamma of bw",0.0042,0.001,0.01)
mean = RooRealVar("mean","mean of gaussian",bumean,bumean-0.1,bumean+0.1);

nSig = RooRealVar("nSig","nSig",2E4,0.,5.0E6)
nBkg = RooRealVar("nBkg","nBkg",1E5,0.,5.0E6)
cheb = RooChebychev("cheb","Background",mass_ref_c_kkk,aset)
gauss = RooGaussian("gauss","gaussian PDF ",mass_ref_c_kkk,mean,sigma)
#signal = RooVoigtian("signal","signal",mass_ref_c_kkk,mean,gamma,sigma)
signal = gauss

B_1     = RooRealVar ( "B_{1}"    , "B_1 "   , 0.3  , -20   , 100   )
B_2     = RooRealVar ( "B_{2}"    , "B_2"    , 0.3  , -20   , 100   )
B_3     = RooRealVar ( "B_{3}"    , "B_3"    , 0.3  , -20   , 100   )
B_4     = RooRealVar ( "B_{4}"    , "B_4"    , 0.3  , -20   , 100   )

bkg    = RooChebychev("pdfB" , "pdfB"    , mass_ref_c_kkk   , RooArgList(aset))

tot = RooAddPdf("tot","g+cheb",RooArgList(signal,bkg),RooArgList(nSig,nBkg))

nfits = 0

mean.setConstant(True)
gamma.setConstant(True)
rPhifit = tot.fitTo(splotData,Range(bumin,bumax),RooFit.NumCPU(8),RooFit.Verbose(False))
nfits = nfits + 1

mean.setConstant(True)
gamma.setConstant(False)
rPhifit = tot.fitTo(splotData,Range(bumin,bumax),RooFit.NumCPU(8),RooFit.Verbose(False))
nfits = nfits + 1

mean.setConstant(False)
gamma.setConstant(False)
rPhifit = tot.fitTo(splotData,Range(bumin,bumax),RooFit.NumCPU(8),RooFit.Verbose(False))
nfits = nfits + 1

c = TCanvas("canvas","canvas",1200,800) 
phiFrame = mass_ref_c_kkk.frame(Range(bumin,bumax),Normalization((nSig.getValV() + nBkg.getValV())))
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

linelow.Draw()
lineup.Draw()

#tot.paramOn(phiFrame,RooFit.Layout(0.57,0.99,0.65))


c.SaveAs("bumassSPlot.png")
c.SaveAs("bumassSPlot.root")
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
cD.SaveAs('SigSPlotBu.gif')
cD.SaveAs('SigSPlotBu.root')


sweightSig   = TH1F('sweightSig','sweightSig', 50, 4.0, 6.0)
sweightSig.Sumw2()
sweightSig.SetLineColor(2)    
sweightSig.SetMarkerColor(2); sweightSig.SetMinimum(0.)
dstree.Project('sweightSig','nSig_sw');  

sweightSig.Draw('e0');
cD.SaveAs('SigWeightPlotBu.gif')
cD.SaveAs('SigWeightPlotBu.root')


#S plot hist for bkg
shistBkg   = TH1F('shistBkg','shistBkg', 50, 4.0, 6.0)
shistBkg.Sumw2()
shistBkg.SetLineColor(2)    
shistBkg.SetMarkerColor(2); shistBkg.SetMinimum(0.)
dstree.Project('shistBkg','dimuonditrk_m_rf_c','nBkg_sw');  

shistBkg.Draw('e0');
cD.SaveAs('BkgSPlotBu.gif')
cD.SaveAs('BkgSPlotBu.root')


sweightBkg   = TH1F('sweightBkg','sweightBkg', 50, 4.0, 6.0)
sweightBkg.Sumw2()
sweightBkg.SetLineColor(2)    
sweightBkg.SetMarkerColor(2); sweightBkg.SetMinimum(0.)
dstree.Project('sweightBkg','nBkg_sw');  

sweightBkg.Draw('e0');
cD.SaveAs('BkgWeightPlotBu.gif')
cD.SaveAs('BkgWeightPlotBu.root')


outputFile = TFile("sPlot2018_outtree.root","RECREATE")
outputFile.cd()

dstree.Write()
outputFile.Close()

