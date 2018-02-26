
# coding: utf-8

# In[1]:

import ROOT
from ROOT import TFile,TH1,TH1F,TCanvas,TNtuple,TTreeReader,TTreeReaderValue
from ROOT import RooFit
from ROOT.RooFit import Layout
from ROOT import RooStats
from ROOT import RooAbsData
RooAbsData.setDefaultStorageType ( RooAbsData.Tree )
from array import array
import sys


# In[3]:


from ROOT import RooRealVar,RooAbsPdf,RooChebychev,RooExponential,RooGaussian,RooAbsPdf,RooPlot,RooAddPdf,RooDataHist,RooArgSet,RooArgList
from ROOT import kGreen,kRed,kBlack,kBlue,kDashed,kDotted,kMagenta,RooVoigtian
from ROOT.RooFit import Components,LineColor,LineStyle,Name,Normalization,Range,SelectVars
from ROOT import RooDataSet,RooFormulaVar,RooLinkedList,RooBernstein


# In[4]:

rootfile = "../rootfiles/Prompt_skim_cut_JPsi_Phi_Tree_skim_trigger_0_JPsi_Phi_Tree_skim_merge.root" #mmkk 2017 bcdef Jan 18 run
#rootfile = "../rootfiles/NONPrompt_skim_cut_JPsi_Phi_Tree_skim_trigger_0_JPsi_Phi_Tree_skim_merge.root" #mmkk 2017 bcdef Jan 18 run
inputfile = TFile(rootfile,"READ")
inputfile.ls()


# In[7]:


xTree = (inputfile.Get("JPsiPhiCuts Tree"))


# In[30]:


massmin = 4.0
massmax = 6.0
phimin = 1.01
phimax = 1.03


# In[48]:


mass = RooRealVar("xMass","M(#mu#muKK)[GeV]",massmin,massmax)
mass.setBins(500)
masskk = RooRealVar("phiMass","phiMass",phimin,phimax)
masskk.setBins(80)
massmumu = RooRealVar("jpsiMass","jpsiMass",2.9,3.3)


# In[49]:


alldata = RooDataSet("alldata","alldata",xTree,RooArgSet(masskk,mass,massmumu))


# In[50]:


alldata.numEntries()


# In[51]:


alldata.numEntries()
c = TCanvas("canvas","canvas",1200,800)

massFrame = mass.frame()
alldata.plotOn(massFrame)

massFrame.Draw()
c.SaveAs("testmass.png")

kkFrame = masskk.frame()
alldata.plotOn(kkFrame)

kkFrame.Draw()
c.SaveAs("testkk.png")

mumuFrame = massmumu.frame()
alldata.plotOn(mumuFrame)

mumuFrame.Draw()
c.SaveAs("testmm.png")


# In[11]:


#b0dataNonPromptMass = b0dataNonPrompt.reduce(SelectVars(RooArgSet(mass)))


# In[52]:


sigma1 = RooRealVar("sigma1","width of gaussian1",0.002,0.0005,0.05);
sigma2 = RooRealVar("sigma2","width of gaussian2",0.004,0.004,0.01);

phimean = 1.019

sigma = RooRealVar("sigma","width of gaussian",0.0013)
gamma = RooRealVar("gamma","gamma of bw",0.004253,0.001,0.01)
mean = RooRealVar("mean","mean of voigtian",phimean,phimean-0.005,phimean+0.005);

a0 = RooRealVar("a0","a0",0.001,-1.,1.)
a1 = RooRealVar("a1","a1",0.001,-0.5,0.5)
a2 = RooRealVar("a2","a2",-0.00001,-2.,2.)
a3 = RooRealVar("a3","a3",-0.000001,-0.1,0.1)
a4 = RooRealVar("a4","a4",-0.000001,-2.,2.)
a5 = RooRealVar("a5","a5",-0.000001)
a6 = RooRealVar("a6","a6",-0.000001,-0.01,0.01)

aset = RooArgList(a0,a1,a2)#,a3)

B_1     = RooRealVar ( "B_1"    , "B_1 "    , 0.3  , -20   , 100   )
B_2     = RooRealVar ( "B_2"    , "B_2"    , 0.3  , -20   , 100   )
B_3     = RooRealVar ( "B_3"    , "B_3"    , 0.3  , -20   , 100   )
B_4     = RooRealVar ( "B_4"    , "B_4"    , 0.3  , -20   , 100   )

gaussFrac = RooRealVar("sig1frac","fraction of component 1 in signal",0.3,0.0,1.0)
nSig = RooRealVar("nSig","nSig",100000,0.,10E6)
nBkg = RooRealVar("nBkg","nBkg",55000,0.,10E6)


# In[53]:


signal = RooVoigtian("signal","signal",masskk,mean,gamma,sigma)
bkg    = RooBernstein("pdfB" , "pdfB"    , masskk   , RooArgList(B_1, B_2,B_3,B_4))
tot = RooAddPdf("tot","g+cheb",RooArgList(signal,bkg),RooArgList(nSig,nBkg))

#mean.setValV(phimean)
gamma.setConstant(ROOT.kTRUE)
mean.setConstant(ROOT.kTRUE)


# In[ ]:


rfit = tot.fitTo(alldata,Range(phimin,phimax),RooFit.NumCPU(8))
mean.setConstant(ROOT.kFALSE)
rfit = tot.fitTo(alldata,Range(phimin,phimax),RooFit.NumCPU(8))
gamma.setConstant(ROOT.kFALSE)
rfit = tot.fitTo(alldata,Range(phimin,phimax),RooFit.NumCPU(8))



# In[ ]:


kkFrame = masskk.frame(Range(phimin,phimax))
alldata.plotOn(kkFrame,RooLinkedList())
tot.plotOn(kkFrame)

kkFrame.Draw()
c.SaveAs("testmassPhiFit.png")


cD=TCanvas("cD","cD",750,600);cD.cd()
splot   = RooStats.SPlot ( "sPlot","sPlot", alldata, tot, RooArgList(nSig,nBkg))


# In[ ]:


dstree  = alldata.store().tree()
dstree.GetEntryNumber(88)


# In[ ]:


sPlot_B0_hist   = TH1F('sPlot_B0_hist','sPlot_B0_hist', 200, 4.00, 6.0)


# In[ ]:


sPlot_B0_hist.Sumw2()
sPlot_B0_hist.SetLineColor(2)
sPlot_B0_hist.SetMarkerColor(2);
sPlot_B0_hist.SetMinimum(0.)
dstree.Project('sPlot_B0_hist','xMass','nSig_sw');


# In[ ]:


sPlot_B0_hist.Draw('e0');
cD.SaveAs('b0_Splot_Phi.gif')


# In[ ]:


sys.exit()
xdataPrompt = (alldata.reduce('xM<4.8')).reduce('xM>4.0').reduce("xL<2.0")


# In[ ]:


massmin = 1.020-0.03
massmax = 1.020+0.03
phimean = 1.020
xdataPrompt.numEntries()


# In[ ]:


a0 = RooRealVar("a0","a0",0.001,-1.,1.)
a1 = RooRealVar("a1","a1",0.001,-0.5,0.5)
a2 = RooRealVar("a2","a2",-0.00001,-2.,2.)
a3 = RooRealVar("a3","a3",0.0)#
a4 = RooRealVar("a4","a4",0.0,-0.1,0.1)
a5 = RooRealVar("a5","a5",0.0,-0.025,0.05)
a6 = RooRealVar("a6","a6",0.0,-0.001,0.001)

aset = RooArgList(a0,a1,a2)#,a3,a4,a5)

sigma = RooRealVar("sigma","width of gaussian",0.0013)
gamma = RooRealVar("gamma","gamma of bw",0.004253)#,0.001,0.01)
mean = RooRealVar("mean","mean of gaussian",phimean,phimean-0.005,phimean+0.005);

nSig = RooRealVar("nSig","nSig",1E6,0.,5.0E6)
nBkg = RooRealVar("nBkg","nBkg",5E5,0.,5.0E6)
cheb = RooChebychev("cheb","Background",masskk,aset)
#gauss = RooGaussian("gauss","gaussian PDF ",mass,mean,sigma)
signal = RooVoigtian("signal","signal",masskk,mean,gamma,sigma)

tot = RooAddPdf("tot","g+cheb",RooArgList(signal,cheb),RooArgList(nSig,nBkg))


# In[ ]:


masskk.setBins(100)
mass.setBins(100)
h = xdataPrompt.createHistogram(masskk,mass,20,20)
h1 = h.ProjectionX()
dh1 = RooDataHist("kmass","kmass",RooArgList(masskk),h1)


# In[ ]:


rPhifit = tot.fitTo(xdataPrompt,Range(massmin,massmax))


# In[ ]:


c = TCanvas("canvas","canvas",1200,800)
phiFrame = masskk.frame(Range(massmin,massmax))
dh1.plotOn(phiFrame,RooLinkedList())
tot.plotOn(phiFrame)

phiFrame.Draw()
c.SaveAs("testmassFitPhi2.png")
c.Clear()


# In[ ]:
