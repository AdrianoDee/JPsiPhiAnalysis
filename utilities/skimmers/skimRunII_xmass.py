
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
from ROOT import RooDataSet,RooFormulaVar,RooLinkedList


# In[7]:


no_hlts = 13              


# In[8]:


#rootfile = "/Users/adrianodiflorio/Desktop/mmkk2017/09Jan2017.root"
#rootfile = "/Users/adrianodiflorio/Desktop/mmkk2017/allphi_DataF.root"
rootfile = "/Users/adrianodiflorio/Desktop/mmkk2017/phiJpsiTriggersBCDEF.root"
inputfile = TFile(rootfile,"READ") 
inputfile.ls()
xTupleDir = (inputfile.Get("rootuple")) 
xTupleDir.ls()
#pTuple = (xTupleDir.Get("pTree")) 
xTuple = (xTupleDir.Get("xTree")) 
#jTuple = (xTupleDir.Get("jTree")) 

event = 2

xTuple.SetBranchAddress("event",AddressOf(milk,'price5'))
newfile = TFile("small.root","recreate")
newtree = xTuple.CloneTree()
newtree.CopyEntries(xTuple)

newtree.Write()

sys.exit()

# In[16]:

file = TFile("newFile.root","RECREATE")
canvas = TCanvas("canvas","canvas",1200,1000)
mass = RooRealVar("xM","M(#mu#mu#mu#mu)[GeV]",5.15,5.55)
trigger = RooRealVar("trigger","trigger",0.0,10000)
vProb = RooRealVar("vProb","vProb",-1.0,1.0)
alldata = RooDataSet("alldata","alldata",xTuple,RooArgSet(mass), RooFormulaVar("vProb","vProb","vProb>0.01",RooArgList(vProb)))#,cutFormula)
frame = mass.frame(Range(5.15,5.55))
alldata.plotOn(frame,RooLinkedList())
alldata.Write()
frame.Draw()


# In[ ]:


canvas.SaveAs("testCanvas.eps")

