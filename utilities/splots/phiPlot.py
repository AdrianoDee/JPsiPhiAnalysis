import ROOT
from ROOT import TFile,TH1,TH1F,TCanvas,TNtuple,TTreeReader,TTreeReaderValue
from ROOT import RooFit
from ROOT.RooFit import Layout
from ROOT import RooStats
from ROOT import RooAbsData
RooAbsData.setDefaultStorageType ( RooAbsData.Tree )
from array import array
import sys

from ROOT import RooRealVar,RooAbsPdf,RooChebychev,RooExponential,RooGaussian
from ROOT import RooBernstein,RooAbsPdf,RooPlot,RooAddPdf,RooDataHist,RooArgSet,RooArgList
from ROOT import kGreen,kRed,kBlack,kBlue,kDashed,kDotted,kMagenta,RooVoigtian
from ROOT.RooFit import Components,LineColor,LineStyle,Name,Normalization,Range,SelectVars
from ROOT import RooDataSet,RooFormulaVar,RooLinkedList


# In[4]:


rootfile = "/Users/adrianodiflorio/Desktop/X4140_roots/mumukk_tree.root" 
inputfile = TFile(rootfile,"READ") 
xTuple = (inputfile.Get("outuple")) 


# In[6]:


myReader = TTreeReader("outuple", inputfile)
nentries = xTuple.GetEntries()
print nentries


# In[8]:
phimean = 1.020
phimin = 1.020-0.015
phimax = 1.020+0.015
massbins = (6.0 - 4.0)/0.005
massmin = 5.2
massmax = 5.55
mass = RooRealVar("xM","M(#mu#muKK)[GeV]",5.2,5.55)
mass.setBins(400)
lxy = RooRealVar("xL","l(xy)",0.0,10000.)
hlt = RooRealVar("xHlt","xHlt",0.0,20.0)
masskk = RooRealVar("kkM","kkM",phimin,phimax);
masskk.setBins(int(200))
massmumu = RooRealVar("mumuM","mumuM",2.5,3.5)


# In[9]:


alldata = RooDataSet("alldata","alldata",xTuple,RooArgSet(masskk,mass,lxy,hlt,massmumu))#,cutFormula)


# In[10]:


alldata.numEntries()


# In[ ]:


xdataPrompt = (alldata.reduce('xM>5.15')).reduce('xM<5.55').reduce("xL>=3.0").reduce("kkM<1.035").reduce("kkM>1.005")
xdataPrompt.numEntries()

B_c     = RooRealVar ( "B_c"    , "B_c "    , 0.3  , -20   , 100   )
B_b     = RooRealVar ( "B_b"    , "B_b "    , 0.3  , -20   , 100   )

S1      = RooRealVar ( "S1"     , "Signal"  , 30000 , 1     , 900000    )
S2      = RooRealVar ( "S2"     , "Signal"  , 30000 , 1     , 900000    )
B       = RooRealVar ( "B"      , "B"       , 50000 , 1     , 900000000 )

mean = RooRealVar("mean","mean of gaussian",5.38,5.31,5.41);
sigma1 = RooRealVar("sigma1","width of gaussian1",0.002,0.0005,0.05);
sigma2 = RooRealVar("sigma2","width of gaussian2",0.004,0.004,0.01);

pdfS1   = RooGaussian( "pdfS1"  , "gaus"    , mass   ,mean, sigma1)
pdfS2   = RooGaussian( "pdfS2"  , "gaus"    , mass   ,mean, sigma2)
# pdfB    = RooExponential("pdfB" , "pdfB"    , mbs   , B_c)
pdfB    = RooBernstein("pdfB" , "pdfB"    , mass   , RooArgList(B_c, B_b))

alist1  = RooArgList (pdfS1, pdfS2, pdfB);  alist2 = RooArgList (S1, S2, B);
tot  = RooAddPdf  ("model", "model", alist1, alist2)


# In[ ]:


#xDataPromptKK = xdataPrompt.reduce(SelectVars(RooArgSet(masskk)))
mean.setConstant(True)
rPhifit = tot.fitTo(xdataPrompt,Range(massmin,massmax),RooFit.NumCPU(8))

mean.setConstant(False)
rPhifit = tot.fitTo(xdataPrompt,Range(massmin,massmax),RooFit.NumCPU(8))

# In[ ]:


c = TCanvas("canvas","canvas",1200,800) 
phiFrame = masskk.frame(Range(massmin,massmax))
xdataPrompt.plotOn(phiFrame)
tot.plotOn(phiFrame)

phiFrame.Draw()
c.SaveAs("phiMassSPlotPhi.png")
c.SaveAs("phiMassSPlotPhi.root")
c.Clear()


# In[ ]:


cD=TCanvas("cD","cD",750,600)
cD.cd()
splot   = RooStats.SPlot ("sPlot","sPlot",xdataPrompt, tot, alist2)
dstree  = xdataPrompt.store().tree()



# In[19]:


shist   = TH1F('shist','shist', 200, 0.97, 1.07)


# In[20]:


shist.Sumw2()
shist.SetLineColor(2)    
shist.SetMarkerColor(2); shist.SetMinimum(0.)
dstree.Project('shist','kkM','S1_sw + S2_sw');  


# In[21]:


shist.Draw('e0');
cD.SaveAs('OtherPlotPhi.gif')
cD.SaveAs('OtherPlotPhi.root')


