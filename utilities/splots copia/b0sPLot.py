import ROOT
from ROOT import TFile,TH1,TH1F,TCanvas,TNtuple,TTreeReader,TTreeReaderValue
from ROOT import RooFit
from ROOT.RooFit import Layout
from ROOT import RooStats
from ROOT import RooAbsData,RooBernstein
RooAbsData.setDefaultStorageType ( RooAbsData.Tree )
from array import array
import sys


# In[3]:


from ROOT import RooRealVar,RooAbsPdf,RooChebychev,RooExponential,RooGaussian,RooAbsPdf,RooPlot,RooAddPdf,RooDataHist,RooArgSet,RooArgList
from ROOT import kGreen,kRed,kBlack,kBlue,kDashed,kDotted,kMagenta,RooVoigtian
from ROOT.RooFit import Components,LineColor,LineStyle,Name,Normalization,Range,SelectVars
from ROOT import RooDataSet,RooFormulaVar,RooLinkedList


# In[4]:


rootfile = "/Users/adrianodiflorio/Desktop/X4140_roots/mumukk_tree.root" 
inputfile = TFile(rootfile,"READ") 
xTuple = (inputfile.Get("outuple")) 


# In[5]:


massmin = 1.2
massmax = 5.55


# In[6]:


myReader = TTreeReader("outuple", inputfile)
nentries = xTuple.GetEntries()
print nentries


# In[8]:


massbins = (6.0 - 4.0)/0.005
mass = RooRealVar("xM","M(#mu#muKK)[GeV]",4.0,6.0)
mass.setBins(int(massbins))
lxy = RooRealVar("xL","l(xy)",0.0,10000.)
hlt = RooRealVar("xHlt","xHlt",0.0,20.0)
masskk = RooRealVar("kkM","kkM",0.5,1.5)
massbins = 100
masskk.setBins(int(massbins))
massmumu = RooRealVar("mumuM","mumuM",2.5,3.5)
cutFormula = RooFormulaVar("cutFormula","cutFormula","xHlt!=8.0",RooArgList(hlt))


# In[9]:


alldata = RooDataSet("alldata","alldata",xTuple,RooArgSet(masskk,mass,lxy,hlt,massmumu))#,cutFormula)
datasetfile = TFile("xMassDataset.root","RECREATE") 
datasetfile.cd()
alldata.Write()


# In[10]:


alldata.numEntries()


# In[ ]:


#xb->setRange("alt","x_coarse_bin1,x_coarse_bin3,x_coarse_bin5,x_coarse_bin7,x_coarse_bin9") ;
b0dataNonPrompt = ((alldata.reduce('xHlt!=8')).reduce('xM>5.2')).reduce("xL>3.0")


# In[ ]:


b0dataNonPromptMass = b0dataNonPrompt.reduce(SelectVars(RooArgSet(mass)))
b0dataNonPrompt.numEntries()


# In[ ]:


c = TCanvas("canvas","canvas",1200,800) 
mass.setRange("fitRange",massmin,massmax)
mass.setBins(200)
massFrame = mass.frame(Range(massmin,massmax))
b0dataNonPrompt.plotOn(massFrame)
massFrame.Draw()
c.SaveAs("testmass.png")


# In[13]:

phimean = 1.019

sigma = RooRealVar("sigma","width of gaussian",0.0013)
gamma = RooRealVar("gamma","gamma of bw",0.004253,0.001,0.01)
mean = RooRealVar("mean","mean of gaussian",phimean,phimean-0.005,phimean+0.005);

nSig = RooRealVar("nSig","nSig",1E6,0.,5.0E6)
nBkg = RooRealVar("nBkg","nBkg",5E5,0.,5.0E6)
#cheb = RooChebychev("cheb","Background",masskk,aset)
#gauss = RooGaussian("gauss","gaussian PDF ",mass,mean,sigma)
signal = RooVoigtian("signal","signal",masskk,mean,gamma,sigma)

B_1     = RooRealVar ( "B_1"    , "B_1 "    , 0.3  , -20   , 100   )
B_2     = RooRealVar ( "B_2"    , "B_2"    , 0.3  , -20   , 100   )
B_3     = RooRealVar ( "B_3"    , "B_3"    , 0.3  , -20   , 100   )
B_4     = RooRealVar ( "B_4"    , "B_4"    , 0.3  , -20   , 100   )

bkg    = RooBernstein("pdfB" , "pdfB"    , masskk   , RooArgList(B_1, B_2,B_3,B_4))

tot = RooAddPdf("tot","g+cheb",RooArgList(signal,bkg),RooArgList(nSig,nBkg))

mean.setVal(phimean)
gamma.setConstant(ROOT.kTRUE)
mean.setConstant(ROOT.kTRUE)

rfit = tot.fitTo(b0dataNonPrompt,Range(massmin,massmax),RooFit.NumCPU(8))
mean.setConstant(ROOT.kFALSE)
rfit = tot.fitTo(b0dataNonPrompt,Range(massmin,massmax),RooFit.NumCPU(8))
gamma.setConstant(ROOT.kFALSE)
rfit = tot.fitTo(b0dataNonPrompt,Range(phimean-0.025,phimean+0.025),RooFit.NumCPU(8))

masskkFrame = masskk.frame(Range(phimean-0.025,phimean+0.025))
b0dataNonPrompt.plotOn(massFrame,RooLinkedList())
tot.plotOn(masskkFrame)

massFrame.Draw()
c.SaveAs("testmassFit.png")

cD=TCanvas("cD","cD",750,600)
cD.cd()
splot   = RooStats.SPlot ( "sPlot","sPlot", b0dataNonPrompt, tot, RooArgList(nSig,nBkg))


# In[18]:


dstree  = b0dataNonPrompt.store().tree()


# In[19]:


shist   = TH1F('shist','shist', 500, 1.00, 1.05)


# In[20]:


shist.Sumw2()
shist.SetLineColor(2)    
shist.SetMarkerColor(2); shist.SetMinimum(0.)
dstree.Project('shist','kkM','nSig_sw');  


# In[21]:


shist.Draw('e0');
cD.SaveAs('OtherPlot.gif')



