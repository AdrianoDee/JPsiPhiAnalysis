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
mass = RooRealVar("xM","M(#mu#muKK)[GeV]",4.0,6.0)
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


xdataPrompt = (alldata.reduce('xM<4.8')).reduce('xM>4.0').reduce("xL<2.0").reduce("kkM<1.035").reduce("kkM>1.005")
xdataPrompt.numEntries()


# In[ ]:


a0 = RooRealVar("a0","a0",0.001,-10.,10.)
a1 = RooRealVar("a1","a1",0.001,-5.0,5.0)
a2 = RooRealVar("a2","a2",-0.00001,-2.,2.)
a3 = RooRealVar("a3","a3",0.0,-0.5,0.5)
a4 = RooRealVar("a4","a4",0.0,-0.2,0.2)
a5 = RooRealVar("a5","a5",0.0,-0.025,0.05)
a6 = RooRealVar("a6","a6",0.0,-0.001,0.001)

aset = RooArgList(a0,a1,a2,a3,a4,a5)

sigma = RooRealVar("sigma","width of gaussian",0.0013)
gamma = RooRealVar("gamma","gamma of bw",0.004253,0.001,0.01)
mean = RooRealVar("mean","mean of gaussian",phimean,phimean-0.005,phimean+0.005);

nSig = RooRealVar("nSig","nSig",1E6,0.,5.0E6)
nBkg = RooRealVar("nBkg","nBkg",5E5,0.,5.0E6)
cheb = RooChebychev("cheb","Background",masskk,aset)
#gauss = RooGaussian("gauss","gaussian PDF ",mass,mean,sigma)
signal = RooVoigtian("signal","signal",masskk,mean,gamma,sigma)

B_1     = RooRealVar ( "B_1"    , "B_1 "    , 0.3  , -20   , 100   )
B_2     = RooRealVar ( "B_2"    , "B_2"    , 0.3  , -20   , 100   )
B_3     = RooRealVar ( "B_3"    , "B_3"    , 0.3  , -20   , 100   )
B_4     = RooRealVar ( "B_4"    , "B_4"    , 0.3  , -20   , 100   )

bkg    = RooBernstein("pdfB" , "pdfB"    , masskk   , RooArgList(B_1, B_2, B_3, B_4))

tot = RooAddPdf("tot","g+cheb",RooArgList(signal,bkg),RooArgList(nSig,nBkg))


# In[ ]:


#xDataPromptKK = xdataPrompt.reduce(SelectVars(RooArgSet(masskk)))
mean.setConstant(True)
gamma.setConstant(True)
rPhifit = tot.fitTo(xdataPrompt,Range(phimin,phimax),RooFit.NumCPU(4))


#rPhifit = tot.fitTo(xdataPrompt,Range(phimin,phimax),RooFit.NumCPU(4))

mean.setConstant(False)
gamma.setConstant(False)
rPhifit = tot.fitTo(xdataPrompt,Range(phimin,phimax),RooFit.NumCPU(4))

# In[ ]:


c = TCanvas("canvas","canvas",1200,800) 
phiFrame = masskk.frame(Range(phimin,phimax),Normalization((nSig.getValV() + nBkg.getValV())))
xdataPrompt.plotOn(phiFrame)
tot.plotOn(phiFrame)

phiFrame.Draw()
c.SaveAs("phiMassSPlot.png")
c.SaveAs("phiMassSPlot.root")
c.Clear()


# In[ ]:


cD=TCanvas("cD","cD",750,600)
cD.cd()
splot   = RooStats.SPlot ("sPlot","sPlot",xdataPrompt, tot, RooArgList(nSig,nBkg))
dstree  = xdataPrompt.store().tree()



# In[19]:


shist   = TH1F('shist','shist', 500, 4.0, 5.0)


# In[20]:


shist.Sumw2()
shist.SetLineColor(2)    
shist.SetMarkerColor(2); shist.SetMinimum(0.)
dstree.Project('shist','xM','nSig_sw');  


# In[21]:


shist.Draw('e0');
cD.SaveAs('OtherPlotX.gif')
cD.SaveAs('OtherPlotX.root')

