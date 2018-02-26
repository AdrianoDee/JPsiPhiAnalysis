#usage : python2.7 sidebands.py
# change roofilepath to proof output file
#change also phimass histo name

# coding: utf-8

# In[1]:


import ROOT
from ROOT import TFile,TH1,TH1F,TCanvas,TLine
from ROOT import RooFit

# In[3]:


#rootfile ="../rootfiles/X4140_MuMuKK_KRe_MuRe_NP3.0_Alpha99_CW5.15-5.55_PtJPsi7_PtMu_4_NOPHI.root"
rootfile ="../rootfiles/X4140_MuMuKK_KRe_MuRef_Sidebands5-7_NP3.0_B0Cuts_CW5.15-5.55.root"
#histname = "Xcand_histo_hlt8_cw_nonprompt_cosalpha"
#histname = "Phi_hist_all_cw_all_cosalpha"
histname = "PhiMassHisto"
#histname = "Xcand_histo_DM_any_cw_nonprompt_cosalpha"
inputfile  = TFile(rootfile,"READ")
hist = inputfile.Get(histname)
c = TCanvas("canvas","canvas",1200,800) ;


# In[4]:


from ROOT import RooRealVar,RooAbsPdf,RooChebychev,RooExponential,RooGaussian,RooAbsPdf,RooPlot,RooAddPdf,RooDataHist,RooArgSet,RooArgList
from ROOT import kGreen,kRed,kBlack,kBlue,kDashed,kDotted,kMagenta
from ROOT.RooFit import Components,LineColor,LineStyle,Name,Normalization,Layout,Format,Label,Parameters,Range


# In[5]:


massmin = 1.020-0.030
massmax = 1.020+0.030
phimean = 1.020
massbins = (massmax - massmin)/hist.GetBinWidth(2)
mass = RooRealVar("mass","M(KK)[GeV]",massmin,massmax)
mass.setBins(int(massbins))


# In[6]:


hist.GetXaxis().SetRangeUser(massmin,massmax)
hist.Draw()
plotname = "plots/phiplot_" + histname
c.SaveAs(plotname  + ".png")
c.SaveAs(plotname  + ".eps")
c.SaveAs(plotname  + ".pdf")


# In[7]:


mean = RooRealVar("mean","mean of gaussian",phimean,phimean-0.005,phimean+0.005);
sigma = RooRealVar("sigma1","width of gaussian1",0.001,0.0005,0.05);

a0 = RooRealVar("a0","a0",0.001,-1.,1.)
a1 = RooRealVar("a1","a1",0.001,-0.5,0.5)
a2 = RooRealVar("a2","a2",-0.00001,-2.,2.)
a3 = RooRealVar("a3","a3",-0.000001,-0.1,0.1)
a4 = RooRealVar("a4","a4",0.0,-0.1,0.1)
a5 = RooRealVar("a5","a5",0.0,-0.025,0.05)
a6 = RooRealVar("a6","a6",0.0,-0.001,0.001)

aset = RooArgList(a0,a1,a2,a3,a4,a5)
sFrac = RooRealVar("sFrac","sFrac",0.5,0.,1.0)


# In[8]:


cheb = RooChebychev("cheb","Background",mass,aset)
gauss = RooGaussian("gauss","gaussian PDF ",mass,mean,sigma)

tot = RooAddPdf("tot","g+cheb",gauss,cheb,sFrac)


# In[9]:


masslist = RooArgList(mass)
dh = RooDataHist("dh","dh",masslist,hist)
numEvts = dh.sum(False)
print numEvts


# In[10]:


tot.fitTo(dh)


# In[11]:


massFrame = mass.frame()
massFrame.SetTitle("Phi signal")
dh.plotOn(massFrame)
tot.plotOn(massFrame)
gauss.plotOn(massFrame,LineColor(kGreen),LineStyle(kDashed),Normalization((sFrac.getValV()*numEvts)/(numEvts)))
cheb.plotOn(massFrame,LineColor(kMagenta),LineStyle(kDotted),Normalization(((1.0-sFrac.getValV())*numEvts)/(numEvts)))
tot.paramOn(massFrame,Layout(0.60,0.99,0.75));
massFrame.Draw()


# In[12]:


plotmax = hist.GetMaximum()*1.05
sidesigma = sigma.getValV()
leftlowside = -7.*sidesigma + mean.getValV()
leftupside = -5.*sidesigma + mean.getValV()
rightlowside = +5.*sidesigma + mean.getValV()
rightupside = +7.*sidesigma + mean.getValV()

signallow = -3.*sidesigma + mean.getValV()
signalup = +3.*sidesigma + mean.getValV()

print "Side Sigma"
print "%.5e" % sidesigma
print "Left Side"
print "%.5f - %.5f" % (leftlowside,leftupside)
print mean.getValV()


# In[13]:


linelow = TLine(leftlowside,0.0,leftlowside,plotmax)
linemidlow = TLine(leftupside,0.0,leftupside,plotmax)
linemidup = TLine(rightlowside,0.0,rightlowside,plotmax)
lineup = TLine(rightupside,0.0,rightupside,plotmax)

linesiglow = TLine(signallow,0.0,signallow,plotmax)
linesigup = TLine(signalup,0.0,signalup,plotmax)

linelow.SetLineColor(kRed)
linemidlow.SetLineColor(kRed)
linemidup.SetLineColor(kRed)
lineup.SetLineColor(kRed)

linesiglow.SetLineColor(kGreen)
linesigup.SetLineColor(kGreen)

linelow.SetLineWidth(2)
linemidlow.SetLineWidth(2)
linemidup.SetLineWidth(2)
lineup.SetLineWidth(2)

linesiglow.SetLineWidth(2)
linesigup.SetLineWidth(2)

linelow.Draw()
linemidlow.Draw()
linemidup.Draw()
lineup.Draw()

linesiglow.Draw()
linesigup.Draw()


# In[14]:


plotname = "plots/phi_" + histname
c.SaveAs(plotname  + ".png")
c.SaveAs(plotname  + ".eps")
c.SaveAs(plotname  + ".pdf")


# In[15]:


mass.setRange("signalrange",signallow,signalup)
mass.setRange("sideleftrange",leftlowside,leftupside)
mass.setRange("siderightrange",rightlowside,rightupside)
signalIntegralBkg = cheb.analyticalIntegral(cheb.getAnalyticalIntegral(RooArgSet(mass),RooArgSet(mass)),"signalrange")
leftsideIntegralBkg = cheb.analyticalIntegral(cheb.getAnalyticalIntegral(RooArgSet(mass),RooArgSet(mass)),"sideleftrange")
rightsideIntegralBkg = cheb.analyticalIntegral(cheb.getAnalyticalIntegral(RooArgSet(mass),RooArgSet(mass)),"siderightrange")


# In[16]:


totIntegralBkg= cheb.analyticalIntegral(cheb.getAnalyticalIntegral(RooArgSet(mass),RooArgSet(mass)))


# In[17]:


sigBkgEvts = signalIntegralBkg/totIntegralBkg*((1.0-sFrac.getValV())*numEvts)
sidBkgEvts = (leftsideIntegralBkg+rightsideIntegralBkg)/totIntegralBkg*((1.0-sFrac.getValV())*numEvts)


# In[18]:


print sigBkgEvts


# In[19]:


print sidBkgEvts


# In[20]:


ratio = sigBkgEvts/sidBkgEvts


# In[21]:


inputfile.ls()


# In[22]:


signalB0 = inputfile.Get("B0_Cand_Mass_NoM")
sideB0 = inputfile.Get("B0_Cand_Mass_Sides_NoM")
notrebin = True


# In[23]:


if notrebin:
    signalB0.Rebin(5)
    sideB0.Rebin(5)
    notrebin = False


# In[24]:


sideB0.Scale(ratio)


# In[25]:


cB0 = TCanvas("cB0","cB0",1200,800)
signalB0.SetFillColor(kBlue)
signalB0.SetFillStyle(3002)
sideB0.SetFillColor(kRed)
sideB0.SetFillStyle(3002)
signalB0.Draw()
sideB0.Draw("SAMEBar")
cB0.SaveAs("plots/" + signalB0.GetName() + "_sidebands.png")


# In[26]:


b0SideSubtracted = signalB0.Clone()
b0SideSubtracted.Add(sideB0,-1.0)


# In[27]:


b0SideSubtracted.SetFillColor(kGreen)
b0SideSubtracted.SetFillStyle(3002)
b0SideSubtracted.Draw("Bar")
cB0.SaveAs("plots/" + signalB0.GetName() + "_subtracted.png")


# In[ ]:
