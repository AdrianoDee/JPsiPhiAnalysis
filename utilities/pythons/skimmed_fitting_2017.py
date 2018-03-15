
import ROOT

from ROOT import TLine, TLegend, TPad
from ROOT import TFile,TH1,TH1F,TCanvas,TNtuple,TTreeReader,TTreeReaderValue
from ROOT import RooFit
from ROOT.RooFit import Layout, Components, Title, Name, Normalization, Layout, Format, Label, Parameters, Range, LineColor, LineStyle
from ROOT import RooStats, gPad, RooAbsData, RooAbsReal, RooBinning
from ROOT.RooAbsReal import Relative
from ROOT import kGreen,kRed,kBlack,kBlue,kDashed,kDotted,kMagenta

RooAbsData.setDefaultStorageType ( RooAbsData.Tree )
from array import array
import sys

import argparse
import math
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--debug', action='store_true')
parser.add_argument('--input', type=str, default="2mu2k_tree.root")
parser.add_argument('--binned', action='store_true')
parser.add_argument('--nonprompt', action='store_true')
parser.add_argument('--prompt', action='store_true')
parser.add_argument('--nofit', action='store_false')
parser.add_argument('--nofitkk', action='store_false')
parser.add_argument('--nofitb0', action='store_false')
parser.add_argument('--nofitmm', action='store_false')
parser.add_argument('--noplot', action='store_false')
parser.add_argument('--numcpu',  type=int, default=4)
parser.add_argument('--nsignal', type=float, default=100000.)
parser.add_argument('--ptcuts', type=float, default=None)
parser.add_argument('--sigmaside', type=float, default=0.001)
parser.add_argument('--signalside', type=float, default=3.0)
parser.add_argument('--sidelow', type=float, default=4.0)
parser.add_argument('--sidehigh', type=float, default=6.0)
parser.add_argument('--binwise', type=int, default=None)
#                    help='number of epochs')
#parser.add_argument('--batch_size', type=int, default=64)
#
args = parser.parse_args()

debugging       = args.debug
binnedfit       = args.binned
numcpus         = args.numcpu
sigmaside_kk    = args.sigmaside
signalside      = args.signalside
sidelow         = args.sidelow
sidehigh        = args.sidehigh
region          = "_overall_"
cuts            = "_"

if(args.prompt and args.nonprompt):
    print("Exiting. Choose prompt or nonprompt.")
    sys.exit()
else:
    if args.nonprompt:
        region = "_nprompt_"
    if args.prompt:
        region = "_promt_"

from ROOT import RooRealVar,RooAbsPdf,RooChebychev,RooExponential,RooGaussian,RooAbsPdf,RooPlot,RooAddPdf,RooDataHist,RooArgSet,RooArgList
from ROOT import kGreen,kRed,kBlack,kBlue,kDashed,kDotted,kMagenta,RooVoigtian,RooArgSet
from ROOT.RooFit import Components,LineColor,LineStyle,Name,Normalization,Range,SelectVars
from ROOT import RooDataSet,RooFormulaVar,RooLinkedList,RooBernstein
from ROOT.RooFit import PrintLevel

# In[4]:

filename = args.input

file_2Trk = TFile(filename,"READ")
file_2Trk.ls()
theTree = (file_2Trk.Get("outuple"))

phimin = 1.0
phimax = 1.05

jpsimin = 2.85
jpsimax = 3.35

xmin = 4.0
xmax = 6.0

bZeromin = 5.15
bZeromax = 5.55

ptmin = 0.0
ptmax = 10000.0

tt_mass = RooRealVar("ttM","ttM;m(KK)[GeV]",phimin,phimax)
mm_mass = RooRealVar("mmM","mmM;m(#mu#mu)[GeV]",jpsimin,jpsimax)
mmtt_mass = RooRealVar("xM","xM;m(#mu#muKK)[GeV]",xmin,xmax)

tt_pt = RooRealVar("ttPt","ttPt;p_{t}(KK)[GeV]",0.0,ptmax)
mm_pt = RooRealVar("mmPt","mmPt;p{t}(#mu#mu)[GeV]",0.0,ptmax)
mmtt_pt = RooRealVar("xPt","xPt;p_{t}(#mu#muKK)[GeV]",0.0,ptmax)

lxysig = RooRealVar("xL","l_{xy} sign.;l_{xy} / #sigma(l_{xy})",0.0,1000.0)

massvars    = RooArgSet(tt_mass, mm_mass, mmtt_mass)
ptvars      = RooArgSet(tt_pt, mm_pt, mmtt_pt)
kinvars     = RooArgSet(massvars,ptvars)
extravars   = RooArgSet(lxysig)

theData = RooDataSet("theData","theData",theTree,RooArgSet(kinvars,extravars))


c = TCanvas("canvas","canvas",1200,800)

if args.nonprompt:
    theData = theData.reduce("xL > 3.0")
    region = "_nprompt_"
if args.prompt:
    theData = theData.reduce("xL < 1.5")
    region = "_promt_"

if args.ptcuts is not None:
    theData = theData.reduce("trigp_pT > " + str(args.ptcuts))
    theData = theData.reduce("trign_pT > " + str(args.ptcuts))
    cuts += "P_t_" + str(args.ptcuts) + "_"
#### #### Plotting variables
#### TrakTrak Data

if args.noplot:

    print("TrakTrak data plotting . . .")

    print("All : " + str(theData.numEntries()))
    ttFrame = tt_mass.frame(Title("KK mass"))
    theData.plotOn(ttFrame)
    ttFrame.Draw()
    c.SaveAs("tt_mass" + cuts + ".png")

    #### MuMu Data
    print("MuMu data plotting . . .")

    print("All : " + str(theData.numEntries()))
    mumuFrame = mm_mass.frame(Title("KK mass"))
    theData.plotOn(mumuFrame)
    mumuFrame.Draw()
    c.SaveAs("mm_mass.png")

    #### X Data
    print("TrakTrakMuMu data plotting . . .")

    print("All : " + str(theData.numEntries()))
    mmttFrame = mmtt_mass.frame(Title("#mu#muKK mass"))
    theData.plotOn(mmttFrame)
    mmttFrame.Draw()
    c.SaveAs("mmtt_mass" + region + ".png")

    #### X lxy Data
    print("TrakTrakMuMu data plotting . . .")

    print("All : " + str(theData.numEntries()))
    lxyFrame = lxysig.frame()
    theData.plotOn(lxyFrame)
    lxyFrame.Draw()
    c.SaveAs("lxy" + region + ".png")


    ##Plotting

    hist_mmtt      = (theData.createHistogram(mmtt_mass,mmtt_mass)).ProjectionX("hist_mmtt_mass")
    hist_mm        = (theData.createHistogram(mm_mass,mmtt_mass)).ProjectionX("hist_mm_mass")
    hist_tt        = (theData.createHistogram(tt_mass,mmtt_mass)).ProjectionX("hist_tt_mass")
    hist_tt_pt     = (theData.createHistogram(mm_pt,mmtt_mass,40,100)).ProjectionX("hist_mm_pt")
    hist_mm_pt     = (theData.createHistogram(tt_pt,mmtt_mass,180,100)).ProjectionX("hist_tt_pt")
    hist_mmtt_pt   = (theData.createHistogram(mmtt_pt,mmtt_mass,180,100)).ProjectionX("hist_mmtt_pt")

    ##Masses hists

    hist_mmtt.SetTitle("M(#mu#muKK)" + region + " ;M(#mu#muKK) [GeV]; no. entries")
    hist_mmtt.SetFillColor(ROOT.kRed)
    hist_mmtt.SetLineColor(ROOT.kRed)
    hist_mmtt.SetFillStyle(3002)
    hist_mmtt.Draw()
    c.SaveAs("mmtt_mass_histo" + region + ".png")
    c.SaveAs("mmtt_mass_histo" + region + ".root")

    hist_mm.SetTitle("M(#mu#mu)" + region + " ;M(#mu#mu) [GeV]; no. entries")
    hist_mm.SetFillColor(ROOT.kGreen)
    hist_mm.SetLineColor(ROOT.kGreen)
    hist_mm.SetFillStyle(3002)
    hist_mm.Draw()
    c.SaveAs("mm_mass_histo" + region + ".png")
    c.SaveAs("mm_mass_histo" + region + ".root")

    hist_tt.SetTitle("M(KK)" + region + " ;M(KK) [GeV]; no. entries")
    hist_tt.GetXaxis().SetRangeUser(1.0,1.05)
    hist_tt.SetFillColor(ROOT.kBlue)
    hist_tt.SetLineColor(ROOT.kBlue)
    hist_tt.SetFillStyle(3002)
    hist_tt.Draw()
    c.SaveAs("tt_mass_histo" + region + ".png")
    c.SaveAs("tt_mass_histo" + region + ".root")

    ##PT hists

    hist_tt_pt.SetTitle("p_t(KK)" + region + " ;p_t(KK) [GeV]; no. entries")
    hist_tt_pt.GetXaxis().SetRangeUser(0.0,20.0)
    hist_tt_pt.SetFillColor(ROOT.kBlue)
    hist_tt_pt.SetLineColor(ROOT.kBlue)
    hist_tt_pt.SetFillStyle(3003)
    hist_tt_pt.Draw()
    c.SaveAs("tt_pt_histo" + region + ".png")
    c.SaveAs("tt_pt_histo" + region + ".root")

    hist_mm_pt.SetTitle("p_t(#mu#mu)" + region + " ;p_t(#mu#mu) [GeV]; no. entries")
    hist_mm_pt.GetXaxis().SetRangeUser(0.0,90.0)
    hist_mm_pt.SetFillColor(ROOT.kGreen)
    hist_mm_pt.SetLineColor(ROOT.kGreen)
    hist_mm_pt.SetFillStyle(3003)
    hist_mm_pt.Draw()
    c.SaveAs("mm_pt_histo" + region + ".png")
    c.SaveAs("mm_pt_histo" + region + ".root")

    hist_mmtt_pt.SetTitle("p_t(#mu#muKK)" + region + " ;p_t(#mu#muKK) [GeV]; no. entries")
    hist_mmtt_pt.GetXaxis().SetRangeUser(0.0,90.0)
    hist_mmtt_pt.SetFillColor(ROOT.kRed)
    hist_mmtt_pt.SetLineColor(ROOT.kRed)
    hist_mmtt_pt.SetFillStyle(3003)
    hist_mmtt_pt.Draw()
    c.SaveAs("mmtt_pt_histo" + region + ".png")
    c.SaveAs("mmtt_pt_histo" + region + ".root")

if args.nofit and args.nofitkk:

    ## Phi fitting

    kcanvas = TCanvas("kcanvas","kcanvas",1400,800)

    pullpad = TPad("pullpad","pullpad",0.0,0.05,1.0,0.33)
    plotpad = TPad("histopad","histopad",0.0,0.35,1.0,1.0)
    plotpad.SetFillStyle(4004)
    pullpad.SetFillStyle(4004)
    plotpad.Draw()
    pullpad.Draw()

    traKFitData = theData.Clone("fitTrakData")

    if binnedfit:
        tt_mass.setBins(30)
        traKFitData = theData.binnedClone("binnedTrakData")

    phimean = 1.019
    gammavalue = 0.0012

    fitphimin = 1.01
    fitphimax = 1.03

    kkSigma = RooRealVar("#sigma","#sigma",0.0013)
    kkGamma = RooRealVar("#Gamma","#Gamma",gammavalue,0.001,0.015)
    kkMean = RooRealVar("m_{kk}","m_{kk}",phimean,phimean-0.005,phimean+0.005);

    # B_1     = RooRealVar ( "B_1"    , "B_1 "    , 0.3  , -20   , 100   )
    # B_2     = RooRealVar ( "B_2"    , "B_2"    , 0.3  , -20   , 100   )
    # B_3     = RooRealVar ( "B_3"    , "B_3"    , 0.3  , -20   , 100   )
    # B_4     = RooRealVar ( "B_4"    , "B_4"    , 0.3  , -20   , 100   )
    # B_5     = RooRealVar ( "B_5"    , "B_5"    , 0.3  , -20   , 100   )
    # B_6     = RooRealVar ( "B_6"    , "B_6"    , 0.3  , -20   , 100   )

    a0 = RooRealVar("p0","p0",0.001,-10.,10.)
    a1 = RooRealVar("p1","p1",0.001,-10.,10.)
    a2 = RooRealVar("p2","p2",-0.00001,-10.,10.)
    a3 = RooRealVar("p3","p3",-0.000001,-10.,10.)
    a4 = RooRealVar("p4","p4",-0.000001,-10.,10.)
    a5 = RooRealVar("a5","a5",-0.000001,-10.,10.)
    poliset = RooArgList(a0,a1,a2,a3,a4)

    # gaussFrac = RooRealVar("s","fraction of component 1 in kkSig",0.3,0.0,1.0)
    nSigKK = RooRealVar("nSig","nSig",theData.numEntries()*0.3,0.0,theData.numEntries()*1.5)
    nBkgKK = RooRealVar("nBkg","nBkg",theData.numEntries()*0.7,0.0,theData.numEntries()*1.5)

    kkSig = RooVoigtian("kkSig","kkSig",tt_mass,kkMean,kkGamma,kkSigma)
    #kkSig = RooGaussian("kkSig","kkSig",tt_mass,kkMean,kkGamma)#,kkSigma)
    #kkBkg = RooBernstein("kkBkg" , "kkBkg", tt_mass, RooArgList(B_1, B_2,B_3,B_4))#,B_5) )#,B_6))
    kkBkg = RooChebychev("kkBkg","Background",tt_mass,poliset)
    kkTot = RooAddPdf("kkTot","kkTot",RooArgList(kkSig,kkBkg),RooArgList(nSigKK,nBkgKK))

    kkGamma.setConstant(ROOT.kTRUE)
    kkMean.setConstant(ROOT.kTRUE)

    nfit = 0

    #kkfit = kkTot.fitTo(traKFitData,Range(fitphimin+0.005,fitphimax-0.005),RooFit.PrintLevel(-1), RooFit.NumCPU(7),RooFit.Save())
    #nfit +=1

    if not debugging:

        kkMean.setConstant(ROOT.kFALSE)
        # kkfit = kkTot.fitTo(traKFitData,Range(fitphimin+0.005,fitphimax-0.005),RooFit.PrintLevel(-1), RooFit.NumCPU(7),RooFit.Save())
        # nfit +=1

        kkGamma.setConstant(ROOT.kFALSE)
        # kkfit = kkTot.fitTo(traKFitData,Range(fitphimin+0.005,fitphimax-0.005),RooFit.PrintLevel(-1), RooFit.NumCPU(7),RooFit.Save())
        # nfit +=1

        #kkfit = kkTot.fitTo(traKFitData,Range(fitphimin+0.005,fitphimax-0.005), RooFit.NumCPU(7),RooFit.Save())
        kkfit = kkTot.fitTo(traKFitData,Range(fitphimin,fitphimax),RooFit.PrintLevel(-1), RooFit.NumCPU(numcpus),RooFit.Save())
        nfit +=1
    else:
        nfit +=1

    sigmaside_kk = math.sqrt(kkGamma.getValV()**2 + kkSigma.getValV()**2)
    sigmaside_kk = kkGamma.getValV()
    if debugging:
	       sigmaside_kk = 0.001

    leftlowside = -signalside*sigmaside_kk + kkMean.getValV()
    leftupside = -sidelow*sigmaside_kk + kkMean.getValV()
    rightlowside = +sidehigh*sigmaside_kk + kkMean.getValV()
    rightupside = +sidehigh*sigmaside_kk + kkMean.getValV()

    signallow = -signalside*sigmaside_kk + kkMean.getValV()
    signalup = +signalside*sigmaside_kk + kkMean.getValV()

    tt_mass.setRange("signalrange",signallow,signalup)
    tt_mass.setRange("sideleftrange",leftlowside,leftupside)
    tt_mass.setRange("siderightrange",rightlowside,rightupside)
    signalIntegralBkg = kkBkg.analyticalIntegral(kkBkg.getAnalyticalIntegral(RooArgSet(tt_mass),RooArgSet(tt_mass)),"signalrange")
    leftsideIntegralBkg = kkBkg.analyticalIntegral(kkBkg.getAnalyticalIntegral(RooArgSet(tt_mass),RooArgSet(tt_mass)),"sideleftrange")
    rightsideIntegralBkg = kkBkg.analyticalIntegral(kkBkg.getAnalyticalIntegral(RooArgSet(tt_mass),RooArgSet(tt_mass)),"siderightrange")

    totIntegralBkg = kkBkg.analyticalIntegral(kkBkg.getAnalyticalIntegral(RooArgSet(tt_mass),RooArgSet(tt_mass)))

    sigBkgEvts = signalIntegralBkg/totIntegralBkg*((nBkgKK.getValV()))
    sidBkgEvts = (leftsideIntegralBkg+rightsideIntegralBkg)/totIntegralBkg*((nBkgKK.getValV()))

    leftsidedata = theData.reduce("ttM<" + str(leftupside))
    leftsidedata = leftsidedata.reduce("ttM>" + str(leftlowside))

    rigthsidedata = theData.reduce("ttM<" + str(rightupside))
    rigthsidedata = rigthsidedata.reduce("ttM>" + str(rightlowside))

    signaldata = theData.reduce("ttM<" + str(signalup))
    signaldata = signaldata.reduce("ttM>" + str(signallow))

    signalhist    = (signaldata.createHistogram(mmtt_mass,mmtt_mass,200,100)).ProjectionX("hist_mmtt_mass_signal")
    leftsidehist  = (leftsidedata.createHistogram(mmtt_mass,mmtt_mass,200,100)).ProjectionX("hist_mmtt_mass_left")
    rightsidehist = (rigthsidedata.createHistogram(mmtt_mass,mmtt_mass,200,100)).ProjectionX("hist_mmtt_mass_right")
    theRatio = sigBkgEvts/sidBkgEvts

    print(theRatio)

    print(signaldata.numEntries())
    print(rigthsidedata.numEntries())
    leftsidehist.Scale(theRatio)
    rightsidehist.Scale(theRatio)

    sidehist = leftsidehist.Clone()
    sidehist.Add(rightsidehist)

    plotpad.cd()

    kkFrame = tt_mass.frame(Range(fitphimin,fitphimax),Title("#phi mass"))

    # kkbins = RooBinning(-15,15)
    # kkbins.addUniform(30,fitphimin,fitphimax)
    traKFitData.plotOn(kkFrame,Name("Data"))
    kkTot.plotOn(kkFrame,RooFit.Normalization(1.0/float(nfit)),Name("Pdf"))
    traKFitData.plotOn(kkFrame)
    kkTot.paramOn(kkFrame,RooFit.Layout(0.57,0.99,0.65))

    kkFrame.Draw()

    pullpad.cd()
    hpull = bZeroFrame.pullHist("Data","Pdf")
    pullframe = mmtt_mass.frame(Title("Pull Distribution"),Range(fitphimin,fitphimax))
    pullframe.SetAxisRange (-5.0,5.0,"Y")
    #pullframe.GetXaxis().SetTitleSize(0.04)
    #pullframe.GetYaxis().SetTitleSize(0.03)
    ROOT.gStyle.SetTitleFontSize(0.07)
    pullframe.addPlotable(hpull,"P")
    pullframe.Draw()

    lineup = TLine(fitphimin,4.0,fitphimax,4.0)
    lineup.SetLineColor(kRed)
    lineup.SetLineStyle(kDashed)
    lineup.Draw()

    linedw = TLine(fitphimin,-4.0,fitphimax,-4.0)
    linedw.SetLineColor(kRed)
    linedw.SetLineStyle(kDashed)
    linedw.Draw()

    kcanvas.SaveAs("kk_Phi_fit" + region + ".png")
    kcanvas.SaveAs("kk_Phi_fit" + region + ".root")
    kcanvas.Clear()

    signalhist.SetFillColor(kBlue)
    signalhist.SetName("B_{0}^{s} Candidates Mass - CC")
    signalhist.SetTitle("B_{0}^{s} Candidates Mass - CC; M(KK#mu#mu)[GeV]; candidates/" + str(signalhist.GetBinWidth(2)*1000)+ "MeV")
    signalhist.GetYaxis().SetTitleOffset(1.3)
    signalhist.SetMarkerColor(kBlue)
    signalhist.SetFillStyle(3002)
    signalhist.SetMarkerStyle(ROOT.kFullCircle)
    signalhist.SetMarkerSize(0.5)
    signalhist.SetLineColor(kBlue)

    leftsidehist.SetFillColor(kRed)
    #leftsidehist.SetMarkerColor(kBlack)
    leftsidehist.SetFillStyle(3002)
    leftsidehist.SetMarkerStyle(ROOT.kFullCircle)
    leftsidehist.SetMarkerSize(0.5)
    leftsidehist.SetLineColor(kBlack)

    rightsidehist.SetFillColor(kGreen)
    #rightsidehist.SetMarkerColor(kBlack)
    rightsidehist.SetMarkerStyle(ROOT.kFullCircle)
    rightsidehist.SetMarkerSize(0.5)
    rightsidehist.SetLineColor(kBlack)
    rightsidehist.SetFillStyle(3002)

    signalhist.Draw("EBar")
    sidehist.SetFillColor(kGreen)
    #rightsidehist.SetMarkerColor(kBlack)
    sidehist.SetMarkerStyle(ROOT.kFullCircle)
    sidehist.SetMarkerSize(0.5)
    sidehist.SetLineColor(kBlack)
    sidehist.SetFillStyle(3002)
    #sideCW = leftsidehist.Clone()
    #sideCW.Add(rightsidehist,+1.0)
    #sideCW.SetFillColor(kRed)
    #sideCW.SetFillStyle(3002)
    #sideCW.Scale(2.0)
    #sideCW.Draw("SAMEBar")
    #rightsidehist.Draw("E0SAMEBar")
    #leftsidehist.Draw("E0SAMEBar")
    sidehist.Draw("E0SAMEBar")

    legend = TLegend(0.75,0.45,0.99,0.75)
    legend.AddEntry(signalhist,"Signal region","f")
    #legend.AddEntry(signalhist,"Signal region (-3.0#sigma,+3.0#sigma)","f")
    #legend.AddEntry(rightsidehist,"R-sideband    (+4.0#sigma,+6.0#sigma)","f")
    #legend.AddEntry(leftsidehist,"L-sideband    (-6.0#sigma,-4.0#sigma)","f")
    legend.AddEntry(sidehist,"Sidebands ","f")
    legend.Draw()
    kcanvas.SaveAs(signalhist.GetName() + "_sidebands" + region + ".png")
    kcanvas.SaveAs(signalhist.GetName() + "_sidebands" + region + ".root")

    ROOT.gStyle.SetOptStat(0)
    b0SideSub = signalhist.Clone()
    b0SideSub.SetTitle("B_{0}^{s} Mass - #phi sides subtracted; M(KK#mu#mu)[GeV]; candidates/" + str(signalhist.GetBinWidth(2)*1000)+ "MeV")
    b0SideSub.Add(leftsidehist,-1.0)
    b0SideSub.Add(rightsidehist,-1.0)

    b0SideSub.SetFillColor(4004)
    b0SideSub.SetFillStyle(3002)
    b0SideSub.SetMarkerStyle(ROOT.kFullCircle)
    b0SideSub.SetMarkerColor(kBlack)
    b0SideSub.SetMarkerSize(0.8)
    b0SideSub.SetLineColor(kBlack)
    b0SideSub.Draw("E0")

    linezero = TLine(b0SideSub.GetBinCenter(1),0.0,b0SideSub.GetBinCenter(b0SideSub.GetNbinsX()),0.0)
    linezero.SetLineColor(kRed)
    linezero.SetLineWidth(2)
    linezero.SetLineStyle(kDotted)
    linezero.Draw()
    kcanvas.SaveAs(signalhist.GetName() + "_subtracted" + region + ".png")
    kcanvas.SaveAs(signalhist.GetName() + "_subtracted" + region + ".root")


    splot   = RooStats.SPlot ( "sPlot","sPlot", theData, kkTot, RooArgList(nSigKK,nBkgKK))

    dstree  = theData.store().tree()
    dstree.GetEntryNumber(88)

    sPlot_B0_hist   = TH1F('sPlot_B0_hist','sPlot_B0_hist', 2000, 4.00, 6.0)

    sPlot_B0_hist.Sumw2()
    sPlot_B0_hist.SetLineColor(2)
    sPlot_B0_hist.SetMarkerColor(2);
    sPlot_B0_hist.SetMinimum(0.)
    dstree.Project('sPlot_B0_hist','xM','nSig_sw');

    sPlot_B0_hist.Draw('e0');
    kcanvas.SaveAs('b0_Splot_Phi' + region + '.png')
    kcanvas.SaveAs('b0_Splot_Phi' + region + '.root')

if args.nofit and args.nofitb0:

    bcanvas = TCanvas("bcanvas","bcanvas",1200,800)

    pullpad = TPad("pullpad","pullpad",0.0,0.05,1.0,0.33)
    plotpad = TPad("histopad","histopad",0.0,0.35,1.0,1.0)
    plotpad.SetFillStyle(4004)
    pullpad.SetFillStyle(4004)
    plotpad.Draw()
    pullpad.Draw()

    bZeroFitData = theData.Clone("fitB0Data")

    fitbZeromin = 5.15
    fitbZeromax = 5.55

    if binnedfit:
        mmtt_mass.setBins(50)
        bZeroFitData = theData.binnedClone("binnedTrakData")

    bZeroFitData = (bZeroFitData.reduce("xM<5.55")).reduce("xM>5.15")

    mean = RooRealVar("m_{#mu#mukk}","m_{#mu#mukk}",5.35,5.2,5.4);
    sigma1 = RooRealVar("#sigma_{1}","#sigma_{1}",0.002,0.0005,0.1);
    sigma2 = RooRealVar("#sigma_{2}","#sigma_{2}",0.004,0.0005,0.1);

    c0 = RooRealVar("p0","p0",0.001,-10.,10.)
    c1 = RooRealVar("p1","p1",0.001,-10.,10.)
    c2 = RooRealVar("p2","p2",-0.00001,-10.,10.)
    c3 = RooRealVar("p3","p3",-0.000001,-10.,10.)
    c4 = RooRealVar("p4","p4",-0.000001,-10.,10.)
    c5 = RooRealVar("p5","p5",-0.000001,-10.,10.)
    c6 = RooRealVar("p6","p6",-0.000001,-0.01,0.01)

    alpha = RooRealVar("#alpha","#alpha",-0.1,-10.0,-0.00001)
    polyset = RooArgList(c0,c1,c2,c3,c4,c5)

    gFraMMKK = RooRealVar("f_{gauss}","f_{gauss}",0.3,0.0,1.0)
    nSigMMKK = RooRealVar("n_{sig}","n_{sig}",10000,0.,10E6)
    nBkgMMKK = RooRealVar("n_{bkg}","n_{bkg}",10000,0.,10E6)

    ##pdf_bkg = RooChebychev("pdf_bkg","pdf_bkg",mmtt_mass,polyset)
    pdf_bkg  = RooExponential("pdf_bkg","pdf_bkg",mmtt_mass,alpha)
    pdf_sig1 = RooGaussian("pdf_sig1","pdf_sig1",mmtt_mass,mean,sigma1)
    pdf_sig2 = RooGaussian("pdf_sig2","pdf_sig2",mmtt_mass,mean,sigma2)

    pdf_sig  = RooAddPdf("pdf_sig","pdf_sig",pdf_sig1,pdf_sig2,gFraMMKK)
    pdf_tot = RooAddPdf("pdf_tot","pdf_tot",RooArgList(pdf_sig,pdf_bkg),RooArgList(nSigMMKK,nBkgMMKK))

    nfit = 0

    bZeroFit = pdf_tot.fitTo(bZeroFitData,Range(fitbZeromin,fitbZeromax),RooFit.PrintLevel(-1), RooFit.NumCPU(numcpus),RooFit.Save())
    nfit += 1

    gPad.SetRightMargin(0.1)

    plotpad.cd()

    bZeroFrame = mmtt_mass.frame(Range(fitbZeromin,fitbZeromax),Title("B_0^s mass"))

    bZeroFitData.plotOn(bZeroFrame)
    pdf_tot.plotOn(bZeroFrame,RooFit.Normalization(1.0/float(nfit)),Name("Pdf"))
    bZeroFitData.plotOn(bZeroFrame,Name("Data"))

    pdf_tot.paramOn(bZeroFrame,RooFit.Layout(0.78,0.99,0.99))

    bZeroFrame.Draw()

    pullpad.cd()
    hpull = bZeroFrame.pullHist("Data","Pdf")
    pullframe = mmtt_mass.frame(Title("Pull Distribution"),Range(fitbZeromin,fitbZeromax))
    pullframe.SetAxisRange (-5.0,5.0,"Y")
    #pullframe.GetXaxis().SetTitleSize(0.04)
    #pullframe.GetYaxis().SetTitleSize(0.03)
    ROOT.gStyle.SetTitleFontSize(0.07)
    pullframe.addPlotable(hpull,"P")
    pullframe.Draw()

    lineup = TLine(fitbZeromin,4.0,fitbZeromax,4.0)
    lineup.SetLineColor(kRed)
    lineup.SetLineStyle(kDashed)
    lineup.Draw()

    linedw = TLine(fitbZeromin,-4.0,fitbZeromax,-4.0)
    linedw.SetLineColor(kRed)
    linedw.SetLineStyle(kDashed)
    linedw.Draw()

    bcanvas.SaveAs("b0_fit" + region + ".png")
    bcanvas.SaveAs("b0_fit" + region + ".root")

    bZeroFrameTT = tt_mass.frame(Range(phimin,phimax),Title("#phi from B_0^s mass"))
    bZeroFitData.plotOn(bZeroFrameTT)
    bZeroFrameTT.Draw()
    bcanvas.SaveAs("b0_TT" + region + ".png")
    bcanvas.SaveAs("b0_TT" + region + ".root")

    bZeroFrameMM = mm_mass.frame(Range(jpsimin,jpsimax),Title("J/#Psi from B_0^s mass"))
    bZeroFitData.plotOn(bZeroFrameMM)
    bZeroFrameMM.Draw()
    bcanvas.SaveAs("b0_MM" + region + ".png")
    bcanvas.SaveAs("b0_MM" + region + ".root")

if args.binwise is not None:

    #_,bins,_ = plt.hist([],range=[xmin,xmax],bins=args.binwise)
    step = (xmax - xmin)/(float(args.binwise))

    scalingData = theData.Clone("binwiseData")
    bwcanvas = TCanvas("bwcanvas","bwcanvas",1200,800)

    phimean = 1.019
    gammavalue = 0.0012

    fitphimin = 1.01
    fitphimax = 1.03

    binwSigma = RooRealVar("#sigma","#sigma",0.0013)
    binwGamma = RooRealVar("#Gamma","#Gamma",gammavalue,0.001,0.015)
    binwMean = RooRealVar("m_{kk}","m_{kk}",phimean,phimean-0.005,phimean+0.005);

    # B_1     = RooRealVar ( "B_1"    , "B_1 "    , 0.3  , -20   , 100   )
    # B_2     = RooRealVar ( "B_2"    , "B_2"    , 0.3  , -20   , 100   )
    # B_3     = RooRealVar ( "B_3"    , "B_3"    , 0.3  , -20   , 100   )
    # B_4     = RooRealVar ( "B_4"    , "B_4"    , 0.3  , -20   , 100   )
    # B_5     = RooRealVar ( "B_5"    , "B_5"    , 0.3  , -20   , 100   )
    # B_6     = RooRealVar ( "B_6"    , "B_6"    , 0.3  , -20   , 100   )

    a0 = RooRealVar("p0","p0",0.001,-10.,10.)
    a1 = RooRealVar("p1","p1",0.001,-10.,10.)
    a2 = RooRealVar("p2","p2",-0.00001,-10.,10.)
    a3 = RooRealVar("p3","p3",-0.000001,-10.,10.)
    a4 = RooRealVar("p4","p4",-0.000001,-10.,10.)
    a5 = RooRealVar("a5","a5",-0.000001,-10.,10.)
    poliset = RooArgList(a0,a1,a2,a3)#,a4)

    for i in range(args.binwise-1):

        bwcanvas.Clear()

        pullpad = TPad("pullpad","pullpad",0.0,0.05,1.0,0.33)
        plotpad = TPad("histopad","histopad",0.0,0.35,1.0,1.0)
        plotpad.SetFillStyle(4004)
        pullpad.SetFillStyle(4004)
        plotpad.Draw()
        pullpad.Draw()

        lowedge = step * i + xmin #bins[i]
        upedge  = step * i + step + xmin #bins[i+1]

        scalingData = scalingData.reduce("xM>" + str(lowedge))

        thisData = scalingData.reduce("xM<" + str(upedge))

        # gaussFrac = RooRealVar("s","fraction of component 1 in kkSig",0.3,0.0,1.0)
        nSigBinW = RooRealVar("nSig_{bw}","nSig_{bw}",theData.numEntries()*0.3,0.0,thisData.numEntries()*1.5)
        nBkgBinW = RooRealVar("nBkg_{bw}","nBkg_{bw}",theData.numEntries()*0.7,0.0,thisData.numEntries()*1.5)

        binwSig = RooVoigtian("binwSig","binwSig",tt_mass,binwMean,binwGamma,binwSigma)
        #binwSig = RooGaussian("binwSig","binwSig",tt_mass,binwMean,binwGamma)#,binwSigma)
        #binwBkg = RooBernstein("binwBkg" , "binwBkg", tt_mass, RooArgList(B_1, B_2,B_3,B_4))#,B_5) )#,B_6))
        binwBkg = RooChebychev("binwBkg","Background",tt_mass,poliset)
        binwTot = RooAddPdf("binwTot","binwTot",RooArgList(binwSig,binwBkg),RooArgList(nSigBinW,nBkgBinW))

        # binwGamma.setConstant(ROOT.kTRUE)
        # binwMean.setConstant(ROOT.kTRUE)

        print("Fitting range " + str(lowedge) + " - " + str(upedge) + " : " + str(thisData.numEntries()))

        binwfit = binwTot.fitTo(thisData,Range(fitphimin,fitphimax),RooFit.PrintLevel(-1), RooFit.NumCPU(numcpus),RooFit.Save())

        plotpad.cd()

        binwFrame = tt_mass.frame(Range(fitphimin,fitphimax),Title("#phi binw mass [" + str(lowedge) + "-" + str(upedge) + "]"))

        thisData.plotOn(binwFrame,Name("Data"))
        binwTot.plotOn(binwFrame,RooFit.Normalization(1.0),Name("Pdf"))
        thisData.plotOn(binwFrame)
        binwTot.paramOn(binwFrame,RooFit.Layout(0.75,0.99,0.99))

        binwFrame.Draw()

        pullpad.cd()
        hpull = binwFrame.pullHist("Data","Pdf")
        pullframe = mmtt_mass.frame(Title("Pull Distribution"),Range(fitphimin,fitphimax))
        pullframe.SetAxisRange (-5.0,5.0,"Y")
        #pullframe.GetXaxis().SetTitleSize(0.04)
        #pullframe.GetYaxis().SetTitleSize(0.03)
        ROOT.gStyle.SetTitleFontSize(0.07)
        pullframe.addPlotable(hpull,"P")
        pullframe.Draw()

        lineup = TLine(fitphimin,4.0,fitphimax,4.0)
        lineup.SetLineColor(kRed)
        lineup.SetLineStyle(kDashed)
        lineup.Draw()

        linedw = TLine(fitphimin,-4.0,fitphimax,-4.0)
        linedw.SetLineColor(kRed)
        linedw.SetLineStyle(kDashed)
        linedw.Draw()

        bwcanvas.SaveAs("binwise_phi_xM_"  + str(lowedge) + "_" + str(upedge) + ".png")
        bwcanvas.SaveAs("binwise_phi_xM_"  + str(lowedge) + "_" + str(upedge) + ".root")

if args.nofit and args.nofitmm:

    jcanvas = TCanvas("jcanvas","jcanvas",1400,800)

    jPsiFitData = theData.Clone("fitB0Data")

    fitjPsimin = 3.0
    fitjPsimax = 3.2

    if binnedfit:
        mm_mass.setBins(100)
        jPsiFitData = theData.binnedClone("binnedTrakData")

    mean = RooRealVar("m","m",3.09,3.06,3.1);
    sigma1 = RooRealVar("#sigma_{1}","#sigma_{1}",0.01,0.001,0.1);
    sigma2 = RooRealVar("#sigma_{2}","#sigma_{2}",0.0011,0.001,0.1);

    # c0 = RooRealVar("p0","p0",0.001,-5.,5.)
    # c1 = RooRealVar("p1","p1",0.001,-4.,4.)
    # c2 = RooRealVar("p2","p2",-0.00001,-4.,4.)
    # c3 = RooRealVar("p3","p3",-0.000001,-4.,4.)
    # c4 = RooRealVar("p4","p4",-0.000001,-2.,2.)
    # c5 = RooRealVar("p5","p5",-0.000001)
    # c6 = RooRealVar("p6","p6",-0.000001,-0.01,0.01)
    #
    # polyset = RooArgList(c0,c1,c2,c3,c4)

    gFraMMKK = RooRealVar("f_{gauss}","f_{gauss}",0.3,0.0,1.0)
    # nSigMMKK = RooRealVar("n_{sig}","n_{sig}",10000,0.,10E6)
    # nBkgMMKK = RooRealVar("n_{bkg}","n_{bkg}",10000,0.,10E6)
    #
    # pdf_bkg = RooChebychev("pdf_bkg","pdf_bkg",mm_mass,polyset)
    # pdf_sig1 = RooGaussian("pdf_sig1","pdf_sig1",mm_mass,mean,sigma1)
    # pdf_sig2 = RooGaussian("pdf_sig2","pdf_sig2",mm_mass,mean,sigma2)
    #
    # pdf_sig  = RooAddPdf("pdf_sig","pdf_sig",pdf_sig1,pdf_sig2,gFraMMKK)
    # pdf_tot = RooAddPdf("pdf_tot","pdf_tot",RooArgList(pdf_sig1,pdf_bkg),RooArgList(nSigMMKK,nBkgMMKK))

    pdf_sig1 = RooGaussian("pdf_sig1","pdf_sig1",mm_mass,mean,sigma1)
    pdf_sig2 = RooGaussian("pdf_sig2","pdf_sig2",mm_mass,mean,sigma2)
    pdf_tot  = RooAddPdf("pdf_sig","pdf_sig",pdf_sig1,pdf_sig2,gFraMMKK)

    nfit = 0

    jPsiFit = pdf_tot.fitTo(jPsiFitData,Range(fitjPsimin+0.005,fitjPsimax-0.005),RooFit.PrintLevel(-1), RooFit.NumCPU(numcpus),RooFit.Save())
    nfit += 1

    jPsiFrame = mm_mass.frame(Range(fitjPsimin+0.005,fitjPsimax-0.005),Title("J/#Psi mass"))

    jPsiFitData.plotOn(jPsiFrame)

    pdf_tot.plotOn(jPsiFrame,RooFit.Normalization(1.0/float(nfit)))
    jPsiFitData.plotOn(jPsiFrame)
    pdf_tot.paramOn(jPsiFrame,RooFit.Layout(0.7,0.99,0.99))
    jPsiFrame.Draw()
    jcanvas.SaveAs("mm_fit" + region + ".png")
    jcanvas.SaveAs("mm_fit" + region + ".root")
