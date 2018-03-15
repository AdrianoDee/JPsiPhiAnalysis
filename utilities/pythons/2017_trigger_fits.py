
import ROOT
from ROOT import TLine
from ROOT import TFile,TH1,TH1F,TCanvas,TNtuple,TTreeReader,TTreeReaderValue
from ROOT import RooFit
from ROOT.RooFit import Layout
from ROOT import RooStats, gPad, RooAbsData, RooAbsReal
from ROOT.RooAbsReal import Relative

RooAbsData.setDefaultStorageType ( RooAbsData.Tree )
from array import array
import sys

import argparse
import math

parser = argparse.ArgumentParser()
parser.add_argument('--debug', action='store_true')
parser.add_argument('--binned', action='store_true')
parser.add_argument('--nonprompt', action='store_true')
parser.add_argument('--prompt', action='store_true')
parser.add_argument('--nofit', action='store_true')
parser.add_argument('--nofitkk', action='store_false')
parser.add_argument('--nofitb0', action='store_false')
parser.add_argument('--noplot', action='store_false')
parser.add_argument('--notrig', action='store_false')
parser.add_argument('--noreco', action='store_false')
parser.add_argument('--dosidebands', action='store_false')
parser.add_argument('--ptcuts', type=float, default=None)

#                    help='number of epochs')
#parser.add_argument('--batch_size', type=int, default=64)
#
args = parser.parse_args()

debugging = args.debug
binnedfit = args.binned




region = "_overall_"
cuts = "_"

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

name_2Trk = "2Trak2Trig_tree.root" #mmkk 2017 bcdef Jan 18 run
name_2Trk2Mu = "2Trak2MuTrig_tree.root"
name_2Mu = "2Mu2Trig_tree.root"

file_2Trk = TFile(name_2Trk,"READ")
file_2Trk.ls()
trakTree = (file_2Trk.Get("outuple"))

file_2Trk2Mu = TFile(name_2Trk2Mu,"READ")
file_2Trk2Mu.ls()
ttMuMuTree = (file_2Trk2Mu.Get("outuple"))

file_2Mu = TFile(name_2Mu,"READ")
file_2Mu.ls()
mumuTree = (file_2Mu.Get("outuple"))

phimin = 1.0
phimax = 1.05

jpsimin = 2.85
jpsimax = 3.35

xmin = 4.0
xmax = 6.0

bZeromin = 5.15
bZeromax = 5.55

ptmin = 0.0
ptmax = 100.0

tt_mass = RooRealVar("ttM","ttM",phimin,phimax)
tt_mass_trig = RooRealVar("trigtrigM","M(KK)_{L} ;M(KK)_{L3} [GeV]",phimin,phimax)

mm_mass = RooRealVar("mmM","mmM",jpsimin,jpsimax)
mm_mass_trig = RooRealVar("trigtrigM","M(#mu#mu)_{L3} ;M(#mu#mu)_{L3} [GeV]",jpsimin,jpsimax)

mmtt_mass = RooRealVar("xM","xM",xmin,xmax)
mmtt_mass_trig = RooRealVar("xTrigM","M(KK#mu#mu)_{L3}" + region + " ;M(KK#mu#mu)_{L3}[GeV]",xmin,xmax)

pt_pos = RooRealVar("trigp_pT","trigp_pT",ptmin,ptmax)
pt_neg = RooRealVar("trign_pT","trign_pT",ptmin,ptmax)

lxysig = RooRealVar("lxysig","l(xy) sign.;l(xy) sign.",0.0,50.0)

c = TCanvas("canvas","canvas",1200,800)

trakData = RooDataSet("trakData","trakData",trakTree,RooArgSet(tt_mass,tt_mass_trig,pt_pos,pt_neg))
mumuData = RooDataSet("mumuData","mumuData",mumuTree,RooArgSet(mm_mass,mm_mass_trig))
mmttData = RooDataSet("mmttData","mmttData",ttMuMuTree,RooArgSet(mmtt_mass,mmtt_mass_trig,lxysig))

if args.nonprompt:
    mmttData = mmttData.reduce("lxysig > 3.0")
    region = "_nprompt_"
if args.prompt:
    mmttData = mmttData.reduce("lxysig < 2.0")
    region = "_promt_"

if args.ptcuts is not None:
    trakData = trakData.reduce("trigp_pT > " + str(args.ptcuts))
    trakData = trakData.reduce("trign_pT > " + str(args.ptcuts))
    cuts += "pt_" + str(args.ptcuts) + "_"
#### #### Plotting variables
#### TrakTrak Data

if args.noplot:

    print("TrakTrak data plotting . . .")


    print("All : " + str(trakData.numEntries()))
    ttFrame = tt_mass.frame()
    trakData.plotOn(ttFrame)
    ttFrame.Draw()
    c.SaveAs("tt_mass" + cuts + ".png")

    tt_trig_Frame = tt_mass_trig.frame()
    trakData.plotOn(tt_trig_Frame)
    tt_trig_Frame.Draw()
    c.SaveAs("tt_mass" + cuts + "trigger.png")

    #### MuMu Data
    print("MuMu data plotting . . .")

    print("All : " + str(mumuData.numEntries()))
    mumuFrame = mm_mass.frame()
    mumuData.plotOn(mumuFrame)
    mumuFrame.Draw()
    c.SaveAs("mm_mass.png")

    mumu_trig_Frame = mm_mass_trig.frame()
    mumuData.plotOn(mumu_trig_Frame)
    mumu_trig_Frame.Draw()
    c.SaveAs("mm_mass_trigger.png")
    # plotmax = trakTree.GetMaximum("ttM")

    #### X Data
    print("TrakTrakMuMu data plotting . . .")

    print("All : " + str(mmttData.numEntries()))
    mmttFrame = mmtt_mass.frame()
    mmttData.plotOn(mmttFrame)
    mmttFrame.Draw()
    c.SaveAs("mmtt_mass" + region + ".png")

    mmtt_trig_Frame = mmtt_mass_trig.frame()
    mmttData.plotOn(mmtt_trig_Frame)
    mmtt_trig_Frame.Draw()
    c.SaveAs("mmtt_mass_trigger" + region + ".png")

    #### X lxy Data
    print("TrakTrakMuMu data plotting . . .")

    print("All : " + str(mmttData.numEntries()))
    lxyFrame = lxysig.frame()
    mmttData.plotOn(lxyFrame)
    lxyFrame.Draw()
    c.SaveAs("lxy" + region + ".png")

    mmttcuts = [4.0,4.2,4.3,4.4,4.6,4.8,5.0]
    trigcands = []
    patcands = []

    hist_mmtt      = (mmttData.createHistogram(mmtt_mass,mmtt_mass_trig)).ProjectionX("hist_mmtt")
    hist_mmtt_trig = (mmttData.createHistogram(mmtt_mass,mmtt_mass_trig)).ProjectionY("hist_mmtt_trig")

    mmtt_all      = hist_mmtt.Integral(0, hist_mmtt.GetNbinsX())
    mmtt_trig_all = hist_mmtt_trig.Integral(0, hist_mmtt_trig.GetNbinsX())
    print("####TRIGGER")
    for i in range(len(mmttcuts)) :
        integral = hist_mmtt_trig.Integral(hist_mmtt_trig.FindBin(mmttcuts[i]), hist_mmtt_trig.GetNbinsX())
        print("#No. cands. > " + str(mmttcuts[i]) + " : " + str(integral) + " fraction : " +  str(integral/mmtt_trig_all))
        trigcands.append(integral)
        
    print(trigcands)
    print("####PAT")
    for i in range(len(mmttcuts)) :
        integral = hist_mmtt.Integral(hist_mmtt.FindBin(mmttcuts[i]), hist_mmtt.GetNbinsX())
        print("#No. cands. > " + str(mmttcuts[i]) + " : " + str(integral) + " fraction : " +  str(integral/mmtt_all))
        patcands.append(integral)

    print(patcands)

    hist_mmtt_trig.SetTitle("M(#mu#muKK) - L3 " + region + " ;M(#mu#muKK)_{L3} [GeV]; no. entries")
    hist_mmtt_trig.SetFillColor(ROOT.kRed)
    hist_mmtt_trig.SetLineColor(ROOT.kRed)
    hist_mmtt_trig.SetFillStyle(3002)

    hist_mmtt.SetTitle("M(#mu#muKK)" + region + " ;M(#mu#muKK) [GeV]; no. entries")
    hist_mmtt.SetLineColor(ROOT.kBlue)

    hist_mmtt.Draw()
    hist_mmtt_trig.Draw("SAME")

    gPad.BuildLegend(0.2,0.2,0.4,0.35)
    c.SaveAs("mmtt_mass_histos" + region + ".png")
    c.SaveAs("mmtt_mass_histos" + region + ".root")

    hist_lxy      = (mmttData.createHistogram(lxysig,mmtt_mass)).ProjectionX("hist_lxy")
    hist_lxy.SetTitle("L_{XY}/#sigma_{L_{XY}}; L_{XY}/#sigma_{L_{XY}}; no. entries")
    hist_lxy.SetMarkerStyle(22)
    hist_lxy.SetMarkerSize(1.5)

    hist_lxy.GetYaxis().SetRangeUser(0.0,hist_lxy.GetMaximum()*1.3)
    linePrompt    = TLine(1.5,0.0,1.5,hist_lxy.GetMaximum())
    linePrompt.SetLineColor(ROOT.kBlue)
    linePrompt.SetLineWidth(3)
    lineNonPrompt = TLine(3.0,0.0,3.0,hist_lxy.GetMaximum())
    lineNonPrompt.SetLineColor(ROOT.kRed)
    lineNonPrompt.SetLineWidth(3)

    promptfraction = hist_lxy.Integral(hist_lxy.FindBin(0.0), hist_lxy.FindBin(1.5))
    nonpromptfraction = hist_lxy.Integral(hist_lxy.FindBin(3.0), hist_lxy.GetNbinsX())
    allfraction = hist_lxy.Integral(hist_lxy.FindBin(0.0), hist_lxy.GetNbinsX())
    hist_lxy.Draw("P")
    linePrompt.Draw()
    lineNonPrompt.Draw()

    print("###Fractions:")
    print("Prompt:      " + str(100.0*promptfraction/(allfraction)) + " %")
    print("NonPrompt:   " + str(100.0*nonpromptfraction/(allfraction)) + " %")
    c.SaveAs("mmtt_xly_histos.png")
    c.SaveAs("mmtt_xly_histos.root")


    hist_tt      = (trakData.createHistogram(tt_mass,tt_mass_trig,100,100,"","TT")).ProjectionX("hist_tt")
    hist_tt_trig = (trakData.createHistogram(tt_mass,tt_mass_trig,100,100,"","TTTrig")).ProjectionY("hist_tt_trig")

    hist_tt_trig.SetTitle("M(KK) - L3  ;M(KK)_{L3} [GeV]; no. entries")
    hist_tt_trig.SetFillColor(ROOT.kRed)
    hist_tt_trig.SetLineColor(ROOT.kRed)
    hist_tt_trig.SetFillStyle(3002)

    hist_tt.SetTitle("M(KK) ;M(KK) [GeV]; no. entries" )
    hist_tt.SetLineColor(ROOT.kBlue)

    hist_tt_trig.Draw()
    hist_tt.Draw("SAME")

    gPad.BuildLegend(0.4,0.5,0.6,0.65)
    c.SaveAs("tt_mass_histos.png")
    c.SaveAs("tt_mass_histos.root")

    a = 0.95
    b = 1.30
    C = 1.15
    d = 1.05

    trigRange1 = hist_tt_trig.Integral(hist_tt_trig.FindBin(a), hist_tt_trig.FindBin(b))
    trigRange2 = hist_tt_trig.Integral(hist_tt_trig.FindBin(a), hist_tt_trig.FindBin(C))
    trigRange3 = hist_tt_trig.Integral(hist_tt_trig.FindBin(a), hist_tt_trig.FindBin(d))

    patRange1 = hist_tt.Integral(hist_tt.FindBin(a), hist_tt.FindBin(b))
    patRange2 = hist_tt.Integral(hist_tt.FindBin(a), hist_tt.FindBin(C))
    patRange3 = hist_tt.Integral(hist_tt.FindBin(a), hist_tt.FindBin(d))

    print("####PAT")
    print("#No. cands. [" + str(a) + " - " + str(b) + "] : " + str(trigRange1) + " fraction : " +  str(trigRange1/trigRange1))
    print("#No. cands. [" + str(a) + " - " + str(C) + "] : " + str(trigRange2) + " fraction : " +  str(trigRange2/trigRange1))
    print("#No. cands. [" + str(a) + " - " + str(d) + "] : " + str(trigRange3) + " fraction : " +  str(trigRange3/trigRange1))

    print("####TRIGGER")
    print("#No. cands. [" + str(a) + " - " + str(b) + "] : " + str(patRange1) + " fraction : " +  str(patRange1/patRange1))
    print("#No. cands. [" + str(a) + " - " + str(C) + "] : " + str(patRange2) + " fraction : " +  str(patRange2/patRange1))
    print("#No. cands. [" + str(a) + " - " + str(d) + "] : " + str(patRange3) + " fraction : " +  str(patRange3/patRange1))


    hist_mm      = (mumuData.createHistogram(mm_mass,mm_mass_trig)).ProjectionX("hist_mm")
    hist_mm_trig = (mumuData.createHistogram(mm_mass,mm_mass_trig)).ProjectionY("hist_mm_trig")

    hist_mm_trig.SetTitle("M(#mu#mu)- L3  ;M(#mu#mu)_{L3} [GeV]; no. entries")
    hist_mm_trig.SetFillColor(ROOT.kRed)
    hist_mm_trig.SetLineColor(ROOT.kRed)
    hist_mm_trig.SetFillStyle(3002)

    hist_mm.SetTitle("M(#mu#mu) ;M(#mu#mu)[GeV]; no. entries")
    hist_mm.SetLineColor(ROOT.kBlue)


    hist_mm.Draw()
    hist_mm_trig.Draw("SAME")

    gPad.BuildLegend(0.65,0.5,0.85,0.65)
    c.SaveAs("mm_mass_histos.png")
    c.SaveAs("mm_mass_histos.root")

if args.nofit:
    sys.exit()

#### #### Done plotting

#### #### Fitting
### Trak trak
#

traKFitData = trakData.Clone("fitTrakData")

if binnedfit:
    tt_mass.setBins(100)
    traKFitData = trakData.binnedClone("binnedTrakData")

phimean = 1.019
gammavalue = 0.01

fitphimin = 1.0
fitphimax = 1.04


if args.nofitkk and args.noreco:


    ## RECO

    kkSigma = RooRealVar("#sigma","#sigma",0.0013)
    kkGamma = RooRealVar("#Gamma","#Gamma",gammavalue,0.001,0.015)
    kkMean = RooRealVar("mean","mean",phimean,phimean-0.007,phimean+0.007);

    # B_1     = RooRealVar ( "B_1"    , "B_1 "    , 0.3  , -20   , 100   )
    # B_2     = RooRealVar ( "B_2"    , "B_2"    , 0.3  , -20   , 100   )
    # B_3     = RooRealVar ( "B_3"    , "B_3"    , 0.3  , -20   , 100   )
    # B_4     = RooRealVar ( "B_4"    , "B_4"    , 0.3  , -20   , 100   )
    # B_5     = RooRealVar ( "B_5"    , "B_5"    , 0.3  , -20   , 100   )
    # B_6     = RooRealVar ( "B_6"    , "B_6"    , 0.3  , -20   , 100   )

    a0 = RooRealVar("p0","p0",0.001,-2.,2.)
    a1 = RooRealVar("p1","p1",0.001,-2.,2.)
    a2 = RooRealVar("p2","p2",-0.00001,-2.,2.)
    a3 = RooRealVar("p3","p3",-0.000001,-2,2.)
    a4 = RooRealVar("p4","p4",-0.000001,-2.,2.)
    # a5 = RooRealVar("a5","a5",-0.000001,-2.,2.)
    poliset = RooArgList(a0,a1,a2,a3,a4)

    # gaussFrac = RooRealVar("s","fraction of component 1 in kkSig",0.3,0.0,1.0)
    nSigKK = RooRealVar("nSig","nSig",1E5,0.,10E6)
    nBkgKK = RooRealVar("nBkg","nBkg",5E5,0.,10E6)

    kkSig = RooVoigtian("kkSig","kkSig",tt_mass,kkMean,kkGamma,kkSigma)
    #kkBkg = RooBernstein("kkBkg" , "kkBkg", tt_mass, RooArgList(B_1, B_2,B_3,B_4))#,B_5) )#,B_6))
    kkBkg = RooChebychev("kkBkg","Background",tt_mass,poliset)
    kkTot = RooAddPdf("kkTot","kkTot",RooArgList(kkSig,kkBkg),RooArgList(nSigKK,nBkgKK))

    kkGamma.setConstant(ROOT.kTRUE)
    kkMean.setConstant(ROOT.kTRUE)

    nfit = 0

    # kkfit = kkTot.fitTo(traKFitData,Range(fitphimin+0.005,fitphimax-0.005),RooFit.PrintLevel(-1), RooFit.NumCPU(7),RooFit.Save())
    # nfit +=1

    if not debugging:

        kkMean.setConstant(ROOT.kFALSE)
        # kkfit = kkTot.fitTo(traKFitData,Range(fitphimin+0.005,fitphimax-0.005),RooFit.PrintLevel(-1), RooFit.NumCPU(7),RooFit.Save())
        # nfit +=1

        kkGamma.setConstant(ROOT.kFALSE)
        # kkfit = kkTot.fitTo(traKFitData,Range(fitphimin+0.005,fitphimax-0.005),RooFit.PrintLevel(-1), RooFit.NumCPU(7),RooFit.Save())
        # nfit +=1

        kkfit = kkTot.fitTo(traKFitData,Range(fitphimin+0.005,fitphimax-0.005), RooFit.NumCPU(7),RooFit.Save())
        nfit +=1

    sigmaside_kk = math.sqrt(kkGamma.getValV()**2 + kkSigma.getValV()**2)

    kkFrame = tt_mass.frame(Range(fitphimin+0.005,fitphimax-0.005))

    leftlowside_kk = -6.*sigmaside_kk + kkMean.getValV()
    leftupside_kk = -4.*sigmaside_kk + kkMean.getValV()
    rightlowside_kk = +4.*sigmaside_kk + kkMean.getValV()
    rightupside_kk = +6.*sigmaside_kk + kkMean.getValV()

    signallow = -3.*sigmaside_kk + kkMean.getValV()
    signalup = +3.*sigmaside_kk + kkMean.getValV()

    trakData.plotOn(kkFrame)

    kkTot.plotOn(kkFrame,RooFit.Normalization(1.0/float(nfit)))
    traKFitData.plotOn(kkFrame)
    kkTot.paramOn(kkFrame,RooFit.Layout(0.57,0.99,0.65))

    kkFrame.Draw()

    c.SaveAs("kk_Phi_fit.png")
    c.SaveAs("kk_Phi_fit.root")



if args.nofitkk and args.notrig:

    ##L3

    kkSigma = RooRealVar("#sigma","#sigma",0.0013)
    kkGamma = RooRealVar("#Gamma","#Gamma",gammavalue,0.001,0.01)
    kkMean = RooRealVar("mean","mean",phimean,phimean-0.007,phimean+0.007);

    # B_1     = RooRealVar ( "B_1"    , "B_1 "    , 0.3  , -20   , 100   )
    # B_2     = RooRealVar ( "B_2"    , "B_2"    , 0.3  , -20   , 100   )
    # B_3     = RooRealVar ( "B_3"    , "B_3"    , 0.3  , -20   , 100   )
    # B_4     = RooRealVar ( "B_4"    , "B_4"    , 0.3  , -20   , 100   )
    # B_5     = RooRealVar ( "B_5"    , "B_5"    , 0.3  , -20   , 100   )
    # B_6     = RooRealVar ( "B_6"    , "B_6"    , 0.3  , -20   , 100   )

    a0 = RooRealVar("p0","p0",0.001,-2.,2.)
    a1 = RooRealVar("p1","p1",0.001,-2.,2.)
    a2 = RooRealVar("p2","p2",-0.00001,-2.,2.)
    a3 = RooRealVar("p3","p3",-0.000001,-2,2.)
    a4 = RooRealVar("p4","p4",-0.000001,-2.,2.)
    # a5_t = RooRealVar("a5_t","a5_t",a5.getValV(),-2.,2.)
    poliset = RooArgList(a0,a1,a2,a3,a4)

    # gaussFrac = RooRealVar("s","fraction of component 1 in kkSig",0.3,0.0,1.0)
    nSig = RooRealVar("nSig","nSig",2E5,0.,10E6)
    nBkg = RooRealVar("nBkg","nBkg",2E5,0.,10E6)

    kkSig = RooVoigtian("kkSig","kkSig",tt_mass_trig,kkMean,kkGamma,kkSigma)
    # kkBkg = RooBernstein("kkBkg" , "kkBkg", tt_mass_trig   , RooArgList(B_1, B_2,B_3))#,B_4,B_5) )#,B_6))
    kkBkg = RooChebychev("kkBkg","kkBkg",tt_mass_trig,poliset)
    kkTot = RooAddPdf("kkTot","kkTot",RooArgList(kkSig,kkBkg),RooArgList(nSig,nBkg))

    kkGamma.setConstant(ROOT.kTRUE)
    kkMean.setConstant(ROOT.kTRUE)

    nfit = 0

    # kkFit = kkTot.fitTo(traKFitData,Range(fitphimin+0.005,fitphimax-0.005),RooFit.PrintLevel(-1), RooFit.NumCPU(7),RooFit.Save())
    # nfit +=1

    if not debugging:

        kkMean.setConstant(ROOT.kFALSE)
        # kkFit = kkTot.fitTo(traKFitData,Range(fitphimin+0.005,fitphimax-0.005),RooFit.PrintLevel(-1), RooFit.NumCPU(7),RooFit.Save())
        # nfit +=1

        kkGamma.setConstant(ROOT.kFALSE)
        # kkFit = kkTot.fitTo(traKFitData,Range(fitphimin+0.005,fitphimax-0.005),RooFit.PrintLevel(-1), RooFit.NumCPU(7),RooFit.Save())
        # nfit +=1

        kkFit = kkTot.fitTo(traKFitData,Range(fitphimin+0.005,fitphimax-0.005), RooFit.NumCPU(7),RooFit.Save())
        nfit +=1

    sigmaside_kk = math.sqrt(kkGamma.getValV()**2 + kkSigma.getValV()**2)

    kkFrame = tt_mass_trig.frame(Range(fitphimin+0.005,fitphimax-0.005))

    leftlowside_kk = -6.*sigmaside_kk + kkMean.getValV()
    leftupside_kk = -4.*sigmaside_kk + kkMean.getValV()
    rightlowside_kk = +4.*sigmaside_kk + kkMean.getValV()
    rightupside_kk = +6.*sigmaside_kk + kkMean.getValV()

    signallow = -3.*sigmaside_kk + kkMean.getValV()
    signalup = +3.*sigmaside_kk + kkMean.getValV()

    trakData.plotOn(kkFrame)

    kkTot.plotOn(kkFrame,RooFit.Normalization(1.0/float(nfit)))
    traKFitData.plotOn(kkFrame)
    kkTot.paramOn(kkFrame,RooFit.Layout(0.57,0.99,0.65))

    kkFrame.Draw()

    c.SaveAs("kkTrig_Phi_Fit.png")
    c.SaveAs("kkTrig_Phi_Fit.root")

########################################################################

### B0s
##Reco

if args.nofitb0 and args.noreco:

    bcanvas = TCanvas("bcanvas","bcanvas",1400,800)

    bZeroFitData = mmttData.Clone("fitB0Data")

    fitbZeromin = 5.15
    fitbZeromax = 5.55

    if binnedfit:
        mmtt_mass.setBins(100)
        bZeroFitData = mmttData.binnedClone("binnedTrakData")

    bZeroFitData = (bZeroFitData.reduce("xM<5.55")).reduce("xM>5.15")

    mean = RooRealVar("m_{L3}","m_{L3}",5.38,5.31,5.41);
    sigma1 = RooRealVar("#sigma_{1}","#sigma_{1}",0.002,0.0005,0.05);
    sigma2 = RooRealVar("#sigma_{2}","#sigma_{2}",0.004,0.004,0.01);

    c0 = RooRealVar("p0","p0",0.001,-5.,5.)
    c1 = RooRealVar("p1","p1",0.001,-4.,4.)
    c2 = RooRealVar("p2","p2",-0.00001,-4.,4.)
    c3 = RooRealVar("p3","p3",-0.000001,-4.,4.)
    c4 = RooRealVar("p4","p4",-0.000001,-2.,2.)
    c5 = RooRealVar("p5","p5",-0.000001)
    c6 = RooRealVar("p6","p6",-0.000001,-0.01,0.01)

    polyset = RooArgList(c0,c1,c2,c3,c4)

    gFraMMKK = RooRealVar("f_{gauss}","f_{gauss}",0.3,0.0,1.0)
    nSigMMKK = RooRealVar("n_{sig}","n_{sig}",10000,0.,10E6)
    nBkgMMKK = RooRealVar("n_{bkg}","n_{bkg}",10000,0.,10E6)

    pdf_bkg = RooChebychev("pdf_bkg","pdf_bkg",mmtt_mass,polyset)
    pdf_sig1 = RooGaussian("pdf_sig1","pdf_sig1",mmtt_mass,mean,sigma1)
    pdf_sig2 = RooGaussian("pdf_sig2","pdf_sig2",mmtt_mass,mean,sigma2)

    pdf_sig  = RooAddPdf("pdf_sig","pdf_sig",pdf_sig1,pdf_sig2,gFraMMKK)
    pdf_tot = RooAddPdf("pdf_tot","pdf_tot",RooArgList(pdf_sig1,pdf_bkg),RooArgList(nSigMMKK,nBkgMMKK))

    nfit = 0

    bZeroFit = pdf_tot.fitTo(bZeroFitData,Range(fitbZeromin+0.005,fitbZeromax-0.005),RooFit.PrintLevel(-1), RooFit.NumCPU(7),RooFit.Save())
    nfit += 1

    bZeroFrame = mmtt_mass.frame(Range(fitbZeromin+0.005,fitbZeromax-0.005))

    bZeroFitData.plotOn(bZeroFrame)

    pdf_tot.plotOn(bZeroFrame,RooFit.Normalization(1.0/float(nfit)))
    bZeroFitData.plotOn(bZeroFrame)
    pdf_tot.paramOn(bZeroFrame,RooFit.Layout(0.7,0.99,0.99))

    bZeroFrame.Draw()

    gPad.SetRightMargin(0.1)
    bcanvas.SaveAs("b0_fit" + region + ".png")
    bcanvas.SaveAs("b0_fit" + region + ".root")


if args.nofitb0 and args.notrig:

    bTrigcanvas = TCanvas("bTrigcanvas","bTrigcanvas",1400,800)

    bZeroTriggerFitData = mmttData.Clone("fitB0TriggerData")

    fitbZeroTriggermin = 5.15
    fitbZeroTriggermax = 5.55

    if binnedfit:
        mmtt_mass_trig.setBins(100)
        bZeroTriggerFitData = mmttData.binnedClone("binnedTrakData")

    bZeroTriggerFitData = (bZeroTriggerFitData.reduce("xTrigM<5.55")).reduce("xTrigM>5.15")

    mmkkMean = RooRealVar("m_{L3}","m_{L3}",5.38,5.31,5.41);
    mmkkSigma1 = RooRealVar("#sigma_{1}","#sigma_{1}",0.002,0.0005,0.05);
    mmkkSigma2 = RooRealVar("#sigma_{2}","#sigma_{2}",0.001,0.001,0.02);

    c0 = RooRealVar("p0","p0",0.001,-100.,100.)
    c1 = RooRealVar("p1","p1",0.001,-10.,10.)
    c2 = RooRealVar("p2","p2",-0.00001,-10.,10.)
    c3 = RooRealVar("p3","p3",-0.000001,-10.,10.)
    c4 = RooRealVar("p4","p4",-0.000001,-10.,10.)
    c5 = RooRealVar("p5","p5",-0.000001)
    c6 = RooRealVar("p6","p6",-0.000001,-0.01,0.01)

    polyset = RooArgList(c0,c1,c2,c3,c4)

    gFrac = RooRealVar("f_{gauss}","f_{gauss}",0.05,0.0,1.0)
    nSig = RooRealVar("n_{sig}","n_{sig}",10000,0.,10E6)
    bBkg = RooRealVar("n_{bkg}","n_{bkg}",10000,0.,10E6)

    pdf_bkg = RooChebychev("pdf_bkg","pdf_bkg",mmtt_mass_trig,polyset)
    pdf_sig1 = RooGaussian("pdf_sig1","pdf_sig1",mmtt_mass_trig,mmkkMean,mmkkSigma1)
    pdf_sig2 = RooGaussian("pdf_sig2","pdf_sig2",mmtt_mass_trig,mmkkMean,mmkkSigma2)

    pdf_sig  = RooAddPdf("pdf_sig","pdf_sig",pdf_sig1,pdf_sig2,gFrac)
    pdf_tot = RooAddPdf("pdf_tot","pdf_tot",RooArgList(pdf_sig1,pdf_bkg),RooArgList(nSig,bBkg))

    nfit = 0

    bZeroTriggerFit = pdf_tot.fitTo(bZeroTriggerFitData,Range(fitbZeroTriggermin+0.005,fitbZeroTriggermax-0.005),RooFit.PrintLevel(-1), RooFit.NumCPU(7),RooFit.Save())
    nfit += 1

    bZeroTriggerFrame = mmtt_mass_trig.frame(Range(fitbZeroTriggermin+0.005,fitbZeroTriggermax-0.005))

    bZeroTriggerFitData.plotOn(bZeroTriggerFrame)

    sFrac = (nSig.getValV())/(bBkg.getValV()+nSig.getValV())
    bFrac = 1.0 - sFrac
    fSig2 = (1.0-gFrac.getValV()) * sFrac
    fSig1 = (gFrac.getValV())# * sFrac

    pdf_tot.plotOn(bZeroTriggerFrame,RooFit.Normalization(1.0/float(nfit)))
    # pdf_sig2.plotOn(bZeroTriggerFrame,LineColor(kRed),LineStyle(kDashed),RooFit.NormRange("full"),RooFit.Range("full")) #Normalization((float((fSig2))/(float(nfit))))#,int(bBkg.getValV()+nSig.getValV())) )
    # pdf_sig1.plotOn(bZeroTriggerFrame,LineColor(kGreen),LineStyle(kDashed),Normalization((float(fSig2)/float(nfit))))
    # pdf_bkg.plotOn(bZeroTriggerFrame,LineColor(kMagenta),LineStyle(kDotted),Normalization(float(bFrac)/float(nfit)))

    #bZeroTriggerFitData.plotOn(bZeroTriggerFrame)
    pdf_tot.paramOn(bZeroTriggerFrame,RooFit.Layout(0.7,0.99,0.9))

    bZeroTriggerFrame.Draw()

    gPad.SetRightMargin(0.1)
    bTrigcanvas.SaveAs("b0Trig_fit" + region + ".png")
    bTrigcanvas.SaveAs("b0Trig_fit" + region + ".root")

if args.dosidebands and args.noreco:
    print "ciao"
