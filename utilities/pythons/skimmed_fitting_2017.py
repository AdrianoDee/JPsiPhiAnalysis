
import ROOT
from ROOT import TLine
from ROOT import TFile,TH1,TH1F,TCanvas,TNtuple,TTreeReader,TTreeReaderValue
from ROOT import RooFit
from ROOT.RooFit import Layout
from ROOT import RooStats, gPad, RooAbsData, RooAbsReal, RooBinning
from ROOT.RooAbsReal import Relative

RooAbsData.setDefaultStorageType ( RooAbsData.Tree )
from array import array
import sys

import argparse
import math

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
# parser.add_argument('--notrig', action='store_false')
# parser.add_argument('--noreco', action='store_false')
# parser.add_argument('--dosidebands', action='store_false')
parser.add_argument('--ptcuts', type=float, default=None)

#                    help='number of epochs')
#parser.add_argument('--batch_size', type=int, default=64)
#
args = parser.parse_args()

debugging = args.debug
binnedfit = args.binned
numcpus = args.numcpu

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
ptmax = 100.0

tt_mass = RooRealVar("ttM","ttM",phimin,phimax)
mm_mass = RooRealVar("mmM","mmM",jpsimin,jpsimax)
mmtt_mass = RooRealVar("xM","xM",xmin,xmax)

tt_pt = RooRealVar("ttPt","ttPt",0.0,ptmax)
mm_pt = RooRealVar("mmPt","mmPt",0.0,ptmax)
mmtt_pt = RooRealVar("xPt","xPt",0.0,ptmax)

lxysig = RooRealVar("xL","l(xy) sign.;l(xy) sign.",0.0,1000.0)

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
    ttFrame = tt_mass.frame()
    theData.plotOn(ttFrame)
    ttFrame.Draw()
    c.SaveAs("tt_mass" + cuts + ".png")

    #### MuMu Data
    print("MuMu data plotting . . .")

    print("All : " + str(theData.numEntries()))
    mumuFrame = mm_mass.frame()
    theData.plotOn(mumuFrame)
    mumuFrame.Draw()
    c.SaveAs("mm_mass.png")

    #### X Data
    print("TrakTrakMuMu data plotting . . .")

    print("All : " + str(theData.numEntries()))
    mmttFrame = mmtt_mass.frame()
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
    hist_tt_pt     = (theData.createHistogram(mm_pt,mmtt_mass)).ProjectionX("hist_mm_pt")
    hist_mm_pt     = (theData.createHistogram(tt_pt,mmtt_mass)).ProjectionX("hist_tt_pt")
    hist_mmtt_pt   = (theData.createHistogram(mmtt_pt,mmtt_mass)).ProjectionX("hist_mmtt_pt")

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
    hist_tt.SetFillColor(ROOT.kBlue)
    hist_tt.SetLineColor(ROOT.kBlue)
    hist_tt.SetFillStyle(3002)
    hist_tt.Draw()
    c.SaveAs("tt_mass_histo" + region + ".png")
    c.SaveAs("tt_mass_histo" + region + ".root")

    ##PT hists

    hist_tt_pt.SetTitle("P_t(KK)" + region + " ;M(KK) [GeV]; no. entries")
    hist_tt_pt.SetFillColor(ROOT.kBlue)
    hist_tt_pt.SetLineColor(ROOT.kBlue)
    hist_tt_pt.SetFillStyle(3003)
    hist_tt_pt.Draw()
    c.SaveAs("tt_pt_histo" + region + ".png")
    c.SaveAs("tt_pt_histo" + region + ".root")

    hist_mm_pt.SetTitle("P_t(#mu#mu)" + region + " ;Pt(#mu#mu) [GeV]; no. entries")
    hist_mm_pt.SetFillColor(ROOT.kGreen)
    hist_mm_pt.SetLineColor(ROOT.kGreen)
    hist_mm_pt.SetFillStyle(3003)
    hist_mm_pt.Draw()
    c.SaveAs("mm_pt_histo" + region + ".png")
    c.SaveAs("mm_pt_histo" + region + ".root")

    hist_tt_pt.SetTitle("P_t(#mu#muKK)" + region + " ;Pt(#mu#muKK) [GeV]; no. entries")
    hist_tt_pt.SetFillColor(ROOT.kRed)
    hist_tt_pt.SetLineColor(ROOT.kRed)
    hist_tt_pt.SetFillStyle(3003)
    hist_tt_pt.Draw()
    c.SaveAs("tt_pt_histo" + region + ".png")
    c.SaveAs("tt_pt_histo" + region + ".root")

if args.nofit and args.nofitkk:

    ## Phi fitting

    traKFitData = theData.Clone("fitTrakData")

    if binnedfit:
        tt_mass.setBins(30)
        traKFitData = theData.binnedClone("binnedTrakData")

    phimean = 1.019
    gammavalue = 0.01

    fitphimin = 1.006
    fitphimax = 1.034

    kkSigma = RooRealVar("#sigma","#sigma",0.0013)
    kkGamma = RooRealVar("#Gamma","#Gamma",gammavalue,0.001,0.015)
    kkMean = RooRealVar("mean","mean",phimean,phimean-0.007,phimean+0.007);

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
        kkfit = kkTot.fitTo(traKFitData,Range(fitphimin+0.005,fitphimax-0.005),RooFit.PrintLevel(-1), RooFit.NumCPU(numcpus),RooFit.Save())
	nfit +=1

    sigmaside_kk = math.sqrt(kkGamma.getValV()**2 + kkSigma.getValV()**2)

    kkFrame = tt_mass.frame(Range(fitphimin+0.005,fitphimax-0.005))

    leftlowside_kk = -6.*sigmaside_kk + kkMean.getValV()
    leftupside_kk = -4.*sigmaside_kk + kkMean.getValV()
    rightlowside_kk = +4.*sigmaside_kk + kkMean.getValV()
    rightupside_kk = +6.*sigmaside_kk + kkMean.getValV()

    signallow = -3.*sigmaside_kk + kkMean.getValV()
    signalup = +3.*sigmaside_kk + kkMean.getValV()



    kkbins = RooBinning(-15,15)
    kkbins.addUniform(30,fitphimin+0.005,fitphimax-0.005)
    traKFitData.plotOn(kkFrame,RooFit.Binning(kkbins))
    kkTot.plotOn(kkFrame,RooFit.Normalization(1.0/float(nfit)))
    traKFitData.plotOn(kkFrame,RooFit.Binning(kkbins))
    kkTot.paramOn(kkFrame,RooFit.Layout(0.57,0.99,0.65))

    kkFrame.Draw()

    c.SaveAs("kk_Phi_fit" + region + ".png")
    c.SaveAs("kk_Phi_fit" + region + ".root")

if args.nofit and args.nofitb0:

    bcanvas = TCanvas("bcanvas","bcanvas",1400,800)

    bZeroFitData = theData.Clone("fitB0Data")

    fitbZeromin = 5.15
    fitbZeromax = 5.55

    if binnedfit:
        mmtt_mass.setBins(50)
        bZeroFitData = theData.binnedClone("binnedTrakData")

    bZeroFitData = (bZeroFitData.reduce("xM<5.55")).reduce("xM>5.15")

    mean = RooRealVar("m_{L3}","m_{L3}",5.35,5.2,5.4);
    sigma1 = RooRealVar("#sigma_{1}","#sigma_{1}",0.002,0.0005,0.05);
    sigma2 = RooRealVar("#sigma_{2}","#sigma_{2}",0.004,0.004,0.01);

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
    pdf_tot = RooAddPdf("pdf_tot","pdf_tot",RooArgList(pdf_sig1,pdf_bkg),RooArgList(nSigMMKK,nBkgMMKK))

    nfit = 0

    bZeroFit = pdf_tot.fitTo(bZeroFitData,Range(fitbZeromin+0.005,fitbZeromax-0.005),RooFit.PrintLevel(-1), RooFit.NumCPU(numcpus),RooFit.Save())
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

    bZeroFrameTT = tt_mass.frame(Range(phimin,phimax))
    bZeroFitData.plotOn(bZeroFrameTT)
    bZeroFrameTT.Draw()
    bcanvas.SaveAs("b0_TT" + region + ".png")
    bcanvas.SaveAs("b0_TT" + region + ".root")

    bZeroFrameMM = mm_mass.frame(Range(jpsimin,jpsimax))
    bZeroFitData.plotOn(bZeroFrameMM)
    bZeroFrameMM.Draw()
    bcanvas.SaveAs("b0_MM" + region + ".png")
    bcanvas.SaveAs("b0_MM" + region + ".root")


if args.nofit and args.nofitb0:

    jcanvas = TCanvas("jcanvas","jcanvas",1400,800)

    jPsiFitData = theData.Clone("fitB0Data")

    fitjPsimin = 3.0
    fitjPsimax = 3.2

    if binnedfit:
        mm_mass.setBins(100)
        jPsiFitData = theData.binnedClone("binnedTrakData")

    mean = RooRealVar("m","m",3.09,3.06,3.1);
    sigma = RooRealVar("#sigma","#sigma",0.01,0.001,0.1);
    # sigma2 = RooRealVar("#sigma_{2}","#sigma_{2}",0.004,0.004,0.01);

    # c0 = RooRealVar("p0","p0",0.001,-5.,5.)
    # c1 = RooRealVar("p1","p1",0.001,-4.,4.)
    # c2 = RooRealVar("p2","p2",-0.00001,-4.,4.)
    # c3 = RooRealVar("p3","p3",-0.000001,-4.,4.)
    # c4 = RooRealVar("p4","p4",-0.000001,-2.,2.)
    # c5 = RooRealVar("p5","p5",-0.000001)
    # c6 = RooRealVar("p6","p6",-0.000001,-0.01,0.01)
    #
    # polyset = RooArgList(c0,c1,c2,c3,c4)

    # gFraMMKK = RooRealVar("f_{gauss}","f_{gauss}",0.3,0.0,1.0)
    # nSigMMKK = RooRealVar("n_{sig}","n_{sig}",10000,0.,10E6)
    # nBkgMMKK = RooRealVar("n_{bkg}","n_{bkg}",10000,0.,10E6)
    #
    # pdf_bkg = RooChebychev("pdf_bkg","pdf_bkg",mm_mass,polyset)
    # pdf_sig1 = RooGaussian("pdf_sig1","pdf_sig1",mm_mass,mean,sigma1)
    # pdf_sig2 = RooGaussian("pdf_sig2","pdf_sig2",mm_mass,mean,sigma2)
    #
    # pdf_sig  = RooAddPdf("pdf_sig","pdf_sig",pdf_sig1,pdf_sig2,gFraMMKK)
    # pdf_tot = RooAddPdf("pdf_tot","pdf_tot",RooArgList(pdf_sig1,pdf_bkg),RooArgList(nSigMMKK,nBkgMMKK))

    pdf_tot = RooGaussian("pdf_sig1","pdf_sig1",mm_mass,mean,sigma)

    nfit = 0

    jPsiFit = pdf_tot.fitTo(jPsiFitData,Range(fitjPsimin+0.005,fitjPsimax-0.005),RooFit.PrintLevel(-1), RooFit.NumCPU(numcpus),RooFit.Save())
    nfit += 1

    jPsiFrame = mm_mass.frame(Range(fitjPsimin+0.005,fitjPsimax-0.005))

    jPsiFitData.plotOn(jPsiFrame)

    pdf_tot.plotOn(jPsiFrame,RooFit.Normalization(1.0/float(nfit)))
    jPsiFitData.plotOn(jPsiFrame)
    pdf_tot.paramOn(jPsiFrame,RooFit.Layout(0.7,0.99,0.99))
    jPsiFrame.Draw()
    jcanvas.SaveAs("mm_fit" + region + ".png")
    jcanvas.SaveAs("mm_fit" + region + ".root")
