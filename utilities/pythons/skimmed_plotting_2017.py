
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
parser.add_argument('--input', type=str, default="2mu2k_tree.root")
# parser.add_argument('--binned', action='store_true')
parser.add_argument('--nonprompt', action='store_true')
parser.add_argument('--prompt', action='store_true')
# parser.add_argument('--nofit', action='store_true')
# parser.add_argument('--nofitkk', action='store_false')
# parser.add_argument('--nofitb0', action='store_false')
parser.add_argument('--noplot', action='store_false')
# parser.add_argument('--notrig', action='store_false')
# parser.add_argument('--noreco', action='store_false')
# parser.add_argument('--dosidebands', action='store_false')
parser.add_argument('--ptcuts', type=float, default=None)

#                    help='number of epochs')
#parser.add_argument('--batch_size', type=int, default=64)
#
args = parser.parse_args()

debugging = args.debug

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

phimin = 0.9
phimax = 1.3

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

lxysig = RooRealVar("xL","l(xy) sign.;l(xy) sign.",0.0,50.0)

massvars    = RooArgList(tt_mass, mm_mass, mmtt_mass)
ptvars      = RooArgList(tt_pt, mm_pt, mmtt_pt)
extravars   = RooArgList(lxysig)

theData = RooDataSet("theData","theData",theTre,RooArgSet(massvars,ptvars,extravars))


c = TCanvas("canvas","canvas",1200,800)

if args.nonprompt:
    theData = theData.reduce("xL > 3.0")
    region = "_nprompt_"
if args.prompt:
    theData = theData.reduce("xL < 2.0")
    region = "_promt_"

if args.ptcuts is not None:
    theData = theData.reduce("trigp_pT > " + str(args.ptcuts))
    theData = theData.reduce("trign_pT > " + str(args.ptcuts))
    cuts += "pt_" + str(args.ptcuts) + "_"
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

    hist_tt_pt.SetTitle("Pt(KK)" + region + " ;M(KK) [GeV]; no. entries")
    hist_tt_pt.SetFillColor(ROOT.kBlue)
    hist_tt_pt.SetLineColor(ROOT.kBlue)
    hist_tt_pt.SetFillStyle(3003)
    hist_tt_pt.Draw()
    c.SaveAs("tt_pt_histo" + region + ".png")
    c.SaveAs("tt_pt_histo" + region + ".root")

    hist_mm_pt.SetTitle("Pt(#mu#mu)" + region + " ;Pt(#mu#mu) [GeV]; no. entries")
    hist_mm_pt.SetFillColor(ROOT.kGreen)
    hist_mm_pt.SetLineColor(ROOT.kGreen)
    hist_mm_pt.SetFillStyle(3003)
    hist_mm_pt.Draw()
    c.SaveAs("mm_pt_histo" + region + ".png")
    c.SaveAs("mm_pt_histo" + region + ".root")

    hist_tt_pt.SetTitle("Pt(#mu#muKK)" + region + " ;Pt(#mu#muKK) [GeV]; no. entries")
    hist_tt_pt.SetFillColor(ROOT.kRed)
    hist_tt_pt.SetLineColor(ROOT.kRed)
    hist_tt_pt.SetFillStyle(3003)
    hist_tt_pt.Draw()
    c.SaveAs("tt_pt_histo" + region + ".png")
    c.SaveAs("tt_pt_histo" + region + ".root")
