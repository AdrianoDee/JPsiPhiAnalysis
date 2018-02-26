#!/usr/bin/env python
import sys, os, os.path, re
import commands,string,getopt
import argparse

from math import *

from ROOT import TString
from ROOT import TFile, TCanvas
import ROOT

#usage="\n Usage: python JPsiByRun.py <options> \n Options: \n --trigger= \t\t trigger path like: Jpsi6p5 Jpsi6p5_Barrel \n --period \t\t reco period like: CertMay6  May10 \n  --selection \t\t selection like: Barrel Barrel_nocuts Barrel_vtx Barrel_muon Barrel_trigger"

x4140jsonfile = "/afs/cern.ch/work/a/adiflori/CMSSW_5_3_22/src/X4140/MuMuKKPAT/test/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON_MuonPhys.txt"

parser = argparse.ArgumentParser()
parser.add_argument('--json', type=str, default=x4140jsonfile,
                    help='input json file for brilcalc')
parser.add_argument('--readbril', type=str, default=None,help='brilcalc lumi result to use as input')
parser.add_argument('--histoname', type=str, default="JPsi_vs_run",help='input no. of JPsi vs run histogram name')
parser.add_argument('--input', type=str, default="JPsiCount.root",help='input root file name')
parser.add_argument('--hlts', nargs='+', help='Select different hlts suffixes (also a list, e.g. --hlts HLT8 HLT4)', required=False,default=None)
parser.add_argument('--outbril', type=str, default="all_lumis.txt",help='brilcalc output txt')

args = parser.parse_args()

filename  = args.input
jsonfile  = args.json
histoname = args.histoname
outbril   = args.outbril
inpbril   = args.readbril

if args.readbril is None:
    bril_cmd      = 'brilcalc lumi -i ' + jsonfile +' > ' + outbril
    bril_cmd_test = 'brilcalc lumi --help'
    inpbril = outbril
else:
    bril_cmd      = 'echo' #dummy comand
    bril_cmd_test = 'echo' #dummy comand

fileIn   = ROOT.TFile(filename)
binmap = {}

# ###############################
# Settin a bit of style
def tdrStyle():
  ROOT.gStyle.SetCanvasBorderMode(0)
  ROOT.gStyle.SetCanvasColor(0)
  ROOT.gStyle.SetCanvasDefH(600) #Height of canvas
  ROOT.gStyle.SetCanvasDefW(600) #Width of canvas
  ROOT.gStyle.SetCanvasDefX(0)   #POsition on screen
  ROOT.gStyle.SetCanvasDefY(0)
  # For the Pad:
  ROOT.gStyle.SetPadBorderMode(0)
  # ROOT.gStyle.SetPadBorderSize(Width_t size = 1)
  ROOT.gStyle.SetPadColor(0)
  ROOT.gStyle.SetPadGridX(False)
  ROOT.gStyle.SetPadGridY(False)
  ROOT.gStyle.SetGridColor(0)
  ROOT.gStyle.SetGridStyle(3)
  ROOT.gStyle.SetGridWidth(1)
  # For the frame:
  ROOT.gStyle.SetFrameBorderMode(0)
  ROOT.gStyle.SetFrameBorderSize(1)
  ROOT.gStyle.SetFrameFillColor(0)
  ROOT.gStyle.SetFrameFillStyle(0)
  ROOT.gStyle.SetFrameLineColor(1)
  ROOT.gStyle.SetFrameLineStyle(1)
  ROOT.gStyle.SetFrameLineWidth(1)
  # For the histo:
  ROOT.gStyle.SetHistLineColor(1)
  ROOT.gStyle.SetHistLineStyle(0)
  ROOT.gStyle.SetHistLineWidth(1)
  ROOT.gStyle.SetEndErrorSize(2)
  ROOT.gStyle.SetMarkerStyle(20)
  #For the date:
  ROOT.gStyle.SetOptDate(0)
  # For the statistics box:
  ROOT.gStyle.SetOptFile(0)
  #ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptStat("eou")
  ROOT.gStyle.SetStatColor(0)
  ROOT.gStyle.SetStatFont(42)
  ROOT.gStyle.SetStatFontSize(0.04)#/---> ROOT.gStyle.SetStatFontSize(0.025)
  ROOT.gStyle.SetStatTextColor(1)
  ROOT.gStyle.SetStatFormat("6.4g")
  ROOT.gStyle.SetStatBorderSize(1)
  ROOT.gStyle.SetStatH(0.15)
  ROOT.gStyle.SetStatW(0.3)#/---> ROOT.gStyle.SetStatW(0.15)
  # Margins:
  ROOT.gStyle.SetPadTopMargin(0.05)
  ROOT.gStyle.SetPadBottomMargin(0.13)
  ROOT.gStyle.SetPadLeftMargin(0.16)
  ROOT.gStyle.SetPadRightMargin(0.04)
  # For the Global title:
  ROOT.gStyle.SetOptTitle(0)
  # For the axis titles:
  ROOT.gStyle.SetTitleColor(1, "XYZ")
  ROOT.gStyle.SetTitleFont(42, "XYZ")
  ROOT.gStyle.SetTitleSize(0.06, "XYZ")
  ROOT.gStyle.SetTitleXOffset(0.9)
  ROOT.gStyle.SetTitleYOffset(1.25)
  # For the axis labels:
  ROOT.gStyle.SetLabelColor(1, "XYZ")
  ROOT.gStyle.SetLabelFont(42, "XYZ")
  ROOT.gStyle.SetLabelOffset(0.007, "XYZ")
  ROOT.gStyle.SetLabelSize(0.05, "XYZ")
  # For the axis:
  ROOT.gStyle.SetAxisColor(1, "XYZ")
  ROOT.gStyle.SetStripDecimals(True)
  ROOT.gStyle.SetTickLength(0.03, "XYZ")
  ROOT.gStyle.SetNdivisions(510, "XYZ")
  ROOT.gStyle.SetPadTickX(1)  # To get tick marks on the opposite side of the frame
  ROOT.gStyle.SetPadTickY(1)
  # Postscript options:
  ROOT.gStyle.SetPaperSize(20.,20.)
  ## OVERRIDES
  ROOT.gStyle.SetHistMinimumZero(1)
  ROOT.gStyle.SetErrorX(0.5)
  ROOT.gStyle.SetOptStat(False)
  # Done
  ROOT.gROOT.ForceStyle()

# ###############################
def testBrilEnv():
  status_cmd,lumiout=commands.getstatusoutput(bril_cmd_test)
  if status_cmd != 0 :
         print "\n ERROR with Environment for brilcalc ==> %s \n"%lumiout
         sys.exit(1)

# ###############################
def ReBinRunHisto(histo,name):

  ## find the new binning based on runs with muons
  nbinX = histo.GetNbinsX()
  nbinOutX =0
  first = True
  bin_first = int(round(histo.GetBinCenter(1),0))
  bin_last = int(round(histo.GetBinCenter(nbinX),0))
  for ix in range(nbinX):
    run = round(histo.GetBinCenter(ix),0)
    irun = int(run)
    if histo.GetBinContent(ix)>0 :  # skip runs without JPsi events
       if (first):
           bin_first = int(round(histo.GetBinCenter(ix),0))
           first = False
       bin_last = int(round(histo.GetBinCenter(ix),0))
       nbinOutX = nbinOutX + 1
       binmap[ix]=histo.GetBinContent(ix)

  print 'First Run = %s  Last Run = %s'%(bin_first,bin_last)
  histo_lumi = ROOT.TH1F(name, 'Nb. of JPsi vs run', nbinOutX,bin_first,bin_last )
  ixnew = 1
  ## fill the histogram with the new binning
  #for ix in range(nbinX):
  for ix, val in binmap.iteritems():
    if val>0 : # skip runs without JPsi events
       run = round(histo.GetBinCenter(ix),0)
       irun = int(run)
       #    print "====> run %i is now in bin %s"%(irun,ixnew)
       runcontent = val #histo.GetBinContent(ix)
       histo_lumi.SetBinContent(ixnew,runcontent)
       histo_lumi.SetBinError(ixnew,histo.GetBinError(ix))
       histo_lumi.GetXaxis().SetBinLabel(ixnew,str(irun))
       ixnew = ixnew + 1
  return histo_lumi

# ###############################
def GetMuPerLumiHisto(histo,hlt):
  nbinX = histo.GetNbinsX()
  histo.GetNbinsX()
#
  histo.Sumw2()
#
  histo_lumi = histo.Clone("lumiJPsirate" + hlt)

  status_cmd,lumiout=commands.getstatusoutput(bril_cmd)

  if status_cmd != 0 :
     print "ERROR %s"%lumiout
     sys.exit(2)


  #for ix in range(1,nbinX):
  for ix, val in binmap.iteritems():
   nmuperlumi = 0.
   lumi=999999.
   #run = round(histo.GetBinCenter(ix),0)
   run = round(histo.GetBinCenter(ix),0)
   irun = int(run) - 1
   #if histo.GetBinContent(ix)>0 : # for runs with JPsi
   if val>0:
      #runcontent = histo.GetBinContent(ix)
      runcontent = val
      grep_cmd='cat ' + str(inpbril) + ' | grep '+str(irun)+': | grep -v Warning | grep -v WARNING | cut -d"|" -f7'
      status_cmd,lumiout=commands.getstatusoutput(grep_cmd)

      if lumiout == '':
          print irun
          continue

      lumiclean = string.strip(lumiout)
      lumiclean = lumiout.strip(" ")
      lumiclean = lumiclean.strip("\n")
      lumiclean = lumiclean.strip("\b")
      lumiclean = lumiclean.strip("\t")
      #
    #   print(lumiclean)
    #   print(status_cmd)
      if status_cmd != 0:
         print "ERROR %s for run %i"%(lumiout,irun)
      else:
         #print(lumiclean)
         lumi =float(string.strip(lumiclean)) / 1000.

         nmuperlumi = runcontent/lumi
         err=sqrt(runcontent)/lumi

         print("run=%i #JPsi=%i lumi=%s  #JPsi per lumi=%s +/- %s"%(irun, runcontent,lumi,nmuperlumi,err))

   histo_lumi.SetBinContent(ix,nmuperlumi)
   histo_lumi.SetBinError(ix,histo.GetBinError(ix)/lumi)
  return histo_lumi
# ###############################
# ###############################
def printRunHisto(histo,title):
    histo.SetXTitle("Run")
    histo.SetYTitle(title)
    histo.SetMarkerStyle(20)
    histo.Draw("E1")


################################
#
################################
if __name__ == '__main__':
  ROOT.gBenchmark.Start('run info')
  tdrStyle()

  if args.hlts is not None:
      hlts = args.hlts
  else:
      hlts = [""]

  hfile = TFile( 'rootfiles/jpsiLumiNormRunII.root', 'RECREATE', 'ROOT file with histograms' )

  for hlt in hlts:
      fileIn.cd()
      histo = fileIn.Get(histoname + hlt)
      ## #mu vs run
      c1 = TCanvas("c1","c1",1200, 900)
      c1.SetLogy(1)
      printRunHisto(histo,'Nb. JPsi ')
      c1.SaveAs("plots/JPsivsrun%s.eps"%(hlt))


      c2 = TCanvas("c2","c2",1200, 900)
      c2.SetLogy(1)
      histo_rebin=ReBinRunHisto(histo,'JPsi_rebin' + hlt)
      printRunHisto(histo_rebin,'Nb. JPsi ')
      c2.SaveAs("plots/JPsivsrunrebin%s.eps"%(hlt))

      testBrilEnv()
      histo_lumi=GetMuPerLumiHisto(histo,hlt)
      cl = TCanvas("cl","cl",1200, 900)
      cl.SetLogy(0)
      printRunHisto(histo_lumi,'# JPsi (from X)/lumi (mb)' + hlt)
      cl.SaveAs("plots/JPsivsrun_lumi%s.eps"%(hlt))

      clr = TCanvas("clr","clr",1200, 900)
      clr.SetLogy(0)
      histo_lumi_rebin=ReBinRunHisto(histo_lumi,'JPsilumi_rebin' + hlt)
      printRunHisto(histo_lumi_rebin,'# JPsi (from X)/lumi (mb)')
      clr.SaveAs("plots/JPsivsrunrebin_lumi%s.eps"%(hlt))

      hfile.cd()
      histo.Write()
      histo_rebin.Write()
      histo_lumi.Write()
      histo_lumi_rebin.Write()


  ROOT.gBenchmark.Show( 'run info' )
