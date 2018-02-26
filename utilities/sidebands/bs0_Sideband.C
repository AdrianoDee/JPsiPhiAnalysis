////////////////////////////////////////////////////
/// signal fit functions:	                 ///
/// GAUSSIAN 		          		 ///
/// VOIGTIAN 			 		 ///
/// BREIT-WIGNER				 /// 	 	
///				                 ///
/// background fit functions:	                 ///
/// ARGUS 					 ///
/// CHEBYCHEV 				         ///
/// POLYNOMIAL 					 ///	
/// EXPONENTIAL 	  			 ///
/// 						 ///
/// Binned Maximum Likelihood		         ///
////////////////////////////////////////////////////


#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooGaussian.h" 
#include "RooVoigtian.h"
#include "RooBreitWigner.h"
#include "RooArgusBG.h" 
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TFile.h"
#include "TBox.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLatex.h"

#include <TROOT.h>
#include <TH1.h>
#include <TF1.h>
#include <TF2.h>
#include <TFormula.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TString.h>
#include <TLine.h>
#include <TPad.h>
#include <TMath.h>
#include <TColor.h>
#include <iostream>

gROOT->Reset();
gROOT->Clear();

using namespace RooFit;

TStyle *myStyle = new TStyle("myStyle","myStyle");

void bs0_Sideband() 

{
  gROOT->SetStyle("Plain");
  gStyle->SetCanvasColor(0);
  gStyle->SetOptStat(10);

  Double_t chiSquare;  

  TFile *f1 = TFile::Open("Y4140_HLT4_Ours.root","read");
  //TFile *f1 = TFile::Open("Y4140_HLT4_Theirs.root","read");
  TH1F* hist = (TH1F*)f1->Get("Bs0_mass_Bs0NP");

  RooRealVar x("x","x",-10,10);
  RooDataHist* Bs0 = new RooDataHist("Bs0",hist->GetTitle(),RooArgSet(x),Import(*hist,kFALSE));
  RooPlot* xframe = x.frame(Title(""));
  xframe->SetTitle("B_{s}^{0} Invariant Mass");
  xframe->SetTitleOffset(1.32,"y");
  xframe->SetYTitle("Candidates / 4MeV/c^{2}");
  xframe->SetTitleOffset(1.26,"x");
  xframe->SetXTitle("M(J/#psi#phi) [Gev/c^{2}]");
  Bs0->plotOn(xframe);


  RooRealVar mean("mean","mean of gaussian",5.35,5.3,5.4);
  RooRealVar sigma1("#sigma_{1}","width of gaussian1", 0.029); // 0.1,0.01,0.5 // 0.029,0.0027,0.1 // 0.016,0.001,0.045
  RooGaussian gauss1("gauss1","Signal",x,mean,sigma1);
  RooRealVar sigma2("#sigma_{2}","width of gaussian2",0.02); // 0.1,0.02,0.5 // 0.052,0.047,0.1 // 0.02,0.001,0.090
  RooGaussian gauss2("gauss2","Signal",x,mean,sigma2); 
  RooRealVar sigfrac("sigfrac","fraction of the signal",0.5,0.,1.);
  RooAddPdf signal("signal","gauus1+gauss2",RooArgList(gauss1,gauss2),sigfrac); 
  RooRealVar a1("a1","a1",-0.1,-10.,10.); 
  RooRealVar a2("a2","a2",-0.1,-10.,10.); 
  RooChebychev background("cheby","Background",x,RooArgSet(a1,a2)); 
  RooRealVar nSig("nSig","Number of Signal Candidates",2e+5,1.,1e+6);
  RooRealVar nBkg("nBkg","Bkg Component",120e+3,1.,1e+6);


  RooAddPdf* totalPDF = new RooAddPdf("totalPDF","totalPDF",RooArgList(signal,background),RooArgList(nSig,nBkg));
  totalPDF->paramOn(xframe,Parameters(RooArgSet(mean,sigma1,sigma2,nSig,nBkg,sigfrac)),Layout(0.59,0.94,0.9));
 
  x.setRange("Signal",-2.0,2.0);
  totalPDF->fitTo(*Bs0,Range("Signal"),Save());
  totalPDF->plotOn(xframe,Range("Signal"),NormRange("Signal"));
 
  x.setRange("Sideband_Left",-6.0,-4.0);
  totalPDF->fitTo(*Bs0,Range("Sideband_Left"),Save());
  totalPDF->plotOn(xframe,Range("Sideband_Left"),NormRange("Sideband_Left"));

  x.SetRange("Sideband_Right",4.0,6.0); 
  totalPDF->fitTo(*Bs0,Range("Sideband_Right"),Save());
  totalPDF->plotOn(xframe,Range("Sideband_Right"),NormRange("Sideband_Right"));
 
  //totalPDF->fitTo(*Bs0,Range("Sideband_Right"),Save());
  //RooFitResult* SB2 = totalPDF.fitTo(*Bs0,Range("Sideband_Right"),Save());

  

 
  TCanvas* myc = new TCanvas("Sideband Subtraction","Sideband Subtraction",700,700);
  myc->SetFrameFillColor(0);
  myc->cd();
  gPad->SetPad(0.,0.3,1.,1.);
  xframe->Draw();
  myc->SaveAs("Sideband.png");

	


}



