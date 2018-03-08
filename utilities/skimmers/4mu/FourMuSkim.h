//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Mar  8 12:32:35 2018 by ROOT version 6.10/09
// from TTree OniaPhiTree/Tree of Onia and Phi
// found on file: /lustre/cms/store/user/adiflori/Charmonium/4mu_miniaod_17Nov2017_BCDEF_2017.root
//////////////////////////////////////////////////////////

#ifndef FourMuSkim_h
#define FourMuSkim_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include "TLorentzVector.h"



class FourMuSkim : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> run = {fReader, "run"};
   TTreeReaderValue<Int_t> event = {fReader, "event"};
   TTreeReaderValue<Int_t> numPrimaryVertices = {fReader, "numPrimaryVertices"};
   TTreeReaderValue<Int_t> trigger = {fReader, "trigger"};
   TTreeReaderValue<TLorentzVector> doubledimuon_p4 = {fReader, "doubledimuon_p4"};
   TTreeReaderValue<TLorentzVector> lowdim_p4 = {fReader, "lowdim_p4"};
   TTreeReaderValue<TLorentzVector> higdim_p4 = {fReader, "higdim_p4"};
   TTreeReaderValue<TLorentzVector> muonLowN_p4 = {fReader, "muonLowN_p4"};
   TTreeReaderValue<TLorentzVector> muonHighN_p4 = {fReader, "muonHighN_p4"};
   TTreeReaderValue<TLorentzVector> muonHighP_p4 = {fReader, "muonHighP_p4"};
   TTreeReaderValue<TLorentzVector> doubledimuon_rf_p4 = {fReader, "doubledimuon_rf_p4"};
   TTreeReaderValue<TLorentzVector> lowdim_rf_p4 = {fReader, "lowdim_rf_p4"};
   TTreeReaderValue<TLorentzVector> higdim_rf_p4 = {fReader, "higdim_rf_p4"};
   TTreeReaderValue<TLorentzVector> muonLowN_rf_p4 = {fReader, "muonLowN_rf_p4"};
   TTreeReaderValue<TLorentzVector> muonHighN_rf_p4 = {fReader, "muonHighN_rf_p4"};
   TTreeReaderValue<TLorentzVector> muonHighP_rf_p4 = {fReader, "muonHighP_rf_p4"};
   TTreeReaderValue<Int_t> doubledimuon_rf_bindx = {fReader, "doubledimuon_rf_bindx"};
   TTreeReaderValue<Double_t> higdim_vProb = {fReader, "higdim_vProb"};
   TTreeReaderValue<Double_t> higdim_vNChi2 = {fReader, "higdim_vNChi2"};
   TTreeReaderValue<Double_t> higdim_DCA = {fReader, "higdim_DCA"};
   TTreeReaderValue<Double_t> higdim_ctauPV = {fReader, "higdim_ctauPV"};
   TTreeReaderValue<Double_t> higdim_ctauErrPV = {fReader, "higdim_ctauErrPV"};
   TTreeReaderValue<Double_t> higdim_cosAlpha = {fReader, "higdim_cosAlpha"};
   TTreeReaderValue<Int_t> higdim_triggerMatch = {fReader, "higdim_triggerMatch"};
   TTreeReaderValue<Double_t> lowdim_vProb = {fReader, "lowdim_vProb"};
   TTreeReaderValue<Double_t> lowdim_vNChi2 = {fReader, "lowdim_vNChi2"};
   TTreeReaderValue<Double_t> lowdim_DCA = {fReader, "lowdim_DCA"};
   TTreeReaderValue<Double_t> lowdim_ctauPV = {fReader, "lowdim_ctauPV"};
   TTreeReaderValue<Double_t> lowdim_ctauErrPV = {fReader, "lowdim_ctauErrPV"};
   TTreeReaderValue<Double_t> lowdim_cosAlpha = {fReader, "lowdim_cosAlpha"};
   TTreeReaderValue<Int_t> lowdim_triggerMatch = {fReader, "lowdim_triggerMatch"};
   TTreeReaderValue<Double_t> doubledimuon_vProb = {fReader, "doubledimuon_vProb"};
   TTreeReaderValue<Double_t> doubledimuon_vChi2 = {fReader, "doubledimuon_vChi2"};
   TTreeReaderValue<Double_t> doubledimuon_cosAlpha = {fReader, "doubledimuon_cosAlpha"};
   TTreeReaderValue<Double_t> doubledimuon_ctauPV = {fReader, "doubledimuon_ctauPV"};
   TTreeReaderValue<Double_t> doubledimuon_ctauErrPV = {fReader, "doubledimuon_ctauErrPV"};
   TTreeReaderValue<Int_t> doubledimuon_charge = {fReader, "doubledimuon_charge"};
   TTreeReaderValue<Double_t> highDiMM_fit = {fReader, "highDiMM_fit"};
   TTreeReaderValue<Double_t> highDiMPx_fit = {fReader, "highDiMPx_fit"};
   TTreeReaderValue<Double_t> highDiMPy_fit = {fReader, "highDiMPy_fit"};
   TTreeReaderValue<Double_t> highDiMPz_fit = {fReader, "highDiMPz_fit"};
   TTreeReaderValue<Bool_t> muonHighP_isLoose = {fReader, "muonHighP_isLoose"};
   TTreeReaderValue<Bool_t> muonHighP_isSoft = {fReader, "muonHighP_isSoft"};
   TTreeReaderValue<Bool_t> muonHighP_isMedium = {fReader, "muonHighP_isMedium"};
   TTreeReaderValue<Bool_t> muonHighP_isHighPt = {fReader, "muonHighP_isHighPt"};
   TTreeReaderValue<Bool_t> muonHighP_isTracker = {fReader, "muonHighP_isTracker"};
   TTreeReaderValue<Bool_t> muonHighP_isGlobal = {fReader, "muonHighP_isGlobal"};
   TTreeReaderValue<Bool_t> muonHighN_isLoose = {fReader, "muonHighN_isLoose"};
   TTreeReaderValue<Bool_t> muonHighN_isSoft = {fReader, "muonHighN_isSoft"};
   TTreeReaderValue<Bool_t> muonHighN_isMedium = {fReader, "muonHighN_isMedium"};
   TTreeReaderValue<Bool_t> muonHighN_isHighPt = {fReader, "muonHighN_isHighPt"};
   TTreeReaderValue<Bool_t> muonHighN_isTracker = {fReader, "muonHighN_isTracker"};
   TTreeReaderValue<Bool_t> muonHighN_isGlobal = {fReader, "muonHighN_isGlobal"};
   TTreeReaderValue<UInt_t> muonHighP_type = {fReader, "muonHighP_type"};
   TTreeReaderValue<UInt_t> muonHighN_type = {fReader, "muonHighN_type"};
   TTreeReaderValue<Bool_t> muonLowP_isLoose = {fReader, "muonLowP_isLoose"};
   TTreeReaderValue<Bool_t> muonLowP_isSoft = {fReader, "muonLowP_isSoft"};
   TTreeReaderValue<Bool_t> muonLowP_isMedium = {fReader, "muonLowP_isMedium"};
   TTreeReaderValue<Bool_t> muonLowP_isHighPt = {fReader, "muonLowP_isHighPt"};
   TTreeReaderValue<Bool_t> muonLowP_isTracker = {fReader, "muonLowP_isTracker"};
   TTreeReaderValue<Bool_t> muonLowP_isGlobal = {fReader, "muonLowP_isGlobal"};
   TTreeReaderValue<Bool_t> muonLowN_isLoose = {fReader, "muonLowN_isLoose"};
   TTreeReaderValue<Bool_t> muonLowN_isSoft = {fReader, "muonLowN_isSoft"};
   TTreeReaderValue<Bool_t> muonLowN_isMedium = {fReader, "muonLowN_isMedium"};
   TTreeReaderValue<Bool_t> muonLowN_isHighPt = {fReader, "muonLowN_isHighPt"};
   TTreeReaderValue<Bool_t> muonLowN_isTracker = {fReader, "muonLowN_isTracker"};
   TTreeReaderValue<Bool_t> muonLowN_isGlobal = {fReader, "muonLowN_isGlobal"};
   TTreeReaderValue<UInt_t> muonLowP_type = {fReader, "muonLowP_type"};
   TTreeReaderValue<UInt_t> muonLowN_type = {fReader, "muonLowN_type"};
   TTreeReaderValue<Bool_t> isBestCandidate = {fReader, "isBestCandidate"};


   FourMuSkim(TTree * /*tree*/ =0) { }
   virtual ~FourMuSkim() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(FourMuSkim,0);

};

#endif

#ifdef FourMuSkim_cxx
void FourMuSkim::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t FourMuSkim::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef FourMuSkim_cxx
