//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Feb 27 17:26:52 2018 by ROOT version 6.10/09
// from TTree DiMuon PiK Tree/Tree of DiMuon and DiTrak
// found on file: 2mukpi_BCDEF_2017.root
//////////////////////////////////////////////////////////

#ifndef TwoMuPiKSkim_h
#define TwoMuPiKSkim_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include "TLorentzVector.h"



class TwoMuPiKSkim : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> run = {fReader, "run"};
   TTreeReaderValue<Int_t> event = {fReader, "event"};
   TTreeReaderValue<Int_t> numPrimaryVertices = {fReader, "numPrimaryVertices"};
   TTreeReaderValue<Int_t> trigger = {fReader, "trigger"};
   TTreeReaderValue<Int_t> noXCandidates = {fReader, "noXCandidates"};
   TTreeReaderValue<TLorentzVector> dimuonditrk_p4 = {fReader, "dimuonditrk_p4"};
   TTreeReaderValue<TLorentzVector> ditrak_p4 = {fReader, "ditrak_p4"};
   TTreeReaderValue<TLorentzVector> dimuon_p4 = {fReader, "dimuon_p4"};
   TTreeReaderValue<TLorentzVector> muonp_p4 = {fReader, "muonp_p4"};
   TTreeReaderValue<TLorentzVector> muonn_p4 = {fReader, "muonn_p4"};
   TTreeReaderValue<TLorentzVector> kaonp_p4 = {fReader, "kaonp_p4"};
   TTreeReaderValue<TLorentzVector> kaonn_p4 = {fReader, "kaonn_p4"};
   TTreeReaderValue<TLorentzVector> dimuonditrk_rf_p4 = {fReader, "dimuonditrk_rf_p4"};
   TTreeReaderValue<TLorentzVector> ditrak_rf_p4 = {fReader, "ditrak_rf_p4"};
   TTreeReaderValue<TLorentzVector> dimuon_rf_p4 = {fReader, "dimuon_rf_p4"};
   TTreeReaderValue<TLorentzVector> muonp_rf_p4 = {fReader, "muonp_rf_p4"};
   TTreeReaderValue<TLorentzVector> muonn_rf_p4 = {fReader, "muonn_rf_p4"};
   TTreeReaderValue<TLorentzVector> kaonp_rf_p4 = {fReader, "kaonp_rf_p4"};
   TTreeReaderValue<TLorentzVector> kaonn_rf_p4 = {fReader, "kaonn_rf_p4"};
   TTreeReaderValue<Double_t> dimuon_vProb = {fReader, "dimuon_vProb"};
   TTreeReaderValue<Double_t> dimuon_vNChi2 = {fReader, "dimuon_vNChi2"};
   TTreeReaderValue<Double_t> dimuon_DCA = {fReader, "dimuon_DCA"};
   TTreeReaderValue<Double_t> dimuon_ctauPV = {fReader, "dimuon_ctauPV"};
   TTreeReaderValue<Double_t> dimuon_ctauErrPV = {fReader, "dimuon_ctauErrPV"};
   TTreeReaderValue<Double_t> dimuon_cosAlpha = {fReader, "dimuon_cosAlpha"};
   TTreeReaderValue<Int_t> dimuon_triggerMatch = {fReader, "dimuon_triggerMatch"};
   TTreeReaderValue<Double_t> dimuonditrk_vProb = {fReader, "dimuonditrk_vProb"};
   TTreeReaderValue<Double_t> dimuonditrk_vChi2 = {fReader, "dimuonditrk_vChi2"};
   TTreeReaderValue<Double_t> dimuonditrk_cosAlpha = {fReader, "dimuonditrk_cosAlpha"};
   TTreeReaderValue<Double_t> dimuonditrk_ctauPV = {fReader, "dimuonditrk_ctauPV"};
   TTreeReaderValue<Double_t> dimuonditrk_ctauErrPV = {fReader, "dimuonditrk_ctauErrPV"};
   TTreeReaderValue<Int_t> dimuonditrk_charge = {fReader, "dimuonditrk_charge"};
   TTreeReaderValue<Bool_t> muonP_isLoose = {fReader, "muonP_isLoose"};
   TTreeReaderValue<Bool_t> muonP_isSoft = {fReader, "muonP_isSoft"};
   TTreeReaderValue<Bool_t> muonP_isMedium = {fReader, "muonP_isMedium"};
   TTreeReaderValue<Bool_t> muonP_isHighPt = {fReader, "muonP_isHighPt"};
   TTreeReaderValue<Bool_t> muonP_isTracker = {fReader, "muonP_isTracker"};
   TTreeReaderValue<Bool_t> muonP_isGlobal = {fReader, "muonP_isGlobal"};
   TTreeReaderValue<Bool_t> muonN_isLoose = {fReader, "muonN_isLoose"};
   TTreeReaderValue<Bool_t> muonN_isSoft = {fReader, "muonN_isSoft"};
   TTreeReaderValue<Bool_t> muonN_isMedium = {fReader, "muonN_isMedium"};
   TTreeReaderValue<Bool_t> muonN_isHighPt = {fReader, "muonN_isHighPt"};
   TTreeReaderValue<Bool_t> muonN_isTracker = {fReader, "muonN_isTracker"};
   TTreeReaderValue<Bool_t> muonN_isGlobal = {fReader, "muonN_isGlobal"};
   TTreeReaderValue<UInt_t> muonP_type = {fReader, "muonP_type"};
   TTreeReaderValue<UInt_t> muonN_type = {fReader, "muonN_type"};
   TTreeReaderValue<Bool_t> isBestCandidate = {fReader, "isBestCandidate"};


   TwoMuPiKSkim(TTree * /*tree*/ =0) { }
   virtual ~TwoMuPiKSkim() { }
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

   ClassDef(TwoMuPiKSkim,0);

};

#endif

#ifdef TwoMuPiKSkim_cxx
void TwoMuPiKSkim::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t TwoMuPiKSkim::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef TwoMuPiKSkim_cxx
