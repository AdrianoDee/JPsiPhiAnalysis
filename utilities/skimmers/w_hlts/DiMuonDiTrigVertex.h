//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 23 16:52:11 2018 by ROOT version 6.10/09
// from TTree dimuonditrkTree/Tree of DiMuonDiTrak
// found on file: /lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2k2Trig_Charmonium_Run2017F-17Nov2017-v1_MINIAOD_305388-309000__20180322_150955/merge.root
//////////////////////////////////////////////////////////

#ifndef DiMuonDiTrigVertex_h
#define DiMuonDiTrigVertex_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include "TLorentzVector.h"

#include <TNtuple.h>
#include <TString.h>
#include <TSelector.h>

#include <TProof.h>
#include <TProofOutputFile.h>



class DiMuonDiTrigVertex : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   TNtuple *outTuple;

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<UInt_t> run = {fReader, "run"};
   TTreeReaderValue<ULong64_t> event = {fReader, "event"};
   TTreeReaderValue<UInt_t> lumiblock = {fReader, "lumiblock"};
   TTreeReaderValue<UInt_t> ndimuonditrk = {fReader, "ndimuonditrk"};
   TTreeReaderValue<UInt_t> trigger = {fReader, "trigger"};
   TTreeReaderValue<Int_t> charge = {fReader, "charge"};
   TTreeReaderValue<Bool_t> isBest = {fReader, "isBest"};
   TTreeReaderValue<TLorentzVector> dimuonditrak_p4 = {fReader, "dimuonditrak_p4"};
   TTreeReaderValue<TLorentzVector> dimuon_p4 = {fReader, "dimuon_p4"};
   TTreeReaderValue<TLorentzVector> muonP_p4 = {fReader, "muonP_p4"};
   TTreeReaderValue<TLorentzVector> muonN_p4 = {fReader, "muonN_p4"};
   TTreeReaderValue<TLorentzVector> ditrak_p4 = {fReader, "ditrak_p4"};
   TTreeReaderValue<TLorentzVector> trakP_p4 = {fReader, "trakP_p4"};
   TTreeReaderValue<TLorentzVector> trakN_p4 = {fReader, "trakN_p4"};
   TTreeReaderValue<TLorentzVector> dimuonditrkTrigger_p4 = {fReader, "dimuonditrkTrigger_p4"};
   TTreeReaderValue<TLorentzVector> ditrakTrigger_p4 = {fReader, "ditrakTrigger_p4"};
   TTreeReaderValue<TLorentzVector> dimuonTrigger_p4 = {fReader, "dimuonTrigger_p4"};
   TTreeReaderValue<Double_t> dimuonditrk_vProb = {fReader, "dimuonditrk_vProb"};
   TTreeReaderValue<Double_t> dimuonditrk_vChi2 = {fReader, "dimuonditrk_vChi2"};
   TTreeReaderValue<Double_t> dimuonditrk_cosAlpha = {fReader, "dimuonditrk_cosAlpha"};
   TTreeReaderValue<Double_t> dimuonditrk_ctauPV = {fReader, "dimuonditrk_ctauPV"};
   TTreeReaderValue<Double_t> dimuonditrk_ctauErrPV = {fReader, "dimuonditrk_ctauErrPV"};
   TTreeReaderValue<Int_t> dimuonditrk_charge = {fReader, "dimuonditrk_charge"};
   TTreeReaderValue<Float_t> dimuonditrk_charge = {fReader, "dimuonditrk_lxy"};
   TTreeReaderValue<Float_t> dimuonditrk_lxyErr = {fReader, "dimuonditrk_lxyErr"};
   TTreeReaderValue<Float_t> dimuonditrk_MassErr = {fReader, "dimuonditrk_MassErr"};
   TTreeReaderValue<Double_t> dimuon_vProb = {fReader, "dimuon_vProb"};
   TTreeReaderValue<Double_t> dimuon_vNChi2 = {fReader, "dimuon_vNChi2"};
   TTreeReaderValue<Double_t> dimuon_DCA = {fReader, "dimuon_DCA"};
   TTreeReaderValue<Double_t> dimuon_ctauPV = {fReader, "dimuon_ctauPV"};
   TTreeReaderValue<Double_t> dimuon_ctauErrPV = {fReader, "dimuon_ctauErrPV"};
   TTreeReaderValue<Double_t> dimuon_cosAlpha = {fReader, "dimuon_cosAlpha"};
   TTreeReaderValue<Int_t> muonP_tMatch = {fReader, "muonP_tMatch"};
   TTreeReaderValue<Int_t> muonN_tMatch = {fReader, "muonN_tMatch"};
   TTreeReaderValue<Int_t> trakP_tMatch = {fReader, "trakP_tMatch"};
   TTreeReaderValue<Int_t> trakN_tMatch = {fReader, "trakN_tMatch"};
   TTreeReaderValue<UInt_t> numPrimaryVertices = {fReader, "numPrimaryVertices"};


   DiMuonDiTrigVertex(TTree * /*tree*/ =0) { }
   virtual ~DiMuonDiTrigVertex() { }
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

   TProofOutputFile *OutFile;
   TFile            *fOut;

   ClassDef(DiMuonDiTrigVertex,0);

};

#endif

#ifdef DiMuonDiTrigVertex_cxx
void DiMuonDiTrigVertex::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t DiMuonDiTrigVertex::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef DiMuonDiTrigVertex_cxx
