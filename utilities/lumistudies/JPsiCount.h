//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 14 18:56:57 2018 by ROOT version 6.10/09
// from TTree dimuonTree/Tree of DiMuon
// found on file: crab_miniaod_2mu2k_Charmonium_Run2016D-07Aug17-v1_MINIAOD___20180311_091533/merge.root
//////////////////////////////////////////////////////////

#ifndef JPsiCount_h
#define JPsiCount_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include "TLorentzVector.h"



class JPsiCount : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<UInt_t> run = {fReader, "run"};
   TTreeReaderValue<ULong64_t> event = {fReader, "event"};
   TTreeReaderValue<UInt_t> lumiblock = {fReader, "lumiblock"};
   TTreeReaderValue<UInt_t> ndimuon = {fReader, "ndimuon"};
   TTreeReaderValue<UInt_t> nmuons = {fReader, "nmuons"};
   TTreeReaderValue<UInt_t> trigger = {fReader, "trigger"};
   TTreeReaderValue<Int_t> charge = {fReader, "charge"};
   TTreeReaderValue<Int_t> tMatch = {fReader, "tMatch"};
   TTreeReaderValue<TLorentzVector> dimuon_p4 = {fReader, "dimuon_p4"};
   TTreeReaderValue<TLorentzVector> muonP_p4 = {fReader, "muonP_p4"};
   TTreeReaderValue<TLorentzVector> muonN_p4 = {fReader, "muonN_p4"};
   TTreeReaderValue<Float_t> MassErr = {fReader, "MassErr"};
   TTreeReaderValue<Float_t> vProb = {fReader, "vProb"};
   TTreeReaderValue<Float_t> DCA = {fReader, "DCA"};
   TTreeReaderValue<Float_t> ppdlPV = {fReader, "ppdlPV"};
   TTreeReaderValue<Float_t> ppdlErrPV = {fReader, "ppdlErrPV"};
   TTreeReaderValue<Float_t> ppdlBS = {fReader, "ppdlBS"};
   TTreeReaderValue<Float_t> ppdlErrBS = {fReader, "ppdlErrBS"};
   TTreeReaderValue<Float_t> cosAlpha = {fReader, "cosAlpha"};
   TTreeReaderValue<Float_t> lxyPV = {fReader, "lxyPV"};
   TTreeReaderValue<Float_t> lxyBS = {fReader, "lxyBS"};
   TTreeReaderValue<UInt_t> numPrimaryVertices = {fReader, "numPrimaryVertices"};


   JPsiCount(TTree * /*tree*/ =0) { }
   virtual ~JPsiCount() { }
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

   ClassDef(JPsiCount,0);

};

#endif

#ifdef JPsiCount_cxx
void JPsiCount::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t JPsiCount::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef JPsiCount_cxx
