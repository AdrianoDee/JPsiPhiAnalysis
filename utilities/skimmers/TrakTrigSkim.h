//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb 28 10:51:08 2018 by ROOT version 6.10/09
// from TTree DiTrakDiTrigTree/Tree of ditrakditrig
// found on file: /lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2k2Trig_Charmonium_Run2017F-17Nov2017-v1_MINIAOD_305388-309000__20180227_203341/merge.root
//////////////////////////////////////////////////////////

#ifndef TrakTrigSkim_h
#define TrakTrigSkim_h

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

#include <bitset>

class TrakTrigSkim : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   Float_t JPsi_mass = 0.0, Phi_mass = 0.0, Phi_mean = 0.0, Phi_sigma = 0.0;
   TNtuple *outTuple;

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<UInt_t> run = {fReader, "run"};
   TTreeReaderValue<ULong64_t> event = {fReader, "event"};
   TTreeReaderValue<UInt_t> lumiblock = {fReader, "lumiblock"};
   TTreeReaderValue<UInt_t> nditrak = {fReader, "nditrak"};
   TTreeReaderValue<UInt_t> ntraks = {fReader, "ntraks"};
   TTreeReaderValue<UInt_t> trigger = {fReader, "trigger"};
   TTreeReaderValue<Int_t> charge = {fReader, "charge"};
   TTreeReaderValue<Int_t> tMatchOne = {fReader, "tMatchOne"};
   TTreeReaderValue<Int_t> tMatchTwo = {fReader, "tMatchTwo"};
   TTreeReaderValue<TLorentzVector> ditrak_p4 = {fReader, "ditrak_p4"};
   TTreeReaderValue<TLorentzVector> ditrig_p4 = {fReader, "ditrig_p4"};
   TTreeReaderValue<TLorentzVector> trigP_p4 = {fReader, "trigP_p4"};
   TTreeReaderValue<TLorentzVector> trigN_p4 = {fReader, "trigN_p4"};
   TTreeReaderValue<UInt_t> numPrimaryVertices = {fReader, "numPrimaryVertices"};


   TrakTrigSkim(TTree * /*tree*/ =0) { }
   virtual ~TrakTrigSkim() { }
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

   ClassDef(TrakTrigSkim,0);

};

#endif

#ifdef TrakTrigSkim_cxx
void TrakTrigSkim::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t TrakTrigSkim::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef TrakTrigSkim_cxx
