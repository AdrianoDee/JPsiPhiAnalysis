
//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul  2 00:35:06 2019 by ROOT version 6.12/07
// from TTree dimuonTree/Tree of DiMuon
// found on file: 0000.root
//////////////////////////////////////////////////////////

#ifndef Dimuon_h
#define Dimuon_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include "TLorentzVector.h"



class Dimuon : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   TTree *outTree;

   Float_t out_run, out_event, out_lumiblock, out_ndimuon, out_nmuons;
    Float_t out_trigger, out_charge, out_tMatch, out_dimuon_p4, out_highMuon_p4;
    Float_t out_lowMuon_p4, out_dimuon_m, out_dimuon_pt, out_dimuon_eta, out_dimuon_phi;
    Float_t out_dimuon_p, out_MassErr, out_vProb, out_DCA, out_ppdlPV;
    Float_t out_ppdlErrPV, out_ppdlBS, out_ppdlErrBS, out_cosAlpha, out_lxyPV;
    Float_t out_lxyBS, out_numPrimaryVertices,

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
   TTreeReaderValue<TLorentzVector> highMuon_p4 = {fReader, "highMuon_p4"};
   TTreeReaderValue<TLorentzVector> lowMuon_p4 = {fReader, "lowMuon_p4"};
   TTreeReaderValue<Float_t> dimuon_m = {fReader, "dimuon_m"};
   TTreeReaderValue<Float_t> dimuon_pt = {fReader, "dimuon_pt"};
   TTreeReaderValue<Float_t> dimuon_eta = {fReader, "dimuon_eta"};
   TTreeReaderValue<Float_t> dimuon_phi = {fReader, "dimuon_phi"};
   TTreeReaderValue<Float_t> dimuon_p = {fReader, "dimuon_p"};
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


   Dimuon(TTree * /*tree*/ =0) { }
   virtual ~Dimuon() { }
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

   ClassDef(Dimuon,0);

   TProofOutputFile *OutFile;
   TFile            *fOut;

};

#endif

#ifdef Dimuon_cxx
void Dimuon::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t Dimuon::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef Dimuon_cxx
~
