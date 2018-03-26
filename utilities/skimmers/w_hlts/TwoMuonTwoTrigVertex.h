//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar 26 12:23:05 2018 by ROOT version 5.34/28
// from TTree dimuonditrkTree/Tree of DiMuonDiTrak
// found on file: /lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2k2Trig_Charmonium_Run2017F-17Nov2017-v1_MINIAOD_305388-309000__20180323_184238/180323_174253/0000/rootuple-2017-ditraktrigger_1.root
//////////////////////////////////////////////////////////

#ifndef TwoMuonTwoTrigVertex_h
#define TwoMuonTwoTrigVertex_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.
#include "/mnt/build/jenkins/workspace/root-release/BUILDTYPE/Release/COMPILER/native/LABEL/slc6/sources/root_v5.34.28/root/math/physics/inc/TLorentzVector.h"

// Fixed size dimensions of array or collections stored in the TTree if any.

class TwoMuonTwoTrigVertex : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   UInt_t          run;
   ULong64_t       event;
   UInt_t          lumiblock;
   UInt_t          ndimuonditrk;
   UInt_t          trigger;
   Int_t           charge;
   Bool_t          isBest;
   TLorentzVector  *dimuonditrak_p4;
   TLorentzVector  *dimuon_p4;
   TLorentzVector  *muonP_p4;
   TLorentzVector  *muonN_p4;
   TLorentzVector  *ditrak_p4;
   TLorentzVector  *trakP_p4;
   TLorentzVector  *trakN_p4;
   TLorentzVector  *dimuonditrkTrigger_p4;
   TLorentzVector  *ditrakTrigger_p4;
   TLorentzVector  *dimuonTrigger_p4;
   Double_t        dimuonditrk_vProb;
   Double_t        dimuonditrk_vChi2;
   Double_t        dimuonditrk_cosAlpha;
   Double_t        dimuonditrk_ctauPV;
   Double_t        dimuonditrk_ctauErrPV;
   Int_t           dimuonditrk_charge;
   Float_t         dimuonditrk_lxy;
   Float_t         dimuonditrk_lxyErr;
   Float_t         dimuonditrk_MassErr;
   Double_t        dimuon_vProb;
   Double_t        dimuon_vNChi2;
   Double_t        dimuon_DCA;
   Double_t        dimuon_ctauPV;
   Double_t        dimuon_ctauErrPV;
   Double_t        dimuon_cosAlpha;
   Int_t           muonP_tMatch;
   Int_t           muonN_tMatch;
   Int_t           trakP_tMatch;
   Int_t           trakN_tMatch;
   UInt_t          numPrimaryVertices;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumiblock;   //!
   TBranch        *b_ndimuonditrk;   //!
   TBranch        *b_trigger;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_isBest;   //!
   TBranch        *b_dimuonditrak_p4;   //!
   TBranch        *b_dimuon_p4;   //!
   TBranch        *b_muonP_p4;   //!
   TBranch        *b_muonN_p4;   //!
   TBranch        *b_ditrak_p4;   //!
   TBranch        *b_trakP_p4;   //!
   TBranch        *b_trakN_p4;   //!
   TBranch        *b_dimuonditrkTrigger_p4;   //!
   TBranch        *b_ditrakTrigger_p4;   //!
   TBranch        *b_dimuonTrigger_p4;   //!
   TBranch        *b_dimuonditrk_vProb;   //!
   TBranch        *b_dimuonditrk_vChi2;   //!
   TBranch        *b_dimuonditrk_cosAlpha;   //!
   TBranch        *b_dimuonditrk_ctauPV;   //!
   TBranch        *b_dimuonditrk_ctauErrPV;   //!
   TBranch        *b_dimuonditrk_charge;   //!
   TBranch        *b_dimuonditrk_charge;   //!
   TBranch        *b_dimuonditrk_lxyErr;   //!
   TBranch        *b_dimuonditrk_MassErr;   //!
   TBranch        *b_dimuon_vProb;   //!
   TBranch        *b_dimuon_vNChi2;   //!
   TBranch        *b_dimuon_DCA;   //!
   TBranch        *b_dimuon_ctauPV;   //!
   TBranch        *b_dimuon_ctauErrPV;   //!
   TBranch        *b_dimuon_cosAlpha;   //!
   TBranch        *b_muonP_tMatch;   //!
   TBranch        *b_muonN_tMatch;   //!
   TBranch        *b_trakP_tMatch;   //!
   TBranch        *b_trakN_tMatch;   //!
   TBranch        *b_numPrimaryVertices;   //!

   TwoMuonTwoTrigVertex(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~TwoMuonTwoTrigVertex() { }
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

   ClassDef(TwoMuonTwoTrigVertex,0);
};

#endif

#ifdef TwoMuonTwoTrigVertex_cxx
void TwoMuonTwoTrigVertex::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   dimuonditrak_p4 = 0;
   dimuon_p4 = 0;
   muonP_p4 = 0;
   muonN_p4 = 0;
   ditrak_p4 = 0;
   trakP_p4 = 0;
   trakN_p4 = 0;
   dimuonditrkTrigger_p4 = 0;
   ditrakTrigger_p4 = 0;
   dimuonTrigger_p4 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
   fChain->SetBranchAddress("ndimuonditrk", &ndimuonditrk, &b_ndimuonditrk);
   fChain->SetBranchAddress("trigger", &trigger, &b_trigger);
   fChain->SetBranchAddress("charge", &charge, &b_charge);
   fChain->SetBranchAddress("isBest", &isBest, &b_isBest);
   fChain->SetBranchAddress("dimuonditrak_p4", &dimuonditrak_p4, &b_dimuonditrak_p4);
   fChain->SetBranchAddress("dimuon_p4", &dimuon_p4, &b_dimuon_p4);
   fChain->SetBranchAddress("muonP_p4", &muonP_p4, &b_muonP_p4);
   fChain->SetBranchAddress("muonN_p4", &muonN_p4, &b_muonN_p4);
   fChain->SetBranchAddress("ditrak_p4", &ditrak_p4, &b_ditrak_p4);
   fChain->SetBranchAddress("trakP_p4", &trakP_p4, &b_trakP_p4);
   fChain->SetBranchAddress("trakN_p4", &trakN_p4, &b_trakN_p4);
   fChain->SetBranchAddress("dimuonditrkTrigger_p4", &dimuonditrkTrigger_p4, &b_dimuonditrkTrigger_p4);
   fChain->SetBranchAddress("ditrakTrigger_p4", &ditrakTrigger_p4, &b_ditrakTrigger_p4);
   fChain->SetBranchAddress("dimuonTrigger_p4", &dimuonTrigger_p4, &b_dimuonTrigger_p4);
   fChain->SetBranchAddress("dimuonditrk_vProb", &dimuonditrk_vProb, &b_dimuonditrk_vProb);
   fChain->SetBranchAddress("dimuonditrk_vChi2", &dimuonditrk_vChi2, &b_dimuonditrk_vChi2);
   fChain->SetBranchAddress("dimuonditrk_cosAlpha", &dimuonditrk_cosAlpha, &b_dimuonditrk_cosAlpha);
   fChain->SetBranchAddress("dimuonditrk_ctauPV", &dimuonditrk_ctauPV, &b_dimuonditrk_ctauPV);
   fChain->SetBranchAddress("dimuonditrk_ctauErrPV", &dimuonditrk_ctauErrPV, &b_dimuonditrk_ctauErrPV);
   fChain->SetBranchAddress("dimuonditrk_charge", &dimuonditrk_charge, &b_dimuonditrk_charge);
   fChain->SetBranchAddress("dimuonditrk_lxy", &dimuonditrk_lxy, &b_dimuonditrk_charge);
   fChain->SetBranchAddress("dimuonditrk_lxyErr", &dimuonditrk_lxyErr, &b_dimuonditrk_lxyErr);
   fChain->SetBranchAddress("dimuonditrk_MassErr", &dimuonditrk_MassErr, &b_dimuonditrk_MassErr);
   fChain->SetBranchAddress("dimuon_vProb", &dimuon_vProb, &b_dimuon_vProb);
   fChain->SetBranchAddress("dimuon_vNChi2", &dimuon_vNChi2, &b_dimuon_vNChi2);
   fChain->SetBranchAddress("dimuon_DCA", &dimuon_DCA, &b_dimuon_DCA);
   fChain->SetBranchAddress("dimuon_ctauPV", &dimuon_ctauPV, &b_dimuon_ctauPV);
   fChain->SetBranchAddress("dimuon_ctauErrPV", &dimuon_ctauErrPV, &b_dimuon_ctauErrPV);
   fChain->SetBranchAddress("dimuon_cosAlpha", &dimuon_cosAlpha, &b_dimuon_cosAlpha);
   fChain->SetBranchAddress("muonP_tMatch", &muonP_tMatch, &b_muonP_tMatch);
   fChain->SetBranchAddress("muonN_tMatch", &muonN_tMatch, &b_muonN_tMatch);
   fChain->SetBranchAddress("trakP_tMatch", &trakP_tMatch, &b_trakP_tMatch);
   fChain->SetBranchAddress("trakN_tMatch", &trakN_tMatch, &b_trakN_tMatch);
   fChain->SetBranchAddress("numPrimaryVertices", &numPrimaryVertices, &b_numPrimaryVertices);
}

Bool_t TwoMuonTwoTrigVertex::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef TwoMuonTwoTrigVertex_cxx
