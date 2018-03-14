//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 14 18:46:03 2018 by ROOT version 5.34/28
// from TTree dimuonTree/Tree of DiMuon
// found on file: crab_miniaod_2mu2k_Charmonium_Run2016D-07Aug17-v1_MINIAOD___20180311_091533/merge.root
//////////////////////////////////////////////////////////

#ifndef JPsiCount_h
#define JPsiCount_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.
#include "/mnt/build/jenkins/workspace/root-release/BUILDTYPE/Release/COMPILER/native/LABEL/slc6/sources/root_v5.34.28/root/math/physics/inc/TLorentzVector.h"

// Fixed size dimensions of array or collections stored in the TTree if any.

class JPsiCount : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   UInt_t          run;
   ULong64_t       event;
   UInt_t          lumiblock;
   UInt_t          ndimuon;
   UInt_t          nmuons;
   UInt_t          trigger;
   Int_t           charge;
   Int_t           tMatch;
   TLorentzVector  *dimuon_p4;
   TLorentzVector  *muonP_p4;
   TLorentzVector  *muonN_p4;
   Float_t         MassErr;
   Float_t         vProb;
   Float_t         DCA;
   Float_t         ppdlPV;
   Float_t         ppdlErrPV;
   Float_t         ppdlBS;
   Float_t         ppdlErrBS;
   Float_t         cosAlpha;
   Float_t         lxyPV;
   Float_t         lxyBS;
   UInt_t          numPrimaryVertices;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumiblock;   //!
   TBranch        *b_ndimuon;   //!
   TBranch        *b_nmuons;   //!
   TBranch        *b_trigger;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_tMatch;   //!
   TBranch        *b_dimuon_p4;   //!
   TBranch        *b_muonP_p4;   //!
   TBranch        *b_muonN_p4;   //!
   TBranch        *b_MassErr;   //!
   TBranch        *b_vProb;   //!
   TBranch        *b_DCA;   //!
   TBranch        *b_ppdlPV;   //!
   TBranch        *b_ppdlErrPV;   //!
   TBranch        *b_ppdlBS;   //!
   TBranch        *b_ppdlErrBS;   //!
   TBranch        *b_cosAlpha;   //!
   TBranch        *b_lxyPV;   //!
   TBranch        *b_lxyBS;   //!
   TBranch        *b_numPrimaryVertices;   //!

   JPsiCount(TTree * /*tree*/ =0) : fChain(0) { }
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
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   dimuon_p4 = 0;
   muonP_p4 = 0;
   muonN_p4 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
   fChain->SetBranchAddress("ndimuon", &ndimuon, &b_ndimuon);
   fChain->SetBranchAddress("nmuons", &nmuons, &b_nmuons);
   fChain->SetBranchAddress("trigger", &trigger, &b_trigger);
   fChain->SetBranchAddress("charge", &charge, &b_charge);
   fChain->SetBranchAddress("tMatch", &tMatch, &b_tMatch);
   fChain->SetBranchAddress("dimuon_p4", &dimuon_p4, &b_dimuon_p4);
   fChain->SetBranchAddress("muonP_p4", &muonP_p4, &b_muonP_p4);
   fChain->SetBranchAddress("muonN_p4", &muonN_p4, &b_muonN_p4);
   fChain->SetBranchAddress("MassErr", &MassErr, &b_MassErr);
   fChain->SetBranchAddress("vProb", &vProb, &b_vProb);
   fChain->SetBranchAddress("DCA", &DCA, &b_DCA);
   fChain->SetBranchAddress("ppdlPV", &ppdlPV, &b_ppdlPV);
   fChain->SetBranchAddress("ppdlErrPV", &ppdlErrPV, &b_ppdlErrPV);
   fChain->SetBranchAddress("ppdlBS", &ppdlBS, &b_ppdlBS);
   fChain->SetBranchAddress("ppdlErrBS", &ppdlErrBS, &b_ppdlErrBS);
   fChain->SetBranchAddress("cosAlpha", &cosAlpha, &b_cosAlpha);
   fChain->SetBranchAddress("lxyPV", &lxyPV, &b_lxyPV);
   fChain->SetBranchAddress("lxyBS", &lxyBS, &b_lxyBS);
   fChain->SetBranchAddress("numPrimaryVertices", &numPrimaryVertices, &b_numPrimaryVertices);
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
