//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Feb 27 16:36:42 2018 by ROOT version 5.34/28
// from TTree DiMuon KPi Tree/Tree of DiMuon and DiTrak
// found on file: 2mukpi_BCDEF_2017.root
//////////////////////////////////////////////////////////

#ifndef TwoMuKPiSkim_h
#define TwoMuKPiSkim_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.
#include "TLorentzVector.h"

// Fixed size dimensions of array or collections stored in the TTree if any.

class TwoMuKPiSkim : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           numPrimaryVertices;
   Int_t           trigger;
   Int_t           noXCandidates;
   TLorentzVector  *dimuonditrk_p4;
   TLorentzVector  *ditrak_p4;
   TLorentzVector  *dimuon_p4;
   TLorentzVector  *muonp_p4;
   TLorentzVector  *muonn_p4;
   TLorentzVector  *kaonp_p4;
   TLorentzVector  *kaonn_p4;
   TLorentzVector  *dimuonditrk_rf_p4;
   TLorentzVector  *ditrak_rf_p4;
   TLorentzVector  *dimuon_rf_p4;
   TLorentzVector  *muonp_rf_p4;
   TLorentzVector  *muonn_rf_p4;
   TLorentzVector  *kaonp_rf_p4;
   TLorentzVector  *kaonn_rf_p4;
   Double_t        dimuon_vProb;
   Double_t        dimuon_vNChi2;
   Double_t        dimuon_DCA;
   Double_t        dimuon_ctauPV;
   Double_t        dimuon_ctauErrPV;
   Double_t        dimuon_cosAlpha;
   Int_t           dimuon_triggerMatch;
   Double_t        dimuonditrk_vProb;
   Double_t        dimuonditrk_vChi2;
   Double_t        dimuonditrk_cosAlpha;
   Double_t        dimuonditrk_ctauPV;
   Double_t        dimuonditrk_ctauErrPV;
   Int_t           dimuonditrk_charge;
   Bool_t          muonP_isLoose;
   Bool_t          muonP_isSoft;
   Bool_t          muonP_isMedium;
   Bool_t          muonP_isHighPt;
   Bool_t          muonP_isTracker;
   Bool_t          muonP_isGlobal;
   Bool_t          muonN_isLoose;
   Bool_t          muonN_isSoft;
   Bool_t          muonN_isMedium;
   Bool_t          muonN_isHighPt;
   Bool_t          muonN_isTracker;
   Bool_t          muonN_isGlobal;
   UInt_t          muonP_type;
   UInt_t          muonN_type;
   Bool_t          isBestCandidate;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_numPrimaryVertices;   //!
   TBranch        *b_trigger;   //!
   TBranch        *b_noXCandidates;   //!
   TBranch        *b_dimuonditrk_p4;   //!
   TBranch        *b_ditrak_p4;   //!
   TBranch        *b_dimuon_p4;   //!
   TBranch        *b_muonp_p4;   //!
   TBranch        *b_muonn_p4;   //!
   TBranch        *b_kaonp_p4;   //!
   TBranch        *b_kaonn_p4;   //!
   TBranch        *b_dimuonditrk_rf_p4;   //!
   TBranch        *b_ditrak_rf_p4;   //!
   TBranch        *b_dimuon_rf_p4;   //!
   TBranch        *b_muonp_rf_p4;   //!
   TBranch        *b_muonn_rf_p4;   //!
   TBranch        *b_kaonp_rf_p4;   //!
   TBranch        *b_kaonn_rf_p4;   //!
   TBranch        *b_dimuon_vProb;   //!
   TBranch        *b_dimuon_vNChi2;   //!
   TBranch        *b_dimuon_DCA;   //!
   TBranch        *b_dimuon_ctauPV;   //!
   TBranch        *b_dimuon_ctauErrPV;   //!
   TBranch        *b_dimuon_cosAlpha;   //!
   TBranch        *b_dimuon_triggerMatch;   //!
   TBranch        *b_dimuonditrk_vProb;   //!
   TBranch        *b_dimuonditrk_vChi2;   //!
   TBranch        *b_dimuonditrk_cosAlpha;   //!
   TBranch        *b_dimuonditrk_ctauPV;   //!
   TBranch        *b_dimuonditrk_ctauErrPV;   //!
   TBranch        *b_dimuonditrk_charge;   //!
   TBranch        *b_muonP_isLoose;   //!
   TBranch        *b_muonP_isSoft;   //!
   TBranch        *b_muonP_isMedium;   //!
   TBranch        *b_muonP_isHighPt;   //!
   TBranch        *b_muonP_isTracker;   //!
   TBranch        *b_muonP_isGlobal;   //!
   TBranch        *b_muonN_isLoose;   //!
   TBranch        *b_muonN_isSoft;   //!
   TBranch        *b_muonN_isMedium;   //!
   TBranch        *b_muonN_isHighPt;   //!
   TBranch        *b_muonN_isTracker;   //!
   TBranch        *b_muonN_isGlobal;   //!
   TBranch        *b_muonP_type;   //!
   TBranch        *b_muonN_type;   //!
   TBranch        *b_isBestCandidate;   //!

   TwoMuKPiSkim(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~TwoMuKPiSkim() { }
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

   ClassDef(TwoMuKPiSkim,0);
};

#endif

#ifdef TwoMuKPiSkim_cxx
void TwoMuKPiSkim::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   dimuonditrk_p4 = 0;
   ditrak_p4 = 0;
   dimuon_p4 = 0;
   muonp_p4 = 0;
   muonn_p4 = 0;
   kaonp_p4 = 0;
   kaonn_p4 = 0;
   dimuonditrk_rf_p4 = 0;
   ditrak_rf_p4 = 0;
   dimuon_rf_p4 = 0;
   muonp_rf_p4 = 0;
   muonn_rf_p4 = 0;
   kaonp_rf_p4 = 0;
   kaonn_rf_p4 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("numPrimaryVertices", &numPrimaryVertices, &b_numPrimaryVertices);
   fChain->SetBranchAddress("trigger", &trigger, &b_trigger);
   fChain->SetBranchAddress("noXCandidates", &noXCandidates, &b_noXCandidates);
   fChain->SetBranchAddress("dimuonditrk_p4", &dimuonditrk_p4, &b_dimuonditrk_p4);
   fChain->SetBranchAddress("ditrak_p4", &ditrak_p4, &b_ditrak_p4);
   fChain->SetBranchAddress("dimuon_p4", &dimuon_p4, &b_dimuon_p4);
   fChain->SetBranchAddress("muonp_p4", &muonp_p4, &b_muonp_p4);
   fChain->SetBranchAddress("muonn_p4", &muonn_p4, &b_muonn_p4);
   fChain->SetBranchAddress("kaonp_p4", &kaonp_p4, &b_kaonp_p4);
   fChain->SetBranchAddress("kaonn_p4", &kaonn_p4, &b_kaonn_p4);
   fChain->SetBranchAddress("dimuonditrk_rf_p4", &dimuonditrk_rf_p4, &b_dimuonditrk_rf_p4);
   fChain->SetBranchAddress("ditrak_rf_p4", &ditrak_rf_p4, &b_ditrak_rf_p4);
   fChain->SetBranchAddress("dimuon_rf_p4", &dimuon_rf_p4, &b_dimuon_rf_p4);
   fChain->SetBranchAddress("muonp_rf_p4", &muonp_rf_p4, &b_muonp_rf_p4);
   fChain->SetBranchAddress("muonn_rf_p4", &muonn_rf_p4, &b_muonn_rf_p4);
   fChain->SetBranchAddress("kaonp_rf_p4", &kaonp_rf_p4, &b_kaonp_rf_p4);
   fChain->SetBranchAddress("kaonn_rf_p4", &kaonn_rf_p4, &b_kaonn_rf_p4);
   fChain->SetBranchAddress("dimuon_vProb", &dimuon_vProb, &b_dimuon_vProb);
   fChain->SetBranchAddress("dimuon_vNChi2", &dimuon_vNChi2, &b_dimuon_vNChi2);
   fChain->SetBranchAddress("dimuon_DCA", &dimuon_DCA, &b_dimuon_DCA);
   fChain->SetBranchAddress("dimuon_ctauPV", &dimuon_ctauPV, &b_dimuon_ctauPV);
   fChain->SetBranchAddress("dimuon_ctauErrPV", &dimuon_ctauErrPV, &b_dimuon_ctauErrPV);
   fChain->SetBranchAddress("dimuon_cosAlpha", &dimuon_cosAlpha, &b_dimuon_cosAlpha);
   fChain->SetBranchAddress("dimuon_triggerMatch", &dimuon_triggerMatch, &b_dimuon_triggerMatch);
   fChain->SetBranchAddress("dimuonditrk_vProb", &dimuonditrk_vProb, &b_dimuonditrk_vProb);
   fChain->SetBranchAddress("dimuonditrk_vChi2", &dimuonditrk_vChi2, &b_dimuonditrk_vChi2);
   fChain->SetBranchAddress("dimuonditrk_cosAlpha", &dimuonditrk_cosAlpha, &b_dimuonditrk_cosAlpha);
   fChain->SetBranchAddress("dimuonditrk_ctauPV", &dimuonditrk_ctauPV, &b_dimuonditrk_ctauPV);
   fChain->SetBranchAddress("dimuonditrk_ctauErrPV", &dimuonditrk_ctauErrPV, &b_dimuonditrk_ctauErrPV);
   fChain->SetBranchAddress("dimuonditrk_charge", &dimuonditrk_charge, &b_dimuonditrk_charge);
   fChain->SetBranchAddress("muonP_isLoose", &muonP_isLoose, &b_muonP_isLoose);
   fChain->SetBranchAddress("muonP_isSoft", &muonP_isSoft, &b_muonP_isSoft);
   fChain->SetBranchAddress("muonP_isMedium", &muonP_isMedium, &b_muonP_isMedium);
   fChain->SetBranchAddress("muonP_isHighPt", &muonP_isHighPt, &b_muonP_isHighPt);
   fChain->SetBranchAddress("muonP_isTracker", &muonP_isTracker, &b_muonP_isTracker);
   fChain->SetBranchAddress("muonP_isGlobal", &muonP_isGlobal, &b_muonP_isGlobal);
   fChain->SetBranchAddress("muonN_isLoose", &muonN_isLoose, &b_muonN_isLoose);
   fChain->SetBranchAddress("muonN_isSoft", &muonN_isSoft, &b_muonN_isSoft);
   fChain->SetBranchAddress("muonN_isMedium", &muonN_isMedium, &b_muonN_isMedium);
   fChain->SetBranchAddress("muonN_isHighPt", &muonN_isHighPt, &b_muonN_isHighPt);
   fChain->SetBranchAddress("muonN_isTracker", &muonN_isTracker, &b_muonN_isTracker);
   fChain->SetBranchAddress("muonN_isGlobal", &muonN_isGlobal, &b_muonN_isGlobal);
   fChain->SetBranchAddress("muonP_type", &muonP_type, &b_muonP_type);
   fChain->SetBranchAddress("muonN_type", &muonN_type, &b_muonN_type);
   fChain->SetBranchAddress("isBestCandidate", &isBestCandidate, &b_isBestCandidate);
}

Bool_t TwoMuKPiSkim::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef TwoMuKPiSkim_cxx
