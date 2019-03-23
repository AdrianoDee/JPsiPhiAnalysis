//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Mar 23 16:50:33 2019 by ROOT version 6.12/07
// from TTree JPsiPhiTree/Tree of DiMuon and DiTrack
// found on file: rootuple-2018-dimuondiTrack_fivedataD2018_2_0_1.root
//////////////////////////////////////////////////////////

#ifndef FourTracks_h
#define FourTracks_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include "TLorentzVector.h"



class FourTracks : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> run = {fReader, "run"};
   TTreeReaderValue<Int_t> event = {fReader, "event"};
   TTreeReaderValue<Int_t> lumi = {fReader, "lumi"};
   TTreeReaderValue<Int_t> numPrimaryVertices = {fReader, "numPrimaryVertices"};
   TTreeReaderValue<Int_t> trigger = {fReader, "trigger"};
   TTreeReaderValue<Int_t> noXCandidates = {fReader, "noXCandidates"};
   TTreeReaderValue<Int_t> dimuon_id = {fReader, "dimuon_id"};
   TTreeReaderValue<Int_t> p_id = {fReader, "p_id"};
   TTreeReaderValue<Int_t> m_id = {fReader, "m_id"};
   TTreeReaderValue<TLorentzVector> dimuonditrk_p4 = {fReader, "dimuonditrk_p4"};
   TTreeReaderValue<TLorentzVector> ditrack_p4 = {fReader, "ditrack_p4"};
   TTreeReaderValue<TLorentzVector> dimuon_p4 = {fReader, "dimuon_p4"};
   TTreeReaderValue<TLorentzVector> lowMuon_p4 = {fReader, "lowMuon_p4"};
   TTreeReaderValue<TLorentzVector> highMuon_p4 = {fReader, "highMuon_p4"};
   TTreeReaderValue<TLorentzVector> highKaon_p4 = {fReader, "highKaon_p4"};
   TTreeReaderValue<TLorentzVector> lowKaon_p4 = {fReader, "lowKaon_p4"};
   TTreeReaderValue<TLorentzVector> dimuonditrk_rf_p4 = {fReader, "dimuonditrk_rf_p4"};
   TTreeReaderValue<TLorentzVector> ditrack_rf_p4 = {fReader, "ditrack_rf_p4"};
   TTreeReaderValue<TLorentzVector> dimuon_rf_p4 = {fReader, "dimuon_rf_p4"};
   TTreeReaderValue<TLorentzVector> lowMuon_rf_p4 = {fReader, "lowMuon_rf_p4"};
   TTreeReaderValue<TLorentzVector> highMuon_rf_p4 = {fReader, "highMuon_rf_p4"};
   TTreeReaderValue<TLorentzVector> kaonp_rf_p4 = {fReader, "kaonp_rf_p4"};
   TTreeReaderValue<TLorentzVector> kaonn_rf_p4 = {fReader, "kaonn_rf_p4"};
   TTreeReaderValue<Double_t> dimuonditrk_m = {fReader, "dimuonditrk_m"};
   TTreeReaderValue<Double_t> dimuonditrk_m_rf = {fReader, "dimuonditrk_m_rf"};
   TTreeReaderValue<Double_t> dimuonditrk_m_rf_d_c = {fReader, "dimuonditrk_m_rf_d_c"};
   TTreeReaderValue<Double_t> dimuonditrk_m_rf_c = {fReader, "dimuonditrk_m_rf_c"};
   TTreeReaderValue<Double_t> dimuonditrk_pt = {fReader, "dimuonditrk_pt"};
   TTreeReaderValue<Double_t> dimuonditrk_eta = {fReader, "dimuonditrk_eta"};
   TTreeReaderValue<Double_t> dimuonditrk_phi = {fReader, "dimuonditrk_phi"};
   TTreeReaderValue<Double_t> dimuonditrk_y = {fReader, "dimuonditrk_y"};
   TTreeReaderValue<Double_t> dimuonditrk_vx = {fReader, "dimuonditrk_vx"};
   TTreeReaderValue<Double_t> dimuonditrk_vy = {fReader, "dimuonditrk_vy"};
   TTreeReaderValue<Double_t> dimuonditrk_vz = {fReader, "dimuonditrk_vz"};
   TTreeReaderValue<Double_t> pv_x = {fReader, "pv_x"};
   TTreeReaderValue<Double_t> pv_y = {fReader, "pv_y"};
   TTreeReaderValue<Double_t> pv_z = {fReader, "pv_z"};
   TTreeReaderValue<Double_t> dimuon_m = {fReader, "dimuon_m"};
   TTreeReaderValue<Double_t> dimuon_pt = {fReader, "dimuon_pt"};
   TTreeReaderValue<Double_t> ditrack_m = {fReader, "ditrack_m"};
   TTreeReaderValue<Double_t> ditrack_pt = {fReader, "ditrack_pt"};
   TTreeReaderValue<Double_t> highKaon_pt = {fReader, "highKaon_pt"};
   TTreeReaderValue<Double_t> lowKaon_pt = {fReader, "lowKaon_pt"};
   TTreeReaderValue<Double_t> highMuon_pt = {fReader, "highMuon_pt"};
   TTreeReaderValue<Double_t> lowMuon_pt = {fReader, "lowMuon_pt"};
   TTreeReaderValue<Double_t> highKaon_eta = {fReader, "highKaon_eta"};
   TTreeReaderValue<Double_t> lowKaon_eta = {fReader, "lowKaon_eta"};
   TTreeReaderValue<Double_t> highMuon_eta = {fReader, "highMuon_eta"};
   TTreeReaderValue<Double_t> lowMuon_eta = {fReader, "lowMuon_eta"};
   TTreeReaderValue<Double_t> highKaon_phi = {fReader, "highKaon_phi"};
   TTreeReaderValue<Double_t> lowKaon_phi = {fReader, "lowKaon_phi"};
   TTreeReaderValue<Double_t> highMuon_phi = {fReader, "highMuon_phi"};
   TTreeReaderValue<Double_t> lowMuon_phi = {fReader, "lowMuon_phi"};
   TTreeReaderValue<Double_t> highKaon_dxy = {fReader, "highKaon_dxy"};
   TTreeReaderValue<Double_t> lowKaon_dxy = {fReader, "lowKaon_dxy"};
   TTreeReaderValue<Double_t> highMuon_dxy = {fReader, "highMuon_dxy"};
   TTreeReaderValue<Double_t> lowMuon_dxy = {fReader, "lowMuon_dxy"};
   TTreeReaderValue<Double_t> highKaon_dz = {fReader, "highKaon_dz"};
   TTreeReaderValue<Double_t> lowKaon_dz = {fReader, "lowKaon_dz"};
   TTreeReaderValue<Double_t> highMuon_dz = {fReader, "highMuon_dz"};
   TTreeReaderValue<Double_t> lowMuon_dz = {fReader, "lowMuon_dz"};
   TTreeReaderValue<Double_t> highKaonMuonDR = {fReader, "highKaonMuonDR"};
   TTreeReaderValue<Double_t> highKaonMuonDP = {fReader, "highKaonMuonDP"};
   TTreeReaderValue<Double_t> highKaonMuonDPt = {fReader, "highKaonMuonDPt"};
   TTreeReaderValue<Double_t> lowKaonMuonDR = {fReader, "lowKaonMuonDR"};
   TTreeReaderValue<Double_t> lowKaonMuonDP = {fReader, "lowKaonMuonDP"};
   TTreeReaderValue<Double_t> lowKaonMuonDPt = {fReader, "lowKaonMuonDPt"};
   TTreeReaderValue<Double_t> dimuonditrk_refPK_mass = {fReader, "dimuonditrk_refPK_mass"};
   TTreeReaderValue<Double_t> dimuonditrk_refKP_mass = {fReader, "dimuonditrk_refKP_mass"};
   TTreeReaderValue<Double_t> dimuonditrk_refPP_mass = {fReader, "dimuonditrk_refPP_mass"};
   TTreeReaderValue<Double_t> dimuonditrk_refPK_vChi2 = {fReader, "dimuonditrk_refPK_vChi2"};
   TTreeReaderValue<Double_t> dimuonditrk_refKP_vChi2 = {fReader, "dimuonditrk_refKP_vChi2"};
   TTreeReaderValue<Double_t> dimuonditrk_refPP_vChi2 = {fReader, "dimuonditrk_refPP_vChi2"};
   TTreeReaderValue<Double_t> dimuonditrk_refPK_nDof = {fReader, "dimuonditrk_refPK_nDof"};
   TTreeReaderValue<Double_t> dimuonditrk_refKP_nDof = {fReader, "dimuonditrk_refKP_nDof"};
   TTreeReaderValue<Double_t> dimuonditrk_refPP_nDof = {fReader, "dimuonditrk_refPP_nDof"};
   TTreeReaderValue<Double_t> dimuonditrk_refPK_vProb = {fReader, "dimuonditrk_refPK_vProb"};
   TTreeReaderValue<Double_t> dimuonditrk_refKP_vProb = {fReader, "dimuonditrk_refKP_vProb"};
   TTreeReaderValue<Double_t> dimuonditrk_refPP_vProb = {fReader, "dimuonditrk_refPP_vProb"};
   TTreeReaderValue<Double_t> dimuon_vProb = {fReader, "dimuon_vProb"};
   TTreeReaderValue<Double_t> dimuon_vChi2 = {fReader, "dimuon_vChi2"};
   TTreeReaderValue<Double_t> dimuon_DCA = {fReader, "dimuon_DCA"};
   TTreeReaderValue<Double_t> dimuon_ctauPV = {fReader, "dimuon_ctauPV"};
   TTreeReaderValue<Double_t> dimuon_ctauErrPV = {fReader, "dimuon_ctauErrPV"};
   TTreeReaderValue<Double_t> dimuon_cosAlpha = {fReader, "dimuon_cosAlpha"};
   TTreeReaderValue<Int_t> dimuon_triggerMatch = {fReader, "dimuon_triggerMatch"};
   TTreeReaderValue<Double_t> dimuonditrk_vProb = {fReader, "dimuonditrk_vProb"};
   TTreeReaderValue<Double_t> dimuonditrk_vChi2 = {fReader, "dimuonditrk_vChi2"};
   TTreeReaderValue<Double_t> dimuonditrk_nDof = {fReader, "dimuonditrk_nDof"};
   TTreeReaderValue<Int_t> dimuonditrk_charge = {fReader, "dimuonditrk_charge"};
   TTreeReaderValue<Double_t> dimuonditrk_cosAlpha = {fReader, "dimuonditrk_cosAlpha"};
   TTreeReaderValue<Double_t> dimuonditrk_ctauPV = {fReader, "dimuonditrk_ctauPV"};
   TTreeReaderValue<Double_t> dimuonditrk_ctauErrPV = {fReader, "dimuonditrk_ctauErrPV"};
   TTreeReaderValue<Double_t> dimuonditrk_tPFromPV = {fReader, "dimuonditrk_tPFromPV"};
   TTreeReaderValue<Double_t> dimuonditrk_tMFromPV = {fReader, "dimuonditrk_tMFromPV"};
   TTreeReaderValue<Double_t> dimuonditrk_cosAlphaCA = {fReader, "dimuonditrk_cosAlphaCA"};
   TTreeReaderValue<Double_t> dimuonditrk_ctauPVCA = {fReader, "dimuonditrk_ctauPVCA"};
   TTreeReaderValue<Double_t> dimuonditrk_ctauErrPVCA = {fReader, "dimuonditrk_ctauErrPVCA"};
   TTreeReaderValue<Double_t> dimuonditrk_tPFromPVCA = {fReader, "dimuonditrk_tPFromPVCA"};
   TTreeReaderValue<Double_t> dimuonditrk_tMFromPVCA = {fReader, "dimuonditrk_tMFromPVCA"};
   TTreeReaderValue<Double_t> dimuonditrk_cosAlphaDZ = {fReader, "dimuonditrk_cosAlphaDZ"};
   TTreeReaderValue<Double_t> dimuonditrk_ctauPVDZ = {fReader, "dimuonditrk_ctauPVDZ"};
   TTreeReaderValue<Double_t> dimuonditrk_ctauErrPVDZ = {fReader, "dimuonditrk_ctauErrPVDZ"};
   TTreeReaderValue<Double_t> dimuonditrk_tPFromPVDZ = {fReader, "dimuonditrk_tPFromPVDZ"};
   TTreeReaderValue<Double_t> dimuonditrk_tMFromPVDZ = {fReader, "dimuonditrk_tMFromPVDZ"};
   TTreeReaderValue<Double_t> dimuonditrk_cosAlphaBS = {fReader, "dimuonditrk_cosAlphaBS"};
   TTreeReaderValue<Double_t> dimuonditrk_ctauPVBS = {fReader, "dimuonditrk_ctauPVBS"};
   TTreeReaderValue<Double_t> dimuonditrk_ctauErrPVBS = {fReader, "dimuonditrk_ctauErrPVBS"};
   TTreeReaderValue<Double_t> dimuonditrk_dca_m1m2 = {fReader, "dimuonditrk_dca_m1m2"};
   TTreeReaderValue<Double_t> dimuonditrk_dca_m1t1 = {fReader, "dimuonditrk_dca_m1t1"};
   TTreeReaderValue<Double_t> dimuonditrk_dca_m1t2 = {fReader, "dimuonditrk_dca_m1t2"};
   TTreeReaderValue<Double_t> dimuonditrk_dca_m2t1 = {fReader, "dimuonditrk_dca_m2t1"};
   TTreeReaderValue<Double_t> dimuonditrk_dca_m2t2 = {fReader, "dimuonditrk_dca_m2t2"};
   TTreeReaderValue<Double_t> dimuonditrk_dca_t1t2 = {fReader, "dimuonditrk_dca_t1t2"};
   TTreeReaderValue<Double_t> dimuonditrk_rf_vProb = {fReader, "dimuonditrk_rf_vProb"};
   TTreeReaderValue<Double_t> dimuonditrk_rf_vChi2 = {fReader, "dimuonditrk_rf_vChi2"};
   TTreeReaderValue<Double_t> dimuonditrk_rf_nDof = {fReader, "dimuonditrk_rf_nDof"};
   TTreeReaderValue<Double_t> dimuonditrk_rf_cosAlpha = {fReader, "dimuonditrk_rf_cosAlpha"};
   TTreeReaderValue<Double_t> dimuonditrk_rf_ctauPV = {fReader, "dimuonditrk_rf_ctauPV"};
   TTreeReaderValue<Double_t> dimuonditrk_rf_ctauErrPV = {fReader, "dimuonditrk_rf_ctauErrPV"};
   TTreeReaderValue<Double_t> dimuonditrk_rf_c_vProb = {fReader, "dimuonditrk_rf_c_vProb"};
   TTreeReaderValue<Double_t> dimuonditrk_rf_c_vChi2 = {fReader, "dimuonditrk_rf_c_vChi2"};
   TTreeReaderValue<Double_t> dimuonditrk_rf_c_nDof = {fReader, "dimuonditrk_rf_c_nDof"};
   TTreeReaderValue<Double_t> dimuonditrk_rf_c_cosAlpha = {fReader, "dimuonditrk_rf_c_cosAlpha"};
   TTreeReaderValue<Double_t> dimuonditrk_rf_c_ctauPV = {fReader, "dimuonditrk_rf_c_ctauPV"};
   TTreeReaderValue<Double_t> dimuonditrk_rf_c_ctauErrPV = {fReader, "dimuonditrk_rf_c_ctauErrPV"};
   TTreeReaderValue<Int_t> highKaonMatch = {fReader, "highKaonMatch"};
   TTreeReaderValue<Int_t> lowKaonMatch = {fReader, "lowKaonMatch"};
   TTreeReaderValue<Int_t> lowMuonMatch = {fReader, "lowMuonMatch"};
   TTreeReaderValue<Int_t> highMuonMatch = {fReader, "highMuonMatch"};
   TTreeReaderValue<Bool_t> lowMuon_isTight = {fReader, "lowMuon_isTight"};
   TTreeReaderValue<Bool_t> lowMuon_isLoose = {fReader, "lowMuon_isLoose"};
   TTreeReaderValue<Bool_t> lowMuon_isSoft = {fReader, "lowMuon_isSoft"};
   TTreeReaderValue<Bool_t> lowMuon_isMedium = {fReader, "lowMuon_isMedium"};
   TTreeReaderValue<Bool_t> lowMuon_isHighPt = {fReader, "lowMuon_isHighPt"};
   TTreeReaderValue<Bool_t> lowMuon_isTracker = {fReader, "lowMuon_isTracker"};
   TTreeReaderValue<Bool_t> lowMuon_isGlobal = {fReader, "lowMuon_isGlobal"};
   TTreeReaderValue<Int_t> lowMuon_NPixelHits = {fReader, "lowMuon_NPixelHits"};
   TTreeReaderValue<Int_t> lowMuon_NStripHits = {fReader, "lowMuon_NStripHits"};
   TTreeReaderValue<Int_t> lowMuon_NTrackhits = {fReader, "lowMuon_NTrackhits"};
   TTreeReaderValue<Int_t> lowMuon_NBPixHits = {fReader, "lowMuon_NBPixHits"};
   TTreeReaderValue<Int_t> lowMuon_NPixLayers = {fReader, "lowMuon_NPixLayers"};
   TTreeReaderValue<Int_t> lowMuon_NTraLayers = {fReader, "lowMuon_NTraLayers"};
   TTreeReaderValue<Int_t> lowMuon_NStrLayers = {fReader, "lowMuon_NStrLayers"};
   TTreeReaderValue<Int_t> lowMuon_NBPixLayers = {fReader, "lowMuon_NBPixLayers"};
   TTreeReaderValue<Bool_t> highMuon_isTight = {fReader, "highMuon_isTight"};
   TTreeReaderValue<Bool_t> highMuon_isLoose = {fReader, "highMuon_isLoose"};
   TTreeReaderValue<Bool_t> highMuon_isSoft = {fReader, "highMuon_isSoft"};
   TTreeReaderValue<Bool_t> highMuon_isMedium = {fReader, "highMuon_isMedium"};
   TTreeReaderValue<Bool_t> highMuon_isHighPt = {fReader, "highMuon_isHighPt"};
   TTreeReaderValue<Bool_t> highMuon_isTracker = {fReader, "highMuon_isTracker"};
   TTreeReaderValue<Bool_t> highMuon_isGlobal = {fReader, "highMuon_isGlobal"};
   TTreeReaderValue<Int_t> highMuon_NPixelHits = {fReader, "highMuon_NPixelHits"};
   TTreeReaderValue<Int_t> highMuon_NStripHits = {fReader, "highMuon_NStripHits"};
   TTreeReaderValue<Int_t> highMuon_NTrackhits = {fReader, "highMuon_NTrackhits"};
   TTreeReaderValue<Int_t> highMuon_NBPixHits = {fReader, "highMuon_NBPixHits"};
   TTreeReaderValue<Int_t> highMuon_NPixLayers = {fReader, "highMuon_NPixLayers"};
   TTreeReaderValue<Int_t> highMuon_NTraLayers = {fReader, "highMuon_NTraLayers"};
   TTreeReaderValue<Int_t> highMuon_NStrLayers = {fReader, "highMuon_NStrLayers"};
   TTreeReaderValue<Int_t> highMuon_NBPixLayers = {fReader, "highMuon_NBPixLayers"};
   TTreeReaderValue<UInt_t> lowMuon_type = {fReader, "lowMuon_type"};
   TTreeReaderValue<UInt_t> highMuon_type = {fReader, "highMuon_type"};
   TTreeReaderValue<Int_t> highKaon_NPixelHits = {fReader, "highKaon_NPixelHits"};
   TTreeReaderValue<Int_t> highKaon_NStripHits = {fReader, "highKaon_NStripHits"};
   TTreeReaderValue<Int_t> highKaon_NTrackhits = {fReader, "highKaon_NTrackhits"};
   TTreeReaderValue<Int_t> highKaon_NBPixHits = {fReader, "highKaon_NBPixHits"};
   TTreeReaderValue<Int_t> highKaon_NPixLayers = {fReader, "highKaon_NPixLayers"};
   TTreeReaderValue<Int_t> highKaon_NTraLayers = {fReader, "highKaon_NTraLayers"};
   TTreeReaderValue<Int_t> highKaon_NStrLayers = {fReader, "highKaon_NStrLayers"};
   TTreeReaderValue<Int_t> highKaon_NBPixLayers = {fReader, "highKaon_NBPixLayers"};
   TTreeReaderValue<Int_t> lowKaon_NPixelHits = {fReader, "lowKaon_NPixelHits"};
   TTreeReaderValue<Int_t> lowKaon_NStripHits = {fReader, "lowKaon_NStripHits"};
   TTreeReaderValue<Int_t> lowKaon_NTrackhits = {fReader, "lowKaon_NTrackhits"};
   TTreeReaderValue<Int_t> lowKaon_NBPixHits = {fReader, "lowKaon_NBPixHits"};
   TTreeReaderValue<Int_t> lowKaon_NPixLayers = {fReader, "lowKaon_NPixLayers"};
   TTreeReaderValue<Int_t> lowKaon_NTraLayers = {fReader, "lowKaon_NTraLayers"};
   TTreeReaderValue<Int_t> lowKaon_NStrLayers = {fReader, "lowKaon_NStrLayers"};
   TTreeReaderValue<Int_t> lowKaon_NBPixLayers = {fReader, "lowKaon_NBPixLayers"};


   FourTracks(TTree * /*tree*/ =0) { }
   virtual ~FourTracks() { }
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

   ClassDef(FourTracks,0);

};

#endif

#ifdef FourTracks_cxx
void FourTracks::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t FourTracks::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef FourTracks_cxx
