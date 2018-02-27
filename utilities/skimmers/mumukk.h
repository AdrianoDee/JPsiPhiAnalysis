//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jun 18 14:15:53 2017 by ROOT version 5.34/26
// from TTree X_data/X(4140) Data
// found on file: MuOniaParked_Run2012C_MuMuKKPAT_merged0.root
//////////////////////////////////////////////////////////

#ifndef mumukk_h
#define mumukk_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

#include <TSystem.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TBranch.h>
//#include <TCint.h>
#include <TRandom.h>
#include <TMath.h>
#include <TDirectory.h>
#include "TEnv.h"
#include <TString.h>
#include <TSelector.h>
#include <TProof.h>
#include <TProofOutputFile.h>
#include <TLorentzVector.h>
#include "TPoint.h"
#include <TH1.h>
#include <TH2.h>
#include <TH2F.h>
#include <TF1.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <map>

// Header file for the classes stored in the TTree if any.
//#include "/build/bellenot/source/root_v5.34.26/root/cint/cint/lib/prec_stl/vector"
//#include "/build/bellenot/source/root_v5.34.26/root/cint/cint/lib/prec_stl/vector"
//#include "/build/bellenot/source/root_v5.34.26/root/cint/cint/lib/prec_stl/vector"
//#include "/build/bellenot/source/root_v5.34.26/root/cint/cint/lib/prec_stl/vector"
//#include "/build/bellenot/source/root_v5.34.26/root/cint/cint/lib/prec_stl/vector"
//#include "/build/bellenot/source/root_v5.34.26/root/cint/cint/lib/prec_stl/vector"

// Fixed size dimensions of array or collections stored in the TTree if any.

class mumukk : public TSelector {
public :
  ////////////////////////////////////////////////////////////
  //INPUT TREE
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   // Declaration of leaf types
   vector<unsigned int> *TrigRes;
   vector<string>  *TrigNames;
   vector<string>  *MatchTriggerNames;
   vector<unsigned int> *L1TrigRes;
   UInt_t          evtNum;
   UInt_t          runNum;
   UInt_t          lumiNum;
   UInt_t          priVtx_n;
   Float_t         priVtx_X;
   Float_t         priVtx_Y;
   Float_t         priVtx_Z;
   Float_t         priVtx_XE;
   Float_t         priVtx_YE;
   Float_t         priVtx_ZE;
   Float_t         priVtx_NormChi2;
   Float_t         priVtx_Chi2;
   Float_t         priVtx_CL;
   UInt_t          priVtx_tracks;
   Float_t         priVtx_tracksPtSq;
   vector<float>   *trNotRef;
   vector<float>   *trRef;
   vector<float>   *trackPx;
   vector<float>   *trackPy;
   vector<float>   *trackPz;
   vector<float>   *trackEnergy;
   vector<int>     *trackNDF;
   vector<int>     *trackPhits;
   vector<int>     *trackShits;
   vector<float>   *trackChi2;
   vector<float>   *trackD0;
   vector<float>   *trackD0Err;
   vector<float>   *trackCharge;
   vector<int>     *TrackHighPurity;
   vector<int>     *TrackTight;
   vector<float>   *trackfHits;
   vector<bool>    *trackFirstBarrel;
   vector<bool>    *trackFirstEndCap;
   vector<float>   *trackDzVtx;
   vector<float>   *trackDxyVtx;
   vector<double>  *tr_nsigdedx;
   vector<float>   *tr_dedx;
   vector<float>   *tr_dedxMass;
   vector<float>   *tr_theo;
   vector<float>   *tr_sigma;
   vector<float>   *tr_dedx_byHits;
   vector<float>   *tr_dedxErr_byHits;
   vector<int>     *tr_saturMeas_byHits;
   vector<int>     *tr_Meas_byHits;
   UInt_t          nMu;
   vector<float>   *muPx;
   vector<float>   *muPy;
   vector<float>   *muPz;
   vector<float>   *muCharge;
   vector<float>   *muD0;
   vector<float>   *muDz;
   vector<float>   *muChi2;
   vector<int>     *muNDF;
   vector<int>     *muPhits;
   vector<int>     *muShits;
   vector<int>     *muLayersTr;
   vector<int>     *muLayersPix;
   vector<float>   *muD0E;
   vector<float>   *muDzVtxErr;
   vector<unsigned int> *muKey;
   vector<bool>    *muIsGlobal;
   vector<bool>    *muIsPF;
   vector<int>     *muGlMuHits;
   vector<float>   *muGlChi2;
   vector<int>     *muGlNDF;
   vector<int>     *muGlMatchedStation;
   vector<float>   *muGlDzVtx;
   vector<float>   *muGlDxyVtx;
   vector<int>     *nMatchedStations;
   vector<int>     *muType;
   vector<int>     *muQual;
   vector<int>     *muTrack;
   vector<int>     *muNOverlap;
   vector<int>     *muNSharingSegWith;
   vector<float>   *mufHits;
   vector<bool>    *muFirstBarrel;
   vector<bool>    *muFirstEndCap;
   vector<float>   *muDzVtx;
   vector<float>   *muDxyVtx;
   UInt_t          nMuMu;
   vector<float>   *MuMuMass;
   vector<float>   *MuMuPx;
   vector<float>   *MuMuPy;
   vector<float>   *MuMuPz;
   vector<float>   *MuMuVtx_CL;
   vector<float>   *MuMuVtx_Chi2;
   vector<float>   *MuMuDecayVtx_X;
   vector<float>   *MuMuDecayVtx_Y;
   vector<float>   *MuMuDecayVtx_Z;
   vector<float>   *MuMuDecayVtx_XE;
   vector<float>   *MuMuDecayVtx_YE;
   vector<float>   *MuMuDecayVtx_ZE;
   vector<int>     *mu1Idx;
   vector<int>     *mu2Idx;
   vector<int>     *Xmu1Idx;
   vector<int>     *Xmu2Idx;
   vector<int>     *X_ka1Idx;
   vector<int>     *X_ka2Idx;
   vector<float>   *mu1Px_MuMu;
   vector<float>   *mu1Py_MuMu;
   vector<float>   *mu1Pz_MuMu;
   vector<float>   *mu1Chi2_MuMu;
   vector<int>     *mu1NDF_MuMu;
   vector<float>   *mu2Px_MuMu;
   vector<float>   *mu2Py_MuMu;
   vector<float>   *mu2Pz_MuMu;
   vector<float>   *mu2Chi2_MuMu;
   vector<int>     *mu2NDF_MuMu;
   vector<int>     *MuMuType;
   vector<bool>    *MuMuMuonTrigMatch;
   vector<int>     *PriVtxMuMuCorr_n;
   vector<float>   *PriVtxMuMuCorr_X;
   vector<float>   *PriVtxMuMuCorr_Y;
   vector<float>   *PriVtxMuMuCorr_Z;
   vector<double>  *PriVtxMuMuCorr_EX;
   vector<double>  *PriVtxMuMuCorr_EY;
   vector<double>  *PriVtxMuMuCorr_EZ;
   vector<float>   *PriVtxMuMuCorr_Chi2;
   vector<float>   *PriVtxMuMuCorr_CL;
   vector<int>     *PriVtxMuMuCorr_tracks;
   vector<int>     *nTrk_afterMuMu;
   UInt_t          nKK;
   vector<float>   *KKMass;
   vector<float>   *KKPx;
   vector<float>   *KKPy;
   vector<float>   *KKPz;
   vector<float>   *KKVtx_CL;
   vector<float>   *KKVtx_Chi2;
   vector<float>   *KKDecayVtx_X;
   vector<float>   *KKDecayVtx_Y;
   vector<float>   *KKDecayVtx_Z;
   vector<float>   *KKDecayVtx_XE;
   vector<float>   *KKDecayVtx_YE;
   vector<float>   *KKDecayVtx_ZE;
   vector<int>     *ka1Idx;
   vector<int>     *ka2Idx;
   vector<float>   *ka1Px_KK;
   vector<float>   *ka1Py_KK;
   vector<float>   *ka1Pz_KK;
   vector<float>   *ka1Chi2_KK;
   vector<int>     *ka1NDF_KK;
   vector<float>   *ka2Px_KK;
   vector<float>   *ka2Py_KK;
   vector<float>   *ka2Pz_KK;
   vector<float>   *ka2Chi2_KK;
   vector<int>     *ka2NDF_KK;
   vector<float>   *DR_MuMu_K1;
   vector<float>   *DR_MuMu_K2;
   vector<float>   *DR_MuMuKK_K1;
   vector<float>   *DR_MuMuKK_K2;
   UInt_t          nX;
   UInt_t          nX_pre0;
   UInt_t          nX_pre1;
   UInt_t          nX_pre2;
   UInt_t          nX_pre3;
   UInt_t          nX_pre4;
   UInt_t          nX_pre5;
   UInt_t          nX_pre6;
   UInt_t          nX_pre7;
   UInt_t          nX_pre8;
   UInt_t          nX_pre9;
   UInt_t          nX_pre10;
   UInt_t          nX_pre11;
   UInt_t          nX_pre12;
   UInt_t          nX_pre13;
   UInt_t          nX_pre14;
   UInt_t          nX_pre15;
   vector<float>   *XMass;
   vector<float>   *XPx;
   vector<float>   *XPy;
   vector<float>   *XPz;
   vector<double>  *XPxE;
   vector<double>  *XPyE;
   vector<double>  *XPzE;
   vector<float>   *XVtx_CL;
   vector<float>   *XVtx_Chi2;
   vector<float>   *XDecayVtx_X;
   vector<float>   *XDecayVtx_Y;
   vector<float>   *XDecayVtx_Z;
   vector<double>  *XDecayVtx_XE;
   vector<double>  *XDecayVtx_YE;
   vector<double>  *XDecayVtx_ZE;
   vector<double>  *XCosAlphaBS;
   vector<double>  *XCosAlpha3DBS;
   vector<double>  *XCTauBS;
   vector<double>  *XCTauBSE;
   vector<double>  *XLxyBS;
   vector<double>  *XLxyBSE;
   vector<double>  *XLxyzBS;
   vector<double>  *XLxyzBSE;
   vector<double>  *XCosAlphaPV;
   vector<double>  *XCosAlpha3DPV;
   vector<double>  *XCTauPV;
   vector<double>  *XCTauPVE;
   vector<double>  *XLxyPV;
   vector<double>  *XLxyPVE;
   vector<double>  *XLxyzPV;
   vector<double>  *XLxyzPVE;
   vector<int>     *PriVtx_XCosAlpha_n;
   vector<float>   *PriVtx_XCosAlpha_X;
   vector<float>   *PriVtx_XCosAlpha_Y;
   vector<float>   *PriVtx_XCosAlpha_Z;
   vector<double>  *PriVtx_XCosAlpha_EX;
   vector<double>  *PriVtx_XCosAlpha_EY;
   vector<double>  *PriVtx_XCosAlpha_EZ;
   vector<float>   *PriVtx_XCosAlpha_Chi2;
   vector<float>   *PriVtx_XCosAlpha_CL;
   vector<int>     *PriVtx_XCosAlpha_tracks;
   vector<double>  *XCosAlphaPVCosAlpha;
   vector<double>  *XCosAlpha3DPVCosAlpha;
   vector<double>  *XCTauPVCosAlpha;
   vector<double>  *XCTauPVCosAlphaE;
   vector<double>  *XLxyPVCosAlpha;
   vector<double>  *XLxyPVCosAlphaE;
   vector<double>  *XLxyzPVCosAlpha;
   vector<double>  *XLxyzPVCosAlphaE;
   vector<int>     *PriVtx_XCosAlpha3D_n;
   vector<float>   *PriVtx_XCosAlpha3D_X;
   vector<float>   *PriVtx_XCosAlpha3D_Y;
   vector<float>   *PriVtx_XCosAlpha3D_Z;
   vector<double>  *PriVtx_XCosAlpha3D_EX;
   vector<double>  *PriVtx_XCosAlpha3D_EY;
   vector<double>  *PriVtx_XCosAlpha3D_EZ;
   vector<float>   *PriVtx_XCosAlpha3D_Chi2;
   vector<float>   *PriVtx_XCosAlpha3D_CL;
   vector<int>     *PriVtx_XCosAlpha3D_tracks;
   vector<double>  *XCosAlphaPVCosAlpha3D;
   vector<double>  *XCosAlpha3DPVCosAlpha3D;
   vector<double>  *XCTauPVCosAlpha3D;
   vector<double>  *XCTauPVCosAlpha3DE;
   vector<double>  *XLxyPVCosAlpha3D;
   vector<double>  *XLxyPVCosAlpha3DE;
   vector<double>  *XLxyzPVCosAlpha3D;
   vector<double>  *XLxyzPVCosAlpha3DE;
   vector<float>   *XLessPV_tracksPtSq;
   vector<float>   *XLessPV_4tracksPtSq;
   vector<int>     *PriVtxXLess_n;
   vector<float>   *PriVtxXLess_X;
   vector<float>   *PriVtxXLess_Y;
   vector<float>   *PriVtxXLess_Z;
   vector<double>  *PriVtxXLess_EX;
   vector<double>  *PriVtxXLess_EY;
   vector<double>  *PriVtxXLess_EZ;
   vector<float>   *PriVtxXLess_Chi2;
   vector<float>   *PriVtxXLess_CL;
   vector<int>     *PriVtxXLess_tracks;
   vector<double>  *XCosAlphaXLessPV;
   vector<double>  *XCosAlpha3DXLessPV;
   vector<double>  *XCTauXLessPV;
   vector<double>  *XCTauXLessPVE;
   vector<double>  *XLxyXLessPV;
   vector<double>  *XLxyXLessPVE;
   vector<double>  *XLxyzXLessPV;
   vector<double>  *XLxyzXLessPVE;
   vector<int>     *PriVtxXLess_XCosAlpha_n;
   vector<float>   *PriVtxXLess_XCosAlpha_X;
   vector<float>   *PriVtxXLess_XCosAlpha_Y;
   vector<float>   *PriVtxXLess_XCosAlpha_Z;
   vector<double>  *PriVtxXLess_XCosAlpha_EX;
   vector<double>  *PriVtxXLess_XCosAlpha_EY;
   vector<double>  *PriVtxXLess_XCosAlpha_EZ;
   vector<float>   *PriVtxXLess_XCosAlpha_Chi2;
   vector<float>   *PriVtxXLess_XCosAlpha_CL;
   vector<int>     *PriVtxXLess_XCosAlpha_tracks;
   vector<double>  *XCosAlphaXLessPVCosAlpha;
   vector<double>  *XCosAlpha3DXLessPVCosAlpha;
   vector<double>  *XCTauXLessPVCosAlpha;
   vector<double>  *XCTauXLessPVCosAlphaE;
   vector<double>  *XLxyXLessPVCosAlpha;
   vector<double>  *XLxyXLessPVCosAlphaE;
   vector<double>  *XLxyzXLessPVCosAlpha;
   vector<double>  *XLxyzXLessPVCosAlphaE;
   vector<int>     *PriVtxXLess_XCosAlpha3D_n;
   vector<float>   *PriVtxXLess_XCosAlpha3D_X;
   vector<float>   *PriVtxXLess_XCosAlpha3D_Y;
   vector<float>   *PriVtxXLess_XCosAlpha3D_Z;
   vector<double>  *PriVtxXLess_XCosAlpha3D_EX;
   vector<double>  *PriVtxXLess_XCosAlpha3D_EY;
   vector<double>  *PriVtxXLess_XCosAlpha3D_EZ;
   vector<float>   *PriVtxXLess_XCosAlpha3D_Chi2;
   vector<float>   *PriVtxXLess_XCosAlpha3D_CL;
   vector<int>     *PriVtxXLess_XCosAlpha3D_tracks;
   vector<double>  *XCosAlphaXLessPVCosAlpha3D;
   vector<double>  *XCosAlpha3DXLessPVCosAlpha3D;
   vector<double>  *XCTauXLessPVCosAlpha3D;
   vector<double>  *XCTauXLessPVCosAlpha3DE;
   vector<double>  *XLxyXLessPVCosAlpha3D;
   vector<double>  *XLxyXLessPVCosAlpha3DE;
   vector<double>  *XLxyzXLessPVCosAlpha3D;
   vector<double>  *XLxyzXLessPVCosAlpha3DE;
   vector<int>     *PriVtxXCorr_n;
   vector<float>   *PriVtxXCorr_X;
   vector<float>   *PriVtxXCorr_Y;
   vector<float>   *PriVtxXCorr_Z;
   vector<double>  *PriVtxXCorr_EX;
   vector<double>  *PriVtxXCorr_EY;
   vector<double>  *PriVtxXCorr_EZ;
   vector<float>   *PriVtxXCorr_Chi2;
   vector<float>   *PriVtxXCorr_CL;
   vector<int>     *PriVtxXCorr_tracks;
   vector<double>  *XCosAlphaPVX;
   vector<double>  *XCTauPVX;
   vector<double>  *XCTauPVXE;
   vector<double>  *XLxyPVX;
   vector<double>  *XLxyzPVX;
   vector<float>   *XCTauPVX_3D;
   vector<float>   *XCTauPVX_3D_err;
   vector<float>   *kaon1_dxy_PV;
   vector<float>   *kaon1_dz_PV;
   vector<float>   *kaon2_dxy_PV;
   vector<float>   *kaon2_dz_PV;
   vector<float>   *kaon1_dxy_BS;
   vector<float>   *kaon1_dz_BS;
   vector<float>   *kaon2_dxy_BS;
   vector<float>   *kaon2_dz_BS;
   vector<float>   *kaon1_dxy_XLessPV;
   vector<float>   *kaon1_dz_XLessPV;
   vector<float>   *kaon2_dxy_XLessPV;
   vector<float>   *kaon2_dz_XLessPV;
   vector<float>   *kaon1_dxyE;
   vector<float>   *kaon1_dzE;
   vector<float>   *kaon2_dxyE;
   vector<float>   *kaon2_dzE;
   vector<int>     *XMuMuIdx;
   vector<int>     *XKaon1Idx;
   vector<int>     *XKaon2Idx;
   vector<bool>    *Kaon1FromPV;
   vector<bool>    *Kaon2FromPV;
   vector<float>   *Muon1Px_MuMuKK;
   vector<float>   *Muon1Py_MuMuKK;
   vector<float>   *Muon1Pz_MuMuKK;
   vector<float>   *Muon1E_MuMuKK;
   vector<float>   *Muon2Px_MuMuKK;
   vector<float>   *Muon2Py_MuMuKK;
   vector<float>   *Muon2Pz_MuMuKK;
   vector<float>   *Muon2E_MuMuKK;
   vector<float>   *Kaon1Px_MuMuKK;
   vector<float>   *Kaon1Py_MuMuKK;
   vector<float>   *Kaon1Pz_MuMuKK;
   vector<float>   *Kaon1E_MuMuKK;
   vector<double>  *kaon1_nsigdedx;
   vector<float>   *kaon1_dedx;
   vector<float>   *kaon1_dedxMass;
   vector<float>   *kaon1_theo;
   vector<float>   *kaon1_sigma;
   vector<float>   *kaon1_dedx_byHits;
   vector<float>   *kaon1_dedxErr_byHits;
   vector<int>     *kaon1_saturMeas_byHits;
   vector<int>     *kaon1_Meas_byHits;
   vector<float>   *Kaon2Px_MuMuKK;
   vector<float>   *Kaon2Py_MuMuKK;
   vector<float>   *Kaon2Pz_MuMuKK;
   vector<float>   *Kaon2E_MuMuKK;
   vector<double>  *kaon2_nsigdedx;
   vector<float>   *kaon2_dedx;
   vector<float>   *kaon2_dedxMass;
   vector<float>   *kaon2_theo;
   vector<float>   *kaon2_sigma;
   vector<float>   *kaon2_dedx_byHits;
   vector<float>   *kaon2_dedxErr_byHits;
   vector<int>     *kaon2_saturMeas_byHits;
   vector<int>     *kaon2_Meas_byHits;

   // List of branches
   TBranch        *b_TrigRes;   //!
   TBranch        *b_TrigNames;   //!
   TBranch        *b_MatchTriggerNames;   //!
   TBranch        *b_L1TrigRes;   //!
   TBranch        *b_evtNum;   //!
   TBranch        *b_runNum;   //!
   TBranch        *b_lumiNum;   //!
   TBranch        *b_priVtx_n;   //!
   TBranch        *b_priVtx_X;   //!
   TBranch        *b_priVtx_Y;   //!
   TBranch        *b_priVtx_Z;   //!
   TBranch        *b_priVtx_XE;   //!
   TBranch        *b_priVtx_YE;   //!
   TBranch        *b_priVtx_ZE;   //!
   TBranch        *b_priVtx_NormChi2;   //!
   TBranch        *b_priVtx_Chi2;   //!
   TBranch        *b_priVtx_CL;   //!
   TBranch        *b_priVtx_tracks;   //!
   TBranch        *b_priVtx_tracksPtSq;   //!
   TBranch        *b_trNotRef;   //!
   TBranch        *b_trRef;   //!
   TBranch        *b_trackPx;   //!
   TBranch        *b_trackPy;   //!
   TBranch        *b_trackPz;   //!
   TBranch        *b_trackEnergy;   //!
   TBranch        *b_trackNDF;   //!
   TBranch        *b_trackPhits;   //!
   TBranch        *b_trackShits;   //!
   TBranch        *b_trackChi2;   //!
   TBranch        *b_trackD0;   //!
   TBranch        *b_trackD0Err;   //!
   TBranch        *b_trackCharge;   //!
   TBranch        *b_TrackHighPurity;   //!
   TBranch        *b_TrackTight;   //!
   TBranch        *b_trackfHits;   //!
   TBranch        *b_trackFirstBarrel;   //!
   TBranch        *b_trackFirstEndCap;   //!
   TBranch        *b_trackDzVtx;   //!
   TBranch        *b_trackDxyVtx;   //!
   TBranch        *b_tr_nsigdedx;   //!
   TBranch        *b_tr_dedx;   //!
   TBranch        *b_tr_dedxMass;   //!
   TBranch        *b_tr_theo;   //!
   TBranch        *b_tr_sigma;   //!
   TBranch        *b_tr_dedx_byHits;   //!
   TBranch        *b_tr_dedxErr_byHits;   //!
   TBranch        *b_tr_saturMeas_byHits;   //!
   TBranch        *b_tr_Meas_byHits;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_muPx;   //!
   TBranch        *b_muPy;   //!
   TBranch        *b_muPz;   //!
   TBranch        *b_muCharge;   //!
   TBranch        *b_muD0;   //!
   TBranch        *b_muDz;   //!
   TBranch        *b_muChi2;   //!
   TBranch        *b_muNDF;   //!
   TBranch        *b_muPhits;   //!
   TBranch        *b_muShits;   //!
   TBranch        *b_muLayersTr;   //!
   TBranch        *b_muLayersPix;   //!
   TBranch        *b_muD0E;   //!
   TBranch        *b_muDzVtxErr;   //!
   TBranch        *b_muKey;   //!
   TBranch        *b_muIsGlobal;   //!
   TBranch        *b_muIsPF;   //!
   TBranch        *b_muGlMuHits;   //!
   TBranch        *b_muGlChi2;   //!
   TBranch        *b_muGlNDF;   //!
   TBranch        *b_muGlMatchedStation;   //!
   TBranch        *b_muGlDzVtx;   //!
   TBranch        *b_muGlDxyVtx;   //!
   TBranch        *b_nMatchedStations;   //!
   TBranch        *b_muType;   //!
   TBranch        *b_muQual;   //!
   TBranch        *b_muTrack;   //!
   TBranch        *b_muNOverlap;   //!
   TBranch        *b_muNSharingSegWith;   //!
   TBranch        *b_mufHits;   //!
   TBranch        *b_muFirstBarrel;   //!
   TBranch        *b_muFirstEndCap;   //!
   TBranch        *b_muDzVtx;   //!
   TBranch        *b_muDxyVtx;   //!
   TBranch        *b_nMuMu;   //!
   TBranch        *b_MuMuMass;   //!
   TBranch        *b_MuMuPx;   //!
   TBranch        *b_MuMuPy;   //!
   TBranch        *b_MuMuPz;   //!
   TBranch        *b_MuMuVtx_CL;   //!
   TBranch        *b_MuMuVtx_Chi2;   //!
   TBranch        *b_MuMuDecayVtx_X;   //!
   TBranch        *b_MuMuDecayVtx_Y;   //!
   TBranch        *b_MuMuDecayVtx_Z;   //!
   TBranch        *b_MuMuDecayVtx_XE;   //!
   TBranch        *b_MuMuDecayVtx_YE;   //!
   TBranch        *b_MuMuDecayVtx_ZE;   //!
   TBranch        *b_mu1Idx;   //!
   TBranch        *b_mu2Idx;   //!
   TBranch        *b_mu1Px_MuMu;   //!
   TBranch        *b_mu1Py_MuMu;   //!
   TBranch        *b_mu1Pz_MuMu;   //!
   TBranch        *b_mu1Chi2_MuMu;   //!
   TBranch        *b_mu1NDF_MuMu;   //!
   TBranch        *b_mu2Px_MuMu;   //!
   TBranch        *b_mu2Py_MuMu;   //!
   TBranch        *b_mu2Pz_MuMu;   //!
   TBranch        *b_mu2Chi2_MuMu;   //!
   TBranch        *b_mu2NDF_MuMu;   //!
   TBranch        *b_MuMuType;   //!
   TBranch        *b_MuMuMuonTrigMatch;   //!
   TBranch        *b_PriVtxMuMuCorr_n;   //!
   TBranch        *b_PriVtxMuMuCorr_X;   //!
   TBranch        *b_PriVtxMuMuCorr_Y;   //!
   TBranch        *b_PriVtxMuMuCorr_Z;   //!
   TBranch        *b_PriVtxMuMuCorr_EX;   //!
   TBranch        *b_PriVtxMuMuCorr_EY;   //!
   TBranch        *b_PriVtxMuMuCorr_EZ;   //!
   TBranch        *b_PriVtxMuMuCorr_Chi2;   //!
   TBranch        *b_PriVtxMuMuCorr_CL;   //!
   TBranch        *b_PriVtxMuMuCorr_tracks;   //!
   TBranch        *b_nTrk_afterMuMu;   //!
   TBranch        *b_nKK;   //!
   TBranch        *b_KKMass;   //!
   TBranch        *b_KKPx;   //!
   TBranch        *b_KKPy;   //!
   TBranch        *b_KKPz;   //!
   TBranch        *b_KKVtx_CL;   //!
   TBranch        *b_KKVtx_Chi2;   //!
   TBranch        *b_KKDecayVtx_X;   //!
   TBranch        *b_KKDecayVtx_Y;   //!
   TBranch        *b_KKDecayVtx_Z;   //!
   TBranch        *b_KKDecayVtx_XE;   //!
   TBranch        *b_KKDecayVtx_YE;   //!
   TBranch        *b_KKDecayVtx_ZE;   //!
   TBranch        *b_ka1Idx;   //!
   TBranch        *b_ka2Idx;   //!
   TBranch        *b_ka1Px_KK;   //!
   TBranch        *b_ka1Py_KK;   //!
   TBranch        *b_ka1Pz_KK;   //!
   TBranch        *b_ka1Chi2_KK;   //!
   TBranch        *b_ka1NDF_KK;   //!
   TBranch        *b_ka2Px_KK;   //!
   TBranch        *b_ka2Py_KK;   //!
   TBranch        *b_ka2Pz_KK;   //!
   TBranch        *b_ka2Chi2_KK;   //!
   TBranch        *b_ka2NDF_KK;   //!
   TBranch        *b_DR_MuMu_K1;   //!
   TBranch        *b_DR_MuMu_K2;   //!
   TBranch        *b_DR_MuMuKK_K1;   //!
   TBranch        *b_DR_MuMuKK_K2;   //!
   TBranch        *b_nX;   //!
   TBranch        *b_nX_pre0;   //!
   TBranch        *b_nX_pre1;   //!
   TBranch        *b_nX_pre2;   //!
   TBranch        *b_nX_pre3;   //!
   TBranch        *b_nX_pre4;   //!
   TBranch        *b_nX_pre5;   //!
   TBranch        *b_nX_pre6;   //!
   TBranch        *b_nX_pre7;   //!
   TBranch        *b_nX_pre8;   //!
   TBranch        *b_nX_pre9;   //!
   TBranch        *b_nX_pre10;   //!
   TBranch        *b_nX_pre11;   //!
   TBranch        *b_nX_pre12;   //!
   TBranch        *b_nX_pre13;   //!
   TBranch        *b_nX_pre14;   //!
   TBranch        *b_nX_pre15;   //!
   TBranch        *b_XMass;   //!
   TBranch        *b_XPx;   //!
   TBranch        *b_XPy;   //!
   TBranch        *b_XPz;   //!
   TBranch        *b_XPxE;   //!
   TBranch        *b_XPyE;   //!
   TBranch        *b_XPzE;   //!
   TBranch        *b_XVtx_CL;   //!
   TBranch        *b_XVtx_Chi2;   //!
   TBranch        *b_XDecayVtx_X;   //!
   TBranch        *b_XDecayVtx_Y;   //!
   TBranch        *b_XDecayVtx_Z;   //!
   TBranch        *b_XDecayVtx_XE;   //!
   TBranch        *b_XDecayVtx_YE;   //!
   TBranch        *b_XDecayVtx_ZE;   //!
   TBranch        *b_XCosAlphaBS;   //!
   TBranch        *b_XCosAlpha3DBS;   //!
   TBranch        *b_XCTauBS;   //!
   TBranch        *b_XCTauBSE;   //!
   TBranch        *b_XLxyBS;   //!
   TBranch        *b_XLxyBSE;   //!
   TBranch        *b_XLxyzBS;   //!
   TBranch        *b_XLxyzBSE;   //!
   TBranch        *b_XCosAlphaPV;   //!
   TBranch        *b_XCosAlpha3DPV;   //!
   TBranch        *b_XCTauPV;   //!
   TBranch        *b_XCTauPVE;   //!
   TBranch        *b_XLxyPV;   //!
   TBranch        *b_XLxyPVE;   //!
   TBranch        *b_XLxyzPV;   //!
   TBranch        *b_XLxyzPVE;   //!
   TBranch        *b_PriVtx_XCosAlpha_n;   //!
   TBranch        *b_PriVtx_XCosAlpha_X;   //!
   TBranch        *b_PriVtx_XCosAlpha_Y;   //!
   TBranch        *b_PriVtx_XCosAlpha_Z;   //!
   TBranch        *b_PriVtx_XCosAlpha_EX;   //!
   TBranch        *b_PriVtx_XCosAlpha_EY;   //!
   TBranch        *b_PriVtx_XCosAlpha_EZ;   //!
   TBranch        *b_PriVtx_XCosAlpha_Chi2;   //!
   TBranch        *b_PriVtx_XCosAlpha_CL;   //!
   TBranch        *b_PriVtx_XCosAlpha_tracks;   //!
   TBranch        *b_XCosAlphaPVCosAlpha;   //!
   TBranch        *b_XCosAlpha3DPVCosAlpha;   //!
   TBranch        *b_XCTauPVCosAlpha;   //!
   TBranch        *b_XCTauPVCosAlphaE;   //!
   TBranch        *b_XLxyPVCosAlpha;   //!
   TBranch        *b_XLxyPVCosAlphaE;   //!
   TBranch        *b_XLxyzPVCosAlpha;   //!
   TBranch        *b_XLxyzPVCosAlphaE;   //!
   TBranch        *b_PriVtx_XCosAlpha3D_n;   //!
   TBranch        *b_PriVtx_XCosAlpha3D_X;   //!
   TBranch        *b_PriVtx_XCosAlpha3D_Y;   //!
   TBranch        *b_PriVtx_XCosAlpha3D_Z;   //!
   TBranch        *b_PriVtx_XCosAlpha3D_EX;   //!
   TBranch        *b_PriVtx_XCosAlpha3D_EY;   //!
   TBranch        *b_PriVtx_XCosAlpha3D_EZ;   //!
   TBranch        *b_PriVtx_XCosAlpha3D_Chi2;   //!
   TBranch        *b_PriVtx_XCosAlpha3D_CL;   //!
   TBranch        *b_PriVtx_XCosAlpha3D_tracks;   //!
   TBranch        *b_XCosAlphaPVCosAlpha3D;   //!
   TBranch        *b_XCosAlpha3DPVCosAlpha3D;   //!
   TBranch        *b_XCTauPVCosAlpha3D;   //!
   TBranch        *b_XCTauPVCosAlpha3DE;   //!
   TBranch        *b_XLxyPVCosAlpha3D;   //!
   TBranch        *b_XLxyPVCosAlpha3DE;   //!
   TBranch        *b_XLxyzPVCosAlpha3D;   //!
   TBranch        *b_XLxyzPVCosAlpha3DE;   //!
   TBranch        *b_XLessPV_tracksPtSq;   //!
   TBranch        *b_XLessPV_4tracksPtSq;   //!
   TBranch        *b_PriVtxXLess_n;   //!
   TBranch        *b_PriVtxXLess_X;   //!
   TBranch        *b_PriVtxXLess_Y;   //!
   TBranch        *b_PriVtxXLess_Z;   //!
   TBranch        *b_PriVtxXLess_EX;   //!
   TBranch        *b_PriVtxXLess_EY;   //!
   TBranch        *b_PriVtxXLess_EZ;   //!
   TBranch        *b_PriVtxXLess_Chi2;   //!
   TBranch        *b_PriVtxXLess_CL;   //!
   TBranch        *b_PriVtxXLess_tracks;   //!
   TBranch        *b_XCosAlphaXLessPV;   //!
   TBranch        *b_XCosAlpha3DXLessPV;   //!
   TBranch        *b_XCTauXLessPV;   //!
   TBranch        *b_XCTauXLessPVE;   //!
   TBranch        *b_XLxyXLessPV;   //!
   TBranch        *b_XLxyXLessPVE;   //!
   TBranch        *b_XLxyzXLessPV;   //!
   TBranch        *b_XLxyzXLessPVE;   //!
   TBranch        *b_PriVtxXLess_XCosAlpha_n;   //!
   TBranch        *b_PriVtxXLess_XCosAlpha_X;   //!
   TBranch        *b_PriVtxXLess_XCosAlpha_Y;   //!
   TBranch        *b_PriVtxXLess_XCosAlpha_Z;   //!
   TBranch        *b_PriVtxXLess_XCosAlpha_EX;   //!
   TBranch        *b_PriVtxXLess_XCosAlpha_EY;   //!
   TBranch        *b_PriVtxXLess_XCosAlpha_EZ;   //!
   TBranch        *b_PriVtxXLess_XCosAlpha_Chi2;   //!
   TBranch        *b_PriVtxXLess_XCosAlpha_CL;   //!
   TBranch        *b_PriVtxXLess_XCosAlpha_tracks;   //!
   TBranch        *b_XCosAlphaXLessPVCosAlpha;   //!
   TBranch        *b_XCosAlpha3DXLessPVCosAlpha;   //!
   TBranch        *b_XCTauXLessPVCosAlpha;   //!
   TBranch        *b_XCTauXLessPVCosAlphaE;   //!
   TBranch        *b_XLxyXLessPVCosAlpha;   //!
   TBranch        *b_XLxyXLessPVCosAlphaE;   //!
   TBranch        *b_XLxyzXLessPVCosAlpha;   //!
   TBranch        *b_XLxyzXLessPVCosAlphaE;   //!
   TBranch        *b_PriVtxXLess_XCosAlpha3D_n;   //!
   TBranch        *b_PriVtxXLess_XCosAlpha3D_X;   //!
   TBranch        *b_PriVtxXLess_XCosAlpha3D_Y;   //!
   TBranch        *b_PriVtxXLess_XCosAlpha3D_Z;   //!
   TBranch        *b_PriVtxXLess_XCosAlpha3D_EX;   //!
   TBranch        *b_PriVtxXLess_XCosAlpha3D_EY;   //!
   TBranch        *b_PriVtxXLess_XCosAlpha3D_EZ;   //!
   TBranch        *b_PriVtxXLess_XCosAlpha3D_Chi2;   //!
   TBranch        *b_PriVtxXLess_XCosAlpha3D_CL;   //!
   TBranch        *b_PriVtxXLess_XCosAlpha3D_tracks;   //!
   TBranch        *b_XCosAlphaXLessPVCosAlpha3D;   //!
   TBranch        *b_XCosAlpha3DXLessPVCosAlpha3D;   //!
   TBranch        *b_XCTauXLessPVCosAlpha3D;   //!
   TBranch        *b_XCTauXLessPVCosAlpha3DE;   //!
   TBranch        *b_XLxyXLessPVCosAlpha3D;   //!
   TBranch        *b_XLxyXLessPVCosAlpha3DE;   //!
   TBranch        *b_XLxyzXLessPVCosAlpha3D;   //!
   TBranch        *b_XLxyzXLessPVCosAlpha3DE;   //!
   TBranch        *b_PriVtxXCorr_n;   //!
   TBranch        *b_PriVtxXCorr_X;   //!
   TBranch        *b_PriVtxXCorr_Y;   //!
   TBranch        *b_PriVtxXCorr_Z;   //!
   TBranch        *b_PriVtxXCorr_EX;   //!
   TBranch        *b_PriVtxXCorr_EY;   //!
   TBranch        *b_PriVtxXCorr_EZ;   //!
   TBranch        *b_PriVtxXCorr_Chi2;   //!
   TBranch        *b_PriVtxXCorr_CL;   //!
   TBranch        *b_PriVtxXCorr_tracks;   //!
   TBranch        *b_XCosAlphaPVX;   //!
   TBranch        *b_XCTauPVX;   //!
   TBranch        *b_XCTauPVXE;   //!
   TBranch        *b_XLxyPVX;   //!
   TBranch        *b_XLxyzPVX;   //!
   TBranch        *b_XCTauPVX_3D;   //!
   TBranch        *b_XCTauPVX_3D_err;   //!
   TBranch        *b_kaon1_dxy_PV;   //!
   TBranch        *b_kaon1_dz_PV;   //!
   TBranch        *b_kaon2_dxy_PV;   //!
   TBranch        *b_kaon2_dz_PV;   //!
   TBranch        *b_kaon1_dxy_BS;   //!
   TBranch        *b_kaon1_dz_BS;   //!
   TBranch        *b_kaon2_dxy_BS;   //!
   TBranch        *b_kaon2_dz_BS;   //!
   TBranch        *b_kaon1_dxy_XLessPV;   //!
   TBranch        *b_kaon1_dz_XLessPV;   //!
   TBranch        *b_kaon2_dxy_XLessPV;   //!
   TBranch        *b_kaon2_dz_XLessPV;   //!
   TBranch        *b_kaon1_dxyE;   //!
   TBranch        *b_kaon1_dzE;   //!
   TBranch        *b_kaon2_dxyE;   //!
   TBranch        *b_kaon2_dzE;   //!
   TBranch        *b_XMuMuIdx;   //!
   TBranch        *b_XKaon1Idx;   //!
   TBranch        *b_XKaon2Idx;   //!
   TBranch        *b_Kaon1FromPV;   //!
   TBranch        *b_Kaon2FromPV;   //!
   TBranch        *b_Muon1Px_MuMuKK;   //!
   TBranch        *b_Muon1Py_MuMuKK;   //!
   TBranch        *b_Muon1Pz_MuMuKK;   //!
   TBranch        *b_Muon1E_MuMuKK;   //!
   TBranch        *b_Muon2Px_MuMuKK;   //!
   TBranch        *b_Muon2Py_MuMuKK;   //!
   TBranch        *b_Muon2Pz_MuMuKK;   //!
   TBranch        *b_Muon2E_MuMuKK;   //!
   TBranch        *b_Kaon1Px_MuMuKK;   //!
   TBranch        *b_Kaon1Py_MuMuKK;   //!
   TBranch        *b_Kaon1Pz_MuMuKK;   //!
   TBranch        *b_Kaon1E_MuMuKK;   //!
   TBranch        *b_kaon1_nsigdedx;   //!
   TBranch        *b_kaon1_dedx;   //!
   TBranch        *b_kaon1_dedxMass;   //!
   TBranch        *b_kaon1_theo;   //!
   TBranch        *b_kaon1_sigma;   //!
   TBranch        *b_kaon1_dedx_byHits;   //!
   TBranch        *b_kaon1_dedxErr_byHits;   //!
   TBranch        *b_kaon1_saturMeas_byHits;   //!
   TBranch        *b_kaon1_Meas_byHits;   //!
   TBranch        *b_Kaon2Px_MuMuKK;   //!
   TBranch        *b_Kaon2Py_MuMuKK;   //!
   TBranch        *b_Kaon2Pz_MuMuKK;   //!
   TBranch        *b_Kaon2E_MuMuKK;   //!
   TBranch        *b_kaon2_nsigdedx;   //!
   TBranch        *b_kaon2_dedx;   //!
   TBranch        *b_kaon2_dedxMass;   //!
   TBranch        *b_kaon2_theo;   //!
   TBranch        *b_kaon2_sigma;   //!
   TBranch        *b_kaon2_dedx_byHits;   //!
   TBranch        *b_kaon2_dedxErr_byHits;   //!
   TBranch        *b_kaon2_saturMeas_byHits;   //!
   TBranch        *b_kaon2_Meas_byHits;   //!
   ////////////////////////////////////////////////////////////
   mumukk(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~mumukk() { }
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

   // TTree *outTree;
   TNtuple *outTuple;
   //
   Float_t          run_out;
   Float_t          evt_out;
   Float_t          lum_out;
   Float_t         X_mass;
   Float_t         kk_mass;
   Float_t         mumu_mass;
   Float_t         X_LFly;
   Float_t         X_pt;
   Float_t         X_eta;
   Float_t         X_vtx;
   Float_t         X_cosAlpha;
   Float_t         X_hlt;
  //
  // TBranch*      X_mass_b;
  // TBranch*      kk_mass_b;
  // TBranch*      mumu_mass_b;
  // TBranch*      X_LFly_b;
  // TBranch*      X_pt_b;
  // TBranch*      X_eta_b;
  // TBranch*      X_vtx_b;
  // TBranch*      X_cosAlpha_b;
  // TBranch*      X_hlt_b;
  // TBranch*      X_run_b;
  // TBranch*      X_evt_b;
  // TBranch*      X_lum_b;

   double JPsi_mass;
   double Phi_mass,Phi_sigma,Phi_mean;
   double Bs0_Low_Mass;
   double Bs0_High_Mass;
   double Y_High_Mass;

  //  int muonQual[4] = {1,3,4,12};


   ClassDef(mumukk,0);
};

#endif

#ifdef mumukk_cxx
void mumukk::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   TrigRes = 0;
   TrigNames = 0;
   MatchTriggerNames = 0;
   L1TrigRes = 0;
   trNotRef = 0;
   trRef = 0;
   trackPx = 0;
   trackPy = 0;
   trackPz = 0;
   trackEnergy = 0;
   trackNDF = 0;
   trackPhits = 0;
   trackShits = 0;
   trackChi2 = 0;
   trackD0 = 0;
   trackD0Err = 0;
   trackCharge = 0;
   TrackHighPurity = 0;
   TrackTight = 0;
   trackfHits = 0;
   trackFirstBarrel = 0;
   trackFirstEndCap = 0;
   trackDzVtx = 0;
   trackDxyVtx = 0;
   tr_nsigdedx = 0;
   tr_dedx = 0;
   tr_dedxMass = 0;
   tr_theo = 0;
   tr_sigma = 0;
   tr_dedx_byHits = 0;
   tr_dedxErr_byHits = 0;
   tr_saturMeas_byHits = 0;
   tr_Meas_byHits = 0;
   muPx = 0;
   muPy = 0;
   muPz = 0;
   muCharge = 0;
   muD0 = 0;
   muDz = 0;
   muChi2 = 0;
   muNDF = 0;
   muPhits = 0;
   muShits = 0;
   muLayersTr = 0;
   muLayersPix = 0;
   muD0E = 0;
   muDzVtxErr = 0;
   muKey = 0;
   muIsGlobal = 0;
   muIsPF = 0;
   muGlMuHits = 0;
   muGlChi2 = 0;
   muGlNDF = 0;
   muGlMatchedStation = 0;
   muGlDzVtx = 0;
   muGlDxyVtx = 0;
   nMatchedStations = 0;
   muType = 0;
   muQual = 0;
   muTrack = 0;
   muNOverlap = 0;
   muNSharingSegWith = 0;
   mufHits = 0;
   muFirstBarrel = 0;
   muFirstEndCap = 0;
   muDzVtx = 0;
   muDxyVtx = 0;
   MuMuMass = 0;
   MuMuPx = 0;
   MuMuPy = 0;
   MuMuPz = 0;
   MuMuVtx_CL = 0;
   MuMuVtx_Chi2 = 0;
   MuMuDecayVtx_X = 0;
   MuMuDecayVtx_Y = 0;
   MuMuDecayVtx_Z = 0;
   MuMuDecayVtx_XE = 0;
   MuMuDecayVtx_YE = 0;
   MuMuDecayVtx_ZE = 0;
   mu1Idx = 0;
   mu2Idx = 0;
   mu1Px_MuMu = 0;
   mu1Py_MuMu = 0;
   mu1Pz_MuMu = 0;
   mu1Chi2_MuMu = 0;
   mu1NDF_MuMu = 0;
   mu2Px_MuMu = 0;
   mu2Py_MuMu = 0;
   mu2Pz_MuMu = 0;
   mu2Chi2_MuMu = 0;
   mu2NDF_MuMu = 0;
   MuMuType = 0;
   MuMuMuonTrigMatch = 0;
   PriVtxMuMuCorr_n = 0;
   PriVtxMuMuCorr_X = 0;
   PriVtxMuMuCorr_Y = 0;
   PriVtxMuMuCorr_Z = 0;
   PriVtxMuMuCorr_EX = 0;
   PriVtxMuMuCorr_EY = 0;
   PriVtxMuMuCorr_EZ = 0;
   PriVtxMuMuCorr_Chi2 = 0;
   PriVtxMuMuCorr_CL = 0;
   PriVtxMuMuCorr_tracks = 0;
   nTrk_afterMuMu = 0;
   KKMass = 0;
   KKPx = 0;
   KKPy = 0;
   KKPz = 0;
   KKVtx_CL = 0;
   KKVtx_Chi2 = 0;
   KKDecayVtx_X = 0;
   KKDecayVtx_Y = 0;
   KKDecayVtx_Z = 0;
   KKDecayVtx_XE = 0;
   KKDecayVtx_YE = 0;
   KKDecayVtx_ZE = 0;
   ka1Idx = 0;
   ka2Idx = 0;
   ka1Px_KK = 0;
   ka1Py_KK = 0;
   ka1Pz_KK = 0;
   ka1Chi2_KK = 0;
   ka1NDF_KK = 0;
   ka2Px_KK = 0;
   ka2Py_KK = 0;
   ka2Pz_KK = 0;
   ka2Chi2_KK = 0;
   ka2NDF_KK = 0;
   DR_MuMu_K1 = 0;
   DR_MuMu_K2 = 0;
   DR_MuMuKK_K1 = 0;
   DR_MuMuKK_K2 = 0;
   XMass = 0;
   XPx = 0;
   XPy = 0;
   XPz = 0;
   XPxE = 0;
   XPyE = 0;
   XPzE = 0;
   XVtx_CL = 0;
   XVtx_Chi2 = 0;
   XDecayVtx_X = 0;
   XDecayVtx_Y = 0;
   XDecayVtx_Z = 0;
   XDecayVtx_XE = 0;
   XDecayVtx_YE = 0;
   XDecayVtx_ZE = 0;
   XCosAlphaBS = 0;
   XCosAlpha3DBS = 0;
   XCTauBS = 0;
   XCTauBSE = 0;
   XLxyBS = 0;
   XLxyBSE = 0;
   XLxyzBS = 0;
   XLxyzBSE = 0;
   XCosAlphaPV = 0;
   XCosAlpha3DPV = 0;
   XCTauPV = 0;
   XCTauPVE = 0;
   XLxyPV = 0;
   XLxyPVE = 0;
   XLxyzPV = 0;
   XLxyzPVE = 0;
   PriVtx_XCosAlpha_n = 0;
   PriVtx_XCosAlpha_X = 0;
   PriVtx_XCosAlpha_Y = 0;
   PriVtx_XCosAlpha_Z = 0;
   PriVtx_XCosAlpha_EX = 0;
   PriVtx_XCosAlpha_EY = 0;
   PriVtx_XCosAlpha_EZ = 0;
   PriVtx_XCosAlpha_Chi2 = 0;
   PriVtx_XCosAlpha_CL = 0;
   PriVtx_XCosAlpha_tracks = 0;
   XCosAlphaPVCosAlpha = 0;
   XCosAlpha3DPVCosAlpha = 0;
   XCTauPVCosAlpha = 0;
   XCTauPVCosAlphaE = 0;
   XLxyPVCosAlpha = 0;
   XLxyPVCosAlphaE = 0;
   XLxyzPVCosAlpha = 0;
   XLxyzPVCosAlphaE = 0;
   PriVtx_XCosAlpha3D_n = 0;
   PriVtx_XCosAlpha3D_X = 0;
   PriVtx_XCosAlpha3D_Y = 0;
   PriVtx_XCosAlpha3D_Z = 0;
   PriVtx_XCosAlpha3D_EX = 0;
   PriVtx_XCosAlpha3D_EY = 0;
   PriVtx_XCosAlpha3D_EZ = 0;
   PriVtx_XCosAlpha3D_Chi2 = 0;
   PriVtx_XCosAlpha3D_CL = 0;
   PriVtx_XCosAlpha3D_tracks = 0;
   XCosAlphaPVCosAlpha3D = 0;
   XCosAlpha3DPVCosAlpha3D = 0;
   XCTauPVCosAlpha3D = 0;
   XCTauPVCosAlpha3DE = 0;
   XLxyPVCosAlpha3D = 0;
   XLxyPVCosAlpha3DE = 0;
   XLxyzPVCosAlpha3D = 0;
   XLxyzPVCosAlpha3DE = 0;
   XLessPV_tracksPtSq = 0;
   XLessPV_4tracksPtSq = 0;
   PriVtxXLess_n = 0;
   PriVtxXLess_X = 0;
   PriVtxXLess_Y = 0;
   PriVtxXLess_Z = 0;
   PriVtxXLess_EX = 0;
   PriVtxXLess_EY = 0;
   PriVtxXLess_EZ = 0;
   PriVtxXLess_Chi2 = 0;
   PriVtxXLess_CL = 0;
   PriVtxXLess_tracks = 0;
   XCosAlphaXLessPV = 0;
   XCosAlpha3DXLessPV = 0;
   XCTauXLessPV = 0;
   XCTauXLessPVE = 0;
   XLxyXLessPV = 0;
   XLxyXLessPVE = 0;
   XLxyzXLessPV = 0;
   XLxyzXLessPVE = 0;
   PriVtxXLess_XCosAlpha_n = 0;
   PriVtxXLess_XCosAlpha_X = 0;
   PriVtxXLess_XCosAlpha_Y = 0;
   PriVtxXLess_XCosAlpha_Z = 0;
   PriVtxXLess_XCosAlpha_EX = 0;
   PriVtxXLess_XCosAlpha_EY = 0;
   PriVtxXLess_XCosAlpha_EZ = 0;
   PriVtxXLess_XCosAlpha_Chi2 = 0;
   PriVtxXLess_XCosAlpha_CL = 0;
   PriVtxXLess_XCosAlpha_tracks = 0;
   XCosAlphaXLessPVCosAlpha = 0;
   XCosAlpha3DXLessPVCosAlpha = 0;
   XCTauXLessPVCosAlpha = 0;
   XCTauXLessPVCosAlphaE = 0;
   XLxyXLessPVCosAlpha = 0;
   XLxyXLessPVCosAlphaE = 0;
   XLxyzXLessPVCosAlpha = 0;
   XLxyzXLessPVCosAlphaE = 0;
   PriVtxXLess_XCosAlpha3D_n = 0;
   PriVtxXLess_XCosAlpha3D_X = 0;
   PriVtxXLess_XCosAlpha3D_Y = 0;
   PriVtxXLess_XCosAlpha3D_Z = 0;
   PriVtxXLess_XCosAlpha3D_EX = 0;
   PriVtxXLess_XCosAlpha3D_EY = 0;
   PriVtxXLess_XCosAlpha3D_EZ = 0;
   PriVtxXLess_XCosAlpha3D_Chi2 = 0;
   PriVtxXLess_XCosAlpha3D_CL = 0;
   PriVtxXLess_XCosAlpha3D_tracks = 0;
   XCosAlphaXLessPVCosAlpha3D = 0;
   XCosAlpha3DXLessPVCosAlpha3D = 0;
   XCTauXLessPVCosAlpha3D = 0;
   XCTauXLessPVCosAlpha3DE = 0;
   XLxyXLessPVCosAlpha3D = 0;
   XLxyXLessPVCosAlpha3DE = 0;
   XLxyzXLessPVCosAlpha3D = 0;
   XLxyzXLessPVCosAlpha3DE = 0;
   PriVtxXCorr_n = 0;
   PriVtxXCorr_X = 0;
   PriVtxXCorr_Y = 0;
   PriVtxXCorr_Z = 0;
   PriVtxXCorr_EX = 0;
   PriVtxXCorr_EY = 0;
   PriVtxXCorr_EZ = 0;
   PriVtxXCorr_Chi2 = 0;
   PriVtxXCorr_CL = 0;
   PriVtxXCorr_tracks = 0;
   XCosAlphaPVX = 0;
   XCTauPVX = 0;
   XCTauPVXE = 0;
   XLxyPVX = 0;
   XLxyzPVX = 0;
   XCTauPVX_3D = 0;
   XCTauPVX_3D_err = 0;
   kaon1_dxy_PV = 0;
   kaon1_dz_PV = 0;
   kaon2_dxy_PV = 0;
   kaon2_dz_PV = 0;
   kaon1_dxy_BS = 0;
   kaon1_dz_BS = 0;
   kaon2_dxy_BS = 0;
   kaon2_dz_BS = 0;
   kaon1_dxy_XLessPV = 0;
   kaon1_dz_XLessPV = 0;
   kaon2_dxy_XLessPV = 0;
   kaon2_dz_XLessPV = 0;
   kaon1_dxyE = 0;
   kaon1_dzE = 0;
   kaon2_dxyE = 0;
   kaon2_dzE = 0;
   XMuMuIdx = 0;
   XKaon1Idx = 0;
   XKaon2Idx = 0;
   Kaon1FromPV = 0;
   Kaon2FromPV = 0;
   Muon1Px_MuMuKK = 0;
   Muon1Py_MuMuKK = 0;
   Muon1Pz_MuMuKK = 0;
   Muon1E_MuMuKK = 0;
   Muon2Px_MuMuKK = 0;
   Muon2Py_MuMuKK = 0;
   Muon2Pz_MuMuKK = 0;
   Muon2E_MuMuKK = 0;
   Kaon1Px_MuMuKK = 0;
   Kaon1Py_MuMuKK = 0;
   Kaon1Pz_MuMuKK = 0;
   Kaon1E_MuMuKK = 0;
   kaon1_nsigdedx = 0;
   kaon1_dedx = 0;
   kaon1_dedxMass = 0;
   kaon1_theo = 0;
   kaon1_sigma = 0;
   kaon1_dedx_byHits = 0;
   kaon1_dedxErr_byHits = 0;
   kaon1_saturMeas_byHits = 0;
   kaon1_Meas_byHits = 0;
   Kaon2Px_MuMuKK = 0;
   Kaon2Py_MuMuKK = 0;
   Kaon2Pz_MuMuKK = 0;
   Kaon2E_MuMuKK = 0;
   kaon2_nsigdedx = 0;
   kaon2_dedx = 0;
   kaon2_dedxMass = 0;
   kaon2_theo = 0;
   kaon2_sigma = 0;
   kaon2_dedx_byHits = 0;
   kaon2_dedxErr_byHits = 0;
   kaon2_saturMeas_byHits = 0;
   kaon2_Meas_byHits = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("TrigRes", &TrigRes, &b_TrigRes);
   fChain->SetBranchAddress("TrigNames", &TrigNames, &b_TrigNames);
   fChain->SetBranchAddress("MatchTriggerNames", &MatchTriggerNames, &b_MatchTriggerNames);
   fChain->SetBranchAddress("L1TrigRes", &L1TrigRes, &b_L1TrigRes);
   fChain->SetBranchAddress("evtNum", &evtNum, &b_evtNum);
   fChain->SetBranchAddress("runNum", &runNum, &b_runNum);
   fChain->SetBranchAddress("lumiNum", &lumiNum, &b_lumiNum);
   fChain->SetBranchAddress("priVtx_n", &priVtx_n, &b_priVtx_n);
   fChain->SetBranchAddress("priVtx_X", &priVtx_X, &b_priVtx_X);
   fChain->SetBranchAddress("priVtx_Y", &priVtx_Y, &b_priVtx_Y);
   fChain->SetBranchAddress("priVtx_Z", &priVtx_Z, &b_priVtx_Z);
   fChain->SetBranchAddress("priVtx_XE", &priVtx_XE, &b_priVtx_XE);
   fChain->SetBranchAddress("priVtx_YE", &priVtx_YE, &b_priVtx_YE);
   fChain->SetBranchAddress("priVtx_ZE", &priVtx_ZE, &b_priVtx_ZE);
   fChain->SetBranchAddress("priVtx_NormChi2", &priVtx_NormChi2, &b_priVtx_NormChi2);
   fChain->SetBranchAddress("priVtx_Chi2", &priVtx_Chi2, &b_priVtx_Chi2);
   fChain->SetBranchAddress("priVtx_CL", &priVtx_CL, &b_priVtx_CL);
   fChain->SetBranchAddress("priVtx_tracks", &priVtx_tracks, &b_priVtx_tracks);
   fChain->SetBranchAddress("priVtx_tracksPtSq", &priVtx_tracksPtSq, &b_priVtx_tracksPtSq);
   fChain->SetBranchAddress("trNotRef", &trNotRef, &b_trNotRef);
   fChain->SetBranchAddress("trRef", &trRef, &b_trRef);
   fChain->SetBranchAddress("trackPx", &trackPx, &b_trackPx);
   fChain->SetBranchAddress("trackPy", &trackPy, &b_trackPy);
   fChain->SetBranchAddress("trackPz", &trackPz, &b_trackPz);
   fChain->SetBranchAddress("trackEnergy", &trackEnergy, &b_trackEnergy);
   fChain->SetBranchAddress("trackNDF", &trackNDF, &b_trackNDF);
   fChain->SetBranchAddress("trackPhits", &trackPhits, &b_trackPhits);
   fChain->SetBranchAddress("trackShits", &trackShits, &b_trackShits);
   fChain->SetBranchAddress("trackChi2", &trackChi2, &b_trackChi2);
   fChain->SetBranchAddress("trackD0", &trackD0, &b_trackD0);
   fChain->SetBranchAddress("trackD0Err", &trackD0Err, &b_trackD0Err);
   fChain->SetBranchAddress("trackCharge", &trackCharge, &b_trackCharge);
   fChain->SetBranchAddress("TrackHighPurity", &TrackHighPurity, &b_TrackHighPurity);
   fChain->SetBranchAddress("TrackTight", &TrackTight, &b_TrackTight);
   fChain->SetBranchAddress("trackfHits", &trackfHits, &b_trackfHits);
   fChain->SetBranchAddress("trackFirstBarrel", &trackFirstBarrel, &b_trackFirstBarrel);
   fChain->SetBranchAddress("trackFirstEndCap", &trackFirstEndCap, &b_trackFirstEndCap);
   fChain->SetBranchAddress("trackDzVtx", &trackDzVtx, &b_trackDzVtx);
   fChain->SetBranchAddress("trackDxyVtx", &trackDxyVtx, &b_trackDxyVtx);
   fChain->SetBranchAddress("tr_nsigdedx", &tr_nsigdedx, &b_tr_nsigdedx);
   fChain->SetBranchAddress("tr_dedx", &tr_dedx, &b_tr_dedx);
   fChain->SetBranchAddress("tr_dedxMass", &tr_dedxMass, &b_tr_dedxMass);
   fChain->SetBranchAddress("tr_theo", &tr_theo, &b_tr_theo);
   fChain->SetBranchAddress("tr_sigma", &tr_sigma, &b_tr_sigma);
   fChain->SetBranchAddress("tr_dedx_byHits", &tr_dedx_byHits, &b_tr_dedx_byHits);
   fChain->SetBranchAddress("tr_dedxErr_byHits", &tr_dedxErr_byHits, &b_tr_dedxErr_byHits);
   fChain->SetBranchAddress("tr_saturMeas_byHits", &tr_saturMeas_byHits, &b_tr_saturMeas_byHits);
   fChain->SetBranchAddress("tr_Meas_byHits", &tr_Meas_byHits, &b_tr_Meas_byHits);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("muPx", &muPx, &b_muPx);
   fChain->SetBranchAddress("muPy", &muPy, &b_muPy);
   fChain->SetBranchAddress("muPz", &muPz, &b_muPz);
   fChain->SetBranchAddress("muCharge", &muCharge, &b_muCharge);
   fChain->SetBranchAddress("muD0", &muD0, &b_muD0);
   fChain->SetBranchAddress("muDz", &muDz, &b_muDz);
   fChain->SetBranchAddress("muChi2", &muChi2, &b_muChi2);
   fChain->SetBranchAddress("muNDF", &muNDF, &b_muNDF);
   fChain->SetBranchAddress("muPhits", &muPhits, &b_muPhits);
   fChain->SetBranchAddress("muShits", &muShits, &b_muShits);
   fChain->SetBranchAddress("muLayersTr", &muLayersTr, &b_muLayersTr);
   fChain->SetBranchAddress("muLayersPix", &muLayersPix, &b_muLayersPix);
   fChain->SetBranchAddress("muD0E", &muD0E, &b_muD0E);
   fChain->SetBranchAddress("muDzVtxErr", &muDzVtxErr, &b_muDzVtxErr);
   fChain->SetBranchAddress("muKey", &muKey, &b_muKey);
   fChain->SetBranchAddress("muIsGlobal", &muIsGlobal, &b_muIsGlobal);
   fChain->SetBranchAddress("muIsPF", &muIsPF, &b_muIsPF);
   fChain->SetBranchAddress("muGlMuHits", &muGlMuHits, &b_muGlMuHits);
   fChain->SetBranchAddress("muGlChi2", &muGlChi2, &b_muGlChi2);
   fChain->SetBranchAddress("muGlNDF", &muGlNDF, &b_muGlNDF);
   fChain->SetBranchAddress("muGlMatchedStation", &muGlMatchedStation, &b_muGlMatchedStation);
   fChain->SetBranchAddress("muGlDzVtx", &muGlDzVtx, &b_muGlDzVtx);
   fChain->SetBranchAddress("muGlDxyVtx", &muGlDxyVtx, &b_muGlDxyVtx);
   fChain->SetBranchAddress("nMatchedStations", &nMatchedStations, &b_nMatchedStations);
   fChain->SetBranchAddress("muType", &muType, &b_muType);
   fChain->SetBranchAddress("muQual", &muQual, &b_muQual);
   fChain->SetBranchAddress("muTrack", &muTrack, &b_muTrack);
   fChain->SetBranchAddress("muNOverlap", &muNOverlap, &b_muNOverlap);
   fChain->SetBranchAddress("muNSharingSegWith", &muNSharingSegWith, &b_muNSharingSegWith);
   fChain->SetBranchAddress("mufHits", &mufHits, &b_mufHits);
   fChain->SetBranchAddress("muFirstBarrel", &muFirstBarrel, &b_muFirstBarrel);
   fChain->SetBranchAddress("muFirstEndCap", &muFirstEndCap, &b_muFirstEndCap);
   fChain->SetBranchAddress("muDzVtx", &muDzVtx, &b_muDzVtx);
   fChain->SetBranchAddress("muDxyVtx", &muDxyVtx, &b_muDxyVtx);
   fChain->SetBranchAddress("nMuMu", &nMuMu, &b_nMuMu);
   fChain->SetBranchAddress("MuMuMass", &MuMuMass, &b_MuMuMass);
   fChain->SetBranchAddress("MuMuPx", &MuMuPx, &b_MuMuPx);
   fChain->SetBranchAddress("MuMuPy", &MuMuPy, &b_MuMuPy);
   fChain->SetBranchAddress("MuMuPz", &MuMuPz, &b_MuMuPz);
   fChain->SetBranchAddress("MuMuVtx_CL", &MuMuVtx_CL, &b_MuMuVtx_CL);
   fChain->SetBranchAddress("MuMuVtx_Chi2", &MuMuVtx_Chi2, &b_MuMuVtx_Chi2);
   fChain->SetBranchAddress("MuMuDecayVtx_X", &MuMuDecayVtx_X, &b_MuMuDecayVtx_X);
   fChain->SetBranchAddress("MuMuDecayVtx_Y", &MuMuDecayVtx_Y, &b_MuMuDecayVtx_Y);
   fChain->SetBranchAddress("MuMuDecayVtx_Z", &MuMuDecayVtx_Z, &b_MuMuDecayVtx_Z);
   fChain->SetBranchAddress("MuMuDecayVtx_XE", &MuMuDecayVtx_XE, &b_MuMuDecayVtx_XE);
   fChain->SetBranchAddress("MuMuDecayVtx_YE", &MuMuDecayVtx_YE, &b_MuMuDecayVtx_YE);
   fChain->SetBranchAddress("MuMuDecayVtx_ZE", &MuMuDecayVtx_ZE, &b_MuMuDecayVtx_ZE);
   fChain->SetBranchAddress("mu1Idx", &mu1Idx, &b_mu1Idx);
   fChain->SetBranchAddress("mu2Idx", &mu2Idx, &b_mu2Idx);
   fChain->SetBranchAddress("mu1Px_MuMu", &mu1Px_MuMu, &b_mu1Px_MuMu);
   fChain->SetBranchAddress("mu1Py_MuMu", &mu1Py_MuMu, &b_mu1Py_MuMu);
   fChain->SetBranchAddress("mu1Pz_MuMu", &mu1Pz_MuMu, &b_mu1Pz_MuMu);
   fChain->SetBranchAddress("mu1Chi2_MuMu", &mu1Chi2_MuMu, &b_mu1Chi2_MuMu);
   fChain->SetBranchAddress("mu1NDF_MuMu", &mu1NDF_MuMu, &b_mu1NDF_MuMu);
   fChain->SetBranchAddress("mu2Px_MuMu", &mu2Px_MuMu, &b_mu2Px_MuMu);
   fChain->SetBranchAddress("mu2Py_MuMu", &mu2Py_MuMu, &b_mu2Py_MuMu);
   fChain->SetBranchAddress("mu2Pz_MuMu", &mu2Pz_MuMu, &b_mu2Pz_MuMu);
   fChain->SetBranchAddress("mu2Chi2_MuMu", &mu2Chi2_MuMu, &b_mu2Chi2_MuMu);
   fChain->SetBranchAddress("mu2NDF_MuMu", &mu2NDF_MuMu, &b_mu2NDF_MuMu);
   fChain->SetBranchAddress("MuMuType", &MuMuType, &b_MuMuType);
   fChain->SetBranchAddress("MuMuMuonTrigMatch", &MuMuMuonTrigMatch, &b_MuMuMuonTrigMatch);
   fChain->SetBranchAddress("PriVtxMuMuCorr_n", &PriVtxMuMuCorr_n, &b_PriVtxMuMuCorr_n);
   fChain->SetBranchAddress("PriVtxMuMuCorr_X", &PriVtxMuMuCorr_X, &b_PriVtxMuMuCorr_X);
   fChain->SetBranchAddress("PriVtxMuMuCorr_Y", &PriVtxMuMuCorr_Y, &b_PriVtxMuMuCorr_Y);
   fChain->SetBranchAddress("PriVtxMuMuCorr_Z", &PriVtxMuMuCorr_Z, &b_PriVtxMuMuCorr_Z);
   fChain->SetBranchAddress("PriVtxMuMuCorr_EX", &PriVtxMuMuCorr_EX, &b_PriVtxMuMuCorr_EX);
   fChain->SetBranchAddress("PriVtxMuMuCorr_EY", &PriVtxMuMuCorr_EY, &b_PriVtxMuMuCorr_EY);
   fChain->SetBranchAddress("PriVtxMuMuCorr_EZ", &PriVtxMuMuCorr_EZ, &b_PriVtxMuMuCorr_EZ);
   fChain->SetBranchAddress("PriVtxMuMuCorr_Chi2", &PriVtxMuMuCorr_Chi2, &b_PriVtxMuMuCorr_Chi2);
   fChain->SetBranchAddress("PriVtxMuMuCorr_CL", &PriVtxMuMuCorr_CL, &b_PriVtxMuMuCorr_CL);
   fChain->SetBranchAddress("PriVtxMuMuCorr_tracks", &PriVtxMuMuCorr_tracks, &b_PriVtxMuMuCorr_tracks);
   fChain->SetBranchAddress("nTrk_afterMuMu", &nTrk_afterMuMu, &b_nTrk_afterMuMu);
   fChain->SetBranchAddress("nKK", &nKK, &b_nKK);
   fChain->SetBranchAddress("KKMass", &KKMass, &b_KKMass);
   fChain->SetBranchAddress("KKPx", &KKPx, &b_KKPx);
   fChain->SetBranchAddress("KKPy", &KKPy, &b_KKPy);
   fChain->SetBranchAddress("KKPz", &KKPz, &b_KKPz);
   fChain->SetBranchAddress("KKVtx_CL", &KKVtx_CL, &b_KKVtx_CL);
   fChain->SetBranchAddress("KKVtx_Chi2", &KKVtx_Chi2, &b_KKVtx_Chi2);
   fChain->SetBranchAddress("KKDecayVtx_X", &KKDecayVtx_X, &b_KKDecayVtx_X);
   fChain->SetBranchAddress("KKDecayVtx_Y", &KKDecayVtx_Y, &b_KKDecayVtx_Y);
   fChain->SetBranchAddress("KKDecayVtx_Z", &KKDecayVtx_Z, &b_KKDecayVtx_Z);
   fChain->SetBranchAddress("KKDecayVtx_XE", &KKDecayVtx_XE, &b_KKDecayVtx_XE);
   fChain->SetBranchAddress("KKDecayVtx_YE", &KKDecayVtx_YE, &b_KKDecayVtx_YE);
   fChain->SetBranchAddress("KKDecayVtx_ZE", &KKDecayVtx_ZE, &b_KKDecayVtx_ZE);
   fChain->SetBranchAddress("ka1Idx", &ka1Idx, &b_ka1Idx);
   fChain->SetBranchAddress("ka2Idx", &ka2Idx, &b_ka2Idx);
   fChain->SetBranchAddress("ka1Px_KK", &ka1Px_KK, &b_ka1Px_KK);
   fChain->SetBranchAddress("ka1Py_KK", &ka1Py_KK, &b_ka1Py_KK);
   fChain->SetBranchAddress("ka1Pz_KK", &ka1Pz_KK, &b_ka1Pz_KK);
   fChain->SetBranchAddress("ka1Chi2_KK", &ka1Chi2_KK, &b_ka1Chi2_KK);
   fChain->SetBranchAddress("ka1NDF_KK", &ka1NDF_KK, &b_ka1NDF_KK);
   fChain->SetBranchAddress("ka2Px_KK", &ka2Px_KK, &b_ka2Px_KK);
   fChain->SetBranchAddress("ka2Py_KK", &ka2Py_KK, &b_ka2Py_KK);
   fChain->SetBranchAddress("ka2Pz_KK", &ka2Pz_KK, &b_ka2Pz_KK);
   fChain->SetBranchAddress("ka2Chi2_KK", &ka2Chi2_KK, &b_ka2Chi2_KK);
   fChain->SetBranchAddress("ka2NDF_KK", &ka2NDF_KK, &b_ka2NDF_KK);
   fChain->SetBranchAddress("DR_MuMu_K1", &DR_MuMu_K1, &b_DR_MuMu_K1);
   fChain->SetBranchAddress("DR_MuMu_K2", &DR_MuMu_K2, &b_DR_MuMu_K2);
   fChain->SetBranchAddress("DR_MuMuKK_K1", &DR_MuMuKK_K1, &b_DR_MuMuKK_K1);
   fChain->SetBranchAddress("DR_MuMuKK_K2", &DR_MuMuKK_K2, &b_DR_MuMuKK_K2);
   fChain->SetBranchAddress("nX", &nX, &b_nX);
   fChain->SetBranchAddress("nX_pre0", &nX_pre0, &b_nX_pre0);
   fChain->SetBranchAddress("nX_pre1", &nX_pre1, &b_nX_pre1);
   fChain->SetBranchAddress("nX_pre2", &nX_pre2, &b_nX_pre2);
   fChain->SetBranchAddress("nX_pre3", &nX_pre3, &b_nX_pre3);
   fChain->SetBranchAddress("nX_pre4", &nX_pre4, &b_nX_pre4);
   fChain->SetBranchAddress("nX_pre5", &nX_pre5, &b_nX_pre5);
   fChain->SetBranchAddress("nX_pre6", &nX_pre6, &b_nX_pre6);
   fChain->SetBranchAddress("nX_pre7", &nX_pre7, &b_nX_pre7);
   fChain->SetBranchAddress("nX_pre8", &nX_pre8, &b_nX_pre8);
   fChain->SetBranchAddress("nX_pre9", &nX_pre9, &b_nX_pre9);
   fChain->SetBranchAddress("nX_pre10", &nX_pre10, &b_nX_pre10);
   fChain->SetBranchAddress("nX_pre11", &nX_pre11, &b_nX_pre11);
   fChain->SetBranchAddress("nX_pre12", &nX_pre12, &b_nX_pre12);
   fChain->SetBranchAddress("nX_pre13", &nX_pre13, &b_nX_pre13);
   fChain->SetBranchAddress("nX_pre14", &nX_pre14, &b_nX_pre14);
   fChain->SetBranchAddress("nX_pre15", &nX_pre15, &b_nX_pre15);
   fChain->SetBranchAddress("XMass", &XMass, &b_XMass);
   fChain->SetBranchAddress("XPx", &XPx, &b_XPx);
   fChain->SetBranchAddress("XPy", &XPy, &b_XPy);
   fChain->SetBranchAddress("XPz", &XPz, &b_XPz);
   fChain->SetBranchAddress("XPxE", &XPxE, &b_XPxE);
   fChain->SetBranchAddress("XPyE", &XPyE, &b_XPyE);
   fChain->SetBranchAddress("XPzE", &XPzE, &b_XPzE);
   fChain->SetBranchAddress("XVtx_CL", &XVtx_CL, &b_XVtx_CL);
   fChain->SetBranchAddress("XVtx_Chi2", &XVtx_Chi2, &b_XVtx_Chi2);
   fChain->SetBranchAddress("XDecayVtx_X", &XDecayVtx_X, &b_XDecayVtx_X);
   fChain->SetBranchAddress("XDecayVtx_Y", &XDecayVtx_Y, &b_XDecayVtx_Y);
   fChain->SetBranchAddress("XDecayVtx_Z", &XDecayVtx_Z, &b_XDecayVtx_Z);
   fChain->SetBranchAddress("XDecayVtx_XE", &XDecayVtx_XE, &b_XDecayVtx_XE);
   fChain->SetBranchAddress("XDecayVtx_YE", &XDecayVtx_YE, &b_XDecayVtx_YE);
   fChain->SetBranchAddress("XDecayVtx_ZE", &XDecayVtx_ZE, &b_XDecayVtx_ZE);
   fChain->SetBranchAddress("XCosAlphaBS", &XCosAlphaBS, &b_XCosAlphaBS);
   fChain->SetBranchAddress("XCosAlpha3DBS", &XCosAlpha3DBS, &b_XCosAlpha3DBS);
   fChain->SetBranchAddress("XCTauBS", &XCTauBS, &b_XCTauBS);
   fChain->SetBranchAddress("XCTauBSE", &XCTauBSE, &b_XCTauBSE);
   fChain->SetBranchAddress("XLxyBS", &XLxyBS, &b_XLxyBS);
   fChain->SetBranchAddress("XLxyBSE", &XLxyBSE, &b_XLxyBSE);
   fChain->SetBranchAddress("XLxyzBS", &XLxyzBS, &b_XLxyzBS);
   fChain->SetBranchAddress("XLxyzBSE", &XLxyzBSE, &b_XLxyzBSE);
   fChain->SetBranchAddress("XCosAlphaPV", &XCosAlphaPV, &b_XCosAlphaPV);
   fChain->SetBranchAddress("XCosAlpha3DPV", &XCosAlpha3DPV, &b_XCosAlpha3DPV);
   fChain->SetBranchAddress("XCTauPV", &XCTauPV, &b_XCTauPV);
   fChain->SetBranchAddress("XCTauPVE", &XCTauPVE, &b_XCTauPVE);
   fChain->SetBranchAddress("XLxyPV", &XLxyPV, &b_XLxyPV);
   fChain->SetBranchAddress("XLxyPVE", &XLxyPVE, &b_XLxyPVE);
   fChain->SetBranchAddress("XLxyzPV", &XLxyzPV, &b_XLxyzPV);
   fChain->SetBranchAddress("XLxyzPVE", &XLxyzPVE, &b_XLxyzPVE);
   fChain->SetBranchAddress("PriVtx_XCosAlpha_n", &PriVtx_XCosAlpha_n, &b_PriVtx_XCosAlpha_n);
   fChain->SetBranchAddress("PriVtx_XCosAlpha_X", &PriVtx_XCosAlpha_X, &b_PriVtx_XCosAlpha_X);
   fChain->SetBranchAddress("PriVtx_XCosAlpha_Y", &PriVtx_XCosAlpha_Y, &b_PriVtx_XCosAlpha_Y);
   fChain->SetBranchAddress("PriVtx_XCosAlpha_Z", &PriVtx_XCosAlpha_Z, &b_PriVtx_XCosAlpha_Z);
   fChain->SetBranchAddress("PriVtx_XCosAlpha_EX", &PriVtx_XCosAlpha_EX, &b_PriVtx_XCosAlpha_EX);
   fChain->SetBranchAddress("PriVtx_XCosAlpha_EY", &PriVtx_XCosAlpha_EY, &b_PriVtx_XCosAlpha_EY);
   fChain->SetBranchAddress("PriVtx_XCosAlpha_EZ", &PriVtx_XCosAlpha_EZ, &b_PriVtx_XCosAlpha_EZ);
   fChain->SetBranchAddress("PriVtx_XCosAlpha_Chi2", &PriVtx_XCosAlpha_Chi2, &b_PriVtx_XCosAlpha_Chi2);
   fChain->SetBranchAddress("PriVtx_XCosAlpha_CL", &PriVtx_XCosAlpha_CL, &b_PriVtx_XCosAlpha_CL);
   fChain->SetBranchAddress("PriVtx_XCosAlpha_tracks", &PriVtx_XCosAlpha_tracks, &b_PriVtx_XCosAlpha_tracks);
   fChain->SetBranchAddress("XCosAlphaPVCosAlpha", &XCosAlphaPVCosAlpha, &b_XCosAlphaPVCosAlpha);
   fChain->SetBranchAddress("XCosAlpha3DPVCosAlpha", &XCosAlpha3DPVCosAlpha, &b_XCosAlpha3DPVCosAlpha);
   fChain->SetBranchAddress("XCTauPVCosAlpha", &XCTauPVCosAlpha, &b_XCTauPVCosAlpha);
   fChain->SetBranchAddress("XCTauPVCosAlphaE", &XCTauPVCosAlphaE, &b_XCTauPVCosAlphaE);
   fChain->SetBranchAddress("XLxyPVCosAlpha", &XLxyPVCosAlpha, &b_XLxyPVCosAlpha);
   fChain->SetBranchAddress("XLxyPVCosAlphaE", &XLxyPVCosAlphaE, &b_XLxyPVCosAlphaE);
   fChain->SetBranchAddress("XLxyzPVCosAlpha", &XLxyzPVCosAlpha, &b_XLxyzPVCosAlpha);
   fChain->SetBranchAddress("XLxyzPVCosAlphaE", &XLxyzPVCosAlphaE, &b_XLxyzPVCosAlphaE);
   fChain->SetBranchAddress("PriVtx_XCosAlpha3D_n", &PriVtx_XCosAlpha3D_n, &b_PriVtx_XCosAlpha3D_n);
   fChain->SetBranchAddress("PriVtx_XCosAlpha3D_X", &PriVtx_XCosAlpha3D_X, &b_PriVtx_XCosAlpha3D_X);
   fChain->SetBranchAddress("PriVtx_XCosAlpha3D_Y", &PriVtx_XCosAlpha3D_Y, &b_PriVtx_XCosAlpha3D_Y);
   fChain->SetBranchAddress("PriVtx_XCosAlpha3D_Z", &PriVtx_XCosAlpha3D_Z, &b_PriVtx_XCosAlpha3D_Z);
   fChain->SetBranchAddress("PriVtx_XCosAlpha3D_EX", &PriVtx_XCosAlpha3D_EX, &b_PriVtx_XCosAlpha3D_EX);
   fChain->SetBranchAddress("PriVtx_XCosAlpha3D_EY", &PriVtx_XCosAlpha3D_EY, &b_PriVtx_XCosAlpha3D_EY);
   fChain->SetBranchAddress("PriVtx_XCosAlpha3D_EZ", &PriVtx_XCosAlpha3D_EZ, &b_PriVtx_XCosAlpha3D_EZ);
   fChain->SetBranchAddress("PriVtx_XCosAlpha3D_Chi2", &PriVtx_XCosAlpha3D_Chi2, &b_PriVtx_XCosAlpha3D_Chi2);
   fChain->SetBranchAddress("PriVtx_XCosAlpha3D_CL", &PriVtx_XCosAlpha3D_CL, &b_PriVtx_XCosAlpha3D_CL);
   fChain->SetBranchAddress("PriVtx_XCosAlpha3D_tracks", &PriVtx_XCosAlpha3D_tracks, &b_PriVtx_XCosAlpha3D_tracks);
   fChain->SetBranchAddress("XCosAlphaPVCosAlpha3D", &XCosAlphaPVCosAlpha3D, &b_XCosAlphaPVCosAlpha3D);
   fChain->SetBranchAddress("XCosAlpha3DPVCosAlpha3D", &XCosAlpha3DPVCosAlpha3D, &b_XCosAlpha3DPVCosAlpha3D);
   fChain->SetBranchAddress("XCTauPVCosAlpha3D", &XCTauPVCosAlpha3D, &b_XCTauPVCosAlpha3D);
   fChain->SetBranchAddress("XCTauPVCosAlpha3DE", &XCTauPVCosAlpha3DE, &b_XCTauPVCosAlpha3DE);
   fChain->SetBranchAddress("XLxyPVCosAlpha3D", &XLxyPVCosAlpha3D, &b_XLxyPVCosAlpha3D);
   fChain->SetBranchAddress("XLxyPVCosAlpha3DE", &XLxyPVCosAlpha3DE, &b_XLxyPVCosAlpha3DE);
   fChain->SetBranchAddress("XLxyzPVCosAlpha3D", &XLxyzPVCosAlpha3D, &b_XLxyzPVCosAlpha3D);
   fChain->SetBranchAddress("XLxyzPVCosAlpha3DE", &XLxyzPVCosAlpha3DE, &b_XLxyzPVCosAlpha3DE);
   fChain->SetBranchAddress("XLessPV_tracksPtSq", &XLessPV_tracksPtSq, &b_XLessPV_tracksPtSq);
   fChain->SetBranchAddress("XLessPV_4tracksPtSq", &XLessPV_4tracksPtSq, &b_XLessPV_4tracksPtSq);
   fChain->SetBranchAddress("PriVtxXLess_n", &PriVtxXLess_n, &b_PriVtxXLess_n);
   fChain->SetBranchAddress("PriVtxXLess_X", &PriVtxXLess_X, &b_PriVtxXLess_X);
   fChain->SetBranchAddress("PriVtxXLess_Y", &PriVtxXLess_Y, &b_PriVtxXLess_Y);
   fChain->SetBranchAddress("PriVtxXLess_Z", &PriVtxXLess_Z, &b_PriVtxXLess_Z);
   fChain->SetBranchAddress("PriVtxXLess_EX", &PriVtxXLess_EX, &b_PriVtxXLess_EX);
   fChain->SetBranchAddress("PriVtxXLess_EY", &PriVtxXLess_EY, &b_PriVtxXLess_EY);
   fChain->SetBranchAddress("PriVtxXLess_EZ", &PriVtxXLess_EZ, &b_PriVtxXLess_EZ);
   fChain->SetBranchAddress("PriVtxXLess_Chi2", &PriVtxXLess_Chi2, &b_PriVtxXLess_Chi2);
   fChain->SetBranchAddress("PriVtxXLess_CL", &PriVtxXLess_CL, &b_PriVtxXLess_CL);
   fChain->SetBranchAddress("PriVtxXLess_tracks", &PriVtxXLess_tracks, &b_PriVtxXLess_tracks);
   fChain->SetBranchAddress("XCosAlphaXLessPV", &XCosAlphaXLessPV, &b_XCosAlphaXLessPV);
   fChain->SetBranchAddress("XCosAlpha3DXLessPV", &XCosAlpha3DXLessPV, &b_XCosAlpha3DXLessPV);
   fChain->SetBranchAddress("XCTauXLessPV", &XCTauXLessPV, &b_XCTauXLessPV);
   fChain->SetBranchAddress("XCTauXLessPVE", &XCTauXLessPVE, &b_XCTauXLessPVE);
   fChain->SetBranchAddress("XLxyXLessPV", &XLxyXLessPV, &b_XLxyXLessPV);
   fChain->SetBranchAddress("XLxyXLessPVE", &XLxyXLessPVE, &b_XLxyXLessPVE);
   fChain->SetBranchAddress("XLxyzXLessPV", &XLxyzXLessPV, &b_XLxyzXLessPV);
   fChain->SetBranchAddress("XLxyzXLessPVE", &XLxyzXLessPVE, &b_XLxyzXLessPVE);
   fChain->SetBranchAddress("PriVtxXLess_XCosAlpha_n", &PriVtxXLess_XCosAlpha_n, &b_PriVtxXLess_XCosAlpha_n);
   fChain->SetBranchAddress("PriVtxXLess_XCosAlpha_X", &PriVtxXLess_XCosAlpha_X, &b_PriVtxXLess_XCosAlpha_X);
   fChain->SetBranchAddress("PriVtxXLess_XCosAlpha_Y", &PriVtxXLess_XCosAlpha_Y, &b_PriVtxXLess_XCosAlpha_Y);
   fChain->SetBranchAddress("PriVtxXLess_XCosAlpha_Z", &PriVtxXLess_XCosAlpha_Z, &b_PriVtxXLess_XCosAlpha_Z);
   fChain->SetBranchAddress("PriVtxXLess_XCosAlpha_EX", &PriVtxXLess_XCosAlpha_EX, &b_PriVtxXLess_XCosAlpha_EX);
   fChain->SetBranchAddress("PriVtxXLess_XCosAlpha_EY", &PriVtxXLess_XCosAlpha_EY, &b_PriVtxXLess_XCosAlpha_EY);
   fChain->SetBranchAddress("PriVtxXLess_XCosAlpha_EZ", &PriVtxXLess_XCosAlpha_EZ, &b_PriVtxXLess_XCosAlpha_EZ);
   fChain->SetBranchAddress("PriVtxXLess_XCosAlpha_Chi2", &PriVtxXLess_XCosAlpha_Chi2, &b_PriVtxXLess_XCosAlpha_Chi2);
   fChain->SetBranchAddress("PriVtxXLess_XCosAlpha_CL", &PriVtxXLess_XCosAlpha_CL, &b_PriVtxXLess_XCosAlpha_CL);
   fChain->SetBranchAddress("PriVtxXLess_XCosAlpha_tracks", &PriVtxXLess_XCosAlpha_tracks, &b_PriVtxXLess_XCosAlpha_tracks);
   fChain->SetBranchAddress("XCosAlphaXLessPVCosAlpha", &XCosAlphaXLessPVCosAlpha, &b_XCosAlphaXLessPVCosAlpha);
   fChain->SetBranchAddress("XCosAlpha3DXLessPVCosAlpha", &XCosAlpha3DXLessPVCosAlpha, &b_XCosAlpha3DXLessPVCosAlpha);
   fChain->SetBranchAddress("XCTauXLessPVCosAlpha", &XCTauXLessPVCosAlpha, &b_XCTauXLessPVCosAlpha);
   fChain->SetBranchAddress("XCTauXLessPVCosAlphaE", &XCTauXLessPVCosAlphaE, &b_XCTauXLessPVCosAlphaE);
   fChain->SetBranchAddress("XLxyXLessPVCosAlpha", &XLxyXLessPVCosAlpha, &b_XLxyXLessPVCosAlpha);
   fChain->SetBranchAddress("XLxyXLessPVCosAlphaE", &XLxyXLessPVCosAlphaE, &b_XLxyXLessPVCosAlphaE);
   fChain->SetBranchAddress("XLxyzXLessPVCosAlpha", &XLxyzXLessPVCosAlpha, &b_XLxyzXLessPVCosAlpha);
   fChain->SetBranchAddress("XLxyzXLessPVCosAlphaE", &XLxyzXLessPVCosAlphaE, &b_XLxyzXLessPVCosAlphaE);
   fChain->SetBranchAddress("PriVtxXLess_XCosAlpha3D_n", &PriVtxXLess_XCosAlpha3D_n, &b_PriVtxXLess_XCosAlpha3D_n);
   fChain->SetBranchAddress("PriVtxXLess_XCosAlpha3D_X", &PriVtxXLess_XCosAlpha3D_X, &b_PriVtxXLess_XCosAlpha3D_X);
   fChain->SetBranchAddress("PriVtxXLess_XCosAlpha3D_Y", &PriVtxXLess_XCosAlpha3D_Y, &b_PriVtxXLess_XCosAlpha3D_Y);
   fChain->SetBranchAddress("PriVtxXLess_XCosAlpha3D_Z", &PriVtxXLess_XCosAlpha3D_Z, &b_PriVtxXLess_XCosAlpha3D_Z);
   fChain->SetBranchAddress("PriVtxXLess_XCosAlpha3D_EX", &PriVtxXLess_XCosAlpha3D_EX, &b_PriVtxXLess_XCosAlpha3D_EX);
   fChain->SetBranchAddress("PriVtxXLess_XCosAlpha3D_EY", &PriVtxXLess_XCosAlpha3D_EY, &b_PriVtxXLess_XCosAlpha3D_EY);
   fChain->SetBranchAddress("PriVtxXLess_XCosAlpha3D_EZ", &PriVtxXLess_XCosAlpha3D_EZ, &b_PriVtxXLess_XCosAlpha3D_EZ);
   fChain->SetBranchAddress("PriVtxXLess_XCosAlpha3D_Chi2", &PriVtxXLess_XCosAlpha3D_Chi2, &b_PriVtxXLess_XCosAlpha3D_Chi2);
   fChain->SetBranchAddress("PriVtxXLess_XCosAlpha3D_CL", &PriVtxXLess_XCosAlpha3D_CL, &b_PriVtxXLess_XCosAlpha3D_CL);
   fChain->SetBranchAddress("PriVtxXLess_XCosAlpha3D_tracks", &PriVtxXLess_XCosAlpha3D_tracks, &b_PriVtxXLess_XCosAlpha3D_tracks);
   fChain->SetBranchAddress("XCosAlphaXLessPVCosAlpha3D", &XCosAlphaXLessPVCosAlpha3D, &b_XCosAlphaXLessPVCosAlpha3D);
   fChain->SetBranchAddress("XCosAlpha3DXLessPVCosAlpha3D", &XCosAlpha3DXLessPVCosAlpha3D, &b_XCosAlpha3DXLessPVCosAlpha3D);
   fChain->SetBranchAddress("XCTauXLessPVCosAlpha3D", &XCTauXLessPVCosAlpha3D, &b_XCTauXLessPVCosAlpha3D);
   fChain->SetBranchAddress("XCTauXLessPVCosAlpha3DE", &XCTauXLessPVCosAlpha3DE, &b_XCTauXLessPVCosAlpha3DE);
   fChain->SetBranchAddress("XLxyXLessPVCosAlpha3D", &XLxyXLessPVCosAlpha3D, &b_XLxyXLessPVCosAlpha3D);
   fChain->SetBranchAddress("XLxyXLessPVCosAlpha3DE", &XLxyXLessPVCosAlpha3DE, &b_XLxyXLessPVCosAlpha3DE);
   fChain->SetBranchAddress("XLxyzXLessPVCosAlpha3D", &XLxyzXLessPVCosAlpha3D, &b_XLxyzXLessPVCosAlpha3D);
   fChain->SetBranchAddress("XLxyzXLessPVCosAlpha3DE", &XLxyzXLessPVCosAlpha3DE, &b_XLxyzXLessPVCosAlpha3DE);
   fChain->SetBranchAddress("PriVtxXCorr_n", &PriVtxXCorr_n, &b_PriVtxXCorr_n);
   fChain->SetBranchAddress("PriVtxXCorr_X", &PriVtxXCorr_X, &b_PriVtxXCorr_X);
   fChain->SetBranchAddress("PriVtxXCorr_Y", &PriVtxXCorr_Y, &b_PriVtxXCorr_Y);
   fChain->SetBranchAddress("PriVtxXCorr_Z", &PriVtxXCorr_Z, &b_PriVtxXCorr_Z);
   fChain->SetBranchAddress("PriVtxXCorr_EX", &PriVtxXCorr_EX, &b_PriVtxXCorr_EX);
   fChain->SetBranchAddress("PriVtxXCorr_EY", &PriVtxXCorr_EY, &b_PriVtxXCorr_EY);
   fChain->SetBranchAddress("PriVtxXCorr_EZ", &PriVtxXCorr_EZ, &b_PriVtxXCorr_EZ);
   fChain->SetBranchAddress("PriVtxXCorr_Chi2", &PriVtxXCorr_Chi2, &b_PriVtxXCorr_Chi2);
   fChain->SetBranchAddress("PriVtxXCorr_CL", &PriVtxXCorr_CL, &b_PriVtxXCorr_CL);
   fChain->SetBranchAddress("PriVtxXCorr_tracks", &PriVtxXCorr_tracks, &b_PriVtxXCorr_tracks);
   fChain->SetBranchAddress("XCosAlphaPVX", &XCosAlphaPVX, &b_XCosAlphaPVX);
   fChain->SetBranchAddress("XCTauPVX", &XCTauPVX, &b_XCTauPVX);
   fChain->SetBranchAddress("XCTauPVXE", &XCTauPVXE, &b_XCTauPVXE);
   fChain->SetBranchAddress("XLxyPVX", &XLxyPVX, &b_XLxyPVX);
   fChain->SetBranchAddress("XLxyzPVX", &XLxyzPVX, &b_XLxyzPVX);
   fChain->SetBranchAddress("XCTauPVX_3D", &XCTauPVX_3D, &b_XCTauPVX_3D);
   fChain->SetBranchAddress("XCTauPVX_3D_err", &XCTauPVX_3D_err, &b_XCTauPVX_3D_err);
   fChain->SetBranchAddress("kaon1_dxy_PV", &kaon1_dxy_PV, &b_kaon1_dxy_PV);
   fChain->SetBranchAddress("kaon1_dz_PV", &kaon1_dz_PV, &b_kaon1_dz_PV);
   fChain->SetBranchAddress("kaon2_dxy_PV", &kaon2_dxy_PV, &b_kaon2_dxy_PV);
   fChain->SetBranchAddress("kaon2_dz_PV", &kaon2_dz_PV, &b_kaon2_dz_PV);
   fChain->SetBranchAddress("kaon1_dxy_BS", &kaon1_dxy_BS, &b_kaon1_dxy_BS);
   fChain->SetBranchAddress("kaon1_dz_BS", &kaon1_dz_BS, &b_kaon1_dz_BS);
   fChain->SetBranchAddress("kaon2_dxy_BS", &kaon2_dxy_BS, &b_kaon2_dxy_BS);
   fChain->SetBranchAddress("kaon2_dz_BS", &kaon2_dz_BS, &b_kaon2_dz_BS);
   fChain->SetBranchAddress("kaon1_dxy_XLessPV", &kaon1_dxy_XLessPV, &b_kaon1_dxy_XLessPV);
   fChain->SetBranchAddress("kaon1_dz_XLessPV", &kaon1_dz_XLessPV, &b_kaon1_dz_XLessPV);
   fChain->SetBranchAddress("kaon2_dxy_XLessPV", &kaon2_dxy_XLessPV, &b_kaon2_dxy_XLessPV);
   fChain->SetBranchAddress("kaon2_dz_XLessPV", &kaon2_dz_XLessPV, &b_kaon2_dz_XLessPV);
   fChain->SetBranchAddress("kaon1_dxyE", &kaon1_dxyE, &b_kaon1_dxyE);
   fChain->SetBranchAddress("kaon1_dzE", &kaon1_dzE, &b_kaon1_dzE);
   fChain->SetBranchAddress("kaon2_dxyE", &kaon2_dxyE, &b_kaon2_dxyE);
   fChain->SetBranchAddress("kaon2_dzE", &kaon2_dzE, &b_kaon2_dzE);
   fChain->SetBranchAddress("XMuMuIdx", &XMuMuIdx, &b_XMuMuIdx);
   fChain->SetBranchAddress("XKaon1Idx", &XKaon1Idx, &b_XKaon1Idx);
   fChain->SetBranchAddress("XKaon2Idx", &XKaon2Idx, &b_XKaon2Idx);
   fChain->SetBranchAddress("Kaon1FromPV", &Kaon1FromPV, &b_Kaon1FromPV);
   fChain->SetBranchAddress("Kaon2FromPV", &Kaon2FromPV, &b_Kaon2FromPV);
   fChain->SetBranchAddress("Muon1Px_MuMuKK", &Muon1Px_MuMuKK, &b_Muon1Px_MuMuKK);
   fChain->SetBranchAddress("Muon1Py_MuMuKK", &Muon1Py_MuMuKK, &b_Muon1Py_MuMuKK);
   fChain->SetBranchAddress("Muon1Pz_MuMuKK", &Muon1Pz_MuMuKK, &b_Muon1Pz_MuMuKK);
   fChain->SetBranchAddress("Muon1E_MuMuKK", &Muon1E_MuMuKK, &b_Muon1E_MuMuKK);
   fChain->SetBranchAddress("Muon2Px_MuMuKK", &Muon2Px_MuMuKK, &b_Muon2Px_MuMuKK);
   fChain->SetBranchAddress("Muon2Py_MuMuKK", &Muon2Py_MuMuKK, &b_Muon2Py_MuMuKK);
   fChain->SetBranchAddress("Muon2Pz_MuMuKK", &Muon2Pz_MuMuKK, &b_Muon2Pz_MuMuKK);
   fChain->SetBranchAddress("Muon2E_MuMuKK", &Muon2E_MuMuKK, &b_Muon2E_MuMuKK);
   fChain->SetBranchAddress("Kaon1Px_MuMuKK", &Kaon1Px_MuMuKK, &b_Kaon1Px_MuMuKK);
   fChain->SetBranchAddress("Kaon1Py_MuMuKK", &Kaon1Py_MuMuKK, &b_Kaon1Py_MuMuKK);
   fChain->SetBranchAddress("Kaon1Pz_MuMuKK", &Kaon1Pz_MuMuKK, &b_Kaon1Pz_MuMuKK);
   fChain->SetBranchAddress("Kaon1E_MuMuKK", &Kaon1E_MuMuKK, &b_Kaon1E_MuMuKK);
   fChain->SetBranchAddress("kaon1_nsigdedx", &kaon1_nsigdedx, &b_kaon1_nsigdedx);
   fChain->SetBranchAddress("kaon1_dedx", &kaon1_dedx, &b_kaon1_dedx);
   fChain->SetBranchAddress("kaon1_dedxMass", &kaon1_dedxMass, &b_kaon1_dedxMass);
   fChain->SetBranchAddress("kaon1_theo", &kaon1_theo, &b_kaon1_theo);
   fChain->SetBranchAddress("kaon1_sigma", &kaon1_sigma, &b_kaon1_sigma);
   fChain->SetBranchAddress("kaon1_dedx_byHits", &kaon1_dedx_byHits, &b_kaon1_dedx_byHits);
   fChain->SetBranchAddress("kaon1_dedxErr_byHits", &kaon1_dedxErr_byHits, &b_kaon1_dedxErr_byHits);
   fChain->SetBranchAddress("kaon1_saturMeas_byHits", &kaon1_saturMeas_byHits, &b_kaon1_saturMeas_byHits);
   fChain->SetBranchAddress("kaon1_Meas_byHits", &kaon1_Meas_byHits, &b_kaon1_Meas_byHits);
   fChain->SetBranchAddress("Kaon2Px_MuMuKK", &Kaon2Px_MuMuKK, &b_Kaon2Px_MuMuKK);
   fChain->SetBranchAddress("Kaon2Py_MuMuKK", &Kaon2Py_MuMuKK, &b_Kaon2Py_MuMuKK);
   fChain->SetBranchAddress("Kaon2Pz_MuMuKK", &Kaon2Pz_MuMuKK, &b_Kaon2Pz_MuMuKK);
   fChain->SetBranchAddress("Kaon2E_MuMuKK", &Kaon2E_MuMuKK, &b_Kaon2E_MuMuKK);
   fChain->SetBranchAddress("kaon2_nsigdedx", &kaon2_nsigdedx, &b_kaon2_nsigdedx);
   fChain->SetBranchAddress("kaon2_dedx", &kaon2_dedx, &b_kaon2_dedx);
   fChain->SetBranchAddress("kaon2_dedxMass", &kaon2_dedxMass, &b_kaon2_dedxMass);
   fChain->SetBranchAddress("kaon2_theo", &kaon2_theo, &b_kaon2_theo);
   fChain->SetBranchAddress("kaon2_sigma", &kaon2_sigma, &b_kaon2_sigma);
   fChain->SetBranchAddress("kaon2_dedx_byHits", &kaon2_dedx_byHits, &b_kaon2_dedx_byHits);
   fChain->SetBranchAddress("kaon2_dedxErr_byHits", &kaon2_dedxErr_byHits, &b_kaon2_dedxErr_byHits);
   fChain->SetBranchAddress("kaon2_saturMeas_byHits", &kaon2_saturMeas_byHits, &b_kaon2_saturMeas_byHits);
   fChain->SetBranchAddress("kaon2_Meas_byHits", &kaon2_Meas_byHits, &b_kaon2_Meas_byHits);


     // outTree = new TTree("outTree","outTree");

    //  //cw values
    // X_mass = 0;
    // kk_mass = 0;
    // mumu_mass = 0;
    // X_LFly = 0;
    // X_pt = 0;
    // X_eta = 0;
    // X_vtx = 0;
    // X_cosAlpha = 0;
        //
        // outTree->SetBranchAddress("X_mass_out", &X_mass, &X_mass_b);
        // outTree->SetBranchAddress("kk_mass_out", &kk_mass, &kk_mass_b);
        // outTree->SetBranchAddress("mumu_mass_out", &mumu_mass, &mumu_mass_b);
        // outTree->SetBranchAddress("X_LFly_out", &X_LFly, &X_LFly_b);
        // outTree->SetBranchAddress("X_eta_out", &X_eta, &X_eta_b);
        // outTree->SetBranchAddress("X_vtx_out", &X_vtx, &X_vtx_b);
        // outTree->SetBranchAddress("X_cosAlpha_out", &X_cosAlpha, &X_cosAlpha_b);
        //


}

Bool_t mumukk::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif
