// -*- C++ -*-
//
// Package:    MuMuKKPAT
// Class:      MuMuKKPAT
//
/**\class MuMuKKPAT MuMuKKPAT.cc myAnalyzers/MuMuKKPAT/src/MuMuKKPAT.cc

   Description: <one line class summary>
   Make rootTuple for JPsiKK reconstruction

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Keith Ulmer   keith.ulmer@colorado.edu
//
//

#ifndef _MuMuKKPAT_h
#define _MuMuKKPAT_h

// system include files
#include <memory>

/// user include files
#include "../interface/VertexReProducer.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "RecoVertex/V0Producer/interface/V0Producer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "TFile.h"
#include "TTree.h"

#include <vector>
#include <string>
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"

#include "DataFormats/Math/interface/LorentzVector.h"
//#include "DataFormats/PatCandidates/interface/GenericParticle.h" // for namespace pat

///
/// class decleration
///

using std::vector;
using namespace edm;
using namespace reco;
using namespace std;

class MuMuKKPAT : public edm::EDAnalyzer {
public:
  explicit MuMuKKPAT(const edm::ParameterSet&);
  ~MuMuKKPAT();

private:
  virtual void beginJob() ;
  virtual void beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  InvariantMassFromVertex massCalculator;
  const reco::DeDxDataValueMap *energyLoss;
  Int_t iexception_dedx;
  /// dE/dx hits
  edm::ValueMap<reco::DeDxData> dEdxTrack, dEdxTrack_Kaon;
  bool isAbHadron(int pdgID);
  bool isAMixedbHadron(int pdgID, int momPdgID);
  std::pair<int, float> findCandMCInfo(reco::GenParticleRef genCand);
  virtual double getSigmaOfLogdEdx(double logde);
  virtual float  getEnergyLoss(const reco::TrackRef & track);
  virtual double nsigmaofdedx(const reco::TrackRef & track, double & theo, double & sigma);
  virtual double getLogdEdx(double bg);
  virtual double GetMass(const reco::TrackRef & track);
  bool isSameMuon(const reco::Muon &mu1, const reco::Muon &mu2) const ;
  template<typename T> bool isBetterMuon(const T &mu1, const T &mu2) const ;

  /// ----------member data ---------------------------
  std::string proccessName_;
  HLTConfigProvider hltConfig_;

  edm::InputTag hlTriggerResults_;
  std::map<std::string,int> *HLTTrig; /// HLT trigger prescale for accepted paths

  edm::InputTag inputGEN_;
  std::string vtxSample;
  bool doData, doMC;
  int  MCParticle;
  bool MCExclusiveDecay;
  int  MCMother, MCDaughtersN;
  bool doMuMuMassConst;
  bool skipJPsi, skipPsi2S;
  int MuMinPixHits, MuMinSiHits;
  double MuMaxNormChi;
  double MuMaxD0;
  bool sharedFraction;
  int    TrMinSiHits;
  double TrMinPt;
  double TrMaxNormChi2;
  std::vector<string> TriggersForMatching_, FiltersForMatching_;
  bool resolveAmbiguity_;
  int  MatchingTriggerResult[50];
  bool   addMuMulessPrimaryVertex_;
  double MuMuMinMass, MuMuMaxMass, JPsiMinMass, JPsiMaxMass;
  double KKMinMass, KKMaxMass, PhiMinMass, PhiMaxMass;
  double JPsiPhiMaxXMass, JPsiPhiMinB0Mass, JPsiPhiMaxB0Mass;
  double MuMuTrackMaxDR, XTrackMaxDR;
  bool   UseXDR ;
  double MuMuKKMinB0Mass, MuMuKKMaxB0Mass, MuMuKKMaxXMass;
  bool   addXlessPrimaryVertex_;
  bool   Debug_;
  std::string DeDxEstimator_, m_dEdxDiscrimTag, m_dEdxDiscrimTag_kaon ;
  TTree* mumukktree;
  unsigned int        runNum, evtNum, lumiNum;
  std::vector<unsigned int>* trigRes;
  std::vector<std::string>*  trigNames;
  std::vector<unsigned int>* L1TT;
  std::vector<std::string>*  MatchTriggerNames;

  /// counters for X
  unsigned int          nMu, nMuMu, nX, nKK;
  unsigned int          nX_pre0, nX_pre1, nX_pre2, nX_pre3, nX_pre4, nX_pre5, nX_pre6, nX_pre7, nX_pre8, nX_pre9, nX_pre10, nX_pre11, nX_pre12, nX_pre13, nX_pre14, nX_pre15;
  int                   n_pV;

  Vertex pV, bS;
  reco::TrackCollection tracks;
  std::vector<Vertex> corrPVs, mumuLessPVs, xLessPvs, xCosAlphaPVs, xCosAlphaXLessPVs, xCosAlpha3DPVs, xCosAlpha3DXLessPVs;
  std::vector<math::XYZTLorentzVector> mumu_p4, muP_p4,muNeg_p4;
  std::vector < CompositeCandidate > ref_Jpsi, ref_mupos, ref_muneg, ref_Phi, ref_kaonpos, ref_kaonneg;
  std::vector < math::XYZTLorentzVector > Jpsi_p4, mupos_p4, muneg_p4, Phi_p4, kpos_p4, kneg_p4;

  // float                 priVtx_X, priVtx_Y, priVtx_Z, priVtx_XE, priVtx_YE, priVtx_ZE, priVtx_NormChi2, priVtx_Chi2, priVtx_CL;
  int                   priVtx_tracks;
  float                 tracksPtSq_pV;
  /// Indices
  std::vector<int>           *mu1Idx, *mu2Idx;
  std::vector<int>           *MuMuType;
  std::vector<int>           *ka1Idx, *ka2Idx;
  std::vector<int>           *X_MuMuIdx, *X_ka1Idx, *X_ka2Idx;
  /// MC Analysis
  // Gen Primary Vertex
  unsigned int          n_genEvtVtx;
  std::vector<float>         *genEvtVtx_X, *genEvtVtx_Y, *genEvtVtx_Z;
  std::vector<int>           *genEvtVtx_particles;
  std::vector<int>           *n_XAncestors;
  unsigned int          nMCAll, nMCX, nMCXVtx;
  std::vector<int>           *MCPdgIdAll, *MCDanNumAll;

  // Gen Primary Vertex

  std::vector<float>         *PriVtxGen_X, *PriVtxGen_Y, *PriVtxGen_Z ;
  std::vector<double>        *PriVtxGen_EX, *PriVtxGen_EY, *PriVtxGen_EZ ;
  std::vector<float>	        *PriVtxGen_Chi2, *PriVtxGen_CL, *PriVtxGen_Ndof;
  std::vector<int>           *PriVtxGen_tracks ;
  std::vector<float>         *MCJPsiPx, *MCJPsiPy, *MCJPsiPz;
  std::vector<float>         *MCmupPx, *MCmupPy, *MCmupPz;
  std::vector<float>         *MCmumPx, *MCmumPy, *MCmumPz;
  std::vector<float>         *MCPhiPx, *MCPhiPy, *MCPhiPz;
  std::vector<float>         *MCkpPx, *MCkpPy, *MCkpPz;
  std::vector<float>         *MCkmPx, *MCkmPy, *MCkmPz;
  //std::vector<float>         *MCpionPx, *MCpionPy, *MCpionPz;
  //std::vector<float>         *MCkaonPx, *MCkaonPy, *MCkaonPz;
  //std::vector<int>           *MCpionCh, *MCkaonCh;
  std::vector<float>         *MCPx, *MCPy, *MCPz;

  /// Generic Muons
  std::vector<float>         *muPx, *muPy, *muPz, *muCharge;
  std::vector<int>           *muPhits, *muShits, *muLayersTr, *muLayersPix;
  std::vector<float>	        *muD0, *muD0E, *muDz, *muChi2 ;
  std::vector<int>           *muNDF;
  std::vector<float>         *mufHits;
  std::vector<bool>          *muFirstBarrel, *muFirstEndCap;
  std::vector<float>	        *muDzVtx, *muDxyVtx, *muDzVtxErr ;
  std::vector<unsigned int>	*muKey;
  std::vector<bool> 	        *muIsGlobal, *muIsPF ;
  std::vector<int>           *muGlMuHits;
  std::vector<float>         *muGlChi2;
  std::vector<int>           *muGlNDF, *muGlMatchedStation;
  std::vector<float>         *muGlDzVtx, *muGlDxyVtx;
  std::vector<int>           *nMatchedStations;
  std::vector<int>           *muType, *muQual, *muTrack, *muNOverlap, *muNSharingSegWith;
  /// Generic tracks
  std::vector<float>         *trNotRef, *trRef;
  std::vector<float>         *trPx, *trPy, *trPz, *trE;
  std::vector<int>           *trNDF, *trPhits, *trShits;
  std::vector<float>         *trChi2;
  std::vector<float>         *trD0, *trD0E, *trCharge;
  std::vector<float>         *trfHits;
  std::vector<bool>          *trFirstBarrel, *trFirstEndCap;
  std::vector<float>         *trDzVtx, *trDxyVtx;
  std::vector<int>           *trQualityHighPurity, *trQualityTight;
  std::vector<double>        *tr_nsigdedx;
  std::vector<float>         *tr_dedx, *tr_dedxMass, *tr_theo, *tr_sigma;
  std::vector<float>         *tr_dedx_byHits, *tr_dedxErr_byHits ;
  std::vector<int>           *tr_saturMeas_byHits, *tr_Meas_byHits ;

  /// MuMu cand & KK cand
  std::vector<float>         *MuMuMass, *MuMuPx, *MuMuPy, *MuMuPz;
  std::vector<float>         *MuMuVtx_CL, *MuMuVtx_Chi2;
  std::vector<float>         *MuMuDecayVtx_X, *MuMuDecayVtx_Y, *MuMuDecayVtx_Z, *MuMuDecayVtx_XE, *MuMuDecayVtx_YE, *MuMuDecayVtx_ZE;
  std::vector<bool>          *MuMuMuonTrigMatch;
  std::vector<float>         *KKMass, *KKPx, *KKPy, *KKPz;
  std::vector<float>         *KKVtx_CL, *KKVtx_Chi2;
  std::vector<float>         *KKDecayVtx_X, *KKDecayVtx_Y, *KKDecayVtx_Z, *KKDecayVtx_XE, *KKDecayVtx_YE, *KKDecayVtx_ZE;
  /// Muons after JPsi (MuMu) fit & Kaons after Phi (KK) fit
  std::vector<float>         *muPos_MuMu_Px, *muPos_MuMu_Py, *muPos_MuMu_Pz ;
  std::vector<float>         *muPos_MuMu_Chi2 ;
  std::vector<int>           *muPos_MuMu_NDF ;
  std::vector<float>         *muNeg_MuMu_Px, *muNeg_MuMu_Py, *muNeg_MuMu_Pz ;
  std::vector<float>         *muNeg_MuMu_Chi2 ;
  std::vector<int>           *muNeg_MuMu_NDF ;
  std::vector<float>         *kaonPos_KK_Px, *kaonPos_KK_Py, *kaonPos_KK_Pz ;
  std::vector<float>         *kaonPos_KK_Chi2 ;
  std::vector<int>           *kaonPos_KK_NDF ;
  std::vector<float>         *kaonNeg_KK_Px, *kaonNeg_KK_Py, *kaonNeg_KK_Pz ;
  std::vector<float>         *kaonNeg_KK_Chi2 ;
  std::vector<int>           *kaonNeg_KK_NDF ;
  std::vector<float>         *DR_MuMu_K1, *DR_MuMu_K2, *DR_MuMuKK_K1, *DR_MuMuKK_K2;
  /// Primary Vertex with "MuMu correction"
  std::vector<int>           *PriVtxMuMuCorr_n, mumuLessPvs_n;
  std::vector<float>         *PriVtxMuMuCorr_X, *PriVtxMuMuCorr_Y, *PriVtxMuMuCorr_Z ;
  std::vector<double>        *PriVtxMuMuCorr_EX, *PriVtxMuMuCorr_EY, *PriVtxMuMuCorr_EZ ;
  std::vector<float>	        *PriVtxMuMuCorr_Chi2, *PriVtxMuMuCorr_CL;
  std::vector<int>           *PriVtxMuMuCorr_tracks ;
  std::vector<int>           *nTrk ;
  /// X candidates
  std::vector<float>         *xMass, *xVtx_CL, *xVtx_Chi2;
  std::vector<float>         *xPx, *xPy, *xPz ;
  std::vector<double>        *xPxE, *xPyE, *xPzE ;
  std::vector<float>         *xDecayVtx_X, *xDecayVtx_Y, *xDecayVtx_Z ;
  std::vector<double>        *xDecayVtx_XE, *xDecayVtx_YE, *xDecayVtx_ZE ;
  /// Muons and tracks after B0 cand fit
  std::vector<float>         *mu1Px_MuMuKK, *mu1Py_MuMuKK, *mu1Pz_MuMuKK, *mu1E_MuMuKK ;
  std::vector<float>         *mu2Px_MuMuKK, *mu2Py_MuMuKK, *mu2Pz_MuMuKK, *mu2E_MuMuKK ;
  std::vector<float>         *k1Px_MuMuKK, *k1Py_MuMuKK, *k1Pz_MuMuKK, *k1E_MuMuKK ;
  std::vector<double>        *kaonPos_nsigdedx;
  std::vector<float>         *kaonPos_dedx, *kaonPos_dedxMass, *kaonPos_theo, *kaonPos_sigma ;
  std::vector<float>         *kaonPos_dedx_byHits, *kaonPos_dedxErr_byHits ;
  std::vector<int>           *kaonPos_saturMeas_byHits, *kaonPos_Meas_byHits ;
  std::vector<float>         *k2Px_MuMuKK, *k2Py_MuMuKK, *k2Pz_MuMuKK, *k2E_MuMuKK ;
  std::vector<double>        *kaonNeg_nsigdedx;
  std::vector<float>         *kaonNeg_dedx, *kaonNeg_dedxMass, *kaonNeg_theo, *kaonNeg_sigma ;
  std::vector<float>         *kaonNeg_dedx_byHits, *kaonNeg_dedxErr_byHits ;
  std::vector<int>           *kaonNeg_saturMeas_byHits, *kaonNeg_Meas_byHits ;
  /// Primary Vertex with largest B0_cos(alpha)
  std::vector<int>           *PriVtx_XCosAlpha_n;
  std::vector<float>         *PriVtx_XCosAlpha_X, *PriVtx_XCosAlpha_Y, *PriVtx_XCosAlpha_Z ;
  std::vector<double>        *PriVtx_XCosAlpha_EX, *PriVtx_XCosAlpha_EY, *PriVtx_XCosAlpha_EZ ;
  std::vector<float>	        *PriVtx_XCosAlpha_Chi2, *PriVtx_XCosAlpha_CL;
  std::vector<int>           *PriVtx_XCosAlpha_tracks ;

  std::vector<int>           *PriVtx_XCosAlpha3D_n;
  std::vector<float>         *PriVtx_XCosAlpha3D_X, *PriVtx_XCosAlpha3D_Y, *PriVtx_XCosAlpha3D_Z ;
  std::vector<double>        *PriVtx_XCosAlpha3D_EX, *PriVtx_XCosAlpha3D_EY, *PriVtx_XCosAlpha3D_EZ ;
  std::vector<float>	        *PriVtx_XCosAlpha3D_Chi2, *PriVtx_XCosAlpha3D_CL;
  std::vector<int>           *PriVtx_XCosAlpha3D_tracks ;
  std::vector<float>         *XLessPV_tracksPtSq, *XLessPV_4tracksPtSq ;
  std::vector<int>           *PriVtxXLess_n;
  std::vector<float>         *PriVtxXLess_X, *PriVtxXLess_Y, *PriVtxXLess_Z ;
  std::vector<double>        *PriVtxXLess_EX, *PriVtxXLess_EY, *PriVtxXLess_EZ ;
  std::vector<float>	        *PriVtxXLess_Chi2, *PriVtxXLess_CL;
  std::vector<int>           *PriVtxXLess_tracks ;
  std::vector<int>           *PriVtxXLess_XCosAlpha_n;
  std::vector<float>         *PriVtxXLess_XCosAlpha_X, *PriVtxXLess_XCosAlpha_Y, *PriVtxXLess_XCosAlpha_Z ;
  std::vector<double>        *PriVtxXLess_XCosAlpha_EX, *PriVtxXLess_XCosAlpha_EY, *PriVtxXLess_XCosAlpha_EZ ;
  std::vector<float>	        *PriVtxXLess_XCosAlpha_Chi2, *PriVtxXLess_XCosAlpha_CL;
  std::vector<int>           *PriVtxXLess_XCosAlpha_tracks ;
  std::vector<int>           *PriVtxXLess_XCosAlpha3D_n;
  std::vector<float>         *PriVtxXLess_XCosAlpha3D_X, *PriVtxXLess_XCosAlpha3D_Y, *PriVtxXLess_XCosAlpha3D_Z ;
  std::vector<double>        *PriVtxXLess_XCosAlpha3D_EX, *PriVtxXLess_XCosAlpha3D_EY, *PriVtxXLess_XCosAlpha3D_EZ ;
  std::vector<float>	        *PriVtxXLess_XCosAlpha3D_Chi2, *PriVtxXLess_XCosAlpha3D_CL;
  std::vector<int>           *PriVtxXLess_XCosAlpha3D_tracks ;
  /// Primary Vertex with "B0 correction"
  std::vector<int>           *PriVtxXCorr_n;
  std::vector<float>         *PriVtxXCorr_X, *PriVtxXCorr_Y, *PriVtxXCorr_Z;
  std::vector<double>        *PriVtxXCorr_EX, *PriVtxXCorr_EY, *PriVtxXCorr_EZ;
  std::vector<float>	        *PriVtxXCorr_Chi2, *PriVtxXCorr_CL;
  std::vector<int>           *PriVtxXCorr_tracks;
  /// Lifetimes variables for B0
  std::vector<double>        *xCosAlphaBS, *xCosAlpha3DBS, *xCTauBS, *xCTauBSE, *xLxyBS, *xLxyBSE, *xLxyzBS, *xLxyzBSE ;
  std::vector<double>        *xCosAlphaPV, *xCosAlpha3DPV, *xCTauPV, *xCTauPVE, *xLxyPV, *xLxyPVE, *xLxyzPV, *xLxyzPVE ;
  std::vector<double>        *xCosAlphaPVCosAlpha, *xCosAlpha3DPVCosAlpha, *xCTauPVCosAlpha, *xCTauPVCosAlphaE, *xLxyPVCosAlpha, *xLxyPVCosAlphaE, *xLxyzPVCosAlpha, *xLxyzPVCosAlphaE ;
  std::vector<double>        *xCosAlphaPVCosAlpha3D, *xCosAlpha3DPVCosAlpha3D, *xCTauPVCosAlpha3D, *xCTauPVCosAlpha3DE, *xLxyPVCosAlpha3D, *xLxyPVCosAlpha3DE, *xLxyzPVCosAlpha3D, *xLxyzPVCosAlpha3DE ;
  std::vector<double>        *xCosAlphaXLessPV, *xCosAlpha3DXLessPV, *xCTauXLessPV, *xCTauXLessPVE, *xLxyXLessPV, *xLxyXLessPVE, *xLxyzXLessPV, *xLxyzXLessPVE ;
  std::vector<double>        *xCosAlphaXLessPVCosAlpha, *xCosAlpha3DXLessPVCosAlpha, *xCTauXLessPVCosAlpha, *xCTauXLessPVCosAlphaE, *xLxyXLessPVCosAlpha, *xLxyXLessPVCosAlphaE, *xLxyzXLessPVCosAlpha, *xLxyzXLessPVCosAlphaE ;
  std::vector<double>        *xCosAlphaXLessPVCosAlpha3D, *xCosAlpha3DXLessPVCosAlpha3D, *xCTauXLessPVCosAlpha3D, *xCTauXLessPVCosAlpha3DE, *xLxyXLessPVCosAlpha3D, *xLxyXLessPVCosAlpha3DE, *xLxyzXLessPVCosAlpha3D, *xLxyzXLessPVCosAlpha3DE ;
  std::vector<double>        *xCosAlphaPVX, *xCTauPVX, *xCTauPVXE, *xLxyPVX, *xLxyPVXE, *xLxyzPVX, *xLxyzPVXE ;
  std::vector<float>	        *xCTauPVX_3D, *xCTauPVX_3D_err;
  /// dxy, dz, dxyE, dzE for kaons from PV, BS, B0LessPV
  std::vector<float>         *kaonPos_dxy_PV, *kaonPos_dz_PV, *kaonNeg_dxy_PV, *kaonNeg_dz_PV;
  std::vector<float>         *kaonPos_dxy_BS, *kaonPos_dz_BS, *kaonNeg_dxy_BS, *kaonNeg_dz_BS;
  std::vector<float>         *kaonPos_dxy_XLessPV, *kaonPos_dz_XLessPV, *kaonNeg_dxy_XLessPV, *kaonNeg_dz_XLessPV;
  std::vector<float>         *kaonPos_dxyE, *kaonPos_dzE, *kaonNeg_dxyE, *kaonNeg_dzE;

  std::vector<bool>          *kaonPosFromPV, *kaonNegFromPV;
};

#endif

// rsync -vut --existing interface/MuMuPiKPAT.h semrat@lxplus.cern.ch:/afs/cern.ch/user/s/semrat/scratch0/CMSSW_5_3_22/src/X4140/MuMuKKPAT/interface/MuMuPiKPAT.h
