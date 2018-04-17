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
  vector<string> TriggersForMatching_, FiltersForMatching_;
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
  TTree* X_One_Tree_;
  unsigned int        runNum, evtNum, lumiNum;
  vector<unsigned int>* trigRes;
  vector<std::string>*  trigNames;
  vector<unsigned int>* L1TT;
  vector<std::string>*  MatchTriggerNames;

  /// counters for X 
  unsigned int          nMu, nMuMu, nX, nKK; 
  unsigned int          nX_pre0, nX_pre1, nX_pre2, nX_pre3, nX_pre4, nX_pre5, nX_pre6, nX_pre7, nX_pre8, nX_pre9, nX_pre10, nX_pre11, nX_pre12, nX_pre13, nX_pre14, nX_pre15; 
  int                   priVtx_n;
  float                 priVtx_X, priVtx_Y, priVtx_Z, priVtx_XE, priVtx_YE, priVtx_ZE, priVtx_NormChi2, priVtx_Chi2, priVtx_CL;
  int                   priVtx_tracks;
  float                 priVtx_tracksPtSq;
  /// Indices
  vector<int>           *mu1Idx, *mu2Idx;
  vector<int>           *MuMuType;
  vector<int>           *ka1Idx, *ka2Idx;
  vector<int>           *X_MuMuIdx, *X_ka1Idx, *X_ka2Idx; 
  /// MC Analysis
  // Gen Primary Vertex
  unsigned int          n_genEvtVtx;
  vector<float>         *genEvtVtx_X, *genEvtVtx_Y, *genEvtVtx_Z; 
  vector<int>           *genEvtVtx_particles;
  vector<int>           *n_XAncestors; 
  unsigned int          nMCAll, nMCX, nMCXVtx; 
  vector<int>           *MCPdgIdAll, *MCDanNumAll;
  // Gen Primary Vertex 
  vector<float>         *PriVtxGen_X, *PriVtxGen_Y, *PriVtxGen_Z ; 
  vector<double>        *PriVtxGen_EX, *PriVtxGen_EY, *PriVtxGen_EZ ;
  vector<float>	        *PriVtxGen_Chi2, *PriVtxGen_CL, *PriVtxGen_Ndof;
  vector<int>           *PriVtxGen_tracks ;
  vector<float>         *MCJPsiPx, *MCJPsiPy, *MCJPsiPz;
  vector<float>         *MCmupPx, *MCmupPy, *MCmupPz;
  vector<float>         *MCmumPx, *MCmumPy, *MCmumPz;
  vector<float>         *MCPhiPx, *MCPhiPy, *MCPhiPz;
  vector<float>         *MCkpPx, *MCkpPy, *MCkpPz;
  vector<float>         *MCkmPx, *MCkmPy, *MCkmPz;
  //vector<float>         *MCpionPx, *MCpionPy, *MCpionPz;
  //vector<float>         *MCkaonPx, *MCkaonPy, *MCkaonPz;
  //vector<int>           *MCpionCh, *MCkaonCh;
  vector<float>         *MCPx, *MCPy, *MCPz;
  /// Generic Muons
  vector<float>         *muPx, *muPy, *muPz, *muCharge;
  vector<int>           *muPhits, *muShits, *muLayersTr, *muLayersPix;
  vector<float>	        *muD0, *muD0E, *muDz, *muChi2 ;
  vector<int>           *muNDF;
  vector<float>         *mufHits;
  vector<bool>          *muFirstBarrel, *muFirstEndCap;
  vector<float>	        *muDzVtx, *muDxyVtx, *muDzVtxErr ;   
  vector<unsigned int>	*muKey;
  vector<bool> 	        *muIsGlobal, *muIsPF ;
  vector<int>           *muGlMuHits;
  vector<float>         *muGlChi2;
  vector<int>           *muGlNDF, *muGlMatchedStation;
  vector<float>         *muGlDzVtx, *muGlDxyVtx;
  vector<int>           *nMatchedStations;
  vector<int>           *muType, *muQual, *muTrack, *muNOverlap, *muNSharingSegWith;
  /// Generic tracks
  vector<float>         *trNotRef, *trRef;
  vector<float>         *trPx, *trPy, *trPz, *trE;
  vector<int>           *trNDF, *trPhits, *trShits;
  vector<float>         *trChi2;
  vector<float>         *trD0, *trD0E, *trCharge;
  vector<float>         *trfHits;
  vector<bool>          *trFirstBarrel, *trFirstEndCap;
  vector<float>         *trDzVtx, *trDxyVtx;
  vector<int>           *trQualityHighPurity, *trQualityTight;
  vector<double>        *tr_nsigdedx;
  vector<float>         *tr_dedx, *tr_dedxMass, *tr_theo, *tr_sigma;
  vector<float>         *tr_dedx_byHits, *tr_dedxErr_byHits ;
  vector<int>           *tr_saturMeas_byHits, *tr_Meas_byHits ;
  /// MuMu cand & KK cand 
  vector<float>         *MuMuMass, *MuMuPx, *MuMuPy, *MuMuPz;
  vector<float>         *MuMuVtx_CL, *MuMuVtx_Chi2;
  vector<float>         *MuMuDecayVtx_X, *MuMuDecayVtx_Y, *MuMuDecayVtx_Z, *MuMuDecayVtx_XE, *MuMuDecayVtx_YE, *MuMuDecayVtx_ZE;
  vector<bool>          *MuMuMuonTrigMatch;
  vector<float>         *KKMass, *KKPx, *KKPy, *KKPz;
  vector<float>         *KKVtx_CL, *KKVtx_Chi2;
  vector<float>         *KKDecayVtx_X, *KKDecayVtx_Y, *KKDecayVtx_Z, *KKDecayVtx_XE, *KKDecayVtx_YE, *KKDecayVtx_ZE;
  /// Muons after JPsi (MuMu) fit & Kaons after Phi (KK) fit
  vector<float>         *mu1_MuMu_Px, *mu1_MuMu_Py, *mu1_MuMu_Pz ;
  vector<float>         *mu1_MuMu_Chi2 ;
  vector<int>           *mu1_MuMu_NDF ;
  vector<float>         *mu2_MuMu_Px, *mu2_MuMu_Py, *mu2_MuMu_Pz ;
  vector<float>         *mu2_MuMu_Chi2 ;
  vector<int>           *mu2_MuMu_NDF ;
  vector<float>         *ka1_KK_Px, *ka1_KK_Py, *ka1_KK_Pz ;
  vector<float>         *ka1_KK_Chi2 ;
  vector<int>           *ka1_KK_NDF ;
  vector<float>         *ka2_KK_Px, *ka2_KK_Py, *ka2_KK_Pz ;
  vector<float>         *ka2_KK_Chi2 ;
  vector<int>           *ka2_KK_NDF ;
  vector<float>         *DR_MuMu_K1, *DR_MuMu_K2, *DR_MuMuKK_K1, *DR_MuMuKK_K2;
  /// Primary Vertex with "MuMu correction"
  vector<int>           *PriVtxMuMuCorr_n;
  vector<float>         *PriVtxMuMuCorr_X, *PriVtxMuMuCorr_Y, *PriVtxMuMuCorr_Z ; 
  vector<double>        *PriVtxMuMuCorr_EX, *PriVtxMuMuCorr_EY, *PriVtxMuMuCorr_EZ ;
  vector<float>	        *PriVtxMuMuCorr_Chi2, *PriVtxMuMuCorr_CL;
  vector<int>           *PriVtxMuMuCorr_tracks ;
  vector<int>           *nTrk ;
  /// X candidates 
  vector<float>         *xMass, *xVtx_CL, *xVtx_Chi2; 
  vector<float>         *xPx, *xPy, *xPz ; 
  vector<double>        *xPxE, *xPyE, *xPzE ; 
  vector<float>         *xDecayVtx_X, *xDecayVtx_Y, *xDecayVtx_Z ; 
  vector<double>        *xDecayVtx_XE, *xDecayVtx_YE, *xDecayVtx_ZE ;
  /// Muons and tracks after B0 cand fit  
  vector<float>         *mu1Px_MuMuKK, *mu1Py_MuMuKK, *mu1Pz_MuMuKK, *mu1E_MuMuKK ;
  vector<float>         *mu2Px_MuMuKK, *mu2Py_MuMuKK, *mu2Pz_MuMuKK, *mu2E_MuMuKK ;
  vector<float>         *k1Px_MuMuKK, *k1Py_MuMuKK, *k1Pz_MuMuKK, *k1E_MuMuKK ;  
  vector<double>        *kaon1_nsigdedx; 
  vector<float>         *kaon1_dedx, *kaon1_dedxMass, *kaon1_theo, *kaon1_sigma ;
  vector<float>         *kaon1_dedx_byHits, *kaon1_dedxErr_byHits ; 
  vector<int>           *kaon1_saturMeas_byHits, *kaon1_Meas_byHits ; 
  vector<float>         *k2Px_MuMuKK, *k2Py_MuMuKK, *k2Pz_MuMuKK, *k2E_MuMuKK ;
  vector<double>        *kaon2_nsigdedx; 
  vector<float>         *kaon2_dedx, *kaon2_dedxMass, *kaon2_theo, *kaon2_sigma ; 
  vector<float>         *kaon2_dedx_byHits, *kaon2_dedxErr_byHits ; 
  vector<int>           *kaon2_saturMeas_byHits, *kaon2_Meas_byHits ; 
  /// Primary Vertex with largest B0_cos(alpha) 
  vector<int>           *PriVtx_XCosAlpha_n; 
  vector<float>         *PriVtx_XCosAlpha_X, *PriVtx_XCosAlpha_Y, *PriVtx_XCosAlpha_Z ; 
  vector<double>        *PriVtx_XCosAlpha_EX, *PriVtx_XCosAlpha_EY, *PriVtx_XCosAlpha_EZ ; 
  vector<float>	        *PriVtx_XCosAlpha_Chi2, *PriVtx_XCosAlpha_CL; 
  vector<int>           *PriVtx_XCosAlpha_tracks ; 

  vector<int>           *PriVtx_XCosAlpha3D_n; 
  vector<float>         *PriVtx_XCosAlpha3D_X, *PriVtx_XCosAlpha3D_Y, *PriVtx_XCosAlpha3D_Z ; 
  vector<double>        *PriVtx_XCosAlpha3D_EX, *PriVtx_XCosAlpha3D_EY, *PriVtx_XCosAlpha3D_EZ ; 
  vector<float>	        *PriVtx_XCosAlpha3D_Chi2, *PriVtx_XCosAlpha3D_CL;
  vector<int>           *PriVtx_XCosAlpha3D_tracks ;
  vector<float>         *XLessPV_tracksPtSq, *XLessPV_4tracksPtSq ;
  vector<int>           *PriVtxXLess_n;
  vector<float>         *PriVtxXLess_X, *PriVtxXLess_Y, *PriVtxXLess_Z ; 
  vector<double>        *PriVtxXLess_EX, *PriVtxXLess_EY, *PriVtxXLess_EZ ;
  vector<float>	        *PriVtxXLess_Chi2, *PriVtxXLess_CL;
  vector<int>           *PriVtxXLess_tracks ;
  vector<int>           *PriVtxXLess_XCosAlpha_n;
  vector<float>         *PriVtxXLess_XCosAlpha_X, *PriVtxXLess_XCosAlpha_Y, *PriVtxXLess_XCosAlpha_Z ; 
  vector<double>        *PriVtxXLess_XCosAlpha_EX, *PriVtxXLess_XCosAlpha_EY, *PriVtxXLess_XCosAlpha_EZ ;
  vector<float>	        *PriVtxXLess_XCosAlpha_Chi2, *PriVtxXLess_XCosAlpha_CL;
  vector<int>           *PriVtxXLess_XCosAlpha_tracks ;
  vector<int>           *PriVtxXLess_XCosAlpha3D_n;
  vector<float>         *PriVtxXLess_XCosAlpha3D_X, *PriVtxXLess_XCosAlpha3D_Y, *PriVtxXLess_XCosAlpha3D_Z ; 
  vector<double>        *PriVtxXLess_XCosAlpha3D_EX, *PriVtxXLess_XCosAlpha3D_EY, *PriVtxXLess_XCosAlpha3D_EZ ;
  vector<float>	        *PriVtxXLess_XCosAlpha3D_Chi2, *PriVtxXLess_XCosAlpha3D_CL;
  vector<int>           *PriVtxXLess_XCosAlpha3D_tracks ;
  /// Primary Vertex with "B0 correction"
  vector<int>           *PriVtxXCorr_n; 
  vector<float>         *PriVtxXCorr_X, *PriVtxXCorr_Y, *PriVtxXCorr_Z; 
  vector<double>        *PriVtxXCorr_EX, *PriVtxXCorr_EY, *PriVtxXCorr_EZ; 
  vector<float>	        *PriVtxXCorr_Chi2, *PriVtxXCorr_CL; 
  vector<int>           *PriVtxXCorr_tracks; 
  /// Lifetimes variables for B0 
  vector<double>        *xCosAlphaBS, *xCosAlpha3DBS, *xCTauBS, *xCTauBSE, *xLxyBS, *xLxyBSE, *xLxyzBS, *xLxyzBSE ;
  vector<double>        *xCosAlphaPV, *xCosAlpha3DPV, *xCTauPV, *xCTauPVE, *xLxyPV, *xLxyPVE, *xLxyzPV, *xLxyzPVE ;
  vector<double>        *xCosAlphaPVCosAlpha, *xCosAlpha3DPVCosAlpha, *xCTauPVCosAlpha, *xCTauPVCosAlphaE, *xLxyPVCosAlpha, *xLxyPVCosAlphaE, *xLxyzPVCosAlpha, *xLxyzPVCosAlphaE ;
  vector<double>        *xCosAlphaPVCosAlpha3D, *xCosAlpha3DPVCosAlpha3D, *xCTauPVCosAlpha3D, *xCTauPVCosAlpha3DE, *xLxyPVCosAlpha3D, *xLxyPVCosAlpha3DE, *xLxyzPVCosAlpha3D, *xLxyzPVCosAlpha3DE ;
  vector<double>        *xCosAlphaXLessPV, *xCosAlpha3DXLessPV, *xCTauXLessPV, *xCTauXLessPVE, *xLxyXLessPV, *xLxyXLessPVE, *xLxyzXLessPV, *xLxyzXLessPVE ;
  vector<double>        *xCosAlphaXLessPVCosAlpha, *xCosAlpha3DXLessPVCosAlpha, *xCTauXLessPVCosAlpha, *xCTauXLessPVCosAlphaE, *xLxyXLessPVCosAlpha, *xLxyXLessPVCosAlphaE, *xLxyzXLessPVCosAlpha, *xLxyzXLessPVCosAlphaE ;
  vector<double>        *xCosAlphaXLessPVCosAlpha3D, *xCosAlpha3DXLessPVCosAlpha3D, *xCTauXLessPVCosAlpha3D, *xCTauXLessPVCosAlpha3DE, *xLxyXLessPVCosAlpha3D, *xLxyXLessPVCosAlpha3DE, *xLxyzXLessPVCosAlpha3D, *xLxyzXLessPVCosAlpha3DE ;
  vector<double>        *xCosAlphaPVX, *xCTauPVX, *xCTauPVXE, *xLxyPVX, *xLxyPVXE, *xLxyzPVX, *xLxyzPVXE ;
  vector<float>	        *xCTauPVX_3D, *xCTauPVX_3D_err;
  /// dxy, dz, dxyE, dzE for kaons from PV, BS, B0LessPV
  vector<float>         *kaon1_dxy_PV, *kaon1_dz_PV, *kaon2_dxy_PV, *kaon2_dz_PV;
  vector<float>         *kaon1_dxy_BS, *kaon1_dz_BS, *kaon2_dxy_BS, *kaon2_dz_BS;
  vector<float>         *kaon1_dxy_XLessPV, *kaon1_dz_XLessPV, *kaon2_dxy_XLessPV, *kaon2_dz_XLessPV;
  vector<float>         *kaon1_dxyE, *kaon1_dzE, *kaon2_dxyE, *kaon2_dzE; 

  vector<bool>          *Kaon1FromPV, *Kaon2FromPV;
};

#endif

// rsync -vut --existing interface/MuMuPiKPAT.h semrat@lxplus.cern.ch:/afs/cern.ch/user/s/semrat/scratch0/CMSSW_5_3_22/src/X4140/MuMuKKPAT/interface/MuMuPiKPAT.h
