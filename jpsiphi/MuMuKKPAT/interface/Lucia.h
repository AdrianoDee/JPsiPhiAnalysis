// -*- C++ -*-
//
// Package:    JPsiPhiKPAT
// Class:      JPsiPhiKPAT
// 
/**\class JPsiPhiKPAT JPsiPhiKPAT.cc myAnalyzers/JPsiPhiKPAT/src/JPsiPhiKPAT.cc

 Description: <one line class summary>
Make rootTuple for JPsiPiPi reconstruction

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Keith Ulmer   keith.ulmer@colorado.edu
//
//

#ifndef _JPsiPhiKPAT_h
#define _JPsiPhiKPAT_h

// system include files
#include <memory>

// user include files
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
//
// class decleration
//

using std::vector;
using namespace edm;
using namespace reco;
using namespace std;

class JPsiPhiKPAT : public edm::EDAnalyzer {
public:
  explicit JPsiPhiKPAT(const edm::ParameterSet&);
  ~JPsiPhiKPAT();
  
private:
  virtual void beginJob() ;
  virtual void beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
	
  const reco::DeDxDataValueMap * energyLoss;
  Int_t iexception_dedx;
  
  virtual double getSigmaOfLogdEdx(double logde);
  virtual float  getEnergyLoss(const reco::TrackRef & track);
  virtual double nsigmaofdedx(const reco::TrackRef & track, double & theo, double & sigma);
  virtual double getLogdEdx(double bg);
  virtual double GetMass(const reco::TrackRef & track);
  
  // ----------member data ---------------------------
  std::string proccessName_;
  HLTConfigProvider hltConfig_;

  edm::InputTag hlTriggerResults_;
  edm::InputTag inputGEN_;
  std::string vtxSample;
  bool doData;
  bool doMC;
  int  MCParticle;

  bool doJPsiMassConst;

  int MuMinPixHits;
  int MuMinSiHits;
  double MuMaxNormChi;
  double MuMaxD0;
  double MuMinPt;

  double JPsiMinMass;
  double JPsiMaxMass;

  double PhiMinMass;
  double PhiMaxMass;

  int    TrMinPixHits;
  int    TrMinSiHits;
  double TrMinPt;
  double TrMaxNormChi2;

  double JPsiTrackMaxDR;
  double BchTrackMaxDR;
  bool   UseBchDR;

  double JPsiPhiKMinMass;
  double JPsiPhiKMaxMass;
	
  bool resolveAmbiguity_; 
  bool addBchLessPrimaryVertex_;
  vector<string>      TriggersForMatching_;
  int  MatchingTriggerResult[50];
  bool Debug_;

  TTree* Y_One_Tree_;

  unsigned int        runNum, evtNum, lumiNum;

  unsigned int        doubleMuonsEvents;
  //
  vector<unsigned int>* trigRes;
  vector<std::string>*  trigNames;
  //
  vector<unsigned int>* L1TT;
  vector<std::string>*  MatchTriggerNames;

  //// MC Analysis
  unsigned int        nMCB; 
  vector<float>       *MCjpsiPx, *MCjpsiPy, *MCjpsiPz;
  vector<float>       *MCphiPx, *MCphiPy, *MCphiPz;
  vector<float>       *MCkaonPx, *MCkaonPy, *MCkaonPz;
  ////////////////

  unsigned int        nB, nJPsi, nPhi, nMu, nTrk;
  float               priVtxX, priVtxY, priVtxZ, priVtxXE, priVtxYE, priVtxZE, priVtxNormChi2, priVtxChi2, priVtxCL;

  //// Indices
  vector<int>         *jpsiIdx, *phiIdx;
  vector<int>         *k1Idx, *k2Idx, *k3Idx;
  vector<int>         *mu1Idx, *mu2Idx;

  //// Generic Muons
  vector<float>       *muPx, *muPy, *muPz, *muChi2;
  vector<int>         *muNDF, *muCharge;
  vector<float>       *muD0, *muDz, *muDzVtx, *muDxyVtx;
  vector<int>         *muType, *muQual;
  vector<float>       *muGlChi2; 
  vector<int>         *muGlNDF;
  vector<float>       *mufHits; 
  vector<bool>        *muFirstBarrel, *muFirstEndCap;
	
  vector<int>         *muPhits, *muShits, *muGlMuonHits;


  //// Generic tracks
  vector<float>       *trPx, *trPy, *trPz, *trE;
  
  vector<int>         *trNDF, *trPhits, *trShits;
  vector<float>       *trChi2;
  vector<float>       *trfHits;
  vector<bool>        *trFirstBarrel, *trFirstEndCap;
  vector<float>       *trDzVtx, *trDxyVtx;
  vector<float>       *trD0, *trD0E, *trDz;
  vector<int>         *trCharge;

  vector<double>      *tr_nsigdedx;
  vector<float>       *tr_dedx, *tr_dedxMass, *tr_theo, *tr_sigma;

  //// J/Psi
  vector<float>       *jpsiMass, *jpsiPx, *jpsiPy, *jpsiPz;
  vector<float>       *jpsiVtxCL, *jpsiVtxChi2;
  vector<float>       *jpsiDecayVtxX, *jpsiDecayVtxY, *jpsiDecayVtxZ;
  vector<float>       *jpsiDecayVtxXE, *jpsiDecayVtxYE, *jpsiDecayVtxZE;
  vector<float>       *jpsiCTauPV, *jpsiCTauPVE, *jpsiFL, *jpsiFLE;
  vector<bool>        *JPsiMuonTrigMatch;

  //// Muons after J/Psi fit
  vector<float>       *mu1Px, *mu1Py, *mu1Pz, *mu1Chi2;
  vector<int>         *mu1NDF;
  vector<float>       *mu2Px, *mu2Py, *mu2Pz, *mu2Chi2;
  vector<int>         *mu2NDF;

  //// Phi
  vector<float>       *phiMass, *phiPx, *phiPy, *phiPz;
  vector<float>       *phiVtxCL, *phiVtxChi2;
  vector<float>       *phiDecayVtxX, *phiDecayVtxY, *phiDecayVtxZ;
  vector<float>       *phiDecayVtxXE, *phiDecayVtxYE, *phiDecayVtxZE;

  //// Kaons after Phi fit
  vector<float>       *k1Px, *k1Py, *k1Pz, *k1Chi2;
  vector<int>         *k1NDF;
  vector<float>       *k2Px, *k2Py, *k2Pz, *k2Chi2;
  vector<int>         *k2NDF;

  //// B charged Cand 
  vector<float>       *bMass, *bVtxCL, *bVtxChi2;
  vector<int>         *bCharge;
  vector<float>       *bPx, *bPy, *bPz;
  vector<double>      *bPxE, *bPyE, *bPzE;
  vector<float>       *bDecayVtxX, *bDecayVtxY, *bDecayVtxZ;
  vector<double>      *bDecayVtxXE, *bDecayVtxYE, *bDecayVtxZE;

  vector<float>       *PriVtxBCorrX, *PriVtxBCorrY, *PriVtxBCorrZ;
  vector<double>      *PriVtxBCorrXE, *PriVtxBCorrYE, *PriVtxBCorrZE;
  vector<float>	      *PriVtxBCorrChi2, *PriVtxBCorrCL;

  vector<double>      *bLxyPV, *bCosAlphaPV, *bCTauPV, *bLxyBS, *bCosAlphaBS, *bCTauBS;
  vector<double>      *bCTauPVE, *bCTauBSE, *bLxyPVE, *bLxyBSE;

  //// Muons and tracks after Bch Cand fit
  vector<float>       *k1fitPx, *k1fitPy, *k1fitPz, *k1fitE;
  vector<float>       *k2fitPx, *k2fitPy, *k2fitPz, *k2fitE;
  vector<float>       *k3fitPx, *k3fitPy, *k3fitPz, *k3fitE;
  vector<float>       *mu1fitPx, *mu1fitPy, *mu1fitPz, *mu1fitE;
  vector<float>       *mu2fitPx, *mu2fitPy, *mu2fitPz, *mu2fitE;

};

#endif
