// -*- C++ -*-
//
// Package:    MuMuProducerPAT
// Class:      MuMuProducerPAT
//
/**\class MuMuProducerPAT MuMuProducerPAT.cc myProducers/MuMuProducerPAT/src/MuMuProducerPAT.cc

   Description: <one line class summary>
   Make rootTuple for JPsiKK reconstruction

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Keith Ulmer   keith.ulmer@colorado.edu
//
//

#ifndef _MuMuProducerPAT_h
#define _MuMuProducerPAT_h

// system include files
#include <memory>

/// user include files
#include "../interface/VertexReProducer.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
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
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
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

template<typename T>
struct GreaterByVProb {
  bool operator()( const T & t1, const T & t2 ) const {
    return t1.userFloat("vProb") > t2.userFloat("vProb");
  }
};

class MuMuProducerPAT : public edm::EDProducer {
public:
  explicit MuMuProducerPAT(const edm::ParameterSet&);
  ~MuMuProducerPAT();

private:
  virtual void beginJob() ;
  virtual void beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup);
  virtual void produce(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  UInt_t isTriggerMatched(pat::CompositeCandidate *diMuon_cand);

  InvariantMassFromVertex massCalculator;


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
  std::string vtxSample_;
  std::vector<float> jspiMassCuts_,psiMassCuts_;
  bool doData_, doMC_,addCommonVertex_,resolveAmbiguity_,addMCTruth_;
  int  MCParticle_;
  bool MCExclusiveDecay_;
  int  MCMother_, MCDaughtersN_;

  int MuMinPixHits_, MuMinSiHits_;
  double MuMaxNormChi_;
  double MuMaxD0_;
  bool resolveAmbiguity_,TriggerCut_;

  std::vector<string> HLTFileters_, FiltersForMatching_;
  bool Debug_;
  
};

#endif

// rsync -vut --existing interface/MuMuPiKPAT.h semrat@lxplus.cern.ch:/afs/cern.ch/user/s/semrat/scratch0/CMSSW_5_3_22/src/X4140/MuMuProducerPAT/interface/MuMuPiKPAT.h
