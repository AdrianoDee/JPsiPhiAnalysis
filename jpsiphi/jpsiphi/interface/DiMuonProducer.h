#ifndef JpsiPhiAnalysis_DiMuonProducer_DiMuonProducerPAT_h
#define JpsiPhiAnalysis_DiMuonProducer_DiMuonProducerPAT_h


// system include files
#include <memory>

// FW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/Utils/interface/PtComparator.h"

// DataFormat includes
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include <CommonTools/UtilAlgos/interface/StringCutObjectSelector.h>
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"

template<typename T>
struct GreaterByVProb {
  bool operator()( const T & t1, const T & t2 ) const {
    return t1.userFloat("vProb") > t2.userFloat("vProb");
  }
};


//
// class decleration
//

class DiMuonProducerPAT : public edm::EDProducer {
 public:
  explicit DiMuonProducerPAT(const edm::ParameterSet&);
  ~DiMuonProducerPAT() override;

 private:
  void beginJob() override ;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endJob() override ;
  bool isAbHadron(int pdgID);
  bool isAMixedbHadron(int pdgID, int momPdgID);
  UInt_t isTriggerMatched(pat::CompositeCandidate *diMuon_cand);
  const pat::TriggerObjectStandAlone BestTriggerMuon(const pat::Muon& m);
  bool MatchByDRDPt(const pat::Muon t1, const pat::TriggerObjectStandAlone t2);
  UInt_t isTriggerMatched(const pat::Muon& m);
  float DeltaR(const pat::Muon t1, const pat::TriggerObjectStandAlone t2);
  float DeltaPt(const pat::Muon t1, const pat::TriggerObjectStandAlone t2);
  const pat::CompositeCandidate makeMuMuTriggerCand(const pat::TriggerObjectStandAlone& muonP, const pat::TriggerObjectStandAlone& muonN);
  std::pair<int, float> findJpsiMCInfo(reco::GenParticleRef genJpsi);

  // ----------member data ---------------------------
 private:

  edm::EDGetTokenT<edm::View<pat::Muon>> muons_;
  double muonPtCut_;
  edm::EDGetTokenT<reco::BeamSpot> thebeamspot_;
  edm::EDGetTokenT<reco::VertexCollection> thePVs_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> TriggerCollection_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  StringCutObjectSelector<reco::Candidate, true> dimuonSelection_;
  bool addCommonVertex_, addMuonlessPrimaryVertex_;
  bool resolveAmbiguity_;
  bool addMCTruth_;
  GreaterByVProb<pat::CompositeCandidate> vPComparator_;
  std::vector<std::string> HLTFilters_;

  InvariantMassFromVertex massCalculator;

  float maxDeltaR, maxDPtRel;

};

#endif
