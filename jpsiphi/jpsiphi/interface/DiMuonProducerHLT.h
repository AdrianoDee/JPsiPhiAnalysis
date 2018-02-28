#ifndef JpsiPhiAnalysis_DiMuonProducer_DiMuonProducerHLTPAT_h
#define JpsiPhiAnalysis_DiMuonProducer_DiMuonProducerHLTPAT_h


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

class DiMuonProducerHLTPAT : public edm::EDProducer {
 public:
  explicit DiMuonProducerHLTPAT(const edm::ParameterSet&);
  ~DiMuonProducerHLTPAT() override;

 private:
  void beginJob() override ;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endJob() override ;
  bool isAbHadron(int pdgID);
  bool isAMixedbHadron(int pdgID, int momPdgID);

  UInt_t isTriggerMatched(pat::CompositeCandidate *diMuon_cand);
  bool isTriggerMatched(pat::Muon *m);
  UInt_t isTriggerMatched(const pat::TriggerObjectStandAlone *t);

  float DeltaR(const pat::Muon m, const pat::TriggerObjectStandAlone t);

  const pat::CompositeCandidate makeMuMuTriggerCand(const pat::TriggerObjectStandAlone& muonP, const pat::TriggerObjectStandAlone& muonN);

  const pat::TriggerObjectStandAlone& BestTriggerMuon(const pat::Muon* m);
  std::pair<int, float> findJpsiMCInfo(reco::GenParticleRef genJpsi);

  // ----------member data ---------------------------
 private:

  edm::EDGetTokenT<edm::View<pat::Muon>> muons_;
  edm::EDGetTokenT<reco::BeamSpot> thebeamspot_;
  edm::EDGetTokenT<reco::VertexCollection> thePVs_;
  edm::EDGetTokenT<reco::TrackCollection> revtxtrks_;
  edm::EDGetTokenT<reco::BeamSpot> revtxbs_;
  StringCutObjectSelector<pat::Muon> higherPuritySelection_;
  StringCutObjectSelector<pat::Muon> lowerPuritySelection_;
  StringCutObjectSelector<reco::Candidate, true> dimuonSelection_;
  bool addCommonVertex_, addMuonlessPrimaryVertex_;
  bool resolveAmbiguity_;
  bool addMCTruth_;
  GreaterByVProb<pat::CompositeCandidate> vPComparator_;
  std::vector<std::string> HLTFilters_;

  InvariantMassFromVertex massCalculator;

};

#endif
