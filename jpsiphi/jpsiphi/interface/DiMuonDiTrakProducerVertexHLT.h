/**
   \file
   Declaration of DiMuonDiTrakProducerVertexHLT

   \author Alberto Sanchez-Hernandez
   \date 2 Mar 2014
*/

#ifndef __DiMuonDiTrakProducerVertexHLT_h_
#define __DiMuonDiTrakProducerVertexHLT_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include <TLorentzVector.h>
#include <vector>

/**
   Create a HF candidate by mathing DiMuon(chi,psi,etc.) and a track (K, pi, etc.)
 */

class DiMuonDiTrakProducerVertexHLT : public edm::EDProducer {

 public:
  explicit DiMuonDiTrakProducerVertexHLT(const edm::ParameterSet& ps);

 private:

  void produce(edm::Event& event, const edm::EventSetup& esetup) override;

  void endJob() override;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> DiMuonCollection_;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate>> TrakCollection_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> TriggerCollection_;

  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;

  std::vector<double> DiMuonMassCuts_;
  std::vector<double> TrakTrakMassCuts_;
  std::vector<double> DiMuonDiTrakMassCuts_;
  std::vector<double> MassTraks_;
  std::vector<double> MaxDeltaRPt_;

  bool OnlyBest_;
  std::string product_name_;

  std::vector<std::string>  HLTs_;
  std::vector<std::string>  HLTFilters_;

  reco::Candidate::LorentzVector convertVector(const math::XYZTLorentzVectorF& v);
  bool IsTheSame(const pat::PackedCandidate& tk, const pat::Muon& mu);
  const pat::CompositeCandidate makeDiMuonTTCandidate(const pat::CompositeCandidate& DiMuon,
						    const pat::CompositeCandidate& tt);

  const pat::CompositeCandidate makeTTCandidate(const pat::PackedCandidate& trak1,
                                                const pat::PackedCandidate& trak2);
  const pat::CompositeCandidate makeTTTriggerCandidate(const pat::TriggerObjectStandAlone& t1,
						    const pat::TriggerObjectStandAlone& t2);
  const pat::CompositeCandidate makeTTTriggerMixedCandidate(
                                                          const pat::PackedCandidate& trakP,
                                                          const pat::TriggerObjectStandAlone& trakN)

  bool MatchByDRDPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2);
  float DeltaR(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2);

  UInt_t isTriggerMatched(const pat::TriggerObjectStandAlone& t);

  int candidates;
  int nevents;
  int ndimuon;
  int nreco;


};

#endif // __DiMuonDiTrakProducerVertexHLT_h_
