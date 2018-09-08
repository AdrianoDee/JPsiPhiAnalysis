/**
   \file
   Declaration of DiMuonDiTrakProducerAOD

   \author Alberto Sanchez-Hernandez
   \date 2 Mar 2014
*/

#ifndef __DiMuonDiTrakProducerAOD_h_
#define __DiMuonDiTrakProducerAOD_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "FWCore/Common/interface/TriggerNames.h"

#include <TLorentzVector.h>
#include <vector>

/**
   Create a HF candidate by mathing DiMuon(chi,psi,etc.) and a track (K, pi, etc.)
 */

class DiMuonDiTrakProducerAOD : public edm::EDProducer {

 public:
  explicit DiMuonDiTrakProducerAOD(const edm::ParameterSet& ps);

 private:

  void produce(edm::Event& event, const edm::EventSetup& esetup) override;

  void endJob() override;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> DiMuonCollection_;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate>> TrakCollection_;
  edm::EDGetTokenT<trigger::TriggerEvent > TriggerEvent_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  std::vector<double> DiMuonMassCuts_;
  std::vector<double> TrakTrakMassCuts_;
  std::vector<double> DiMuonDiTrakMassCuts_;
  std::vector<double> MassTraks_;
  bool OnlyBest_;
  std::string product_name_;
  std::vector<std::string>  HLTFilters_;
  bool isMC_;

  reco::Candidate::LorentzVector convertVector(const math::XYZTLorentzVectorF& v);
  bool IsTheSame(const pat::PackedCandidate& tk, const pat::Muon& mu);
  pat::CompositeCandidate makeDiMuonTTCandidate(const pat::CompositeCandidate& DiMuon,
						    const pat::CompositeCandidate& tt);
  pat::CompositeCandidate makeTTCandidate(const pat::PackedCandidate& trak1,
                                                const pat::PackedCandidate& trak2);

  bool MatchByDRDPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2);
  float DeltaR(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2);

  float maxDeltaR;
  float maxDPtRel;

  int candidates;
  int nevents;
  int ndimuon;
  int nreco;
};

#endif // __DiMuonDiTrakProducerAOD_h_
