/**
   \file
   Declaration of DiTrackHLTProducer

   \author Alberto Sanchez-Hernandez
   \date 2 Mar 2014
*/

#ifndef __DiTrackHLTProducer_h_
#define __DiTrackHLTProducer_h_

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

#include "CommonTools/UtilAlgos/interface/DeltaR.h"
#include "CommonTools/UtilAlgos/interface/MatchByDRDPt.h"
#include <TLorentzVector.h>
#include <vector>
#include <cmath>

/**
   Create a HF candidate by mathing DiMuon(chi,psi,etc.) and a track (K, pi, etc.)
 */

class DiTrackHLTProducer : public edm::EDProducer {

 public:
  explicit DiTrackHLTProducer(const edm::ParameterSet& ps);

 private:

  void produce(edm::Event& event, const edm::EventSetup& esetup) override;

  void endJob() override;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate>> TrakCollection_;
  std::vector<double> TrakTrakMassCuts_;
  std::vector<double> MassTraks_;
  bool OnlyBest_;
  std::string product_name_;
  std::string HLTFilters_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObj_;

  reco::Candidate::LorentzVector convertVector(const math::XYZTLorentzVectorF& v);
  bool IsTheSame(const pat::PackedCandidate& tk, const pat::Muon& mu);
  const pat::CompositeCandidate makeDiMuonTTCandidate(const pat::CompositeCandidate& DiMuon,
						    const pat::CompositeCandidate& tt);
  const pat::CompositeCandidate makeTTCandidate(const pat::PackedCandidate& trak1,
                                                const pat::PackedCandidate& trak2);

  bool MatchByDRDPt(pat::CompositeCandidate t1, pat::TriggerObjectStandAlone t2);
  float DeltaR(pat::CompositeCandidate t1, pat::TriggerObjectStandAlone t2);

  int candidates;
  int nevents;
  int ndimuon;
  int nreco;
  float maxDeltaR;
  float maxDPtRel;
};

#endif // __DiTrackHLTProducer_h_
