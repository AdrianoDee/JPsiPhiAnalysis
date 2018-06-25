#ifndef __DoubleDiMuonProducer_h_
#define __DoubleDiMuonProducer_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include <TLorentzVector.h>
#include <vector>

/**
   Create a HF candidate by mathing Onia(chi,psi,etc.) and a track (K, pi, etc.)
 */

class DoubleDiMuonProducer : public edm::EDProducer {

 public:
  explicit DoubleDiMuonProducer(const edm::ParameterSet& ps);

 private:

  void produce(edm::Event& event, const edm::EventSetup& esetup) override;

  void endJob() override;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> HighDiMuonCollection_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> LowDiMuonCollection_;
  std::vector<double> HighDiMuonMassCuts_;
  std::vector<double> LowDiMuonMassCuts_;
  std::vector<double> DoubleDiMuonMassCuts_;
  bool addMCTruth_;

  reco::Candidate::LorentzVector convertVector(const math::XYZTLorentzVectorF& v);
  std::tuple<int, float, float> findJpsiMCInfo(reco::GenParticleRef genJpsi);
  const pat::CompositeCandidate makeCandidate(const pat::CompositeCandidate& l,
    const pat::CompositeCandidate& h);

  int candidates;
  int nevents;
  int nLdM;
  int nHdM;

};

#endif
