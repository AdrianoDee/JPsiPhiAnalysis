/**
   \file
   Declaration of OniaPFPFProducer

   \author Alberto Sanchez-Hernandez
   \date 2 Mar 2014
*/

#ifndef __OniaPFPFProducer_h_
#define __OniaPFPFProducer_h_

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
#include <TLorentzVector.h>
#include <vector>

/**
   Create a HF candidate by mathing Onia(chi,psi,etc.) and a track (K, pi, etc.)
 */

class OniaPFPFProducer : public edm::EDProducer {

 public:
  explicit OniaPFPFProducer(const edm::ParameterSet& ps);

 private:

  void produce(edm::Event& event, const edm::EventSetup& esetup) override;

  void endJob() override;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> OniaCollection_;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate>> PFCandCollection_;
  std::vector<double> OniaMassCuts_;
  std::vector<double> TrakTrakMassCuts_;
  std::vector<double> OniaPFPFMassCuts_;
  std::vector<double> MassTraks_;
  bool OnlyBest_;

  reco::Candidate::LorentzVector convertVector(const math::XYZTLorentzVectorF& v);
  bool IsTheSame(const pat::PackedCandidate& tk, const pat::Muon& mu);
  const pat::CompositeCandidate makeOniaTTCandidate(const pat::CompositeCandidate& onia,
						    const pat::CompositeCandidate& tt);
  const pat::CompositeCandidate makeTTCandidate(const pat::PackedCandidate& trak1,
                                                const pat::PackedCandidate& trak2);
  int candidates;
  int nevents;
  int nonia;
  int nreco;
};

#endif // __OniaPFPFProducer_h_
