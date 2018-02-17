/**
   \file
   Declaration of OniaTrakTrakProducer

   \author Alberto Sanchez-Hernandez
   \date 2 Mar 2014
*/

#ifndef __PsiPhiFourMuonsProducer_h_
#define __PsiPhiFourMuonsProducer_h_

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

class PsiPhiFourMuonsProducer : public edm::EDProducer {

 public:
  explicit PsiPhiFourMuonsProducer(const edm::ParameterSet& ps);

 private:

  void produce(edm::Event& event, const edm::EventSetup& esetup) override;

  void endJob() override;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> PsiCollection_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> PhiCollection_;
  std::vector<double> JPsiMassCuts_;
  std::vector<double> PhiMassCuts_;
  std::vector<double> FourOniaMassCuts_;

  reco::Candidate::LorentzVector convertVector(const math::XYZTLorentzVectorF& v);

  const pat::CompositeCandidate makeCandidate(const pat::CompositeCandidate& phi,
    const pat::CompositeCandidate& jpsi);

  int candidates;
  int nevents;
  int nPhi;
  int nJps;

};

#endif
