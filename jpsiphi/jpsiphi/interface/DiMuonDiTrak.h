#ifndef JpsiPhiAnalysis_DiMuonDiTrak_h
#define JpsiPhiAnalysis_DiMuonDiTrak_h


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
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "FWCore/Common/interface/TriggerNames.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"


class DiMuonDiTrakPAT : public edm::EDProducer {
 public:
  explicit DiMuonDiTrakPAT(const edm::ParameterSet&);
  ~DiMuonDiTrakPAT() override;

 private:
  void beginJob() override ;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endJob() override ;
  const pat::CompositeCandidate makeTTCandidate(const pat::PackedCandidate& trakP,
    const pat::PackedCandidate& trakN);


  // ----------member data ---------------------------
 private:

  edm::EDGetTokenT<pat::CompositeCandidateCollection> ditraks_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuons_;
  std::vector<double> DiMuonDiTrakMassCuts_;

};

#endif
