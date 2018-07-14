// system include files
#include <memory>

// FW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// DataFormat includes
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"


class DiMuonFilter : public edm::EDProducer {
 public:
  explicit DiMuonFilter(const edm::ParameterSet&);
  ~DiMuonFilter() override {};
  UInt_t isTriggerMatched(const pat::CompositeCandidate *);
 private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  edm::EDGetTokenT<std::vector<pat::CompositeCandidate>> theOnias_;
  StringCutObjectSelector<reco::Candidate, true> SingleMuonSelection_;
  StringCutObjectSelector<reco::Candidate, true> DiMuonSelection_;
  bool do_trigger_match_;
  std::vector<std::string> HLTFilters_;
};

DiMuonFilter::DiMuonFilter(const edm::ParameterSet& iConfig):
  theOnias_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("OniaTag"))),
  SingleMuonSelection_(iConfig.existsAs<std::string>("singlemuonSelection") ? iConfig.getParameter<std::string>("singlemuonSelection") : ""),
  DiMuonSelection_(iConfig.existsAs<std::string>("dimuonSelection") ? iConfig.getParameter<std::string>("dimuonSelection") : ""),
  do_trigger_match_(iConfig.getParameter<bool>("do_trigger_match")),
  HLTFilters_(iConfig.getParameter<std::vector<std::string>>("HLTFilters"))
{
  produces<pat::CompositeCandidateCollection>();
}

UInt_t DiMuonFilter::isTriggerMatched(const pat::CompositeCandidate *diMuon_cand) {
  const pat::Muon* highMuon = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("highMuon"));
  const pat::Muon* lowMuon = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("lowMuon"));
  UInt_t matched = 0;  // if no list is given, is not matched

// if matched a given trigger, set the bit, in the same order as listed
  for (unsigned int iTr = 0; iTr<HLTFilters_.size(); iTr++ ) {
     // std::cout << HLTFilters_[iTr] << std::endl;
     const pat::TriggerObjectStandAloneCollection mu1HLTMatches = highMuon->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
     const pat::TriggerObjectStandAloneCollection mu2HLTMatches = lowMuon->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
     if (!mu1HLTMatches.empty() && !mu2HLTMatches.empty()) matched += (1<<iTr);
     // if (!mu1HLTMatches.empty() || !mu2HLTMatches.empty()) std::cout << " MMM " << std::endl;
  }

  // const pat::TriggerObjectStandAloneCollection highMuonCollection = highMuon->triggerObjectMatches();
  // const pat::TriggerObjectStandAloneCollection lowMuonCollection = lowMuon->triggerObjectMatches();
  //
  // for ( size_t i = 0; i < highMuonCollection.size(); ++i )
  //   for ( size_t j = 0; j < highMuon->triggerObjectMatch(i)->filterLabels().size(); ++j )
  //     std::cout << (highMuon->triggerObjectMatch(i)->filterLabels())[j] << std::endl;

  // std::cout << "Triggers matched : " << matched << std::endl;
  // std::cout << "Sizes : " << highMuonCollection.size() << " - " << lowMuonCollection.size() << std::endl;
  return matched;
}

// ------------ method called to produce the data  ------------
void DiMuonFilter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  std::unique_ptr<pat::CompositeCandidateCollection> mumuOutput(new pat::CompositeCandidateCollection);
  edm::Handle<pat::CompositeCandidateCollection> onias_;
  iEvent.getByToken(theOnias_, onias_);
  if (onias_.isValid() && !onias_->empty()) {

    const pat::CompositeCandidate *ionia = nullptr;
    for (size_t ii = 0, nn=onias_->size(); ii < nn; ii++ ) {
       ionia = &(onias_->at(ii));

       if (ionia && DiMuonSelection_(*ionia) &&
           SingleMuonSelection_(*ionia->daughter("highMuon")) &&
           SingleMuonSelection_(*ionia->daughter("lowMuon")) &&
           ( !do_trigger_match_ || isTriggerMatched(ionia))
          ) mumuOutput->push_back(*ionia);
    }
  }
  iEvent.put(std::move(mumuOutput));
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonFilter);
