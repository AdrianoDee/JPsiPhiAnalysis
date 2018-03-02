#include "../interface/TriggersByFilters.h"

TriggersByFilters::TriggersByFilters(const edm::ParameterSet& iConfig):
triggers_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("TriggerInput"))),
triggerResults(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
HLTFilters_(iConfig.getParameter<std::vector<std::string>>("Filters"))
{
  produces<pat::CompositeCandidateCollection>();
}


TriggersByFilters::~TriggersByFilters()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


void
TriggersByFilters::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  typedef Candidate::LorentzVector LorentzVector;

  std::unique_ptr<pat::TriggerObjectStandAloneCollection> trigCollection(new pat::TriggerObjectStandAloneCollection);

  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> trigs;
  iEvent.getByToken(triggers_,trigs);

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( triggerResults , triggerResults_handle);

  const edm::TriggerNames & names = iEvent.triggerNames( *triggerResults_handle );

  for ( size_t iTrigObj = 0; iTrigObj < trigs->size(); ++iTrigObj ) {

    pat::TriggerObjectStandAlone unPackedTrigger( trigs->at( iTrigObj ) );

    unPackedTrigger.unpackPathNames( names );
    unPackedTrigger.unpackFilterLabels(iEvent,*triggerResults_handle);

    bool filtered = false;

    for (size_t i = 0; i < HLTFilters_.size(); i++)
      if(unPackedTrigger.hasFilterLabel(HLTFilters_[i]))
        filtered = true;

    if(filtered)
      trigCollection->push_back(unPackedTrigger);

  }

  iEvent.put(std::move(trigCollection));

}


// ------------ method called once each job just before starting event loop  ------------
void
TriggersByFilters::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
TriggersByFilters::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggersByFilters);
