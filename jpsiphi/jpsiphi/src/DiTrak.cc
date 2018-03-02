#include "../interface/DiTrak.h"

DiTrakPAT::DiTrakPAT(const edm::ParameterSet& iConfig):
traks_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Traks"))),
ditrakMassCuts_(iConfig.getParameter<std::vector<double>>("DiTrakCuts")),
massTraks_(iConfig.getParameter<std::vector<double>>("TraksMasses")),
{
  produces<pat::CompositeCandidateCollection>();
}


DiTrakPAT::~DiTrakPAT()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


void
DiTrakPAT::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  typedef Candidate::LorentzVector LorentzVector;

  std::unique_ptr<pat::CompositeCandidateCollection> trakCollection(new pat::CompositeCandidateCollection);

  edm::Handle<std::vector<pat::PackedCandidate> > traks;
  iEvent.getByToken(traks_,traks);

  float TrakTrakMassMax_ = ditrakMassCuts_[1];
  float TrakTrakMassMin_ = ditrakMassCuts_[0];

  for (size_t i = 0; i < traks->size(); i++)
  {
    auto posTrack = traks->at(i);

    if(posTrack.charge() <= 0 ) continue;
    if(posTrack.pt()<0.5) continue;
    if(fabs(posTrack.pdgId())!=211) continue;


    for (size_t j = 0; j < traks->size(); j++){

      if (i == j) continue;

      auto negTrack = traks->at(j);

      if(negTrack.charge() >= 0 ) continue;
      if(negTrack.pt()<0.5) continue;
      if(fabs(negTrack.pdgId())!=211) continue;

      pat::CompositeCandidate TTCand = makeTTCandidate(posTrack,negTrack);

      if ( TTCand.mass() < TrakTrakMassMax_ && TTCand.mass() > TrakTrakMassMin_ ) {

        trakCollection->push_back(TTCand);

      }
    } // loop over second track
  }

  iEvent.put(std::move(oniaOutput));

}


// ------------ method called once each job just before starting event loop  ------------
void
DiTrakPAT::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
DiTrakPAT::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiTrakPAT);
