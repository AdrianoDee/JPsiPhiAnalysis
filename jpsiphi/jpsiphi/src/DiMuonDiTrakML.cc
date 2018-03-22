#include "../interface/DiMuonDiTrakML.h"

DiMuonDiTrakML::DiMuonDiTrakML(const edm::ParameterSet& iConfig):
muons_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("Muons"))),
traks_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("Track"))),
// thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpot"))),
// thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PrimaryVertex"))),
// DiMuonMassCuts_(iConfig.getParameter<std::vector<double>>("DiMuonCuts")),
// DiTrakMassCuts_(iConfig.getParameter<std::vector<double>>("DiTrakCuts")),
// DiMuonDiTrakMassCuts_(iConfig.getParameter<std::vector<double>>("DiMuonDiTrakCuts")),
// massCands_(iConfig.getParameter<std::vector<double>>("CandsMasses"))
{
  // produces<pat::CompositeCandidateCollection>();
}


DiMuonDiTrakML::~DiMuonDiTrakML()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


void
DiMuonDiTrakML::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(traks_,tracks);

  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByToken(muons_,muons);

  for(reco::TrackCollection::const_iterator itTrack = tracks->begin();itTrack != tracks->end(); ++itTrack )
    std::cout<<itTrack->recHitsSize()<<std::endl;


  // std::sort(mmttCollection->begin(),mmttCollection->end(),vPComparator_);
  // iEvent.put(std::move(mmttCollection));

}


// ------------ method called once each job just before starting event loop  ------------
void
DiMuonDiTrakML::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
DiMuonDiTrakML::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonDiTrakML);
