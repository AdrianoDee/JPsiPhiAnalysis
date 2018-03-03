#include "../interface/DiMuonDiTrak.h"

DiMuonDiTrakPAT::DiMuonDiTrakPAT(const edm::ParameterSet& iConfig):
dimuons_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("DiMuons"))),
ditraks_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("DiTraks"))),
DiMuonDiTrakMassCuts_(iConfig.getParameter<std::vector<double>>("DiMuonDiTrakCuts"))
{
  produces<pat::CompositeCandidateCollection>();
  muon_mass = 0.1056583715;
}


DiMuonDiTrakPAT::~DiMuonDiTrakPAT()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


void
DiMuonDiTrakPAT::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;

  typedef Candidate::LorentzVector LorentzVector;

  std::unique_ptr<pat::CompositeCandidateCollection> mumuCollection(new pat::CompositeCandidateCollection);

  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByToken(muons_,muons);

  float DiMuonDiTrakMassMax_ = DiMuonDiTrakMassCuts_[1];
  float DiMuonDiTrakMassMin_ = DiMuonDiTrakMassCuts_[0];

  for(View<pat::Muon>::const_iterator mNeg = muons->begin(), itend = muons->end(); mNeg != itend; ++mNeg){

    if(mNeg->charge()>=0.0) continue;

    for(View<pat::Muon>::const_iterator mPos = muons->begin(), itend = muons->end(); mPos != itend; ++mPos)
    {
      if(mNeg == mPos) continue;

      if(mPos->charge()<=0.0) continue;

      if (!(mNeg->track().isNonnull() && mPos->track().isNonnull())) continue;

      pat::CompositeCandidate mumucand;

      mumucand.addDaughter(*mNeg,"muonN");
      mumucand.addDaughter(*mPos,"muonP");

      LorentzVector mumu = mNeg->p4() + mPos->p4();
      TLorentzVector mu1, mu2,mumuP4;

      mu1.SetXYZM(mNeg->track()->px(),mNeg->track()->py(),mNeg->track()->pz(),muon_mass);
      mu2.SetXYZM(mPos->track()->px(),mPos->track()->py(),mPos->track()->pz(),muon_mass);

      mumuP4=mu1+mu2;
      mumucand.setP4(mumu);
      mumucand.setCharge(mNeg->charge()+mPos->charge());

      if ( mumucand.mass() < DiMuonDiTrakMassMax_ && mumucand.mass() > DiMuonDiTrakMassMin_ )
        mumuCollection->push_back(mumucand);

    }
  }

  iEvent.put(std::move(mumuCollection));

}


// ------------ method called once each job just before starting event loop  ------------
void
DiMuonDiTrakPAT::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
DiMuonDiTrakPAT::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonDiTrakPAT);
