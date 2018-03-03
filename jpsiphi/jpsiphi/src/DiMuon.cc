#include "../interface/DiMuon.h"

DiMuonPAT::DiMuonPAT(const edm::ParameterSet& iConfig):
muons_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("Muons"))),
dimuonMassCuts_(iConfig.getParameter<std::vector<double>>("DiMuonCuts"))
{
  produces<pat::CompositeCandidateCollection>();
}


DiMuonPAT::~DiMuonPAT()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


void
DiMuonPAT::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;

  typedef Candidate::LorentzVector LorentzVector;

  std::unique_ptr<pat::CompositeCandidateCollection> mumuCollection(new pat::CompositeCandidateCollection);

  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByToken(muons_,muons);

  float DiMuonMassMax_ = dimuonMassCuts_[1];
  float DiMuonMassMin_ = dimuonMassCuts_[0];

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

      if ( mumucand.mass() < DiMuonMassMax_ && mumucand.mass() > DiMuonMassMin_ )
        mumuCollection->push_back(TTCand);

    }
  }

  iEvent.put(std::move(mumuCollection));

}


// ------------ method called once each job just before starting event loop  ------------
void
DiMuonPAT::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
DiMuonPAT::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonPAT);
