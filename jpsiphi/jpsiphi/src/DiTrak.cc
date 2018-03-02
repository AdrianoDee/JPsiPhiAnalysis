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

const pat::CompositeCandidate DiTrakPAT::makeTTCandidate(
                                          const pat::PackedCandidate& trakP,
                                          const pat::PackedCandidate& trakN
                                         ){

  pat::CompositeCandidate TTCand;
  TTCand.addDaughter(trakP,"trakP");
  TTCand.addDaughter(trakN,"trakN");
  TTCand.setCharge(trakP.charge()+trakN.charge());

  double m_kaon1 = MassTraks_[0];
  math::XYZVector mom_kaon1 = trakP.momentum();
  double e_kaon1 = sqrt(m_kaon1*m_kaon1 + mom_kaon1.Mag2());
  math::XYZTLorentzVector p4_kaon1 = math::XYZTLorentzVector(mom_kaon1.X(),mom_kaon1.Y(),mom_kaon1.Z(),e_kaon1);
  double m_kaon2 = MassTraks_[1];
  math::XYZVector mom_kaon2 = trakN.momentum();
  double e_kaon2 = sqrt(m_kaon2*m_kaon2 + mom_kaon2.Mag2());
  math::XYZTLorentzVector p4_kaon2 = math::XYZTLorentzVector(mom_kaon2.X(),mom_kaon2.Y(),mom_kaon2.Z(),e_kaon2);
  reco::Candidate::LorentzVector vTT = p4_kaon1 + p4_kaon2;
  TTCand.setP4(vTT);

  return TTCand;
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

  iEvent.put(std::move(trakCollection));

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
