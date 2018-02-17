#include "jpsiphi/jpsiphi/interface/DoubleDiMuonProducer.h"

DoubleDiMuonProducer::DoubleDiMuonProducer(const edm::ParameterSet& ps):
  HighDiMuonCollection_(consumes<pat::CompositeCandidateCollection>(ps.getParameter<edm::InputTag>("HighDiMuonCollection"))),
  LowDiMuonCollection_(consumes<pat::CompositeCandidateCollection>(ps.getParameter<edm::InputTag>("LowDiMuonCollection"))),
  HighDiMuonMassCuts_(ps.getParameter<std::vector<double>>("HighDiMuonMassCuts")),
  LowDiMuonMassCuts_(ps.getParameter<std::vector<double>>("LowDiMuonMassCuts")),
  DoubleDiMuonMassCuts_(ps.getParameter<std::vector<double>>("DoubleDiMuonMassCuts"))
{
  produces<pat::CompositeCandidateCollection>("DoubleDiMuonCandidates");
  candidates = 0;
  nevents = 0;
  nLdM = 0;
  nHdM = 0;
}

void DoubleDiMuonProducer::produce(edm::Event& event, const edm::EventSetup& esetup){

  std::unique_ptr<pat::CompositeCandidateCollection> DoubleDiMuonCandColl(new pat::CompositeCandidateCollection);

  edm::Handle<pat::CompositeCandidateCollection> highDiMuon;
  event.getByToken(HighDiMuonCollection_,highDiMuon);

  edm::Handle<pat::CompositeCandidateCollection> lowDiMuon;
  event.getByToken(LowDiMuonCollection_,lowDiMuon);

  float HighDiMuonMassMax_ = HighDiMuonMassCuts_[1];
  float HighDiMuonMassMin_ = HighDiMuonMassCuts_[0];
  float LowDiMuonMassMax_ = LowDiMuonMassCuts_[1];
  float LowDiMuonMassMin_ = LowDiMuonMassCuts_[0];

  float DoubleDiMuonMassMax_ = DoubleDiMuonMassCuts_[1];
  float DoubleDiMuonMassMin_ = DoubleDiMuonMassCuts_[0];

// Note: Dimuon cand are sorted by decreasing vertex probability then first is associated with "best" dimuon
  //Looking for J/Psi

  for (pat::CompositeCandidateCollection::const_iterator highCand = highDiMuon->begin(); highCand != highDiMuon->end(); ++highCand){

     if ( highCand->mass() < HighDiMuonMassMax_  && highCand->mass() > HighDiMuonMassMin_ ) {
       const pat::Muon *jPsiMu1 = dynamic_cast<const pat::Muon*>(highCand->daughter("muon1"));
       const pat::Muon *jPsiMu2 = dynamic_cast<const pat::Muon*>(highCand->daughter("muon2"));


       for (pat::CompositeCandidateCollection::const_iterator lowCand = lowDiMuon->begin(); lowCand != lowDiMuon->end(); ++lowCand){

          if ( lowCand->mass() < LowDiMuonMassMax_  && lowCand->mass() > LowDiMuonMassMin_ ) {

            const pat::Muon *phiMu1 = dynamic_cast<const pat::Muon*>(lowCand->daughter("muon1"));
            const pat::Muon *phiMu2 = dynamic_cast<const pat::Muon*>(lowCand->daughter("muon2"));

            if( phiMu1 == phiMu2 || phiMu1 == jPsiMu1 || phiMu1 == jPsiMu2 ) continue;
            if( phiMu2 == jPsiMu1 || phiMu2 == jPsiMu2 ) continue;
            if( jPsiMu1 == jPsiMu2 ) continue;

            pat::CompositeCandidate DoubleDiMuonCandidate = makeCandidate(*lowCand, *highCand);

            if(DoubleDiMuonCandidate.charge() != 0.0) continue;

            if ( DoubleDiMuonCandidate.mass() < DoubleDiMuonMassMax_ && DoubleDiMuonCandidate.mass() > DoubleDiMuonMassMin_)
              {
                candidates++;
                DoubleDiMuonCandColl->push_back(DoubleDiMuonCandidate);
              }
            }
          }
        }
      }
     // if (OnlyBest_) break;

     if ( !(highDiMuon->empty()) )  nLdM++;
     if ( !(lowDiMuon->empty()) )  nHdM++;

     event.put(std::move(DoubleDiMuonCandColl),"DoubleDiMuonCandidates");
     nevents++;
  }


void DoubleDiMuonProducer::endJob(){
  std::cout << "###########################" << std::endl;
  std::cout << "DoubleDiMuon Candidate producer report:" << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << nevents << " Events" << std::endl;
  std::cout << "Events with LowDiMuon  candidates " << nLdM << std::endl;
  std::cout << "Events with HighDiMuon candidates " << nHdM << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << candidates << " DoubleDiMuon candidates." << std::endl;
  std::cout << "###########################" << std::endl;
}


const pat::CompositeCandidate DoubleDiMuonProducer::makeCandidate(const pat::CompositeCandidate& lowDiMuon,
  const pat::CompositeCandidate& higDiMuon){
    pat::CompositeCandidate xCand;
    xCand.addDaughter(lowDiMuon,"lowdimuon");
    xCand.addDaughter(higDiMuon,"higdimuon");
    reco::Candidate::LorentzVector vX = lowDiMuon.p4() + higDiMuon.p4();
    xCand.setP4(vX);
    return xCand;
  }

reco::Candidate::LorentzVector DoubleDiMuonProducer::convertVector(const math::XYZTLorentzVectorF& v){

  return reco::Candidate::LorentzVector(v.x(),v.y(), v.z(), v.t());
}
//define this as a plug-in
DEFINE_FWK_MODULE(DoubleDiMuonProducer);
