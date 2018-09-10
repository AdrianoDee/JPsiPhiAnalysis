#include "../interface/DoubleDiMuonProducer.h"

DoubleDiMuonProducer::DoubleDiMuonProducer(const edm::ParameterSet& ps):
  HighDiMuonCollection_(consumes<pat::CompositeCandidateCollection>(ps.getParameter<edm::InputTag>("HighDiMuonCollection"))),
  LowDiMuonCollection_(consumes<pat::CompositeCandidateCollection>(ps.getParameter<edm::InputTag>("LowDiMuonCollection"))),
  HighDiMuonMassCuts_(ps.getParameter<std::vector<double>>("HighDiMuonMassCuts")),
  LowDiMuonMassCuts_(ps.getParameter<std::vector<double>>("LowDiMuonMassCuts")),
  DoubleDiMuonMassCuts_(ps.getParameter<std::vector<double>>("DoubleDiMuonMassCuts")),
  addMCTruth_(ps.getParameter<bool>("AddMCTruth"))
{
  produces<pat::CompositeCandidateCollection>("DoubleDiMuonCandidates");
  candidates = 0;
  nevents = 0;
  nLdM = 0;
  nHdM = 0;
}

std::tuple<int, float, float>
DoubleDiMuonProducer::findJpsiMCInfo(reco::GenParticleRef genJpsi) {

  // std::cout << "findJpsiMCInfo 1 " << std::endl;
  int momJpsiID = 0;
  float trueLife = -99.;
  float isPrompt = -99.;
  if (genJpsi->numberOfMothers()>0) {

    // std::cout << "findJpsiMCInfo 1 " << std::endl;

    TVector3 trueVtx(0.0,0.0,0.0);
    TVector3 trueP(0.0,0.0,0.0);
    TVector3 trueVtxMom(0.0,0.0,0.0);

    trueVtx.SetXYZ(genJpsi->vertex().x(),genJpsi->vertex().y(),genJpsi->vertex().z());
    trueP.SetXYZ(genJpsi->momentum().x(),genJpsi->momentum().y(),genJpsi->momentum().z());

    reco::GenParticleRef Jpsimom = genJpsi->motherRef();       // find mothers
    // std::cout << "findJpsiMCInfo 1 " << std::endl;
    if (Jpsimom.isNull()) {
      std::tuple<int, float, float> result = std::make_tuple(momJpsiID, trueLife,isPrompt);
      return result;
    } else
    {
    momJpsiID = Jpsimom->pdgId();
    isPrompt = Jpsimom->isPromptDecayed();
    trueVtxMom.SetXYZ(Jpsimom->vertex().x(),Jpsimom->vertex().y(),Jpsimom->vertex().z());
    TVector3 vdiff = trueVtx - trueVtxMom;
    trueLife = vdiff.Perp()*genJpsi->mass()/trueP.Perp();
  }
}
  std::tuple<int,float,float> result = std::make_tuple(momJpsiID, trueLife,isPrompt);
  return result;

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

                if (addMCTruth_) {
                  reco::GenParticleRef genMu1 = phiMu1->genParticleRef();
                  reco::GenParticleRef genMu2 = phiMu2->genParticleRef();
                  // reco::GenParticleRef genKaon1 = posTrack.genParticleRef();
                  // reco::GenParticleRef genKaon2 = negTrack.genParticleRef();

                  if (genMu1.isNonnull() && genMu2.isNonnull()) {
                    if (genMu1->numberOfMothers()>0 && genMu2->numberOfMothers()>0){
                      reco::GenParticleRef mumu_mom1 = genMu1->motherRef();
                      reco::GenParticleRef mumu_mom2 = genMu2->motherRef();

                      if (mumu_mom1.isNonnull() && (mumu_mom1 == mumu_mom2)) {

                        std::tuple<int,float,float> MCinfo = findJpsiMCInfo(mumu_mom1);
                        DoubleDiMuonCandidate.addUserInt("phiGenPdgId",mumu_mom1->pdgId());
                        DoubleDiMuonCandidate.addUserFloat("phiPpdlTrue",std::get<1>(MCinfo));
                        DoubleDiMuonCandidate.addUserInt("xGenPdgId",std::get<0>(MCinfo));
                        DoubleDiMuonCandidate.addUserFloat("xGenIsPrompt",std::get<2>(MCinfo));
                      } else {
                        DoubleDiMuonCandidate.addUserInt("phiGenPdgId",0.0);
                        DoubleDiMuonCandidate.addUserFloat("phiPpdlTrue",-99.0);
                        DoubleDiMuonCandidate.addUserInt("xGenPdgId",0.0);
                        DoubleDiMuonCandidate.addUserFloat("xGenIsPrompt",-99.0);
                      }

                    }
                 } else {
                   DoubleDiMuonCandidate.addUserInt("phiGenPdgId",0.0);
                   DoubleDiMuonCandidate.addUserFloat("phiPpdlTrue",-99.0);
                   DoubleDiMuonCandidate.addUserInt("xGenPdgId",0.0);
                   DoubleDiMuonCandidate.addUserFloat("xGenIsPrompt",-99.0);
                  }
                }

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
  std::cout << "#########################################" << std::endl;
  std::cout << "DoubleDiMuon Candidate producer report:" << std::endl;
  std::cout << "#########################################" << std::endl;
  std::cout << "Found " << nevents << " Events" << std::endl;
  std::cout << "Events with LowDiMuon  candidates " << nLdM << std::endl;
  std::cout << "Events with HighDiMuon candidates " << nHdM << std::endl;
  std::cout << "#########################################" << std::endl;
  std::cout << "Found " << candidates << " DoubleDiMuon candidates." << std::endl;
  std::cout << "#########################################" << std::endl;
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
