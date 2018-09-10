#include "../interface/DiMuonIncoherentProducer.h"

DiMuonIncoherentProducer::DiMuonIncoherentProducer(const edm::ParameterSet& ps):
  DiMuonCollection_(consumes<pat::CompositeCandidateCollection>(ps.getParameter<edm::InputTag>("DiMuon"))),
  TrakCollection_(consumes<std::vector<pat::PackedCandidate>>(ps.getParameter<edm::InputTag>("PFCandidates"))),
  DiMuonMassCuts_(ps.getParameter<std::vector<double>>("DiMuonMassCuts")),
  TrakTrakMassCuts_(ps.getParameter<std::vector<double>>("TrakTrakMassCuts")),
  DiMuonDiTrakMassCuts_(ps.getParameter<std::vector<double>>("DiMuonDiTrakMassCuts")),
  MassTraks_(ps.getParameter<std::vector<double>>("MassTraks")),
  OnlyBest_(ps.getParameter<bool>("OnlyBest")),
  product_name_(ps.getParameter<std::string>("Product"))
{
  produces<pat::CompositeCandidateCollection>(product_name_);
  candidates = 0;
  nevents = 0;
  ndimuon = 0;
  nreco = 0;
}

void DiMuonIncoherentProducer::produce(edm::Event& event, const edm::EventSetup& esetup){

  std::unique_ptr<pat::CompositeCandidateCollection> DiMuonTTCandColl(new pat::CompositeCandidateCollection);

  edm::Handle<pat::CompositeCandidateCollection> dimuon;
  event.getByToken(DiMuonCollection_,dimuon);

  edm::Handle<std::vector<pat::PackedCandidate> > trak;
  event.getByToken(TrakCollection_,trak);

  uint ncombo = 0;
  float DiMuonMassMax_ = DiMuonMassCuts_[1];
  float DiMuonMassMin_ = DiMuonMassCuts_[0];
  float TrakTrakMassMax_ = TrakTrakMassCuts_[1];
  float TrakTrakMassMin_ = TrakTrakMassCuts_[0];
  float DiMuonDiTrakMassMax_ = DiMuonDiTrakMassCuts_[1];
  float DiMuonDiTrakMassMin_ = DiMuonDiTrakMassCuts_[0];

// Note: Dimuon cand are sorted by decreasing vertex probability then first is associated with "best" dimuon
  for (pat::CompositeCandidateCollection::const_iterator dimuonCand = dimuon->begin(); dimuonCand != dimuon->end(); ++dimuonCand){
     if ( dimuonCand->mass() < DiMuonMassMax_  && dimuonCand->mass() > DiMuonMassMin_ ) {
       const pat::Muon *pmu1 = dynamic_cast<const pat::Muon*>(dimuonCand->daughter("muon1"));
       const pat::Muon *pmu2 = dynamic_cast<const pat::Muon*>(dimuonCand->daughter("muon2"));

// loop on track candidates, make DiMuonT candidate, positive charge
       for (std::vector<pat::PackedCandidate>::const_iterator firstTrack = trak->begin(), trakend=trak->end(); firstTrack!= trakend; ++firstTrack){

         if(firstTrack->charge()==0) continue;
         if(firstTrack->pt()<0.5) continue;
	       if(fabs(firstTrack->pdgId())!=211) continue;
	       if(!(firstTrack->trackHighPurity())) continue;

         if ( IsTheSame(*firstTrack,*pmu1) || IsTheSame(*firstTrack,*pmu2)) continue;

// loop over second track candidate, negative charge
         for (std::vector<pat::PackedCandidate>::const_iterator secondTrack = firstTrack + 1; secondTrack!= trakend; ++secondTrack){

           if(secondTrack->charge()==0) continue;
           if(secondTrack->pt()<0.5) continue;
  	       if(fabs(secondTrack->pdgId())!=211) continue;
  	       if(!(secondTrack->trackHighPurity())) continue;

           if(secondTrack->charge() * firstTrack->charge() <= 0) continue;

           if (secondTrack == firstTrack) continue;
           if ( IsTheSame(*secondTrack,*pmu1) || IsTheSame(*secondTrack,*pmu2)) continue;

           pat::CompositeCandidate TTCand = makeTTCandidate(*firstTrack, *secondTrack);

           if ( TTCand.mass() < TrakTrakMassMax_ && TTCand.mass() > TrakTrakMassMin_ ) {

           pat::CompositeCandidate DiMuonTTCand = makeDiMuonTTCandidate(*dimuonCand, *&TTCand);

           if ( DiMuonTTCand.mass() < DiMuonDiTrakMassMax_ && DiMuonTTCand.mass() > DiMuonDiTrakMassMin_) {

             DiMuonTTCandColl->push_back(DiMuonTTCand);
             candidates++;
             ncombo++;
           }
        }

         }
         } // loop over second track
       }   // loop on track candidates
       if (OnlyBest_) break;
     }

  if ( ncombo != DiMuonTTCandColl->size() ) std::cout <<"ncombo ("<<ncombo<< ") != DiMuonTT ("<<DiMuonTTCandColl->size()<<")"<< std::endl;
  if ( !dimuon->empty() )  ndimuon++;
  if ( ncombo > 0 ) nreco++;
  event.put(std::move(DiMuonTTCandColl),product_name_);
  nevents++;
}

void DiMuonIncoherentProducer::endJob(){
  std::cout << "#########################################" << std::endl;
  std::cout << "DiMuonDiTrak Candidate producer report:" << std::endl;
  std::cout << "#########################################" << std::endl;
  std::cout << "Found " << nevents << " Events" << std::endl;
  std::cout << "Events with DiMuon candidates " << ndimuon << std::endl;
  std::cout << "Events with DiMuonDiTrak candidates " << nreco << std::endl;
  std::cout << "#########################################" << std::endl;
  std::cout << "Found " << candidates << " DiMuonDiTrak candidates." << std::endl;
  std::cout << "#########################################" << std::endl;
}

bool DiMuonIncoherentProducer::IsTheSame(const pat::PackedCandidate& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

const pat::CompositeCandidate DiMuonIncoherentProducer::makeDiMuonTTCandidate(
                                          const pat::CompositeCandidate& dimuon,
				          const pat::CompositeCandidate& tt
                                         ){

  pat::CompositeCandidate DiMuonTCand;
  DiMuonTCand.addDaughter(dimuon,"dimuon");
  DiMuonTCand.addDaughter(tt,"ditrak");
  DiMuonTCand.setVertex(dimuon.vertex());
  DiMuonTCand.setCharge(tt.charge());

  reco::Candidate::LorentzVector vDiMuonT = dimuon.p4() + tt.p4();
  DiMuonTCand.setP4(vDiMuonT);

  return DiMuonTCand;

}

const pat::CompositeCandidate DiMuonIncoherentProducer::makeTTCandidate(
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


reco::Candidate::LorentzVector DiMuonIncoherentProducer::convertVector(const math::XYZTLorentzVectorF& v){

  return reco::Candidate::LorentzVector(v.x(),v.y(), v.z(), v.t());
}
//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonIncoherentProducer);
