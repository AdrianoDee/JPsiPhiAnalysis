#include "../interface/DiTrackHLTProducer.h"


float DiTrackHLTProducer::DeltaR(pat::CompositeCandidate t1, pat::TriggerObjectStandAlone t2)
{
   float p1 = t1.phi();
   float p2 = t2.phi();
   float e1 = t1.eta();
   float e2 = t2.eta();
   auto dp=std::abs(p1-p2); if (dp>Float(M_PI)) dp-=Float(2*M_PI);

   return sqrt((e1-e2)*(e1-e2) + dp*dp);
}

bool DiTrackHLTProducer::MatchByDRDPt(pat::CompositeCandidate t1, pat::TriggerObjectStandAlone t2)
{
  return (fabs(t1.pt()-t2.pt())/t2.pt()<maxDPtRel &&
	DeltaR(t1,t2) < maxDPtRel);
}

DiTrackHLTProducer::DiTrackHLTProducer(const edm::ParameterSet& ps):
  TrakCollection_(consumes<std::vector<pat::PackedCandidate>>(ps.getParameter<edm::InputTag>("PFCandidates"))),
  TrakTrakMassCuts_(ps.getParameter<std::vector<double>>("TrakTrakMassCuts")),
  MassTraks_(ps.getParameter<std::vector<double>>("MassTraks")),
  OnlyBest_(ps.getParameter<bool>("OnlyBest")),
  product_name_(ps.getParameter<std::string>("Product")),
  HLTFilters_(ps.getParameter<std::vector<std::string>>("HLTFilters")),
  triggerObj_(consumes<std::vector<pat::TriggerObjectStandAlone>>("TriggerInput"))
{

  produces<pat::CompositeCandidateCollection>(product_name_);
  produces<pat::TriggerObjectStandAloneCollection>(product_name_);
  candidates = 0;
  nevents = 0;
  ndimuon = 0;
  nreco = 0;
  maxDeltaR = 0.1;
  maxDPtRel = 10.0;
}

void DiTrackHLTProducer::produce(edm::Event& event, const edm::EventSetup& esetup){

  std::unique_ptr<pat::CompositeCandidateCollection> DiTrackColl(new pat::CompositeCandidateCollection);
  std::unique_ptr<pat::TriggerObjectStandAloneCollection> DiTriggColl(new pat::TriggerObjectStandAloneCollection);

  edm::Handle<std::vector<pat::PackedCandidate> > trakColl;
  event.getByToken(TrakCollection_,trak);

  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerColl;
  event.getByToken(triggerObj_,triggerColl);

  uint ncombo = 0;

  float
  float TrakTrakMassMax_ = TrakTrakMassCuts_[1];
  float TrakTrakMassMin_ = TrakTrakMassCuts_[0];

  pat::TriggerObjectStandAloneCollection filteredColl, matchedColl;

  //Filtering

  for (std::vector<pat::TriggerObjectStandAlone>::const_iterator trigger = triggerColl->begin(), triggerEnd=triggerColl->end(); trigger!= triggerEnd; ++triggerColl)
  {
    std::vector< std::string > thisFilters trigger->filterLabels();
    std::vector< std::string > matchFilters;

    std::set_intersection(thisFilters.begin(),thisFilters.end(),HLTFilters_.begin(),HLTFilters_.end(),back_inserter(matchFilters))

    if(interSection.size()>0)
      filteredColl.push_back(*trigger);

  }

  //Matching

  for (std::vector<pat::PackedCandidate>::const_iterator trak = trakColl->begin(), trakend=trakColl->end(); trak!= trakend; ++trakColl)
  {
    bool matched = false;
    for (std::vector<pat::TriggerObjectStandAlone>::const_iterator trigger = filteredColl->begin(), triggerEnd=filteredColl->end(); trigger!= triggerEnd; ++filteredColl)
    {
      if(MatchByDRDPt(*trak,*trigger))
      {
        if(matched)
        {
          if(DeltaR(*trak,matchedColl.back()) > DeltaR(*trak,*trigger))
          {
            matchedColl.pop_back();
            matchedColl.push_back(*trigger);
          }
        }
        matched = true;
      }
    }
  }
//
//   for (std::vector<pat::TriggerObjectStandAlone>::const_iterator trigger = triggerColl->begin(), triggerEnd=triggerColl->end(); trigger!= triggerEnd; ++triggerColl)
//   {
//
//   }
//
//   for (std::vector<pat::TriggerObjectStandAlone>::const_iterator trigger = filteredColl->begin(), triggerEnd=filteredColl->end(); trigger!= triggerEnd; ++filteredColl)
//   {
//     if(posTrack->charge()==0) continue;
//
//
//   }
//
//
// // Note: Dimuon cand are sorted by decreasing vertex probability then first is associated with "best" dimuon
//   for (pat::CompositeCandidateCollection::const_iterator dimuonCand = dimuon->begin(); dimuonCand != dimuon->end(); ++dimuonCand){
//      if ( dimuonCand->mass() < DiMuonMassMax_  && dimuonCand->mass() > DiMuonMassMin_ ) {
//        const pat::Muon *pmu1 = dynamic_cast<const pat::Muon*>(dimuonCand->daughter("muon1"));
//        const pat::Muon *pmu2 = dynamic_cast<const pat::Muon*>(dimuonCand->daughter("muon2"));
//
// // loop on track candidates, make DiMuonT candidate, positive charge
//        for (std::vector<pat::PackedCandidate>::const_iterator posTrack = trak->begin(), trakend=trak->end(); posTrack!= trakend; ++posTrack){
//
//          if(posTrack->charge()==0) continue;
//          if(posTrack->pt()<0.5) continue;
// 	       if(fabs(posTrack->pdgId())!=211) continue;
// 	       if(!(posTrack->trackHighPurity())) continue;
//
//          if ( IsTheSame(*posTrack,*pmu1) || IsTheSame(*posTrack,*pmu2) || posTrack->charge() < 0 ) continue;
//
// // loop over second track candidate, negative charge
//          for (std::vector<pat::PackedCandidate>::const_iterator negTrack = trak->begin(); negTrack!= trakend; ++negTrack){
//
//            if(negTrack->charge()==0) continue;
//            if(negTrack->pt()<0.5) continue;
//   	       if(fabs(negTrack->pdgId())!=211) continue;
//   	       if(!(negTrack->trackHighPurity())) continue;
//
//            if (negTrack == posTrack) continue;
//            if ( IsTheSame(*negTrack,*pmu1) || IsTheSame(*negTrack,*pmu2) || negTrack->charge() > 0 ) continue;
//
//            pat::CompositeCandidate TTCand = makeTTCandidate(*posTrack, *negTrack);
//
//            if ( TTCand.mass() < TrakTrakMassMax_ && TTCand.mass() > TrakTrakMassMin_ ) {
//
//            pat::CompositeCandidate DiMuonTTCand = makeDiMuonTTCandidate(*dimuonCand, *&TTCand);
//
//            if ( DiMuonTTCand.mass() < DiMuonDiTrakMassMax_ && DiMuonTTCand.mass() > DiMuonDiTrakMassMin_) {
//
//              DiMuonTTCandColl->push_back(DiMuonTTCand);
//              candidates++;
//              ncombo++;
//            }
//         }
//
//          }
//          } // loop over second track
//        }   // loop on track candidates
//        if (OnlyBest_) break;
//      }

  if ( ncombo != DiTrackColl->size() ) std::cout <<"ncombo ("<<ncombo<< ") != DiMuonTT ("<<DiTrackColl->size()<<")"<< std::endl;
  if ( ncombo > 0 ) nreco++;
  event.put(std::move(DiTrackColl),product_name_);

  nevents++;
}

void DiTrackHLTProducer::endJob(){
  std::cout << "###########################" << std::endl;
  std::cout << "DiMuonDiTrak Candidate producer report:" << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << nevents << " Events" << std::endl;
  std::cout << "Events with DiMuon candidates " << ndimuon << std::endl;
  std::cout << "Events with DiMuonDiTrak candidates " << nreco << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << candidates << " DiMuonDiTrak candidates." << std::endl;
  std::cout << "###########################" << std::endl;
}

bool DiTrackHLTProducer::IsTheSame(const pat::PackedCandidate& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

const pat::CompositeCandidate DiTrackHLTProducer::makeDiMuonTTCandidate(
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

const pat::CompositeCandidate DiTrackHLTProducer::makeTTCandidate(
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


reco::Candidate::LorentzVector DiTrackHLTProducer::convertVector(const math::XYZTLorentzVectorF& v){

  return reco::Candidate::LorentzVector(v.x(),v.y(), v.z(), v.t());
}
//define this as a plug-in
DEFINE_FWK_MODULE(DiTrackHLTProducer);
