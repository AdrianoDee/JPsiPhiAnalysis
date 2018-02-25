#include "../interface/DiTrackHLTProducer.h"


float DiTrackHLTProducer::DeltaR(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2)
{
   float p1 = t1.phi();
   float p2 = t2.phi();
   float e1 = t1.eta();
   float e2 = t2.eta();
   auto dp=std::abs(p1-p2); if (dp>float(M_PI)) dp-=float(2*M_PI);

   return sqrt((e1-e2)*(e1-e2) + dp*dp);
}

bool DiTrackHLTProducer::MatchByDRDPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2)
{
  return (fabs(t1.pt()-t2.pt())/t2.pt()<maxDPtRel &&
	DeltaR(t1,t2) < maxDPtRel);
}

DiTrackHLTProducer::DiTrackHLTProducer(const edm::ParameterSet& ps):
  TrakCollection_(consumes<std::vector<pat::PackedCandidate>>(ps.getParameter<edm::InputTag>("PFCandidates"))),
  TriggerCollection_(consumes<std::vector<pat::TriggerObjectStandAlone>>(ps.getParameter<edm::InputTag>("TriggerInput"))),
  TrakTrakMassCuts_(ps.getParameter<std::vector<double>>("TrakTrakMassCuts")),
  MassTraks_(ps.getParameter<std::vector<double>>("MassTraks")),
  OnlyBest_(ps.getParameter<bool>("OnlyBest")),
  TTCandidate_name_(ps.getParameter<std::string>("TTCandidate_name")),
  TTTrigger_name_(ps.getParameter<std::string>("TTTrigger_name")),
  HLTFilters_(ps.getParameter<std::vector<std::string>>("HLTFilters")),
{

  produces<pat::CompositeCandidateCollection>(TTCandidate_name_);
  produces<pat::TriggerObjectStandAloneCollection>(TTTrigger_name_);
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
  event.getByToken(TrakCollection_,trakColl);

  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerColl;
  event.getByToken(TriggerCollection_,triggerColl);

  uint ncombo = 0;

  float TrakTrakMassMax_ = TrakTrakMassCuts_[1];
  float TrakTrakMassMin_ = TrakTrakMassCuts_[0];

  pat::TriggerObjectStandAloneCollection filteredColl, matchedColl;
  std::vector< pat::PackedCandidate> filteredTracks;
  //Filtering

  for (std::vector<pat::TriggerObjectStandAlone>::const_iterator trigger = triggerColl->begin(), triggerEnd=triggerColl->end(); trigger!= triggerEnd; ++trigger)
  {
    std::vector< std::string > thisFilters = trigger->filterLabels();
    std::vector< std::string > matchFilters;

    std::set_intersection(thisFilters.begin(),thisFilters.end(),HLTFilters_.begin(),HLTFilters_.end(),back_inserter(matchFilters));

    if(matchFilters.size()>0)
      filteredColl.push_back(*trigger);

  }

  //Matching

  for (std::vector<pat::PackedCandidate>::const_iterator trak = trakColl->begin(), trakend=trakColl->end(); trak!= trakend; ++trak)
  {
    bool matched = false;
    for (std::vector<pat::TriggerObjectStandAlone>::const_iterator trigger = filteredColl.begin(), triggerEnd=filteredColl.end(); trigger!= triggerEnd; ++trigger)
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

        if(!matched)
          filteredTracks.push_back(*trak);

        matched = true;
      }
    }
  }

  std::cout << matchedColl.size() << " vs " << filteredTracks.size() << std::endl;
  // for (std::vector<pat::PackedCandidate>::const_iterator posTrack = filteredTracks.begin(), trakend=filteredTracks.end(); posTrack!= trakend; ++posTrack)
  for (size_t i = 0; i < filteredTracks.size(); i++) {
  {
           auto posTrack = filteredTracks[i];
           if(posTrack.charge()==0) continue;
           if(posTrack.pt()<0.5) continue;
  	       if(fabs(posTrack.pdgId())!=211) continue;
  	       if(!(posTrack.trackHighPurity())) continue;

           if ( IsTheSame(*posTrack,*pmu1) || IsTheSame(*posTrack,*pmu2) || posTrack.charge() < 0 ) continue;

  // loop over second track candidate, negative charge
           // for (std::vector<pat::PackedCandidate>::const_iterator negTrack = trak->begin(); negTrack!= trakend; ++negTrack){
           for (size_t j = 0; j < filteredTracks.size(); j++) {
             auto negTrack = filteredTracks[j];
             if(negTrack.charge()==0) continue;
             if(negTrack.pt()<0.5) continue;
    	       if(fabs(negTrack.pdgId())!=211) continue;
    	       if(!(negTrack.trackHighPurity())) continue;

             if (i == j) continue;
             if ( IsTheSame(negTrack,*pmu1) || IsTheSame(negTrack,*pmu2) || negTrack.charge() > 0 ) continue;

             pat::CompositeCandidate TTCand = makeTTCandidate(posTrack,negTrack);
             pat::CompositeCandidate TTTrigger = makeTTTriggerCandidate(matchedColl[i],matchedColl[j])
             if ( TTCand.mass() < TrakTrakMassMax_ && TTCand.mass() > TrakTrakMassMin_ ) {

               DiTrackColl->push_back(TTCand);
               DiTriggColl->push_back(TTTrigger);

             }
           } // loop over second track
         }   // loop on track candidates

  // if ( ncombo != DiTrackColl->size() ) std::cout <<"ncombo ("<<ncombo<< ") != DiMuonTT ("<<DiTrackColl->size()<<")"<< std::endl;
  // if ( ncombo > 0 ) nreco++;
  event.put(std::move(DiTrackColl),TTCandidate_name_);
  event.put(std::move(DiTriggColl),TTTrigger_name_);

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

const pat::CompositeCandidate DiTrackHLTProducer::makeTTTriggerCandidate(
                                          const pat::TriggerObjectStandAlone& trakP,
                                          const pat::TriggerObjectStandAlone& trakN
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
