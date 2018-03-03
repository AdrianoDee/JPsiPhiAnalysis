#include "../interface/DiMuonDiTrakProducerHLT.h"




DiMuonDiTrakProducerHLT::DiMuonDiTrakProducerHLT(const edm::ParameterSet& iConfig):
  DiMuonCollection_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("DiMuon"))),
  TrakCollection_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("PFCandidates"))),
  TriggerCollection_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("TriggerInput"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  DiMuonMassCuts_(iConfig.getParameter<std::vector<double>>("DiMuonMassCuts")),
  TrakTrakMassCuts_(iConfig.getParameter<std::vector<double>>("TrakTrakMassCuts")),
  DiMuonDiTrakMassCuts_(iConfig.getParameter<std::vector<double>>("DiMuonDiTrakMassCuts")),
  MassTraks_(iConfig.getParameter<std::vector<double>>("MassTraks")),
  MaxDeltaRPt_(iConfig.getParameter<std::vector<double>>("MaxDeltaRPt")),
  OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  product_name_(iConfig.getParameter<std::string>("Product")),
  HLTs_(iConfig.getParameter<std::vector<std::string>>("HLTs")),
  HLTFilters_(iConfig.getParameter<std::vector<std::string>>("Filters"))
{
  produces<pat::CompositeCandidateCollection>(product_name_);
  candidates = 0;
  nevents = 0;
  ndimuon = 0;
  nreco = 0;
}

void DiMuonDiTrakProducerHLT::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

  std::unique_ptr<pat::CompositeCandidateCollection> DiMuonTTCandColl(new pat::CompositeCandidateCollection);

  edm::Handle<pat::CompositeCandidateCollection> dimuon;
  iEvent.getByToken(DiMuonCollection_,dimuon);

  edm::Handle<std::vector<pat::PackedCandidate> > trakCollection;
  iEvent.getByToken(TrakCollection_,trakCollection);

  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerColl;
  iEvent.getByToken(TriggerCollection_,triggerColl);

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( triggerResults_Label , triggerResults_handle);


  const edm::TriggerNames & names = iEvent.triggerNames( *triggerResults_handle );

  uint ncombo = 0;
  float DiMuonMassMax_ = DiMuonMassCuts_[1];
  float DiMuonMassMin_ = DiMuonMassCuts_[0];
  float TrakTrakMassMax_ = TrakTrakMassCuts_[1];
  float TrakTrakMassMin_ = TrakTrakMassCuts_[0];
  float DiMuonDiTrakMassMax_ = DiMuonDiTrakMassCuts_[1];
  float DiMuonDiTrakMassMin_ = DiMuonDiTrakMassCuts_[0];


  pat::TriggerObjectStandAloneCollection filteredColl, matchedColl;
  std::vector< pat::PackedCandidate> filteredTracks;
  std::vector < UInt_t > filterResults;

  //Get the filtered trigger objects
  for ( size_t iTrigObj = 0; iTrigObj < triggerColl->size(); ++iTrigObj ) {

    pat::TriggerObjectStandAlone unPackedTrigger( triggerColl->at( iTrigObj ) );

    unPackedTrigger.unpackPathNames( names );
    unPackedTrigger.unpackFilterLabels(iEvent,*triggerResults_handle);

    bool filtered = false;
    UInt_t thisFilter = 0;

    for (size_t i = 0; i < HLTFilters_.size(); i++)
      if(unPackedTrigger.hasFilterLabel(HLTFilters_[i]))
        {
          thisFilter += (1<<i);
          filtered = true;
        }

    if(filtered)
    {
      filteredColl.push_back(unPackedTrigger);
      filterResults.push_back(thisFilter);
    }
  }

  //Matching

  for (std::vector<pat::PackedCandidate>::const_iterator trak = trakCollection->begin(), trakend=trakCollection->end(); trak!= trakend; ++trak)
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
          {
            filteredTracks.push_back(*trak);
            matchedColl.push_back(*trigger);
          }

        matched = true;
      }
    }
  }


// Note: Dimuon cand are sorted by decreasing vertex probability then first is associated with "best" dimuon
  for (pat::CompositeCandidateCollection::const_iterator dimuonCand = dimuon->begin(); dimuonCand != dimuon->end(); ++dimuonCand){
     if ( dimuonCand->mass() < DiMuonMassMax_  && dimuonCand->mass() > DiMuonMassMin_ ) {
       const pat::Muon *pmuon = dynamic_cast<const pat::Muon*>(dimuonCand->daughter("muonP"));
       const pat::Muon *nmuon = dynamic_cast<const pat::Muon*>(dimuonCand->daughter("muonN"));

// loop on track candidates, make DiMuonT candidate, positive charge
       for (size_t i = 0; i < filteredTracks.size(); i++)
       {
         auto posTrack = filteredTracks[i];
         if(posTrack.charge()<=0) continue;
         if(posTrack.pt()<0.5) continue;
	       if(fabs(posTrack.pdgId())!=211) continue;
	       if(!(posTrack.trackHighPurity())) continue;

         if ( IsTheSame(posTrack,*pmuon) || IsTheSame(posTrack,*nmuon)) continue;

// loop over second track candidate, negative charge
         for (size_t j = 0; j < filteredTracks.size(); j++)
         {

           if (i == j) continue;
           auto negTrack = filteredTracks[j];

           if(negTrack.charge()>=0) continue;
           if(negTrack.pt()<0.5) continue;
  	       if(fabs(negTrack.pdgId())!=211) continue;
  	       if(!(negTrack.trackHighPurity())) continue;

           if ( IsTheSame(negTrack,*pmuon) || IsTheSame(negTrack,*nmuon)) continue;

           pat::CompositeCandidate TTCand = makeTTCandidate(posTrack, negTrack);
           pat::CompositeCandidate TTTrigger = makeTTTriggerCandidate(matchedColl[i],matchedColl[j]);

           if ( TTCand.mass() < TrakTrakMassMax_ && TTCand.mass() > TrakTrakMassMin_ ) {

           pat::CompositeCandidate DiMuonTTCand = makeDiMuonTTCandidate(*dimuonCand, *&TTCand);

           const pat::CompositeCandidate* dimuonTriggerCand = dynamic_cast <const pat::CompositeCandidate *> (dimuonCand->daughter("mumuTrigger"));// ->userData<CompositeCandidate>()
           pat::CompositeCandidate DiMuonTTTriggerCand = makeDiMuonTTCandidate(*dimuonTriggerCand, *&TTTrigger);

           DiMuonTTCand.addDaughter(DiMuonTTTriggerCand,"dimuonTTTrigger");
           DiMuonTTCand.addUserInt("trakMatchP",isTriggerMatched(matchedColl[i]));
           DiMuonTTCand.addUserInt("trakMatchN",isTriggerMatched(matchedColl[j]));

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
  iEvent.put(std::move(DiMuonTTCandColl),product_name_);
  nevents++;
}

void DiMuonDiTrakProducerHLT::endJob(){
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

bool DiMuonDiTrakProducerHLT::IsTheSame(const pat::PackedCandidate& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

UInt_t DiMuonDiTrakProducerHLT::isTriggerMatched(const pat::TriggerObjectStandAlone& t) {

  UInt_t matched = 0;

  for (unsigned int iTr = 0; iTr<HLTFilters_.size(); iTr++ )
    if(t.hasFilterLabel(HLTFilters_[iTr]))  matched += (1<<iTr);

  return matched;
}

const pat::CompositeCandidate DiMuonDiTrakProducerHLT::makeDiMuonTTCandidate(
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


const pat::CompositeCandidate DiMuonDiTrakProducerHLT::makeTTCandidate(
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

float DiMuonDiTrakProducerHLT::DeltaR(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2)
{
   float p1 = t1.phi();
   float p2 = t2.phi();
   float e1 = t1.eta();
   float e2 = t2.eta();
   auto dp=std::abs(p1-p2); if (dp>float(M_PI)) dp-=float(2*M_PI);

   return sqrt((e1-e2)*(e1-e2) + dp*dp);
}

bool DiMuonDiTrakProducerHLT::MatchByDRDPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2)
{
  return (fabs(t1.pt()-t2.pt())/t2.pt()<MaxDeltaRPt_[1] &&
	DeltaR(t1,t2) < MaxDeltaRPt_[0]);
}

const pat::CompositeCandidate DiMuonDiTrakProducerHLT::makeTTTriggerCandidate(
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


reco::Candidate::LorentzVector DiMuonDiTrakProducerHLT::convertVector(const math::XYZTLorentzVectorF& v){

  return reco::Candidate::LorentzVector(v.x(),v.y(), v.z(), v.t());
}
//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonDiTrakProducerHLT);
