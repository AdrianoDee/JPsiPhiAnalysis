#include "../interface/DiTrackHLTAnalyzer.h"


UInt_t DiTrackHLTAnalyzer::isTriggerMatched(const pat::CompositeCandidate *diTrig_Candidate) {
  const pat::TriggerObjectStandAlone* trig1 = dynamic_cast<const pat::TriggerObjectStandAlone*>(diTrig_Candidate->daughter("trakP"));
  const pat::TriggerObjectStandAlone* trig2 = dynamic_cast<const pat::TriggerObjectStandAlone*>(diTrig_Candidate->daughter("trakN"));
  UInt_t matched = 0;  // if no list is given, is not matched

  // if matched a given trigger, set the bit, in the same order as listed
  for (unsigned int iTr = 0; iTr<HLTFilters_.size(); iTr++ ) {

    if(std::find((trig1->filterLabels()).begin(),(trig1->filterLabels()).end(),HLTFilters_[iTr])!=(trig1->filterLabels()).end())
      if(std::find((trig2->filterLabels()).begin(),(trig2->filterLabels()).end(),HLTFilters_[iTr])!=(trig2->filterLabels()).end())
        matched += (1<<iTr);

  }

  return matched;
}

float DiTrackHLTAnalyzer::DeltaR(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2)
{
   float p1 = t1.phi();
   float p2 = t2.phi();
   float e1 = t1.eta();
   float e2 = t2.eta();
   auto dp=std::abs(p1-p2); if (dp>float(M_PI)) dp-=float(2*M_PI);

   return sqrt((e1-e2)*(e1-e2) + dp*dp);
}

bool DiTrackHLTAnalyzer::MatchByDRDPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2)
{
  return (fabs(t1.pt()-t2.pt())/t2.pt()<maxDPtRel &&
	DeltaR(t1,t2) < maxDPtRel);
}

DiTrackHLTAnalyzer::DiTrackHLTAnalyzer(const edm::ParameterSet& iConfig):
  TrakCollection_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("PFCandidates"))),
  TriggerCollection_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("TriggerInput"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"))),
  TrakTrakMassCuts_(iConfig.getParameter<std::vector<double>>("TrakTrakMassCuts")),
  MassTraks_(iConfig.getParameter<std::vector<double>>("MassTraks")),
  OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  TTTrigger_name_(iConfig.getParameter<std::string>("TTTrigger_name")),
  HLTFilters_(iConfig.getParameter<std::vector<std::string>>("HLTFilters"))
{

  edm::Service < TFileService > fs;
  ditrak_tree = fs->make < TTree > ("DiTrakDiTrigTree", "Tree of ditrakditrig");

  ditrak_tree->Branch("run",      &run,      "run/i");
  ditrak_tree->Branch("event",    &event,    "event/l");
  ditrak_tree->Branch("lumiblock",&lumiblock,"lumiblock/i");


  ditrak_tree->Branch("nditrak",    &nditrak,    "nditrak/i");
  ditrak_tree->Branch("ntraks",   &ntraks,   "ntraks/i");
  ditrak_tree->Branch("trigger",  &trigger,  "trigger/i");
  ditrak_tree->Branch("charge",   &charge,   "charge/I");

  ditrak_tree->Branch("tMatch",    &tMatch,      "tMatch/I");

  // ditrak_tree->Branch("ditrak_p4", "TLorentzVector", &ditrak_p4);
  // ditrak_tree->Branch("trakP_p4",  "TLorentzVector", &trakP_p4);
  // ditrak_tree->Branch("trakN_p4",  "TLorentzVector", &trakN_p4);

  ditrak_tree->Branch("ditrig_p4", "TLorentzVector", &ditrig_p4);
  ditrak_tree->Branch("trigP_p4",  "TLorentzVector", &trigP_p4);
  ditrak_tree->Branch("trigN_p4",  "TLorentzVector", &trigN_p4);

  // ditrak_tree->Branch("MassErr",   &MassErr,    "MassErr/F");
  // ditrak_tree->Branch("vProb",     &vProb,      "vProb/F");
  // ditrak_tree->Branch("DCA",       &DCA,        "DCA/F");
  // ditrak_tree->Branch("ppdlPV",    &ppdlPV,     "ppdlPV/F");
  // ditrak_tree->Branch("ppdlErrPV", &ppdlErrPV,  "ppdlErrPV/F");
  // ditrak_tree->Branch("ppdlBS",    &ppdlBS,     "ppdlBS/F");
  // ditrak_tree->Branch("ppdlErrBS", &ppdlErrBS,  "ppdlErrBS/F");
  // ditrak_tree->Branch("cosAlpha",  &cosAlpha,   "cosAlpha/F");
  // ditrak_tree->Branch("lxyPV",     &lxyPV,      "lxyPV/F");
  // ditrak_tree->Branch("lxyBS",     &lxyBS,      "lxyBS/F");

  ditrak_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/i");


  candidates = 0;
  nevents = 0;
  ndimuon = 0;
  nreco = 0;
  maxDeltaR = 0.1;
  maxDPtRel = 10.0;
}

UInt_t DiTrackHLTAnalyzer::getTriggerBits(const edm::Event& iEvent ) {

  UInt_t trigger = 0;

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( triggerResults_Label , triggerResults_handle);

  if (triggerResults_handle.isValid()) {
     const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
     unsigned int NTRIGGERS = HLTs_.size();

     for (unsigned int i = 0; i < NTRIGGERS; i++) {
        for (int version = 1; version < 20; version++) {
           std::stringstream ss;
           ss << HLTs_[i] << "_v" << version;
           unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label());
           if (bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
              trigger += (1<<i);
              break;
           }
        }
     }
   } else std::cout << "*** NO triggerResults found " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;

   return trigger;
}

void DiTrackHLTAnalyzer::analyze(edm::Event& iEvent, const edm::EventSetup& iSetup){

  run       = iEvent.id().run();
  event     = iEvent.id().event();
  lumiblock = iEvent.id().luminosityBlock();

  // std::unique_ptr<pat::CompositeCandidateCollection> DiTrackColl(new pat::CompositeCandidateCollection);
  // std::unique_ptr<pat::TriggerObjectStandAloneCollection> DiTriggColl(new pat::CompositeCandidateCollection);

  edm::Handle<std::vector<pat::PackedCandidate> > trakColl;
  iEvent.getByToken(TrakCollection_,trakColl);

  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerColl;
  iEvent.getByToken(TriggerCollection_,triggerColl);

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( triggerResults_Label , triggerResults_handle);

  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  numPrimaryVertices = 0;
  if (primaryVertices_handle.isValid()) numPrimaryVertices = (int) primaryVertices_handle->size();
  trigger = getTriggerBits(iEvent);

  nditrak  = 0;
  ntraks = 0;

  ditrak_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  trakP_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  trakN_p4.SetPtEtaPhiM(0.,0.,0.,0.);

  ditrig_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  trigP_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  trigN_p4.SetPtEtaPhiM(0.,0.,0.,0.);

  uint ncombo = 0;

  float TrakTrakMassMax_ = TrakTrakMassCuts_[1];
  float TrakTrakMassMin_ = TrakTrakMassCuts_[0];

  pat::TriggerObjectStandAloneCollection filteredColl, matchedColl;
  std::vector< pat::PackedCandidate> filteredTracks;

  //Filtering

  for ( size_t iTrigObj = 0; iTrigObj < triggerColl->size(); ++iTrigObj ) {

    pat::TriggerObjectStandAlone unPackedTrigger( triggerColl->at( iTrigObj ) );

    const edm::TriggerNames & names = iEvent.triggerNames( *triggerResults_handle );

    unPackedTrigger.unpackPathNames( names );
    unPackedTrigger.unpackFilterLabels(iEvent,*triggerResults_handle);

    bool filtered = false;

    for (size_t i = 0; i < HLTFilters_.size(); i++)
      if(unPackedTrigger.hasFilterLabel(HLTFilters_[i]))
        filtered = true;

    if(filtered)
      filteredColl.push_back(unPackedTrigger);

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
          {
            filteredTracks.push_back(*trak);
            matchedColl.push_back(*trigger);
          }

        matched = true;
      }
    }
  }

  // std::cout << matchedColl.size() << " vs " << filteredTracks.size() << std::endl;
  // for (std::vector<pat::PackedCandidate>::const_iterator posTrack = filteredTracks.begin(), trakend=filteredTracks.end(); posTrack!= trakend; ++posTrack)
  for (size_t i = 0; i < filteredTracks.size(); i++)
  {
           auto posTrack = filteredTracks[i];
           if(posTrack.charge() <= 0 ) continue;
           if(posTrack.pt()<0.5) continue;
  	       if(fabs(posTrack.pdgId())!=211) continue;
  	       if(!(posTrack.trackHighPurity())) continue;



  // loop over second track candidate, negative charge
           // for (std::vector<pat::PackedCandidate>::const_iterator negTrack = trak->begin(); negTrack!= trakend; ++negTrack){
           for (size_t j = 0; j < filteredTracks.size(); j++) {

             if (i == j) continue;

             auto negTrack = filteredTracks[j];
             if(negTrack.charge() >= 0 ) continue;
             if(negTrack.pt()<0.5) continue;
    	       if(fabs(negTrack.pdgId())!=211) continue;
    	       if(!(negTrack.trackHighPurity())) continue;

             pat::CompositeCandidate TTCand = makeTTCandidate(posTrack,negTrack);
             pat::CompositeCandidate TTTrigger = makeTTTriggerCandidate(matchedColl[i],matchedColl[j]);

             if ( TTCand.mass() < TrakTrakMassMax_ && TTCand.mass() > TrakTrakMassMin_ ) {

               ditrak_p4.SetPtEtaPhiM(TTCand.pt(),TTCand.eta(),TTCand.phi(),TTCand.mass());

               charge = TTTrigger.charge();

               ditrig_p4.SetPtEtaPhiM(TTTrigger.pt(),TTTrigger.eta(),TTTrigger.phi(),TTTrigger.mass());

               reco::Candidate::LorentzVector vP = TTTrigger.daughter("trigP")->p4();
               reco::Candidate::LorentzVector vM = TTTrigger.daughter("trigN")->p4();

               trigP_p4.SetPtEtaPhiM(vP.pt(),vP.eta(),vP.phi(),vP.mass());
               trigN_p4.SetPtEtaPhiM(vM.pt(),vM.eta(),vM.phi(),vM.mass());

               tMatch = DiTrackHLTAnalyzer::isTriggerMatched(&TTTrigger);

               ditrak_tree->Fill();

               // DiTriggColl->push_back(TTTrigger);

             }
           } // loop over second track
         }   // loop on track candidates

  candidates = DiTriggColl->size();

  // if ( ncombo != DiTrackColl->size() ) std::cout <<"ncombo ("<<ncombo<< ") != DiMuonTT ("<<DiTrackColl->size()<<")"<< std::endl;
  // if ( ncombo > 0 ) nreco++;
  // iEvent.put(std::move(DiTrackColl),TTCandidate_name_);

}

void DiTrackHLTAnalyzer::endJob(){
  std::cout << "###########################" << std::endl;
  std::cout << "DiTrak(DiTrig) Candidate Analyzer report:" << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << candidates << " DiTrak (DiTrig) candidates." << std::endl;
  std::cout << "###########################" << std::endl;
}

void DiTrackHLTAnalyzer::beginRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a run  ------------
void DiTrackHLTAnalyzer::endRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when starting to processes a luminosity block  ------------
void DiTrackHLTAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void DiTrackHLTAnalyzer::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}


const pat::CompositeCandidate DiTrackHLTAnalyzer::makeTTCandidate(
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

const pat::CompositeCandidate DiTrackHLTAnalyzer::makeTTTriggerCandidate(
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


reco::Candidate::LorentzVector DiTrackHLTAnalyzer::convertVector(const math::XYZTLorentzVectorF& v){

  return reco::Candidate::LorentzVector(v.x(),v.y(), v.z(), v.t());
}

void DiTrackHLTAnalyzer::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiTrackHLTAnalyzer);
