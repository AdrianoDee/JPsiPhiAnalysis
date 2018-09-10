// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "FWCore/Common/interface/TriggerNames.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"

//
// class declaration
//

class DiTrakHLT:public edm::EDAnalyzer {
      public:
	explicit DiTrakHLT(const edm::ParameterSet &);
	~DiTrakHLT() override;

	static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

      private:
        UInt_t getTriggerBits(const edm::Event& iEvent, const edm::Handle< edm::TriggerResults >& triggerResults_handle);
        bool   isAncestor(const reco::Candidate *, const reco::Candidate *);
        const  reco::Candidate* GetAncestor(const reco::Candidate *);
        UInt_t isTriggerMatched(const pat::CompositeCandidate *diTrig_Candidate);

	void beginJob() override;
	void analyze(const edm::Event &, const edm::EventSetup &) override;
	void endJob() override;

	void beginRun(edm::Run const &, edm::EventSetup const &) override;
	void endRun(edm::Run const &, edm::EventSetup const &) override;
	void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) override;
	void endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) override;

	// ----------member data ---------------------------
	std::string file_name;
	// edm::EDGetTokenT<pat::CompositeCandidateCollection> diTrak_label;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate>> TrakCollection_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> TriggerCollection_;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  std::vector<double> ditrakMassCuts_;
  std::vector<double> MassTraks_;
	bool isMC_;
  bool OnlyBest_;
  std::vector<std::string>  HLTs_;
  std::vector<std::string>  HLTFilters_;

  reco::Candidate::LorentzVector convertVector(const math::XYZTLorentzVectorF& v);
  bool IsTheSame(const pat::PackedCandidate& t1, const pat::PackedCandidate& t2);

  const pat::CompositeCandidate makeTTTriggerCandidate(const pat::TriggerObjectStandAlone& t1,
						    const pat::TriggerObjectStandAlone& t2);
  const pat::CompositeCandidate makeTTCandidate(const pat::PackedCandidate& t1,
                                                const pat::PackedCandidate& t2);

  bool MatchByDRDPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2);
  float DeltaR(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2);

  int candidates;
  int nevents;
  int ndimuon;
  int nreco;
  float maxDeltaR;
  float maxDPtRel;

	UInt_t    run;
	ULong64_t event;
  UInt_t    lumiblock;

  UInt_t    trigger;
  UInt_t    tMatchOne,tMatchTwo;

  UInt_t charge;

	TLorentzVector ditrak_p4;
	TLorentzVector trakP_p4;
	TLorentzVector trakN_p4;

  TLorentzVector ditrig_p4;
  TLorentzVector trigP_p4;
  TLorentzVector trigN_p4;

	UInt_t numPrimaryVertices;

	TTree *ditrak_tree;

  UInt_t nditrak, ntraks;


};

UInt_t DiTrakHLT::isTriggerMatched(const pat::CompositeCandidate *diTrig_Candidate) {
  const pat::TriggerObjectStandAlone* trig1 = dynamic_cast<const pat::TriggerObjectStandAlone*>(diTrig_Candidate->daughter("trigP"));
  const pat::TriggerObjectStandAlone* trig2 = dynamic_cast<const pat::TriggerObjectStandAlone*>(diTrig_Candidate->daughter("trigN"));
  UInt_t matched = 0;  // if no list is given, is not matched

  // if matched a given trigger, set the bit, in the same order as listed
  for (unsigned int iTr = 0; iTr<HLTFilters_.size(); iTr++ ) {

    if(std::find((trig1->filterLabels()).begin(),(trig1->filterLabels()).end(),HLTFilters_[iTr])!=(trig1->filterLabels()).end())
      if(std::find((trig2->filterLabels()).begin(),(trig2->filterLabels()).end(),HLTFilters_[iTr])!=(trig2->filterLabels()).end())
        matched += (1<<iTr);

  }

  return matched;
}

float DiTrakHLT::DeltaR(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2)
{
   float p1 = t1.phi();
   float p2 = t2.phi();
   float e1 = t1.eta();
   float e2 = t2.eta();
   auto dp=std::abs(p1-p2); if (dp>float(M_PI)) dp-=float(2*M_PI);

   return sqrt((e1-e2)*(e1-e2) + dp*dp);
}

bool DiTrakHLT::MatchByDRDPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2)
{
  return (fabs(t1.pt()-t2.pt())/t2.pt()<maxDPtRel &&
	DeltaR(t1,t2) < maxDeltaR);
}

const pat::CompositeCandidate DiTrakHLT::makeTTCandidate(
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

const pat::CompositeCandidate DiTrakHLT::makeTTTriggerCandidate(
                                          const pat::TriggerObjectStandAlone& trakP,
                                          const pat::TriggerObjectStandAlone& trakN
                                         ){

  pat::CompositeCandidate TTCand;
  TTCand.addDaughter(trakP,"trigP");
  TTCand.addDaughter(trakN,"trigN");
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


//
// constructors and destructor
//

DiTrakHLT::DiTrakHLT(const edm::ParameterSet & iConfig):
TrakCollection_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("PFCandidates"))),
TriggerCollection_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("TriggerInput"))),
primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"))),
triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
ditrakMassCuts_(iConfig.getParameter<std::vector<double>>("TrakTrakMassCuts")),
MassTraks_(iConfig.getParameter<std::vector<double>>("MassTraks")),
isMC_(iConfig.getParameter<bool>("isMC")),
OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
HLTs_(iConfig.getParameter<std::vector<std::string>>("HLTs")),
HLTFilters_(iConfig.getParameter<std::vector<std::string>>("Filters"))
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

  ditrak_tree->Branch("tMatchOne",    &tMatchOne,      "tMatchOne/I");
  ditrak_tree->Branch("tMatchTwo",    &tMatchTwo,      "tMatchTwo/I");
  ditrak_tree->Branch("ditrak_p4", "TLorentzVector", &ditrak_p4);
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
  maxDeltaR = 0.01;
  maxDPtRel = 2.0;

}

DiTrakHLT::~DiTrakHLT() {}

/* Grab Trigger information. Save it in variable trigger, trigger is an uint between 0 and 256, in binary it is:
   (pass 2)(pass 1)(pass 0)
   ex. 7 = pass 0, 1 and 2
   ex. 6 = pass 1, 2
   ex. 1 = pass 0
*/

UInt_t DiTrakHLT::getTriggerBits(const edm::Event& iEvent, const edm::Handle< edm::TriggerResults >& triggerResults_handle) {

  UInt_t trigger = 0;
  const edm::TriggerNames & names = iEvent.triggerNames( *triggerResults_handle );

     unsigned int NTRIGGERS = HLTs_.size();

     for (unsigned int i = 0; i < NTRIGGERS; i++) {
        for (int version = 1; version < 20; version++) {
           std::stringstream ss;
           ss << HLTs_[i] << "_v" << version;
           unsigned int bit = names.triggerIndex(edm::InputTag(ss.str()).label());
           if (bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
              trigger += (1<<i);
              break;
           }
        }
     }

   return trigger;
}

// ------------ method called for each event  ------------
void DiTrakHLT::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  // edm::Handle<pat::CompositeCandidateCollection> ditraks;
  // iEvent.getByToken(diTrak_label,ditraks);

  // edm::Handle<pat::CompositeCandidateCollection> ditrigs;
  // iEvent.getByToken(diTrig_label,ditrigs);

  edm::Handle<std::vector<pat::PackedCandidate> > trakColl;
  iEvent.getByToken(TrakCollection_,trakColl);

  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerColl;
  iEvent.getByToken(TriggerCollection_,triggerColl);

  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  run       = iEvent.id().run();
  event     = iEvent.id().event();
  lumiblock = iEvent.id().luminosityBlock();

  numPrimaryVertices = 0;
  if (primaryVertices_handle.isValid()) numPrimaryVertices = (int) primaryVertices_handle->size();

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( triggerResults_Label , triggerResults_handle);

  const edm::TriggerNames & names = iEvent.triggerNames( *triggerResults_handle );

  trigger = 0;

  if (triggerResults_handle.isValid())
    trigger = getTriggerBits(iEvent,triggerResults_handle);
  else std::cout << "*** NO triggerResults found " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;


  nditrak  = 0;
  ntraks = 0;

  ditrak_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  trakP_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  trakN_p4.SetPtEtaPhiM(0.,0.,0.,0.);

  ditrig_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  trigP_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  trigN_p4.SetPtEtaPhiM(0.,0.,0.,0.);

  float TrakTrakMassMax_ = ditrakMassCuts_[1];
  float TrakTrakMassMin_ = ditrakMassCuts_[0];

  pat::TriggerObjectStandAloneCollection filteredColl, matchedColl;
  std::vector< pat::PackedCandidate> filteredTracks;
  std::vector < UInt_t > filterResults;

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

               tMatchOne = filterResults[i];
               tMatchTwo = filterResults[j];
               ditrak_tree->Fill();

               candidates++;
               // DiTriggColl->push_back(TTTrigger);

             }
           } // loop over second track
         }   // loop on track candidates
}

// ------------ method called once each job just before starting event loop  ------------
void DiTrakHLT::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void DiTrakHLT::endJob() {}

// ------------ method called when starting to processes a run  ------------
void DiTrakHLT::beginRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a run  ------------
void DiTrakHLT::endRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when starting to processes a luminosity block  ------------
void DiTrakHLT::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void DiTrakHLT::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DiTrakHLT::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiTrakHLT);
