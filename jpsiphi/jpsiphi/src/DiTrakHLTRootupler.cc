// -*- C++ -*-
//
// Package:    DiTrakHLTRootupler
// Class:      DiTrakHLTRootupler
//
// Description: ditrak(mu+ mu-)  producer
//
// Author:  Adriano Di Florio
//    based on : Alberto Sanchez Hernandez Onia2MuMu code
//

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

class DiTrakHLTRootupler:public edm::EDAnalyzer {
      public:
	explicit DiTrakHLTRootupler(const edm::ParameterSet &);
	~DiTrakHLTRootupler() override;

	static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

      private:
        UInt_t getTriggerBits(const edm::Event &);
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
  edm::EDGetTokenT<pat::CompositeCandidateCollection> diTrig_label;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  std::vector<double> ditrakMassCuts_;
	bool isMC_;
  bool OnlyBest_;
  std::vector<std::string>  HLTs_;
  std::vector<std::string>  HLTFilters_;

	UInt_t    run;
	ULong64_t event;
  UInt_t    lumiblock;

  UInt_t    trigger;
  UInt_t    tMatch;

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

UInt_t DiTrakHLTRootupler::isTriggerMatched(const pat::CompositeCandidate *diTrig_Candidate) {
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

//
// constructors and destructor
//

DiTrakHLTRootupler::DiTrakHLTRootupler(const edm::ParameterSet & iConfig):
// diTrak_label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("ditraks"))),
diTrig_label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("ditrigs"))),
primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"))),
triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
ditrakMassCuts_(iConfig.getParameter<std::vector<double>>("TrakTrakMassCuts")),
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

}

DiTrakHLTRootupler::~DiTrakHLTRootupler() {}

/* Grab Trigger information. Save it in variable trigger, trigger is an uint between 0 and 256, in binary it is:
   (pass 2)(pass 1)(pass 0)
   ex. 7 = pass 0, 1 and 2
   ex. 6 = pass 1, 2
   ex. 1 = pass 0
*/

UInt_t DiTrakHLTRootupler::getTriggerBits(const edm::Event& iEvent ) {

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

// ------------ method called for each event  ------------
void DiTrakHLTRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  // edm::Handle<pat::CompositeCandidateCollection> ditraks;
  // iEvent.getByToken(diTrak_label,ditraks);

  edm::Handle<pat::CompositeCandidateCollection> ditrigs;
  iEvent.getByToken(diTrig_label,ditrigs);

  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  run       = iEvent.id().run();
  event     = iEvent.id().event();
  lumiblock = iEvent.id().luminosityBlock();

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

  float ditrakMassMax_ = ditrakMassCuts_[1];
  float ditrakMassMin_ = ditrakMassCuts_[0];

  bool already_stored = false;

  if ( ditrigs.isValid() && !ditrigs->empty()) {
  // if ( ditraks.isValid() && !ditraks->empty()) {
    // for ( pat::CompositeCandidateCollection::const_iterator ditrakCand = ditraks->begin(); ditrakCand != ditraks->end(); ++ditrakCand ) {
    for (size_t i = 0; i < (*ditrigs).size(); ++i) {
      std::cout<< "Valids"<<std::endl;
      // const pat::CompositeCandidate ditrakCand = (*ditraks)[i];
      const pat::CompositeCandidate ditrigCand = (*ditrigs)[i];
      pat::CompositeCandidate ditrakCand = ditrigCand.daughter("TrakTrak");
      
      if (ditrakCand.mass() > ditrakMassMin_ && ditrakCand.mass() < ditrakMassMax_ && ditrakCand.charge() == 0) {

        ditrak_p4.SetPtEtaPhiM(ditrakCand.pt(),ditrakCand.eta(),ditrakCand.phi(),ditrakCand.mass());

        reco::Candidate::LorentzVector vP = ditrakCand.daughter("trakP")->p4();
        reco::Candidate::LorentzVector vM = ditrakCand.daughter("trakN")->p4();

        trakP_p4.SetPtEtaPhiM(vP.pt(),vP.eta(),vP.phi(),vP.mass());
        trakN_p4.SetPtEtaPhiM(vM.pt(),vM.eta(),vM.phi(),vM.mass());

        charge = ditrakCand.charge();

        ditrig_p4.SetPtEtaPhiM(ditrigCand.pt(),ditrigCand.eta(),ditrigCand.phi(),ditrigCand.mass());

        vP = ditrigCand.daughter("trigP")->p4();
        vM = ditrigCand.daughter("trigN")->p4();

        trigP_p4.SetPtEtaPhiM(vP.pt(),vP.eta(),vP.phi(),vP.mass());
        trigN_p4.SetPtEtaPhiM(vM.pt(),vM.eta(),vM.phi(),vM.mass());

        tMatch = DiTrakHLTRootupler::isTriggerMatched(&ditrigCand);

        // MassErr = -1.0;
        // if (ditrakCand.hasUserFloat("MassErr"))    MassErr = ditrakCand.userFloat("MassErr");
        // vProb = ditrakCand.userFloat("vProb");
        // DCA = -1.;
        // if (ditrakCand.hasUserFloat("DCA"))  DCA = ditrakCand.userFloat("DCA");
        // ppdlPV = ditrakCand.userFloat("ppdlPV");
        // ppdlErrPV = ditrakCand.userFloat("ppdlErrPV");
        // ppdlBS = ditrakCand.userFloat("ppdlBS");
        // ppdlErrBS = ditrakCand.userFloat("ppdlErrBS");
        // cosAlpha = ditrakCand.userFloat("cosAlpha");
        // tMatch = ditrakCand.userInt("isTriggerMatched");
        // charge = ditrakCand.charge();
        // TVector3 pperp(ditrakCand.px(),ditrakCand.py(),0);
        // lxyPV = ppdlPV * pperp.Perp() / ditrakCand.mass();
        // lxyBS = ppdlBS * pperp.Perp() / ditrakCand.mass();
        nditrak++;
        if (OnlyBest_) break;
        else {
          ditrak_tree->Fill();   // be aware, we are storing all combinations
          already_stored = true;
        }
      }
    }
  } //..else {
    //std::cout << "DiTrakHLTRootupler: (" << run << "," << event << ") -> ";


  if ( !already_stored ) {  // we have to make sure, we are not double storing an combination
    if ( !isMC_ ) {
      if ( nditrak > 0 ) ditrak_tree->Fill();   // if not MC filter out
    } else ditrak_tree->Fill();
  }
}

// ------------ method called once each job just before starting event loop  ------------
void DiTrakHLTRootupler::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void DiTrakHLTRootupler::endJob() {}

// ------------ method called when starting to processes a run  ------------
void DiTrakHLTRootupler::beginRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a run  ------------
void DiTrakHLTRootupler::endRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when starting to processes a luminosity block  ------------
void DiTrakHLTRootupler::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void DiTrakHLTRootupler::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DiTrakHLTRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiTrakHLTRootupler);
