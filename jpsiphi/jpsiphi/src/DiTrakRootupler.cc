// -*- C++ -*-
//
// Package:    DiTrakRootupler
// Class:      DiTrakRootupler
//
// Description: Ditrak(trk+ trk-)  producer
//
// Author:  Adriano Di Florio

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

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"

//
// class declaration
//

class DiTrakRootupler:public edm::EDAnalyzer {
      public:
	explicit DiTrakRootupler(const edm::ParameterSet &);
	~DiTrakRootupler() override;

	static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

      private:
        UInt_t getTriggerBits(const edm::Event &);
        bool   isAncestor(const reco::Candidate *, const reco::Candidate *);
        const  reco::Candidate* GetAncestor(const reco::Candidate *);

	void beginJob() override;
	void analyze(const edm::Event &, const edm::EventSetup &) override;
	void endJob() override;

	void beginRun(edm::Run const &, edm::EventSetup const &) override;
	void endRun(edm::Run const &, edm::EventSetup const &) override;
	void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) override;
	void endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) override;

	// ----------member data ---------------------------
	std::string file_name;
	edm::EDGetTokenT<pat::CompositeCandidateCollection> ditrak_Label;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggers_;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;

  bool addTrigger_;
  bool OnlyBest_;
  std::vector<std::string>  HLTs_;
  std::vector<std::string>  HLTFilters_;

	UInt_t    run;
	ULong64_t event;
  UInt_t    lumiblock;
  UInt_t    nditrak;
  UInt_t    trigger;
  Int_t     charge;

	TLorentzVector ditrak_p4;
	TLorentzVector trakP_p4;
	TLorentzVector trakN_p4;

  std::vector < Float_t > trigs_pt;
  std::vector < Float_t > trigs_eta;
  std::vector < Float_t > trigs_phi;
  std::vector < Float_t > trigs_m;
  std::vector < UInt_t > trigs_filters;

  Bool_t isBest;

  Float_t MassErr;
  Float_t vProb;
  Float_t DCA;
  Float_t ctauPV;
  Float_t ctauErrPV;
  Float_t cosAlpha;
  Float_t lxyPV;
  Float_t lxyErrPV;

	UInt_t numPrimaryVertices;

	TTree *ditrak_tree;

};

//
// constructors and destructor
//

DiTrakRootupler::DiTrakRootupler(const edm::ParameterSet & iConfig):
ditrak_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("ditraks"))),
triggers_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("TriggerInput"))),
primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"))),
triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
addTrigger_(iConfig.getParameter<bool>("AddTriggers")),
OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
HLTs_(iConfig.getParameter<std::vector<std::string>>("HLTs")),
HLTFilters_(iConfig.getParameter<std::vector<std::string>>("Filters"))
{
  edm::Service < TFileService > fs;
  ditrak_tree = fs->make < TTree > ("ditrakTree", "Tree of DiTrak");

  ditrak_tree->Branch("run",      &run,      "run/i");
  ditrak_tree->Branch("event",    &event,    "event/l");
  ditrak_tree->Branch("lumiblock",&lumiblock,"lumiblock/i");

  ditrak_tree->Branch("nditrak",  &nditrak,    "nditrak/i");
  ditrak_tree->Branch("trigger",  &trigger,  "trigger/i");
  ditrak_tree->Branch("charge",   &charge,   "charge/I");

  ditrak_tree->Branch("isBest",   &isBest,   "isBest/O");

  if(addTrigger_)
  {
    ditrak_tree->Branch("trigs_pt",   &trigs_pt);
    ditrak_tree->Branch("trigs_eta",   &trigs_eta);
    ditrak_tree->Branch("trigs_phi",   &trigs_phi);
    ditrak_tree->Branch("trigs_m",   &trigs_m);
    ditrak_tree->Branch("trigs_filters", &trigs_filters);
  }
  ditrak_tree->Branch("ditrak_p4", "TLorentzVector", &ditrak_p4);
  ditrak_tree->Branch("trakP_p4",  "TLorentzVector", &trakP_p4);
  ditrak_tree->Branch("trakN_p4",  "TLorentzVector", &trakN_p4);

  ditrak_tree->Branch("MassErr",   &MassErr,    "MassErr/F");
  ditrak_tree->Branch("vProb",     &vProb,      "vProb/F");
  ditrak_tree->Branch("DCA",       &DCA,        "DCA/F");
  ditrak_tree->Branch("ctauPV",    &ctauPV,     "ctauPV/F");
  ditrak_tree->Branch("ctauErrPV", &ctauErrPV,  "ctauErrPV/F");
  ditrak_tree->Branch("cosAlpha",  &cosAlpha,   "cosAlpha/F");
  ditrak_tree->Branch("lxy",       &lxyPV,      "lxy/F");
  ditrak_tree->Branch("lxyErrPV",    &lxyErrPV,      "lxyErr/F");

  ditrak_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/i");


}

DiTrakRootupler::~DiTrakRootupler() {}

//
// member functions
//

const reco::Candidate* DiTrakRootupler::GetAncestor(const reco::Candidate* p) {
   if (p->numberOfMothers()) {
      if  ((p->mother(0))->pdgId() == p->pdgId()) return GetAncestor(p->mother(0));
      else return p->mother(0);
   }
   return p;
}


/* Grab Trigger information. Save it in variable trigger, trigger is an uint between 0 and 256, in binary it is:
   (pass 2)(pass 1)(pass 0)
   ex. 7 = pass 0, 1 and 2
   ex. 6 = pass 1, 2
   ex. 1 = pass 0
*/


UInt_t DiTrakRootupler::getTriggerBits(const edm::Event& iEvent ) {

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
void DiTrakRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> trigs;
  iEvent.getByToken(triggers_,trigs);

  edm::Handle<pat::CompositeCandidateCollection> ditraks;
  iEvent.getByToken(ditrak_Label,ditraks);

  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( triggerResults_Label , triggerResults_handle);

  run       = iEvent.id().run();
  event     = iEvent.id().event();
  lumiblock = iEvent.id().luminosityBlock();

  numPrimaryVertices = 0;
  if (primaryVertices_handle.isValid()) numPrimaryVertices = (int) primaryVertices_handle->size();
  trigger = getTriggerBits(iEvent);

  const edm::TriggerNames & names = iEvent.triggerNames( *triggerResults_handle );

  nditrak  = ditraks->size();

  ditrak_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  trakP_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  trakN_p4.SetPtEtaPhiM(0.,0.,0.,0.);

  isBest = true;

  trigs_filters.clear();
  trigs_pt.clear();
  trigs_eta.clear();
  trigs_phi.clear();
  trigs_m.clear();

  if(addTrigger_)
  {
    for ( size_t iTrigObj = 0; iTrigObj < trigs->size(); ++iTrigObj ) {


      pat::TriggerObjectStandAlone unPackedTrigger( trigs->at( iTrigObj ) );

      if(unPackedTrigger.charge()==0) continue;

      unPackedTrigger.unpackPathNames( names );
      unPackedTrigger.unpackFilterLabels(iEvent,*triggerResults_handle);

      bool filtered = false;
      UInt_t thisFilter = 0;

      for (size_t i = 0; i < HLTFilters_.size(); i++)
      {
        if(unPackedTrigger.hasFilterLabel(HLTFilters_[i]))
          {
            thisFilter += (1<<i);
            filtered = true;
          }
      }

      if(filtered)
      {
        trigs_filters.push_back(thisFilter);
        trigs_pt.push_back(unPackedTrigger.pt());
        trigs_eta.push_back(unPackedTrigger.eta());
        trigs_phi.push_back(unPackedTrigger.phi());
        trigs_m.push_back(unPackedTrigger.mass());
      }


    }
}
  if ( ditraks.isValid() && !ditraks->empty()) {
    for ( pat::CompositeCandidateCollection::const_iterator ditrakCand = ditraks->begin(); ditrakCand != ditraks->end(); ++ditrakCand ) {

        ditrak_p4.SetPtEtaPhiM(ditrakCand->pt(),ditrakCand->eta(),ditrakCand->phi(),ditrakCand->mass());

        reco::Candidate::LorentzVector vP = ditrakCand->daughter("trakP")->p4();
        reco::Candidate::LorentzVector vM = ditrakCand->daughter("trakN")->p4();

        trakP_p4.SetPtEtaPhiM(vP.pt(),vP.eta(),vP.phi(),vP.mass());
        trakN_p4.SetPtEtaPhiM(vM.pt(),vM.eta(),vM.phi(),vM.mass());

        MassErr = -1.0;
        if (ditrakCand->hasUserFloat("MassErr"))
          MassErr = ditrakCand->userFloat("MassErr");
        vProb = ditrakCand->userFloat("vProb");

        DCA = -1.;
        if (ditrakCand->hasUserFloat("DCA"))
          DCA = ditrakCand->userFloat("DCA");

        ctauPV = ditrakCand->userFloat("ctauPV");
        ctauErrPV = ditrakCand->userFloat("ctauErrPV");
        lxyPV = ditrakCand->userFloat("lxy");
        lxyErrPV = ditrakCand->userFloat("lErrxy");

        cosAlpha = ditrakCand->userFloat("cosAlpha");

        charge = ditrakCand->charge();

        ditrak_tree->Fill();   // be aware, we are storing all combinations

        if (OnlyBest_) break;
        isBest = false;
    }
  }


}

// ------------ method called once each job just before starting event loop  ------------
void DiTrakRootupler::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void DiTrakRootupler::endJob() {}

// ------------ method called when starting to processes a run  ------------
void DiTrakRootupler::beginRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a run  ------------
void DiTrakRootupler::endRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when starting to processes a luminosity block  ------------
void DiTrakRootupler::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void DiTrakRootupler::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DiTrakRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiTrakRootupler);
