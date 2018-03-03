// -*- C++ -*-
//
// Package:    DiMuonRootupler
// Class:      DiMuonRootupler
//
// Description: Dimuon(mu+ mu-)  producer
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

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"

//
// class declaration
//

class DiMuonRootupler:public edm::EDAnalyzer {
      public:
	explicit DiMuonRootupler(const edm::ParameterSet &);
	~DiMuonRootupler() override;

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
	edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuon_Label;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;

  bool OnlyBest_;
  std::vector<std::string>  HLTs_;

	UInt_t    run;
	ULong64_t event;
  UInt_t    lumiblock;
  UInt_t    ndimuon;
  UInt_t    trigger;
  Int_t     charge;

	TLorentzVector dimuon_p4;
	TLorentzVector muonP_p4;
	TLorentzVector muonN_p4;

  Float_t MassErr;
  Float_t vProb;
  Float_t DCA;
  Float_t ctauPV;
  Float_t ctauErrPV;
  Float_t cosAlpha;
  Float_t lxyPV;
  Float_t lxyErrPV;

  Bool_t isBest;

	UInt_t numPrimaryVertices;

	TTree *dimuon_tree;

};

//
// constructors and destructor
//

DiMuonRootupler::DiMuonRootupler(const edm::ParameterSet & iConfig):
dimuon_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("dimuons"))),
primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"))),
triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
HLTs_(iConfig.getParameter<std::vector<std::string>>("HLTs"))
{
  edm::Service < TFileService > fs;
  dimuon_tree = fs->make < TTree > ("dimuonTree", "Tree of DiMuon");

  dimuon_tree->Branch("run",      &run,      "run/i");
  dimuon_tree->Branch("event",    &event,    "event/l");
  dimuon_tree->Branch("lumiblock",&lumiblock,"lumiblock/i");

  dimuon_tree->Branch("ndimuon",    &ndimuon,    "ndimuon/i");
  dimuon_tree->Branch("trigger",  &trigger,  "trigger/i");
  dimuon_tree->Branch("charge",   &charge,   "charge/I");

  dimuon_tree->Branch("isBest",   &isBest,   "isBest/O");

  dimuon_tree->Branch("dimuon_p4", "TLorentzVector", &dimuon_p4);
  dimuon_tree->Branch("muonP_p4",  "TLorentzVector", &muonP_p4);
  dimuon_tree->Branch("muonN_p4",  "TLorentzVector", &muonN_p4);

  dimuon_tree->Branch("MassErr",   &MassErr,    "MassErr/F");
  dimuon_tree->Branch("vProb",     &vProb,      "vProb/F");
  dimuon_tree->Branch("DCA",       &DCA,        "DCA/F");
  dimuon_tree->Branch("ctauPV",    &ctauPV,     "ctauPV/F");
  dimuon_tree->Branch("ctauErrPV", &ctauErrPV,  "ctauErrPV/F");
  dimuon_tree->Branch("cosAlpha",  &cosAlpha,   "cosAlpha/F");
  dimuon_tree->Branch("lxy",       &lxyPV,      "lxy/F");
  dimuon_tree->Branch("lxyErrPV",    &lxyErrPV,      "lxyErr/F");

  dimuon_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/i");


}

DiMuonRootupler::~DiMuonRootupler() {}

//
// member functions
//

const reco::Candidate* DiMuonRootupler::GetAncestor(const reco::Candidate* p) {
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


UInt_t DiMuonRootupler::getTriggerBits(const edm::Event& iEvent ) {

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
void DiMuonRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  iEvent.getByToken(dimuon_Label,dimuons);

  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  run       = iEvent.id().run();
  event     = iEvent.id().event();
  lumiblock = iEvent.id().luminosityBlock();

  numPrimaryVertices = 0;
  if (primaryVertices_handle.isValid()) numPrimaryVertices = (int) primaryVertices_handle->size();
  trigger = getTriggerBits(iEvent);

  ndimuon  = 0;

  dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  muonP_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  muonN_p4.SetPtEtaPhiM(0.,0.,0.,0.);

  isBest = false;

  if ( dimuons.isValid() && !dimuons->empty()) {
    for ( pat::CompositeCandidateCollection::const_iterator dimuonCand = dimuons->begin(); dimuonCand != dimuons->end(); ++dimuonCand ) {

        dimuon_p4.SetPtEtaPhiM(dimuonCand->pt(),dimuonCand->eta(),dimuonCand->phi(),dimuonCand->mass());

        reco::Candidate::LorentzVector vP = dimuonCand->daughter("muonP")->p4();
        reco::Candidate::LorentzVector vM = dimuonCand->daughter("muonN")->p4();

        muonP_p4.SetPtEtaPhiM(vP.pt(),vP.eta(),vP.phi(),vP.mass());
        muonN_p4.SetPtEtaPhiM(vM.pt(),vM.eta(),vM.phi(),vM.mass());

        MassErr = -1.0;
        if (dimuonCand->hasUserFloat("MassErr"))
          MassErr = dimuonCand->userFloat("MassErr");
        vProb = dimuonCand->userFloat("vProb");

        DCA = -1.;
        if (dimuonCand->hasUserFloat("DCA"))
          DCA = dimuonCand->userFloat("DCA");

        ctauPV = dimuonCand->userFloat("catuPV");
        ctauErrPV = dimuonCand->userFloat("catuErrPV");
        lxyPV = dimuonCand->userFloat("lxy");
        lxyErrPV = dimuonCand->userFloat("lErrxy");

        cosAlpha = dimuonCand->userFloat("cosAlpha");

        charge = dimuonCand->charge();

        dimuon_tree->Fill();   // be aware, we are storing all combinations
        isBest = false;

        if (OnlyBest_) break;
    }
  }

}

// ------------ method called once each job just before starting event loop  ------------
void DiMuonRootupler::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void DiMuonRootupler::endJob() {}

// ------------ method called when starting to processes a run  ------------
void DiMuonRootupler::beginRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a run  ------------
void DiMuonRootupler::endRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when starting to processes a luminosity block  ------------
void DiMuonRootupler::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void DiMuonRootupler::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DiMuonRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonRootupler);
