// -*- C++ -*-
//
// Package:    DiMuonDiTrakRootupler
// Class:      DiMuonDiTrakRootupler
//
// Description: DiMuonDiTrak  rootupler
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

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"

//
// class declaration
//

class DiMuonDiTrakRootupler:public edm::EDAnalyzer {
      public:
	explicit DiMuonDiTrakRootupler(const edm::ParameterSet &);
	~DiMuonDiTrakRootupler() override;

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
	edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuonditrk_Label;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggers_;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;

  bool OnlyBest_;
  std::vector<std::string>  HLTs_;
  std::vector<std::string>  HLTFilters_;

	UInt_t    run;
	ULong64_t event;
  UInt_t    lumiblock;
  UInt_t    ndimuonditrk;
  UInt_t    trigger;
  Int_t     charge;

	TLorentzVector dimuonditrak_p4;
  TLorentzVector ditrak_p4;
  TLorentzVector dimuon_p4;

  std::vector < TLorentzVector > trigs_p4;

  std::vector < Float_t > trigs_pt;
  std::vector < Float_t > trigs_eta;
  std::vector < Float_t > trigs_phi;
  std::vector < Float_t > trigs_m;
  std::vector < UInt_t > trigs_filters;

	TLorentzVector muonP_p4;
	TLorentzVector muonN_p4;
  TLorentzVector trakP_p4;
  TLorentzVector trakN_p4;

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

	TTree *dimuonditrk_tree;

};

//
// constructors and destructor
//

DiMuonDiTrakRootupler::DiMuonDiTrakRootupler(const edm::ParameterSet & iConfig):
dimuonditrk_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("dimuonditrks"))),
triggers_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("TriggerInput"))),
primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"))),
triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
HLTs_(iConfig.getParameter<std::vector<std::string>>("HLTs")),
HLTFilters_(iConfig.getParameter<std::vector<std::string>>("Filters"))
{
  edm::Service < TFileService > fs;
  dimuonditrk_tree = fs->make < TTree > ("dimuonditrkTree", "Tree of DiMuonDiTrak");

  dimuonditrk_tree->Branch("run",      &run,      "run/i");
  dimuonditrk_tree->Branch("event",    &event,    "event/l");
  dimuonditrk_tree->Branch("lumiblock",&lumiblock,"lumiblock/i");

  dimuonditrk_tree->Branch("ndimuonditrk",    &ndimuonditrk,    "ndimuonditrk/i");
  dimuonditrk_tree->Branch("trigger",  &trigger,  "trigger/i");
  dimuonditrk_tree->Branch("charge",   &charge,   "charge/I");

  dimuonditrk_tree->Branch("isBest",   &isBest,   "isBest/O");

  dimuonditrk_tree->Branch("dimuonditrak_p4", "TLorentzVector", &dimuonditrak_p4);

  dimuonditrk_tree->Branch("dimuon_p4", "TLorentzVector", &dimuon_p4);
  dimuonditrk_tree->Branch("muonP_p4",  "TLorentzVector", &muonP_p4);
  dimuonditrk_tree->Branch("muonN_p4",  "TLorentzVector", &muonN_p4);

  dimuonditrk_tree->Branch("trigs_p4", &trigs_p4);
  dimuonditrk_tree->Branch("trigs_pt",   &trigs_pt);
  dimuonditrk_tree->Branch("trigs_eta",   &trigs_eta);
  dimuonditrk_tree->Branch("trigs_phi",   &trigs_phi);
  dimuonditrk_tree->Branch("trigs_m",   &trigs_m);
  dimuonditrk_tree->Branch("trigs_filters", &trigs_filters);

  dimuonditrk_tree->Branch("ditrak_p4", "TLorentzVector", &ditrak_p4);
  dimuonditrk_tree->Branch("trakP_p4",  "TLorentzVector", &trakP_p4);
  dimuonditrk_tree->Branch("trakN_p4",  "TLorentzVector", &trakN_p4);

  dimuonditrk_tree->Branch("MassErr",   &MassErr,    "MassErr/F");
  dimuonditrk_tree->Branch("vProb",     &vProb,      "vProb/F");
  dimuonditrk_tree->Branch("DCA",       &DCA,        "DCA/F");
  dimuonditrk_tree->Branch("ctauPV",    &ctauPV,     "ctauPV/F");
  dimuonditrk_tree->Branch("ctauErrPV", &ctauErrPV,  "ctauErrPV/F");
  dimuonditrk_tree->Branch("cosAlpha",  &cosAlpha,   "cosAlpha/F");
  dimuonditrk_tree->Branch("lxy",       &lxyPV,      "lxy/F");
  dimuonditrk_tree->Branch("lxyErrPV",    &lxyErrPV,      "lxyErr/F");

  dimuonditrk_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/i");


}

DiMuonDiTrakRootupler::~DiMuonDiTrakRootupler() {}

//
// member functions
//

UInt_t DiMuonDiTrakRootupler::getTriggerBits(const edm::Event& iEvent ) {

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
void DiMuonDiTrakRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> trigs;
  iEvent.getByToken(triggers_,trigs);

  edm::Handle<pat::CompositeCandidateCollection> dimuonditrks;
  iEvent.getByToken(dimuonditrk_Label,dimuonditrks);

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

  ndimuonditrk  = dimuonditrks->size();

  dimuonditrak_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  muonP_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  muonN_p4.SetPtEtaPhiM(0.,0.,0.,0.);

  isBest = true;

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

    if(!filtered) continue;

    trigs_filters.push_back(thisFilter);
    trigs_pt.push_back(unPackedTrigger.pt());
    trigs_eta.push_back(unPackedTrigger.eta());
    trigs_phi.push_back(unPackedTrigger.phi());
    trigs_m.push_back(unPackedTrigger.mass());

  }

  if ( dimuonditrks.isValid() && !dimuonditrks->empty()) {
    for ( pat::CompositeCandidateCollection::const_iterator dimuonditrkCand = dimuonditrks->begin(); dimuonditrkCand != dimuonditrks->end(); ++dimuonditrkCand ) {

        dimuonditrak_p4.SetPtEtaPhiM(dimuonditrkCand->pt(),dimuonditrkCand->eta(),dimuonditrkCand->phi(),dimuonditrkCand->mass());

        const pat::CompositeCandidate *dimuon_cand = dynamic_cast <const pat::CompositeCandidate *>(dimuonditrkCand->daughter("dimuon"));
        const pat::CompositeCandidate *ditrak_cand = dynamic_cast <const pat::CompositeCandidate *>(dimuonditrkCand->daughter("ditrak"));

        dimuon_p4.SetPtEtaPhiM(dimuon_cand->pt(),dimuon_cand->eta(),dimuon_cand->phi(),dimuon_cand->mass());
        ditrak_p4.SetPtEtaPhiM(ditrak_cand->pt(),ditrak_cand->eta(),ditrak_cand->phi(),ditrak_cand->mass());

        reco::Candidate::LorentzVector vP = dimuon_cand->daughter("muonP")->p4();
        reco::Candidate::LorentzVector vM = dimuon_cand->daughter("muonN")->p4();

        muonP_p4.SetPtEtaPhiM(vP.pt(),vP.eta(),vP.phi(),vP.mass());
        muonN_p4.SetPtEtaPhiM(vM.pt(),vM.eta(),vM.phi(),vM.mass());

        reco::Candidate::LorentzVector tP = ditrak_cand->daughter("trakP")->p4();
        reco::Candidate::LorentzVector tM = ditrak_cand->daughter("trakN")->p4();

        trakP_p4.SetPtEtaPhiM(tP.pt(),tP.eta(),tP.phi(),tP.mass());
        trakN_p4.SetPtEtaPhiM(tM.pt(),tM.eta(),tM.phi(),tM.mass());

        MassErr = -1.0;
        if (dimuonditrkCand->hasUserFloat("MassErr"))
          MassErr = dimuonditrkCand->userFloat("MassErr");
        vProb = dimuonditrkCand->userFloat("vProb");

        DCA = -1.;
        if (dimuonditrkCand->hasUserFloat("DCA"))
          DCA = dimuonditrkCand->userFloat("DCA");

        ctauPV = dimuonditrkCand->userFloat("ctauPV");
        ctauErrPV = dimuonditrkCand->userFloat("ctauErrPV");
        lxyPV = dimuonditrkCand->userFloat("lxy");
        lxyErrPV = dimuonditrkCand->userFloat("lErrxy");

        cosAlpha = dimuonditrkCand->userFloat("cosAlpha");

        charge = dimuonditrkCand->charge();

        dimuonditrk_tree->Fill();
        isBest = false;
        if (OnlyBest_) break;
    }
  }


}

// ------------ method called once each job just before starting event loop  ------------
void DiMuonDiTrakRootupler::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void DiMuonDiTrakRootupler::endJob() {}

// ------------ method called when starting to processes a run  ------------
void DiMuonDiTrakRootupler::beginRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a run  ------------
void DiMuonDiTrakRootupler::endRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when starting to processes a luminosity block  ------------
void DiMuonDiTrakRootupler::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void DiMuonDiTrakRootupler::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DiMuonDiTrakRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonDiTrakRootupler);
