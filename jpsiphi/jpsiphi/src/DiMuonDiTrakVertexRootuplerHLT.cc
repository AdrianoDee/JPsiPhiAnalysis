// -*- C++ -*-
//
// Package:    DiMuonDiTrakVertexRootuplerHLT
// Class:      DiMuonDiTrakVertexRootuplerHLT
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

class DiMuonDiTrakVertexRootuplerHLT:public edm::EDAnalyzer {
      public:
	explicit DiMuonDiTrakVertexRootuplerHLT(const edm::ParameterSet &);
	~DiMuonDiTrakVertexRootuplerHLT() override;

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

  TLorentzVector dimuonditrkTrigger_p4;
  TLorentzVector dimuonTrigger_p4;
  TLorentzVector ditrakTrigger_p4;

	TLorentzVector muonP_p4;
	TLorentzVector muonN_p4;
  TLorentzVector trakP_p4;
  TLorentzVector trakN_p4;

  Int_t dimuonditrk_charge;

  Double_t dimuon_vProb, dimuon_vNChi2, dimuon_DCA, dimuon_ctauPV, dimuon_ctauErrPV, dimuon_cosAlpha;

  Double_t dimuonditrk_vProb, dimuonditrk_vChi2, dimuonditrk_cosAlpha, dimuonditrk_ctauPV, dimuonditrk_ctauErrPV;
  Double_t dimuonditrk_MassErr,dimuonditrk_lxy, dimuonditrk_lxyErr;

  Int_t noXCandidates;

  UInt_t muonN_tMatch, muonP_tMatch, trakN_tMatch, trakP_tMatch;

  Bool_t isBest;

	UInt_t numPrimaryVertices;

	TTree *dimuonditrk_tree;

};

//
// constructors and destructor
//

DiMuonDiTrakVertexRootuplerHLT::DiMuonDiTrakVertexRootuplerHLT(const edm::ParameterSet & iConfig):
dimuonditrk_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("dimuonditrks"))),
primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"))),
triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
HLTs_(iConfig.getParameter<std::vector<std::string>>("HLTs")),
HLTFilters_(iConfig.getParameter<std::vector<std::string>>("Filters"))
{
  edm::Service < TFileService > fs;
  dimuonditrk_tree = fs->make < TTree > ("dimuonditrkTree", "Tree of DiMuonDiTrak");

  dimuonditrk_tree->Branch("dimuonditrak_p4", "TLorentzVector", &dimuonditrak_p4);
  dimuonditrk_tree->Branch("dimuon_p4", "TLorentzVector", &dimuon_p4);
  dimuonditrk_tree->Branch("muonP_p4",  "TLorentzVector", &muonP_p4);
  dimuonditrk_tree->Branch("muonN_p4",  "TLorentzVector", &muonN_p4);
  dimuonditrk_tree->Branch("ditrak_p4", "TLorentzVector", &ditrak_p4);
  dimuonditrk_tree->Branch("trakP_p4",  "TLorentzVector", &trakP_p4);
  dimuonditrk_tree->Branch("trakN_p4",  "TLorentzVector", &trakN_p4);

  dimuonditrk_tree->Branch("dimuonditrkTrigger_p4",   "TLorentzVector", &dimuonditrkTrigger_p4);
  dimuonditrk_tree->Branch("ditrakTrigger_p4",     "TLorentzVector", &ditrakTrigger_p4);
  dimuonditrk_tree->Branch("dimuonTrigger_p4",     "TLorentzVector", &dimuonTrigger_p4);

  
  dimuonditrk_tree->Branch("run",      &run,      "run/i");
  dimuonditrk_tree->Branch("event",    &event,    "event/l");
  dimuonditrk_tree->Branch("lumiblock",&lumiblock,"lumiblock/i");



  dimuonditrk_tree->Branch("ndimuonditrk",    &ndimuonditrk,    "ndimuonditrk/i");
  dimuonditrk_tree->Branch("trigger",  &trigger,  "trigger/i");

  dimuonditrk_tree->Branch("isBest",   &isBest,   "isBest/O");

  //2mu+2Trk vertexing
  dimuonditrk_tree->Branch("dimuonditrk_vProb",      &dimuonditrk_vProb,        "dimuonditrk_vProb/D");
  dimuonditrk_tree->Branch("dimuonditrk_vChi2",      &dimuonditrk_vChi2,        "dimuonditrk_vChi2/D");
  dimuonditrk_tree->Branch("dimuonditrk_cosAlpha",   &dimuonditrk_cosAlpha,     "dimuonditrk_cosAlpha/D");
  dimuonditrk_tree->Branch("dimuonditrk_ctauPV",     &dimuonditrk_ctauPV,       "dimuonditrk_ctauPV/D");
  dimuonditrk_tree->Branch("dimuonditrk_ctauErrPV",  &dimuonditrk_ctauErrPV,    "dimuonditrk_ctauErrPV/D");
  dimuonditrk_tree->Branch("dimuonditrk_charge",     &dimuonditrk_charge,       "dimuonditrk_charge/I");
  dimuonditrk_tree->Branch("dimuonditrk_lxy",        &dimuonditrk_lxy,      "dimuonditrk_lxy/F");
  dimuonditrk_tree->Branch("dimuonditrk_lxyErr",     &dimuonditrk_lxyErr,      "dimuonditrk_lxyErr/F");
  dimuonditrk_tree->Branch("dimuonditrk_MassErr",     &dimuonditrk_MassErr,      "dimuonditrk_MassErr/F");

  dimuonditrk_tree->Branch("dimuon_vProb",        &dimuon_vProb,        "dimuon_vProb/D");
  dimuonditrk_tree->Branch("dimuon_vNChi2",       &dimuon_vNChi2,        "dimuon_vNChi2/D");
  dimuonditrk_tree->Branch("dimuon_DCA",          &dimuon_DCA,          "dimuon_DCA/D");
  dimuonditrk_tree->Branch("dimuon_ctauPV",       &dimuon_ctauPV,       "dimuon_ctauPV/D");
  dimuonditrk_tree->Branch("dimuon_ctauErrPV",    &dimuon_ctauErrPV,    "dimuon_ctauErrPV/D");
  dimuonditrk_tree->Branch("dimuon_cosAlpha",     &dimuon_cosAlpha,     "dimuon_cosAlpha/D");

  dimuonditrk_tree->Branch("muonP_tMatch", &muonP_tMatch, "muonP_tMatch/I");
  dimuonditrk_tree->Branch("muonN_tMatch", &muonN_tMatch, "muonN_tMatch/I");
  dimuonditrk_tree->Branch("trakP_tMatch", &trakP_tMatch, "trakP_tMatch/I");
  dimuonditrk_tree->Branch("trakN_tMatch", &trakN_tMatch, "trakN_tMatch/I");

  dimuonditrk_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/i");


}

DiMuonDiTrakVertexRootuplerHLT::~DiMuonDiTrakVertexRootuplerHLT() {}

//
// member functions
//

UInt_t DiMuonDiTrakVertexRootuplerHLT::getTriggerBits(const edm::Event& iEvent ) {

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
void DiMuonDiTrakVertexRootuplerHLT::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  edm::Handle<pat::CompositeCandidateCollection> dimuonditrks;
  iEvent.getByToken(dimuonditrk_Label,dimuonditrks);

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( triggerResults_Label , triggerResults_handle);

  edm::Handle<std::vector<reco::Vertex >> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  run       = iEvent.id().run();
  event     = iEvent.id().event();
  lumiblock = iEvent.id().luminosityBlock();

  numPrimaryVertices = 0;
  if (primaryVertices_handle.isValid()) numPrimaryVertices = (int) primaryVertices_handle->size();
  trigger = getTriggerBits(iEvent);

  const edm::TriggerNames & names = iEvent.triggerNames( *triggerResults_handle );

  ndimuonditrk  = dimuonditrks->size();

  dimuonditrak_p4.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
  dimuon_p4.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
  muonP_p4.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
  muonN_p4.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
  ditrak_p4.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
  trakP_p4.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
  trakN_p4.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
  dimuonditrkTrigger_p4.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
  ditrakTrigger_p4.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
  dimuonTrigger_p4.SetPtEtaPhiM(0.0,0.0,0.0,0.0);

  isBest = true;

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

        const pat::CompositeCandidate *dimuonditrkTrigger_cand = dynamic_cast <const pat::CompositeCandidate *>(dimuonditrkCand->daughter("dimuonTTTrigger"));
        const pat::CompositeCandidate *dimuonTrigger_cand = dynamic_cast <const pat::CompositeCandidate *>(dimuonditrkTrigger_cand->daughter("dimuon"));
        const pat::CompositeCandidate *ditrakTrigger_cand = dynamic_cast <const pat::CompositeCandidate *>(dimuonditrkTrigger_cand->daughter("ditrak"));

        dimuonditrkTrigger_p4.SetPtEtaPhiM(dimuonditrkTrigger_cand->pt(),dimuonditrkTrigger_cand->eta(),dimuonditrkTrigger_cand->phi(),dimuonditrkTrigger_cand->mass());
        dimuonTrigger_p4.SetPtEtaPhiM(dimuonTrigger_cand->pt(),dimuonTrigger_cand->eta(),dimuonTrigger_cand->phi(),dimuonTrigger_cand->mass());
        ditrakTrigger_p4.SetPtEtaPhiM(ditrakTrigger_cand->pt(),ditrakTrigger_cand->eta(),ditrakTrigger_cand->phi(),ditrakTrigger_cand->mass());

        dimuon_vProb        = dimuon_cand->userFloat("vProb");
        dimuon_vNChi2        = dimuon_cand->userFloat("vNChi2");
        dimuon_DCA          = dimuon_cand->userFloat("DCA");
        dimuon_ctauPV       = dimuon_cand->userFloat("ppdlPV");
        dimuon_ctauErrPV    = dimuon_cand->userFloat("ppdlErrPV");
        dimuon_cosAlpha     = dimuon_cand->userFloat("cosAlpha");

        muonP_tMatch  = dimuon_cand->userInt("tMatchP");
        muonN_tMatch  = dimuon_cand->userInt("tMatchN");
        trakP_tMatch  = dimuonditrkCand->userInt("trakMatchP");
        trakN_tMatch  = dimuonditrkCand->userInt("trakMatchN");

        dimuonditrk_MassErr = -1.0;
        if (dimuonditrkCand->hasUserFloat("MassErr"))
          dimuonditrk_MassErr = dimuonditrkCand->userFloat("MassErr");

        dimuonditrk_vProb     = dimuonditrkCand->userFloat("vProb");
        dimuonditrk_vChi2     = dimuonditrkCand->userFloat("vNChi2");
        dimuonditrk_cosAlpha  = dimuonditrkCand->userFloat("cosAlpha");
        dimuonditrk_ctauPV    = dimuonditrkCand->userFloat("ctauPV");
        dimuonditrk_ctauErrPV = dimuonditrkCand->userFloat("ctauErrPV");
        dimuonditrk_ctauErrPV = dimuonditrkCand->userFloat("ctauErrPV");
        dimuonditrk_lxyErr    = dimuonditrkCand->userFloat("lErrxy");
        dimuonditrk_lxy       = dimuonditrkCand->userFloat("lxy");
        dimuonditrk_charge    = dimuonditrkCand->charge();

        dimuonditrk_tree->Fill();
        isBest = false;
        if (OnlyBest_) break;
    }
  }


}

// ------------ method called once each job just before starting event loop  ------------
void DiMuonDiTrakVertexRootuplerHLT::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void DiMuonDiTrakVertexRootuplerHLT::endJob() {}

// ------------ method called when starting to processes a run  ------------
void DiMuonDiTrakVertexRootuplerHLT::beginRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a run  ------------
void DiMuonDiTrakVertexRootuplerHLT::endRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when starting to processes a luminosity block  ------------
void DiMuonDiTrakVertexRootuplerHLT::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void DiMuonDiTrakVertexRootuplerHLT::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DiMuonDiTrakVertexRootuplerHLT::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonDiTrakVertexRootuplerHLT);
