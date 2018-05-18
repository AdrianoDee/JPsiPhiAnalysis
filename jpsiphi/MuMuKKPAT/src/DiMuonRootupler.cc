// -*- C++ -*-
//
// Package:    DiMuonRootupler
// Class:      DiMuonRootupler
//
// Description: dimuonCand(mu+ mu-)  producer
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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

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
	std::string dimuon_Label;
  std::string primaryVertices_Label;
  edm::InputTag hlTriggerResults;
  int  pdgid_;
  std::vector<double> DimuonMassCuts_;
	bool isMC_;
  bool OnlyBest_;
  bool OnlyGen_;
  std::vector<std::string>  HLTs_, Filters_;

	UInt_t    run;
	ULong64_t event;
  UInt_t    lumiblock;
  UInt_t    ndimuon;
  UInt_t    nmuons;
  UInt_t    trigger;
  UInt_t    tMatch;
  Int_t     charge;
  Int_t     dimuonType;
  Int_t isTriggerMatched;

	TLorentzVector dimuon_p4;
	TLorentzVector muonP_p4;
	TLorentzVector muonN_p4;
  TLorentzVector dimuon_unref_p4;
  TLorentzVector muonP_ref_p4;
  TLorentzVector muonN_ref_p4;

  Float_t invalidMu, deltaRMuMu, mass, mass_unref;
  Float_t cosAlpha, cosAlpha3D, ctau, ctauErr, lxy, lxyErr, lxyz, lxyzErr;
  Float_t cosAlphaBS, cosAlpha3DBS, ctauBS, ctauErrBS, lxyBS, lxyErrBS, lxyzBS, lxyzErrBS;
  Float_t posMuonDzVtx, posMuonDxyVtx, negMuonDzVtx, negMuonDxyVtx;
  Float_t muPos_Chi2, muPos_NDF, muNeg_Chi2, muNeg_NDF, vProb, Chi2, NDof;

  Int_t nMatchedStationsPos, nOverlapMusPos, nSharingSegWithPos, nMatchedStationsNeg, nOverlapMusNeg;
  Int_t nSharingSegWithNeg, posMuonTrackType, posMuonType, negMuonTrackType, negMuonType;

	UInt_t numPrimaryVertices;

	TTree *dimuon_tree;

  Int_t mother_pdgId;
  Int_t dimuon_pdgId;
	TLorentzVector gen_dimuon_p4;
	TLorentzVector gen_muonP_p4;
	TLorentzVector gen_muonM_p4;

  std::string  genCands_;

};

//
// constructors and destructor
//

DiMuonRootupler::DiMuonRootupler(const edm::ParameterSet & iConfig):
dimuon_Label(iConfig.getUntrackedParameter<std::string>("dimuons")),
primaryVertices_Label(iConfig.getUntrackedParameter<std::string>("primaryVertices",std::string("offlinePrimaryVertices"))),
hlTriggerResults(iConfig.getUntrackedParameter<edm::InputTag>("HLTriggerResults",edm::InputTag("TriggerResults::HLT")) ),
pdgid_(iConfig.getParameter<uint32_t>("dimuon_pdgid")),
DimuonMassCuts_(iConfig.getParameter<std::vector<double>>("dimuon_mass_cuts")),
isMC_(iConfig.getParameter<bool>("isMC")),
OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
HLTs_(iConfig.getUntrackedParameter<std::vector<std::string> >("TriggersForMatching")),
Filters_(iConfig.getUntrackedParameter<std::vector<std::string> >("FiltersForMatching"))
{
  edm::Service < TFileService > fs;
  dimuon_tree = fs->make < TTree > ("dimuonTree", "Tree of dimuonCand");

  dimuon_tree->Branch("run",      &run,      "run/i");
  dimuon_tree->Branch("event",    &event,    "event/l");
  dimuon_tree->Branch("lumiblock",&lumiblock,"lumiblock/i");

  if (!OnlyGen_) {
    dimuon_tree->Branch("ndimuon",    &ndimuon,    "ndimuon/i");
    dimuon_tree->Branch("nmuons",   &nmuons,   "nmuons/i");
    dimuon_tree->Branch("trigger",  &trigger,  "trigger/i");
    dimuon_tree->Branch("charge",   &charge,   "charge/I");

    dimuon_tree->Branch("deltaRMuMu",   &deltaRMuMu,   "deltaRMuMu/D");

    dimuon_tree->Branch("mass",   &mass,   "mass/D");
    dimuon_tree->Branch("mass_unref",   &mass_unref,   "mass_unref/D");

    dimuon_tree->Branch("dimuon_p4", "TLorentzVector", &dimuon_p4);
    dimuon_tree->Branch("muonP_p4",  "TLorentzVector", &muonP_p4);
    dimuon_tree->Branch("muonN_p4",  "TLorentzVector", &muonN_p4);

    dimuon_tree->Branch("dimuon_unref_p4", "TLorentzVector", &dimuon_unref_p4);
    dimuon_tree->Branch("muonP_ref_p4", "TLorentzVector", &muonP_ref_p4);
    dimuon_tree->Branch("muonN_ref_p4", "TLorentzVector", &muonN_ref_p4);

    dimuon_tree->Branch("nMatchedStationsPos",   &nMatchedStationsPos,   "nMatchedStationsPos/I");
    dimuon_tree->Branch("nOverlapMusPos",  &nOverlapMusPos,  "nOverlapMusPos/I");
    dimuon_tree->Branch("nSharingSegWithPos",   &nSharingSegWithPos,   "nSharingSegWithPos/I");

    dimuon_tree->Branch("nMatchedStationsNeg",   &nMatchedStationsNeg,   "nMatchedStationsNeg/I");
    dimuon_tree->Branch("nOverlapMusNeg",  &nOverlapMusNeg,  "nOverlapMusNeg/I");
    dimuon_tree->Branch("nSharingSegWithNeg",   &nSharingSegWithNeg,   "nSharingSegWithNeg/I");

    dimuon_tree->Branch("posMuonTrackType",  &posMuonTrackType,  "posMuonTrackType/I");
    dimuon_tree->Branch("posMuonType",   &posMuonType,   "posMuonType/I");

    dimuon_tree->Branch("negMuonTrackType",  &negMuonTrackType,  "negMuonTrackType/I");
    dimuon_tree->Branch("negMuonType",   &negMuonType,   "negMuonType/I");

    dimuon_tree->Branch("posMuonDzVtx",  &posMuonDzVtx,  "posMuonDzVtx/D");
    dimuon_tree->Branch("posMuonDxyVtx",   &posMuonDxyVtx,   "posMuonDxyVtx/D");

    dimuon_tree->Branch("negMuonDzVtx",  &negMuonDzVtx,  "negMuonDzVtx/D");
    dimuon_tree->Branch("negMuonDxyVtx",   &negMuonDxyVtx,   "negMuonDxyVtx/D");

    dimuon_tree->Branch("muPos_Chi2",  &muPos_Chi2,  "muPos_Chi2/D");
    dimuon_tree->Branch("muPos_NDF",   &muPos_NDF,   "muPos_NDF/D");
    dimuon_tree->Branch("muNeg_Chi2",  &muNeg_Chi2,  "muNeg_Chi2/D");
    dimuon_tree->Branch("muNeg_NDF",   &muNeg_NDF,   "muNeg_NDF/D");
    dimuon_tree->Branch("vProb",  &vProb,  "vProb/D");
    dimuon_tree->Branch("Chi2",   &Chi2,   "Chi2/D");
    dimuon_tree->Branch("NDof",   &NDof,   "NDof/D");

    dimuon_tree->Branch("cosAlpha",  &cosAlpha,  "cosAlpha/D");
    dimuon_tree->Branch("cosAlpha3D",   &cosAlpha3D,   "cosAlpha3D/D");
    dimuon_tree->Branch("ctau",  &ctau,  "ctau/D");
    dimuon_tree->Branch("ctauErr",   &ctauErr,   "ctauErr/D");
    dimuon_tree->Branch("lxy",  &lxy,  "lxy/D");
    dimuon_tree->Branch("lxyErr",   &lxyErr,   "lxyErr/D");
    dimuon_tree->Branch("lxyz",   &lxyz,   "lxyz/D");
    dimuon_tree->Branch("lxyzErr",   &lxyzErr,   "lxyzErr/D");

    dimuon_tree->Branch("cosAlphaBS",  &cosAlphaBS,  "cosAlphaBS/D");
    dimuon_tree->Branch("cosAlpha3DBS",   &cosAlpha3DBS,   "cosAlpha3DBS/D");
    dimuon_tree->Branch("ctauBS",  &ctauBS,  "ctauBS/D");
    dimuon_tree->Branch("ctauErrBS",   &ctauErrBS,   "ctauErrBS/D");
    dimuon_tree->Branch("lxyBS",  &lxyBS,  "lxyBS/D");
    dimuon_tree->Branch("lxyErrBS",   &lxyErrBS,   "lxyErrBS/D");
    dimuon_tree->Branch("lxyzBS",   &lxyzBS,   "lxyzBS/D");
    dimuon_tree->Branch("lxyzErrBS",   &lxyzErrBS,   "lxyzErrBS/D");

    dimuon_tree->Branch("invalidMu",   &invalidMu,   "invalidMu/D");
    dimuon_tree->Branch("dimuonType",   &dimuonType,   "lxyzErrBS/I");
    dimuon_tree->Branch("isTriggerMatched",   &isTriggerMatched,   "lxyzErrBS/I");

    dimuon_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/i");
  }

  if (isMC_ || OnlyGen_) {
     std::cout << "DiMuonRootupler::DiMuonRootupler: dimuonCand id " << pdgid_ << std::endl;
     dimuon_tree->Branch("mother_pdgId",  &mother_pdgId,     "mother_pdgId/I");
     dimuon_tree->Branch("dimuon_pdgId",  &dimuon_pdgId,     "dimuon_pdgId/I");
     dimuon_tree->Branch("gen_dimuon_p4", "TLorentzVector",  &gen_dimuon_p4);
     dimuon_tree->Branch("gen_muonP_p4",  "TLorentzVector",  &gen_muonP_p4);
     dimuon_tree->Branch("gen_muonN_p4",  "TLorentzVector",  &gen_muonM_p4);
  }
  genCands_ = "prunedGenParticles";
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

//Check recursively if any ancestor of particle is the given one
bool DiMuonRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   for (size_t i=0; i< particle->numberOfMothers(); i++) {
      if (isAncestor(ancestor, particle->mother(i))) return true;
   }
   return false;
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
  iEvent.getByLabel( hlTriggerResults , triggerResults_handle);

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

  using namespace reco;
  using namespace pat;

  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  iEvent.getByLabel(dimuon_Label,dimuons);

  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByLabel(primaryVertices_Label, primaryVertices_handle);

  run       = iEvent.id().run();
  event     = iEvent.id().event();
  lumiblock = iEvent.id().luminosityBlock();

  numPrimaryVertices = 0;
  if (primaryVertices_handle.isValid()) numPrimaryVertices = (int) primaryVertices_handle->size();
  trigger = getTriggerBits(iEvent);

  dimuon_pdgId = 0;
  mother_pdgId = 0;
  ndimuon  = 0;
  nmuons = 0;

  dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  muonP_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  muonN_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  dimuon_unref_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  muonP_ref_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  muonN_ref_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muonP_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muonM_p4.SetPtEtaPhiM(0.,0.,0.,0.);

  // Pruned particles are the one containing "important" stuff
  edm::Handle<reco::GenParticleCollection> pruned;
  iEvent.getByLabel(genCands_, pruned);

  // // Packed particles are all the status 1. The navigation to pruned is possible (the other direction should be made by hand)
  // edm::Handle<pat::PackedGenParticleCollection> packed;
  // iEvent.getByLabel(packCands_,  packed);
  //
  // if ( (isMC_ || OnlyGen_) && packed.isValid() && pruned.isValid() ) {
  //   for (size_t i=0; i<pruned->size(); i++) {
  //     const reco::Candidate *adimuon = &(*pruned)[i];
  //     if ( (abs(adimuon->pdgId()) == pdgid_) && (adimuon->status() == 2) ) {
  //       int foundit = 1;
  //       dimuon_pdgId = adimuon->pdgId();
  //       for ( size_t j=0; j<packed->size(); j++ ) { //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
  //         const reco::Candidate * motherInPrunedCollection = (*packed)[j].mother(0);
  //         const reco::Candidate * d = &(*packed)[j];
  //         if ( motherInPrunedCollection != nullptr && (d->pdgId() ==  13 ) && isAncestor(adimuon , motherInPrunedCollection) ) {
  //           gen_muonM_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
  //           foundit++;
  //         }
  //         if ( motherInPrunedCollection != nullptr && (d->pdgId() == -13 ) && isAncestor(adimuon , motherInPrunedCollection) ) {
  //           gen_muonP_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
  //           foundit++;
  //         }
  //         if ( foundit == 3 ) break;
  //       }
  //       if ( foundit == 3 ) {
  //         gen_dimuon_p4 = gen_muonM_p4 + gen_muonP_p4;   // this should take into account FSR
  //         mother_pdgId  = GetAncestor(adimuon)->pdgId();
  //         break;
  //       } else dimuon_pdgId = 0;
  //     }  // if ( p_id
  //   } // for (size
  //   if ( ! dimuon_pdgId ) std::cout << "DiMuonRootupler: does not found the given decay " << run << "," << event << std::endl; // sanity check
  // }  // end if isMC

  float DimuonMassMax_ = DimuonMassCuts_[1];
  float DimuonMassMin_ = DimuonMassCuts_[0];

  bool already_stored = false;
  if ( ! OnlyGen_ ) { // we will look for dimuons, then for muons
    if ( dimuons.isValid() && !dimuons->empty()) {
      for ( pat::CompositeCandidateCollection::const_iterator dimuonCand = dimuons->begin(); dimuonCand != dimuons->end(); ++dimuonCand ) {
        if (dimuonCand->mass() > DimuonMassMin_ && dimuonCand->mass() < DimuonMassMax_ && dimuonCand->charge() == 0) {
          dimuon_p4.SetPtEtaPhiM(dimuonCand->pt(),dimuonCand->eta(),dimuonCand->phi(),dimuonCand->mass());

          reco::Candidate::LorentzVector vP = dimuonCand->daughter("posMuon")->p4();
          reco::Candidate::LorentzVector vM = dimuonCand->daughter("negMuon")->p4();

          muonP_p4.SetPtEtaPhiM(vP.pt(),vP.eta(),vP.phi(),vP.mass());
          muonN_p4.SetPtEtaPhiM(vM.pt(),vM.eta(),vM.phi(),vM.mass());

          vProb = dimuonCand->userFloat("vProb");
          // DCA = -1.;
          // if (dimuonCand->hasUserFloat("DCA"))  DCA = dimuonCand->userFloat("DCA");

          deltaRMuMu = dimuonCand->userFloat("deltaR");
          mass = dimuon_p4.M();
          mass_unref = dimuonCand->userFloat("mumuP4");

          reco::Candidate::LorentzVector unrefMuMu = dimuonCand->daughter("mumuCandidate")->p4();

          dimuon_unref_p4.SetPtEtaPhiM(unrefMuMu.pt(),unrefMuMu.eta(),unrefMuMu.phi(),unrefMuMu.mass());

          vM = dimuonCand->daughter("ref_muonNeg")->p4();
          vP = dimuonCand->daughter("ref_muonPos")->p4();

          muonP_ref_p4.SetPtEtaPhiM(vP.pt(),vP.eta(),vP.phi(),vP.mass());
          muonN_ref_p4.SetPtEtaPhiM(vM.pt(),vM.eta(),vM.phi(),vM.mass());

          nMatchedStationsPos = dimuonCand->userInt("nMatchedStationsPos");
          nOverlapMusPos =  dimuonCand->userInt("nOverlapMusPos");
          nSharingSegWithPos =  dimuonCand->userInt("nSharingSegWithPos");

          nMatchedStationsNeg = dimuonCand->userInt("nMatchedStationsNeg");
          nOverlapMusNeg =  dimuonCand->userInt("nOverlapMusNeg");
          nSharingSegWithNeg =  dimuonCand->userInt("nSharingSegWithNeg");

          posMuonTrackType =  dimuonCand->userInt("posMuonTrackType");
          posMuonType = dimuonCand->userInt("posMuonType");

          posMuonDzVtx =  dimuonCand->userFloat("posMuonDzVtx");
          posMuonDxyVtx = dimuonCand->userFloat("posMuonDxyVtx");

          negMuonTrackType =  dimuonCand->userInt("negMuonTrackType");
          negMuonType = dimuonCand->userInt("negMuonType");

          negMuonDzVtx = dimuonCand->userFloat("negMuonDzVtx");
          negMuonDxyVtx = dimuonCand->userFloat("negMuonDxyVtx");
          muPos_Chi2 = dimuonCand->userFloat("muPos_Chi2");
          muPos_NDF = dimuonCand->userFloat("muPos_NDF");
          muNeg_Chi2 = dimuonCand->userFloat("muNeg_Chi2");
          muNeg_NDF = dimuonCand->userFloat("muNeg_NDF");

          vProb = dimuonCand->userFloat("vProb");
          Chi2 = dimuonCand->userFloat("Chi2");
          NDof = dimuonCand->userFloat("NDof");

          cosAlpha = dimuonCand->userFloat("cosAlpha");
          cosAlpha3D = dimuonCand->userFloat("cosAlpha3D");
          ctau = dimuonCand->userFloat("ctau");
          ctauErr = dimuonCand->userFloat("ctauErr");
          lxy = dimuonCand->userFloat("lxy");
          lxyErr = dimuonCand->userFloat("lxyErr");
          lxyz = dimuonCand->userFloat("lxyz");
          lxyzErr = dimuonCand->userFloat("lxyzErr");

          cosAlphaBS =  dimuonCand->userFloat("cosAlphaBS");
          cosAlpha3DBS =  dimuonCand->userFloat("cosAlpha3DBS");
          ctauBS =  dimuonCand->userFloat("ctauBS");
          ctauErrBS = dimuonCand->userFloat("ctauErrBS");
          lxyBS = dimuonCand->userFloat("lxyBS");
          lxyErrBS =  dimuonCand->userFloat("lxyErrBS");
          lxyzBS =  dimuonCand->userFloat("lxyzBS");
          lxyzErrBS = dimuonCand->userFloat("lxyzErrBS");

          invalidMu = dimuonCand->userFloat("isEventWithInvalidMu");
          dimuonType = dimuonCand->userInt("dimuonType");
          isTriggerMatched = dimuonCand->userInt("isTriggerMatched");

          // pat_ref_JPsi.addUserData("muonlessPV");

          charge = dimuonCand->charge();

          ndimuon++;
          if (OnlyBest_) break;
          else {
            dimuon_tree->Fill();   // be aware, we are storing all combinations
            already_stored = true;
          }
        }
      }
    } //..else {
      //std::cout << "DiMuonRootupler: (" << run << "," << event << ") -> ";

  }  // !OnlyGen_

  if ( !already_stored ) {  // we have to make sure, we are not double storing an combination
    if ( !isMC_ ) {
      if ( ndimuon > 0 ) dimuon_tree->Fill();   // if not MC filter out
    } else dimuon_tree->Fill();
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
