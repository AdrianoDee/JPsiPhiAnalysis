// -*- C++ -*-
//
// Package:    DiTrackRootupler
// Class:      DiTrackRootupler
//
// Description: DiTrackCand(mu+ mu-)  producer
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

class DiTrackRootupler:public edm::EDAnalyzer {
      public:
	explicit DiTrackRootupler(const edm::ParameterSet &);
	~DiTrackRootupler() override;

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
	edm::InputTag DiTrack_Label;
  std::string primaryVertices_Label;
  edm::InputTag hlTriggerResults;
  int  pdgid_;
  std::vector<double> DiTrackMassCuts_;
	bool isMC_;
  bool OnlyBest_;
  bool OnlyGen_;
  std::vector<std::string>  HLTs_, Filters_;

	UInt_t    run;
	ULong64_t event;
  UInt_t    lumiblock;
  UInt_t    nDiTrack;
  UInt_t    ntracks;
  UInt_t    trigger;
  UInt_t    tMatch;
  Int_t     charge;
  Int_t     DiTrackType;
  Int_t isTriggerMatched;

	TLorentzVector DiTrack_p4;
	TLorentzVector trackP_p4;
	TLorentzVector trackN_p4;
  TLorentzVector DiTrack_unref_p4;
  TLorentzVector trackP_ref_p4;
  TLorentzVector trackN_ref_p4;

  Double_t invalidMu, deltaRMuMu, mass, mass_unref;
  Float_t cosAlpha, cosAlpha3D, ctau, ctauErr, lxy, lxyErr, lxyz, lxyzErr;
  Float_t cosAlphaBS, cosAlpha3DBS, ctauBS, ctauErrBS, lxyBS, lxyErrBS, lxyzBS, lxyzErrBS;
  Float_t posTrackDzVtx, posTrackDxyVtx, negTrackDzVtx, negTrackDxyVtx;
  Float_t trPos_Chi2, trPos_NDF, trNeg_Chi2, trNeg_NDF, vProb, Chi2, NDof,SS;

  Int_t nMatchedStationsPos, nOverlapMusPos, nSharingSegWithPos, nMatchedStationsNeg, nOverlapMusNeg;
  Int_t nSharingSegWithNeg, posTrackTrackType, posTrackType, negTrackTrackType, negTrackType;

	UInt_t numPrimaryVertices;

	TTree *DiTrack_tree;

  Int_t mother_pdgId;
  Int_t DiTrack_pdgId;
	TLorentzVector gen_DiTrack_p4;
	TLorentzVector gen_trackP_p4;
	TLorentzVector gen_trackM_p4;

  std::string  genCands_;

};

//
// constructors and destructor
//

DiTrackRootupler::DiTrackRootupler(const edm::ParameterSet & iConfig):
DiTrack_Label(iConfig.getUntrackedParameter<edm::InputTag>("DiTracks")),
primaryVertices_Label(iConfig.getUntrackedParameter<std::string>("primaryVertices",std::string("offlinePrimaryVertices"))),
hlTriggerResults(iConfig.getUntrackedParameter<edm::InputTag>("HLTriggerResults",edm::InputTag("TriggerResults::HLT")) ),
pdgid_(iConfig.getParameter<uint32_t>("DiTrack_pdgid")),
DiTrackMassCuts_(iConfig.getParameter<std::vector<double>>("DiTrack_mass_cuts")),
isMC_(iConfig.getParameter<bool>("isMC")),
OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
HLTs_(iConfig.getUntrackedParameter<std::vector<std::string> >("TriggersForMatching")),
Filters_(iConfig.getUntrackedParameter<std::vector<std::string> >("FiltersForMatching"))
{
  edm::Service < TFileService > fs;
  DiTrack_tree = fs->make < TTree > ("DiTrackTree", "Tree of DiTrackCand");

  DiTrack_tree->Branch("run",      &run,      "run/i");
  DiTrack_tree->Branch("event",    &event,    "event/l");
  DiTrack_tree->Branch("lumiblock",&lumiblock,"lumiblock/i");

  if (!OnlyGen_) {
    DiTrack_tree->Branch("nDiTrack",    &nDiTrack,    "nDiTrack/i");
    DiTrack_tree->Branch("ntracks",   &ntracks,   "ntracks/i");
    DiTrack_tree->Branch("trigger",  &trigger,  "trigger/i");
    DiTrack_tree->Branch("charge",   &charge,   "charge/I");

    DiTrack_tree->Branch("deltaRMuMu",   &deltaRMuMu,   "deltaRMuMu/D");

    DiTrack_tree->Branch("mass",   &mass,   "mass/D");
    DiTrack_tree->Branch("mass_unref",   &mass_unref,   "mass_unref/D");

    DiTrack_tree->Branch("DiTrack_p4", "TLorentzVector", &DiTrack_p4);
    DiTrack_tree->Branch("trackP_p4",  "TLorentzVector", &trackP_p4);
    DiTrack_tree->Branch("trackN_p4",  "TLorentzVector", &trackN_p4);

    DiTrack_tree->Branch("DiTrack_unref_p4", "TLorentzVector", &DiTrack_unref_p4);
    DiTrack_tree->Branch("trackP_ref_p4", "TLorentzVector", &trackP_ref_p4);
    DiTrack_tree->Branch("trackN_ref_p4", "TLorentzVector", &trackN_ref_p4);

    DiTrack_tree->Branch("nMatchedStationsPos",   &nMatchedStationsPos,   "nMatchedStationsPos/I");
    DiTrack_tree->Branch("nOverlapMusPos",  &nOverlapMusPos,  "nOverlapMusPos/I");
    DiTrack_tree->Branch("nSharingSegWithPos",   &nSharingSegWithPos,   "nSharingSegWithPos/I");

    DiTrack_tree->Branch("nMatchedStationsNeg",   &nMatchedStationsNeg,   "nMatchedStationsNeg/I");
    DiTrack_tree->Branch("nOverlapMusNeg",  &nOverlapMusNeg,  "nOverlapMusNeg/I");
    DiTrack_tree->Branch("nSharingSegWithNeg",   &nSharingSegWithNeg,   "nSharingSegWithNeg/I");

    DiTrack_tree->Branch("posTrackTrackType",  &posTrackTrackType,  "posTrackTrackType/I");
    DiTrack_tree->Branch("posTrackType",   &posTrackType,   "posTrackType/I");

    DiTrack_tree->Branch("negTrackTrackType",  &negTrackTrackType,  "negTrackTrackType/I");
    DiTrack_tree->Branch("negTrackType",   &negTrackType,   "negTrackType/I");

    DiTrack_tree->Branch("posTrackDzVtx",  &posTrackDzVtx,  "posTrackDzVtx/D");
    DiTrack_tree->Branch("posTrackDxyVtx",   &posTrackDxyVtx,   "posTrackDxyVtx/D");

    DiTrack_tree->Branch("negTrackDzVtx",  &negTrackDzVtx,  "negTrackDzVtx/D");
    DiTrack_tree->Branch("negTrackDxyVtx",   &negTrackDxyVtx,   "negTrackDxyVtx/D");

    DiTrack_tree->Branch("trPos_Chi2",  &trPos_Chi2,  "trPos_Chi2/D");
    DiTrack_tree->Branch("trPos_NDF",   &trPos_NDF,   "trPos_NDF/D");
    DiTrack_tree->Branch("trNeg_Chi2",  &trNeg_Chi2,  "trNeg_Chi2/D");
    DiTrack_tree->Branch("trNeg_NDF",   &trNeg_NDF,   "trNeg_NDF/D");
    DiTrack_tree->Branch("vProb",  &vProb,  "vProb/D");
    DiTrack_tree->Branch("Chi2",   &Chi2,   "Chi2/D");
    DiTrack_tree->Branch("NDof",   &NDof,   "NDof/D");
    DiTrack_tree->Branch("SS",   &SS,   "SS/D");

    DiTrack_tree->Branch("cosAlpha",  &cosAlpha,  "cosAlpha/D");
    DiTrack_tree->Branch("cosAlpha3D",   &cosAlpha3D,   "cosAlpha3D/D");
    DiTrack_tree->Branch("ctau",  &ctau,  "ctau/D");
    DiTrack_tree->Branch("ctauErr",   &ctauErr,   "ctauErr/D");
    DiTrack_tree->Branch("lxy",  &lxy,  "lxy/D");
    DiTrack_tree->Branch("lxyErr",   &lxyErr,   "lxyErr/D");
    DiTrack_tree->Branch("lxyz",   &lxyz,   "lxyz/D");
    DiTrack_tree->Branch("lxyzErr",   &lxyzErr,   "lxyzErr/D");

    DiTrack_tree->Branch("cosAlphaBS",  &cosAlphaBS,  "cosAlphaBS/D");
    DiTrack_tree->Branch("cosAlpha3DBS",   &cosAlpha3DBS,   "cosAlpha3DBS/D");
    DiTrack_tree->Branch("ctauBS",  &ctauBS,  "ctauBS/D");
    DiTrack_tree->Branch("ctauErrBS",   &ctauErrBS,   "ctauErrBS/D");
    DiTrack_tree->Branch("lxyBS",  &lxyBS,  "lxyBS/D");
    DiTrack_tree->Branch("lxyErrBS",   &lxyErrBS,   "lxyErrBS/D");
    DiTrack_tree->Branch("lxyzBS",   &lxyzBS,   "lxyzBS/D");
    DiTrack_tree->Branch("lxyzErrBS",   &lxyzErrBS,   "lxyzErrBS/D");

    DiTrack_tree->Branch("invalidMu",   &invalidMu,   "invalidMu/D");
    DiTrack_tree->Branch("DiTrackType",   &DiTrackType,   "lxyzErrBS/I");
    DiTrack_tree->Branch("isTriggerMatched",   &isTriggerMatched,   "lxyzErrBS/I");

    DiTrack_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/i");
  }

  if (isMC_ || OnlyGen_) {
     std::cout << "DiTrackRootupler::DiTrackRootupler: DiTrackCand id " << pdgid_ << std::endl;
     DiTrack_tree->Branch("mother_pdgId",  &mother_pdgId,     "mother_pdgId/I");
     DiTrack_tree->Branch("DiTrack_pdgId",  &DiTrack_pdgId,     "DiTrack_pdgId/I");
     DiTrack_tree->Branch("gen_DiTrack_p4", "TLorentzVector",  &gen_DiTrack_p4);
     DiTrack_tree->Branch("gen_trackP_p4",  "TLorentzVector",  &gen_trackP_p4);
     DiTrack_tree->Branch("gen_trackN_p4",  "TLorentzVector",  &gen_trackM_p4);
  }
  genCands_ = "prunedGenParticles";
}

DiTrackRootupler::~DiTrackRootupler() {}

//
// member functions
//

const reco::Candidate* DiTrackRootupler::GetAncestor(const reco::Candidate* p) {
   if (p->numberOfMothers()) {
      if  ((p->mother(0))->pdgId() == p->pdgId()) return GetAncestor(p->mother(0));
      else return p->mother(0);
   }
   return p;
}

//Check recursively if any ancestor of particle is the given one
bool DiTrackRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
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

UInt_t DiTrackRootupler::getTriggerBits(const edm::Event& iEvent ) {

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
void DiTrackRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  using namespace reco;
  using namespace pat;

  edm::Handle<pat::CompositeCandidateCollection> DiTracks;
  iEvent.getByLabel(DiTrack_Label,DiTracks);

  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByLabel(primaryVertices_Label, primaryVertices_handle);

  run       = iEvent.id().run();
  event     = iEvent.id().event();
  lumiblock = iEvent.id().luminosityBlock();

  numPrimaryVertices = 0;
  if (primaryVertices_handle.isValid()) numPrimaryVertices = (int) primaryVertices_handle->size();
  trigger = getTriggerBits(iEvent);

  DiTrack_pdgId = 0;
  mother_pdgId = 0;
  nDiTrack  = 0;
  ntracks = 0;

  DiTrack_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  trackP_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  trackN_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  DiTrack_unref_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  trackP_ref_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  trackN_ref_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_DiTrack_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_trackP_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_trackM_p4.SetPtEtaPhiM(0.,0.,0.,0.);

  // Pruned particles are the one containing "important" stuff
  edm::Handle<reco::GenParticleCollection> pruned;
  iEvent.getByLabel(genCands_, pruned);

  // // Packed particles are all the status 1. The navigation to pruned is possible (the other direction should be made by hand)
  // edm::Handle<pat::PackedGenParticleCollection> packed;
  // iEvent.getByLabel(packCands_,  packed);
  //
  // if ( (isMC_ || OnlyGen_) && packed.isValid() && pruned.isValid() ) {
  //   for (size_t i=0; i<pruned->size(); i++) {
  //     const reco::Candidate *aDiTrack = &(*pruned)[i];
  //     if ( (abs(aDiTrack->pdgId()) == pdgid_) && (aDiTrack->status() == 2) ) {
  //       int foundit = 1;
  //       DiTrack_pdgId = aDiTrack->pdgId();
  //       for ( size_t j=0; j<packed->size(); j++ ) { //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
  //         const reco::Candidate * motherInPrunedCollection = (*packed)[j].mother(0);
  //         const reco::Candidate * d = &(*packed)[j];
  //         if ( motherInPrunedCollection != nullptr && (d->pdgId() ==  13 ) && isAncestor(aDiTrack , motherInPrunedCollection) ) {
  //           gen_trackM_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
  //           foundit++;
  //         }
  //         if ( motherInPrunedCollection != nullptr && (d->pdgId() == -13 ) && isAncestor(aDiTrack , motherInPrunedCollection) ) {
  //           gen_trackP_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
  //           foundit++;
  //         }
  //         if ( foundit == 3 ) break;
  //       }
  //       if ( foundit == 3 ) {
  //         gen_DiTrack_p4 = gen_trackM_p4 + gen_trackP_p4;   // this should take into account FSR
  //         mother_pdgId  = GetAncestor(aDiTrack)->pdgId();
  //         break;
  //       } else DiTrack_pdgId = 0;
  //     }  // if ( p_id
  //   } // for (size
  //   if ( ! DiTrack_pdgId ) std::cout << "DiTrackRootupler: does not found the given decay " << run << "," << event << std::endl; // sanity check
  // }  // end if isMC

  float DiTrackMassMax_ = DiTrackMassCuts_[1];
  float DiTrackMassMin_ = DiTrackMassCuts_[0];

  bool already_stored = false;
  if ( ! OnlyGen_ ) { // we will look for DiTracks, then for tracks
    if ( DiTracks.isValid() && !DiTracks->empty()) {
      std::cout << "DiTrack Valid!" << std::endl;

      for ( pat::CompositeCandidateCollection::const_iterator DiTrackCand = DiTracks->begin(); DiTrackCand != DiTracks->end(); ++DiTrackCand ) {
        std::cout << "DiTrack Size:" << DiTracks->size() <<std::endl;
        if (DiTrackCand->mass() > DiTrackMassMin_ && DiTrackCand->mass() < DiTrackMassMax_ && DiTrackCand->charge() == 0) {
          DiTrack_p4.SetPtEtaPhiM(DiTrackCand->pt(),DiTrackCand->eta(),DiTrackCand->phi(),DiTrackCand->mass());

          reco::Candidate::LorentzVector vP = DiTrackCand->daughter("trackPos")->p4();
          reco::Candidate::LorentzVector vM = DiTrackCand->daughter("trackNeg")->p4();

          trackP_p4.SetPtEtaPhiM(vP.pt(),vP.eta(),vP.phi(),vP.mass());
          trackN_p4.SetPtEtaPhiM(vM.pt(),vM.eta(),vM.phi(),vM.mass());

          vProb = DiTrackCand->userFloat("vProb");
          // DCA = -1.;
          // if (DiTrackCand->hasUserFloat("DCA"))  DCA = DiTrackCand->userFloat("DCA");

          deltaRMuMu = DiTrackCand->userFloat("deltaR");
          mass = DiTrack_p4.M();
          // mass_unref = DiTrackCand->userFloat("mumuP4");

          reco::Candidate::LorentzVector unrefMuMu = DiTrackCand->daughter("trktrkCandidate")->p4();

          DiTrack_unref_p4.SetPtEtaPhiM(unrefMuMu.pt(),unrefMuMu.eta(),unrefMuMu.phi(),unrefMuMu.mass());

          vP = DiTrackCand->daughter("ref_kaonPos")->p4();
          vM = DiTrackCand->daughter("ref_kaonNeg")->p4();

          trackP_ref_p4.SetPtEtaPhiM(vP.pt(),vP.eta(),vP.phi(),vP.mass());
          trackN_ref_p4.SetPtEtaPhiM(vM.pt(),vM.eta(),vM.phi(),vM.mass());

          posTrackDzVtx =  DiTrackCand->userFloat("trackPosDzVtx");
          posTrackDxyVtx = DiTrackCand->userFloat("trackPosDxyVtx");

          negTrackDzVtx = DiTrackCand->userFloat("trackNegDzVtx");
          negTrackDxyVtx = DiTrackCand->userFloat("trackNegDxyVtx");

          trPos_Chi2 = DiTrackCand->userFloat("trPos_Chi2");
          trPos_NDF = DiTrackCand->userFloat("trPos_NDF");
          trNeg_Chi2 = DiTrackCand->userFloat("trNeg_Chi2");
          trNeg_NDF = DiTrackCand->userFloat("trNeg_NDF");

          vProb = DiTrackCand->userFloat("vProb");
          Chi2 = DiTrackCand->userFloat("Chi2");
          NDof = DiTrackCand->userFloat("NDof");
          SS   = DiTrackCand->userFloat("SS");

          cosAlpha = DiTrackCand->userFloat("cosAlpha");
          cosAlpha3D = DiTrackCand->userFloat("cosAlpha3D");
          ctau = DiTrackCand->userFloat("ctau");
          ctauErr = DiTrackCand->userFloat("ctauErr");
          lxy = DiTrackCand->userFloat("lxy");
          lxyErr = DiTrackCand->userFloat("lxyErr");
          lxyz = DiTrackCand->userFloat("lxyz");
          lxyzErr = DiTrackCand->userFloat("lxyzErr");

          cosAlphaBS =  DiTrackCand->userFloat("cosAlphaBS");
          cosAlpha3DBS =  DiTrackCand->userFloat("cosAlpha3DBS");
          ctauBS =  DiTrackCand->userFloat("ctauBS");
          ctauErrBS = DiTrackCand->userFloat("ctauErrBS");
          lxyBS = DiTrackCand->userFloat("lxyBS");
          lxyErrBS =  DiTrackCand->userFloat("lxyErrBS");
          lxyzBS =  DiTrackCand->userFloat("lxyzBS");
          lxyzErrBS = DiTrackCand->userFloat("lxyzErrBS");

          invalidMu = DiTrackCand->userFloat("isEventWithInvalidMu");
          DiTrackType = DiTrackCand->userInt("DiTrackType");
          isTriggerMatched = DiTrackCand->userInt("isTriggerMatched");

          // pat_ref_JPsi.addUserData("tracklessPV");

          charge = DiTrackCand->charge();

          nDiTrack++;

          DiTrack_tree->Fill();   // be aware, we are storing all combinations
          already_stored = true;


          if (OnlyBest_ && already_stored) break;
        }
      }
    } //..else {
      //std::cout << "DiTrackRootupler: (" << run << "," << event << ") -> ";

  }  // !OnlyGen_

  if ( !already_stored ) {  // we have to make sure, we are not double storing an combination
    if ( !isMC_ ) {
      if ( nDiTrack > 0 ) DiTrack_tree->Fill();   // if not MC filter out
    } else DiTrack_tree->Fill();
  }
}

// ------------ method called once each job just before starting event loop  ------------
void DiTrackRootupler::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void DiTrackRootupler::endJob() {}

// ------------ method called when starting to processes a run  ------------
void DiTrackRootupler::beginRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a run  ------------
void DiTrackRootupler::endRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when starting to processes a luminosity block  ------------
void DiTrackRootupler::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void DiTrackRootupler::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DiTrackRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiTrackRootupler);
