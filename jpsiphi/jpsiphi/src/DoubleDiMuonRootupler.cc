/*
   Package:    DoubleDiMuonRootupler
   Class:      DoubleDiMuonRootupler

   Description: make rootuple of DiMuon DiTrack combination

   Original Author:  Adriano Di Florio
   Created:  based on Alberto Sanchez Hernandez RootuplePsiDiTrak

*/

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/UserData.h"

#include "TLorentzVector.h"
#include "TTree.h"
#include <vector>
#include <sstream>

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

//
// class declaration
//


class DoubleDiMuonRootupler : public edm::EDAnalyzer {
   public:
      explicit DoubleDiMuonRootupler(const edm::ParameterSet&);
      ~DoubleDiMuonRootupler() override;

      bool isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle);
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      void beginJob() override ;
      void analyze(const edm::Event&, const edm::EventSetup&) override;
      void endJob() override ;

      void beginRun(edm::Run const&, edm::EventSetup const&) override;
      void endRun(edm::Run const&, edm::EventSetup const&) override;
      void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      UInt_t isTriggerMatched(pat::CompositeCandidate *diMuon_cand);

  // ----------member data ---------------------------
  std::string file_name;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> doubledimuon_cand_Label;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> doubledimuon_rf_cand_Label;
  edm::EDGetTokenT<reco::BeamSpot> thebeamspot_;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  int  doubledimuon_pdgid_, higdim_pdgid_, lowdim_pdgid_;
  bool isMC_,OnlyBest_;
  std::vector<std::string>                            HLTs_;
  std::vector<std::string>                            HLTFilters_;

  UInt_t run,event,numPrimaryVertices, trigger;

  TLorentzVector doubledimuon_p4;
  TLorentzVector higdim_p4;
  TLorentzVector lowdim_p4;
  TLorentzVector muonLowP_p4;
  TLorentzVector muonHighN_p4;
  TLorentzVector muonHighP_p4;
  TLorentzVector muonLowN_p4;

  TLorentzVector doubledimuon_rf_p4;
  TLorentzVector higdim_rf_p4;
  TLorentzVector lowdim_rf_p4;
  TLorentzVector muonLowP_rf_p4;
  TLorentzVector muonHighN_rf_p4;
  TLorentzVector muonHighP_rf_p4;
  TLorentzVector muonLowN_rf_p4;
  TLorentzVector doubledimuon_not_rf_p4;
  TLorentzVector higdim_not_rf_p4;
  TLorentzVector lowdim_not_rf_p4;

  Int_t    doubledimuon_charge, higdim_triggerMatch, lowdim_triggerMatch, doubledimuon_rf_bindx;
  Double_t doubledimuon_vProb,  doubledimuon_vChi2, doubledimuon_cosAlpha, doubledimuon_ctauPV, doubledimuon_ctauErrPV;
  Double_t doubledimuon_rf_vProb,  doubledimuon_rf_vChi2, doubledimuon_rf_cosAlpha, doubledimuon_rf_ctauPV, doubledimuon_rf_ctauErrPV;

  Int_t    higdim_triggerMatch_rf, lowdim_triggerMatch_rf;
  Double_t higdim_vProb_rf, higdim_vChi2_rf, higdim_DCA_rf, higdim_ctauPV_rf, higdim_ctauErrPV_rf, higdim_cosAlpha_rf;
  Double_t lowdim_vProb_rf, lowdim_vChi2_rf, lowdim_DCA_rf, lowdim_ctauPV_rf, lowdim_ctauErrPV_rf, lowdim_cosAlpha_rf;

  //Double_t track_d0, track_d0Err, track_dz, track_dxy;
  Double_t higdim_vProb, higdim_vChi2, higdim_DCA, higdim_ctauPV, higdim_ctauErrPV, higdim_cosAlpha;
  Double_t lowdim_vProb, lowdim_vChi2, lowdim_DCA, lowdim_ctauPV, lowdim_ctauErrPV, lowdim_cosAlpha;
  Double_t  highDiMM_fit, highDiMPx_fit, highDiMPy_fit, highDiMPz_fit;

  Bool_t muonHighP_isLoose, muonHighP_isSoft, muonHighP_isMedium, muonHighP_isHighPt;
  Bool_t muonHighN_isLoose, muonHighN_isSoft, muonHighN_isMedium, muonHighN_isHighPt;
  Bool_t muonLowP_isLoose, muonLowP_isSoft, muonLowP_isMedium, muonLowP_isHighPt;
  Bool_t muonLowN_isLoose, muonLowN_isSoft, muonLowN_isMedium, muonLowN_isHighPt;

  Bool_t muonHighP_isTracker, muonHighP_isGlobal, muonHighN_isTracker, muonHighN_isGlobal;
  Bool_t muonLowP_isTracker, muonLowP_isGlobal, muonLowN_isTracker, muonLowN_isGlobal;

  Bool_t isBestCandidate;

  UInt_t muonHighP_type, muonHighN_type, muonLowP_type, muonLowN_type;
  Int_t          gen_doubledimuon_pdgId;
  TLorentzVector gen_doubledimuon_p4;
  TLorentzVector gen_higdim_p4;
  TLorentzVector gen_lowdim_p4;
  TLorentzVector gen_muonLowP_p4;
  TLorentzVector gen_muonHighN_p4;
  TLorentzVector gen_muonHighP_p4;
  TLorentzVector gen_muonLowN_p4;

  TTree* doubledimuon_tree, *doubledimuon_tree_rf;
  edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;
};

//
// constants, enums and typedefs
//

UInt_t DoubleDiMuonRootupler::isTriggerMatched(pat::CompositeCandidate *diMuon_cand) {
  const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muon1"));
  const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muon2"));
  UInt_t matched = 0;  // if no list is given, is not matched

  // if matched a given trigger, set the bit, in the same order as listed
  for (unsigned int iTr = 0; iTr<HLTFilters_.size(); iTr++ ) {
    // std::cout << HLTFilters_[iTr] << std::endl;
    const pat::TriggerObjectStandAloneCollection mu1HLTMatches = muon1->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
    const pat::TriggerObjectStandAloneCollection mu2HLTMatches = muon2->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
    if (!mu1HLTMatches.empty() && !mu2HLTMatches.empty()) matched += (1<<iTr);
    // if (!mu1HLTMatches.empty() && !mu2HLTMatches.empty()) std::cout << std::endl << HLTFilters_[iTr] << std::endl;
  }

  return matched;
}

//
// constructors and destructor
//
DoubleDiMuonRootupler::DoubleDiMuonRootupler(const edm::ParameterSet& iConfig):
        doubledimuon_cand_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("doubledimuon_cand"))),
        doubledimuon_rf_cand_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("doubledimuon_rf_cand"))),
        thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"))),
        primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
        triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
	      isMC_(iConfig.getParameter<bool>("isMC")),
        OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
        HLTs_(iConfig.getParameter<std::vector<std::string>>("HLTs")),
        HLTFilters_(iConfig.getParameter<std::vector<std::string>>("filters"))
{
	      edm::Service<TFileService> fs;
        doubledimuon_tree = fs->make<TTree>("OniaPhiTree","Tree of Onia and Phi");

        doubledimuon_tree->Branch("run",                &run,                "run/I");
        doubledimuon_tree->Branch("event",              &event,              "event/I");
        doubledimuon_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");
        doubledimuon_tree->Branch("trigger",            &trigger,            "trigger/I");

        doubledimuon_tree->Branch("doubledimuon_p4",   "TLorentzVector", &doubledimuon_not_rf_p4);
        doubledimuon_tree->Branch("lowdim_p4",     "TLorentzVector", &lowdim_not_rf_p4);
        doubledimuon_tree->Branch("higdim_p4",     "TLorentzVector", &higdim_not_rf_p4);
        doubledimuon_tree->Branch("muonLowN_p4",   "TLorentzVector", &muonLowN_p4);
        doubledimuon_tree->Branch("muonHighN_p4",   "TLorentzVector", &muonHighN_p4);
        doubledimuon_tree->Branch("muonHighP_p4",   "TLorentzVector", &muonHighP_p4);
        doubledimuon_tree->Branch("muonLowN_p4",   "TLorentzVector", &muonLowN_p4);

        doubledimuon_tree->Branch("doubledimuon_rf_p4",   "TLorentzVector", &doubledimuon_rf_p4);
        doubledimuon_tree->Branch("lowdim_rf_p4",     "TLorentzVector", &lowdim_rf_p4);
        doubledimuon_tree->Branch("higdim_rf_p4",     "TLorentzVector", &higdim_rf_p4);
        doubledimuon_tree->Branch("muonLowN_rf_p4",   "TLorentzVector", &muonLowN_rf_p4);
        doubledimuon_tree->Branch("muonHighN_rf_p4",   "TLorentzVector", &muonHighN_rf_p4);
        doubledimuon_tree->Branch("muonHighP_rf_p4",   "TLorentzVector", &muonHighP_rf_p4);
        doubledimuon_tree->Branch("muonLowN_rf_p4",   "TLorentzVector", &muonLowN_rf_p4);

        doubledimuon_tree->Branch("doubledimuon_rf_bindx", &doubledimuon_rf_bindx, "doubledimuon_rf_bindx/I");

        doubledimuon_tree->Branch("higdim_vProb",        &higdim_vProb,        "higdim_vProb/D");
        doubledimuon_tree->Branch("higdim_vNChi2",       &higdim_vChi2,        "higdim_vNChi2/D");
        doubledimuon_tree->Branch("higdim_DCA",          &higdim_DCA,          "higdim_DCA/D");
        doubledimuon_tree->Branch("higdim_ctauPV",       &higdim_ctauPV,       "higdim_ctauPV/D");
        doubledimuon_tree->Branch("higdim_ctauErrPV",    &higdim_ctauErrPV,    "higdim_ctauErrPV/D");
        doubledimuon_tree->Branch("higdim_cosAlpha",     &higdim_cosAlpha,     "higdim_cosAlpha/D");
        doubledimuon_tree->Branch("higdim_triggerMatch", &higdim_triggerMatch, "higdim_triggerMatch/I");

        doubledimuon_tree->Branch("lowdim_vProb",        &lowdim_vProb,        "lowdim_vProb/D");
        doubledimuon_tree->Branch("lowdim_vNChi2",       &lowdim_vChi2,        "lowdim_vNChi2/D");
        doubledimuon_tree->Branch("lowdim_DCA",          &lowdim_DCA,          "lowdim_DCA/D");
        doubledimuon_tree->Branch("lowdim_ctauPV",       &lowdim_ctauPV,       "lowdim_ctauPV/D");
        doubledimuon_tree->Branch("lowdim_ctauErrPV",    &lowdim_ctauErrPV,    "lowdim_ctauErrPV/D");
        doubledimuon_tree->Branch("lowdim_cosAlpha",     &lowdim_cosAlpha,     "lowdim_cosAlpha/D");
        doubledimuon_tree->Branch("lowdim_triggerMatch", &lowdim_triggerMatch, "lowdim_triggerMatch/I");

        doubledimuon_tree->Branch("doubledimuon_vProb",      &doubledimuon_vProb,        "doubledimuon_vProb/D");
        doubledimuon_tree->Branch("doubledimuon_vChi2",      &doubledimuon_vChi2,        "doubledimuon_vChi2/D");
        doubledimuon_tree->Branch("doubledimuon_cosAlpha",   &doubledimuon_cosAlpha,     "doubledimuon_cosAlpha/D");
        doubledimuon_tree->Branch("doubledimuon_ctauPV",     &doubledimuon_ctauPV,       "doubledimuon_ctauPV/D");
        doubledimuon_tree->Branch("doubledimuon_ctauErrPV",  &doubledimuon_ctauErrPV,    "doubledimuon_ctauErrPV/D");
        doubledimuon_tree->Branch("doubledimuon_charge",     &doubledimuon_charge,       "doubledimuon_charge/I");

        doubledimuon_tree->Branch("highDiMM_fit",  &highDiMM_fit,    "highDiMM_fit/D");
        doubledimuon_tree->Branch("highDiMPx_fit",  &highDiMPx_fit,    "highDiMPx_fit/D");
        doubledimuon_tree->Branch("highDiMPy_fit",  &highDiMPy_fit,    "highDiMPy_fit/D");
        doubledimuon_tree->Branch("highDiMPz_fit",  &highDiMPz_fit,    "highDiMPz_fit/D");


        doubledimuon_tree->Branch("muonHighP_isLoose",        &muonHighP_isLoose,        "muonHighP_isLoose/O");
        doubledimuon_tree->Branch("muonHighP_isSoft",        &muonHighP_isSoft,        "muonHighP_isSoft/O");
        doubledimuon_tree->Branch("muonHighP_isMedium",        &muonHighP_isMedium,        "muonHighP_isMedium/O");
        doubledimuon_tree->Branch("muonHighP_isHighPt",        &muonHighP_isHighPt,        "muonHighP_isHighPt/O");

        doubledimuon_tree->Branch("muonHighP_isTracker",        &muonHighP_isTracker,        "muonHighP_isTracker/O");
        doubledimuon_tree->Branch("muonHighP_isGlobal",        &muonHighP_isGlobal,        "muonHighP_isGlobal/O");

        doubledimuon_tree->Branch("muonHighN_isLoose",        &muonHighN_isLoose,        "muonHighN_isLoose/O");
        doubledimuon_tree->Branch("muonHighN_isSoft",        &muonHighN_isSoft,        "muonHighN_isSoft/O");
        doubledimuon_tree->Branch("muonHighN_isMedium",        &muonHighN_isMedium,        "muonHighN_isMedium/O");
        doubledimuon_tree->Branch("muonHighN_isHighPt",        &muonHighN_isHighPt,        "muonHighN_isHighPt/O");

        doubledimuon_tree->Branch("muonHighN_isTracker",        &muonHighN_isTracker,        "muonHighN_isTracker/O");
        doubledimuon_tree->Branch("muonHighN_isGlobal",        &muonHighN_isGlobal,        "muonHighN_isGlobal/O");

        doubledimuon_tree->Branch("muonHighP_type",     &muonHighP_type,       "muonHighP_type/i");
        doubledimuon_tree->Branch("muonHighN_type",     &muonHighN_type,       "muonHighN_type/i");

        doubledimuon_tree->Branch("muonLowP_isLoose",        &muonLowP_isLoose,        "muonLowP_isLoose/O");
        doubledimuon_tree->Branch("muonLowP_isSoft",        &muonLowP_isSoft,        "muonLowP_isSoft/O");
        doubledimuon_tree->Branch("muonLowP_isMedium",        &muonLowP_isMedium,        "muonLowP_isMedium/O");
        doubledimuon_tree->Branch("muonLowP_isHighPt",        &muonLowP_isHighPt,        "muonLowP_isHighPt/O");

        doubledimuon_tree->Branch("muonLowP_isTracker",        &muonLowP_isTracker,        "muonLowP_isTracker/O");
        doubledimuon_tree->Branch("muonLowP_isGlobal",        &muonLowP_isGlobal,        "muonLowP_isGlobal/O");

        doubledimuon_tree->Branch("muonLowN_isLoose",        &muonLowN_isLoose,        "muonLowN_isLoose/O");
        doubledimuon_tree->Branch("muonLowN_isSoft",        &muonLowN_isSoft,        "muonLowN_isSoft/O");
        doubledimuon_tree->Branch("muonLowN_isMedium",        &muonLowN_isMedium,        "muonLowN_isMedium/O");
        doubledimuon_tree->Branch("muonLowN_isHighPt",        &muonLowN_isHighPt,        "muonLowN_isHighPt/O");

        doubledimuon_tree->Branch("muonLowN_isTracker",        &muonLowN_isTracker,        "muonLowN_isTracker/O");
        doubledimuon_tree->Branch("muonLowN_isGlobal",        &muonLowN_isGlobal,        "muonLowN_isGlobal/O");

        doubledimuon_tree->Branch("muonLowP_type",     &muonLowP_type,       "muonLowP_type/i");
        doubledimuon_tree->Branch("muonLowN_type",     &muonLowN_type,       "muonLowN_type/i");

        doubledimuon_tree->Branch("isBestCandidate",        &isBestCandidate,        "isBestCandidate/O");

	if(isMC_)
	  {
            doubledimuon_tree->Branch("gen_doubledimuon_pdgId", &gen_doubledimuon_pdgId, "gen_doubledimuon_pdgId/I");
      	    doubledimuon_tree->Branch("gen_doubledimuon_p4",    "TLorentzVector", &gen_doubledimuon_p4);
      	    doubledimuon_tree->Branch("gen_higdim_p4",      "TLorentzVector", &gen_higdim_p4);
      	    doubledimuon_tree->Branch("gen_lowdim_p4",      "TLorentzVector", &gen_lowdim_p4);
            doubledimuon_tree->Branch("gen_muonLowP_p4",    "TLorentzVector", &gen_muonLowP_p4);
            doubledimuon_tree->Branch("gen_muonHighN_p4",    "TLorentzVector", &gen_muonHighN_p4);
            doubledimuon_tree->Branch("gen_muonHighP_p4",    "TLorentzVector", &gen_muonHighP_p4);
            doubledimuon_tree->Branch("gen_muonLowN_p4",    "TLorentzVector", &gen_muonLowN_p4);
	  }
        genCands_ = consumes<reco::GenParticleCollection>((edm::InputTag)"prunedGenParticles");
        packCands_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");
}

DoubleDiMuonRootupler::~DoubleDiMuonRootupler() {}

//
// member functions
//

bool DoubleDiMuonRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   for (size_t i=0; i< particle->numberOfMothers(); i++) {
      if (isAncestor(ancestor, particle->mother(i))) return true;
   }
   return false;
}


// ------------ method called for each event  ------------
void DoubleDiMuonRootupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
//  using namespace edm;
  using namespace std;

  edm::Handle<std::vector<pat::CompositeCandidate>> doubledimuon_cand_handle;
  iEvent.getByToken(doubledimuon_cand_Label, doubledimuon_cand_handle);

  edm::Handle<std::vector<pat::CompositeCandidate>> doubledimuon_rf_cand_handle;
  iEvent.getByToken(doubledimuon_rf_cand_Label, doubledimuon_rf_cand_handle);

  edm::Handle<std::vector<reco::Vertex >> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( triggerResults_Label , triggerResults_handle);

  numPrimaryVertices = primaryVertices_handle->size();
  run = iEvent.id().run();
  event = iEvent.id().event();

  reco::Vertex thePrimaryV;
  reco::Vertex theBeamSpotV;

  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(thebeamspot_,theBeamSpot);
  reco::BeamSpot bs = *theBeamSpot;

  if ( primaryVertices_handle->begin() != primaryVertices_handle->end() ) {
    thePrimaryV = reco::Vertex(*(primaryVertices_handle->begin()));
  }
  else {
    thePrimaryV = reco::Vertex(bs.position(), bs.covariance3D());
  }

	// grab Trigger information
	// save it in variable trigger, trigger is an int between 0 and 7, in binary it is:
	// (pass 10)(pass 8)(pass 0)
	// ex. 7 = pass 0, 8 and 10
	// ex. 6 = pass 8, 10
        // ex. 1 = pass 0

  trigger = 0;

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

// grabbing doubledimuon information
  if (!doubledimuon_cand_handle.isValid()) std::cout<< "No doubledimuon information " << run << "," << event <<std::endl;
  if (!doubledimuon_rf_cand_handle.isValid()) std::cout<< "No doubledimuon_rf information " << run << "," << event <<std::endl;
// get rf information. Notice we are just keeping combinations with succesfull vertex fit
  if (doubledimuon_rf_cand_handle.isValid() && doubledimuon_cand_handle.isValid()) {



    pat::CompositeCandidate doubledimuon_rf_cand, doubledimuon_cand, *higdim_cand, *lowdim_cand, *higdim_cand_rf, *lowdim_cand_rf;

    for (unsigned int i=0; i< doubledimuon_rf_cand_handle->size(); i++){

      doubledimuon_rf_cand      = doubledimuon_rf_cand_handle->at(i);
      doubledimuon_rf_bindx     = doubledimuon_rf_cand.userInt("bIndex");

      if (doubledimuon_rf_bindx<0 || doubledimuon_rf_bindx>(int) doubledimuon_cand_handle->size()) {
        std::cout << "Incorrect index for oniatt combination " << run << "," << event <<"," << doubledimuon_rf_bindx << std::endl;
        continue;
      }

      doubledimuon_vProb     = doubledimuon_rf_cand.userFloat("vProb");
      doubledimuon_vChi2     = doubledimuon_rf_cand.userFloat("vChi2");
      doubledimuon_cosAlpha  = doubledimuon_rf_cand.userFloat("cosAlpha");
      doubledimuon_ctauPV    = doubledimuon_rf_cand.userFloat("ctauPV");
      doubledimuon_ctauErrPV = doubledimuon_rf_cand.userFloat("ctauErrPV");
      highDiMM_fit = doubledimuon_rf_cand.userFloat("highDiMM_fit");
      highDiMPx_fit = doubledimuon_rf_cand.userFloat("highDiMPx_fit");
      highDiMPy_fit = doubledimuon_rf_cand.userFloat("highDiMPy_fit");
      highDiMPz_fit = doubledimuon_rf_cand.userFloat("highDiMPz_fit");

      doubledimuon_rf_p4.SetPtEtaPhiM(doubledimuon_rf_cand.pt(),doubledimuon_rf_cand.eta(),doubledimuon_rf_cand.phi(),doubledimuon_rf_cand.mass());
      higdim_rf_p4.SetPtEtaPhiM(doubledimuon_rf_cand.daughter("higdimuon")->pt(),doubledimuon_rf_cand.daughter("higdimuon")->eta(),
                              doubledimuon_rf_cand.daughter("higdimuon")->phi(),doubledimuon_rf_cand.daughter("higdimuon")->mass());
      lowdim_rf_p4.SetPtEtaPhiM(doubledimuon_rf_cand.daughter("lowdimuon")->pt(),doubledimuon_rf_cand.daughter("lowdimuon")->eta(),
                              doubledimuon_rf_cand.daughter("lowdimuon")->phi(),doubledimuon_rf_cand.daughter("lowdimuon")->mass());


      higdim_cand_rf = dynamic_cast <pat::CompositeCandidate *>(doubledimuon_rf_cand.daughter("higdimuon"));
      lowdim_cand_rf = dynamic_cast <pat::CompositeCandidate *>(doubledimuon_rf_cand.daughter("lowdimuon"));

      reco::Candidate::LorentzVector vJpsiP = higdim_cand_rf->daughter("muon1")->p4();
      reco::Candidate::LorentzVector vJpsiM = higdim_cand_rf->daughter("muon2")->p4();

      if (higdim_cand_rf->daughter("muon1")->charge() < 0) {
         vJpsiP = higdim_cand_rf->daughter("muon2")->p4();
         vJpsiM = higdim_cand_rf->daughter("muon1")->p4();
      }

      muonLowN_rf_p4.SetPtEtaPhiM(vJpsiP.pt(), vJpsiP.eta(), vJpsiP.phi(), vJpsiP.mass());
      muonHighN_rf_p4.SetPtEtaPhiM(vJpsiM.pt(), vJpsiM.eta(), vJpsiM.phi(), vJpsiM.mass());

      reco::Candidate::LorentzVector vPhiP = lowdim_cand_rf->daughter("muon1")->p4();
      reco::Candidate::LorentzVector vPhiM = lowdim_cand_rf->daughter("muon2")->p4();

      if (lowdim_cand_rf->daughter("muon1")->charge() < 0) {
         vPhiP = lowdim_cand_rf->daughter("muon2")->p4();
         vPhiM = lowdim_cand_rf->daughter("muon1")->p4();
      }

      pat::CompositeCandidate doubledimuon_not_rf_cand = doubledimuon_cand_handle->at(doubledimuon_rf_bindx);

      higdim_cand = dynamic_cast <pat::CompositeCandidate *>(doubledimuon_not_rf_cand.daughter("higdimuon"));
      lowdim_cand = dynamic_cast <pat::CompositeCandidate *>(doubledimuon_not_rf_cand.daughter("lowdimuon"));

      doubledimuon_not_rf_p4.SetPtEtaPhiM(doubledimuon_not_rf_cand.pt(),doubledimuon_not_rf_cand.eta(),doubledimuon_not_rf_cand.phi(),doubledimuon_not_rf_cand.mass());
      higdim_not_rf_p4.SetPtEtaPhiM(doubledimuon_not_rf_cand.daughter("higdimuon")->pt(),doubledimuon_not_rf_cand.daughter("higdimuon")->eta(),
                              doubledimuon_not_rf_cand.daughter("higdimuon")->phi(),doubledimuon_not_rf_cand.daughter("higdimuon")->mass());
      lowdim_not_rf_p4.SetPtEtaPhiM(doubledimuon_not_rf_cand.daughter("lowdimuon")->pt(),doubledimuon_not_rf_cand.daughter("lowdimuon")->eta(),
                              doubledimuon_not_rf_cand.daughter("lowdimuon")->phi(),doubledimuon_not_rf_cand.daughter("lowdimuon")->mass());

      higdim_vProb        = higdim_cand->userFloat("vProb");
      higdim_vChi2        = higdim_cand->userFloat("vNChi2");
      higdim_DCA          = higdim_cand->userFloat("DCA");
      higdim_ctauPV       = higdim_cand->userFloat("ppdlPV");
      higdim_ctauErrPV    = higdim_cand->userFloat("ppdlErrPV");
      higdim_cosAlpha     = higdim_cand->userFloat("cosAlpha");
      higdim_triggerMatch = DoubleDiMuonRootupler::isTriggerMatched(higdim_cand);

      lowdim_vProb        = lowdim_cand->userFloat("vProb");
      lowdim_vChi2        = lowdim_cand->userFloat("vNChi2");
      lowdim_DCA          = lowdim_cand->userFloat("DCA");
      lowdim_ctauPV       = lowdim_cand->userFloat("ppdlPV");
      lowdim_ctauErrPV    = lowdim_cand->userFloat("ppdlErrPV");
      lowdim_cosAlpha     = lowdim_cand->userFloat("cosAlpha");
      lowdim_triggerMatch = DoubleDiMuonRootupler::isTriggerMatched(lowdim_cand);

      const pat::Muon *highPatMuonP,  *highPatMuonN, *lowPatMuonP, *lowPatMuonN;

      if (higdim_cand->daughter("muon1")->charge() < 0) {
         vJpsiP = higdim_cand->daughter("muon2")->p4();
         vJpsiM = higdim_cand->daughter("muon1")->p4();
         highPatMuonN = dynamic_cast<const pat::Muon*>(higdim_cand->daughter("muon1"));
         highPatMuonP = dynamic_cast<const pat::Muon*>(higdim_cand->daughter("muon2"));
      } else
      {
        highPatMuonP = dynamic_cast<const pat::Muon*>(higdim_cand->daughter("muon1"));
        highPatMuonN = dynamic_cast<const pat::Muon*>(higdim_cand->daughter("muon2"));
      }

      muonHighP_p4.SetPtEtaPhiM(vJpsiP.pt(), vJpsiP.eta(), vJpsiP.phi(), vJpsiP.mass());
      muonHighN_p4.SetPtEtaPhiM(vJpsiM.pt(), vJpsiM.eta(), vJpsiM.phi(), vJpsiM.mass());

      muonHighP_isLoose   =  highPatMuonP->isLooseMuon();
      muonHighP_isSoft   =  highPatMuonP->isSoftMuon(thePrimaryV);
      muonHighP_isMedium   = highPatMuonP->isMediumMuon();
      muonHighP_isHighPt   = highPatMuonP->isHighPtMuon(thePrimaryV);
      muonHighP_isTracker   = highPatMuonP->isTrackerMuon();
      muonHighP_isGlobal   = highPatMuonP->isGlobalMuon();
      muonHighN_isLoose   = highPatMuonN->isLooseMuon();
      muonHighN_isSoft   = highPatMuonN->isSoftMuon(thePrimaryV);
      muonHighN_isMedium   = highPatMuonN->isMediumMuon();
      muonHighN_isHighPt   = highPatMuonN->isHighPtMuon(thePrimaryV);
      muonHighN_isTracker   = highPatMuonN->isTrackerMuon();
      muonHighN_isGlobal   = highPatMuonN->isGlobalMuon();
      muonHighP_type   =  highPatMuonP->type();
      muonHighN_type   = highPatMuonN->type();

      //double kmass = 0.4936770;
      doubledimuon_p4.SetPtEtaPhiM(doubledimuon_cand.pt(),doubledimuon_cand.eta(),doubledimuon_cand.phi(),doubledimuon_cand.mass());
      higdim_p4.SetPtEtaPhiM(higdim_cand->pt(),higdim_cand->eta(),higdim_cand->phi(),higdim_cand->mass());
      lowdim_p4.SetPtEtaPhiM(lowdim_cand->pt(), lowdim_cand->eta(), lowdim_cand->phi(), lowdim_cand->mass());

      vPhiP = lowdim_cand->daughter("muon1")->p4();
      vPhiM = lowdim_cand->daughter("muon2")->p4();

      if (lowdim_cand->daughter("muon1")->charge() < 0) {
         vPhiP = lowdim_cand->daughter("muon2")->p4();
         vPhiM = lowdim_cand->daughter("muon1")->p4();
         lowPatMuonN = dynamic_cast<const pat::Muon*>(lowdim_cand->daughter("muon1"));
         lowPatMuonP = dynamic_cast<const pat::Muon*>(lowdim_cand->daughter("muon2"));
      } else
      {
        lowPatMuonP = dynamic_cast<const pat::Muon*>(lowdim_cand->daughter("muon1"));
        lowPatMuonN = dynamic_cast<const pat::Muon*>(lowdim_cand->daughter("muon2"));
      }

      muonLowP_isLoose   =  lowPatMuonP->isLooseMuon();
      muonLowP_isSoft   =  lowPatMuonP->isSoftMuon(thePrimaryV);
      muonLowP_isMedium   = lowPatMuonP->isMediumMuon();
      muonLowP_isHighPt   = lowPatMuonP->isHighPtMuon(thePrimaryV);
      muonLowP_isTracker   = lowPatMuonP->isTrackerMuon();
      muonLowP_isGlobal   = lowPatMuonP->isGlobalMuon();
      muonLowN_isLoose   = lowPatMuonN->isLooseMuon();
      muonLowN_isSoft   = lowPatMuonN->isSoftMuon(thePrimaryV);
      muonLowN_isMedium   = lowPatMuonN->isMediumMuon();
      muonLowN_isHighPt   = lowPatMuonN->isHighPtMuon(thePrimaryV);
      muonLowN_isTracker   = lowPatMuonN->isTrackerMuon();
      muonLowN_isGlobal   = lowPatMuonN->isGlobalMuon();
      muonLowP_type   =  lowPatMuonP->type();
      muonLowN_type   = lowPatMuonN->type();

      muonLowP_p4.SetPtEtaPhiM(vPhiP.pt(), vPhiP.eta(), vPhiP.phi(), vPhiP.mass());
      muonLowN_p4.SetPtEtaPhiM(vPhiM.pt(), vPhiM.eta(), vPhiM.phi(), vPhiM.mass());

      doubledimuon_tree->Fill();

      if (OnlyBest_) break;
      else if(i==0)
      isBestCandidate = false;

      }

  }

}

// ------------ method called once each job just before starting event loop  ------------
void DoubleDiMuonRootupler::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void DoubleDiMuonRootupler::endJob() {}

// ------------ method called when starting to processes a run  ------------
void DoubleDiMuonRootupler::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void DoubleDiMuonRootupler::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void DoubleDiMuonRootupler::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void DoubleDiMuonRootupler::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DoubleDiMuonRootupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DoubleDiMuonRootupler);
