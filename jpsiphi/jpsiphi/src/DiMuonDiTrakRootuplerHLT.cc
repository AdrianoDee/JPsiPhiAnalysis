/*
   Package:    DiMuonDiTrakRootuplerHLT
   Class:      DiMuonDiTrakRootuplerHLT

   Description: make rootuple of DiMuon-DiTrack combination

   Original Author: Adriano Di Florio
   Created:  based on Alberto Sanchez Hernandez PsiTrkTrk Code

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
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

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

class DiMuonDiTrakRootuplerHLT : public edm::EDAnalyzer {
   public:
      explicit DiMuonDiTrakRootuplerHLT(const edm::ParameterSet&);
      ~DiMuonDiTrakRootuplerHLT() override;

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
  edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuonditrk_cand_Label;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuonditrk_rf_cand_Label;
  edm::EDGetTokenT<reco::BeamSpot> thebeamspot_;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  int  dimuonditrk_pdgid_, dimuon_pdgid_, trak_pdgid_;
  bool isMC_,OnlyBest_;
  std::vector<std::string>  HLTs_;
  std::vector<std::string>  HLTFilters_;
  std::string treeName_;

  UInt_t run,event,numPrimaryVertices, trigger;

  TLorentzVector dimuonditrk_p4;
  TLorentzVector dimuon_p4;
  TLorentzVector ditrak_p4;
  TLorentzVector muonp_p4;
  TLorentzVector muonn_p4;
  TLorentzVector kaonp_p4;
  TLorentzVector kaonn_p4;

  TLorentzVector dimuonditrkTrigger_p4;
  TLorentzVector dimuonTrigger_p4;
  TLorentzVector ditrakTrigger_p4;

  TLorentzVector dimuonditrk_rf_p4;
  TLorentzVector dimuonditrk_not_rf_p4;
  TLorentzVector dimuon_rf_p4, dimuon_not_rf_p4;
  TLorentzVector ditrak_rf_p4, ditrak_not_rf_p4;
  TLorentzVector muonp_rf_p4;
  TLorentzVector muonn_rf_p4;
  TLorentzVector kaonp_rf_p4;
  TLorentzVector kaonn_rf_p4;

  Int_t dimuonditrk_charge;

  UInt_t dimuon_triggerMatch, dimuon_triggerMatch_rf;
  UInt_t muonN_tMatch, muonP_tMatch, trakN_tMatch, trakP_tMatch;

  Double_t dimuonditrk_vProb,  dimuonditrk_vChi2, dimuonditrk_cosAlpha, dimuonditrk_ctauPV, dimuonditrk_ctauErrPV;

  Double_t dimuon_vProb, dimuon_vChi2, dimuon_DCA, dimuon_ctauPV, dimuon_ctauErrPV, dimuon_cosAlpha;

  Bool_t muonP_isLoose, muonP_isSoft, muonP_isMedium, muonP_isHighPt;
  Bool_t muonN_isLoose, muonN_isSoft, muonN_isMedium, muonN_isHighPt;

  Bool_t muonP_isTracker, muonP_isGlobal, muonN_isTracker, muonN_isGlobal;
  UInt_t muonP_type, muonN_type;

  Bool_t muonP_rf_isLoose, muonP_rf_isSoft, muonP_rf_isMedium, muonP_rf_isHighPt;
  Bool_t muonN_rf_isLoose, muonN_rf_isSoft, muonN_rf_isMedium, muonN_rf_isHighPt;

  Bool_t muonP_rf_isTracker, muonP_rf_isGlobal, muonN_rf_isTracker, muonN_rf_isGlobal;
  UInt_t muonP_rf_type, muonN_rf_type;

  Double_t track_KP_d0, track_KP_d0Err, track_KP_dz, track_KP_dxy;
  Int_t track_KP_nvsh, track_KP_nvph;

  Int_t track_KN_nvsh, track_KN_nvph;

  Int_t dimuonditrk_rf_bindx;

  Int_t noXCandidates;

  Bool_t isBestCandidate;

  Int_t          gen_dimuonditrk_pdgId;
  TLorentzVector gen_dimuonditrk_p4;
  TLorentzVector gen_dimuon_p4;
  TLorentzVector gen_ditrak_p4;
  TLorentzVector gen_muonp_p4;
  TLorentzVector gen_muonn_p4;
  TLorentzVector gen_kaonp_p4;
  TLorentzVector gen_kaonn_p4;

  TTree* dimuonditrk_tree, *dimuonditrk_tree_rf;
  edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;
};

UInt_t DiMuonDiTrakRootuplerHLT::isTriggerMatched(pat::CompositeCandidate *diMuon_cand) {
  const pat::Muon* muonN = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muonN"));
  const pat::Muon* muonP = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muonP"));
  UInt_t matched = 0;  // if no list is given, is not matched

  // if matched a given trigger, set the bit, in the same order as listed
  for (unsigned int iTr = 0; iTr<HLTFilters_.size(); iTr++ ) {
    // std::cout << HLTFilters_[iTr] << std::endl;
    const pat::TriggerObjectStandAloneCollection mu1HLTMatches = muonN->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
    const pat::TriggerObjectStandAloneCollection mu2HLTMatches = muonP->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
    if (!mu1HLTMatches.empty() && !mu2HLTMatches.empty()) matched += (1<<iTr);
    // if (!mu1HLTMatches.empty() && !mu2HLTMatches.empty()) std::cout << std::endl << HLTFilters_[iTr] << std::endl;
  }

  return matched;
}

//
// constants, enums and typedefs
//

//
// static data member definitions
//
static const Double_t psi1SMass =  3.09691;


//
// constructors and destructor
//
DiMuonDiTrakRootuplerHLT::DiMuonDiTrakRootuplerHLT(const edm::ParameterSet& iConfig):
        dimuonditrk_cand_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("dimuonditrk_cand"))),
        dimuonditrk_rf_cand_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("dimuonditrk_rf_cand"))),
        thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"))),
        primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
        triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
	      isMC_(iConfig.getParameter<bool>("isMC")),
        OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
        HLTs_(iConfig.getParameter<std::vector<std::string>>("HLTs")),
        HLTFilters_(iConfig.getParameter<std::vector<std::string>>("Filters")),
        treeName_(iConfig.getParameter<std::string>("TreeName"))
{
	      edm::Service<TFileService> fs;
        dimuonditrk_tree = fs->make<TTree>(treeName_.data(),"Tree of DiMuon and DiTrak");

        dimuonditrk_tree->Branch("run",                &run,                "run/I");
        dimuonditrk_tree->Branch("event",              &event,              "event/I");
        dimuonditrk_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");
        dimuonditrk_tree->Branch("trigger",            &trigger,            "trigger/I");
        dimuonditrk_tree->Branch("noXCandidates",      &noXCandidates,      "noXCandidates/I");

        //p4s
        dimuonditrk_tree->Branch("dimuonditrk_p4",   "TLorentzVector", &dimuonditrk_p4);
        dimuonditrk_tree->Branch("ditrak_p4",     "TLorentzVector", &ditrak_p4);
        dimuonditrk_tree->Branch("dimuon_p4",     "TLorentzVector", &dimuon_p4);
        dimuonditrk_tree->Branch("muonp_p4",   "TLorentzVector", &muonp_p4);
        dimuonditrk_tree->Branch("muonn_p4",   "TLorentzVector", &muonn_p4);
        dimuonditrk_tree->Branch("kaonp_p4",   "TLorentzVector", &kaonp_p4);
        dimuonditrk_tree->Branch("kaonn_p4",   "TLorentzVector", &kaonn_p4);

        //Trigger
        dimuonditrk_tree->Branch("dimuonditrkTrigger_p4",   "TLorentzVector", &dimuonditrkTrigger_p4);
        dimuonditrk_tree->Branch("ditrakTrigger_p4",     "TLorentzVector", &ditrakTrigger_p4);
        dimuonditrk_tree->Branch("dimuonTrigger_p4",     "TLorentzVector", &dimuonTrigger_p4);

        //refitted p4s
        dimuonditrk_tree->Branch("dimuonditrk_rf_p4",   "TLorentzVector", &dimuonditrk_rf_p4);
        dimuonditrk_tree->Branch("ditrak_rf_p4",     "TLorentzVector", &ditrak_rf_p4);
        dimuonditrk_tree->Branch("dimuon_rf_p4",     "TLorentzVector", &dimuon_rf_p4);
        dimuonditrk_tree->Branch("muonp_rf_p4",   "TLorentzVector", &muonp_rf_p4);
        dimuonditrk_tree->Branch("muonn_rf_p4",   "TLorentzVector", &muonn_rf_p4);
        dimuonditrk_tree->Branch("kaonp_rf_p4",   "TLorentzVector", &kaonp_rf_p4);
        dimuonditrk_tree->Branch("kaonn_rf_p4",   "TLorentzVector", &kaonn_rf_p4);

        // //2mu vertexing
        dimuonditrk_tree->Branch("dimuon_vProb",        &dimuon_vProb,        "dimuon_vProb/D");
        dimuonditrk_tree->Branch("dimuon_vNChi2",       &dimuon_vChi2,        "dimuon_vNChi2/D");
        dimuonditrk_tree->Branch("dimuon_DCA",          &dimuon_DCA,          "dimuon_DCA/D");
        dimuonditrk_tree->Branch("dimuon_ctauPV",       &dimuon_ctauPV,       "dimuon_ctauPV/D");
        dimuonditrk_tree->Branch("dimuon_ctauErrPV",    &dimuon_ctauErrPV,    "dimuon_ctauErrPV/D");
        dimuonditrk_tree->Branch("dimuon_cosAlpha",     &dimuon_cosAlpha,     "dimuon_cosAlpha/D");
        dimuonditrk_tree->Branch("dimuon_triggerMatch", &dimuon_triggerMatch, "dimuon_triggerMatch/I");
        //
        //2mu+2Trk vertexing
        dimuonditrk_tree->Branch("dimuonditrk_vProb",      &dimuonditrk_vProb,        "dimuonditrk_vProb/D");
        dimuonditrk_tree->Branch("dimuonditrk_vChi2",      &dimuonditrk_vChi2,        "dimuonditrk_vChi2/D");
        dimuonditrk_tree->Branch("dimuonditrk_cosAlpha",   &dimuonditrk_cosAlpha,     "dimuonditrk_cosAlpha/D");
        dimuonditrk_tree->Branch("dimuonditrk_ctauPV",     &dimuonditrk_ctauPV,       "dimuonditrk_ctauPV/D");
        dimuonditrk_tree->Branch("dimuonditrk_ctauErrPV",  &dimuonditrk_ctauErrPV,    "dimuonditrk_ctauErrPV/D");
        dimuonditrk_tree->Branch("dimuonditrk_charge",     &dimuonditrk_charge,       "dimuonditrk_charge/I");

        // //Muon flags
        dimuonditrk_tree->Branch("muonP_isLoose",        &muonP_isLoose,        "muonP_isLoose/O");
        dimuonditrk_tree->Branch("muonP_isSoft",        &muonP_isSoft,        "muonP_isSoft/O");
        dimuonditrk_tree->Branch("muonP_isMedium",        &muonP_isMedium,        "muonP_isMedium/O");
        dimuonditrk_tree->Branch("muonP_isHighPt",        &muonP_isHighPt,        "muonP_isHighPt/O");
        //
        dimuonditrk_tree->Branch("muonP_isTracker",        &muonP_isTracker,        "muonP_isTracker/O");
        dimuonditrk_tree->Branch("muonP_isGlobal",        &muonP_isGlobal,        "muonP_isGlobal/O");

        dimuonditrk_tree->Branch("muonN_isLoose",        &muonN_isLoose,        "muonN_isLoose/O");
        dimuonditrk_tree->Branch("muonN_isSoft",        &muonN_isSoft,        "muonN_isSoft/O");
        dimuonditrk_tree->Branch("muonN_isMedium",        &muonN_isMedium,        "muonN_isMedium/O");
        dimuonditrk_tree->Branch("muonN_isHighPt",        &muonN_isHighPt,        "muonN_isHighPt/O");

        dimuonditrk_tree->Branch("muonN_isTracker",        &muonN_isTracker,        "muonN_isTracker/O");
        dimuonditrk_tree->Branch("muonN_isGlobal",        &muonN_isGlobal,        "muonN_isGlobal/O");

        dimuonditrk_tree->Branch("muonP_type",     &muonP_type,       "muonP_type/i");
        dimuonditrk_tree->Branch("muonN_type",     &muonN_type,       "muonN_type/i");

        dimuonditrk_tree->Branch("isBestCandidate",        &isBestCandidate,        "isBestCandidate/O");

        dimuonditrk_tree->Branch("muonN_tMatch", &muonN_tMatch, "muonN_tMatch/I");
        dimuonditrk_tree->Branch("muonP_tMatch", &muonP_tMatch, "muonP_tMatch/I");
        dimuonditrk_tree->Branch("trakN_tMatch", &trakN_tMatch, "trakN_tMatch/I");
        dimuonditrk_tree->Branch("trakP_tMatch", &trakP_tMatch, "trakP_tMatch/I");
}

DiMuonDiTrakRootuplerHLT::~DiMuonDiTrakRootuplerHLT() {}

//
// member functions
//

bool DiMuonDiTrakRootuplerHLT::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   for (size_t i=0; i< particle->numberOfMothers(); i++) {
      if (isAncestor(ancestor, particle->mother(i))) return true;
   }
   return false;
}


// ------------ method called for each event  ------------
void DiMuonDiTrakRootuplerHLT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
//  using namespace edm;
  using namespace std;

  edm::Handle<std::vector<pat::CompositeCandidate>> dimuonditrk_cand_handle;
  iEvent.getByToken(dimuonditrk_cand_Label, dimuonditrk_cand_handle);

  edm::Handle<std::vector<pat::CompositeCandidate>> dimuonditrk_rf_cand_handle;
  iEvent.getByToken(dimuonditrk_rf_cand_Label, dimuonditrk_rf_cand_handle);

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

  isBestCandidate = true;

// grabbing dimuontt information
  if (!dimuonditrk_cand_handle.isValid()) std::cout<< "No dimuontt information " << run << "," << event <<std::endl;
  if (!dimuonditrk_rf_cand_handle.isValid()) std::cout<< "No dimuonditrk_rf information " << run << "," << event <<std::endl;
// get rf information. Notice we are just keeping combinations with succesfull vertex fit
  if (dimuonditrk_rf_cand_handle.isValid() && dimuonditrk_cand_handle.isValid()) {

    pat::CompositeCandidate dimuonditrk_rf_cand, dimuonditrk_cand, *dimuon_cand, *ditrak_cand, *dimuon_cand_rf, *ditrak_cand_rf;
    pat::CompositeCandidate *dimuonditrkTrigger_cand, *dimuonTrigger_cand, *ditrakTrigger_cand;

    noXCandidates = (Int_t)(dimuonditrk_rf_cand_handle->size());
    //Refitted Handle
    for (unsigned int i=0; i< dimuonditrk_rf_cand_handle->size(); i++){


      dimuonditrk_rf_cand   = dimuonditrk_rf_cand_handle->at(i);
      dimuonditrk_rf_bindx = dimuonditrk_rf_cand.userInt("bIndex");

      if (dimuonditrk_rf_bindx<0 || dimuonditrk_rf_bindx>(int) dimuonditrk_cand_handle->size()) {
        std::cout << "Incorrect index for dimuontt combination " << run << "," << event <<"," << dimuonditrk_rf_bindx << std::endl;
        continue;
      }

      dimuonditrk_vProb     = dimuonditrk_rf_cand.userFloat("vProb");
      dimuonditrk_vChi2     = dimuonditrk_rf_cand.userFloat("vChi2");
      dimuonditrk_cosAlpha  = dimuonditrk_rf_cand.userFloat("cosAlpha");
      dimuonditrk_ctauPV    = dimuonditrk_rf_cand.userFloat("ctauPV");
      dimuonditrk_ctauErrPV = dimuonditrk_rf_cand.userFloat("ctauErrPV");
      dimuonditrk_charge    = dimuonditrk_cand.charge();

      dimuonditrk_rf_p4.SetPtEtaPhiM(dimuonditrk_rf_cand.pt(),dimuonditrk_rf_cand.eta(),dimuonditrk_rf_cand.phi(),dimuonditrk_rf_cand.mass());
      dimuon_rf_p4.SetPtEtaPhiM(dimuonditrk_rf_cand.daughter("dimuon")->pt(),dimuonditrk_rf_cand.daughter("dimuon")->eta(),
                              dimuonditrk_rf_cand.daughter("dimuon")->phi(),dimuonditrk_rf_cand.daughter("dimuon")->mass());
      ditrak_rf_p4.SetPtEtaPhiM(dimuonditrk_rf_cand.daughter("ditrak")->pt(),dimuonditrk_rf_cand.daughter("ditrak")->eta(),
                              dimuonditrk_rf_cand.daughter("ditrak")->phi(),dimuonditrk_rf_cand.daughter("ditrak")->mass());

      dimuon_cand_rf = dynamic_cast <pat::CompositeCandidate *>(dimuonditrk_rf_cand.daughter("dimuon"));
      ditrak_cand_rf = dynamic_cast <pat::CompositeCandidate *>(dimuonditrk_rf_cand.daughter("ditrak"));

      reco::Candidate::LorentzVector vN = dimuon_cand_rf->daughter("muonN")->p4();
      reco::Candidate::LorentzVector VP = dimuon_cand_rf->daughter("muonP")->p4();

      muonp_rf_p4.SetPtEtaPhiM(vN.pt(), vN.eta(), vN.phi(), vN.mass());
      muonn_rf_p4.SetPtEtaPhiM(VP.pt(), VP.eta(), VP.phi(), VP.mass());

      reco::Candidate::LorentzVector kP = ditrak_cand_rf->daughter("trakP")->p4();
      reco::Candidate::LorentzVector kM = ditrak_cand_rf->daughter("trakN")->p4();

      kaonp_rf_p4.SetPtEtaPhiM(kP.pt(), kP.eta(), kP.phi(), kP.mass());
      kaonn_rf_p4.SetPtEtaPhiM(kM.pt(), kM.eta(), kM.phi(), kM.mass());

      //unref corresponding

      dimuonditrk_cand = dimuonditrk_cand_handle->at(dimuonditrk_rf_bindx);
      dimuon_cand = dynamic_cast <pat::CompositeCandidate *>(dimuonditrk_cand.daughter("dimuon"));
      ditrak_cand = dynamic_cast <pat::CompositeCandidate *>(dimuonditrk_cand.daughter("ditrak"));

      dimuonditrkTrigger_cand = dynamic_cast <pat::CompositeCandidate *>(dimuonditrk_cand.daughter("dimuonTTTrigger"));
      dimuonTrigger_cand = dynamic_cast <pat::CompositeCandidate *>(dimuonditrkTrigger_cand->daughter("dimuon"));
      ditrakTrigger_cand = dynamic_cast <pat::CompositeCandidate *>(dimuonditrkTrigger_cand->daughter("ditrak"));

      const pat::Muon *muonP, *muonN;

      if (dimuon_cand_rf->daughter("muonN")->charge() < 0) {
         vN = dimuon_cand->daughter("muonP")->p4();
         VP = dimuon_cand->daughter("muonN")->p4();
         muonN = dynamic_cast<const pat::Muon*>(dimuon_cand->daughter("muonN"));
         muonP = dynamic_cast<const pat::Muon*>(dimuon_cand->daughter("muonP"));
      } else
      {
        muonP = dynamic_cast<const pat::Muon*>(dimuon_cand->daughter("muonN"));
        muonN = dynamic_cast<const pat::Muon*>(dimuon_cand->daughter("muonP"));
      }

      muonP_isLoose    =  muonP->isLooseMuon();
      muonP_isSoft     =  muonP->isSoftMuon(thePrimaryV);
      muonP_isMedium   = muonP->isMediumMuon();
      muonP_isHighPt   = muonP->isHighPtMuon(thePrimaryV);
      muonP_isTracker  = muonP->isTrackerMuon();
      muonP_isGlobal   = muonP->isGlobalMuon();
      muonN_isLoose    = muonN->isLooseMuon();
      muonN_isSoft     = muonN->isSoftMuon(thePrimaryV);
      muonN_isMedium   = muonN->isMediumMuon();
      muonN_isHighPt   = muonN->isHighPtMuon(thePrimaryV);
      muonN_isTracker  = muonN->isTrackerMuon();
      muonN_isGlobal   = muonN->isGlobalMuon();
      muonP_type       = muonP->type();
      muonN_type       = muonN->type();

      muonp_p4.SetPtEtaPhiM(vN.pt(), vN.eta(), vN.phi(), vN.mass());
      muonn_p4.SetPtEtaPhiM(VP.pt(), VP.eta(), VP.phi(), VP.mass());

      kP = ditrak_cand->daughter("trakP")->p4();
      kM = ditrak_cand->daughter("trakN")->p4();

      kaonp_p4.SetPtEtaPhiM(kP.pt(), kP.eta(), kP.phi(), kP.mass());
      kaonn_p4.SetPtEtaPhiM(kM.pt(), kM.eta(), kM.phi(), kM.mass());

      //double kmass = 0.4936770;
      dimuonditrk_p4.SetPtEtaPhiM(dimuonditrk_cand.pt(),dimuonditrk_cand.eta(),dimuonditrk_cand.phi(),dimuonditrk_cand.mass());
      dimuon_p4.SetPtEtaPhiM(dimuon_cand->pt(),dimuon_cand->eta(),dimuon_cand->phi(),dimuon_cand->mass());
      ditrak_p4.SetPtEtaPhiM(ditrak_cand->pt(), ditrak_cand->eta(), ditrak_cand->phi(), ditrak_cand->mass());

      dimuonditrkTrigger_p4.SetPtEtaPhiM(dimuonditrkTrigger_cand->pt(),dimuonditrkTrigger_cand->eta(),dimuonditrkTrigger_cand->phi(),dimuonditrkTrigger_cand->mass());
      dimuonTrigger_p4.SetPtEtaPhiM(dimuonTrigger_cand->pt(),dimuonTrigger_cand->eta(),dimuonTrigger_cand->phi(),dimuonTrigger_cand->mass());
      ditrakTrigger_p4.SetPtEtaPhiM(ditrakTrigger_cand->pt(),ditrakTrigger_cand->eta(),ditrakTrigger_cand->phi(),ditrakTrigger_cand->mass());

      dimuon_vProb        = dimuon_cand->userFloat("vProb");
      dimuon_vChi2        = dimuon_cand->userFloat("vNChi2");
      dimuon_DCA          = dimuon_cand->userFloat("DCA");
      dimuon_ctauPV       = dimuon_cand->userFloat("ppdlPV");
      dimuon_ctauErrPV    = dimuon_cand->userFloat("ppdlErrPV");
      dimuon_cosAlpha     = dimuon_cand->userFloat("cosAlpha");
      dimuon_triggerMatch = DiMuonDiTrakRootuplerHLT::isTriggerMatched(dimuon_cand);

      muonP_tMatch  = dimuon_cand->userInt("tMatchP");
      muonN_tMatch  = dimuon_cand->userInt("tMatchN");
      trakP_tMatch  = dimuonditrk_cand.userInt("trakMatchP");
      trakN_tMatch  = dimuonditrk_cand.userInt("trakMatchN");

      dimuonditrk_tree->Fill();

      if (OnlyBest_) break;
      else if(i==0)
      isBestCandidate = false;

        // dimuontt candidates are sorted by vProb
    }

  }

}

// ------------ method called once each job just before starting event loop  ------------
void DiMuonDiTrakRootuplerHLT::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void DiMuonDiTrakRootuplerHLT::endJob() {}

// ------------ method called when starting to processes a run  ------------
void DiMuonDiTrakRootuplerHLT::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void DiMuonDiTrakRootuplerHLT::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void DiMuonDiTrakRootuplerHLT::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void DiMuonDiTrakRootuplerHLT::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DiMuonDiTrakRootuplerHLT::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonDiTrakRootuplerHLT);
