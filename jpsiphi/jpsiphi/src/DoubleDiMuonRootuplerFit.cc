/*
   Package:    DoubleDiMuonRootuplerFit
   Class:      DoubleDiMuonRootuplerFit

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


class DoubleDiMuonRootuplerFit : public edm::EDAnalyzer {
   public:
      explicit DoubleDiMuonRootuplerFit(const edm::ParameterSet&);
      ~DoubleDiMuonRootuplerFit() override;

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
  bool isMC_,OnlyBest_,OnlyGen_;
  UInt_t motherpdgid_;
  std::vector<std::string>                            HLTs_;
  std::vector<std::string>                            HLTFilters_;

  UInt_t run,event,numPrimaryVertices, trigger;

  TLorentzVector doubledimuon_p4;
  TLorentzVector higdim_p4;
  TLorentzVector lowdim_p4;
  TLorentzVector mHighPhi_p4;
  TLorentzVector mLowJPsi_p4;
  TLorentzVector mHighJPsi_p4;
  TLorentzVector mLowPhi_p4;

  TLorentzVector doubledimuon_rf_p4;
  TLorentzVector higdim_rf_p4;
  TLorentzVector lowdim_rf_p4;
  TLorentzVector mHighPhi_rf_p4;
  TLorentzVector mLowJPsi_rf_p4;
  TLorentzVector mHighJPsi_rf_p4;
  TLorentzVector mLowPhi_rf_p4;
  TLorentzVector doubledimuon_not_rf_p4;
  TLorentzVector higdim_not_rf_p4;
  TLorentzVector lowdim_not_rf_p4;

  Int_t    doubledimuon_charge, higdim_triggerMatch, lowdim_triggerMatch, doubledimuon_rf_bindx;
  Double_t doubledimuon_vProb,  doubledimuon_vChi2, doubledimuon_cosAlpha, doubledimuon_ctauPV, doubledimuon_ctauErrPV;
  Double_t doubledimuon_rf_vProb,  doubledimuon_rf_vChi2, doubledimuon_rf_cosAlpha, doubledimuon_rf_ctauPV, doubledimuon_rf_ctauErrPV;

  Double_t gen_doubledimuon_m,doubledimuon_m,doubledimuon_pt,dimuon_m,dimuon_pt;
  Double_t mHighPhi_pt,mLowPhi_pt,mHighJPsi_pt,mLowJPsi_pt,doubledimuon_nDof,doubledimuon_m_rf;

  Double_t doubledimuon_pdgid, doubledimuon_phipdg, doubledimuon_isprompt, doubledimuon_phippdl;

  Double_t highDiM_m, highDiM_pt, lowDiM_m, lowDiM_pt;
  Int_t    higdim_triggerMatch_rf, lowdim_triggerMatch_rf;
  Double_t higdim_vProb_rf, higdim_vChi2_rf, higdim_DCA_rf, higdim_ctauPV_rf, higdim_ctauErrPV_rf, higdim_cosAlpha_rf;
  Double_t lowdim_vProb_rf, lowdim_vChi2_rf, lowdim_DCA_rf, lowdim_ctauPV_rf, lowdim_ctauErrPV_rf, lowdim_cosAlpha_rf;

  //Double_t track_d0, track_d0Err, track_dz, track_dxy;
  Double_t higdim_vProb, higdim_vChi2, higdim_DCA, higdim_ctauPV, higdim_ctauErrPV, higdim_cosAlpha;
  Double_t lowdim_vProb, lowdim_vChi2, lowdim_DCA, lowdim_ctauPV, lowdim_ctauErrPV, lowdim_cosAlpha;
  Double_t  highDiMM_fit, highDiMPx_fit, highDiMPy_fit, highDiMPz_fit;

  Double_t mHighJPsi_isLoose, mHighJPsi_isSoft, mHighJPsi_isMedium, mHighJPsi_isHighPt;
  Double_t mLowJPsi_isLoose, mLowJPsi_isSoft, mLowJPsi_isMedium, mLowJPsi_isHighPt;
  Double_t mHighPhi_isLoose, mHighPhi_isSoft, mHighPhi_isMedium, mHighPhi_isHighPt;
  Double_t mLowPhi_isLoose, mLowPhi_isSoft, mLowPhi_isMedium, mLowPhi_isHighPt;

  Double_t mHighJPsi_isTracker, mHighJPsi_isGlobal, mLowJPsi_isTracker, mLowJPsi_isGlobal;
  Double_t mHighPhi_isTracker, mHighPhi_isGlobal, mLowPhi_isTracker, mLowPhi_isGlobal;

  Double_t isBestCandidate;

  Double_t mHighJPsi_type, mLowJPsi_type, mHighPhi_type, mLowPhi_type;
  Double_t          gen_doubledimuon_pdgId;
  TLorentzVector gen_doubledimuon_p4;
  TLorentzVector gen_higdim_p4;
  TLorentzVector gen_lowdim_p4;
  TLorentzVector gen_mHighPhi_p4;
  TLorentzVector gen_mLowJPsi_p4;
  TLorentzVector gen_mHighJPsi_p4;
  TLorentzVector gen_mLowPhi_p4;

  TLorentzVector gen_b_p4;
  TLorentzVector gen_b4_p4;
  TLorentzVector gen_d1_p4;
  TLorentzVector gen_d2_p4;
  TLorentzVector gen_gd1_p4;
  TLorentzVector gen_gd2_p4;
  TLorentzVector gen_gd3_p4;
  TLorentzVector gen_gd4_p4;
  TLorentzVector gen_gd5_p4;
  TLorentzVector gen_gd6_p4;

  TTree* fourmuon_tree, *fourmuon_tree_rf;
  edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;
};

//
// constants, enums and typedefs
//

UInt_t DoubleDiMuonRootuplerFit::isTriggerMatched(pat::CompositeCandidate *diMuon_cand) {
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
DoubleDiMuonRootuplerFit::DoubleDiMuonRootuplerFit(const edm::ParameterSet& iConfig):
        doubledimuon_cand_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("doubledimuon_cand"))),
        doubledimuon_rf_cand_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("doubledimuon_rf_cand"))),
        thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"))),
        primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
        triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
	      isMC_(iConfig.getParameter<bool>("isMC")),
        OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
        OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
        motherpdgid_(iConfig.getParameter<uint32_t>("Mother_pdg")),
        HLTs_(iConfig.getParameter<std::vector<std::string>>("HLTs")),
        HLTFilters_(iConfig.getParameter<std::vector<std::string>>("filters"))
{
	      edm::Service<TFileService> fs;
        fourmuon_tree = fs->make<TTree>("FourMuonTree","Tree of JPsi and Phi in 4 Muons");

        fourmuon_tree->Branch("run",                &run,                "run/D");
        fourmuon_tree->Branch("event",              &event,              "event/D");
        fourmuon_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/D");
        fourmuon_tree->Branch("trigger",            &trigger,            "trigger/D");

        fourmuon_tree->Branch("doubledimuon_p4",   "TLorentzVector", &doubledimuon_not_rf_p4);
        fourmuon_tree->Branch("lowdim_p4",     "TLorentzVector", &lowdim_not_rf_p4);
        fourmuon_tree->Branch("higdim_p4",     "TLorentzVector", &higdim_not_rf_p4);
        fourmuon_tree->Branch("mLowPhi_p4",   "TLorentzVector", &mLowPhi_p4);
        fourmuon_tree->Branch("mLowJPsi_p4",   "TLorentzVector", &mLowJPsi_p4);
        fourmuon_tree->Branch("mHighJPsi_p4",   "TLorentzVector", &mHighJPsi_p4);
        fourmuon_tree->Branch("mLowPhi_p4",   "TLorentzVector", &mLowPhi_p4);

        fourmuon_tree->Branch("doubledimuon_rf_p4",   "TLorentzVector", &doubledimuon_rf_p4);
        fourmuon_tree->Branch("lowdim_rf_p4",     "TLorentzVector", &lowdim_rf_p4);
        fourmuon_tree->Branch("higdim_rf_p4",     "TLorentzVector", &higdim_rf_p4);
        fourmuon_tree->Branch("mLowPhi_rf_p4",   "TLorentzVector", &mLowPhi_rf_p4);
        fourmuon_tree->Branch("mLowJPsi_rf_p4",   "TLorentzVector", &mLowJPsi_rf_p4);
        fourmuon_tree->Branch("mHighJPsi_rf_p4",   "TLorentzVector", &mHighJPsi_rf_p4);
        fourmuon_tree->Branch("mLowPhi_rf_p4",   "TLorentzVector", &mLowPhi_rf_p4);

        // fourmuon_tree->Branch("doubledimuon_rf_bindx", &doubledimuon_rf_bindx, "doubledimuon_rf_bindx/D");

        fourmuon_tree->Branch("higdim_vProb",        &higdim_vProb,        "higdim_vProb/D");
        fourmuon_tree->Branch("higdim_vNChi2",       &higdim_vChi2,        "higdim_vNChi2/D");
        fourmuon_tree->Branch("higdim_DCA",          &higdim_DCA,          "higdim_DCA/D");
        fourmuon_tree->Branch("higdim_ctauPV",       &higdim_ctauPV,       "higdim_ctauPV/D");
        fourmuon_tree->Branch("higdim_ctauErrPV",    &higdim_ctauErrPV,    "higdim_ctauErrPV/D");
        fourmuon_tree->Branch("higdim_cosAlpha",     &higdim_cosAlpha,     "higdim_cosAlpha/D");
        fourmuon_tree->Branch("higdim_triggerMatch", &higdim_triggerMatch, "higdim_triggerMatch/D");

        fourmuon_tree->Branch("lowdim_vProb",        &lowdim_vProb,        "lowdim_vProb/D");
        fourmuon_tree->Branch("lowdim_vNChi2",       &lowdim_vChi2,        "lowdim_vNChi2/D");
        fourmuon_tree->Branch("lowdim_DCA",          &lowdim_DCA,          "lowdim_DCA/D");
        fourmuon_tree->Branch("lowdim_ctauPV",       &lowdim_ctauPV,       "lowdim_ctauPV/D");
        fourmuon_tree->Branch("lowdim_ctauErrPV",    &lowdim_ctauErrPV,    "lowdim_ctauErrPV/D");
        fourmuon_tree->Branch("lowdim_cosAlpha",     &lowdim_cosAlpha,     "lowdim_cosAlpha/D");
        fourmuon_tree->Branch("lowdim_triggerMatch", &lowdim_triggerMatch, "lowdim_triggerMatch/D");

        fourmuon_tree->Branch("doubledimuon_vProb",      &doubledimuon_vProb,        "doubledimuon_vProb/D");
        fourmuon_tree->Branch("doubledimuon_vChi2",      &doubledimuon_vChi2,        "doubledimuon_vChi2/D");
        fourmuon_tree->Branch("doubledimuon_nDof",      &doubledimuon_nDof,        "doubledimuon_nDof/D");
        fourmuon_tree->Branch("doubledimuon_cosAlpha",   &doubledimuon_cosAlpha,     "doubledimuon_cosAlpha/D");
        fourmuon_tree->Branch("doubledimuon_ctauPV",     &doubledimuon_ctauPV,       "doubledimuon_ctauPV/D");
        fourmuon_tree->Branch("doubledimuon_ctauErrPV",  &doubledimuon_ctauErrPV,    "doubledimuon_ctauErrPV/D");
        fourmuon_tree->Branch("doubledimuon_charge",     &doubledimuon_charge,       "doubledimuon_charge/D");

        fourmuon_tree->Branch("dimuonditrk_cosAlphaDZ",      &dimuonditrk_cosAlphaDZ,        "dimuonditrk_cosAlphaDZ/D");
        fourmuon_tree->Branch("dimuonditrk_ctauPVDZ",      &dimuonditrk_ctauPVDZ,        "dimuonditrk_ctauPVDZ/D");
        fourmuon_tree->Branch("dimuonditrk_ctauErrPVDZ",      &dimuonditrk_ctauErrPVDZ,        "dimuonditrk_ctauErrPVDZ/D");

        fourmuon_tree->Branch("dimuonditrk_cosAlphaBS",      &dimuonditrk_cosAlphaBS,        "dimuonditrk_cosAlphaBS/D");
        fourmuon_tree->Branch("dimuonditrk_ctauPVBS",      &dimuonditrk_ctauPVBS,        "dimuonditrk_ctauPVBS/D");
        fourmuon_tree->Branch("dimuonditrk_ctauErrPVBS",      &dimuonditrk_ctauErrPVBS,        "dimuonditrk_ctauErrPVBS/D");

        fourmuon_tree->Branch("dimuonditrk_vx",          &dimuonditrk_vx,          "dimuonditrk_vx/D");
        fourmuon_tree->Branch("dimuonditrk_vy",          &dimuonditrk_vy,          "dimuonditrk_vy/D");
        fourmuon_tree->Branch("dimuonditrk_vz",          &dimuonditrk_vz,          "dimuonditrk_vz/D");

        fourmuon_tree->Branch("dimuonditrk_dca_m1m2",      &dimuonditrk_vProb,        "dimuonditrk_dca_m1m2/D");
        fourmuon_tree->Branch("dimuonditrk_dca_m1t1",      &dimuonditrk_vProb,        "dimuonditrk_dca_m1t1/D");
        fourmuon_tree->Branch("dimuonditrk_dca_m1t2",      &dimuonditrk_vProb,        "dimuonditrk_dca_m1t2/D");
        fourmuon_tree->Branch("dimuonditrk_dca_m2t1",      &dimuonditrk_vProb,        "dimuonditrk_dca_m2t1/D");
        fourmuon_tree->Branch("dimuonditrk_dca_m2t2",      &dimuonditrk_vProb,        "dimuonditrk_dca_m2t2/D");
        fourmuon_tree->Branch("dimuonditrk_dca_t1t2",      &dimuonditrk_vProb,        "dimuonditrk_dca_t1t2/D");

        fourmuon_tree->Branch("dimuonditrk_rf_vProb",      &dimuonditrk_rf_vProb,        "dimuonditrk_rf_vProb/D");
        fourmuon_tree->Branch("dimuonditrk_rf_vChi2",      &dimuonditrk_rf_vChi2,        "dimuonditrk_rf_vChi2/D");
        fourmuon_tree->Branch("dimuonditrk_rf_nDof",       &dimuonditrk_rf_nDof,         "dimuonditrk_rf_nDof/D");
        fourmuon_tree->Branch("dimuonditrk_rf_cosAlpha",   &dimuonditrk_rf_cosAlpha,     "dimuonditrk_rf_cosAlpha/D");
        fourmuon_tree->Branch("dimuonditrk_rf_ctauPV",     &dimuonditrk_rf_ctauPV,       "dimuonditrk_rf_ctauPV/D");
        fourmuon_tree->Branch("dimuonditrk_rf_ctauErrPV",  &dimuonditrk_rf_ctauErrPV,    "dimuonditrk_rf_ctauErrPV/D");

        fourmuon_tree->Branch("dimuonditrk_rf_c_vProb",      &dimuonditrk_rf_c_vProb,        "dimuonditrk_rf_c_vProb/D");
        fourmuon_tree->Branch("dimuonditrk_rf_c_vChi2",      &dimuonditrk_rf_c_vChi2,        "dimuonditrk_rf_c_vChi2/D");
        fourmuon_tree->Branch("dimuonditrk_rf_c_nDof",       &dimuonditrk_rf_c_nDof,         "dimuonditrk_rf_c_nDof/D");
        fourmuon_tree->Branch("dimuonditrk_rf_c_cosAlpha",   &dimuonditrk_rf_c_cosAlpha,     "dimuonditrk_rf_c_cosAlpha/D");
        fourmuon_tree->Branch("dimuonditrk_rf_c_ctauPV",     &dimuonditrk_rf_c_ctauPV,       "dimuonditrk_rf_c_ctauPV/D");
        fourmuon_tree->Branch("dimuonditrk_rf_c_ctauErrPV",  &dimuonditrk_rf_c_ctauErrPV,    "dimuonditrk_rf_c_ctauErrPV/D");

        //
        // fourmuon_tree->Branch("highDiMM_fit",  &highDiMM_fit,    "highDiMM_fit/D");
        // fourmuon_tree->Branch("highDiMPx_fit",  &highDiMPx_fit,    "highDiMPx_fit/D");
        // fourmuon_tree->Branch("highDiMPy_fit",  &highDiMPy_fit,    "highDiMPy_fit/D");
        // fourmuon_tree->Branch("highDiMPz_fit",  &highDiMPz_fit,    "highDiMPz_fit/D");

        fourmuon_tree->Branch("gen_doubledimuon_m",        &gen_doubledimuon_m,        "gen_doubledimuon_m/D");
        fourmuon_tree->Branch("doubledimuon_m",       &doubledimuon_m,        "doubledimuon_m/D");
        fourmuon_tree->Branch("doubledimuon_m_rf",       &doubledimuon_m_rf,        "doubledimuon_m_rf/D");
        fourmuon_tree->Branch("doubledimuon_pt",          &doubledimuon_pt,          "doubledimuon_pt/D");
        fourmuon_tree->Branch("highDiM_m",       &highDiM_m,       "highDiMn_m/D");
        fourmuon_tree->Branch("highDiM_pt",    &highDiM_pt,    "highDiM_pt/D");
        fourmuon_tree->Branch("lowDiM_m",     &lowDiM_m,     "lowDiM_m/D");
        fourmuon_tree->Branch("lowDiM_pt",       &lowDiM_pt,        "lowDiM_pt/D");
        fourmuon_tree->Branch("mHighJPsi_pt",          &mHighJPsi_pt,          "mHighJPsi_pt/D");
        fourmuon_tree->Branch("mLowJPsi_pt",       &mLowJPsi_pt,       "mLowJPsi_pt/D");
        fourmuon_tree->Branch("mHighPhi_pt",    &mHighPhi_pt,    "mHighPhi_pt/D");
        fourmuon_tree->Branch("mLowPhi_pt",     &mLowPhi_pt,     "mLowPhi_pt/D");

        fourmuon_tree->Branch("mHighJPsi_eta",        &mHighJPsi_eta,        "mHighJPsi_eta/D");
        fourmuon_tree->Branch("lowKaon_eta",        &lowKaon_eta,        "lowKaon_eta/D");
        fourmuon_tree->Branch("highMuon_eta",        &highMuon_eta,        "highMuon_eta/D");
        fourmuon_tree->Branch("lowMuon_eta",        &lowMuon_eta,        "lowMuon_eta/D");

        fourmuon_tree->Branch("mHighJPsi_phi",        &mHighJPsi_phi,        "mHighJPsi_phi/D");
        fourmuon_tree->Branch("lowKaon_phi",        &lowKaon_phi,        "lowKaon_phi/D");
        fourmuon_tree->Branch("highMuon_phi",        &highMuon_phi,        "highMuon_phi/D");
        fourmuon_tree->Branch("lowMuon_phi",        &lowMuon_phi,        "lowMuon_phi/D");

        fourmuon_tree->Branch("mHighJPsiMatch",     &mHighJPsiMatch,       "mHighJPsiMatch/D");
        fourmuon_tree->Branch("mHighJPsiMatch",     &mHighJPsiMatch,       "mHighJPsiMatch/D");
        fourmuon_tree->Branch("mHighPhiMatch",     &mHighPhiMatch,       "mHighPhiMatch/D");
        fourmuon_tree->Branch("mLowPhiMatch",     &mLowPhiMatch,       "mLowPhiMatch/D");


        //Muon flags
        fourmuon_tree->Branch("mHighJPsi_isLoose",       &mHighJPsi_isLoose,        "mHighJPsi_isLoose/D");
        fourmuon_tree->Branch("mHighJPsi_isSoft",        &mHighJPsi_isSoft,        "mHighJPsi_isSoft/D");
        fourmuon_tree->Branch("mHighJPsi_isMedium",      &mHighJPsi_isMedium,        "mHighJPsi_isMedium/D");
        fourmuon_tree->Branch("mHighJPsi_isHighPt",      &mHighJPsi_isHighPt,        "mHighJPsi_isHighPt/D");

        fourmuon_tree->Branch("mHighJPsi_isTracker",        &mHighJPsi_isTracker,        "mHighJPsi_isTracker/D");
        fourmuon_tree->Branch("mHighJPsi_isGlobal",        &mHighJPsi_isGlobal,        "mHighJPsi_isGlobal/D");

        fourmuon_tree->Branch("mLowJPsi_isLoose",        &mLowJPsi_isLoose,        "mLowJPsi_isLoose/D");
        fourmuon_tree->Branch("mLowJPsi_isSoft",        &mLowJPsi_isSoft,        "mLowJPsi_isSoft/D");
        fourmuon_tree->Branch("mLowJPsi_isMedium",        &mLowJPsi_isMedium,        "mLowJPsi_isMedium/D");
        fourmuon_tree->Branch("mLowJPsi_isHighPt",        &mLowJPsi_isHighPt,        "mLowJPsi_isHighPt/D");

        fourmuon_tree->Branch("mLowJPsi_isTracker",        &mLowJPsi_isTracker,        "mLowJPsi_isTracker/D");
        fourmuon_tree->Branch("mLowJPsi_isGlobal",        &mLowJPsi_isGlobal,        "mLowJPsi_isGlobal/D");

        fourmuon_tree->Branch("mHighJPsi_type",     &mHighJPsi_type,       "mHighJPsi_type/D");
        fourmuon_tree->Branch("mLowJPsi_type",     &mLowJPsi_type,       "mLowJPsi_type/D");

        fourmuon_tree->Branch("mHighPhi_isLoose",        &mHighPhi_isLoose,        "mHighPhi_isLoose/D");
        fourmuon_tree->Branch("mHighPhi_isSoft",        &mHighPhi_isSoft,        "mHighPhi_isSoft/D");
        fourmuon_tree->Branch("mHighPhi_isMedium",        &mHighPhi_isMedium,        "mHighPhi_isMedium/D");
        fourmuon_tree->Branch("mHighPhi_isHighPt",        &mHighPhi_isHighPt,        "mHighPhi_isHighPt/D");

        fourmuon_tree->Branch("mHighPhi_isTracker",        &mHighPhi_isTracker,        "mHighPhi_isTracker/D");
        fourmuon_tree->Branch("mHighPhi_isGlobal",        &mHighPhi_isGlobal,        "mHighPhi_isGlobal/D");

        fourmuon_tree->Branch("mLowPhi_isLoose",        &mLowPhi_isLoose,        "mLowPhi_isLoose/D");
        fourmuon_tree->Branch("mLowPhi_isSoft",        &mLowPhi_isSoft,        "mLowPhi_isSoft/D");
        fourmuon_tree->Branch("mLowPhi_isMedium",        &mLowPhi_isMedium,        "mLowPhi_isMedium/D");
        fourmuon_tree->Branch("mLowPhi_isHighPt",        &mLowPhi_isHighPt,        "mLowPhi_isHighPt/D");

        fourmuon_tree->Branch("mLowPhi_isTracker",        &mLowPhi_isTracker,        "mLowPhi_isTracker/D");
        fourmuon_tree->Branch("mLowPhi_isGlobal",        &mLowPhi_isGlobal,        "mLowPhi_isGlobal/D");

        fourmuon_tree->Branch("mHighPhi_type",     &mHighPhi_type,       "mHighPhi_type/D");
        fourmuon_tree->Branch("mLowPhi_type",     &mLowPhi_type,       "mLowPhi_type/D");

        fourmuon_tree->Branch("mLowPhi_NPixelHits",        &mLowPhi_NPixelHits,        "mLowPhi_NPixelHits/D");
        fourmuon_tree->Branch("mLowPhi_NStripHits",        &mLowPhi_NStripHits,        "mLowPhi_NStripHits/D");
        fourmuon_tree->Branch("mLowPhi_NTrackhits",        &mLowPhi_NTrackhits,        "mLowPhi_NTrackhits/D");
        fourmuon_tree->Branch("mLowPhi_NBPixHits",        &mLowPhi_NBPixHits,        "mLowPhi_NBPixHits/D");

        fourmuon_tree->Branch("mLowPhi_NPixLayers",        &mLowPhi_NPixLayers,        "mLowPhi_NPixLayers/D");
        fourmuon_tree->Branch("mLowPhi_NTraLayers",        &mLowPhi_NTraLayers,        "mLowPhi_NTraLayers/D");
        fourmuon_tree->Branch("mLowPhi_NStrLayers",        &mLowPhi_NStrLayers,        "mLowPhi_NStrLayers/D");
        fourmuon_tree->Branch("mLowPhi_NBPixLayers",        &mLowPhi_NBPixLayers,        "mLowPhi_NBPixLayers/D");

        fourmuon_tree->Branch("mHighPhi_NPixelHits",        &mHighPhi_NPixelHits,        "mHighPhi_NPixelHits/D");
        fourmuon_tree->Branch("mHighPhi_NStripHits",        &mHighPhi_NStripHits,        "mHighPhi_NStripHits/D");
        fourmuon_tree->Branch("mHighPhi_NTrackhits",        &mHighPhi_NTrackhits,        "mHighPhi_NTrackhits/D");
        fourmuon_tree->Branch("mHighPhi_NBPixHits",        &mHighPhi_NBPixHits,        "mHighPhi_NBPixHits/D");

        fourmuon_tree->Branch("mHighPhi_NPixLayers",        &mHighPhi_NPixLayers,        "mHighPhi_NPixLayers/D");
        fourmuon_tree->Branch("mHighPhi_NTraLayers",        &mHighPhi_NTraLayers,        "mHighPhi_NTraLayers/D");
        fourmuon_tree->Branch("mHighPhi_NStrLayers",        &mHighPhi_NStrLayers,        "mHighPhi_NStrLayers/D");
        fourmuon_tree->Branch("mHighPhi_NBPixLayers",        &mHighPhi_NBPixLayers,        "mHighPhi_NBPixLayers/D");

        fourmuon_tree->Branch("mLowJPsi_NPixelHits",        &mLowJPsi_NPixelHits,        "mLowJPsi_NPixelHits/D");
        fourmuon_tree->Branch("mLowJPsi_NStripHits",        &mLowJPsi_NStripHits,        "mLowJPsi_NStripHits/D");
        fourmuon_tree->Branch("mLowJPsi_NTrackhits",        &mLowJPsi_NTrackhits,        "mLowJPsi_NTrackhits/D");
        fourmuon_tree->Branch("mLowJPsi_NBPixHits",        &mLowJPsi_NBPixHits,        "mLowJPsi_NBPixHits/D");

        fourmuon_tree->Branch("mLowJPsi_NPixLayers",        &mLowJPsi_NPixLayers,        "mLowJPsi_NPixLayers/D");
        fourmuon_tree->Branch("mLowJPsi_NTraLayers",        &mLowJPsi_NTraLayers,        "mLowJPsi_NTraLayers/D");
        fourmuon_tree->Branch("mLowJPsi_NStrLayers",        &mLowJPsi_NStrLayers,        "mLowJPsi_NStrLayers/D");
        fourmuon_tree->Branch("mLowJPsi_NBPixLayers",        &mLowJPsi_NBPixLayers,        "mLowJPsi_NBPixLayers/D");

        fourmuon_tree->Branch("mHighJPsi_NPixelHits",        &mHighJPsi_NPixelHits,        "mHighJPsi_NPixelHits/D");
        fourmuon_tree->Branch("mHighJPsi_NStripHits",        &mHighJPsi_NStripHits,        "mHighJPsi_NStripHits/D");
        fourmuon_tree->Branch("mHighJPsi_NTrackhits",        &mHighJPsi_NTrackhits,        "mHighJPsi_NTrackhits/D");
        fourmuon_tree->Branch("mHighJPsi_NBPixHits",        &mHighJPsi_NBPixHits,        "mHighJPsi_NBPixHits/D");

        fourmuon_tree->Branch("mHighJPsi_NPixLayers",        &mHighJPsi_NPixLayers,        "mHighJPsi_NPixLayers/D");
        fourmuon_tree->Branch("mHighJPsi_NTraLayers",        &mHighJPsi_NTraLayers,        "mHighJPsi_NTraLayers/D");
        fourmuon_tree->Branch("mHighJPsi_NStrLayers",        &mHighJPsi_NStrLayers,        "mHighJPsi_NStrLayers/D");
        fourmuon_tree->Branch("mHighJPsi_NBPixLayers",        &mHighJPsi_NBPixLayers,        "mHighJPsi_NBPixLayers/D");


        fourmuon_tree->Branch("mHighPhi_dpt",    &mHighPhi_dpt,    "mHighPhi_dpt/D");
        fourmuon_tree->Branch("mLowPhi_dpt",     &mLowPhi_dpt,     "mLowPhi_dpt/D");
        fourmuon_tree->Branch("mHighPhi_dr",    &mHighPhi_dr,    "mHighPhi_dr/D");
        fourmuon_tree->Branch("mLowPhi_dr",     &mLowPhi_dr,     "mLowPhi_dr/D");

        fourmuon_tree->Branch("mHighJPsi_dpt",    &mHighJPsi_dpt,    "mHighJPsi_dpt/D");
        fourmuon_tree->Branch("mLowJPsi_dpt",     &mLowJPsi_dpt,     "mLowJPsi_dpt/D");
        fourmuon_tree->Branch("mHighJPsi_dr",    &mHighJPsi_dr,    "mHighJPsi_dr/D");
        fourmuon_tree->Branch("mLowJPsi_dr",     &mLowJPsi_dr,     "mLowJPsi_dr/D");

        fourmuon_tree->Branch("isBestCandidate",        &isBestCandidate,        "isBestCandidate/D");

	if(isMC_)
	  {
            fourmuon_tree->Branch("gen_doubledimuon_pdgId", &gen_doubledimuon_pdgId, "gen_doubledimuon_pdgId/D");
      	    fourmuon_tree->Branch("gen_doubledimuon_p4",    "TLorentzVector", &gen_doubledimuon_p4);
      	    fourmuon_tree->Branch("gen_higdim_p4",      "TLorentzVector", &gen_higdim_p4);
      	    fourmuon_tree->Branch("gen_lowdim_p4",      "TLorentzVector", &gen_lowdim_p4);
            fourmuon_tree->Branch("gen_mHighPhi_p4",    "TLorentzVector", &gen_mHighPhi_p4);
            fourmuon_tree->Branch("gen_mLowJPsi_p4",    "TLorentzVector", &gen_mLowJPsi_p4);
            fourmuon_tree->Branch("gen_mHighJPsi_p4",    "TLorentzVector", &gen_mHighJPsi_p4);
            fourmuon_tree->Branch("gen_mLowPhi_p4",    "TLorentzVector", &gen_mLowPhi_p4);

            fourmuon_tree->Branch("gen_doubledimuon_m",  &gen_doubledimuon_m,    "gen_doubledimuon_m/D");

            fourmuon_tree->Branch("gen_b4_p4", "TLorentzVector",  &gen_b4_p4);
            fourmuon_tree->Branch("gen_d1_p4",  "TLorentzVector",  &gen_d1_p4);
            fourmuon_tree->Branch("gen_d2_p4",  "TLorentzVector",  &gen_d2_p4);
            fourmuon_tree->Branch("gen_gd1_p4",  "TLorentzVector",  &gen_gd1_p4);
            fourmuon_tree->Branch("gen_gd2_p4",  "TLorentzVector",  &gen_gd2_p4);
            fourmuon_tree->Branch("gen_gd3_p4", "TLorentzVector",  &gen_gd3_p4);
            fourmuon_tree->Branch("gen_gd4_p4",  "TLorentzVector",  &gen_gd4_p4);
            fourmuon_tree->Branch("gen_gd5_p4",  "TLorentzVector",  &gen_gd5_p4);
            fourmuon_tree->Branch("gen_gd6_p4",  "TLorentzVector",  &gen_gd6_p4);

            fourmuon_tree->Branch("doubledimuon_pdgid",  &doubledimuon_pdgid,    "doubledimuon_pdgid/D");
            fourmuon_tree->Branch("doubledimuon_phipdg",  &doubledimuon_phipdg,    "doubledimuon_phipdg/D");
            fourmuon_tree->Branch("doubledimuon_isprompt",  &doubledimuon_isprompt,    "doubledimuon_isprompt/D");
            fourmuon_tree->Branch("doubledimuon_phippdl",  &doubledimuon_phippdl,    "doubledimuon_phippdl/D");

	  }
    genCands_ = consumes< std::vector <reco::GenParticle> >((edm::InputTag)"prunedGenParticles");
    packCands_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");
}

DoubleDiMuonRootuplerFit::~DoubleDiMuonRootuplerFit() {}

//
// member functions
//

bool DoubleDiMuonRootuplerFit::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   for (size_t i=0; i< particle->numberOfMothers(); i++) {
      if (isAncestor(ancestor, particle->mother(i))) return true;
   }
   return false;
}


// ------------ method called for each event  ------------
void DoubleDiMuonRootuplerFit::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
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
  gen_doubledimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_higdim_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_lowdim_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_mHighJPsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_mLowJPsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_mHighPhi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_mLowPhi_p4.SetPtEtaPhiM(0.,0.,0.,0.);

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

   gen_doubledimuon_pdgId = 0;

   edm::Handle< std::vector <reco::GenParticle> > pruned;
   iEvent.getByToken(genCands_, pruned);

   // Packed particles are all the status 1. The navigation to pruned is possible (the other direction should be made by hand)
   edm::Handle<pat::PackedGenParticleCollection> packed;
   iEvent.getByToken(packCands_,  packed);

   if ( (isMC_ || OnlyGen_) && packed.isValid() && pruned.isValid() ) {
     for (size_t i=0; i<pruned->size(); i++) {
       // std::cout << "Valid"<<std::endl;
       const reco::Candidate *afourmuon = &(*pruned)[i];

       if ( (abs(afourmuon->pdgId()) == motherpdgid_) && (afourmuon->status() == 2))
         gen_b_p4.SetPtEtaPhiM(afourmuon->pt(),afourmuon->eta(),afourmuon->phi(),afourmuon->mass());

       if ( (abs(afourmuon->pdgId()) == motherpdgid_) && (afourmuon->status() == 2) && (afourmuon->numberOfDaughters() > 1) && (afourmuon->numberOfDaughters() < 7) ) {
         //asking for decay (status==2) && two daughters
         bool goToJPsi = false;
         bool goToPhi = false;

         bool muJPsiP = false, muJPsiN = false, muPhiP = false, muPhiN = false;

         int noDaughters = 0, noGDaughters = 0;
         int theJPsi = 0, thePhi = 0, theJPsiMuP = 0, theJPsiMuN = 0, thePhiMuP = 0, thePhiMuN = 0;

         std::vector<const reco::Candidate *> daughters,gdaughters;
         for(size_t j = 0; j < afourmuon->numberOfDaughters(); ++j)
         {
           const reco::Candidate * daughter = afourmuon->daughter(j);
           if(daughter->mother(daughter->numberOfMothers()-1) != afourmuon) continue;
           if(daughter->pdgId() == 443)
           {
             goToJPsi=true;
             theJPsi = j;
           }
           if(daughter->pdgId() == 333)
           {
             thePhi = j;
             goToPhi=true;
           }

           daughters.push_back(daughter);

           ++noDaughters;

         }

         for (size_t j = 0; j < daughters.size(); j++) {

           if(daughters[j]->status() != 2) continue;

           for(size_t k = 0; k <daughters[j]->numberOfDaughters(); ++k)
           {
             const reco::Candidate * gdaughter = daughters[j]->daughter(k);
             if(gdaughter->mother(gdaughter->numberOfMothers()-1) != daughters[j]) continue;
             gdaughters.push_back(gdaughter);

             if(goToPhi && goToJPsi)
             {
               if(gdaughter->pdgId()==-13)
               {
                 theJPsiMuP = j;
                 muJPsiP=true;
               }
               if(gdaughter->pdgId()==13)
               {
                 theJPsiMuN = j;
                 muJPsiN=true;
               }
               if(gdaughter->pdgId()==13)
               {
                 thePhiMuP = j;
                 muPhiP=true;
               }
               if(gdaughter->pdgId()==-13)
               {
                 thePhiMuN = j;
                 muPhiN=true;
               }
             }

             ++noGDaughters;
           }
         }

         if(noDaughters == 2 && noGDaughters > 3 && noGDaughters < 7 && goToJPsi)
         {

           // for (size_t j = 0; j < daughters.size(); j++)
           //   std::cout << "Daughter no. " << j << " - id : " << daughters[j]->pdgId() << std::endl;
           //
           // for (size_t j = 0; j < gdaughters.size(); j++)
           //   std::cout << "GrandDaughter no. " << j << " - id : " << gdaughters[j]->pdgId() << std::endl;

           gen_b4_p4.SetPtEtaPhiM(afourmuon->pt(),afourmuon->eta(),afourmuon->phi(),afourmuon->mass());
           gen_d1_p4.SetPtEtaPhiM(daughters[0]->pt(),daughters[0]->eta(),daughters[0]->phi(),daughters[0]->mass());
           gen_d2_p4.SetPtEtaPhiM(daughters[1]->pt(),daughters[1]->eta(),daughters[1]->phi(),daughters[1]->mass());

           gen_gd1_p4.SetPtEtaPhiM(gdaughters[0]->pt(),gdaughters[0]->eta(),gdaughters[0]->phi(),gdaughters[0]->mass());
           gen_gd2_p4.SetPtEtaPhiM(gdaughters[1]->pt(),gdaughters[1]->eta(),gdaughters[1]->phi(),gdaughters[1]->mass());
           gen_gd3_p4.SetPtEtaPhiM(gdaughters[2]->pt(),gdaughters[2]->eta(),gdaughters[2]->phi(),gdaughters[2]->mass());
           gen_gd4_p4.SetPtEtaPhiM(gdaughters[3]->pt(),gdaughters[3]->eta(),gdaughters[3]->phi(),gdaughters[3]->mass());

           if(noGDaughters > 4)
             gen_gd5_p4.SetPtEtaPhiM(gdaughters[4]->pt(),gdaughters[4]->eta(),gdaughters[4]->phi(),gdaughters[4]->mass());
           if(noGDaughters > 5)
             gen_gd6_p4.SetPtEtaPhiM(gdaughters[5]->pt(),gdaughters[5]->eta(),gdaughters[5]->phi(),gdaughters[5]->mass());
         }

         if(muJPsiP && muJPsiN && muPhiP && muPhiN)
         {

           gen_doubledimuon_p4.SetPtEtaPhiM(afourmuon->pt(),afourmuon->eta(),afourmuon->phi(),afourmuon->mass());
           gen_higdim_p4.SetPtEtaPhiM(daughters[theJPsi]->pt(),daughters[theJPsi]->eta(),daughters[theJPsi]->phi(),daughters[theJPsi]->mass());
           gen_lowdim_p4.SetPtEtaPhiM(daughters[thePhi]->pt(),daughters[thePhi]->eta(),daughters[thePhi]->phi(),daughters[thePhi]->mass());
           gen_mLowJPsi_p4.SetPtEtaPhiM(gdaughters[theJPsiMuN]->pt(),gdaughters[theJPsiMuN]->eta(),gdaughters[theJPsiMuN]->phi(),gdaughters[theJPsiMuN]->mass());
           gen_mHighJPsi_p4.SetPtEtaPhiM(gdaughters[theJPsiMuP]->pt(),gdaughters[theJPsiMuP]->eta(),gdaughters[theJPsiMuP]->phi(),gdaughters[theJPsiMuP]->mass());
           gen_mLowPhi_p4.SetPtEtaPhiM(gdaughters[thePhiMuN]->pt(),gdaughters[thePhiMuN]->eta(),gdaughters[thePhiMuN]->phi(),gdaughters[thePhiMuN]->mass());
           gen_mHighPhi_p4.SetPtEtaPhiM(gdaughters[thePhiMuP]->pt(),gdaughters[thePhiMuP]->eta(),gdaughters[thePhiMuP]->phi(),gdaughters[thePhiMuP]->mass());
           gen_doubledimuon_pdgId = afourmuon->pdgId();
       }
     } // for (size
   }
   }

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
      doubledimuon_nDof      = doubledimuon_rf_cand.userFloat("nDof");
      doubledimuon_cosAlpha  = doubledimuon_rf_cand.userFloat("cosAlpha");
      doubledimuon_ctauPV    = doubledimuon_rf_cand.userFloat("ctauPV");
      doubledimuon_ctauErrPV = doubledimuon_rf_cand.userFloat("ctauErrPV");
      doubledimuon_m_rf      = doubledimuon_rf_cand.mass();

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

      doubledimuon_m   = doubledimuon_cand.mass();
      doubledimuon_m_rf= doubledimuon_rf_cand.mass();
      doubledimuon_pt  = doubledimuon_cand.pt();

      doubledimuon_pdgid    = doubledimuon_rf_cand.userInt("phiGenPdgId");
      doubledimuon_phipdg   = doubledimuon_rf_cand.userFloat("phiPpdlTrue");
      doubledimuon_isprompt = doubledimuon_rf_cand.userInt("xGenPdgId");
      doubledimuon_phippdl  = doubledimuon_rf_cand.userFloat("xGenIsPrompt");

      reco::Candidate::LorentzVector vJpsiP = higdim_cand_rf->daughter("muon1")->p4();
      reco::Candidate::LorentzVector vJpsiM = higdim_cand_rf->daughter("muon2")->p4();

      if (higdim_cand_rf->daughter("muon1")->charge() < 0) {
         vJpsiP = higdim_cand_rf->daughter("muon2")->p4();
         vJpsiM = higdim_cand_rf->daughter("muon1")->p4();
      }

      mLowPhi_rf_p4.SetPtEtaPhiM(vJpsiP.pt(), vJpsiP.eta(), vJpsiP.phi(), vJpsiP.mass());
      mLowJPsi_rf_p4.SetPtEtaPhiM(vJpsiM.pt(), vJpsiM.eta(), vJpsiM.phi(), vJpsiM.mass());

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
      higdim_triggerMatch = DoubleDiMuonRootuplerFit::isTriggerMatched(higdim_cand);

      lowdim_vProb        = lowdim_cand->userFloat("vProb");
      lowdim_vChi2        = lowdim_cand->userFloat("vNChi2");
      lowdim_DCA          = lowdim_cand->userFloat("DCA");
      lowdim_ctauPV       = lowdim_cand->userFloat("ppdlPV");
      lowdim_ctauErrPV    = lowdim_cand->userFloat("ppdlErrPV");
      lowdim_cosAlpha     = lowdim_cand->userFloat("cosAlpha");
      lowdim_triggerMatch = DoubleDiMuonRootuplerFit::isTriggerMatched(lowdim_cand);

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

      mHighJPsi_p4.SetPtEtaPhiM(vJpsiP.pt(), vJpsiP.eta(), vJpsiP.phi(), vJpsiP.mass());
      mLowJPsi_p4.SetPtEtaPhiM(vJpsiM.pt(), vJpsiM.eta(), vJpsiM.phi(), vJpsiM.mass());

      mHighJPsi_isLoose   = (Double_t) highPatMuonP->isLooseMuon();
      mHighJPsi_isSoft    = (Double_t) highPatMuonP->isSoftMuon(thePrimaryV);
      mHighJPsi_isMedium  = (Double_t) highPatMuonP->isMediumMuon();
      mHighJPsi_isHighPt  = (Double_t) highPatMuonP->isHighPtMuon(thePrimaryV);
      mHighJPsi_isTracker = (Double_t) highPatMuonP->isTrackerMuon();
      mHighJPsi_isGlobal  = (Double_t) highPatMuonP->isGlobalMuon();
      mLowJPsi_isLoose    = (Double_t) highPatMuonN->isLooseMuon();
      mLowJPsi_isSoft     = (Double_t) highPatMuonN->isSoftMuon(thePrimaryV);
      mLowJPsi_isMedium   = (Double_t) highPatMuonN->isMediumMuon();
      mLowJPsi_isHighPt   = (Double_t) highPatMuonN->isHighPtMuon(thePrimaryV);
      mLowJPsi_isTracker  = (Double_t) highPatMuonN->isTrackerMuon();
      mLowJPsi_isGlobal   = (Double_t) highPatMuonN->isGlobalMuon();
      mHighJPsi_type      = (Double_t) highPatMuonP->type();
      mLowJPsi_type       = (Double_t) highPatMuonN->type();

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

      highDiM_m        = (Double_t) higdim_cand->mass();
      highDiM_pt       = (Double_t) higdim_cand->pt();
      lowDiM_m         = (Double_t) lowdim_cand->mass();
      lowDiM_pt        = (Double_t) lowdim_cand->pt();

      mHighJPsi_pt     = (Double_t) std::max(vJpsiP.pt(),vJpsiM.pt());
      mLowJPsi_pt      = (Double_t) -std::max(-vJpsiP.pt(),-vJpsiM.pt());
      mHighPhi_pt      = (Double_t) std::max(vPhiP.pt(),vPhiP.pt());
      mLowPhi_pt       = (Double_t) -std::max(-vPhiP.pt(),-vPhiP.pt());

      mHighPhi_isLoose    = (Double_t)  lowPatMuonP->isLooseMuon();
      mHighPhi_isSoft     = (Double_t)  lowPatMuonP->isSoftMuon(thePrimaryV);
      mHighPhi_isMedium   = (Double_t) lowPatMuonP->isMediumMuon();
      mHighPhi_isHighPt   = (Double_t) lowPatMuonP->isHighPtMuon(thePrimaryV);
      mHighPhi_isTracker  = (Double_t) lowPatMuonP->isTrackerMuon();
      mHighPhi_isGlobal   = (Double_t) lowPatMuonP->isGlobalMuon();
      mLowPhi_isLoose     = (Double_t) lowPatMuonN->isLooseMuon();
      mLowPhi_isSoft      = (Double_t) lowPatMuonN->isSoftMuon(thePrimaryV);
      mLowPhi_isMedium    = (Double_t) lowPatMuonN->isMediumMuon();
      mLowPhi_isHighPt    = (Double_t) lowPatMuonN->isHighPtMuon(thePrimaryV);
      mLowPhi_isTracker   = (Double_t) lowPatMuonN->isTrackerMuon();
      mLowPhi_isGlobal    = (Double_t) lowPatMuonN->isGlobalMuon();
      mHighPhi_type       = (Double_t)  lowPatMuonP->type();
      mLowPhi_type        = (Double_t) lowPatMuonN->type();

      mHighPhi_p4.SetPtEtaPhiM(vPhiP.pt(), vPhiP.eta(), vPhiP.phi(), vPhiP.mass());
      mLowPhi_p4.SetPtEtaPhiM(vPhiM.pt(), vPhiM.eta(), vPhiM.phi(), vPhiM.mass());

      fourmuon_tree->Fill();

      if (OnlyBest_) break;
      else if(i==0)
      isBestCandidate = false;

      }

  }

}

// ------------ method called once each job just before starting event loop  ------------
void DoubleDiMuonRootuplerFit::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void DoubleDiMuonRootuplerFit::endJob() {}

// ------------ method called when starting to processes a run  ------------
void DoubleDiMuonRootuplerFit::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void DoubleDiMuonRootuplerFit::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void DoubleDiMuonRootuplerFit::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void DoubleDiMuonRootuplerFit::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DoubleDiMuonRootuplerFit::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DoubleDiMuonRootuplerFit);
