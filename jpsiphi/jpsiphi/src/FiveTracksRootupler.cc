/*
   Package:    FiveTracksRootupler
   Class:      FiveTracksRootupler

   Description: make rootuple of DiMuon-DiTrack combination

   Original Author: Adriano Di Florio
   Created:  based on Alberto Sanchez Hernandez PsiTrkTrk Code

*/

// Gen Particles
// 0 : an empty entry with no meaningful information and therefore to be skipped unconditionally
// 1 : a final-state particle, i.e. a particle that is not decayed further by the generator
// 2 : decayed Standard Model hadron or tau or mu lepton, excepting virtual intermediate states thereof (i.e. the particle must undergo a normal decay, not e.g. a shower branching).
// 3 : a documentation entry
// 4 : an incoming beam particle
// 5-10 : undefined, reserved for future standards
// 11-200: an intermediate (decayed/branched/...) particle that does not fulfill the criteria of status code 2, with a generator-dependent classification of its nature.
// 201- : at the disposal of the user, in particular for event tracking in the detector

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

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
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

class FiveTracksRootupler : public edm::EDAnalyzer {
   public:
      explicit FiveTracksRootupler(const edm::ParameterSet&);
      ~FiveTracksRootupler() override;

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
  edm::EDGetTokenT<pat::CompositeCandidateCollection> fivetracks_cand_Label;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuonditrk_rf_cand_Label;
  edm::EDGetTokenT<reco::BeamSpot> thebeamspot_;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  int  dimuonditrk_pdgid_, dimuon_pdgid_, trak_pdgid_, pdgid_;
  bool isMC_,OnlyBest_,OnlyGen_ ;
  UInt_t motherpdgid_,phipdgid_,jpspdgid_;
  std::vector<std::string>  HLTs_;
  std::vector<std::string>  HLTFilters_;
  std::string treeName_;

  std::vector < Double_t > fiveTracksMass, fiveTracksMassRef, triTrakMass;
  std::vector < Double_t > fiveTracksVProb, fiveTracksVNDof, fiveTracksVChi2;
  std::vector < Double_t > fiveTracksCosAlpha, fiveTracksCTau, fiveTracksCTauErr;
  std::vector < Double_t > psiPrimeSame, psiPrimeMixed;
  std::vector < Double_t > psiPrimeSame_ditrak, psiPrimeMixed_ditrak;
  std::vector < Double_t > dimuonDiTrkOne, dimuonDiTrkTwo, dimuonDiTrkThree;
  std::vector < Double_t > ditrakOne, ditrakTwo, ditrakThree;
  std::vector < Double_t > trackOneMass, trackTwoMass, trackThreeMass;
  std::vector < Double_t > psiPrimeSame_p_m, psiPrimeSame_m_m, psiPrimeMixed_p_m, psiPrimeMixed_m_m;

  Double_t testMass;

  std::vector < TLorentzVector > five_p4, five_ref_p4;
  std::vector < TLorentzVector > psiPrimeSame_p4, psiPrimeMixed_p4, psiPrimeSame_ditrak_p4, psiPrimeMixed_ditrak_p4;
  std::vector < TLorentzVector >  dimuonDiTrkOne_p4, dimuonDiTrkTwo_p4, dimuonDiTrkThree_p4;
  std::vector < TLorentzVector > ditrakOne_p4, ditrakTwo_p4, ditrakThree_p4;
  UInt_t run, event, lumi, numPrimaryVertices, trigger, dimuonditrk_id;

  TLorentzVector dimuonditrk_p4, dimuon_p4, ditrak_p4;
  TLorentzVector muonp_p4;
  TLorentzVector muonn_p4;
  TLorentzVector kaonp_p4;
  TLorentzVector kaonn_p4;

  TLorentzVector dimuonditrk_not_rf_p4;
  TLorentzVector dimuon_rf_p4, dimuon_not_rf_p4;
  TLorentzVector ditrakOne_rf_p4, ditrakOne_not_rf_p4;
  TLorentzVector muonp_rf_p4;
  TLorentzVector muonn_rf_p4;
  TLorentzVector kaonp_rf_p4;
  TLorentzVector kaonn_rf_p4;

  Int_t dimuonditrk_charge;

  Double_t highKaonMatch, lowKaonMatch, highMuonMatch, lowMuonMatch;
  Double_t dimuonditrk_vProb,  dimuonditrk_vChi2, dimuonditrk_cosAlpha, dimuonditrk_ctauPV, dimuonditrk_ctauErrPV;

  Double_t dimuon_vProb, dimuon_vChi2, dimuon_DCA, dimuon_ctauPV, dimuon_ctauErrPV, dimuon_cosAlpha;

  Double_t highTrack_pt, highTrack_eta, highTrack_phi, highTrack_charge;
  Double_t lowTrack_pt, lowTrack_eta, lowTrack_phi, lowTrack_charge;
  Double_t thirdTrack_pt, thirdTrack_eta, thirdTrack_phi, thirdTrack_charge;
  Double_t thirdTrack_dz, thirdTrack_dxy;

  Double_t dimuonDiTrkOne_pt, dimuonDiTrkOne_eta, dimuonDiTrkOne_phi, dimuonDiTrkOne_charge;
  Double_t dimuonDiTrkTwo_pt, dimuonDiTrkTwo_eta, dimuonDiTrkTwo_phi, dimuonDiTrkTwo_charge;
  Double_t dimuonDiTrkThree_pt, dimuonDiTrkThree_eta, dimuonDiTrkThree_phi, dimuonDiTrkThree_charge;

  Double_t psiPrimeSame_pt, psiPrimeSame_eta, psiPrimeSame_phi, psiPrimeSame_n;
  Double_t psiPrimeSame_p_pt, psiPrimeSame_p_eta, psiPrimeSame_p_phi, psiPrimeSame_p_n;
  Double_t psiPrimeSame_m_pt, psiPrimeSame_m_eta, psiPrimeSame_m_phi, psiPrimeSame_m_n;

  Double_t psiPrimeMixed_pt, psiPrimeMixed_eta, psiPrimeMixed_phi, psiPrimeMixed_n;
  Double_t psiPrimeMixed_p_pt, psiPrimeMixed_p_eta, psiPrimeMixed_p_phi, psiPrimeMixed_p_n;
  Double_t psiPrimeMixed_m_pt, psiPrimeMixed_m_eta, psiPrimeMixed_m_phi, psiPrimeMixed_m_n;
  Double_t psiPrimeSame_ditrak_pt, psiPrimeSame_ditrak_eta, psiPrimeSame_ditrak_phi, psiPrimeSame_ditrak_n;
  Double_t psiPrimeMixed_ditrak_pt, psiPrimeMixed_ditrak_eta, psiPrimeMixed_ditrak_phi, psiPrimeMixed_ditrak_n;


  Double_t dimuonditrk_m, dimuonditrk_eta, dimuonditrk_pt , dimuonditrk_phi, dimuonditrk_p;
  Double_t dimuon_m, dimuon_pt, dimuon_eta, dimuon_phi, dimuon_p;
  Double_t ditrak_m, ditrakOne_pt, ditrakOne_eta, ditrakOne_phi, ditrakOne_p;
  Double_t ditrakTwo_pt, ditrakTwo_eta, ditrakTwo_phi, ditrakTwo_p;
  Double_t ditrakThree_pt, ditrakThree_eta, ditrakThree_phi, ditrakThree_p;

  Double_t triTrak_pt, triTrak_eta, triTrak_phi, triTrak_charge;

  Bool_t muonP_isLoose, muonP_isSoft, muonP_isMedium, muonP_isHighPt;
  Bool_t muonN_isLoose, muonN_isSoft, muonN_isMedium, muonN_isHighPt;

  Bool_t muonP_isTracker, muonP_isGlobal, muonN_isTracker, muonN_isGlobal;
  UInt_t muonP_type, muonN_type;

  Bool_t muhighTrakose, muonP_rf_isSoft, muonP_rf_isMedium, muonP_rf_isHighPt;
  Bool_t muonN_rf_isLoose, muonN_rf_isSoft, muonN_rf_isMedium, muonN_rf_isHighPt;

  Bool_t muonP_rf_isTracker, muonP_rf_isGlobal, muonN_rf_isTracker, muonN_rf_isGlobal;
  UInt_t muonP_rf_type, muonN_rf_type;

  Double_t track_KP_d0, track_KP_d0Err, track_KP_dz, track_KP_dxy;
  Int_t track_KP_nvsh, track_KP_nvph;

  UInt_t tPMatch, tNMatch;

  Int_t track_KN_nvsh, track_KN_nvph;

  Int_t dimuonditrk_rf_bindx;

  Int_t noFiveCandidates;

  Bool_t isBestCandidate;

  size_t numMasses;

  Int_t          gen_dimuonditrk_pdgId;
  TLorentzVector gen_dimuonditrk_p4;
  TLorentzVector gen_b_p4;
  TLorentzVector gen_dimuon_p4;
  TLorentzVector gen_ditrak_p4;
  TLorentzVector gen_muonp_p4;
  TLorentzVector gen_muonn_p4;
  TLorentzVector gen_kaonp_p4;
  TLorentzVector gen_kaonn_p4;

  TLorentzVector gen_b4_p4;
  TLorentzVector gen_d1_p4;
  TLorentzVector gen_d2_p4;
  TLorentzVector gen_gd1_p4;
  TLorentzVector gen_gd2_p4;
  TLorentzVector gen_gd3_p4;
  TLorentzVector gen_gd4_p4;
  TLorentzVector gen_gd5_p4;
  TLorentzVector gen_gd6_p4;

  TTree* fivetracks_tree;
  edm::EDGetTokenT< std::vector <reco::GenParticle> > genCands_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;
};

UInt_t FiveTracksRootupler::isTriggerMatched(pat::CompositeCandidate *diMuon_cand) {
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
// constants, enums and typedefs
//

//
// static data member definitions
//
static const Double_t psi1SMass =  3.09691;


//
// constructors and destructor
//
FiveTracksRootupler::FiveTracksRootupler(const edm::ParameterSet& iConfig):
        fivetracks_cand_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("FiveTracksCand"))),
        thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"))),
        primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
        triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
	      isMC_(iConfig.getParameter<bool>("isMC")),
        OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
        OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
        motherpdgid_(iConfig.getParameter<uint32_t>("Mother_pdg")),
        phipdgid_(iConfig.getParameter<uint32_t>("JPsi_pdg")),
        jpspdgid_(iConfig.getParameter<uint32_t>("Phi_pdg")),
        HLTs_(iConfig.getParameter<std::vector<std::string>>("HLTs")),
        HLTFilters_(iConfig.getParameter<std::vector<std::string>>("Filters")),
        treeName_(iConfig.getParameter<std::string>("TreeName"))
{
	      edm::Service<TFileService> fs;
        fivetracks_tree = fs->make<TTree>(treeName_.data(),"Tree of DiMuon and DiTrak");

        fivetracks_tree->Branch("run",                &run,                "run/I");
        fivetracks_tree->Branch("event",              &event,              "event/I");
        fivetracks_tree->Branch("lumi",              &lumi,              "lumi/I");
        fivetracks_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");
        fivetracks_tree->Branch("trigger",            &trigger,            "trigger/I");

        fivetracks_tree->Branch("noFiveCandidates",      &noFiveCandidates,      "noFiveCandidates/I");

        fivetracks_tree->Branch("dimuonditrk_id",      &dimuonditrk_id,      "dimuonditrk_id/I");

        fivetracks_tree->Branch("dimuonditrk_p4",     "TLorentzVector", &dimuonditrk_p4);
        fivetracks_tree->Branch("ditrak_p4",     "TLorentzVector", &ditrak_p4);
        fivetracks_tree->Branch("dimuon_p4",     "TLorentzVector", &dimuon_p4);

        fivetracks_tree->Branch("dimuonditrk_m",       &dimuonditrk_m,        "dimuonditrk_m/D");
        fivetracks_tree->Branch("dimuonditrk_pt",      &dimuonditrk_pt,       "dimuonditrk_pt/D");
        fivetracks_tree->Branch("dimuonditrk_eta",     &dimuonditrk_eta,      "dimuonditrk_eta/D");
        fivetracks_tree->Branch("dimuonditrk_phi",     &dimuonditrk_phi,      "dimuonditrk_phi/D");
        fivetracks_tree->Branch("dimuonditrk_p",       &dimuonditrk_p,      "dimuonditrk_p/D");

        fivetracks_tree->Branch("dimuon_m",      &dimuon_m,     "dimuon_m/D");
        fivetracks_tree->Branch("dimuon_pt",     &dimuon_pt,    "dimuon_pt/D");
        fivetracks_tree->Branch("dimuon_eta",    &dimuon_pt,    "dimuon_eta/D");
        fivetracks_tree->Branch("dimuon_phi",    &dimuon_phi,   "dimuon_phi/D");
        fivetracks_tree->Branch("dimuon_p",    &dimuon_p,   "dimuon_p/D");

        fivetracks_tree->Branch("dimuon_phi",    &dimuon_phi,   "dimuon_phi/D");
        fivetracks_tree->Branch("dimuon_p",    &dimuon_p,   "dimuon_p/D");

        fivetracks_tree->Branch("highKaonMatch",    &highKaonMatch,   "highKaonMatch/D");
        fivetracks_tree->Branch("lowKaonMatch",    &lowKaonMatch,   "lowKaonMatch/D");
        fivetracks_tree->Branch("lowMuonMatch",    &lowMuonMatch,   "lowMuonMatch/D");
        fivetracks_tree->Branch("highMuonMatch",    &highMuonMatch,   "highMuonMatch/D");

        //Di Traks
        fivetracks_tree->Branch("ditrak_m",      &ditrak_m,     "ditrak_m/D"); //the original ditrak (supposed to be Phi->KK)

        fivetracks_tree->Branch("ditrakOne_pt",     &ditrakOne_pt,    "ditrakOne_pt/D");
        fivetracks_tree->Branch("ditrakOne_eta",    &ditrakOne_eta,    "ditrakOne_eta/D");
        fivetracks_tree->Branch("ditrakOne_phi",    &ditrakOne_phi,   "ditrakOne_phi/D");
        fivetracks_tree->Branch("ditrakOne_p",    &ditrakOne_p,   "ditrakOne_p/D");

        fivetracks_tree->Branch("ditrakTwo_pt",     &ditrakTwo_pt,    "ditrakTwo_pt/D");
        fivetracks_tree->Branch("ditrakTwo_eta",    &ditrakTwo_eta,    "ditrakTwo_eta/D");
        fivetracks_tree->Branch("ditrakTwo_phi",    &ditrakTwo_phi,   "ditrakTwo_phi/D");
        fivetracks_tree->Branch("ditrakTwo_p",    &ditrakTwo_p,   "ditrakTwo_p/D");

        fivetracks_tree->Branch("ditrakThree_pt",     &ditrakThree_pt,    "ditrakThree_pt/D");
        fivetracks_tree->Branch("ditrakThree_eta",    &ditrakThree_eta,    "ditrakThree_eta/D");
        fivetracks_tree->Branch("ditrakThree_phi",    &ditrakThree_phi,   "ditrakThree_phi/D");
        fivetracks_tree->Branch("ditrakThree_p",    &ditrakThree_p,   "ditrakThree_p/D");


        //The kinematic doesn't change, only mass
        fivetracks_tree->Branch("highTrack_pt",          &highTrack_pt,          "highTrack_pt/D");
        fivetracks_tree->Branch("highTrack_eta",        &highTrack_eta,        "highTrack_eta/D");
        fivetracks_tree->Branch("highTrack_phi",        &highTrack_phi,        "highTrack_phi/D");
        fivetracks_tree->Branch("highTrack_charge",        &highTrack_charge,        "highTrack_charge/D");

        fivetracks_tree->Branch("lowTrack_pt",          &lowTrack_pt,          "lowTrack_pt/D");
        fivetracks_tree->Branch("lowTrack_eta",        &lowTrack_eta,        "lowTrack_eta/D");
        fivetracks_tree->Branch("lowTrack_phi",        &lowTrack_phi,        "lowTrack_phi/D");
        fivetracks_tree->Branch("lowTrack_charge",        &lowTrack_charge,        "lowTrack_charge/D");

        fivetracks_tree->Branch("thirdTrack_pt",          &thirdTrack_pt,          "thirdTrack_pt/D");
        fivetracks_tree->Branch("thirdTrack_eta",        &thirdTrack_eta,        "thirdTrack_eta/D");
        fivetracks_tree->Branch("thirdTrack_phi",        &thirdTrack_phi,        "thirdTrack_phi/D");
        fivetracks_tree->Branch("thirdTrack_charge",        &thirdTrack_charge,        "thirdTrack_charge/D");
        fivetracks_tree->Branch("thirdTrack_dz",        &thirdTrack_dz,        "thirdTrack_dz/D");
        fivetracks_tree->Branch("thirdTrack_dxy",        &thirdTrack_dxy,        "thirdTrack_dxy/D");

        //J/Psi TrTr system
        fivetracks_tree->Branch("dimuonDiTrkOne_pt",          &dimuonDiTrkOne_pt,          "dimuonDiTrkOne_pt/D");
        fivetracks_tree->Branch("dimuonDiTrkOne_eta",        &dimuonDiTrkOne_eta,        "dimuonDiTrkOne_eta/D");
        fivetracks_tree->Branch("dimuonDiTrkOne_phi",        &dimuonDiTrkOne_phi,        "dimuonDiTrkOne_phi/D");
        fivetracks_tree->Branch("dimuonDiTrkOne_charge",        &dimuonDiTrkOne_charge,        "dimuonDiTrkOne_charge/D");

        fivetracks_tree->Branch("dimuonDiTrkTwo_pt",          &dimuonDiTrkTwo_pt,          "dimuonDiTrkTwo_pt/D");
        fivetracks_tree->Branch("dimuonDiTrkTwo_eta",        &dimuonDiTrkTwo_eta,        "dimuonDiTrkTwo_eta/D");
        fivetracks_tree->Branch("dimuonDiTrkTwo_phi",        &dimuonDiTrkTwo_phi,        "dimuonDiTrkTwo_phi/D");
        fivetracks_tree->Branch("dimuonDiTrkTwo_charge",        &dimuonDiTrkTwo_charge,        "dimuonDiTrkTwo_charge/D");

        fivetracks_tree->Branch("dimuonDiTrkThree_pt",          &dimuonDiTrkThree_pt,          "dimuonDiTrkThree_pt/D");
        fivetracks_tree->Branch("dimuonDiTrkThree_eta",        &dimuonDiTrkThree_eta,        "dimuonDiTrkThree_eta/D");
        fivetracks_tree->Branch("dimuonDiTrkThree_phi",        &dimuonDiTrkThree_phi,        "dimuonDiTrkThree_phi/D");
        fivetracks_tree->Branch("dimuonDiTrkThree_charge",        &dimuonDiTrkThree_charge,        "dimuonDiTrkThree_charge/D");

        //PsiPrime
        //One
        fivetracks_tree->Branch("psiPrimeSame_pt",         &psiPrimeSame_pt,        "psiPrimeSame_pt/D");
        fivetracks_tree->Branch("psiPrimeSame_eta",        &psiPrimeSame_eta,       "psiPrimeSame_eta/D");
        fivetracks_tree->Branch("psiPrimeSame_phi",        &psiPrimeSame_phi,       "psiPrimeSame_phi/D");
        fivetracks_tree->Branch("psiPrimeSame_n",          &psiPrimeSame_n,         "psiPrimeSame_n/D");

        fivetracks_tree->Branch("psiPrimeSame_p_pt",       &psiPrimeSame_p_pt,       "psiPrimeSame_p_pt/D");
        fivetracks_tree->Branch("psiPrimeSame_p_eta",      &psiPrimeSame_p_eta,      "psiPrimeSame_p_eta/D");
        fivetracks_tree->Branch("psiPrimeSame_p_phi",      &psiPrimeSame_p_phi,      "psiPrimeSame_p_phi/D");
        fivetracks_tree->Branch("psiPrimeSame_p_n",        &psiPrimeSame_p_n,        "psiPrimeSame_p_n/D");

        fivetracks_tree->Branch("psiPrimeSame_m_pt",       &psiPrimeSame_m_pt,       "psiPrimeSame_m_pt/D");
        fivetracks_tree->Branch("psiPrimeSame_m_eta",      &psiPrimeSame_m_eta,      "psiPrimeSame_m_eta/D");
        fivetracks_tree->Branch("psiPrimeSame_m_phi",      &psiPrimeSame_m_phi,      "psiPrimeSame_m_phi/D");
        fivetracks_tree->Branch("psiPrimeSame_m_n",        &psiPrimeSame_m_n,        "psiPrimeSame_m_n/D");

        //Two
        fivetracks_tree->Branch("psiPrimeMixed_pt",         &psiPrimeMixed_pt,        "psiPrimeMixed_pt/D");
        fivetracks_tree->Branch("psiPrimeMixed_eta",        &psiPrimeMixed_eta,       "psiPrimeMixed_eta/D");
        fivetracks_tree->Branch("psiPrimeMixed_phi",        &psiPrimeMixed_phi,       "psiPrimeMixed_phi/D");
        fivetracks_tree->Branch("psiPrimeMixed_n",         &psiPrimeMixed_n,          "psiPrimeMixed_n/D");

        fivetracks_tree->Branch("psiPrimeMixed_p_pt",       &psiPrimeMixed_p_pt,       "psiPrimeMixed_p_pt/D");
        fivetracks_tree->Branch("psiPrimeMixed_p_eta",      &psiPrimeMixed_p_eta,      "psiPrimeMixed_p_eta/D");
        fivetracks_tree->Branch("psiPrimeMixed_p_phi",      &psiPrimeMixed_p_phi,      "psiPrimeMixed_p_phi/D");
        fivetracks_tree->Branch("psiPrimeMixed_p_n",        &psiPrimeMixed_p_n,        "psiPrimeMixed_p_n/D");

        fivetracks_tree->Branch("psiPrimeMixed_m_pt",       &psiPrimeMixed_m_pt,       "psiPrimeMixed_m_pt/D");
        fivetracks_tree->Branch("psiPrimeMixed_m_eta",      &psiPrimeMixed_m_eta,      "psiPrimeMixed_m_eta/D");
        fivetracks_tree->Branch("psiPrimeMixed_m_phi",      &psiPrimeMixed_m_phi,      "psiPrimeMixed_m_phi/D");
        fivetracks_tree->Branch("psiPrimeMixed_m_n",        &psiPrimeMixed_m_n,        "psiPrimeMixed_m_n/D");

        //Relative ditraks

        fivetracks_tree->Branch("psiPrimeSame_ditrak_pt",         &psiPrimeSame_ditrak_pt,        "psiPrimeSame_ditrak_pt/D");
        fivetracks_tree->Branch("psiPrimeSame_ditrak_eta",        &psiPrimeSame_ditrak_eta,       "psiPrimeSame_ditrak_eta/D");
        fivetracks_tree->Branch("psiPrimeSame_ditrak_phi",        &psiPrimeSame_ditrak_phi,       "psiPrimeSame_ditrak_phi/D");
        fivetracks_tree->Branch("psiPrimeSame_ditrak_n",          &psiPrimeSame_ditrak_n,         "psiPrimeSame_ditrak_n/D");

        fivetracks_tree->Branch("psiPrimeMixed_ditrak_pt",         &psiPrimeMixed_ditrak_pt,        "psiPrimeMixed_ditrak_pt/D");
        fivetracks_tree->Branch("psiPrimeMixed_ditrak_eta",        &psiPrimeMixed_ditrak_eta,       "psiPrimeMixed_ditrak_eta/D");
        fivetracks_tree->Branch("psiPrimeMixed_ditrak_phi",        &psiPrimeMixed_ditrak_phi,       "psiPrimeMixed_ditrak_phi/D");
        fivetracks_tree->Branch("psiPrimeMixed_ditrak_n",         &psiPrimeMixed_ditrak_n,          "psiPrimeMixed_ditrak_n/D");


        //TriTrak system
        fivetracks_tree->Branch("triTrak_pt",         &triTrak_pt,         "triTrak_pt/D");
        fivetracks_tree->Branch("triTrak_eta",        &triTrak_eta,        "triTrak_eta/D");
        fivetracks_tree->Branch("triTrak_phi",        &triTrak_phi,        "triTrak_phi/D");
        fivetracks_tree->Branch("triTrak_charge",     &triTrak_charge,     "triTrak_charge/D");

        numMasses = 5;

        TLorentzVector zero;
        zero.SetPtEtaPhiM(-10.0,-10.0,-10.0,-10.0);

        for(size_t i = 0; i<numMasses;i++)
        {
          fiveTracksMass.push_back(-1.0);
          fiveTracksMassRef.push_back(-1.0);
          triTrakMass.push_back(-1.0);

          fiveTracksVNDof.push_back(-1.0);
          fiveTracksVChi2.push_back(-1.0);
          fiveTracksVProb.push_back(-1.0);

          fiveTracksCTau.push_back(-1000.0);
          fiveTracksCTauErr.push_back(-1000.0);
          fiveTracksCosAlpha.push_back(-1.1);

          dimuonDiTrkOne.push_back(-1.0);
          dimuonDiTrkTwo.push_back(-1.0);
          dimuonDiTrkThree.push_back(-1.0);

          ditrakOne.push_back(-1.0);
          ditrakTwo.push_back(-1.0);
          ditrakThree.push_back(-1.0);

          trackOneMass.push_back(-1.0);
          trackTwoMass.push_back(-1.0);
          trackThreeMass.push_back(-1.0);

          psiPrimeSame.push_back(-1.0);
          psiPrimeMixed.push_back(-1.0);
          psiPrimeSame_ditrak.push_back(-1.0);
          psiPrimeMixed_ditrak.push_back(-1.0);

          psiPrimeSame_p_m.push_back(-1.0);
          psiPrimeSame_m_m.push_back(-1.0);
          psiPrimeMixed_p_m.push_back(-1.0);
          psiPrimeMixed_m_m.push_back(-1.0);

          five_p4.push_back(zero);
          // five_ref_p4.push_back(zero);

          dimuonDiTrkOne_p4.push_back(zero);
          dimuonDiTrkTwo_p4.push_back(zero);
          dimuonDiTrkThree_p4.push_back(zero);

          ditrakOne_p4.push_back(zero);
          ditrakTwo_p4.push_back(zero);
          ditrakThree_p4.push_back(zero);

          psiPrimeSame_p4.push_back(zero);
          psiPrimeMixed_p4.push_back(zero);

          psiPrimeSame_ditrak_p4.push_back(zero);
          psiPrimeMixed_ditrak_p4.push_back(zero);

        }

        std::vector < std::string > refNames;
        refNames.push_back("kkk"); refNames.push_back("ppk"); refNames.push_back("kpp");
        refNames.push_back("pkp"); refNames.push_back("ppp");


        for(size_t i = 0; i<numMasses;i++)
        {

          std::string name = "mass_" + refNames[i];
          std::string var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&fiveTracksMass[i],var.c_str());
          name = "mass_ref_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&fiveTracksMassRef[i],var.c_str());

          name = "vProb_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&fiveTracksVProb[i],var.c_str());
          name = "nDof_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&fiveTracksVNDof[i],var.c_str());
          name = "vChi2_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&fiveTracksVChi2[i],var.c_str());


          name = "ctau_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&fiveTracksCTau[i],var.c_str());
          name = "ctauErr_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&fiveTracksCTauErr[i],var.c_str());
          name = "cosAlpha_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&fiveTracksCosAlpha[i],var.c_str());


          name = "onePsiPrime_m_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&psiPrimeSame[i],var.c_str());
          name = "twoPsiPrime_m_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&psiPrimeMixed[i],var.c_str());

          name = "onePsiPrime_p_m" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&psiPrimeSame_p_m[i],var.c_str());
          name = "onePsiPrime_m_m" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&psiPrimeSame_m_m[i],var.c_str());


          name = "twoPsiPrime_p_m" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&psiPrimeMixed_p_m[i],var.c_str());
          name = "twoPsiPrime_m_m" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&psiPrimeMixed_m_m[i],var.c_str());

          name = "dimuonDiTrkOne_m_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&dimuonDiTrkOne[i],var.c_str());
          name = "dimuonDiTrkTwo_m_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&dimuonDiTrkTwo[i],var.c_str());
          name = "dimuonDiTrkThree_m_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&dimuonDiTrkThree[i],var.c_str());

          name = "ditrakOne_m_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&ditrakOne[i],var.c_str());
          name = "ditrakTwo_m_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&ditrakTwo[i],var.c_str());
          name = "ditrakThree_m_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&ditrakThree[i],var.c_str());

          name = "trackOne_m_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&trackOneMass[i],var.c_str());
          name = "trackTwo_m_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&trackTwoMass[i],var.c_str());
          name = "trackThree_m_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&trackThreeMass[i],var.c_str());

          name = "five_p4_" + refNames[i];
          fivetracks_tree->Branch(name.c_str(),"TLorentzVector",&five_p4[i]);
          // name = "five_ref_p4_" + refNames[i];
          // fivetracks_tree->Branch(name.c_str(),"TLorentzVector",&five_ref_p4[i]);

          name = "dimuonDiTrkOne_p4_" + refNames[i];
          fivetracks_tree->Branch(name.c_str(),"TLorentzVector",&dimuonDiTrkOne_p4[i]);
          name = "dimuonDiTrkTwo_p4_" + refNames[i];
          fivetracks_tree->Branch(name.c_str(),"TLorentzVector",&dimuonDiTrkTwo_p4[i]);
          name = "dimuonDiTrkThree_p4_" + refNames[i];
          fivetracks_tree->Branch(name.c_str(),"TLorentzVector",&dimuonDiTrkThree_p4[i]);

          name = "ditrakOne_p4_" + refNames[i];
          fivetracks_tree->Branch(name.c_str(),"TLorentzVector",&ditrakOne_p4[i]);
          name = "ditrakTwo_p4_" + refNames[i];
          fivetracks_tree->Branch(name.c_str(),"TLorentzVector",&ditrakTwo_p4[i]);
          name = "ditrakThree_p4_" + refNames[i];
          fivetracks_tree->Branch(name.c_str(),"TLorentzVector",&ditrakThree_p4[i]);

          name = "psiPrimeSame_p4_" + refNames[i];
          fivetracks_tree->Branch(name.c_str(),"TLorentzVector",&psiPrimeSame_p4[i]);
          name = "psiPrimeMixed_p4_" + refNames[i];
          fivetracks_tree->Branch(name.c_str(),"TLorentzVector",&psiPrimeMixed_p4[i]);

          name = "psiPrimeSame_ditrak_p4_" + refNames[i];
          fivetracks_tree->Branch(name.c_str(),"TLorentzVector",&psiPrimeSame_ditrak_p4[i]);
          name = "psiPrimeMixed_ditrak_p4_" + refNames[i];
          fivetracks_tree->Branch(name.c_str(),"TLorentzVector",&psiPrimeMixed_ditrak_p4[i]);

          name = "triTrak_m_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&triTrakMass[i],var.c_str());

        }

        genCands_ = consumes< std::vector <reco::GenParticle> >((edm::InputTag)"prunedGenParticles");
        packCands_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");

}

FiveTracksRootupler::~FiveTracksRootupler() {}

//
// member functions
//

bool FiveTracksRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   for (size_t i=0; i< particle->numberOfMothers(); i++) {
      if (isAncestor(ancestor, particle->mother(i))) return true;
   }
   return false;
}


// ------------ method called for each event  ------------
void FiveTracksRootupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
//  using namespace edm;
  using namespace std;

  edm::Handle<std::vector<pat::CompositeCandidate>> fivetracks_cand_handle;
  iEvent.getByToken(fivetracks_cand_Label, fivetracks_cand_handle);

  edm::Handle<std::vector<reco::Vertex >> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( triggerResults_Label , triggerResults_handle);

  numPrimaryVertices = primaryVertices_handle->size();
  run   = iEvent.id().run();
  event = iEvent.id().event();
  lumi  = iEvent.id().luminosityBlock();

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

edm::Handle< std::vector <reco::GenParticle> > pruned;
iEvent.getByToken(genCands_, pruned);

// Packed particles are all the status 1. The navigation to pruned is possible (the other direction should be made by hand)
edm::Handle<pat::PackedGenParticleCollection> packed;
iEvent.getByToken(packCands_,  packed);

//
// if ( motherInPrunedCollection != nullptr && (d->pdgId() ==  13 ) && isAncestor(aditrkdimu , motherInPrunedCollection) ) {
//   gen_muonn_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
//   foundit++;
// }
// if ( motherInPrunedCollection != nullptr && (d->pdgId() == -13 ) && isAncestor(aditrkdimu , motherInPrunedCollection) ) {
//   gen_muonp_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
//   foundit++;
// }

gen_dimuonditrk_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_ditrak_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_muonp_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_muonn_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_kaonp_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_kaonn_p4.SetPtEtaPhiM(0.,0.,0.,0.);

gen_b_p4.SetPtEtaPhiM(0.,0.,0.,0.);

gen_b4_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_d1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_d2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_gd1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_gd2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_gd3_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_gd4_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_gd5_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_gd6_p4.SetPtEtaPhiM(0.,0.,0.,0.);

gen_dimuonditrk_pdgId = 0;

//std::cout << "Debug  1" << std::endl;
if(!OnlyGen_)
  if (!fivetracks_cand_handle.isValid()) std::cout<< "No five tracks information " << run << "," << event <<std::endl;

  if (fivetracks_cand_handle.isValid()) {


    //std::cout << "Debug  2" << std::endl;
    noFiveCandidates = (Int_t)(fivetracks_cand_handle->size());

    for (unsigned int i=0; i< fivetracks_cand_handle->size(); i++)
    {

      pat::CompositeCandidate five_cand;
      const pat::PackedCandidate *trakOne_cand, *trakTwo_cand, *trakThree_cand;
      const pat::CompositeCandidate *dimuonDiTrkOne_cand, *dimuonDiTrkTwo_cand, *dimuonDiTrkThree_cand, *first_five_ref;
      const pat::CompositeCandidate *dimuonditrk_cand, *dimuon_cand, *ditrakOne_cand;
      const pat::CompositeCandidate *triTrak_cand, *ditrakTwo_cand, *ditrakThree_cand;

      five_cand  = fivetracks_cand_handle->at(i);
      dimuonditrk_id = five_cand.userInt("dimuontt_index");

      std::string name = "fiveCand_" + std::to_string(0);

      first_five_ref     = dynamic_cast<const pat::CompositeCandidate*>(five_cand.daughter(name));

      dimuonditrk_cand = dynamic_cast<const pat::CompositeCandidate*>(five_cand.daughter("dimuonditrak"));
      dimuon_cand = dynamic_cast<const pat::CompositeCandidate*>(dimuonditrk_cand->daughter("dimuon"));
      ditrakOne_cand = dynamic_cast<const pat::CompositeCandidate*>(dimuonditrk_cand->daughter("ditrak"));

      highKaonMatch  = (Double_t)dimuonditrk_cand->userInt("highKaonMatch");
      lowKaonMatch   = (Double_t)dimuonditrk_cand->userInt("lowKaonMatch");
      lowMuonMatch   = (Double_t)dimuon_cand->userInt("highMuonTMatch");
      highMuonMatch  = (Double_t)dimuon_cand->userInt("lowMuonTMatch");

      trakOne_cand = dynamic_cast<const pat::PackedCandidate*>(five_cand.daughter("trakOne"));
      trakTwo_cand = dynamic_cast<const pat::PackedCandidate*>(five_cand.daughter("trakTwo"));
      trakThree_cand = dynamic_cast<const pat::PackedCandidate*>(five_cand.daughter("trakThree"));

      dimuonDiTrkOne_cand = dynamic_cast<const pat::CompositeCandidate*>(first_five_ref->daughter("dimuonDiTrakOne"));
      dimuonDiTrkTwo_cand = dynamic_cast<const pat::CompositeCandidate*>(first_five_ref->daughter("dimuonDiTrakTwo"));
      dimuonDiTrkThree_cand = dynamic_cast<const pat::CompositeCandidate*>(first_five_ref->daughter("dimuonDiTrakThree"));

      ditrakTwo_cand    = dynamic_cast<const pat::CompositeCandidate*>(dimuonDiTrkTwo_cand->daughter("ditrak"));
      ditrakThree_cand  = dynamic_cast<const pat::CompositeCandidate*>(dimuonDiTrkThree_cand->daughter("ditrak"));

      triTrak_cand = dynamic_cast<const pat::CompositeCandidate*>(first_five_ref->daughter("triTrak"));

      dimuonditrk_p4.SetPtEtaPhiM(dimuonditrk_cand->pt(),dimuonditrk_cand->eta(),dimuonditrk_cand->phi(),dimuonditrk_cand->mass());
      dimuonditrk_m = dimuonditrk_cand->mass();
      dimuonditrk_pt = dimuonditrk_cand->pt();
      dimuonditrk_eta = dimuonditrk_cand->eta();
      dimuonditrk_phi = dimuonditrk_cand->phi();
      dimuonditrk_p = dimuonditrk_cand->p();

      dimuon_p4.SetPtEtaPhiM(dimuon_cand->pt(),dimuon_cand->eta(),dimuon_cand->phi(),dimuon_cand->mass());
      dimuon_m = dimuon_cand->mass();
      dimuon_pt = dimuon_cand->pt();
      dimuon_eta = dimuon_cand->eta();
      dimuon_phi = dimuon_cand->phi();
      dimuon_p = dimuon_cand->p();

      ditrak_p4.SetPtEtaPhiM(ditrakOne_cand->pt(),ditrakOne_cand->eta(),ditrakOne_cand->phi(),ditrakOne_cand->mass());
      ditrak_m = ditrakOne_cand->mass();
      ditrakOne_pt = ditrakOne_cand->pt();
      ditrakOne_eta = ditrakOne_cand->eta();
      ditrakOne_phi = ditrakOne_cand->phi();
      ditrakOne_p = ditrakOne_cand->p();

      ditrakTwo_pt = ditrakTwo_cand->pt();
      ditrakTwo_eta = ditrakTwo_cand->eta();
      ditrakTwo_phi = ditrakTwo_cand->phi();
      ditrakTwo_p = ditrakTwo_cand->p();

      ditrakThree_pt = ditrakThree_cand->pt();
      ditrakThree_eta = ditrakThree_cand->eta();
      ditrakThree_phi = ditrakThree_cand->phi();
      ditrakThree_p = ditrakThree_cand->p();

      //trakOne_cand.SetPtEtaPhiM(dimuon_cand->pt(),dimuon_cand->eta(),dimuon_cand->phi(),dimuon_cand->mass());
      highTrack_pt = trakOne_cand->pt();
      highTrack_eta = trakOne_cand->eta();
      highTrack_phi = trakOne_cand->phi();
      highTrack_charge = trakOne_cand->charge();

      lowTrack_pt = trakTwo_cand->pt();
      lowTrack_eta = trakTwo_cand->eta();
      lowTrack_phi = trakTwo_cand->phi();
      lowTrack_charge = trakTwo_cand->charge();

      thirdTrack_pt     = trakThree_cand->pt();
      thirdTrack_eta    = trakThree_cand->eta();
      thirdTrack_phi    = trakThree_cand->phi();
      thirdTrack_charge = trakThree_cand->charge();
      thirdTrack_dz     = trakThree_cand->bestTrack()->dz();
      thirdTrack_dxy    = trakThree_cand->bestTrack()->dxy();

      dimuonDiTrkOne_pt     = dimuonDiTrkOne_cand->pt();
      dimuonDiTrkOne_eta    = dimuonDiTrkOne_cand->eta();
      dimuonDiTrkOne_phi    = dimuonDiTrkOne_cand->phi();
      dimuonDiTrkOne_charge = dimuonDiTrkOne_cand->charge();

      dimuonDiTrkTwo_pt     = dimuonDiTrkTwo_cand->pt();
      dimuonDiTrkTwo_eta    = dimuonDiTrkTwo_cand->eta();
      dimuonDiTrkTwo_phi    = dimuonDiTrkTwo_cand->phi();
      dimuonDiTrkTwo_charge = dimuonDiTrkTwo_cand->charge();

      dimuonDiTrkThree_pt     = dimuonDiTrkThree_cand->pt();
      dimuonDiTrkThree_eta    = dimuonDiTrkThree_cand->eta();
      dimuonDiTrkThree_phi    = dimuonDiTrkThree_cand->phi();
      dimuonDiTrkThree_charge = dimuonDiTrkThree_cand->charge();

      triTrak_pt     = triTrak_cand->pt();
      triTrak_eta    = triTrak_cand->eta();
      triTrak_phi    = triTrak_cand->phi();
      triTrak_charge = triTrak_cand->charge();

      psiPrimeSame_pt  = dimuonDiTrkOne_cand->pt();
      psiPrimeSame_eta = dimuonDiTrkOne_cand->eta();
      psiPrimeSame_phi = dimuonDiTrkOne_cand->phi();
      psiPrimeSame_ditrak_n   = 1.0;
      psiPrimeSame_ditrak_pt  = dimuonDiTrkOne_cand->pt();
      psiPrimeSame_ditrak_eta = dimuonDiTrkOne_cand->eta();
      psiPrimeSame_ditrak_phi = dimuonDiTrkOne_cand->phi();
      //std::cout << "Debug  3" << std::endl;


      for(size_t j = 0; j<numMasses;j++)
      {

        const pat::CompositeCandidate* five_cand_ref, *five_cand_ref_ref, *triTrak_cand_ref;
        const pat::CompositeCandidate *trakOne_cand_ref, *trakTwo_cand_ref, *trakThree_cand_ref;

        const pat::CompositeCandidate *dimuonDiTrkOne_cand_ref, *dimuonDiTrkTwo_cand_ref, *dimuonDiTrkThree_cand_ref;
        const pat::CompositeCandidate *psiPrimeMixed_cand;

        const pat::CompositeCandidate *ditrakOne_cand_ref, *ditrakTwo_cand_ref, *ditrakThree_cand_ref;
        // const pat::CompositeCandidate *psiPrimeMixed_ditrak_cand;

        name = "fiveCand_" + std::to_string(j);

        five_cand_ref     = dynamic_cast<const pat::CompositeCandidate*>(five_cand.daughter(name));
        //
        // if(five_cand_ref->userFloat("has_ref")>0.0)
        // {
        //   five_cand_ref_ref = dynamic_cast<const pat::CompositeCandidate*>(five_cand_ref->daughter("fiveRef"));
        //   five_ref_p4[j].SetPtEtaPhiM(five_cand_ref_ref->pt(),five_cand_ref_ref->eta(),five_cand_ref_ref->phi(),five_cand_ref_ref->mass());
        // }

        fiveTracksMass[j] = five_cand_ref->mass();

        fiveTracksMassRef[j] = five_cand_ref->userFloat("mass_ref");

        fiveTracksVProb[j] = five_cand_ref->userFloat("vProb");
        fiveTracksVNDof[j] = five_cand_ref->userFloat("nDof");
        fiveTracksVChi2[j] = five_cand_ref->userFloat("vChi2");

        fiveTracksCTau[j] = five_cand_ref->userFloat("ctauPV");
        fiveTracksCTauErr[j] = five_cand_ref->userFloat("ctauErrPV");
        fiveTracksCosAlpha[j] = five_cand_ref->userFloat("cosAlpha");

        trakOne_cand_ref = dynamic_cast<const pat::CompositeCandidate*>(five_cand_ref->daughter("trakOne"));
        trakTwo_cand_ref = dynamic_cast<const pat::CompositeCandidate*>(five_cand_ref->daughter("trakTwo"));
        trakThree_cand_ref = dynamic_cast<const pat::CompositeCandidate*>(five_cand_ref->daughter("trakThree"));

        dimuonDiTrkOne_cand_ref = dynamic_cast<const pat::CompositeCandidate*>(five_cand_ref->daughter("dimuonDiTrakOne"));
        dimuonDiTrkTwo_cand_ref = dynamic_cast<const pat::CompositeCandidate*>(five_cand_ref->daughter("dimuonDiTrakTwo"));
        dimuonDiTrkThree_cand_ref = dynamic_cast<const pat::CompositeCandidate*>(five_cand_ref->daughter("dimuonDiTrakThree"));

        ditrakOne_cand_ref = dynamic_cast<const pat::CompositeCandidate*>(dimuonDiTrkOne_cand_ref->daughter("ditrak"));
        ditrakTwo_cand_ref = dynamic_cast<const pat::CompositeCandidate*>(dimuonDiTrkTwo_cand_ref->daughter("ditrak"));
        ditrakThree_cand_ref = dynamic_cast<const pat::CompositeCandidate*>(dimuonDiTrkThree_cand_ref->daughter("ditrak"));

        triTrak_cand_ref = dynamic_cast<const pat::CompositeCandidate*>(five_cand_ref->daughter("triTrak"));

        dimuonDiTrkOne[j]  = dimuonDiTrkOne_cand_ref->mass();
        dimuonDiTrkTwo[j]  = dimuonDiTrkTwo_cand_ref->mass();
        dimuonDiTrkThree[j]  = dimuonDiTrkThree_cand_ref->mass();

        ditrakOne[j]  = ditrakOne_cand_ref->mass();
        ditrakTwo[j]  = ditrakTwo_cand_ref->mass();
        ditrakThree[j]  = ditrakThree_cand_ref->mass();

        trackOneMass[j]  = trakOne_cand_ref->mass();
        trackTwoMass[j]  = trakTwo_cand_ref->mass();
        trackThreeMass[j]  = trakThree_cand_ref->mass();

        five_p4[j].SetPtEtaPhiM(five_cand_ref->pt(),five_cand_ref->eta(),five_cand_ref->phi(),five_cand_ref->mass());

        dimuonDiTrkOne_p4[j].SetPtEtaPhiM(dimuonDiTrkOne_cand_ref->pt(),dimuonDiTrkOne_cand_ref->eta(),dimuonDiTrkOne_cand_ref->phi(),dimuonDiTrkOne_cand_ref->mass());
        dimuonDiTrkTwo_p4[j].SetPtEtaPhiM(dimuonDiTrkTwo_cand_ref->pt(),dimuonDiTrkTwo_cand_ref->eta(),dimuonDiTrkTwo_cand_ref->phi(),dimuonDiTrkTwo_cand_ref->mass());
        dimuonDiTrkThree_p4[j].SetPtEtaPhiM(dimuonDiTrkThree_cand_ref->pt(),dimuonDiTrkThree_cand_ref->eta(),dimuonDiTrkThree_cand_ref->phi(),dimuonDiTrkThree_cand_ref->mass());

        ditrakOne_p4[j].SetPtEtaPhiM(ditrakOne_cand_ref->pt(),ditrakOne_cand_ref->eta(),ditrakOne_cand_ref->phi(),ditrakOne_cand_ref->mass());
        ditrakTwo_p4[j].SetPtEtaPhiM(ditrakTwo_cand_ref->pt(),ditrakTwo_cand_ref->eta(),ditrakTwo_cand_ref->phi(),ditrakTwo_cand_ref->mass());
        ditrakThree_p4[j].SetPtEtaPhiM(ditrakThree_cand_ref->pt(),ditrakThree_cand_ref->eta(),ditrakThree_cand_ref->phi(),ditrakThree_cand_ref->mass());
        //std::cout << "Debug  5" << std::endl;
        //////////////////////////
        ////PSI PRIME SAME
        // My PsiPrime Same come from the original Phi->kk traks + J/Psi
        // Is simply the dimuonditrk refit with ππ πk kπ masses
        // psiPrimeSame_cand == dimuonDiTrkOne_cand_ref;

        const pat::CompositeCandidate *psiPrimeSame_m_cand_ref, *psiPrimeSame_p_cand_ref;

        psiPrimeSame_p4[j].SetPtEtaPhiM(dimuonDiTrkOne_cand_ref->pt(),dimuonDiTrkOne_cand_ref->eta(),dimuonDiTrkOne_cand_ref->phi(),dimuonDiTrkOne_cand_ref->mass());
        psiPrimeSame[j] = dimuonDiTrkOne_cand_ref->mass();

        psiPrimeSame_ditrak_p4[j].SetPtEtaPhiM(ditrakOne_cand_ref->pt(),ditrakOne_cand_ref->eta(),ditrakOne_cand_ref->phi(),ditrakOne_cand_ref->mass());
        psiPrimeSame_ditrak[j]  = ditrakOne_cand_ref->mass();

        if(highTrack_charge>0)
        {
          psiPrimeSame_p_cand_ref  = dynamic_cast<const pat::CompositeCandidate*>(five_cand_ref->daughter("trakOne"));
          psiPrimeSame_m_cand_ref = dynamic_cast<const pat::CompositeCandidate*>(five_cand_ref->daughter("trakTwo"));

          psiPrimeSame_p_n = 1.0;
          psiPrimeSame_m_n = 2.0;

          psiPrimeSame_p_m[j] = trakOne_cand_ref->mass();
          psiPrimeSame_m_m[j] = trakTwo_cand_ref->mass();
        }
        else
        {
          psiPrimeSame_p_cand_ref = dynamic_cast<const pat::CompositeCandidate*>(five_cand_ref->daughter("trakTwo"));
          psiPrimeSame_m_cand_ref  = dynamic_cast<const pat::CompositeCandidate*>(five_cand_ref->daughter("trakOne"));

          psiPrimeSame_p_n = 2.0;
          psiPrimeSame_m_n = 1.0;

          psiPrimeSame_p_m[j] = trakTwo_cand_ref->mass();
          psiPrimeSame_m_m[j] = trakOne_cand_ref->mass();
        }
        //std::cout << "Debug  6" << std::endl;
        psiPrimeSame_p_pt   = psiPrimeSame_p_cand_ref->pt();
        psiPrimeSame_p_eta  = psiPrimeSame_p_cand_ref->eta();
        psiPrimeSame_p_phi  = psiPrimeSame_p_cand_ref->phi();

        psiPrimeSame_m_pt   = psiPrimeSame_m_cand_ref->pt();
        psiPrimeSame_m_eta  = psiPrimeSame_m_cand_ref->eta();
        psiPrimeSame_m_phi  = psiPrimeSame_m_cand_ref->phi();

        //////////////////////////
        ////PSI PRIME MIXED
        // This depends on the third track charge

        const pat::CompositeCandidate *psiPrimeMixed_m_cand_ref, *psiPrimeMixed_p_cand_ref;
        if(dimuonDiTrkTwo_cand_ref->charge()==0)
        {
          psiPrimeMixed_p4[j].SetPtEtaPhiM(dimuonDiTrkTwo_cand_ref->pt(),dimuonDiTrkTwo_cand_ref->eta(),dimuonDiTrkTwo_cand_ref->phi(),dimuonDiTrkTwo_cand_ref->mass());
          psiPrimeMixed[j]  = dimuonDiTrkTwo_cand_ref->mass();

          psiPrimeMixed_n = 2.0;

          if(thirdTrack_charge>0)
          {
            psiPrimeMixed_p_pt   = trakThree_cand_ref->pt();
            psiPrimeMixed_p_eta  = trakThree_cand_ref->eta();
            psiPrimeMixed_p_phi  = trakThree_cand_ref->phi();

            psiPrimeMixed_m_pt   = trakOne_cand_ref->pt();
            psiPrimeMixed_m_eta  = trakOne_cand_ref->eta();
            psiPrimeMixed_m_phi  = trakOne_cand_ref->phi();

            psiPrimeMixed_p_n = 3.0;
            psiPrimeMixed_m_n = 1.0;

            psiPrimeMixed_p_m[j] = trakThree_cand_ref->mass();
            psiPrimeMixed_m_m[j] = trakOne_cand_ref->mass();
          }
          else
          {
            psiPrimeMixed_p_pt   = trakOne_cand_ref->pt();
            psiPrimeMixed_p_eta  = trakOne_cand_ref->eta();
            psiPrimeMixed_p_phi  = trakOne_cand_ref->phi();

            psiPrimeMixed_m_pt   = trakThree_cand_ref->pt();
            psiPrimeMixed_m_eta  = trakThree_cand_ref->eta();
            psiPrimeMixed_m_phi  = trakThree_cand_ref->phi();

            psiPrimeMixed_p_n = 1.0;
            psiPrimeMixed_m_n = 3.0;

            psiPrimeMixed_p_m[j] = trakOne_cand_ref->mass();
            psiPrimeMixed_m_m[j] = trakThree_cand_ref->mass();
          }
        }
        else
        if(dimuonDiTrkThree_cand_ref->charge()==0)
        {

          psiPrimeMixed_p4[j].SetPtEtaPhiM(ditrakThree_cand_ref->pt(),ditrakThree_cand_ref->eta(),ditrakThree_cand_ref->phi(),ditrakThree_cand_ref->mass());
          psiPrimeMixed[j]  = ditrakThree_cand_ref->mass();

          psiPrimeMixed_n = 3.0;

          if(thirdTrack_charge>0)
          {
            psiPrimeMixed_p_pt   = trakThree_cand_ref->pt();
            psiPrimeMixed_p_eta  = trakThree_cand_ref->eta();
            psiPrimeMixed_p_phi  = trakThree_cand_ref->phi();

            psiPrimeMixed_m_pt   = trakTwo_cand_ref->pt();
            psiPrimeMixed_m_eta  = trakTwo_cand_ref->eta();
            psiPrimeMixed_m_phi  = trakTwo_cand_ref->phi();

            psiPrimeMixed_p_n = 3.0;
            psiPrimeMixed_m_n = 2.0;

            psiPrimeMixed_p_m[j] = trakThree_cand_ref->mass();
            psiPrimeMixed_m_m[j] = trakTwo_cand_ref->mass();

          }
          else
          {
            psiPrimeMixed_p_pt   = trakTwo_cand_ref->pt();
            psiPrimeMixed_p_eta  = trakTwo_cand_ref->eta();
            psiPrimeMixed_p_phi  = trakTwo_cand_ref->phi();

            psiPrimeMixed_m_pt   = trakThree_cand_ref->pt();
            psiPrimeMixed_m_eta  = trakThree_cand_ref->eta();
            psiPrimeMixed_m_phi  = trakThree_cand_ref->phi();

            psiPrimeMixed_p_n = 2.0;
            psiPrimeMixed_m_n = 3.0;

            psiPrimeMixed_p_m[j] = trakTwo_cand_ref->mass();
            psiPrimeMixed_m_m[j] = trakThree_cand_ref->mass();
          }

        }
        else
        if(thirdTrack_charge==0) //K0, in this case the PsiPrimeMixed is just a copy of PsiPrimeSame
        {

          psiPrimeMixed_p4[j].SetPtEtaPhiM(dimuonDiTrkOne_cand_ref->pt(),dimuonDiTrkOne_cand_ref->eta(),dimuonDiTrkOne_cand_ref->phi(),dimuonDiTrkOne_cand_ref->mass());
          psiPrimeMixed[j]  = dimuonDiTrkOne_cand_ref->mass();
          psiPrimeMixed_n = 1.0;

          psiPrimeMixed_ditrak_p4[j].SetPtEtaPhiM(ditrakOne_cand_ref->pt(),ditrakOne_cand_ref->eta(),ditrakOne_cand_ref->phi(),ditrakOne_cand_ref->mass());
          psiPrimeMixed_ditrak[j]  = ditrakOne_cand_ref->mass();


          if(highTrack_charge>0)
          {

            psiPrimeMixed_m_n = 2.0;
            psiPrimeMixed_p_n = 1.0;

            psiPrimeMixed_p_pt   = trakOne_cand_ref->pt();
            psiPrimeMixed_p_eta  = trakOne_cand_ref->eta();
            psiPrimeMixed_p_phi  = trakOne_cand_ref->phi();

            psiPrimeMixed_m_pt   = trakTwo_cand_ref->pt();
            psiPrimeMixed_m_eta  = trakTwo_cand_ref->eta();
            psiPrimeMixed_m_phi  = trakTwo_cand_ref->phi();

          }
          else
          {


            psiPrimeMixed_p_n = 1.0;
            psiPrimeMixed_m_n = 2.0;

            psiPrimeMixed_p_pt   = trakTwo_cand_ref->pt();
            psiPrimeMixed_p_eta  = trakTwo_cand_ref->eta();
            psiPrimeMixed_p_phi  = trakTwo_cand_ref->phi();

            psiPrimeMixed_m_pt   = trakOne_cand_ref->pt();
            psiPrimeMixed_m_eta  = trakOne_cand_ref->eta();
            psiPrimeMixed_m_phi  = trakOne_cand_ref->phi();

          }

        }


        psiPrimeMixed_pt  = psiPrimeMixed_p4[j].Pt();
        psiPrimeMixed_eta = psiPrimeMixed_p4[j].Eta();
        psiPrimeMixed_phi = psiPrimeMixed_p4[j].Phi();
        //std::cout << "Debug  10" << std::endl;


        triTrakMass[j] = triTrak_cand_ref->mass();
        //std::cout << "Debug  11" << std::endl;
      }

      fivetracks_tree->Fill();

      if (OnlyBest_) break;
      else if(i==0)
      isBestCandidate = false;

        // dimuontt candidates are sorted by vProb
    }

  }

}

// ------------ method called once each job just before starting event loop  ------------
void FiveTracksRootupler::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void FiveTracksRootupler::endJob() {}

// ------------ method called when starting to processes a run  ------------
void FiveTracksRootupler::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void FiveTracksRootupler::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void FiveTracksRootupler::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void FiveTracksRootupler::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void FiveTracksRootupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FiveTracksRootupler);
