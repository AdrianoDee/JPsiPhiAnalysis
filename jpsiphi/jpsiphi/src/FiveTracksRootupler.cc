/*
   Package:    FiveTracksRootupler
   Class:      FiveTracksRootupler

   Description: make rootuple of DiMuon-DiTrack combination

   Original Author: Adriano Di Florio

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
  edm::EDGetTokenT<pat::CompositeCandidateCollection> FiveTracksCollection_;
  edm::EDGetTokenT<edm::TriggerResults> TriggerResults_;
  edm::EDGetTokenT<reco::VertexCollection> thePVs_;
  int  dimuonditrk_pdgid_, dimuon_pdgid_, track_pdgid_, pdgid_;
  bool IsMC_,OnlyGen_ ;
  std::vector<std::string>  HLTs_;
  std::vector<std::string>  HLTFilters_;
  std::string TreeName_;
  std::vector < Double_t > fiveTracksMass, fiveTracksMassRef, triTrackMass;
  std::vector < Double_t > fiveTracksVProb, fiveTracksVNDof, fiveTracksVChi2;
  std::vector < Double_t > fiveTracksCosAlpha, fiveTracksCTau, fiveTracksCTauErr;
  std::vector < Double_t > psiPrimeSame, psiPrimeMixed;
  std::vector < Double_t > psiPrimeSame_ditrack, psiPrimeMixed_ditrack;
  std::vector < Double_t > dimuonDiTrkOne, dimuonDiTrkTwo, dimuonDiTrkThree;
  std::vector < Double_t > ditrackOne, ditrackTwo, ditrackThree;
  std::vector < Double_t > trackOneMass, trackTwoMass, trackThreeMass;
  std::vector < Double_t > psiPrimeSame_p_m, psiPrimeSame_m_m, psiPrimeMixed_p_m, psiPrimeMixed_m_m;

  Double_t testMass;

  UInt_t run, event, lumi, numPrimaryVertices, trigger;
  UInt_t dimuon_id, dimuonditrk_id, p_id, m_id, t_id;

  TLorentzVector dimuonditrk_p4, dimuon_p4, ditrack_p4, five_p4;
  TLorentzVector lowPion_p4, lowProton_p4, highProton_p4, highPion_p4, thirdProton_p4, thirdPion_p4;
  TLorentzVector lowMuon_p4, highMuon_p4, lowKaon_p4, thirdKaon_p4, highKaon_p4;

  TLorentzVector dimuonditrk_not_rf_p4;
  TLorentzVector ditrackOne_rf_p4, ditrackOne_not_rf_p4;
  TLorentzVector muonp_rf_p4;
  TLorentzVector muonn_rf_p4;

  TLorentzVector dimuonditrk_rf_p4;
  TLorentzVector dimuonditrk_rf_const_p4;
  TLorentzVector dimuon_rf_p4, dimuon_not_rf_p4;
  TLorentzVector ditrack_rf_p4, ditrack_not_rf_p4;
  TLorentzVector lowMuon_rf_p4;
  TLorentzVector highMuon_rf_p4;
  TLorentzVector kaonp_rf_p4;
  TLorentzVector kaonn_rf_p4;

  Double_t five_m, five_m_ref, five_mass_ppk, five_mass_kpp, five_mass_pkp, five_mass_ppp, five_pt, five_eta, five_phi, five_p;
  Double_t five_cosAlpha, five_ctauPV, five_ctauErrPV, five_cosAlphaCA, five_ctauPVCA, five_ctauErrPVCA, five_cosAlphaDZ;
  Double_t five_ctauPVDZ, five_ctauErrPVDZ, five_cosAlphaBS, five_ctauPVBS, five_ctauErrPVBS;
  Double_t five_vProb, five_nDof, five_vChi2, five_vx, five_vy, five_vz;

  Int_t dimuonditrk_charge, five_charge;

  Double_t highTrackMatch, lowTrackMatch, thirdTrackMatch, highMuonMatch, lowMuonMatch;

  Double_t dimuonditrk_vProb, dimuonditrk_vChi2;
  Double_t dimuonditrk_pt, dimuonditrk_eta, dimuonditrk_phi, dimuonditrk_y, dimuonditrk_vx, dimuonditrk_vy, dimuonditrk_vz, dimuonditrk_p;
  Double_t highKaon_pt,lowKaon_pt,highMuon_pt,lowMuon_pt;

  Double_t highTrackMuonDR, highTrackMuonDP, highTrackMuonDPt;
  Double_t lowTrackMuonDR, lowTrackMuonDP, lowTrackMuonDPt;
  Double_t thirdTrackMuonDR, thirdTrackMuonDP, thirdTrackMuonDPt;

  Double_t highMuon_eta, lowMuon_eta, highMuon_phi, lowMuon_phi, dimuonditrk_m;
  Double_t highMuon_dz, lowMuon_dz, highMuon_dxy, lowMuon_dxy;
  Double_t highMuon_charge, lowMuon_charge;

  UInt_t lowMuon_NPixelHits, lowMuon_NStripHits, lowMuon_NTrackhits, lowMuon_NBPixHits, lowMuon_NPixLayers, lowMuon_NTraLayers, lowMuon_NStrLayers, lowMuon_NBPixLayers;
  UInt_t highMuon_NPixelHits, highMuon_NStripHits, highMuon_NTrackhits, highMuon_NBPixHits, highMuon_NPixLayers, highMuon_NTraLayers, highMuon_NStrLayers, highMuon_NBPixLayers;

  Bool_t lowMuon_isLoose, lowMuon_isSoft, lowMuon_isMedium, lowMuon_isHighPt, lowMuon_isTight;
  Bool_t highMuon_isLoose, highMuon_isSoft, highMuon_isMedium, highMuon_isHighPt, highMuon_isTight;

  Bool_t lowMuon_isTracker, lowMuon_isGlobal, highMuon_isTracker, highMuon_isGlobal;
  UInt_t lowMuon_type, highMuon_type;

  Double_t pv_x, pv_y, pv_z, bestPV_Z, bestPV_Y, bestPV_X;
  Double_t bS_Z, bS_X, bS_Y, cosAlphaPV_Z, cosAlphaPV_Y, cosAlphaPV_X;
  Double_t zPV_Z, zPV_X, zPV_Y;

  Double_t dimuonditrk_cosAlpha, dimuonditrk_ctauPV, dimuonditrk_ctauErrPV;
  Double_t dimuonditrk_cosAlphaDZ, dimuonditrk_ctauPVDZ, dimuonditrk_ctauErrPVDZ;
  Double_t dimuonditrk_cosAlphaBS, dimuonditrk_ctauPVBS, dimuonditrk_ctauErrPVBS;
  Double_t dimuonditrk_cosAlphaCA, dimuonditrk_ctauPVCA,dimuonditrk_ctauErrPVCA;

  Double_t tPFromPVBS, tMFromPVBS, tPFromPVCA, tMFromPVCA;
  Double_t tPFromPVDZ, tMFromPVDZ, tPFromPV, tMFromPV;
  Double_t tTFromPV, tTFromPVDZ, tTFromPVBS, tTFromPVCA;

  Double_t dca_m1m2, dca_m1t1, dca_m1t2, dca_m2t1, dca_m2t2, dca_t1t2;
  Double_t dca_m1t3, dca_m2t3, dca_t1t3, dca_t2t3;
  Double_t dimuon_vProb, dimuon_vChi2, dimuon_DCA, dimuon_ctauPV, dimuon_ctauErrPV, dimuon_cosAlpha;

  Double_t highTrack_pt, highTrack_eta, highTrack_phi, highTrack_charge;
  Double_t lowTrack_pt, lowTrack_eta, lowTrack_phi, lowTrack_charge;
  Double_t thirdTrack_pt, thirdTrack_eta, thirdTrack_phi, thirdTrack_charge;
  Double_t thirdTrack_dz, thirdTrack_dxy, lowTrack_dz, lowTrack_dxy, highTrack_dz, highTrack_dxy;

  Double_t dimuonDiTrkOne_pt, dimuonDiTrkOne_eta, dimuonDiTrkOne_phi, dimuonDiTrkOne_charge, dimuonDiTrkOne_p;
  Double_t dimuonDiTrkTwo_pt, dimuonDiTrkTwo_eta, dimuonDiTrkTwo_phi, dimuonDiTrkTwo_charge, dimuonDiTrkTwo_p;
  Double_t dimuonDiTrkThree_pt, dimuonDiTrkThree_eta, dimuonDiTrkThree_phi, dimuonDiTrkThree_charge, dimuonDiTrkThree_p;

  Double_t dimuonditrk_nDof,dimuonditrk_m_rf,dimuonditrk_m_rf_c,dimuonditrk_m_rf_d_c;
  Int_t highTrack_NPixelHits, highTrack_NStripHits, highTrack_NTrackhits, highTrack_NBPixHits, highTrack_NPixLayers;
  Int_t highTrack_NTraLayers, highTrack_NStrLayers, highTrack_NBPixLayers, lowTrack_NPixelHits, lowTrack_NStripHits;
  Int_t lowTrack_NTrackhits, lowTrack_NBPixHits, lowTrack_NPixLayers, lowTrack_NTraLayers, lowTrack_NStrLayers, lowTrack_NBPixLayers;
  Int_t thirdTrack_NPixelHits, thirdTrack_NStripHits, thirdTrack_NTrackhits, thirdTrack_NBPixHits, thirdTrack_NPixLayers;
  Int_t thirdTrack_NTraLayers, thirdTrack_NStrLayers, thirdTrack_NBPixLayers;

  Double_t dimuon_m, dimuon_pt, dimuon_eta, dimuon_phi, dimuon_p;
  Double_t ditrack_m, ditrackOne_pt, ditrackOne_eta, ditrackOne_phi, ditrackOne_p;
  Double_t ditrackTwo_pt, ditrackTwo_eta, ditrackTwo_phi, ditrackTwo_p;
  Double_t ditrackThree_pt, ditrackThree_eta, ditrackThree_phi, ditrackThree_p;

  Double_t triTrack_pt, triTrack_eta, triTrack_phi, triTrack_charge, triTrack_m;

  Bool_t muonP_isLoose, muonP_isSoft, muonP_isMedium, muonP_isHighPt;
  Bool_t muonN_isLoose, muonN_isSoft, muonN_isMedium, muonN_isHighPt;

  Bool_t muonP_isTracker, muonP_isGlobal, muonN_isTracker, muonN_isGlobal;
  UInt_t muonP_type, muonN_type;

  Bool_t muhighTrackose, muonP_rf_isSoft, muonP_rf_isMedium, muonP_rf_isHighPt;
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


  //MC
  TLorentzVector gen_jpsi_p4, gen_phi_p4, gen_five_p4, gen_dimuonditrk_p4;
  TLorentzVector gen_lowMuon_p4, gen_highMuon_p4, gen_highKaon_p4, gen_lowKaon_p4, genThirdTrack_p4;

  Double_t gen_five_prompt, gen_five_pt, gen_five_p, gen_five_eta, gen_five_phi;
  Double_t gen_five_pdg, gen_phi_pdg, gen_jpsi_pdg, gen_dimuonditrk_pdg;
  Double_t gen_lowMuon_pdg, gen_highMuon_pdg, gen_highKaon_pdg, gen_lowKaon_pdg, genThirdTrack_pdg;
  Double_t gen_lowMuon_mompdg, gen_highMuon_mompdg, gen_highKaon_mompdg, gen_lowKaon_mompdg, genThirdTrack_mompdg;
  Double_t gen_lowMuon_status, gen_highMuon_status, gen_highKaon_status, gen_lowKaon_status, genThirdTrack_status;
  Double_t gen_dimuonditrk_prompt, gen_phi_prompt, gen_jpsi_prompt;
  Double_t gen_dimuonditrk_pt, gen_dimuonditrk_p, gen_dimuonditrk_eta;
  Double_t gen_phi_pt, gen_phi_p, gen_phi_eta;
  Double_t gen_jpsi_pt, gen_jpsi_p, gen_jpsi_eta;
  Double_t gen_lowMuon_phi, gen_highMuon_phi, gen_highKaon_phi, gen_lowKaon_phi, genThirdTrack_phi;
  Double_t gen_dimuonditrk_phi, gen_phi_phi, gen_jpsi_phi;
  Double_t gen_lowMuon_p, gen_highMuon_p, gen_highKaon_p, gen_lowKaon_p, genThirdTrack_p;
  Double_t gen_lowMuon_pt, gen_highMuon_pt, gen_highKaon_pt, gen_lowKaon_pt, genThirdTrack_pt;
  Double_t gen_lowMuon_eta, gen_highMuon_eta, gen_highKaon_eta, gen_lowKaon_eta, genThirdTrack_eta;

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
        FiveTracksCollection_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("FiveTracksCand"))),
        TriggerResults_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
        thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PrimaryVertex"))),
	      IsMC_(iConfig.getParameter<bool>("IsMC")),
        OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
        HLTs_(iConfig.getParameter<std::vector<std::string>>("HLTs")),
        HLTFilters_(iConfig.getParameter<std::vector<std::string>>("Filters")),
        TreeName_(iConfig.getParameter<std::string>("TreeName"))
{
	      edm::Service<TFileService> fs;
        fivetracks_tree = fs->make<TTree>(TreeName_.data(),"Tree of DiMuon and Three Tracks");

        fivetracks_tree->Branch("run",                &run,                "run/I");
        fivetracks_tree->Branch("event",              &event,              "event/I");
        fivetracks_tree->Branch("lumi",              &lumi,              "lumi/I");
        fivetracks_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");
        fivetracks_tree->Branch("trigger",            &trigger,            "trigger/I");

        if(!OnlyGen_)
        {
        fivetracks_tree->Branch("noFiveCandidates",      &noFiveCandidates,      "noFiveCandidates/I");

        fivetracks_tree->Branch("dimuonditrk_id",      &dimuonditrk_id,      "dimuonditrk_id/I");
        fivetracks_tree->Branch("dimuon_id",      &dimuon_id,      "dimuon_id/I");
        fivetracks_tree->Branch("p_id",      &p_id,      "p_id/I");
        fivetracks_tree->Branch("m_id",      &m_id,      "m_id/I");
        fivetracks_tree->Branch("t_id",      &t_id,      "m_id/I");

        fivetracks_tree->Branch("five_p4",     "TLorentzVector", &five_p4);
        fivetracks_tree->Branch("dimuonditrk_p4",     "TLorentzVector", &dimuonditrk_p4);
        fivetracks_tree->Branch("ditrack_p4",      "TLorentzVector", &ditrack_p4);
        fivetracks_tree->Branch("dimuon_p4",      "TLorentzVector", &dimuon_p4);

        fivetracks_tree->Branch("lowMuon_p4",    "TLorentzVector", &lowMuon_p4);
        fivetracks_tree->Branch("highMuon_p4",   "TLorentzVector", &highMuon_p4);

        fivetracks_tree->Branch("highKaon_p4",   "TLorentzVector", &highKaon_p4);
        fivetracks_tree->Branch("lowKaon_p4",    "TLorentzVector", &lowKaon_p4);
        fivetracks_tree->Branch("thirdKaon_p4",  "TLorentzVector", &thirdKaon_p4);

        fivetracks_tree->Branch("highPion_p4",   "TLorentzVector", &highPion_p4);
        fivetracks_tree->Branch("lowPion_p4",    "TLorentzVector", &lowPion_p4);
        fivetracks_tree->Branch("thirdPion_p4",  "TLorentzVector", &thirdPion_p4);

        fivetracks_tree->Branch("highProton_p4",   "TLorentzVector", &highProton_p4);
        fivetracks_tree->Branch("lowProton_p4",    "TLorentzVector", &lowProton_p4);
        fivetracks_tree->Branch("thirdProton_p4",  "TLorentzVector", &thirdProton_p4);

        fivetracks_tree->Branch("dimuonditrk_m",       &dimuonditrk_m,        "dimuonditrk_m/D"); // Original DiMuonDiTrack (J/Psi Phi)
        fivetracks_tree->Branch("dimuonditrk_pt",      &dimuonditrk_pt,       "dimuonditrk_pt/D");
        fivetracks_tree->Branch("dimuonditrk_eta",     &dimuonditrk_eta,      "dimuonditrk_eta/D");
        fivetracks_tree->Branch("dimuonditrk_phi",     &dimuonditrk_phi,      "dimuonditrk_phi/D");
        fivetracks_tree->Branch("dimuonditrk_p",       &dimuonditrk_p,      "dimuonditrk_p/D");

        fivetracks_tree->Branch("dimuon_m",      &dimuon_m,     "dimuon_m/D");
        fivetracks_tree->Branch("dimuon_pt",     &dimuon_pt,    "dimuon_pt/D");
        fivetracks_tree->Branch("dimuon_eta",    &dimuon_pt,    "dimuon_eta/D");
        fivetracks_tree->Branch("dimuon_p",      &dimuon_p,     "dimuon_p/D");
        fivetracks_tree->Branch("dimuon_phi",    &dimuon_phi,   "dimuon_phi/D");

        fivetracks_tree->Branch("highTrackMatch",     &highTrackMatch,    "highTrackMatch/D");
        fivetracks_tree->Branch("lowTrackMatch",      &lowTrackMatch,     "lowTrackMatch/D");
        fivetracks_tree->Branch("lowMuonMatch",       &lowMuonMatch,      "lowMuonMatch/D");
        fivetracks_tree->Branch("thirdTrackMatch",    &thirdTrackMatch,   "thirdTrackMatch/D");
        //Di Tracks
        fivetracks_tree->Branch("ditrack_m",      &ditrack_m,     "ditrack_m/D"); //the original ditrack (supposed to be Phi->KK)

        fivetracks_tree->Branch("ditrackOne_pt",     &ditrackOne_pt,    "ditrackOne_pt/D");
        fivetracks_tree->Branch("ditrackOne_eta",    &ditrackOne_eta,    "ditrackOne_eta/D");
        fivetracks_tree->Branch("ditrackOne_phi",    &ditrackOne_phi,   "ditrackOne_phi/D");
        fivetracks_tree->Branch("ditrackOne_p",      &ditrackOne_p,   "ditrackOne_p/D");

        fivetracks_tree->Branch("ditrackTwo_pt",     &ditrackTwo_pt,    "ditrackTwo_pt/D");
        fivetracks_tree->Branch("ditrackTwo_eta",    &ditrackTwo_eta,    "ditrackTwo_eta/D");
        fivetracks_tree->Branch("ditrackTwo_phi",    &ditrackTwo_phi,   "ditrackTwo_phi/D");
        fivetracks_tree->Branch("ditrackTwo_p",      &ditrackTwo_p,   "ditrackTwo_p/D");

        fivetracks_tree->Branch("ditrackThree_pt",     &ditrackThree_pt,    "ditrackThree_pt/D");
        fivetracks_tree->Branch("ditrackThree_eta",    &ditrackThree_eta,    "ditrackThree_eta/D");
        fivetracks_tree->Branch("ditrackThree_phi",    &ditrackThree_phi,   "ditrackThree_phi/D");
        fivetracks_tree->Branch("ditrackThree_p",      &ditrackThree_p,   "ditrackThree_p/D");

        //The kinematic doesn't change, only mass
        fivetracks_tree->Branch("highMuon_pt",         &highMuon_pt,         "highMuon_pt/D");
        fivetracks_tree->Branch("highMuon_eta",        &highMuon_eta,        "highMuon_eta/D");
        fivetracks_tree->Branch("highMuon_phi",        &highMuon_phi,        "highMuon_phi/D");
        fivetracks_tree->Branch("highMuon_charge",     &highMuon_charge,     "highMuon_charge/D");
        fivetracks_tree->Branch("highMuon_dz",         &highMuon_dz,         "highMuon_dz/D");
        fivetracks_tree->Branch("highMuon_dxy",        &highMuon_dxy,        "highMuon_dxy/D");

        fivetracks_tree->Branch("lowMuon_pt",         &lowMuon_pt,         "lowMuon_pt/D");
        fivetracks_tree->Branch("lowMuon_eta",        &lowMuon_eta,        "lowMuon_eta/D");
        fivetracks_tree->Branch("lowMuon_phi",        &lowMuon_phi,        "lowMuon_phi/D");
        fivetracks_tree->Branch("lowMuon_charge",     &lowMuon_charge,     "lowMuon_charge/D");
        fivetracks_tree->Branch("lowMuon_dz",         &lowMuon_dz,         "lowMuon_dz/D");
        fivetracks_tree->Branch("lowMuon_dxy",        &lowMuon_dxy,        "lowMuon_dxy/D");

        //The kinematic doesn't change, only mass
        fivetracks_tree->Branch("highTrack_pt",         &highTrack_pt,         "highTrack_pt/D");
        fivetracks_tree->Branch("highTrack_eta",        &highTrack_eta,        "highTrack_eta/D");
        fivetracks_tree->Branch("highTrack_phi",        &highTrack_phi,        "highTrack_phi/D");
        fivetracks_tree->Branch("highTrack_charge",     &highTrack_charge,     "highTrack_charge/D");
        fivetracks_tree->Branch("highTrack_dz",         &highTrack_dz,         "highTrack_dz/D");
        fivetracks_tree->Branch("highTrack_dxy",        &highTrack_dxy,        "highTrack_dxy/D");

        fivetracks_tree->Branch("lowTrack_pt",           &lowTrack_pt,        "lowTrack_pt/D");
        fivetracks_tree->Branch("lowTrack_eta",          &lowTrack_eta,       "lowTrack_eta/D");
        fivetracks_tree->Branch("lowTrack_phi",          &lowTrack_phi,       "lowTrack_phi/D");
        fivetracks_tree->Branch("lowTrack_charge",       &lowTrack_charge,    "lowTrack_charge/D");
        fivetracks_tree->Branch("lowTrack_dz",           &lowTrack_dz,          "lowTrack_dz/D");
        fivetracks_tree->Branch("lowTrack_dxy",          &lowTrack_dxy,        "lowTrack_dxy/D");

        fivetracks_tree->Branch("thirdTrack_pt",         &thirdTrack_pt,        "thirdTrack_pt/D");
        fivetracks_tree->Branch("thirdTrack_eta",        &thirdTrack_eta,        "thirdTrack_eta/D");
        fivetracks_tree->Branch("thirdTrack_phi",        &thirdTrack_phi,        "thirdTrack_phi/D");
        fivetracks_tree->Branch("thirdTrack_charge",     &thirdTrack_charge,  "thirdTrack_charge/D");
        fivetracks_tree->Branch("thirdTrack_dz",         &thirdTrack_dz,          "thirdTrack_dz/D");
        fivetracks_tree->Branch("thirdTrack_dxy",        &thirdTrack_dxy,        "thirdTrack_dxy/D");

        //J/Psi TrTr system
        fivetracks_tree->Branch("dimuonDiTrkOne_pt",         &dimuonDiTrkOne_pt,          "dimuonDiTrkOne_pt/D");
        fivetracks_tree->Branch("dimuonDiTrkOne_eta",        &dimuonDiTrkOne_eta,        "dimuonDiTrkOne_eta/D");
        fivetracks_tree->Branch("dimuonDiTrkOne_phi",        &dimuonDiTrkOne_phi,        "dimuonDiTrkOne_phi/D");
        fivetracks_tree->Branch("dimuonDiTrkOne_charge",     &dimuonDiTrkOne_charge,        "dimuonDiTrkOne_charge/D");
        fivetracks_tree->Branch("dimuonDiTrkOne_p",          &dimuonDiTrkOne_p,          "dimuonDiTrkOne_p/D");

        fivetracks_tree->Branch("dimuonDiTrkTwo_pt",         &dimuonDiTrkTwo_pt,         "dimuonDiTrkTwo_pt/D");
        fivetracks_tree->Branch("dimuonDiTrkTwo_eta",        &dimuonDiTrkTwo_eta,        "dimuonDiTrkTwo_eta/D");
        fivetracks_tree->Branch("dimuonDiTrkTwo_phi",        &dimuonDiTrkTwo_phi,        "dimuonDiTrkTwo_phi/D");
        fivetracks_tree->Branch("dimuonDiTrkTwo_charge",     &dimuonDiTrkTwo_charge,     "dimuonDiTrkTwo_charge/D");
        fivetracks_tree->Branch("dimuonDiTrkTwo_p",         &dimuonDiTrkTwo_p,         "dimuonDiTrkTwo_p/D");

        fivetracks_tree->Branch("dimuonDiTrkThree_pt",       &dimuonDiTrkThree_pt,         "dimuonDiTrkThree_pt/D");
        fivetracks_tree->Branch("dimuonDiTrkThree_eta",      &dimuonDiTrkThree_eta,        "dimuonDiTrkThree_eta/D");
        fivetracks_tree->Branch("dimuonDiTrkThree_phi",      &dimuonDiTrkThree_phi,        "dimuonDiTrkThree_phi/D");
        fivetracks_tree->Branch("dimuonDiTrkThree_charge",   &dimuonDiTrkThree_charge,     "dimuonDiTrkThree_charge/D");
        fivetracks_tree->Branch("dimuonDiTrkThree_p",        &dimuonDiTrkThree_p,          "dimuonDiTrkThree_p/D");

        fivetracks_tree->Branch("dimuon_vProb",        &dimuon_vProb,        "dimuon_vProb/D");
        fivetracks_tree->Branch("dimuon_vChi2",       &dimuon_vChi2,        "dimuon_vChi2/D");
        fivetracks_tree->Branch("dimuon_DCA",          &dimuon_DCA,          "dimuon_DCA/D");
        fivetracks_tree->Branch("dimuon_ctauPV",       &dimuon_ctauPV,       "dimuon_ctauPV/D");
        fivetracks_tree->Branch("dimuon_ctauErrPV",    &dimuon_ctauErrPV,    "dimuon_ctauErrPV/D");
        fivetracks_tree->Branch("dimuon_cosAlpha",     &dimuon_cosAlpha,     "dimuon_cosAlpha/D");
        // fivetracks_tree->Branch("dimuon_triggerMatch", &dimuon_triggerMatch, "dimuon_triggerMatch/I");

        //TriTrack system
        fivetracks_tree->Branch("triTrack_m",         &triTrack_m,         "triTrack_m/D");
        fivetracks_tree->Branch("triTrack_pt",         &triTrack_pt,         "triTrack_pt/D");
        fivetracks_tree->Branch("triTrack_eta",        &triTrack_eta,        "triTrack_eta/D");
        fivetracks_tree->Branch("triTrack_phi",        &triTrack_phi,        "triTrack_phi/D");
        fivetracks_tree->Branch("triTrack_charge",     &triTrack_charge,     "triTrack_charge/D");

        fivetracks_tree->Branch("dimuonditrk_vProb",      &dimuonditrk_vProb,        "dimuonditrk_vProb/D");
        fivetracks_tree->Branch("dimuonditrk_vChi2",      &dimuonditrk_vChi2,        "dimuonditrk_vChi2/D");
        fivetracks_tree->Branch("dimuonditrk_nDof",       &dimuonditrk_nDof,         "dimuonditrk_nDof/D");
        fivetracks_tree->Branch("dimuonditrk_charge",     &dimuonditrk_charge,       "dimuonditrk_charge/I");

        fivetracks_tree->Branch("dimuonditrk_cosAlpha",      &dimuonditrk_cosAlpha,        "dimuonditrk_cosAlpha/D");
        fivetracks_tree->Branch("dimuonditrk_ctauPV",      &dimuonditrk_ctauPV,        "dimuonditrk_ctauPV/D");
        fivetracks_tree->Branch("dimuonditrk_ctauErrPV",      &dimuonditrk_ctauErrPV,        "dimuonditrk_ctauErrPV/D");

        fivetracks_tree->Branch("dimuonditrk_cosAlphaCA",      &dimuonditrk_cosAlphaCA,        "dimuonditrk_cosAlphaCA/D");
        fivetracks_tree->Branch("dimuonditrk_ctauPVCA",      &dimuonditrk_ctauPVCA,        "dimuonditrk_ctauPVCA/D");
        fivetracks_tree->Branch("dimuonditrk_ctauErrPVCA",      &dimuonditrk_ctauErrPVCA,        "dimuonditrk_ctauErrPVCA/D");

        fivetracks_tree->Branch("dimuonditrk_cosAlphaDZ",      &dimuonditrk_cosAlphaDZ,        "dimuonditrk_cosAlphaDZ/D");
        fivetracks_tree->Branch("dimuonditrk_ctauPVDZ",      &dimuonditrk_ctauPVDZ,        "dimuonditrk_ctauPVDZ/D");
        fivetracks_tree->Branch("dimuonditrk_ctauErrPVDZ",      &dimuonditrk_ctauErrPVDZ,        "dimuonditrk_ctauErrPVDZ/D");

        fivetracks_tree->Branch("dimuonditrk_cosAlphaBS",      &dimuonditrk_cosAlphaBS,        "dimuonditrk_cosAlphaBS/D");
        fivetracks_tree->Branch("dimuonditrk_ctauPVBS",      &dimuonditrk_ctauPVBS,        "dimuonditrk_ctauPVBS/D");
        fivetracks_tree->Branch("dimuonditrk_ctauErrPVBS",      &dimuonditrk_ctauErrPVBS,        "dimuonditrk_ctauErrPVBS/D");

        fivetracks_tree->Branch("dimuonditrk_vx",     &dimuonditrk_vx,      "dimuonditrk_vx/D");
        fivetracks_tree->Branch("dimuonditrk_vy",     &dimuonditrk_vy,      "dimuonditrk_vy/D");
        fivetracks_tree->Branch("dimuonditrk_vz",     &dimuonditrk_vz,      "dimuonditrk_vz/D");

        fivetracks_tree->Branch("dca_m1m2",      &dca_m1m2,        "dca_m1m2/D");
        fivetracks_tree->Branch("dca_m1t1",      &dca_m1t1,        "dca_m1t1/D");
        fivetracks_tree->Branch("dca_m1t2",      &dca_m1t2,        "dca_m1t2/D");
        fivetracks_tree->Branch("dca_m2t1",      &dca_m2t1,        "dca_m2t1/D");
        fivetracks_tree->Branch("dca_m2t2",      &dca_m2t2,        "dca_m2t2/D");
        fivetracks_tree->Branch("dca_t1t2",      &dca_t1t2,        "dca_t1t2/D");
        fivetracks_tree->Branch("dca_m1t3",      &dca_m1t3,        "dca_m1t3/D");
        fivetracks_tree->Branch("dca_m2t3",      &dca_m2t3,        "dca_m2t3/D");
        fivetracks_tree->Branch("dca_t1t3",      &dca_t1t3,        "dca_t1t3/D");
        fivetracks_tree->Branch("dca_t2t3",      &dca_t2t3,        "dca_t2t3/D");

        fivetracks_tree->Branch("highTrackMuonDR",        &highTrackMuonDR,        "highTrackMuonDR/D");
        fivetracks_tree->Branch("highTrackMuonDP",        &highTrackMuonDP,        "highTrackMuonDP/D");
        fivetracks_tree->Branch("highTrackMuonDPt",        &highTrackMuonDPt,        "highTrackMuonDPt/D");

        fivetracks_tree->Branch("lowTrackMuonDR",        &lowTrackMuonDR,        "lowTrackMuonDR/D");
        fivetracks_tree->Branch("lowTrackMuonDP",        &lowTrackMuonDP,        "lowTrackMuonDP/D");
        fivetracks_tree->Branch("lowTrackMuonDPt",        &lowTrackMuonDPt,        "lowTrackMuonDPt/D");

        fivetracks_tree->Branch("thirdTrackMuonDR",        &thirdTrackMuonDR,        "thirdTrackMuonDR/D");
        fivetracks_tree->Branch("thirdTrackMuonDP",        &thirdTrackMuonDP,        "thirdTrackMuonDP/D");
        fivetracks_tree->Branch("thirdTrackMuonDPt",        &thirdTrackMuonDPt,        "thirdTrackMuonDPt/D");

        fivetracks_tree->Branch("tPFromPV",        &tPFromPV,        "tPFromPV/D");
        fivetracks_tree->Branch("tMFromPV",        &tMFromPV,        "tMFromPV/D");
        fivetracks_tree->Branch("tTFromPV",        &tTFromPV,        "tMFTomPV/D");
        fivetracks_tree->Branch("tPFromPVCA",        &tPFromPVCA,        "tPFromPVCA/D");
        fivetracks_tree->Branch("tMFromPVCA",        &tMFromPVCA,        "tMFromPVCA/D");
        fivetracks_tree->Branch("tTFromPVCA",        &tTFromPVCA,        "tMFTomPVCA/D");
        fivetracks_tree->Branch("tPFromPVDZ",      &tPFromPVDZ,        "tPFromPVDZ/D");
        fivetracks_tree->Branch("tMFromPVDZ",      &tMFromPVDZ,        "tMFromPVDZ/D");
        fivetracks_tree->Branch("tTFromPVDZ",      &tTFromPVDZ,        "tTFromPVDZ/D");
        // fivetracks_tree->Branch("tPFromPVBS",      &tPFromPVBS,        "tPFromPVBS/D");
        // fivetracks_tree->Branch("tMFromPVBS",      &tMFromPVBS,        "tMFromPVBS/D");
        // fivetracks_tree->Branch("tTFromPVBS",      &tTFromPVBS,        "tTFromPVBS/D");

        fivetracks_tree->Branch("five_m",             &five_m,               "five_m/D");
        fivetracks_tree->Branch("five_m_ref",         &five_m_ref,           "five_m_ref/D");
        fivetracks_tree->Branch("five_mass_ppk",      &five_mass_ppk,        "five_mass_ppk/D");
        fivetracks_tree->Branch("five_mass_kpp",      &five_mass_kpp,        "five_mass_kpp/D");
        fivetracks_tree->Branch("five_mass_pkp",      &five_mass_pkp,        "five_mass_pkp/D");
        fivetracks_tree->Branch("five_mass_ppp",      &five_mass_ppp,        "five_mass_ppp/D");

        fivetracks_tree->Branch("five_pt",      &five_pt,       "five_pt/D");
        fivetracks_tree->Branch("five_eta",     &five_eta,      "five_eta/D");
        fivetracks_tree->Branch("five_phi",     &five_phi,      "five_phi/D");
        fivetracks_tree->Branch("five_p",       &five_p,        "five_p/D");

        fivetracks_tree->Branch("five_cosAlpha",     &five_cosAlpha,        "five_cosAlpha/D");
        fivetracks_tree->Branch("five_ctauPV",       &five_ctauPV,        "five_ctauPV/D");
        fivetracks_tree->Branch("five_ctauErrPV",    &five_ctauErrPV,        "five_ctauErrPV/D");

        fivetracks_tree->Branch("five_cosAlphaCA",    &five_cosAlphaCA,        "five_cosAlphaCA/D");
        fivetracks_tree->Branch("five_ctauPVCA",      &five_ctauPVCA,        "five_ctauPVCA/D");
        fivetracks_tree->Branch("five_ctauErrPVCA",   &five_ctauErrPVCA,        "five_ctauErrPVCA/D");

        fivetracks_tree->Branch("five_cosAlphaDZ",    &five_cosAlphaDZ,        "five_cosAlphaDZ/D");
        fivetracks_tree->Branch("five_ctauPVDZ",      &five_ctauPVDZ,        "five_ctauPVDZ/D");
        fivetracks_tree->Branch("five_ctauErrPVDZ",   &five_ctauErrPVDZ,        "five_ctauErrPVDZ/D");

        fivetracks_tree->Branch("five_cosAlphaBS",    &five_cosAlphaBS,        "five_cosAlphaBS/D");
        fivetracks_tree->Branch("five_ctauPVBS",      &five_ctauPVBS,        "five_ctauPVBS/D");
        fivetracks_tree->Branch("five_ctauErrPVBS",   &five_ctauErrPVBS,        "five_ctauErrPVBS/D");

        fivetracks_tree->Branch("five_vProb",       &five_vProb,      "five_vProb/D");
        fivetracks_tree->Branch("five_nDof",        &five_nDof,       "five_nDof/D");
        fivetracks_tree->Branch("five_vChi2",       &five_vChi2,      "five_vChi2/D");

        fivetracks_tree->Branch("five_vx",     &five_vx,      "five_vx/D");
        fivetracks_tree->Branch("five_vy",     &five_vy,      "five_vy/D");
        fivetracks_tree->Branch("five_vz",     &five_vz,      "five_vz/D");

        fivetracks_tree->Branch("bestPV_X",     &bestPV_X,      "bestPV_X/D");
        fivetracks_tree->Branch("bestPV_Y",     &bestPV_Y,      "bestPV_Y/D");
        fivetracks_tree->Branch("bestPV_Z",     &bestPV_Z,      "bestPV_Z/D");

        fivetracks_tree->Branch("cosAlphaPV_X",     &cosAlphaPV_X,      "cosAlphaPV_X/D");
        fivetracks_tree->Branch("cosAlphaPV_Y",     &cosAlphaPV_Y,      "cosAlphaPV_Y/D");
        fivetracks_tree->Branch("cosAlphaPV_Z",     &cosAlphaPV_Z,      "cosAlphaPV_Z/D");

        fivetracks_tree->Branch("bS_X",     &bS_X,      "bS_X/D");
        fivetracks_tree->Branch("bS_Y",     &bS_Y,      "bS_Y/D");
        fivetracks_tree->Branch("bS_Z",     &bS_Z,      "bS_Z/D");

        fivetracks_tree->Branch("zPV_X",     &zPV_X,      "zPV_X/D");
        fivetracks_tree->Branch("zPV_Y",     &zPV_Y,      "zPV_Y/D");
        fivetracks_tree->Branch("zPV_Z",     &zPV_Z,      "zPV_Z/D");

        fivetracks_tree->Branch("five_charge",     &five_charge,       "five_charge/I");

        //Muon flags
        fivetracks_tree->Branch("lowMuon_isTight",        &lowMuon_isTight,        "lowMuon_isTight/O");
        fivetracks_tree->Branch("lowMuon_isLoose",        &lowMuon_isLoose,        "lowMuon_isLoose/O");
        fivetracks_tree->Branch("lowMuon_isSoft",        &lowMuon_isSoft,        "lowMuon_isSoft/O");
        fivetracks_tree->Branch("lowMuon_isMedium",        &lowMuon_isMedium,        "lowMuon_isMedium/O");
        fivetracks_tree->Branch("lowMuon_isHighPt",        &lowMuon_isHighPt,        "lowMuon_isHighPt/O");

        fivetracks_tree->Branch("lowMuon_isTracker",        &lowMuon_isTracker,        "lowMuon_isTracker/O");
        fivetracks_tree->Branch("lowMuon_isGlobal",        &lowMuon_isGlobal,        "lowMuon_isGlobal/O");

        fivetracks_tree->Branch("lowMuon_NPixelHits",        &lowMuon_NPixelHits,        "lowMuon_NPixelHits/I");
        fivetracks_tree->Branch("lowMuon_NStripHits",        &lowMuon_NStripHits,        "lowMuon_NStripHits/I");
        fivetracks_tree->Branch("lowMuon_NTrackhits",        &lowMuon_NTrackhits,        "lowMuon_NTrackhits/I");
        fivetracks_tree->Branch("lowMuon_NBPixHits",        &lowMuon_NBPixHits,        "lowMuon_NBPixHits/I");

        fivetracks_tree->Branch("lowMuon_NPixLayers",        &lowMuon_NPixLayers,        "lowMuon_NPixLayers/I");
        fivetracks_tree->Branch("lowMuon_NTraLayers",        &lowMuon_NTraLayers,        "lowMuon_NTraLayers/I");
        fivetracks_tree->Branch("lowMuon_NStrLayers",        &lowMuon_NStrLayers,        "lowMuon_NStrLayers/I");
        fivetracks_tree->Branch("lowMuon_NBPixLayers",        &lowMuon_NBPixLayers,        "lowMuon_NBPixLayers/I");

        fivetracks_tree->Branch("highMuon_isTight",        &highMuon_isTight,        "highMuon_isTight/O");
        fivetracks_tree->Branch("highMuon_isLoose",        &highMuon_isLoose,        "highMuon_isLoose/O");
        fivetracks_tree->Branch("highMuon_isSoft",        &highMuon_isSoft,        "highMuon_isSoft/O");
        fivetracks_tree->Branch("highMuon_isMedium",        &highMuon_isMedium,        "highMuon_isMedium/O");
        fivetracks_tree->Branch("highMuon_isHighPt",        &highMuon_isHighPt,        "highMuon_isHighPt/O");

        fivetracks_tree->Branch("highMuon_isTracker",        &highMuon_isTracker,        "highMuon_isTracker/O");
        fivetracks_tree->Branch("highMuon_isGlobal",        &highMuon_isGlobal,        "highMuon_isGlobal/O");

        fivetracks_tree->Branch("highMuon_NPixelHits",        &highMuon_NPixelHits,        "highMuon_NPixelHits/I");
        fivetracks_tree->Branch("highMuon_NStripHits",        &highMuon_NStripHits,        "highMuon_NStripHits/I");
        fivetracks_tree->Branch("highMuon_NTrackhits",        &highMuon_NTrackhits,        "highMuon_NTrackhits/I");
        fivetracks_tree->Branch("highMuon_NBPixHits",        &highMuon_NBPixHits,        "highMuon_NBPixHits/I");

        fivetracks_tree->Branch("highMuon_NPixLayers",        &highMuon_NPixLayers,        "highMuon_NPixLayers/I");
        fivetracks_tree->Branch("highMuon_NTraLayers",        &highMuon_NTraLayers,        "highMuon_NTraLayers/I");
        fivetracks_tree->Branch("highMuon_NStrLayers",        &highMuon_NStrLayers,        "highMuon_NStrLayers/I");
        fivetracks_tree->Branch("highMuon_NBPixLayers",        &highMuon_NBPixLayers,        "highMuon_NBPixLayers/I");

        fivetracks_tree->Branch("lowMuon_type",     &lowMuon_type,       "lowMuon_type/i");
        fivetracks_tree->Branch("highMuon_type",     &highMuon_type,       "highMuon_type/i");


        //Tracks Flags

        fivetracks_tree->Branch("highTrack_NPixelHits",        &highTrack_NPixelHits,        "highTrack_NPixelHits/I");
        fivetracks_tree->Branch("highTrack_NStripHits",        &highTrack_NStripHits,        "highTrack_NStripHits/I");
        fivetracks_tree->Branch("highTrack_NTrackhits",        &highTrack_NTrackhits,        "highTrack_NTrackhits/I");
        fivetracks_tree->Branch("highTrack_NBPixHits",         &highTrack_NBPixHits,        "highTrack_NBPixHits/I");

        fivetracks_tree->Branch("highTrack_NPixLayers",        &highTrack_NPixLayers,        "highTrack_NPixLayers/I");
        fivetracks_tree->Branch("highTrack_NTraLayers",        &highTrack_NTraLayers,        "highTrack_NTraLayers/I");
        fivetracks_tree->Branch("highTrack_NStrLayers",        &highTrack_NStrLayers,        "highTrack_NStrLayers/I");
        fivetracks_tree->Branch("highTrack_NBPixLayers",       &highTrack_NBPixLayers,        "highTrack_NBPixLayers/I");

        fivetracks_tree->Branch("lowTrack_NPixelHits",        &lowTrack_NPixelHits,        "lowTrack_NPixelHits/I");
        fivetracks_tree->Branch("lowTrack_NStripHits",        &lowTrack_NStripHits,        "lowTrack_NStripHits/I");
        fivetracks_tree->Branch("lowTrack_NTrackhits",        &lowTrack_NTrackhits,        "lowTrack_NTrackhits/I");
        fivetracks_tree->Branch("lowTrack_NBPixHits",         &lowTrack_NBPixHits,        "lowTrack_NBPixHits/I");

        fivetracks_tree->Branch("lowTrack_NPixLayers",        &lowTrack_NPixLayers,        "lowTrack_NPixLayers/I");
        fivetracks_tree->Branch("lowTrack_NTraLayers",        &lowTrack_NTraLayers,        "lowTrack_NTraLayers/I");
        fivetracks_tree->Branch("lowTrack_NStrLayers",        &lowTrack_NStrLayers,        "lowTrack_NStrLayers/I");
        fivetracks_tree->Branch("lowTrack_NBPixLayers",       &lowTrack_NBPixLayers,        "lowTrack_NBPixLayers/I");

        fivetracks_tree->Branch("thirdTrack_NPixelHits",        &thirdTrack_NPixelHits,        "thirdTrack_NPixelHits/I");
        fivetracks_tree->Branch("thirdTrack_NStripHits",        &thirdTrack_NStripHits,        "thirdTrack_NStripHits/I");
        fivetracks_tree->Branch("thirdTrack_NTrackhits",        &thirdTrack_NTrackhits,        "thirdTrack_NTrackhits/I");
        fivetracks_tree->Branch("thirdTrack_NBPixHits",         &thirdTrack_NBPixHits,        "thirdTrack_NBPixHits/I");

        fivetracks_tree->Branch("thirdTrack_NPixLayers",        &thirdTrack_NPixLayers,        "thirdTrack_NPixLayers/I");
        fivetracks_tree->Branch("thirdTrack_NTraLayers",        &thirdTrack_NTraLayers,        "thirdTrack_NTraLayers/I");
        fivetracks_tree->Branch("thirdTrack_NStrLayers",        &thirdTrack_NStrLayers,        "thirdTrack_NStrLayers/I");
        fivetracks_tree->Branch("thirdTrack_NBPixLayers",       &thirdTrack_NBPixLayers,        "thirdTrack_NBPixLayers/I");



      }

      if (IsMC_ || OnlyGen_) {

        fivetracks_tree->Branch("gen_five_p4", "TLorentzVector",  &gen_five_p4);
        fivetracks_tree->Branch("gen_dimuonditrk_p4", "TLorentzVector",  &gen_dimuonditrk_p4);
        fivetracks_tree->Branch("gen_jpsi_p4", "TLorentzVector",  &gen_jpsi_p4);
        fivetracks_tree->Branch("gen_phi_p4", "TLorentzVector",  &gen_phi_p4);

        fivetracks_tree->Branch("gen_highKaon_p4", "TLorentzVector",  &gen_highKaon_p4);
        fivetracks_tree->Branch("gen_lowMuon_p4",  "TLorentzVector",  &gen_lowMuon_p4);
        fivetracks_tree->Branch("gen_highMuon_p4",  "TLorentzVector",  &gen_highMuon_p4);
        fivetracks_tree->Branch("gen_lowKaon_p4",  "TLorentzVector",  &gen_lowKaon_p4);
        fivetracks_tree->Branch("genThirdTrack_p4",  "TLorentzVector",  &genThirdTrack_p4);

        fivetracks_tree->Branch("gen_five_pdg",&gen_five_pdg,"gen_five_pdg/D");
        fivetracks_tree->Branch("gen_dimuonditrk_pdg",&gen_dimuonditrk_pdg,"gen_dimuonditrk_pdg/D");
        fivetracks_tree->Branch("gen_phi_pdg",&gen_phi_pdg,"gen_phi_pdg/D");
        fivetracks_tree->Branch("gen_jpsi_pdg",&gen_jpsi_pdg,"gen_jpsi_pdg/D");

        fivetracks_tree->Branch("gen_lowMuon_pdg",&gen_lowMuon_pdg,"gen_lowMuon_pdg/D");
        fivetracks_tree->Branch("gen_highMuon_pdg",&gen_highMuon_pdg,"gen_highMuon_pdg/D");
        fivetracks_tree->Branch("gen_highKaon_pdg",&gen_highKaon_pdg,"gen_highKaon_pdg/D");
        fivetracks_tree->Branch("gen_lowKaon_pdg",&gen_lowKaon_pdg,"gen_lowKaon_pdg/D");
        fivetracks_tree->Branch("genThirdTrack_pdg",&genThirdTrack_pdg,"genThirdTrack_pdg/D");

        fivetracks_tree->Branch("gen_lowMuon_mompdg",&gen_lowMuon_mompdg,"gen_lowMuon_mompdg/D");
        fivetracks_tree->Branch("gen_highMuon_mompdg",&gen_highMuon_mompdg,"gen_highMuon_mompdg/D");
        fivetracks_tree->Branch("gen_highKaon_mompdg",&gen_highKaon_mompdg,"gen_highKaon_mompdg/D");
        fivetracks_tree->Branch("gen_lowKaon_mompdg",&gen_lowKaon_mompdg,"gen_lowKaon_mompdg/D");
        fivetracks_tree->Branch("genThirdTrack_mompdg",&genThirdTrack_mompdg,"genThirdTrack_mompdg/D");

        fivetracks_tree->Branch("gen_lowMuon_status",&gen_lowMuon_status,"gen_lowMuon_status/D");
        fivetracks_tree->Branch("gen_highMuon_status",&gen_highMuon_status,"gen_highMuon_status/D");
        fivetracks_tree->Branch("gen_highKaon_status",&gen_highKaon_status,"gen_highKaon_status/D");
        fivetracks_tree->Branch("gen_lowKaon_status",&gen_lowKaon_status,"gen_lowKaon_status/D");
        fivetracks_tree->Branch("genThirdTrack_status",&genThirdTrack_status,"genThirdTrack_status/D");

        fivetracks_tree->Branch("gen_lowMuon_p",&gen_lowMuon_p,"gen_lowMuon_p/D");
        fivetracks_tree->Branch("gen_highMuon_p",&gen_highMuon_p,"gen_highMuon_p/D");
        fivetracks_tree->Branch("gen_highKaon_p",&gen_highKaon_p,"gen_highKaon_p/D");
        fivetracks_tree->Branch("gen_lowKaon_p",&gen_lowKaon_p,"gen_lowKaon_p/D");
        fivetracks_tree->Branch("genThirdTrack_p",&genThirdTrack_p,"genThirdTrack_p/D");

        fivetracks_tree->Branch("gen_lowMuon_pt",&gen_lowMuon_pt,"gen_lowMuon_pt/D");
        fivetracks_tree->Branch("gen_highMuon_pt",&gen_highMuon_pt,"gen_highMuon_pt/D");
        fivetracks_tree->Branch("gen_highKaon_pt",&gen_highKaon_pt,"gen_highKaon_pt/D");
        fivetracks_tree->Branch("gen_lowKaon_pt",&gen_lowKaon_pt,"gen_lowKaon_pt/D");
        fivetracks_tree->Branch("genThirdTrack_pt",&genThirdTrack_pt,"genThirdTrack_pt/D");

        fivetracks_tree->Branch("gen_lowMuon_eta",&gen_lowMuon_eta,"gen_lowMuon_eta/D");
        fivetracks_tree->Branch("gen_highMuon_eta",&gen_highMuon_eta,"gen_highMuon_eta/D");
        fivetracks_tree->Branch("gen_highKaon_eta",&gen_highKaon_eta,"gen_highKaon_eta/D");
        fivetracks_tree->Branch("gen_lowKaon_eta",&gen_lowKaon_eta,"gen_lowKaon_eta/D");
        fivetracks_tree->Branch("genThirdTrack_eta",&genThirdTrack_eta,"genThirdTrack_eta/D");

        fivetracks_tree->Branch("gen_lowMuon_phi",&gen_lowMuon_phi,"gen_lowMuon_phi/D");
        fivetracks_tree->Branch("gen_highMuon_phi",&gen_highMuon_phi,"gen_highMuon_phi/D");
        fivetracks_tree->Branch("gen_highKaon_phi",&gen_highKaon_phi,"gen_highKaon_phi/D");
        fivetracks_tree->Branch("gen_lowKaon_phi",&gen_lowKaon_phi,"gen_lowKaon_phi/D");
        fivetracks_tree->Branch("genThirdTrack_phi",&genThirdTrack_phi,"genThirdTrack_phi/D");

        fivetracks_tree->Branch("gen_five_prompt",&gen_five_prompt,"gen_five_prompt/D");
        fivetracks_tree->Branch("gen_five_pt",&gen_five_pt,"gen_five_pt/D");
        fivetracks_tree->Branch("gen_five_p",&gen_five_p,"gen_five_p/D");
        fivetracks_tree->Branch("gen_five_eta",&gen_five_eta,"gen_five_eta/D");
        fivetracks_tree->Branch("gen_five_phi",&gen_five_phi,"gen_five_phi/D");

        fivetracks_tree->Branch("gen_dimuonditrk_prompt",&gen_dimuonditrk_prompt,"gen_dimuonditrk_prompt/D");
        fivetracks_tree->Branch("gen_phi_prompt",&gen_phi_prompt,"gen_phi_prompt/D");
        fivetracks_tree->Branch("gen_jpsi_prompt",&gen_jpsi_prompt,"gen_jpsi_prompt/D");

        fivetracks_tree->Branch("gen_dimuonditrk_pt",&gen_dimuonditrk_pt,"gen_dimuonditrk_pt/D");
        fivetracks_tree->Branch("gen_phi_pt",&gen_phi_pt,"phigen_phi_pt_pt/D");
        fivetracks_tree->Branch("gen_jpsi_pt",&gen_jpsi_pt,"gen_jpsi_pt/D");

        fivetracks_tree->Branch("gen_dimuonditrk_p",&gen_dimuonditrk_p,"gen_dimuonditrk_p/D");
        fivetracks_tree->Branch("gen_phi_p",&gen_phi_p,"phigen_phi_p_p/D");
        fivetracks_tree->Branch("gen_jpsi_p",&gen_jpsi_p,"gen_jpsi_p/D");

        fivetracks_tree->Branch("gen_dimuonditrk_eta",&gen_dimuonditrk_eta,"gen_dimuonditrk_eta/D");
        fivetracks_tree->Branch("gen_phi_eta",&gen_phi_eta,"gen_phi_eta/D");
        fivetracks_tree->Branch("gen_jpsi_eta",&gen_jpsi_eta,"gen_jpsi_eta/D");

        fivetracks_tree->Branch("gen_dimuonditrk_phi",&gen_dimuonditrk_phi,"gen_dimuonditrk_phi/D");
        fivetracks_tree->Branch("gen_phi_phi",&gen_phi_phi,"gen_phi_phi/D");
        fivetracks_tree->Branch("gen_jpsi_phi",&gen_jpsi_phi,"gen_jpsi_phi/D");
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
  iEvent.getByToken(FiveTracksCollection_, fivetracks_cand_handle);

  edm::Handle<std::vector<reco::Vertex >> primaryVertices_handle;
  iEvent.getByToken(thePVs_, primaryVertices_handle);

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( TriggerResults_ , triggerResults_handle);

  numPrimaryVertices = primaryVertices_handle->size();
  run   = iEvent.id().run();
  event = iEvent.id().event();
  lumi  = iEvent.id().luminosityBlock();

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
gen_five_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_phi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_lowMuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_highMuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_highKaon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_lowKaon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
genThirdTrack_p4.SetPtEtaPhiM(0.,0.,0.,0.);

gen_dimuonditrk_pdg = 0;
gen_five_pdg = 0;

//std::cout << "Debug  1" << std::endl;
if(!OnlyGen_)
  if (!fivetracks_cand_handle.isValid()) std::cout<< "No five tracks information " << run << "," << event <<std::endl;

  if (fivetracks_cand_handle.isValid()) {


    //std::cout << "Debug  2" << std::endl;
    noFiveCandidates = (Int_t)(fivetracks_cand_handle->size());

    for (unsigned int i=0; i< fivetracks_cand_handle->size(); i++)
    {

      pat::CompositeCandidate five_cand;
      const pat::PackedCandidate *trackOne_cand, *trackTwo_cand, *trackThree_cand;
      const pat::CompositeCandidate *dimuonDiTrkOne_cand, *dimuonDiTrkTwo_cand, *dimuonDiTrkThree_cand, *first_five_ref;
      const pat::CompositeCandidate *dimuonditrk_cand, *dimuon_cand, *ditrackOne_cand;
      const pat::CompositeCandidate *triTrack_cand, *ditrackTwo_cand, *ditrackThree_cand;

      five_cand  = fivetracks_cand_handle->at(i);
      dimuonditrk_id = five_cand.userInt("dimuontt_index");

      first_five_ref = dynamic_cast<const pat::CompositeCandidate*>(five_cand.daughter("first_five_ref"));
      dimuonditrk_cand = dynamic_cast<const pat::CompositeCandidate*>(five_cand.daughter("dimuonditrack"));

      dimuon_id = dimuonditrk_cand->userInt("dimuon_id");
      p_id = dimuonditrk_cand->userInt("pId");
      m_id = dimuonditrk_cand->userInt("mId");

      dimuon_cand = dynamic_cast<const pat::CompositeCandidate*>(dimuonditrk_cand->daughter("dimuon"));
      ditrackOne_cand = dynamic_cast<const pat::CompositeCandidate*>(dimuonditrk_cand->daughter("ditrack"));

      const reco::Vertex bestPV = *(five_cand.userData<reco::Vertex>("bestPV"));
      const reco::Vertex cosPV  = *(five_cand.userData<reco::Vertex>("cosPV"));
      const reco::Vertex zPV    = *(five_cand.userData<reco::Vertex>("zPV"));
      const reco::Vertex bs     = *(five_cand.userData<reco::Vertex>("bS"));


      bestPV_X = bestPV.position().x();
      bestPV_Y = bestPV.position().y();
      bestPV_Z = bestPV.position().z();

      cosAlphaPV_X = cosPV.position().x();
      cosAlphaPV_Y = cosPV.position().y();
      cosAlphaPV_Z = cosPV.position().z();

      zPV_X = cosPV.position().x();
      zPV_Y = cosPV.position().y();
      zPV_Z = cosPV.position().z();

      bS_X = cosPV.position().x();
      bS_Y = cosPV.position().y();
      bS_Z = cosPV.position().z();

      dimuonditrk_vProb     = dimuonditrk_cand->userFloat("vProb");
      dimuonditrk_vChi2     = dimuonditrk_cand->userFloat("vChi2");
      dimuonditrk_nDof      = dimuonditrk_cand->userFloat("nDof");
      dimuonditrk_charge    = dimuonditrk_cand->charge();

      dimuonditrk_cosAlphaBS = dimuonditrk_cand->userFloat("cosAlphaBS");
      dimuonditrk_ctauPVBS = dimuonditrk_cand->userFloat("ctauPVBS");
      dimuonditrk_ctauErrPVBS = dimuonditrk_cand->userFloat("ctauErrPVBS");

      // tPFromPVBS = dimuonditrk_cand->userFloat("tPFromPVBS");
      // tMFromPVBS = dimuonditrk_cand->userFloat("tMFromPVBS");

      five_pt       = five_cand.pt();
      five_eta      = five_cand.eta();
      five_phi      = five_cand.phi();
      five_p        = five_cand.p();
      five_vProb    = five_cand.userFloat("vProb");
      five_nDof     = five_cand.userFloat("nDof");
      five_vChi2    = five_cand.userFloat("vChi2");
      five_charge   = five_cand.charge();

      five_m           = five_cand.mass();
      five_m_ref       = five_cand.userFloat("mass_ref_0");
      five_mass_ppk    = five_cand.userFloat("mass_ref_1");
      five_mass_kpp    = five_cand.userFloat("mass_ref_2");
      five_mass_pkp    = five_cand.userFloat("mass_ref_3");
      five_mass_ppp    = five_cand.userFloat("mass_ref_4");

      dca_m1t3 = five_cand.userFloat("dca_m1t3");
      dca_m2t3 = five_cand.userFloat("dca_m2t3");
      dca_t1t3 = five_cand.userFloat("dca_t1t3");
      dca_t2t3 = five_cand.userFloat("dca_t2t3");

      tTFromPV        = five_cand.userFloat("tTFromPV");
      tTFromPVDZ      = five_cand.userFloat("tTFromPVDZ");
      // tTFromPVBS      = five_cand.userFloat("tTFromPVBS");
      tTFromPVCA      = five_cand.userFloat("tTFromPVCA");



      dimuonditrk_cosAlpha = dimuonditrk_cand->userFloat("cosAlpha");
      dimuonditrk_ctauPV = dimuonditrk_cand->userFloat("ctauPV");
      dimuonditrk_ctauErrPV = dimuonditrk_cand->userFloat("ctauErrPV");

      tPFromPV = dimuonditrk_cand->userFloat("tPFromPV");
      tMFromPV = dimuonditrk_cand->userFloat("tMFromPV");

      dimuonditrk_cosAlphaDZ = dimuonditrk_cand->userFloat("cosAlphaDZ");
      dimuonditrk_ctauPVDZ = dimuonditrk_cand->userFloat("ctauPVDZ");
      dimuonditrk_ctauErrPVDZ = dimuonditrk_cand->userFloat("ctauErrPVDZ");

      tPFromPVDZ = dimuonditrk_cand->userFloat("tPFromPVDZ");
      tMFromPVDZ = dimuonditrk_cand->userFloat("tMFromPVDZ");

      dimuonditrk_cosAlphaCA = dimuonditrk_cand->userFloat("cosAlphaCA");
      dimuonditrk_ctauPVCA = dimuonditrk_cand->userFloat("ctauPVCA");
      dimuonditrk_ctauErrPVCA = dimuonditrk_cand->userFloat("ctauErrPVCA");

      tPFromPVCA = dimuonditrk_cand->userFloat("tPFromPVCA");
      tMFromPVCA = dimuonditrk_cand->userFloat("tMFromPVCA");


      dca_m1m2 = dimuonditrk_cand->userFloat("dca_m1m2");
      dca_m1t1 = dimuonditrk_cand->userFloat("dca_m1t1");
      dca_m1t2 = dimuonditrk_cand->userFloat("dca_m1t2");
      dca_m2t1 = dimuonditrk_cand->userFloat("dca_m2t1");
      dca_m2t2 = dimuonditrk_cand->userFloat("dca_m2t2");
      dca_t1t2 = dimuonditrk_cand->userFloat("dca_t1t2");

      dimuon_vProb        = dimuon_cand->userFloat("vProb");
      dimuon_vChi2        = dimuon_cand->userFloat("vNChi2");
      dimuon_DCA          = dimuon_cand->userFloat("DCA");
      dimuon_ctauPV       = dimuon_cand->userFloat("ppdlPV");
      dimuon_ctauErrPV    = dimuon_cand->userFloat("ppdlErrPV");
      dimuon_cosAlpha     = dimuon_cand->userFloat("cosAlpha");

      const pat::Muon *lowMuon, *highMuon;

      reco::Candidate::LorentzVector vhighMuon = dimuon_cand->daughter("highMuon")->p4();
      reco::Candidate::LorentzVector vlowMuon = dimuon_cand->daughter("lowMuon")->p4();

      // if (dimuon_cand->daughter("highMuon")->charge() < 0) {
      if (vhighMuon.pt() < vlowMuon.pt()) {
        vhighMuon = dimuon_cand->daughter("lowMuon")->p4();
        vlowMuon = dimuon_cand->daughter("highMuon")->p4();
        highMuon = dynamic_cast<const pat::Muon*>(dimuon_cand->daughter("highMuon"));
        lowMuon = dynamic_cast<const pat::Muon*>(dimuon_cand->daughter("lowMuon"));
      } else
      {
        lowMuon = dynamic_cast<const pat::Muon*>(dimuon_cand->daughter("highMuon"));
        highMuon = dynamic_cast<const pat::Muon*>(dimuon_cand->daughter("lowMuon"));
      }

      highMuon_pt  = highMuon->innerTrack()->pt();
      highMuon_eta  = highMuon->innerTrack()->eta();
      highMuon_phi  = highMuon->innerTrack()->phi();
      highMuon_dz  = highMuon->innerTrack()->dz();
      highMuon_dxy  = highMuon->innerTrack()->dxy();

      lowMuon_pt  = lowMuon->innerTrack()->pt();
      lowMuon_eta  = lowMuon->innerTrack()->eta();
      lowMuon_phi  = lowMuon->innerTrack()->phi();
      lowMuon_dz  = lowMuon->innerTrack()->dz();
      lowMuon_dxy  = lowMuon->innerTrack()->dxy();

      lowMuon_isTight    = lowMuon->isTightMuon(bestPV);
      lowMuon_isLoose    = lowMuon->isLooseMuon();
      lowMuon_isSoft     = lowMuon->isSoftMuon(bestPV);
      lowMuon_isMedium   = lowMuon->isMediumMuon();
      lowMuon_isHighPt   = lowMuon->isHighPtMuon(bestPV);
      lowMuon_isTracker  = lowMuon->isTrackerMuon();
      lowMuon_isGlobal   = lowMuon->isGlobalMuon();
      lowMuon_NPixelHits = lowMuon->innerTrack()->hitPattern().numberOfValidPixelHits();
      lowMuon_NStripHits = lowMuon->innerTrack()->hitPattern().numberOfValidStripHits();
      lowMuon_NTrackhits = lowMuon->innerTrack()->hitPattern().numberOfValidTrackerHits();
      lowMuon_NBPixHits  = lowMuon->innerTrack()->hitPattern().numberOfValidStripHits();
      lowMuon_NPixLayers = lowMuon->innerTrack()->hitPattern().pixelLayersWithMeasurement();
      lowMuon_NTraLayers = lowMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement();
      lowMuon_NStrLayers = lowMuon->innerTrack()->hitPattern().stripLayersWithMeasurement();
      lowMuon_NBPixLayers = lowMuon->innerTrack()->hitPattern().pixelBarrelLayersWithMeasurement();

      highMuon_isTight    = highMuon->isTightMuon(bestPV);
      highMuon_isLoose    = highMuon->isLooseMuon();
      highMuon_isSoft     = highMuon->isSoftMuon(bestPV);
      highMuon_isMedium   = highMuon->isMediumMuon();
      highMuon_isHighPt   = highMuon->isHighPtMuon(bestPV);
      highMuon_isTracker  = highMuon->isTrackerMuon();
      highMuon_isGlobal   = highMuon->isGlobalMuon();
      highMuon_NPixelHits = highMuon->innerTrack()->hitPattern().numberOfValidPixelHits();
      highMuon_NStripHits = highMuon->innerTrack()->hitPattern().numberOfValidStripHits();
      highMuon_NTrackhits = highMuon->innerTrack()->hitPattern().numberOfValidTrackerHits();
      highMuon_NBPixHits  = highMuon->innerTrack()->hitPattern().numberOfValidStripHits();
      highMuon_NPixLayers = highMuon->innerTrack()->hitPattern().pixelLayersWithMeasurement();
      highMuon_NTraLayers = highMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement();
      highMuon_NStrLayers = highMuon->innerTrack()->hitPattern().stripLayersWithMeasurement();
      highMuon_NBPixLayers = highMuon->innerTrack()->hitPattern().pixelBarrelLayersWithMeasurement();

      lowMuon_type       = lowMuon->type();
      highMuon_type       = highMuon->type();

      lowMuon_p4.SetPtEtaPhiM(vhighMuon.pt(), vhighMuon.eta(), vhighMuon.phi(), vhighMuon.mass());
      highMuon_p4.SetPtEtaPhiM(vlowMuon.pt(), vlowMuon.eta(), vlowMuon.phi(), vlowMuon.mass());

      highTrackMatch   = (Double_t)dimuonditrk_cand->userInt("highKaonMatch");
      lowTrackMatch    = (Double_t)dimuonditrk_cand->userInt("lowKaonMatch");
      lowMuonMatch     = (Double_t)dimuon_cand->userInt("highMuonTMatch");
      highMuonMatch    = (Double_t)dimuon_cand->userInt("lowMuonTMatch");
      thirdTrackMatch  = (Double_t)five_cand.userInt("thirdKaonMatch");

      trackOne_cand   = dynamic_cast<const pat::PackedCandidate*>(five_cand.daughter("trackOne"));
      trackTwo_cand   = dynamic_cast<const pat::PackedCandidate*>(five_cand.daughter("trackTwo"));
      trackThree_cand = dynamic_cast<const pat::PackedCandidate*>(five_cand.daughter("trackThree"));

      float pionmass = 0.13957061, protonmass = 0.93827208;

      highKaon_p4.SetPtEtaPhiM(trackOne_cand->pt(), trackOne_cand->eta(), trackOne_cand->phi(), trackOne_cand->mass());
      lowKaon_p4.SetPtEtaPhiM(trackTwo_cand->pt(), trackTwo_cand->eta(), trackTwo_cand->phi(), trackTwo_cand->mass());
      thirdKaon_p4.SetPtEtaPhiM(trackThree_cand->pt(), trackThree_cand->eta(), trackThree_cand->phi(), trackThree_cand->mass());

      highPion_p4.SetPtEtaPhiM(trackOne_cand->pt(), trackOne_cand->eta(), trackOne_cand->phi(), pionmass);
      lowPion_p4.SetPtEtaPhiM(trackTwo_cand->pt(), trackTwo_cand->eta(), trackTwo_cand->phi(), pionmass);
      thirdPion_p4.SetPtEtaPhiM(trackThree_cand->pt(), trackThree_cand->eta(), trackThree_cand->phi(), pionmass);

      highProton_p4.SetPtEtaPhiM(trackOne_cand->pt(), trackOne_cand->eta(), trackOne_cand->phi(), protonmass);
      lowProton_p4.SetPtEtaPhiM(trackTwo_cand->pt(), trackTwo_cand->eta(), trackTwo_cand->phi(),protonmass);
      thirdProton_p4.SetPtEtaPhiM(trackThree_cand->pt(), trackThree_cand->eta(), trackThree_cand->phi(), protonmass);

      lowMuon_p4.SetPtEtaPhiM(vhighMuon.pt(), vhighMuon.eta(), vhighMuon.phi(), vhighMuon.mass());
      highMuon_p4.SetPtEtaPhiM(vlowMuon.pt(), vlowMuon.eta(), vlowMuon.phi(), vlowMuon.mass());

      dimuonDiTrkOne_cand = dynamic_cast<const pat::CompositeCandidate*>(first_five_ref->daughter("dimuonDiTrackOne"));
      dimuonDiTrkTwo_cand = dynamic_cast<const pat::CompositeCandidate*>(first_five_ref->daughter("dimuonDiTrackTwo"));
      dimuonDiTrkThree_cand = dynamic_cast<const pat::CompositeCandidate*>(first_five_ref->daughter("dimuonDiTrackThree"));

      ditrackTwo_cand    = dynamic_cast<const pat::CompositeCandidate*>(dimuonDiTrkTwo_cand->daughter("ditrack"));
      ditrackThree_cand  = dynamic_cast<const pat::CompositeCandidate*>(dimuonDiTrkThree_cand->daughter("ditrack"));

      triTrack_cand = dynamic_cast<const pat::CompositeCandidate*>(first_five_ref->daughter("triTrack"));

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

      ditrack_p4.SetPtEtaPhiM(ditrackOne_cand->pt(),ditrackOne_cand->eta(),ditrackOne_cand->phi(),ditrackOne_cand->mass());
      ditrack_m = ditrackOne_cand->mass();
      ditrackOne_pt = ditrackOne_cand->pt();
      ditrackOne_eta = ditrackOne_cand->eta();
      ditrackOne_phi = ditrackOne_cand->phi();
      ditrackOne_p = ditrackOne_cand->p();

      ditrackTwo_pt = ditrackTwo_cand->pt();
      ditrackTwo_eta = ditrackTwo_cand->eta();
      ditrackTwo_phi = ditrackTwo_cand->phi();
      ditrackTwo_p = ditrackTwo_cand->p();

      ditrackThree_pt = ditrackThree_cand->pt();
      ditrackThree_eta = ditrackThree_cand->eta();
      ditrackThree_phi = ditrackThree_cand->phi();
      ditrackThree_p = ditrackThree_cand->p();

      //trackOne_cand->SetPtEtaPhiM(dimuon_cand->pt(),dimuon_cand->eta(),dimuon_cand->phi(),dimuon_cand->mass());

      highTrackMuonDR = dimuonditrk_cand->userFloat("highKaonMuonDR");
      highTrackMuonDP = dimuonditrk_cand->userFloat("highKaonMuonDP");
      highTrackMuonDPt = dimuonditrk_cand->userFloat("highKaonMuonDPt");

      lowTrackMuonDR  = dimuonditrk_cand->userFloat("lowKaonMuonDR");
      lowTrackMuonDP  = dimuonditrk_cand->userFloat("lowKaonMuonDP");
      lowTrackMuonDPt  = dimuonditrk_cand->userFloat("lowKaonMuonDPt");

      thirdTrackMuonDR   = five_cand.userFloat("thirdTrackMuonDR");
      thirdTrackMuonDP   = five_cand.userFloat("thirdTrackMuonDP");
      thirdTrackMuonDPt  = five_cand.userFloat("thirdTrackMuonDPt");


      dimuonditrk_vx    = dimuonditrk_cand->userFloat("vtxX");
      dimuonditrk_vy    = dimuonditrk_cand->userFloat("vtxY");
      dimuonditrk_vz    = dimuonditrk_cand->userFloat("vtxZ");

      five_vx    = five_cand.userFloat("vtxX");
      five_vy    = five_cand.userFloat("vtxY");
      five_vz    = five_cand.userFloat("vtxZ");

      highTrack_pt      = trackOne_cand->pt();
      highTrack_eta     = trackOne_cand->eta();
      highTrack_phi     = trackOne_cand->phi();
      highTrack_charge  = trackOne_cand->charge();
      thirdTrack_dz     = trackOne_cand->bestTrack()->dz();
      thirdTrack_dxy    = trackOne_cand->bestTrack()->dxy();

      highTrack_NPixelHits = trackOne_cand->bestTrack()->hitPattern().numberOfValidPixelHits();
      highTrack_NStripHits = trackOne_cand->bestTrack()->hitPattern().numberOfValidStripHits();
      highTrack_NTrackhits = trackOne_cand->bestTrack()->hitPattern().numberOfValidTrackerHits();
      highTrack_NBPixHits  = trackOne_cand->bestTrack()->hitPattern().numberOfValidStripHits();
      highTrack_NPixLayers = trackOne_cand->bestTrack()->hitPattern().pixelLayersWithMeasurement();
      highTrack_NTraLayers = trackOne_cand->bestTrack()->hitPattern().trackerLayersWithMeasurement();
      highTrack_NStrLayers = trackOne_cand->bestTrack()->hitPattern().stripLayersWithMeasurement();
      highTrack_NBPixLayers = trackOne_cand->bestTrack()->hitPattern().pixelBarrelLayersWithMeasurement();

      lowTrack_pt = trackTwo_cand->pt();
      lowTrack_eta = trackTwo_cand->eta();
      lowTrack_phi = trackTwo_cand->phi();
      lowTrack_charge = trackTwo_cand->charge();
      lowTrack_dz     = trackTwo_cand->bestTrack()->dz();
      lowTrack_dxy    = trackTwo_cand->bestTrack()->dxy();

      lowTrack_NPixelHits = trackTwo_cand->bestTrack()->hitPattern().numberOfValidPixelHits();
      lowTrack_NStripHits = trackTwo_cand->bestTrack()->hitPattern().numberOfValidStripHits();
      lowTrack_NTrackhits = trackTwo_cand->bestTrack()->hitPattern().numberOfValidTrackerHits();
      lowTrack_NBPixHits  = trackTwo_cand->bestTrack()->hitPattern().numberOfValidStripHits();
      lowTrack_NPixLayers = trackTwo_cand->bestTrack()->hitPattern().pixelLayersWithMeasurement();
      lowTrack_NTraLayers = trackTwo_cand->bestTrack()->hitPattern().trackerLayersWithMeasurement();
      lowTrack_NStrLayers = trackTwo_cand->bestTrack()->hitPattern().stripLayersWithMeasurement();
      lowTrack_NBPixLayers = trackTwo_cand->bestTrack()->hitPattern().pixelBarrelLayersWithMeasurement();


      thirdTrack_pt     = trackThree_cand->pt();
      thirdTrack_eta    = trackThree_cand->eta();
      thirdTrack_phi    = trackThree_cand->phi();
      thirdTrack_charge = trackThree_cand->charge();
      thirdTrack_dz     = trackThree_cand->bestTrack()->dz();
      thirdTrack_dxy    = trackThree_cand->bestTrack()->dxy();

      thirdTrack_NPixelHits = trackThree_cand->bestTrack()->hitPattern().numberOfValidPixelHits();
      thirdTrack_NStripHits = trackThree_cand->bestTrack()->hitPattern().numberOfValidStripHits();
      thirdTrack_NTrackhits = trackThree_cand->bestTrack()->hitPattern().numberOfValidTrackerHits();
      thirdTrack_NBPixHits  = trackThree_cand->bestTrack()->hitPattern().numberOfValidStripHits();
      thirdTrack_NPixLayers = trackThree_cand->bestTrack()->hitPattern().pixelLayersWithMeasurement();
      thirdTrack_NTraLayers = trackThree_cand->bestTrack()->hitPattern().trackerLayersWithMeasurement();
      thirdTrack_NStrLayers = trackThree_cand->bestTrack()->hitPattern().stripLayersWithMeasurement();
      thirdTrack_NBPixLayers = trackThree_cand->bestTrack()->hitPattern().pixelBarrelLayersWithMeasurement();

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

      triTrack_m      = triTrack_cand->mass();
      triTrack_pt     = triTrack_cand->pt();
      triTrack_eta    = triTrack_cand->eta();
      triTrack_phi    = triTrack_cand->phi();
      triTrack_charge = triTrack_cand->charge();

      if(IsMC_ || OnlyGen_)
      {

        gen_five_p4.SetPtEtaPhiM(-1.0,0.0,0.0,-0.01);
        gen_dimuonditrk_p4.SetPtEtaPhiM(-1.0,0.,0.,-0.01);
        gen_jpsi_p4.SetPtEtaPhiM(-1.0,0.0,0.0,-0.01);
        gen_phi_p4.SetPtEtaPhiM(-1.0,0.0,0.0,-0.01);
        gen_highKaon_p4.SetPtEtaPhiM(-1.0,0.0,0.0,-0.01);
        gen_lowMuon_p4.SetPtEtaPhiM(-1.0,0.0,0.0,-0.01);
        gen_highMuon_p4.SetPtEtaPhiM(-1.0,0.0,0.0,-0.01);
        gen_lowKaon_p4.SetPtEtaPhiM(-1.0,0.0,0.0,-0.01);
        genThirdTrack_p4.SetPtEtaPhiM(-1.0,0.0,0.0,-0.01);

        gen_dimuonditrk_pdg = 0.0;
        gen_phi_pdg         = 0.0;
        gen_jpsi_pdg        = 0.0;
        gen_five_pdg        = 0.0;

        gen_lowMuon_pdg     = 0.0;
        gen_highMuon_pdg    = 0.0;
        gen_highKaon_pdg    = 0.0;
        gen_lowKaon_pdg     = 0.0;
        genThirdTrack_pdg     = 0.0;

        gen_lowMuon_mompdg     = 0.0;
        gen_highMuon_mompdg    = 0.0;
        gen_highKaon_mompdg    = 0.0;
        gen_lowKaon_mompdg     = 0.0;

        gen_lowMuon_status     = 0.0;
        gen_highMuon_status    = 0.0;
        gen_highKaon_status    = 0.0;
        gen_lowKaon_status     = 0.0;

        gen_dimuonditrk_prompt = 0.0;
        gen_five_prompt     = 0.0;
        gen_phi_prompt      = 0.0;
        gen_jpsi_prompt     = 0.0;

        gen_dimuonditrk_pt  = 0.0;
        gen_phi_pt          = 0.0;
        gen_jpsi_pt         = 0.0;
        gen_five_pt         = 0.0;

        gen_dimuonditrk_p   = 0.0;
        gen_phi_p           = 0.0;
        gen_jpsi_p          = 0.0;
        gen_five_p          = 0.0;

        gen_dimuonditrk_eta = 0.0;
        gen_phi_eta         = 0.0;
        gen_jpsi_eta        = 0.0;

        gen_dimuonditrk_phi = 0.0;
        gen_phi_phi         = 0.0;
        gen_jpsi_phi        = 0.0;

        gen_lowMuon_pt     = 0.0;
        gen_highMuon_pt    = 0.0;
        gen_highKaon_pt    = 0.0;
        gen_lowKaon_pt     = 0.0;
        genThirdTrack_pt  = 0.0;

        gen_lowMuon_p     = 0.0;
        gen_highMuon_p    = 0.0;
        gen_highKaon_p    = 0.0;
        gen_lowKaon_p     = 0.0;
        genThirdTrack_p  = 0.0;

        gen_lowMuon_eta     = 0.0;
        gen_highMuon_eta    = 0.0;
        gen_highKaon_eta    = 0.0;
        gen_lowKaon_eta     = 0.0;
        genThirdTrack_eta  = 0.0;

        gen_lowMuon_phi     = 0.0;
        gen_highMuon_phi    = 0.0;
        gen_highKaon_phi    = 0.0;
        gen_lowKaon_phi     = 0.0;

        reco::GenParticleRef genhighMuon  = highMuon->genParticleRef();
        reco::GenParticleRef genlowMuon   = lowMuon->genParticleRef();

        const reco::GenParticle *genhighKaon,*genlowKaon,*genThirdTrack;
        reco::GenParticleRef phiMomHigh, phiMomLow, jpsiMomHigh, jpsiMomLow;
        reco::GenParticleRef jpsiMom, phiMom, thirdMom;

        Double_t hasHighGen  = dimuonditrk_cand->userFloat("hasHighGen");
        Double_t hasLowGen   = dimuonditrk_cand->userFloat("hasLowGen");
        Double_t hasThirdGen = five_cand.userFloat("hasThirdGen");

        if(hasHighGen>0.0)
          genhighKaon = dynamic_cast <const reco::GenParticle *>(dimuonditrk_cand->daughter("highKaonGen"));
        if(hasLowGen>0.0)
          genlowKaon = dynamic_cast <const reco::GenParticle *>(dimuonditrk_cand->daughter("lowKaonGen"));
        if(hasThirdGen>0.0)
          genThirdTrack = dynamic_cast <const reco::GenParticle *>(five_cand.daughter("thirdTrackGen"));

        if(hasHighGen>0.0)
        {
          gen_highKaon_p4.SetPtEtaPhiM(genhighKaon->pt(),genhighKaon->eta(),genhighKaon->phi(),genhighKaon->mass());
          if(genhighKaon->numberOfMothers()>0)
            phiMomHigh  = genhighKaon->motherRef();

          gen_highKaon_pdg     = (Double_t)genhighKaon->pdgId();

          if(phiMomHigh.isNonnull() && genhighKaon->numberOfMothers()>0)
            gen_highKaon_mompdg  = phiMomHigh->pdgId();

          gen_highKaon_status  = (Double_t)genhighKaon->status();
          gen_highKaon_pt      = (Double_t)genhighKaon->pt();
          gen_highKaon_p       = (Double_t)genhighKaon->p();
          gen_highKaon_eta     = (Double_t)genhighKaon->eta();
          gen_highKaon_phi     = (Double_t)genhighKaon->phi();

        }

        if(hasLowGen>0.0)
        {
          gen_lowKaon_p4.SetPtEtaPhiM(genlowKaon->pt(),genlowKaon->eta(),genlowKaon->phi(),genlowKaon->mass());

          if(genlowKaon->numberOfMothers()>0)
            phiMomLow  = genlowKaon->motherRef();

          gen_lowKaon_pdg     = (Double_t)genlowKaon->pdgId();

          if(phiMomLow.isNonnull() && genlowKaon->numberOfMothers()>0)
            gen_lowKaon_mompdg  = phiMomLow->pdgId();

          gen_lowKaon_status  = (Double_t)genlowKaon->status();
          gen_lowKaon_pt      = (Double_t)genlowKaon->pt();
          gen_lowKaon_p       = (Double_t)genlowKaon->p();
          gen_lowKaon_eta     = (Double_t)genlowKaon->eta();
          gen_lowKaon_phi     = (Double_t)genlowKaon->phi();
        }

        if(hasThirdGen>0.0)
        {
          genThirdTrack_p4.SetPtEtaPhiM(genThirdTrack->pt(),genThirdTrack->eta(),genThirdTrack->phi(),genThirdTrack->mass());

          if(genThirdTrack->numberOfMothers()>0)
            thirdMom  = genThirdTrack->motherRef();

          genThirdTrack_pdg     = (Double_t)genThirdTrack->pdgId();

          if(phiMomLow.isNonnull() && genThirdTrack->numberOfMothers()>0)
            genThirdTrack_mompdg  = thirdMom->pdgId();

          genThirdTrack_status  = (Double_t)genThirdTrack->status();
          genThirdTrack_pt      = (Double_t)genThirdTrack->pt();
          genThirdTrack_p       = (Double_t)genThirdTrack->p();
          genThirdTrack_eta     = (Double_t)genThirdTrack->eta();
          genThirdTrack_phi     = (Double_t)genThirdTrack->phi();


        }


        if(genhighMuon.isNonnull())
        {
          gen_highMuon_p4.SetPtEtaPhiM(genhighMuon->pt(),genhighMuon->eta(),genhighMuon->phi(),genhighMuon->mass());
          if(genhighMuon->numberOfMothers()>0)
          jpsiMomHigh  = genhighMuon->motherRef();

          gen_highMuon_pdg     = (Double_t)genhighMuon->pdgId();

          if(jpsiMomHigh.isNonnull() && genhighMuon->numberOfMothers()>0)
            gen_highMuon_mompdg  = jpsiMomHigh->pdgId();

          gen_highMuon_status  = (Double_t)genhighMuon->status();
          gen_highMuon_pt      = (Double_t)genhighMuon->pt();
          gen_highMuon_p       = (Double_t)genhighMuon->p();
          gen_highMuon_eta     = (Double_t)genhighMuon->eta();
          gen_highMuon_phi     = (Double_t)genhighMuon->phi();
        }

        if(genlowMuon.isNonnull())
        {
          gen_lowMuon_p4.SetPtEtaPhiM(genlowMuon->pt(),genlowMuon->eta(),genlowMuon->phi(),genlowMuon->mass());
          if(genlowMuon->numberOfMothers()>0)
            jpsiMomLow  = genlowMuon->motherRef();

          gen_lowMuon_pdg     = (Double_t)genlowMuon->pdgId();

          if(jpsiMomLow.isNonnull() && genlowMuon->numberOfMothers()>0)
            gen_lowMuon_mompdg  = jpsiMomLow->pdgId();

          gen_lowMuon_status  = (Double_t)genlowMuon->status();
          gen_lowMuon_pt      = (Double_t)genlowMuon->pt();
          gen_lowMuon_p       = (Double_t)genlowMuon->p();
          gen_lowMuon_eta     = (Double_t)genlowMuon->eta();
          gen_lowMuon_phi     = (Double_t)genlowMuon->phi();
        }

        bool samePhiMom = false, samejpsiMom = false;

        if(phiMomLow.isNonnull() && phiMomHigh.isNonnull())
        {
          samePhiMom = (phiMomHigh == phiMomLow);

          if(samePhiMom)
          {
            gen_phi_p4.SetPtEtaPhiM(phiMomHigh->pt(),phiMomHigh->eta(),phiMomHigh->phi(),phiMomHigh->mass());
            gen_phi_pdg     = phiMomHigh->pdgId();
            gen_phi_prompt  = phiMomHigh->isPromptDecayed();
            gen_phi_p       = phiMomHigh->p();
            gen_phi_pt      = phiMomHigh->pt();
            gen_phi_eta     = phiMomHigh->eta();
            gen_phi_phi     = phiMomHigh->phi();
          }
        }

        if(jpsiMomLow.isNonnull() && jpsiMomHigh.isNonnull())
        {
          samejpsiMom = (jpsiMomHigh == jpsiMomLow);

          if(samejpsiMom)
          {
            gen_jpsi_p4.SetPtEtaPhiM(jpsiMomHigh->pt(),jpsiMomHigh->eta(),jpsiMomHigh->phi(),jpsiMomHigh->mass());
            gen_jpsi_pdg     = jpsiMomHigh->pdgId();
            gen_jpsi_prompt  = jpsiMomHigh->isPromptDecayed();
            gen_jpsi_p       = jpsiMomHigh->p();
            gen_jpsi_pt      = jpsiMomHigh->pt();
            gen_jpsi_eta     = jpsiMomHigh->eta();
            gen_jpsi_phi     = jpsiMomHigh->phi();
          }
        }

        if(samejpsiMom && samePhiMom  && jpsiMomHigh->numberOfMothers()>0 && phiMomHigh->numberOfMothers()>0)
        {

          jpsiMom = jpsiMomHigh->motherRef();
          phiMom  = phiMomHigh->motherRef();

          //check for X->J/Psi Phi Trk
          if(jpsiMom==phiMom && jpsiMom==thirdMom && jpsiMom.isNonnull() && phiMom.isNonnull() && thirdMom.isNonnull())
          {
            gen_five_p4.SetPtEtaPhiM(jpsiMom->pt(),jpsiMom->eta(),jpsiMom->phi(),jpsiMom->mass());
            gen_five_pdg = (Double_t) jpsiMom->pdgId();
            gen_five_prompt = (Double_t) jpsiMom->isPromptDecayed();
            gen_five_p = (Double_t) jpsiMom->p();
            gen_five_pt = (Double_t) jpsiMom->pt();
            gen_five_eta = (Double_t) jpsiMom->eta();
            gen_five_phi = (Double_t) jpsiMom->phi();
          } else if (jpsiMom==phiMom && jpsiMom!=thirdMom && jpsiMom.isNonnull() && phiMom.isNonnull() )
          {

            gen_dimuonditrk_p4.SetPtEtaPhiM(jpsiMom->pt(),jpsiMom->eta(),jpsiMom->phi(),jpsiMom->mass());
            gen_dimuonditrk_pdg = (Double_t) jpsiMom->pdgId();
            gen_dimuonditrk_prompt = (Double_t) jpsiMom->isPromptDecayed();
            gen_dimuonditrk_p = (Double_t) jpsiMom->p();
            gen_dimuonditrk_pt = (Double_t) jpsiMom->pt();
            gen_dimuonditrk_eta = (Double_t) jpsiMom->eta();
            gen_dimuonditrk_phi = (Double_t) jpsiMom->phi();

          }

        }

      } //IsMC || onlyGen


      fivetracks_tree->Fill();

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
