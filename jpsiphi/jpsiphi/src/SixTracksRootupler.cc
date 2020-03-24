/*
   Package:    SixTracksRootupler
   Class:      SixTracksRootupler

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

class SixTracksRootupler : public edm::EDAnalyzer {
   public:
      explicit SixTracksRootupler(const edm::ParameterSet&);
      ~SixTracksRootupler() override;

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
      UInt_t isTriggerMatched(pat::CompositeCandidate *mumu_cand);

  // ----------member data ---------------------------
  std::string file_name;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> SixTracksCollection_;
  edm::EDGetTokenT<edm::TriggerResults> TriggerResults_;
  edm::EDGetTokenT<reco::VertexCollection> thePVs_;
  bool IsMC_,OnlyGen_ ;
  UInt_t motherpdgid_,phipdgid_,jpspdgid_;
  std::vector<std::string>  HLTs_;
  std::vector<std::string>  HLTFilters_;
  std::string TreeName_;

  std::vector < Double_t > fiveTracksMass, fiveTracksMassRef, triTrackMass;
  std::vector < Double_t > fiveTracksVProb, fiveTracksVNDof, fiveTracksVChi2;
  std::vector < Double_t > fiveTracksCosAlpha, fiveTracksCTau, fiveTracksCTauErr;
  std::vector < Double_t > psiPrimeSame, psiPrimeMixed;
  std::vector < Double_t > psiPrimeSame_ditrack, psiPrimeMixed_ditrack;
  std::vector < Double_t > dimuonDiTrkOne, dimuonDiTrkTwo, dimuonDiTrkThree;
  std::vector < Double_t > diTrackOne, diTrackTwo, diTrackThree;
  std::vector < Double_t > trackOneMass, trackTwoMass, trackThreeMass;
  std::vector < Double_t > psiPrimeSame_p_m, psiPrimeSame_m_m, psiPrimeMixed_p_m, psiPrimeMixed_m_m;

  Double_t testMass;


  UInt_t run, event, lumi, numPrimaryVertices, trigger;
  UInt_t mumukk_id, five_id, dimuon_id, p1_id, p2_id, m2_id, m1_id;
  TLorentzVector mumukk_p4, dimuon_p4, ditrack_p4, six_p4, five_p4;
  TLorentzVector lowPion_p4, lowProton_p4, highProton_p4, highPion_p4, thirdProton_p4, thirdPion_p4;
  TLorentzVector lowMuon_p4, highMuon_p4, lowKaon_p4, thirdKaon_p4, highKaon_p4;
  TLorentzVector fourthProton_p4, fourthPion_p4, fourthKaon_p4;


  Double_t five_m, five_m_ref, five_mass_ppk, five_mass_kpp, five_mass_pkp, five_mass_ppp, five_pt, five_eta, five_phi, five_p;
  Double_t five_cosAlpha, five_ctauPV, five_ctauErrPV, five_cosAlphaCA, five_ctauPVCA, five_ctauErrPVCA, five_cosAlphaDZ;
  Double_t five_ctauPVDZ, five_ctauErrPVDZ, five_cosAlphaBS, five_ctauPVBS, five_ctauErrPVBS;
  Double_t five_vProb, five_nDof, five_vChi2, five_vx, five_vy, five_vz;

  Double_t fiveOne_pt, fiveOne_eta, fiveOne_phi, fiveOne_p;
  Double_t fiveTwo_pt, fiveTwo_eta, fiveTwo_phi, fiveTwo_p;
  Double_t fiveThree_pt, fiveThree_eta, fiveThree_phi, fiveThree_p;
  Double_t fiveFour_pt, fiveFour_eta, fiveFour_phi, fiveFour_p;
  Double_t fiveFive_pt, fiveFive_eta, fiveFive_phi, fiveFive_p;

  Int_t mumukk_charge, five_charge;

  Double_t highTrackMatch, lowTrackMatch, thirdTrackMatch, fourthTrackMatch, highMuonMatch, lowMuonMatch;

  Double_t mumukk_vProb, mumukk_vChi2;
  Double_t mumukk_pt, mumukk_eta, mumukk_phi, mumukk_y, mumukk_vx, mumukk_vy, mumukk_vz, mumukk_p;
  Double_t highKaon_pt,lowKaon_pt,highMuon_pt,lowMuon_pt;

  Double_t highTrackMuonDR, highTrackMuonDP, highTrackMuonDPt;
  Double_t lowTrackMuonDR, lowTrackMuonDP, lowTrackMuonDPt;
  Double_t thirdTrackMuonDR, thirdTrackMuonDP, thirdTrackMuonDPt;
  Double_t fourthTrackMuonDR, fourthTrackMuonDP, fourthTrackMuonDPt;

  Double_t highMuon_eta, lowMuon_eta, highMuon_phi, lowMuon_phi, mumukk_m;
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

  Double_t mumukk_cosAlpha, mumukk_ctauPV, mumukk_ctauErrPV;
  Double_t mumukk_cosAlphaDZ, mumukk_ctauPVDZ, mumukk_ctauErrPVDZ;
  Double_t mumukk_cosAlphaBS, mumukk_ctauPVBS, mumukk_ctauErrPVBS;
  Double_t mumukk_cosAlphaCA, mumukk_ctauPVCA,mumukk_ctauErrPVCA;

  Double_t highTrackFromPVBS, lowTrackFromPVBS, highTrackFromPVCA, lowTrackFromPVCA;
  Double_t highTrackFromPVDZ, lowTrackFromPVDZ, highTrackFromPV, lowTrackFromPV;
  Double_t thirdTrackFromPV, thirdTrackFromPVDZ, thirdTrackFromPVBS, thirdTrackFromPVCA;
  Double_t fourthTrackFromPV, fourthTrackFromPVCA, fourthTrackFromPVDZ, fourthTrackFromPVBS;

  Double_t dca_m1m2, dca_m1t1, dca_m1t2, dca_m2t1, dca_m2t2, dca_t1t2;
  Double_t dca_m1t3, dca_m2t3, dca_t1t3, dca_t2t3;
  Double_t dca_m1t4, dca_m2t4, dca_t1t4, dca_t2t4, dca_t3t4;

  Double_t dimuon_vProb, dimuon_vChi2, dimuon_DCA, dimuon_ctauPV, dimuon_ctauErrPV, dimuon_cosAlpha;

  Double_t highTrack_pt, highTrack_eta, highTrack_phi, highTrack_charge;
  Double_t lowTrack_pt, lowTrack_eta, lowTrack_phi, lowTrack_charge;
  Double_t thirdTrack_pt, thirdTrack_eta, thirdTrack_phi, thirdTrack_charge;
  Double_t thirdTrack_dz, thirdTrack_dxy, lowTrack_dz, lowTrack_dxy, highTrack_dz, highTrack_dxy;

  Double_t fourthTrack_pt, fourthTrack_eta, fourthTrack_phi, fourthTrack_charge;
  Double_t fourthTrack_dz, fourthTrack_dxy;

  Double_t dimuonDiTrkOne_pt, dimuonDiTrkOne_eta, dimuonDiTrkOne_phi, dimuonDiTrkOne_charge, dimuonDiTrkOne_p;
  Double_t dimuonDiTrkTwo_pt, dimuonDiTrkTwo_eta, dimuonDiTrkTwo_phi, dimuonDiTrkTwo_charge, dimuonDiTrkTwo_p;
  Double_t dimuonDiTrkThree_pt, dimuonDiTrkThree_eta, dimuonDiTrkThree_phi, dimuonDiTrkThree_charge, dimuonDiTrkThree_p;
  Double_t dimuonDiTrkFour_pt, dimuonDiTrkFour_eta, dimuonDiTrkFour_phi, dimuonDiTrkFour_charge, dimuonDiTrkFour_p;
  Double_t dimuonDiTrkFive_pt, dimuonDiTrkFive_eta, dimuonDiTrkFive_phi, dimuonDiTrkFive_charge, dimuonDiTrkFive_p;
  Double_t dimuonDiTrkSix_pt, dimuonDiTrkSix_eta, dimuonDiTrkSix_phi, dimuonDiTrkSix_charge, dimuonDiTrkSix_p;

  Double_t mumukk_nDof,mumukk_m_rf,mumukk_m_rf_c,mumukk_m_rf_d_c;
  Int_t highTrack_NPixelHits, highTrack_NStripHits, highTrack_NTrackhits, highTrack_NBPixHits, highTrack_NPixLayers;
  Int_t highTrack_NTraLayers, highTrack_NStrLayers, highTrack_NBPixLayers, lowTrack_NPixelHits, lowTrack_NStripHits;
  Int_t lowTrack_NTrackhits, lowTrack_NBPixHits, lowTrack_NPixLayers, lowTrack_NTraLayers, lowTrack_NStrLayers, lowTrack_NBPixLayers;
  Int_t thirdTrack_NPixelHits, thirdTrack_NStripHits, thirdTrack_NTrackhits, thirdTrack_NBPixHits, thirdTrack_NPixLayers;
  Int_t thirdTrack_NTraLayers, thirdTrack_NStrLayers, thirdTrack_NBPixLayers;
  Int_t fourthTrack_NPixelHits, fourthTrack_NStripHits, fourthTrack_NTrackhits, fourthTrack_NBPixHits, fourthTrack_NPixLayers;
  Int_t fourthTrack_NTraLayers, fourthTrack_NStrLayers, fourthTrack_NBPixLayers;

  Double_t dimuon_m, dimuon_pt, dimuon_eta, dimuon_phi, dimuon_p;
  Double_t ditrack_m, diTrackOne_pt, diTrackOne_eta, diTrackOne_phi, diTrackOne_p;
  Double_t diTrackTwo_pt, diTrackTwo_eta, diTrackTwo_phi, diTrackTwo_p;
  Double_t diTrackThree_pt, diTrackThree_eta, diTrackThree_phi, diTrackThree_p;
  Double_t diTrackFour_pt, diTrackFour_eta, diTrackFour_phi, diTrackFour_p;
  Double_t diTrackFive_pt, diTrackFive_eta, diTrackFive_phi, diTrackFive_p;
  Double_t diTrackSix_pt, diTrackSix_eta, diTrackSix_phi, diTrackSix_p;

  Double_t diTrackOne_kk, diTrackOne_pp, diTrackOne_pk, diTrackOne_kp;
  Double_t diTrackTwo_kk, diTrackTwo_pp, diTrackTwo_pk, diTrackTwo_kp;
  Double_t diTrackThree_kk, diTrackThree_pp, diTrackThree_pk, diTrackThree_kp;
  Double_t diTrackFour_kk, diTrackFour_pp, diTrackFour_pk, diTrackFour_kp;
  Double_t diTrackFive_kk, diTrackFive_pp, diTrackFive_pk, diTrackFive_kp;
  Double_t diTrackSix_kk, diTrackSix_pp, diTrackSix_pk, diTrackSix_kp;

  Double_t dimuonDiTrkOne_mmpp, dimuonDiTrkTwo_mmpp, dimuonDiTrkThree_mmpp, dimuonDiTrkFour_mmpp, dimuonDiTrkFive_mmpp, dimuonDiTrkSix_mmpp;
  Double_t dimuonDiTrkOne_mmkk, dimuonDiTrkTwo_mmkk, dimuonDiTrkThree_mmkk, dimuonDiTrkFour_mmkk, dimuonDiTrkFive_mmkk, dimuonDiTrkSix_mmkk;
  Double_t dimuonDiTrkOne_mmpk, dimuonDiTrkTwo_mmpk, dimuonDiTrkThree_mmpk, dimuonDiTrkFour_mmpk, dimuonDiTrkFive_mmpk, dimuonDiTrkSix_mmpk;
  Double_t dimuonDiTrkOne_mmkp, dimuonDiTrkTwo_mmkp, dimuonDiTrkThree_mmkp, dimuonDiTrkFour_mmkp, dimuonDiTrkFive_mmkp, dimuonDiTrkSix_mmkp;

  Double_t triTrackOne_kkk, triTrackOne_kkp, triTrackOne_kpp, triTrackOne_ppp, triTrackOne_pt, triTrackOne_eta, triTrackOne_phi, triTrackOne_charge;
  Double_t triTrackFour_kkk, triTrackFour_kkp, triTrackFour_kpp, triTrackFour_ppp, triTrackFour_pt, triTrackFour_eta, triTrackFour_phi, triTrackFour_charge;
  Double_t triTrackThree_kkk, triTrackThree_kkp, triTrackThree_kpp, triTrackThree_ppp, triTrackThree_pt, triTrackThree_eta, triTrackThree_phi, triTrackThree_charge;
  Double_t triTrackTwo_kkk, triTrackTwo_kkp, triTrackTwo_kpp, triTrackTwo_ppp, triTrackTwo_pt, triTrackTwo_eta, triTrackTwo_phi, triTrackTwo_charge;

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

  Int_t track_KN_nvsh, track_KN_nvph;

  Int_t mumukk_rf_bindx;

  Int_t noSixCandidates;

  Bool_t isBestCandidate;


  //MC
  TLorentzVector gen_jpsi_p4, gen_phi_p4, gen_five_p4, gen_mumukk_p4, gen_fourthtrackp4,gen_six_p4;
  TLorentzVector gen_lowMuon_p4, gen_highMuon_p4, gen_highKaon_p4, gen_lowKaon_p4, gen_thirdtrack_p4;

  Double_t gen_five_prompt, gen_five_pt, gen_five_p, gen_five_eta, gen_five_phi;
  Double_t gen_five_pdg, gen_phi_pdg, gen_jpsi_pdg, gen_mumukk_pdg, gen_six_pdg;
  Double_t gen_lowMuon_pdg, gen_highMuon_pdg, gen_highKaon_pdg, gen_lowKaon_pdg, gen_thirdtrack_pdg;
  Double_t gen_lowMuon_mompdg, gen_highMuon_mompdg, gen_highKaon_mompdg, gen_lowKaon_mompdg, gen_thirdtrack_mompdg;
  Double_t gen_lowMuon_status, gen_highMuon_status, gen_highKaon_status, gen_lowKaon_status, gen_thirdtrack_status;
  Double_t gen_mumukk_prompt, gen_phi_prompt, gen_jpsi_prompt;
  Double_t gen_mumukk_pt, gen_mumukk_p, gen_mumukk_eta;
  Double_t gen_phi_pt, gen_phi_p, gen_phi_eta;
  Double_t gen_jpsi_pt, gen_jpsi_p, gen_jpsi_eta;
  Double_t gen_lowMuon_phi, gen_highMuon_phi, gen_highKaon_phi, gen_lowKaon_phi, gen_thirdtrack_phi;
  Double_t gen_mumukk_phi, gen_phi_phi, gen_jpsi_phi;
  Double_t gen_lowMuon_p, gen_highMuon_p, gen_highKaon_p, gen_lowKaon_p, gen_thirdtrack_p;
  Double_t gen_lowMuon_pt, gen_highMuon_pt, gen_highKaon_pt, gen_lowKaon_pt, gen_thirdtrack_pt;
  Double_t gen_lowMuon_eta, gen_highMuon_eta, gen_highKaon_eta, gen_lowKaon_eta, gen_thirdtrack_eta;

  Double_t gen_fourthtrackpdg, gen_fourthtrackmompdg, gen_fourthtrackstatus, gen_fourthtrackphi;
  Double_t gen_fourthtrackp, gen_fourthtrackpt, gen_fourthtracketa;

  Double_t six_m, six_m_ref, six_mass_ppkk, six_mass_pkpk;
  Double_t six_mass_pppp, six_mass_kpkp, six_mass_kppk, six_mass_kkkk;
  Double_t six_pt, six_eta, six_phi, six_p;
  Double_t six_cosAlpha, six_ctauPV, six_ctauErrPV;
  Double_t six_cosAlphaCA, six_ctauPVCA, six_ctauErrPVCA;
  Double_t six_cosAlphaDZ, six_ctauPVDZ, six_ctauErrPVDZ;
  Double_t six_cosAlphaBS, six_ctauPVBS, six_ctauErrPVBS;
  Double_t six_vProb, six_nDof, six_vChi2;
  Double_t six_vx, six_vy, six_vz, six_charge;

  Double_t swapped, sameSign, sameSign_mmtt;

  TTree* sixtracks_tree;
  edm::EDGetTokenT< std::vector <reco::GenParticle> > genCands_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;
};

UInt_t SixTracksRootupler::isTriggerMatched(pat::CompositeCandidate *mumu_cand) {
  const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(mumu_cand->daughter("muon1"));
  const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(mumu_cand->daughter("muon2"));
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
SixTracksRootupler::SixTracksRootupler(const edm::ParameterSet& iConfig):
        SixTracksCollection_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("SixTracksCand"))),
        TriggerResults_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
        thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PrimaryVertex"))),
        IsMC_(iConfig.getParameter<bool>("IsMC")),
        OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
        HLTs_(iConfig.getParameter<std::vector<std::string>>("HLTs")),
        HLTFilters_(iConfig.getParameter<std::vector<std::string>>("Filters")),
        TreeName_(iConfig.getParameter<std::string>("TreeName"))
{
	      edm::Service<TFileService> fs;
        sixtracks_tree = fs->make<TTree>(TreeName_.data(),"Tree of DiMuon and Four Tracks");

        sixtracks_tree->Branch("run",                &run,                "run/I");
        sixtracks_tree->Branch("event",              &event,              "event/I");
        sixtracks_tree->Branch("lumi",              &lumi,              "lumi/I");
        sixtracks_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");
        sixtracks_tree->Branch("trigger",            &trigger,            "trigger/I");

        if(!OnlyGen_)
        {
        sixtracks_tree->Branch("noSixCandidates",      &noSixCandidates,      "noSixCandidates/I");

        //sixtracks_tree->Branch("mumukk_id",      &mumukk_id,      "mumukk_id/I");
        sixtracks_tree->Branch("five_id",      &five_id,      "five_id/I");
        sixtracks_tree->Branch("dimuon_id",      &dimuon_id,      "dimuon_id/I");

        sixtracks_tree->Branch("p1_id",      &p1_id,      "p1_id/I");
        sixtracks_tree->Branch("m1_id",      &m1_id,      "m1_id/I");
        sixtracks_tree->Branch("p2_id",      &p2_id,      "p2_id/I");
        sixtracks_tree->Branch("m2_id",      &m2_id,      "m2_id/I");

        sixtracks_tree->Branch("swapped",      &swapped,      "swapped/D");
        sixtracks_tree->Branch("sameSign",      &sameSign,      "sameSign/D");
        sixtracks_tree->Branch("sameSign_mmtt",      &sameSign_mmtt,      "sameSign_mmtt/D");

        sixtracks_tree->Branch("six_p4",         "TLorentzVector", &six_p4);
        sixtracks_tree->Branch("five_p4",        "TLorentzVector", &five_p4);
        sixtracks_tree->Branch("mumukk_p4", "TLorentzVector", &mumukk_p4);
        sixtracks_tree->Branch("ditrack_p4",     "TLorentzVector", &ditrack_p4);
        sixtracks_tree->Branch("dimuon_p4",      "TLorentzVector", &dimuon_p4);

        sixtracks_tree->Branch("lowMuon_p4",    "TLorentzVector", &lowMuon_p4);
        sixtracks_tree->Branch("highMuon_p4",   "TLorentzVector", &highMuon_p4);

        sixtracks_tree->Branch("highKaon_p4",   "TLorentzVector", &highKaon_p4);
        sixtracks_tree->Branch("lowKaon_p4",    "TLorentzVector", &lowKaon_p4);
        sixtracks_tree->Branch("thirdKaon_p4",  "TLorentzVector", &thirdKaon_p4);
        sixtracks_tree->Branch("fourthKaon_p4",  "TLorentzVector", &fourthKaon_p4);

        sixtracks_tree->Branch("highPion_p4",   "TLorentzVector", &highPion_p4);
        sixtracks_tree->Branch("lowPion_p4",    "TLorentzVector", &lowPion_p4);
        sixtracks_tree->Branch("thirdPion_p4",  "TLorentzVector", &thirdPion_p4);
        sixtracks_tree->Branch("fourthPion_p4",  "TLorentzVector", &fourthPion_p4);

        sixtracks_tree->Branch("highProton_p4",   "TLorentzVector", &highProton_p4);
        sixtracks_tree->Branch("lowProton_p4",    "TLorentzVector", &lowProton_p4);
        sixtracks_tree->Branch("thirdProton_p4",  "TLorentzVector", &thirdProton_p4);
        sixtracks_tree->Branch("fourthProton_p4",  "TLorentzVector", &fourthProton_p4);

        sixtracks_tree->Branch("mumukk_m",       &mumukk_m,        "mumukk_m/D"); // Original DiMuonDiTrack (J/Psi Phi)
        sixtracks_tree->Branch("mumukk_pt",      &mumukk_pt,       "mumukk_pt/D");
        sixtracks_tree->Branch("mumukk_eta",     &mumukk_eta,      "mumukk_eta/D");
        sixtracks_tree->Branch("mumukk_phi",     &mumukk_phi,      "mumukk_phi/D");
        sixtracks_tree->Branch("mumukk_p",       &mumukk_p,      "mumukk_p/D");

        sixtracks_tree->Branch("dimuon_m",      &dimuon_m,     "dimuon_m/D");
        sixtracks_tree->Branch("dimuon_pt",     &dimuon_pt,    "dimuon_pt/D");
        sixtracks_tree->Branch("dimuon_eta",    &dimuon_pt,    "dimuon_eta/D");
        sixtracks_tree->Branch("dimuon_phi",    &dimuon_phi,   "dimuon_phi/D");
        sixtracks_tree->Branch("dimuon_p",      &dimuon_p,     "dimuon_p/D");

        sixtracks_tree->Branch("highTrackMatch",      &highTrackMatch,    "highTrackMatch/D");
        sixtracks_tree->Branch("lowTrackMatch",       &lowTrackMatch,     "lowTrackMatch/D");
        sixtracks_tree->Branch("lowMuonMatch",        &lowMuonMatch,      "lowMuonMatch/D");
        sixtracks_tree->Branch("highMuonMatch",       &highMuonMatch,     "highMuonMatch/D");
        sixtracks_tree->Branch("thirdTrackMatch",     &thirdTrackMatch,   "thirdTrackMatch/D");
        sixtracks_tree->Branch("fourthTrackMatch",    &fourthTrackMatch,  "fourthTrackMatch/D");

        //Di Tracks
        sixtracks_tree->Branch("ditrack_m",      &ditrack_m,     "ditrack_m/D"); //the original ditrack (supposed to be Phi->KK)

        sixtracks_tree->Branch("diTrackOne_pt",     &diTrackOne_pt,    "diTrackOne_pt/D");
        sixtracks_tree->Branch("diTrackOne_eta",    &diTrackOne_eta,    "diTrackOne_eta/D");
        sixtracks_tree->Branch("diTrackOne_phi",    &diTrackOne_phi,   "diTrackOne_phi/D");
        sixtracks_tree->Branch("diTrackOne_p",      &diTrackOne_p,   "diTrackOne_p/D");

        sixtracks_tree->Branch("diTrackTwo_pt",     &diTrackTwo_pt,    "diTrackTwo_pt/D");
        sixtracks_tree->Branch("diTrackTwo_eta",    &diTrackTwo_eta,    "diTrackTwo_eta/D");
        sixtracks_tree->Branch("diTrackTwo_phi",    &diTrackTwo_phi,   "diTrackTwo_phi/D");
        sixtracks_tree->Branch("diTrackTwo_p",      &diTrackTwo_p,   "diTrackTwo_p/D");

        sixtracks_tree->Branch("diTrackThree_pt",     &diTrackThree_pt,    "diTrackThree_pt/D");
        sixtracks_tree->Branch("diTrackThree_eta",    &diTrackThree_eta,    "diTrackThree_eta/D");
        sixtracks_tree->Branch("diTrackThree_phi",    &diTrackThree_phi,   "diTrackThree_phi/D");
        sixtracks_tree->Branch("diTrackThree_p",      &diTrackThree_p,   "diTrackThree_p/D");

        sixtracks_tree->Branch("diTrackFour_pt",     &diTrackFour_pt,    "diTrackFour_pt/D");
        sixtracks_tree->Branch("diTrackFour_eta",    &diTrackFour_eta,    "diTrackFour_eta/D");
        sixtracks_tree->Branch("diTrackFour_phi",    &diTrackFour_phi,   "diTrackFour_phi/D");
        sixtracks_tree->Branch("diTrackFour_p",      &diTrackFour_p,   "diTrackFour_p/D");

        sixtracks_tree->Branch("diTrackFive_pt",     &diTrackFive_pt,    "diTrackFive_pt/D");
        sixtracks_tree->Branch("diTrackFive_eta",    &diTrackFive_eta,    "diTrackFive_eta/D");
        sixtracks_tree->Branch("diTrackFive_phi",    &diTrackFive_phi,   "diTrackFive_phi/D");
        sixtracks_tree->Branch("diTrackFive_p",      &diTrackFive_p,   "diTrackFive_p/D");

        sixtracks_tree->Branch("diTrackSix_pt",     &diTrackSix_pt,    "diTrackSix_pt/D");
        sixtracks_tree->Branch("diTrackSix_eta",    &diTrackSix_eta,    "diTrackSix_eta/D");
        sixtracks_tree->Branch("diTrackSix_phi",    &diTrackSix_phi,   "diTrackSix_phi/D");
        sixtracks_tree->Branch("diTrackSix_p",      &diTrackSix_p,   "diTrackSix_p/D");

        // sixtracks_tree->Branch("dimuonDiTrkOne_mmpp",      &dimuonDiTrkOne_mmpp,     "dimuonDiTrkOne_mmpp/D");
        // sixtracks_tree->Branch("dimuonDiTrkTwo_mmpp",      &dimuonDiTrkTwo_mmpp,     "dimuonDiTrkTwo_mmpp/D");
        // sixtracks_tree->Branch("dimuonDiTrkThree_mmpp",    &dimuonDiTrkThree_mmpp,   "dimuonDiTrkThree_mmpp/D");
        // sixtracks_tree->Branch("dimuonDiTrkFour_mmpp",     &dimuonDiTrkFour_mmpp,    "dimuonDiTrkFour_mmpp/D");
        // sixtracks_tree->Branch("dimuonDiTrkFive_mmpp",     &dimuonDiTrkFive_mmpp,    "dimuonDiTrkFive_mmpp/D");
        // sixtracks_tree->Branch("dimuonDiTrkSix_mmpp",     &dimuonDiTrkSix_mmpp,    "dimuonDiTrkSix_mmpp/D");
        //
        // sixtracks_tree->Branch("dimuonDiTrkOne_mmkk",      &dimuonDiTrkOne_mmkk,    "dimuonDiTrkOne_mmkk/D");
        // sixtracks_tree->Branch("dimuonDiTrkTwo_mmkk",      &dimuonDiTrkTwo_mmkk,    "dimuonDiTrkTwo_mmkk/D");
        // sixtracks_tree->Branch("dimuonDiTrkThree_mmkk",    &dimuonDiTrkThree_mmkk,  "dimuonDiTrkThree_mmkk/D");
        // sixtracks_tree->Branch("dimuonDiTrkFour_mmkk",     &dimuonDiTrkFour_mmkk,   "dimuonDiTrkFour_mmkk/D");
        // sixtracks_tree->Branch("dimuonDiTrkFive_mmkk",     &dimuonDiTrkFive_mmkk,   "dimuonDiTrkFive_mmkk/D");
        // sixtracks_tree->Branch("dimuonDiTrkSix_mmkk",     &dimuonDiTrkSix_mmkk,   "dimuonDiTrkSix_mmkk/D");
        //
        // sixtracks_tree->Branch("dimuonDiTrkOne_mmpk",      &dimuonDiTrkOne_mmpk,    "dimuonDiTrkOne_mmpk/D");
        // sixtracks_tree->Branch("dimuonDiTrkTwo_mmpk",      &dimuonDiTrkTwo_mmpk,    "dimuonDiTrkTwo_mmpk/D");
        // sixtracks_tree->Branch("dimuonDiTrkThree_mmpk",    &dimuonDiTrkThree_mmpk,  "dimuonDiTrkThree_mmpk/D");
        // sixtracks_tree->Branch("dimuonDiTrkFour_mmpk",     &dimuonDiTrkFour_mmpk,   "dimuonDiTrkFour_mmpk/D");
        // sixtracks_tree->Branch("dimuonDiTrkFive_mmpk",     &dimuonDiTrkFive_mmpk,   "dimuonDiTrkFive_mmpk/D");
        // sixtracks_tree->Branch("dimuonDiTrkSix_mmpk",     &dimuonDiTrkSix_mmpk,   "dimuonDiTrkSix_mmpk/D");
        //
        // sixtracks_tree->Branch("dimuonDiTrkOne_mmkp",      &dimuonDiTrkOne_mmkp,    "dimuonDiTrkOne_mmkp/D");
        // sixtracks_tree->Branch("dimuonDiTrkTwo_mmkp",      &dimuonDiTrkTwo_mmkp,    "dimuonDiTrkTwo_mmkp/D");
        // sixtracks_tree->Branch("dimuonDiTrkThree_mmkp",    &dimuonDiTrkThree_mmkp,  "dimuonDiTrkThree_mmkp/D");
        // sixtracks_tree->Branch("dimuonDiTrkFour_mmkp",     &dimuonDiTrkFour_mmkp,   "dimuonDiTrkFour_mmkp/D");
        // sixtracks_tree->Branch("dimuonDiTrkFive_mmkp",     &dimuonDiTrkFive_mmkp,   "dimuonDiTrkFive_mmkp/D");
        // sixtracks_tree->Branch("dimuonDiTrkSix_mmkp",     &dimuonDiTrkSix_mmkp,   "dimuonDiTrkSix_mmkp/D");

        sixtracks_tree->Branch("diTrackOne_kk",     &diTrackOne_kk,    "diTrackOne_kk/D");
        sixtracks_tree->Branch("diTrackTwo_kk",     &diTrackTwo_kk,    "diTrackTwo_kk/D");
        sixtracks_tree->Branch("diTrackThree_kk",     &diTrackThree_kk,    "diTrackThree_kk/D");
        sixtracks_tree->Branch("diTrackFour_kk",     &diTrackFour_kk,    "diTrackFour_kk/D");
        sixtracks_tree->Branch("diTrackFive_kk",     &diTrackFive_kk,    "diTrackFive_kk/D");
        sixtracks_tree->Branch("diTrackSix_kk",     &diTrackSix_kk,    "diTrackSix_kk/D");

        sixtracks_tree->Branch("diTrackOne_pp",     &diTrackOne_pp,    "diTrackOne_pp/D");
        sixtracks_tree->Branch("diTrackTwo_pp",     &diTrackTwo_pp,    "diTrackTwo_pp/D");
        sixtracks_tree->Branch("diTrackThree_pp",     &diTrackThree_pp,    "diTrackThree_pp/D");
        sixtracks_tree->Branch("diTrackFour_pp",     &diTrackFour_pp,    "diTrackFour_pp/D");
        sixtracks_tree->Branch("diTrackFive_pp",     &diTrackFive_pp,    "diTrackFive_pp/D");
        sixtracks_tree->Branch("diTrackSix_pp",     &diTrackSix_pp,    "diTrackSix_pp/D");

        sixtracks_tree->Branch("diTrackOne_pk",     &diTrackOne_pk,    "diTrackOne_pk/D");
        sixtracks_tree->Branch("diTrackTwo_pk",     &diTrackTwo_pk,    "diTrackTwo_pk/D");
        sixtracks_tree->Branch("diTrackThree_pk",     &diTrackThree_pk,    "diTrackThree_pk/D");
        sixtracks_tree->Branch("diTrackFour_pk",     &diTrackFour_pk,    "diTrackFour_pk/D");
        sixtracks_tree->Branch("diTrackFive_pk",     &diTrackFive_pk,    "diTrackFive_pk/D");
        sixtracks_tree->Branch("diTrackSix_pk",     &diTrackSix_pk,    "diTrackSix_pk/D");

        sixtracks_tree->Branch("diTrackOne_kp",     &diTrackOne_kp,    "diTrackOne_kp/D");
        sixtracks_tree->Branch("diTrackTwo_kp",     &diTrackTwo_kp,    "diTrackTwo_kp/D");
        sixtracks_tree->Branch("diTrackThree_kp",     &diTrackThree_kp,    "diTrackThree_kp/D");
        sixtracks_tree->Branch("diTrackFour_kp",     &diTrackFour_kp,    "diTrackFour_kp/D");
        sixtracks_tree->Branch("diTrackFive_kp",     &diTrackFive_kp,    "diTrackFive_kp/D");
        sixtracks_tree->Branch("diTrackSix_kp",     &diTrackSix_kp,    "diTrackSix_kp/D");

        //The kinematic doesn't change, only mass
        sixtracks_tree->Branch("highMuon_pt",         &highMuon_pt,         "highMuon_pt/D");
        sixtracks_tree->Branch("highMuon_eta",        &highMuon_eta,        "highMuon_eta/D");
        sixtracks_tree->Branch("highMuon_phi",        &highMuon_phi,        "highMuon_phi/D");
        sixtracks_tree->Branch("highMuon_charge",     &highMuon_charge,     "highMuon_charge/D");
        sixtracks_tree->Branch("highMuon_dz",         &highMuon_dz,         "highMuon_dz/D");
        sixtracks_tree->Branch("highMuon_dxy",        &highMuon_dxy,        "highMuon_dxy/D");

        sixtracks_tree->Branch("lowMuon_pt",         &lowMuon_pt,         "lowMuon_pt/D");
        sixtracks_tree->Branch("lowMuon_eta",        &lowMuon_eta,        "lowMuon_eta/D");
        sixtracks_tree->Branch("lowMuon_phi",        &lowMuon_phi,        "lowMuon_phi/D");
        sixtracks_tree->Branch("lowMuon_charge",     &lowMuon_charge,     "lowMuon_charge/D");
        sixtracks_tree->Branch("lowMuon_dz",         &lowMuon_dz,         "lowMuon_dz/D");
        sixtracks_tree->Branch("lowMuon_dxy",        &lowMuon_dxy,        "lowMuon_dxy/D");

        //The kinematic doesn't change, only mass
        sixtracks_tree->Branch("highTrack_pt",         &highTrack_pt,         "highTrack_pt/D");
        sixtracks_tree->Branch("highTrack_eta",        &highTrack_eta,        "highTrack_eta/D");
        sixtracks_tree->Branch("highTrack_phi",        &highTrack_phi,        "highTrack_phi/D");
        sixtracks_tree->Branch("highTrack_charge",     &highTrack_charge,     "highTrack_charge/D");
        sixtracks_tree->Branch("highTrack_dz",         &highTrack_dz,         "highTrack_dz/D");
        sixtracks_tree->Branch("highTrack_dxy",        &highTrack_dxy,        "highTrack_dxy/D");

        sixtracks_tree->Branch("lowTrack_pt",           &lowTrack_pt,        "lowTrack_pt/D");
        sixtracks_tree->Branch("lowTrack_eta",          &lowTrack_eta,       "lowTrack_eta/D");
        sixtracks_tree->Branch("lowTrack_phi",          &lowTrack_phi,       "lowTrack_phi/D");
        sixtracks_tree->Branch("lowTrack_charge",       &lowTrack_charge,    "lowTrack_charge/D");
        sixtracks_tree->Branch("lowTrack_dz",           &lowTrack_dz,          "lowTrack_dz/D");
        sixtracks_tree->Branch("lowTrack_dxy",          &lowTrack_dxy,        "lowTrack_dxy/D");

        sixtracks_tree->Branch("thirdTrack_pt",         &thirdTrack_pt,        "thirdTrack_pt/D");
        sixtracks_tree->Branch("thirdTrack_eta",        &thirdTrack_eta,        "thirdTrack_eta/D");
        sixtracks_tree->Branch("thirdTrack_phi",        &thirdTrack_phi,        "thirdTrack_phi/D");
        sixtracks_tree->Branch("thirdTrack_charge",     &thirdTrack_charge,  "thirdTrack_charge/D");
        sixtracks_tree->Branch("thirdTrack_dz",         &thirdTrack_dz,          "thirdTrack_dz/D");
        sixtracks_tree->Branch("thirdTrack_dxy",        &thirdTrack_dxy,        "thirdTrack_dxy/D");

        sixtracks_tree->Branch("fourthTrack_pt",         &fourthTrack_pt,        "fourthTrack_pt/D");
        sixtracks_tree->Branch("fourthTrack_eta",        &fourthTrack_eta,        "fourthTrack_eta/D");
        sixtracks_tree->Branch("fourthTrack_phi",        &fourthTrack_phi,        "fourthTrack_phi/D");
        sixtracks_tree->Branch("fourthTrack_charge",     &fourthTrack_charge,  "fourthTrack_charge/D");
        sixtracks_tree->Branch("fourthTrack_dz",         &fourthTrack_dz,          "fourthTrack_dz/D");
        sixtracks_tree->Branch("fourthTrack_dxy",        &fourthTrack_dxy,        "fourthTrack_dxy/D");

        //J/Psi TrTr system
        sixtracks_tree->Branch("dimuonDiTrkOne_pt",         &dimuonDiTrkOne_pt,         "dimuonDiTrkOne_pt/D");
        sixtracks_tree->Branch("dimuonDiTrkOne_eta",        &dimuonDiTrkOne_eta,        "dimuonDiTrkOne_eta/D");
        sixtracks_tree->Branch("dimuonDiTrkOne_phi",        &dimuonDiTrkOne_phi,        "dimuonDiTrkOne_phi/D");
        sixtracks_tree->Branch("dimuonDiTrkOne_charge",     &dimuonDiTrkOne_charge,     "dimuonDiTrkOne_charge/D");
        sixtracks_tree->Branch("dimuonDiTrkOne_p",          &dimuonDiTrkOne_p,          "dimuonDiTrkOne_p/D");

        sixtracks_tree->Branch("dimuonDiTrkTwo_pt",         &dimuonDiTrkTwo_pt,         "dimuonDiTrkTwo_pt/D");
        sixtracks_tree->Branch("dimuonDiTrkTwo_eta",        &dimuonDiTrkTwo_eta,        "dimuonDiTrkTwo_eta/D");
        sixtracks_tree->Branch("dimuonDiTrkTwo_phi",        &dimuonDiTrkTwo_phi,        "dimuonDiTrkTwo_phi/D");
        sixtracks_tree->Branch("dimuonDiTrkTwo_charge",     &dimuonDiTrkTwo_charge,     "dimuonDiTrkTwo_charge/D");
        sixtracks_tree->Branch("dimuonDiTrkTwo_p",         &dimuonDiTrkTwo_p,         "dimuonDiTrkTwo_p/D");

        sixtracks_tree->Branch("dimuonDiTrkThree_pt",       &dimuonDiTrkThree_pt,         "dimuonDiTrkThree_pt/D");
        sixtracks_tree->Branch("dimuonDiTrkThree_eta",      &dimuonDiTrkThree_eta,        "dimuonDiTrkThree_eta/D");
        sixtracks_tree->Branch("dimuonDiTrkThree_phi",      &dimuonDiTrkThree_phi,        "dimuonDiTrkThree_phi/D");
        sixtracks_tree->Branch("dimuonDiTrkThree_charge",   &dimuonDiTrkThree_charge,     "dimuonDiTrkThree_charge/D");
        sixtracks_tree->Branch("dimuonDiTrkThree_p",        &dimuonDiTrkThree_p,          "dimuonDiTrkThree_p/D");

        sixtracks_tree->Branch("dimuonDiTrkFour_pt",         &dimuonDiTrkFour_pt,          "dimuonDiTrkFour_pt/D");
        sixtracks_tree->Branch("dimuonDiTrkFour_eta",        &dimuonDiTrkFour_eta,        "dimuonDiTrkFour_eta/D");
        sixtracks_tree->Branch("dimuonDiTrkFour_phi",        &dimuonDiTrkFour_phi,        "dimuonDiTrkFour_phi/D");
        sixtracks_tree->Branch("dimuonDiTrkFour_charge",     &dimuonDiTrkFour_charge,        "dimuonDiTrkFour_charge/D");
        sixtracks_tree->Branch("dimuonDiTrkFour_p",          &dimuonDiTrkFour_p,          "dimuonDiTrkFour_p/D");

        sixtracks_tree->Branch("dimuonDiTrkFive_pt",         &dimuonDiTrkFive_pt,         "dimuonDiTrkFive_pt/D");
        sixtracks_tree->Branch("dimuonDiTrkFive_eta",        &dimuonDiTrkFive_eta,        "dimuonDiTrkFive_eta/D");
        sixtracks_tree->Branch("dimuonDiTrkFive_phi",        &dimuonDiTrkFive_phi,        "dimuonDiTrkFive_phi/D");
        sixtracks_tree->Branch("dimuonDiTrkFive_charge",     &dimuonDiTrkFive_charge,     "dimuonDiTrkFive_charge/D");
        sixtracks_tree->Branch("dimuonDiTrkFive_p",         &dimuonDiTrkFive_p,         "dimuonDiTrkFive_p/D");

        sixtracks_tree->Branch("dimuonDiTrkSix_pt",       &dimuonDiTrkSix_pt,         "dimuonDiTrkSix_pt/D");
        sixtracks_tree->Branch("dimuonDiTrkSix_eta",      &dimuonDiTrkSix_eta,        "dimuonDiTrkSix_eta/D");
        sixtracks_tree->Branch("dimuonDiTrkSix_phi",      &dimuonDiTrkSix_phi,        "dimuonDiTrkSix_phi/D");
        sixtracks_tree->Branch("dimuonDiTrkSix_charge",   &dimuonDiTrkSix_charge,     "dimuonDiTrkSix_charge/D");
        sixtracks_tree->Branch("dimuonDiTrkSix_p",        &dimuonDiTrkSix_p,          "dimuonDiTrkSix_p/D");

        sixtracks_tree->Branch("dimuon_vProb",        &dimuon_vProb,        "dimuon_vProb/D");
        sixtracks_tree->Branch("dimuon_vChi2",       &dimuon_vChi2,        "dimuon_vChi2/D");
        sixtracks_tree->Branch("dimuon_DCA",          &dimuon_DCA,          "dimuon_DCA/D");
        sixtracks_tree->Branch("dimuon_ctauPV",       &dimuon_ctauPV,       "dimuon_ctauPV/D");
        sixtracks_tree->Branch("dimuon_ctauErrPV",    &dimuon_ctauErrPV,    "dimuon_ctauErrPV/D");
        sixtracks_tree->Branch("dimuon_cosAlpha",     &dimuon_cosAlpha,     "dimuon_cosAlpha/D");
        // sixtracks_tree->Branch("dimuon_triggerMatch", &dimuon_triggerMatch, "dimuon_triggerMatch/I");

        //TriTrack system
        sixtracks_tree->Branch("triTrackOne_kkk",         &triTrackOne_kkk,         "triTrackOne_kkk/D");
        sixtracks_tree->Branch("triTrackTwo_kkk",         &triTrackTwo_kkk,         "triTrackTwo_kkk/D");
        sixtracks_tree->Branch("triTrackThree_kkk",         &triTrackThree_kkk,         "triTrackThree_kkk/D");
        sixtracks_tree->Branch("triTrackFour_kkk",         &triTrackFour_kkk,         "triTrackFour_kkk/D");

        sixtracks_tree->Branch("triTrackOne_kkp",         &triTrackOne_kkp,         "triTrackOne_kkp/D");
        sixtracks_tree->Branch("triTrackTwo_kkp",         &triTrackTwo_kkp,         "triTrackTwo_kkp/D");
        sixtracks_tree->Branch("triTrackThree_kkp",         &triTrackThree_kkp,         "triTrackThree_kkp/D");
        sixtracks_tree->Branch("triTrackThree_kkp",         &triTrackThree_kkp,         "triTrackThree_kkp/D");

        sixtracks_tree->Branch("triTrackOne_kpp",         &triTrackOne_kpp,         "triTrackOne_kpp/D");
        sixtracks_tree->Branch("triTrackTwo_kpp",         &triTrackTwo_kpp,         "triTrackTwo_kpp/D");
        sixtracks_tree->Branch("triTrackThree_kpp",         &triTrackThree_kpp,         "triTrackThree_kpp/D");
        sixtracks_tree->Branch("triTrackThree_kpp",         &triTrackThree_kpp,         "triTrackThree_kpp/D");

        sixtracks_tree->Branch("triTrackOne_ppp",         &triTrackOne_ppp,         "triTrackOne_ppp/D");
        sixtracks_tree->Branch("triTrackTwo_ppp",         &triTrackTwo_ppp,         "triTrackTwo_ppp/D");
        sixtracks_tree->Branch("triTrackThree_ppp",         &triTrackThree_ppp,         "triTrackThree_ppp/D");
        sixtracks_tree->Branch("triTrackThree_ppp",         &triTrackThree_ppp,         "triTrackThree_ppp/D");

        sixtracks_tree->Branch("triTrackOne_pt",         &triTrackOne_pt,         "triTrackOne_pt/D");
        sixtracks_tree->Branch("triTrackOne_eta",        &triTrackOne_eta,        "triTrackOne_eta/D");
        sixtracks_tree->Branch("triTrackOne_phi",        &triTrackOne_phi,        "triTrackOne_phi/D");
        sixtracks_tree->Branch("triTrackOne_charge",     &triTrackOne_charge,     "triTrackOne_charge/D");

        sixtracks_tree->Branch("triTrackTwo_pt",         &triTrackTwo_pt,         "triTrackTwo_pt/D");
        sixtracks_tree->Branch("triTrackTwo_eta",        &triTrackTwo_eta,        "triTrackTwo_eta/D");
        sixtracks_tree->Branch("triTrackTwo_phi",        &triTrackTwo_phi,        "triTrackTwo_phi/D");
        sixtracks_tree->Branch("triTrackTwo_charge",     &triTrackTwo_charge,     "triTrackTwo_charge/D");

        sixtracks_tree->Branch("triTrackThree_pt",         &triTrackThree_pt,         "triTrackThree_pt/D");
        sixtracks_tree->Branch("triTrackThree_eta",        &triTrackThree_eta,        "triTrackThree_eta/D");
        sixtracks_tree->Branch("triTrackThree_phi",        &triTrackThree_phi,        "triTrackThree_phi/D");
        sixtracks_tree->Branch("triTrackThree_charge",     &triTrackThree_charge,     "triTrackThree_charge/D");

        sixtracks_tree->Branch("triTrackFour_pt",         &triTrackFour_pt,         "triTrackFour_pt/D");
        sixtracks_tree->Branch("triTrackFour_eta",        &triTrackFour_eta,        "triTrackFour_eta/D");
        sixtracks_tree->Branch("triTrackFour_phi",        &triTrackFour_phi,        "triTrackFour_phi/D");
        sixtracks_tree->Branch("triTrackFour_charge",     &triTrackFour_charge,     "triTrackFour_charge/D");

        //FourTrack System
        sixtracks_tree->Branch("mumukk_vProb",      &mumukk_vProb,        "mumukk_vProb/D");
        sixtracks_tree->Branch("mumukk_vChi2",      &mumukk_vChi2,        "mumukk_vChi2/D");
        sixtracks_tree->Branch("mumukk_nDof",       &mumukk_nDof,         "mumukk_nDof/D");
        sixtracks_tree->Branch("mumukk_charge",     &mumukk_charge,       "mumukk_charge/I");

        sixtracks_tree->Branch("mumukk_cosAlpha",      &mumukk_cosAlpha,        "mumukk_cosAlpha/D");
        sixtracks_tree->Branch("mumukk_ctauPV",      &mumukk_ctauPV,        "mumukk_ctauPV/D");
        sixtracks_tree->Branch("mumukk_ctauErrPV",      &mumukk_ctauErrPV,        "mumukk_ctauErrPV/D");

        sixtracks_tree->Branch("mumukk_cosAlphaCA",      &mumukk_cosAlphaCA,        "mumukk_cosAlphaCA/D");
        sixtracks_tree->Branch("mumukk_ctauPVCA",      &mumukk_ctauPVCA,        "mumukk_ctauPVCA/D");
        sixtracks_tree->Branch("mumukk_ctauErrPVCA",      &mumukk_ctauErrPVCA,        "mumukk_ctauErrPVCA/D");

        sixtracks_tree->Branch("mumukk_cosAlphaDZ",      &mumukk_cosAlphaDZ,        "mumukk_cosAlphaDZ/D");
        sixtracks_tree->Branch("mumukk_ctauPVDZ",      &mumukk_ctauPVDZ,        "mumukk_ctauPVDZ/D");
        sixtracks_tree->Branch("mumukk_ctauErrPVDZ",      &mumukk_ctauErrPVDZ,        "mumukk_ctauErrPVDZ/D");

        sixtracks_tree->Branch("mumukk_cosAlphaBS",      &mumukk_cosAlphaBS,        "mumukk_cosAlphaBS/D");
        sixtracks_tree->Branch("mumukk_ctauPVBS",      &mumukk_ctauPVBS,        "mumukk_ctauPVBS/D");
        sixtracks_tree->Branch("mumukk_ctauErrPVBS",      &mumukk_ctauErrPVBS,        "mumukk_ctauErrPVBS/D");

        sixtracks_tree->Branch("mumukk_vx",     &mumukk_vx,      "mumukk_vx/D");
        sixtracks_tree->Branch("mumukk_vy",     &mumukk_vy,      "mumukk_vy/D");
        sixtracks_tree->Branch("mumukk_vz",     &mumukk_vz,      "mumukk_vz/D");

        //DCAs
        sixtracks_tree->Branch("dca_m1m2",      &dca_m1m2,        "dca_m1m2/D");
        sixtracks_tree->Branch("dca_m1t1",      &dca_m1t1,        "dca_m1t1/D");
        sixtracks_tree->Branch("dca_m1t2",      &dca_m1t2,        "dca_m1t2/D");
        sixtracks_tree->Branch("dca_m2t1",      &dca_m2t1,        "dca_m2t1/D");
        sixtracks_tree->Branch("dca_m2t2",      &dca_m2t2,        "dca_m2t2/D");
        sixtracks_tree->Branch("dca_t1t2",      &dca_t1t2,        "dca_t1t2/D");
        sixtracks_tree->Branch("dca_m1t3",      &dca_m1t3,        "dca_m1t3/D");
        sixtracks_tree->Branch("dca_m2t3",      &dca_m2t3,        "dca_m2t3/D");
        sixtracks_tree->Branch("dca_t1t3",      &dca_t1t3,        "dca_t1t3/D");
        sixtracks_tree->Branch("dca_t2t3",      &dca_t2t3,        "dca_t2t3/D");
        sixtracks_tree->Branch("dca_m1t4",      &dca_m1t4,        "dca_m1t4/D");
        sixtracks_tree->Branch("dca_m2t4",      &dca_m2t4,        "dca_m2t4/D");
        sixtracks_tree->Branch("dca_t1t4",      &dca_t1t4,        "dca_t1t4/D");
        sixtracks_tree->Branch("dca_t2t4",      &dca_t2t4,        "dca_t2t4/D");
        sixtracks_tree->Branch("dca_t3t4",      &dca_t3t4,        "dca_t3t4/D");

        sixtracks_tree->Branch("highTrackMuonDR",       &highTrackMuonDR,        "highTrackMuonDR/D");
        sixtracks_tree->Branch("highTrackMuonDP",       &highTrackMuonDP,        "highTrackMuonDP/D");
        sixtracks_tree->Branch("highTrackMuonDPt",      &highTrackMuonDPt,        "highTrackMuonDPt/D");

        sixtracks_tree->Branch("lowTrackMuonDR",        &lowTrackMuonDR,        "lowTrackMuonDR/D");
        sixtracks_tree->Branch("lowTrackMuonDP",        &lowTrackMuonDP,        "lowTrackMuonDP/D");
        sixtracks_tree->Branch("lowTrackMuonDPt",       &lowTrackMuonDPt,        "lowTrackMuonDPt/D");

        sixtracks_tree->Branch("thirdTrackMuonDR",      &thirdTrackMuonDR,        "thirdTrackMuonDR/D");
        sixtracks_tree->Branch("thirdTrackMuonDP",      &thirdTrackMuonDP,        "thirdTrackMuonDP/D");
        sixtracks_tree->Branch("thirdTrackMuonDPt",     &thirdTrackMuonDPt,        "thirdTrackMuonDPt/D");

        sixtracks_tree->Branch("fourthTrackMuonDR",     &fourthTrackMuonDR,        "fourthTrackMuonDR/D");
        sixtracks_tree->Branch("fourthTrackMuonDP",     &fourthTrackMuonDP,        "fourthTrackMuonDP/D");
        sixtracks_tree->Branch("fourthTrackMuonDPt",    &fourthTrackMuonDPt,        "fourthTrackMuonDPt/D");

        sixtracks_tree->Branch("highTrackFromPV",       &highTrackFromPV,        "highTrackFromPV/D");
        sixtracks_tree->Branch("lowTrackFromPV",        &lowTrackFromPV,        "lowTrackFromPV/D");
        sixtracks_tree->Branch("thirdTrackFromPV",      &thirdTrackFromPV,        "thirdTrackFromPV/D");
        sixtracks_tree->Branch("fourthTrackFromPV",     &fourthTrackFromPV,        "fourthTrackFromPV/D");

        sixtracks_tree->Branch("highTrackFromPVCA",     &highTrackFromPVCA,        "highTrackFromPVCA/D");
        sixtracks_tree->Branch("lowTrackFromPVCA",      &lowTrackFromPVCA,        "lowTrackFromPVCA/D");
        sixtracks_tree->Branch("thirdTrackFromPVCA",    &thirdTrackFromPVCA,      "thirdTrackFromPVCA/D");
        sixtracks_tree->Branch("fourthTrackFromPVCA",   &fourthTrackFromPVCA,     "fourthTrackFromPVCA/D");

        sixtracks_tree->Branch("highTrackFromPVDZ",     &highTrackFromPVDZ,        "highTrackFromPVDZ/D");
        sixtracks_tree->Branch("lowTrackFromPVDZ",      &lowTrackFromPVDZ,        "lowTrackFromPVDZ/D");
        sixtracks_tree->Branch("thirdTrackFromPVDZ",    &thirdTrackFromPVDZ,        "thirdTrackFromPVDZ/D");
        sixtracks_tree->Branch("fourthTrackFromPVDZ",   &fourthTrackFromPVDZ,        "fourthTrackFromPVDZ/D");

        // sixtracks_tree->Branch("highTrackFromPVBS",      &highTrackFromPVBS,        "highTrackFromPVBS/D");
        // sixtracks_tree->Branch("lowTrackFromPVBS",      &lowTrackFromPVBS,        "lowTrackFromPVBS/D");
        // sixtracks_tree->Branch("thirdTrackFromPVBS",      &thirdTrackFromPVBS,        "thirdTrackFromPVBS/D");
        // sixtracks_tree->Branch("fourthTrackFromPVBS",      &fourthTrackFromPVBS,        "fourthTrackFromPVBS/D");


        //Five Tracks System
        sixtracks_tree->Branch("five_m",             &five_m,               "five_m/D");
        sixtracks_tree->Branch("five_m_ref",         &five_m_ref,           "five_m_ref/D");
        sixtracks_tree->Branch("five_mass_ppk",      &five_mass_ppk,        "five_mass_ppk/D");
        sixtracks_tree->Branch("five_mass_kpp",      &five_mass_kpp,        "five_mass_kpp/D");
        sixtracks_tree->Branch("five_mass_pkp",      &five_mass_pkp,        "five_mass_pkp/D");
        sixtracks_tree->Branch("five_mass_ppp",      &five_mass_ppp,        "five_mass_ppp/D");

        sixtracks_tree->Branch("fiveOne_pt",      &fiveOne_pt,       "fiveOne_pt/D");
        sixtracks_tree->Branch("fiveOne_eta",     &fiveOne_eta,      "fiveOne_eta/D");
        sixtracks_tree->Branch("fiveOne_phi",     &fiveOne_phi,      "fiveOne_phi/D");
        sixtracks_tree->Branch("fiveOne_p",       &fiveOne_p,        "fiveOne_p/D");

        sixtracks_tree->Branch("fiveTwo_pt",      &fiveTwo_pt,       "fiveTwo_pt/D");
        sixtracks_tree->Branch("fiveTwo_eta",     &fiveTwo_eta,      "fiveTwo_eta/D");
        sixtracks_tree->Branch("fiveTwo_phi",     &fiveTwo_phi,      "fiveTwo_phi/D");
        sixtracks_tree->Branch("fiveTwo_p",       &fiveTwo_p,        "fiveTwo_p/D");

        sixtracks_tree->Branch("fiveThree_pt",      &fiveThree_pt,       "fiveThree_pt/D");
        sixtracks_tree->Branch("fiveThree_eta",     &fiveThree_eta,      "fiveThree_eta/D");
        sixtracks_tree->Branch("fiveThree_phi",     &fiveThree_phi,      "fiveThree_phi/D");
        sixtracks_tree->Branch("fiveThree_p",       &fiveThree_p,        "fiveThree_p/D");

        sixtracks_tree->Branch("fiveFour_pt",      &fiveFour_pt,       "fiveFour_pt/D");
        sixtracks_tree->Branch("fiveFour_eta",     &fiveFour_eta,      "fiveFour_eta/D");
        sixtracks_tree->Branch("fiveFour_phi",     &fiveFour_phi,      "fiveFour_phi/D");
        sixtracks_tree->Branch("fiveFour_p",       &fiveFour_p,        "fiveFour_p/D");

        sixtracks_tree->Branch("fiveFive_pt",      &fiveFive_pt,       "fiveFive_pt/D");
        sixtracks_tree->Branch("fiveFive_eta",     &fiveFive_eta,      "fiveFive_eta/D");
        sixtracks_tree->Branch("fiveFive_phi",     &fiveFive_phi,      "fiveFive_phi/D");
        sixtracks_tree->Branch("fiveFive_p",       &fiveFive_p,        "fiveFive_p/D");

        sixtracks_tree->Branch("five_cosAlpha",     &five_cosAlpha,        "five_cosAlpha/D");
        sixtracks_tree->Branch("five_ctauPV",       &five_ctauPV,        "five_ctauPV/D");
        sixtracks_tree->Branch("five_ctauErrPV",    &five_ctauErrPV,        "five_ctauErrPV/D");

        sixtracks_tree->Branch("five_cosAlphaCA",    &five_cosAlphaCA,        "five_cosAlphaCA/D");
        sixtracks_tree->Branch("five_ctauPVCA",      &five_ctauPVCA,        "five_ctauPVCA/D");
        sixtracks_tree->Branch("five_ctauErrPVCA",   &five_ctauErrPVCA,        "five_ctauErrPVCA/D");

        sixtracks_tree->Branch("five_cosAlphaDZ",    &five_cosAlphaDZ,        "five_cosAlphaDZ/D");
        sixtracks_tree->Branch("five_ctauPVDZ",      &five_ctauPVDZ,        "five_ctauPVDZ/D");
        sixtracks_tree->Branch("five_ctauErrPVDZ",   &five_ctauErrPVDZ,        "five_ctauErrPVDZ/D");

        sixtracks_tree->Branch("five_cosAlphaBS",    &five_cosAlphaBS,        "five_cosAlphaBS/D");
        sixtracks_tree->Branch("five_ctauPVBS",      &five_ctauPVBS,        "five_ctauPVBS/D");
        sixtracks_tree->Branch("five_ctauErrPVBS",   &five_ctauErrPVBS,        "five_ctauErrPVBS/D");

        sixtracks_tree->Branch("five_vProb",       &five_vProb,      "five_vProb/D");
        sixtracks_tree->Branch("five_nDof",        &five_nDof,       "five_nDof/D");
        sixtracks_tree->Branch("five_vChi2",       &five_vChi2,      "five_vChi2/D");

        sixtracks_tree->Branch("five_vx",     &five_vx,      "five_vx/D");
        sixtracks_tree->Branch("five_vy",     &five_vy,      "five_vy/D");
        sixtracks_tree->Branch("five_vz",     &five_vz,      "five_vz/D");

        sixtracks_tree->Branch("five_charge",     &five_charge,       "five_charge/I");

        //Vertex
        sixtracks_tree->Branch("bestPV_X",     &bestPV_X,      "bestPV_X/D");
        sixtracks_tree->Branch("bestPV_Y",     &bestPV_Y,      "bestPV_Y/D");
        sixtracks_tree->Branch("bestPV_Z",     &bestPV_Z,      "bestPV_Z/D");

        sixtracks_tree->Branch("cosAlphaPV_X",     &cosAlphaPV_X,      "cosAlphaPV_X/D");
        sixtracks_tree->Branch("cosAlphaPV_Y",     &cosAlphaPV_Y,      "cosAlphaPV_Y/D");
        sixtracks_tree->Branch("cosAlphaPV_Z",     &cosAlphaPV_Z,      "cosAlphaPV_Z/D");

        sixtracks_tree->Branch("bS_X",     &bS_X,      "bS_X/D");
        sixtracks_tree->Branch("bS_Y",     &bS_Y,      "bS_Y/D");
        sixtracks_tree->Branch("bS_Z",     &bS_Z,      "bS_Z/D");

        sixtracks_tree->Branch("zPV_X",     &zPV_X,      "zPV_X/D");
        sixtracks_tree->Branch("zPV_Y",     &zPV_Y,      "zPV_Y/D");
        sixtracks_tree->Branch("zPV_Z",     &zPV_Z,      "zPV_Z/D");



        //Muon flags
        sixtracks_tree->Branch("lowMuon_isTight",        &lowMuon_isTight,        "lowMuon_isTight/O");
        sixtracks_tree->Branch("lowMuon_isLoose",        &lowMuon_isLoose,        "lowMuon_isLoose/O");
        sixtracks_tree->Branch("lowMuon_isSoft",        &lowMuon_isSoft,        "lowMuon_isSoft/O");
        sixtracks_tree->Branch("lowMuon_isMedium",        &lowMuon_isMedium,        "lowMuon_isMedium/O");
        sixtracks_tree->Branch("lowMuon_isHighPt",        &lowMuon_isHighPt,        "lowMuon_isHighPt/O");

        sixtracks_tree->Branch("lowMuon_isTracker",        &lowMuon_isTracker,        "lowMuon_isTracker/O");
        sixtracks_tree->Branch("lowMuon_isGlobal",        &lowMuon_isGlobal,        "lowMuon_isGlobal/O");

        sixtracks_tree->Branch("lowMuon_NPixelHits",        &lowMuon_NPixelHits,        "lowMuon_NPixelHits/I");
        sixtracks_tree->Branch("lowMuon_NStripHits",        &lowMuon_NStripHits,        "lowMuon_NStripHits/I");
        sixtracks_tree->Branch("lowMuon_NTrackhits",        &lowMuon_NTrackhits,        "lowMuon_NTrackhits/I");
        sixtracks_tree->Branch("lowMuon_NBPixHits",        &lowMuon_NBPixHits,        "lowMuon_NBPixHits/I");

        sixtracks_tree->Branch("lowMuon_NPixLayers",        &lowMuon_NPixLayers,        "lowMuon_NPixLayers/I");
        sixtracks_tree->Branch("lowMuon_NTraLayers",        &lowMuon_NTraLayers,        "lowMuon_NTraLayers/I");
        sixtracks_tree->Branch("lowMuon_NStrLayers",        &lowMuon_NStrLayers,        "lowMuon_NStrLayers/I");
        sixtracks_tree->Branch("lowMuon_NBPixLayers",        &lowMuon_NBPixLayers,        "lowMuon_NBPixLayers/I");

        sixtracks_tree->Branch("highMuon_isTight",        &highMuon_isTight,        "highMuon_isTight/O");
        sixtracks_tree->Branch("highMuon_isLoose",        &highMuon_isLoose,        "highMuon_isLoose/O");
        sixtracks_tree->Branch("highMuon_isSoft",        &highMuon_isSoft,        "highMuon_isSoft/O");
        sixtracks_tree->Branch("highMuon_isMedium",        &highMuon_isMedium,        "highMuon_isMedium/O");
        sixtracks_tree->Branch("highMuon_isHighPt",        &highMuon_isHighPt,        "highMuon_isHighPt/O");

        sixtracks_tree->Branch("highMuon_isTracker",        &highMuon_isTracker,        "highMuon_isTracker/O");
        sixtracks_tree->Branch("highMuon_isGlobal",        &highMuon_isGlobal,        "highMuon_isGlobal/O");

        sixtracks_tree->Branch("highMuon_NPixelHits",        &highMuon_NPixelHits,        "highMuon_NPixelHits/I");
        sixtracks_tree->Branch("highMuon_NStripHits",        &highMuon_NStripHits,        "highMuon_NStripHits/I");
        sixtracks_tree->Branch("highMuon_NTrackhits",        &highMuon_NTrackhits,        "highMuon_NTrackhits/I");
        sixtracks_tree->Branch("highMuon_NBPixHits",        &highMuon_NBPixHits,        "highMuon_NBPixHits/I");

        sixtracks_tree->Branch("highMuon_NPixLayers",        &highMuon_NPixLayers,        "highMuon_NPixLayers/I");
        sixtracks_tree->Branch("highMuon_NTraLayers",        &highMuon_NTraLayers,        "highMuon_NTraLayers/I");
        sixtracks_tree->Branch("highMuon_NStrLayers",        &highMuon_NStrLayers,        "highMuon_NStrLayers/I");
        sixtracks_tree->Branch("highMuon_NBPixLayers",        &highMuon_NBPixLayers,        "highMuon_NBPixLayers/I");

        sixtracks_tree->Branch("lowMuon_type",     &lowMuon_type,       "lowMuon_type/i");
        sixtracks_tree->Branch("highMuon_type",     &highMuon_type,       "highMuon_type/i");

        //Tracks Flags

        sixtracks_tree->Branch("highTrack_NPixelHits",        &highTrack_NPixelHits,        "highTrack_NPixelHits/I");
        sixtracks_tree->Branch("highTrack_NStripHits",        &highTrack_NStripHits,        "highTrack_NStripHits/I");
        sixtracks_tree->Branch("highTrack_NTrackhits",        &highTrack_NTrackhits,        "highTrack_NTrackhits/I");
        sixtracks_tree->Branch("highTrack_NBPixHits",         &highTrack_NBPixHits,        "highTrack_NBPixHits/I");

        sixtracks_tree->Branch("highTrack_NPixLayers",        &highTrack_NPixLayers,        "highTrack_NPixLayers/I");
        sixtracks_tree->Branch("highTrack_NTraLayers",        &highTrack_NTraLayers,        "highTrack_NTraLayers/I");
        sixtracks_tree->Branch("highTrack_NStrLayers",        &highTrack_NStrLayers,        "highTrack_NStrLayers/I");
        sixtracks_tree->Branch("highTrack_NBPixLayers",       &highTrack_NBPixLayers,        "highTrack_NBPixLayers/I");

        sixtracks_tree->Branch("lowTrack_NPixelHits",        &lowTrack_NPixelHits,        "lowTrack_NPixelHits/I");
        sixtracks_tree->Branch("lowTrack_NStripHits",        &lowTrack_NStripHits,        "lowTrack_NStripHits/I");
        sixtracks_tree->Branch("lowTrack_NTrackhits",        &lowTrack_NTrackhits,        "lowTrack_NTrackhits/I");
        sixtracks_tree->Branch("lowTrack_NBPixHits",         &lowTrack_NBPixHits,        "lowTrack_NBPixHits/I");

        sixtracks_tree->Branch("lowTrack_NPixLayers",        &lowTrack_NPixLayers,        "lowTrack_NPixLayers/I");
        sixtracks_tree->Branch("lowTrack_NTraLayers",        &lowTrack_NTraLayers,        "lowTrack_NTraLayers/I");
        sixtracks_tree->Branch("lowTrack_NStrLayers",        &lowTrack_NStrLayers,        "lowTrack_NStrLayers/I");
        sixtracks_tree->Branch("lowTrack_NBPixLayers",       &lowTrack_NBPixLayers,        "lowTrack_NBPixLayers/I");

        sixtracks_tree->Branch("thirdTrack_NPixelHits",        &thirdTrack_NPixelHits,        "thirdTrack_NPixelHits/I");
        sixtracks_tree->Branch("thirdTrack_NStripHits",        &thirdTrack_NStripHits,        "thirdTrack_NStripHits/I");
        sixtracks_tree->Branch("thirdTrack_NTrackhits",        &thirdTrack_NTrackhits,        "thirdTrack_NTrackhits/I");
        sixtracks_tree->Branch("thirdTrack_NBPixHits",         &thirdTrack_NBPixHits,        "thirdTrack_NBPixHits/I");

        sixtracks_tree->Branch("thirdTrack_NPixLayers",        &thirdTrack_NPixLayers,        "thirdTrack_NPixLayers/I");
        sixtracks_tree->Branch("thirdTrack_NTraLayers",        &thirdTrack_NTraLayers,        "thirdTrack_NTraLayers/I");
        sixtracks_tree->Branch("thirdTrack_NStrLayers",        &thirdTrack_NStrLayers,        "thirdTrack_NStrLayers/I");
        sixtracks_tree->Branch("thirdTrack_NBPixLayers",       &thirdTrack_NBPixLayers,        "thirdTrack_NBPixLayers/I");

        sixtracks_tree->Branch("fourthTrack_NPixelHits",        &fourthTrack_NPixelHits,        "fourthTrack_NPixelHits/I");
        sixtracks_tree->Branch("fourthTrack_NStripHits",        &fourthTrack_NStripHits,        "fourthTrack_NStripHits/I");
        sixtracks_tree->Branch("fourthTrack_NTrackhits",        &fourthTrack_NTrackhits,        "fourthTrack_NTrackhits/I");
        sixtracks_tree->Branch("fourthTrack_NBPixHits",         &fourthTrack_NBPixHits,        "fourthTrack_NBPixHits/I");

        sixtracks_tree->Branch("fourthTrack_NPixLayers",        &fourthTrack_NPixLayers,        "fourthTrack_NPixLayers/I");
        sixtracks_tree->Branch("fourthTrack_NTraLayers",        &fourthTrack_NTraLayers,        "fourthTrack_NTraLayers/I");
        sixtracks_tree->Branch("fourthTrack_NStrLayers",        &fourthTrack_NStrLayers,        "fourthTrack_NStrLayers/I");
        sixtracks_tree->Branch("fourthTrack_NBPixLayers",       &fourthTrack_NBPixLayers,        "fourthTrack_NBPixLayers/I");


        //Six Tracks System
        sixtracks_tree->Branch("six_m_kkpp",             &six_m,               "six_m_kkpp/D");
        sixtracks_tree->Branch("six_m_ref_kkpp",         &six_m_ref,           "six_m_ref_kkpp/D");
        sixtracks_tree->Branch("six_mass_ppkk",      &six_mass_ppkk,        "six_mass_ppkk/D");
        sixtracks_tree->Branch("six_mass_pkpk",      &six_mass_pkpk,        "six_mass_pkpk/D");
        sixtracks_tree->Branch("six_mass_pppp",      &six_mass_pppp,        "six_mass_pppp/D");
        sixtracks_tree->Branch("six_mass_kpkp",      &six_mass_kpkp,        "six_mass_kpkp/D");
        sixtracks_tree->Branch("six_mass_kppk",      &six_mass_kppk,        "six_mass_kppk/D");
        sixtracks_tree->Branch("six_mass_kkkk",      &six_mass_kkkk,        "six_mass_kkkk/D");

        sixtracks_tree->Branch("six_pt",      &six_pt,       "six_pt/D");
        sixtracks_tree->Branch("six_eta",     &six_eta,      "six_eta/D");
        sixtracks_tree->Branch("six_phi",     &six_phi,      "six_phi/D");
        sixtracks_tree->Branch("six_p",       &six_p,        "six_p/D");

        sixtracks_tree->Branch("six_cosAlpha",     &six_cosAlpha,        "six_cosAlpha/D");
        sixtracks_tree->Branch("six_ctauPV",       &six_ctauPV,        "six_ctauPV/D");
        sixtracks_tree->Branch("six_ctauErrPV",    &six_ctauErrPV,        "six_ctauErrPV/D");

        sixtracks_tree->Branch("six_cosAlphaCA",    &six_cosAlphaCA,        "six_cosAlphaCA/D");
        sixtracks_tree->Branch("six_ctauPVCA",      &six_ctauPVCA,        "six_ctauPVCA/D");
        sixtracks_tree->Branch("six_ctauErrPVCA",   &six_ctauErrPVCA,        "six_ctauErrPVCA/D");

        sixtracks_tree->Branch("six_cosAlphaDZ",    &six_cosAlphaDZ,        "six_cosAlphaDZ/D");
        sixtracks_tree->Branch("six_ctauPVDZ",      &six_ctauPVDZ,        "six_ctauPVDZ/D");
        sixtracks_tree->Branch("six_ctauErrPVDZ",   &six_ctauErrPVDZ,        "six_ctauErrPVDZ/D");

        sixtracks_tree->Branch("six_cosAlphaBS",    &six_cosAlphaBS,        "six_cosAlphaBS/D");
        sixtracks_tree->Branch("six_ctauPVBS",      &six_ctauPVBS,        "six_ctauPVBS/D");
        sixtracks_tree->Branch("six_ctauErrPVBS",   &six_ctauErrPVBS,        "six_ctauErrPVBS/D");

        sixtracks_tree->Branch("six_vProb",       &six_vProb,      "six_vProb/D");
        sixtracks_tree->Branch("six_nDof",        &six_nDof,       "six_nDof/D");
        sixtracks_tree->Branch("six_vChi2",       &six_vChi2,      "six_vChi2/D");

        sixtracks_tree->Branch("six_vx",     &six_vx,      "six_vx/D");
        sixtracks_tree->Branch("six_vy",     &six_vy,      "six_vy/D");
        sixtracks_tree->Branch("six_vz",     &six_vz,      "six_vz/D");

        sixtracks_tree->Branch("six_charge",     &six_charge,       "six_charge/I");

      }

      if (IsMC_ || OnlyGen_) {

        sixtracks_tree->Branch("gen_six_p4",    "TLorentzVector",  &gen_six_p4);
        sixtracks_tree->Branch("gen_five_p4",   "TLorentzVector",  &gen_five_p4);
        sixtracks_tree->Branch("gen_mumukk_p4", "TLorentzVector",  &gen_mumukk_p4);
        sixtracks_tree->Branch("gen_jpsi_p4",   "TLorentzVector",  &gen_jpsi_p4);
        sixtracks_tree->Branch("gen_phi_p4",    "TLorentzVector",  &gen_phi_p4);

        sixtracks_tree->Branch("gen_highKaon_p4", "TLorentzVector",  &gen_highKaon_p4);
        sixtracks_tree->Branch("gen_lowMuon_p4",  "TLorentzVector",  &gen_lowMuon_p4);
        sixtracks_tree->Branch("gen_highMuon_p4",  "TLorentzVector",  &gen_highMuon_p4);
        sixtracks_tree->Branch("gen_lowKaon_p4",  "TLorentzVector",  &gen_lowKaon_p4);
        sixtracks_tree->Branch("gen_thirdtrack_p4",  "TLorentzVector",  &gen_thirdtrack_p4);

        sixtracks_tree->Branch("gen_six_pdg",&gen_six_pdg,"gen_six_pdg/D");
        sixtracks_tree->Branch("gen_five_pdg",&gen_five_pdg,"gen_five_pdg/D");
        sixtracks_tree->Branch("gen_mumukk_pdg",&gen_mumukk_pdg,"gen_mumukk_pdg/D");
        sixtracks_tree->Branch("gen_phi_pdg",&gen_phi_pdg,"gen_phi_pdg/D");
        sixtracks_tree->Branch("gen_jpsi_pdg",&gen_jpsi_pdg,"gen_jpsi_pdg/D");

        sixtracks_tree->Branch("gen_lowMuon_pdg",&gen_lowMuon_pdg,"gen_lowMuon_pdg/D");
        sixtracks_tree->Branch("gen_highMuon_pdg",&gen_highMuon_pdg,"gen_highMuon_pdg/D");
        sixtracks_tree->Branch("gen_highKaon_pdg",&gen_highKaon_pdg,"gen_highKaon_pdg/D");
        sixtracks_tree->Branch("gen_lowKaon_pdg",&gen_lowKaon_pdg,"gen_lowKaon_pdg/D");
        sixtracks_tree->Branch("gen_thirdtrack_pdg",&gen_thirdtrack_pdg,"gen_thirdtrack_pdg/D");

        sixtracks_tree->Branch("gen_lowMuon_mompdg",&gen_lowMuon_mompdg,"gen_lowMuon_mompdg/D");
        sixtracks_tree->Branch("gen_highMuon_mompdg",&gen_highMuon_mompdg,"gen_highMuon_mompdg/D");
        sixtracks_tree->Branch("gen_highKaon_mompdg",&gen_highKaon_mompdg,"gen_highKaon_mompdg/D");
        sixtracks_tree->Branch("gen_lowKaon_mompdg",&gen_lowKaon_mompdg,"gen_lowKaon_mompdg/D");
        sixtracks_tree->Branch("gen_thirdtrack_mompdg",&gen_thirdtrack_mompdg,"gen_thirdtrack_mompdg/D");

        sixtracks_tree->Branch("gen_lowMuon_status",&gen_lowMuon_status,"gen_lowMuon_status/D");
        sixtracks_tree->Branch("gen_highMuon_status",&gen_highMuon_status,"gen_highMuon_status/D");
        sixtracks_tree->Branch("gen_highKaon_status",&gen_highKaon_status,"gen_highKaon_status/D");
        sixtracks_tree->Branch("gen_lowKaon_status",&gen_lowKaon_status,"gen_lowKaon_status/D");
        sixtracks_tree->Branch("gen_thirdtrack_status",&gen_thirdtrack_status,"gen_thirdtrack_status/D");

        sixtracks_tree->Branch("gen_lowMuon_p",&gen_lowMuon_p,"gen_lowMuon_p/D");
        sixtracks_tree->Branch("gen_highMuon_p",&gen_highMuon_p,"gen_highMuon_p/D");
        sixtracks_tree->Branch("gen_highKaon_p",&gen_highKaon_p,"gen_highKaon_p/D");
        sixtracks_tree->Branch("gen_lowKaon_p",&gen_lowKaon_p,"gen_lowKaon_p/D");
        sixtracks_tree->Branch("gen_thirdtrack_p",&gen_thirdtrack_p,"gen_thirdtrack_p/D");

        sixtracks_tree->Branch("gen_lowMuon_pt",&gen_lowMuon_pt,"gen_lowMuon_pt/D");
        sixtracks_tree->Branch("gen_highMuon_pt",&gen_highMuon_pt,"gen_highMuon_pt/D");
        sixtracks_tree->Branch("gen_highKaon_pt",&gen_highKaon_pt,"gen_highKaon_pt/D");
        sixtracks_tree->Branch("gen_lowKaon_pt",&gen_lowKaon_pt,"gen_lowKaon_pt/D");
        sixtracks_tree->Branch("gen_thirdtrack_pt",&gen_thirdtrack_pt,"gen_thirdtrack_pt/D");

        sixtracks_tree->Branch("gen_lowMuon_eta",&gen_lowMuon_eta,"gen_lowMuon_eta/D");
        sixtracks_tree->Branch("gen_highMuon_eta",&gen_highMuon_eta,"gen_highMuon_eta/D");
        sixtracks_tree->Branch("gen_highKaon_eta",&gen_highKaon_eta,"gen_highKaon_eta/D");
        sixtracks_tree->Branch("gen_lowKaon_eta",&gen_lowKaon_eta,"gen_lowKaon_eta/D");
        sixtracks_tree->Branch("gen_thirdtrack_eta",&gen_thirdtrack_eta,"gen_thirdtrack_eta/D");

        sixtracks_tree->Branch("gen_lowMuon_phi",&gen_lowMuon_phi,"gen_lowMuon_phi/D");
        sixtracks_tree->Branch("gen_highMuon_phi",&gen_highMuon_phi,"gen_highMuon_phi/D");
        sixtracks_tree->Branch("gen_highKaon_phi",&gen_highKaon_phi,"gen_highKaon_phi/D");
        sixtracks_tree->Branch("gen_lowKaon_phi",&gen_lowKaon_phi,"gen_lowKaon_phi/D");
        sixtracks_tree->Branch("gen_thirdtrack_phi",&gen_thirdtrack_phi,"gen_thirdtrack_phi/D");

        sixtracks_tree->Branch("gen_five_prompt",&gen_five_prompt,"gen_five_prompt/D");
        sixtracks_tree->Branch("gen_five_pt",&gen_five_pt,"gen_five_pt/D");
        sixtracks_tree->Branch("gen_five_p",&gen_five_p,"gen_five_p/D");
        sixtracks_tree->Branch("gen_five_eta",&gen_five_eta,"gen_five_eta/D");
        sixtracks_tree->Branch("gen_five_phi",&gen_five_phi,"gen_five_phi/D");

        sixtracks_tree->Branch("gen_mumukk_prompt",&gen_mumukk_prompt,"gen_mumukk_prompt/D");
        sixtracks_tree->Branch("gen_phi_prompt",&gen_phi_prompt,"gen_phi_prompt/D");
        sixtracks_tree->Branch("gen_jpsi_prompt",&gen_jpsi_prompt,"gen_jpsi_prompt/D");

        sixtracks_tree->Branch("gen_mumukk_pt",&gen_mumukk_pt,"gen_mumukk_pt/D");
        sixtracks_tree->Branch("gen_phi_pt",&gen_phi_pt,"phigen_phi_pt_pt/D");
        sixtracks_tree->Branch("gen_jpsi_pt",&gen_jpsi_pt,"gen_jpsi_pt/D");

        sixtracks_tree->Branch("gen_mumukk_p",&gen_mumukk_p,"gen_mumukk_p/D");
        sixtracks_tree->Branch("gen_phi_p",&gen_phi_p,"phigen_phi_p_p/D");
        sixtracks_tree->Branch("gen_jpsi_p",&gen_jpsi_p,"gen_jpsi_p/D");

        sixtracks_tree->Branch("gen_mumukk_eta",&gen_mumukk_eta,"gen_mumukk_eta/D");
        sixtracks_tree->Branch("gen_phi_eta",&gen_phi_eta,"gen_phi_eta/D");
        sixtracks_tree->Branch("gen_jpsi_eta",&gen_jpsi_eta,"gen_jpsi_eta/D");

        sixtracks_tree->Branch("gen_mumukk_phi",&gen_mumukk_phi,"gen_mumukk_phi/D");
        sixtracks_tree->Branch("gen_phi_phi",&gen_phi_phi,"gen_phi_phi/D");
        sixtracks_tree->Branch("gen_jpsi_phi",&gen_jpsi_phi,"gen_jpsi_phi/D");

        sixtracks_tree->Branch("gen_fourthtrackp4",  "TLorentzVector",  &gen_fourthtrackp4);
        sixtracks_tree->Branch("gen_fourthtrackpdg",&gen_fourthtrackpdg,"gen_fourthtrackpdg/D");
        sixtracks_tree->Branch("gen_fourthtrackmompdg",&gen_fourthtrackmompdg,"gen_fourthtrackmompdg/D");
        sixtracks_tree->Branch("gen_fourthtrackstatus",&gen_fourthtrackstatus,"gen_fourthtrackstatus/D");
        sixtracks_tree->Branch("gen_fourthtrackp",&gen_fourthtrackp,"gen_fourthtrackp/D");
        sixtracks_tree->Branch("gen_fourthtrackpt",&gen_fourthtrackpt,"gen_fourthtrackpt/D");
        sixtracks_tree->Branch("gen_fourthtracketa",&gen_fourthtracketa,"gen_fourthtracketa/D");
        sixtracks_tree->Branch("gen_fourthtrackphi",&gen_fourthtrackphi,"gen_fourthtrackphi/D");

      }

        genCands_ = consumes< std::vector <reco::GenParticle> >((edm::InputTag)"prunedGenParticles");
        packCands_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");

}

SixTracksRootupler::~SixTracksRootupler() {}

//
// member functions
//

bool SixTracksRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   for (size_t i=0; i< particle->numberOfMothers(); i++) {
      if (isAncestor(ancestor, particle->mother(i))) return true;
   }
   return false;
}


// ------------ method called for each event  ------------
void SixTracksRootupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
//  using namespace edm;
  using namespace std;

  edm::Handle<std::vector<pat::CompositeCandidate>> sixtracks_cand_handle;
  iEvent.getByToken(SixTracksCollection_, sixtracks_cand_handle);

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

gen_six_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_mumukk_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_five_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_phi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_lowMuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_highMuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_highKaon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_lowKaon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_thirdtrack_p4.SetPtEtaPhiM(0.,0.,0.,0.);

gen_mumukk_pdg = 0;
gen_five_pdg = 0;

//std::cout << "Debug  1" << std::endl;
if(!OnlyGen_)
  if (!sixtracks_cand_handle.isValid()) std::cout<< "No five tracks information " << run << "," << event <<std::endl;

  if (sixtracks_cand_handle.isValid()) {


    //std::cout << "Debug  2" << std::endl;
    noSixCandidates = (Int_t)(sixtracks_cand_handle->size());

    for (unsigned int i=0; i< sixtracks_cand_handle->size(); i++)
    {

      pat::CompositeCandidate six_cand;

      // const pat::CompositeCandidate *six_ref;

      const pat::PackedCandidate *trackOne_cand, *trackTwo_cand, *trackThree_cand,*trackFour_cand;

      const pat::CompositeCandidate *mumukk_cand;

      const pat::CompositeCandidate *kk_cand;

      const pat::CompositeCandidate *mumu_cand;

      const pat::CompositeCandidate *triTrack_cand;

      const pat::CompositeCandidate *five_cand;

      six_cand  = sixtracks_cand_handle->at(i);

      // six_ref       = dynamic_cast<const pat::CompositeCandidate*>(six_cand.daughter("first_six_ref")); // MUMU KK PP
      five_cand     = dynamic_cast<const pat::CompositeCandidate*>(six_cand.daughter("fiveCand")); // KK
      mumukk_cand   = dynamic_cast<const pat::CompositeCandidate*>(five_cand->daughter("dimuonditrack")); // P N
      mumu_cand     = dynamic_cast<const pat::CompositeCandidate*>(mumukk_cand->daughter("dimuon"));
      kk_cand       = dynamic_cast<const pat::CompositeCandidate*>(mumukk_cand->daughter("ditrack"));

      const pat::Muon *lowMuon, *highMuon;

      reco::Candidate::LorentzVector vhighMuon = mumu_cand->daughter("highMuon")->p4();
      reco::Candidate::LorentzVector vlowMuon = mumu_cand->daughter("lowMuon")->p4();

      // if (mumu_cand->daughter("highMuon")->charge() < 0) {
      if (vhighMuon.pt() < vlowMuon.pt()) {
        vhighMuon = mumu_cand->daughter("lowMuon")->p4();
        vlowMuon = mumu_cand->daughter("highMuon")->p4();
        highMuon = dynamic_cast<const pat::Muon*>(mumu_cand->daughter("highMuon"));
        lowMuon = dynamic_cast<const pat::Muon*>(mumu_cand->daughter("lowMuon"));
      } else
      {
        lowMuon = dynamic_cast<const pat::Muon*>(mumu_cand->daughter("highMuon"));
        highMuon = dynamic_cast<const pat::Muon*>(mumu_cand->daughter("lowMuon"));
      }

      swapped  = six_cand.userFloat("swapped");
      sameSign = six_cand.userFloat("isSameSign");
      sameSign_mmtt = mumukk_cand->userFloat("isSameSign");

      const reco::Vertex bestPV = *(six_cand.userData<reco::Vertex>("bestPV"));
      const reco::Vertex cosPV  = *(six_cand.userData<reco::Vertex>("cosPV"));
      const reco::Vertex zPV    = *(six_cand.userData<reco::Vertex>("zPV"));
      const reco::Vertex bs     = *(six_cand.userData<reco::Vertex>("bS"));

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

      highTrackMatch   = (Double_t)mumukk_cand->userInt("highKaonMatch");
      lowTrackMatch    = (Double_t)mumukk_cand->userInt("lowKaonMatch");
      lowMuonMatch     = (Double_t)mumu_cand->userInt("highMuonTMatch");
      highMuonMatch    = (Double_t)mumu_cand->userInt("lowMuonTMatch");
      thirdTrackMatch  = (Double_t)five_cand->userInt("thirdKaonMatch");
      fourthTrackMatch = (Double_t)six_cand.userInt("fourthKaonMatch");

      trackOne_cand   = dynamic_cast<const pat::PackedCandidate*>(six_cand.daughter("trackOne"));
      trackTwo_cand   = dynamic_cast<const pat::PackedCandidate*>(six_cand.daughter("trackTwo"));
      trackThree_cand = dynamic_cast<const pat::PackedCandidate*>(six_cand.daughter("trackThree"));
      trackFour_cand  = dynamic_cast<const pat::PackedCandidate*>(six_cand.daughter("trackFour"));

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

      fourthKaon_p4.SetPtEtaPhiM(trackFour_cand->pt(), trackFour_cand->eta(), trackFour_cand->phi(), trackFour_cand->mass());
      fourthPion_p4.SetPtEtaPhiM(trackFour_cand->pt(), trackFour_cand->eta(), trackFour_cand->phi(), pionmass);
      fourthProton_p4.SetPtEtaPhiM(trackFour_cand->pt(), trackFour_cand->eta(), trackFour_cand->phi(), protonmass);

      //
      //
      // fiveTwo_cand   = dynamic_cast<const pat::CompositeCandidate*>(six_ref->daughter("fiveTrackTwo")); // P N 4
      // fiveThree_cand = dynamic_cast<const pat::CompositeCandidate*>(six_ref->daughter("fiveTrackThree")); // P 3 4
      // fiveFour_cand  = dynamic_cast<const pat::CompositeCandidate*>(six_ref->daughter("fiveTrackFour")); // N 3 4
      //
      //
      // dimuonDiTrkTwo_cand   = dynamic_cast<const pat::CompositeCandidate*>(five_thr.daughter("dimuonDiTrackOne")); // P 3
      // dimuonDiTrkThree_cand = dynamic_cast<const pat::CompositeCandidate*>(five_thr.daughter("dimuonDiTrackTwo")); // P 4
      // dimuonDiTrkFour_cand  = dynamic_cast<const pat::CompositeCandidate*>(five_fou.daughter("dimuonDiTrackOne"));  // N 3
      // dimuonDiTrkFive_cand  = dynamic_cast<const pat::CompositeCandidate*>(five_fou.daughter("dimuonDiTrackTwo")); // N 4
      // dimuonDiTrkSix_cand   = dynamic_cast<const pat::CompositeCandidate*>(five_fou.daughter("dimuonDiTrackThree")); // 3 4
      //
      //
      // diTrackTwo_cand     = dynamic_cast<const pat::CompositeCandidate*>(four_two.daughter("ditrack"));
      // kk_thr.   = dynamic_cast<const pat::CompositeCandidate*>(four_thr.daughter("ditrack"));
      // kk_fou.    = dynamic_cast<const pat::CompositeCandidate*>(four_fou.daughter("ditrack"));
      // kk_fiv.    = dynamic_cast<const pat::CompositeCandidate*>(four_fiv.daughter("ditrack"));
      // kk_six.     = dynamic_cast<const pat::CompositeCandidate*>(four_six.daughter("ditrack"));

      five_id = six_cand.userInt("five_index");
      p1_id = six_cand.userInt("p1Id");
      m1_id = six_cand.userInt("m1Id");
      p2_id = six_cand.userInt("p2Id");
      m2_id = six_cand.userInt("m2Id");
      dimuon_id = mumukk_cand->userInt("dimuon_id");



      six_m    = six_cand.mass();
      six_m_ref       = six_cand.userFloat("mass_ref");

      six_mass_ppkk    = (highMuon_p4 + lowMuon_p4 + highPion_p4 + lowPion_p4 + thirdKaon_p4 + fourthKaon_p4).M();
      six_mass_pkpk    = (highMuon_p4 + lowMuon_p4 + highPion_p4 + lowKaon_p4 + thirdPion_p4 + fourthPion_p4).M();
      six_mass_kpkp    = (highMuon_p4 + lowMuon_p4 + highKaon_p4 + lowPion_p4 + thirdKaon_p4 + fourthPion_p4).M();
      six_mass_kppk    = (highMuon_p4 + lowMuon_p4 + highKaon_p4 + lowPion_p4 + thirdPion_p4 + fourthKaon_p4).M();

      six_mass_kkkk    = (highMuon_p4 + lowMuon_p4 + highKaon_p4 + lowKaon_p4 + thirdKaon_p4 + fourthKaon_p4).M();
      six_mass_pppp    = (highMuon_p4 + lowMuon_p4 + highPion_p4 + lowPion_p4 + thirdPion_p4 + fourthPion_p4).M();

      six_pt    = six_cand.pt();
      six_eta   = six_cand.eta();
      six_phi   = six_cand.phi();
      six_p     = six_cand.p();

      six_cosAlpha    = six_cand.userFloat("cosAlpha");
      six_ctauPV      = six_cand.userFloat("ctauPV");
      six_ctauErrPV   = six_cand.userFloat("ctauErrPV");
      six_cosAlphaCA  = six_cand.userFloat("cosAlphaCA");
      six_ctauPVCA    = six_cand.userFloat("ctauPVCA");
      six_ctauErrPVCA = six_cand.userFloat("ctauErrPVCA");
      six_cosAlphaDZ  = six_cand.userFloat("cosAlphaDZ");
      six_ctauPVDZ    = six_cand.userFloat("ctauPVDZ");
      six_ctauErrPVDZ = six_cand.userFloat("ctauErrPVDZ");
      six_cosAlphaBS  = six_cand.userFloat("cosAlphaBS");
      six_ctauPVBS    = six_cand.userFloat("ctauPVBS");
      six_ctauErrPVBS = six_cand.userFloat("ctauErrPVBS");

      six_vProb   = six_cand.userFloat("vProb");
      six_nDof    = six_cand.userFloat("nDof");
      six_vChi2   = six_cand.userFloat("vChi2");

      six_vx = six_cand.userFloat("vtxX");
      six_vy = six_cand.userFloat("vtxY");
      six_vz = six_cand.userFloat("vtxZ");
      six_charge = six_cand.charge();


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

      mumukk_vProb     = mumukk_cand->userFloat("vProb");
      mumukk_vChi2     = mumukk_cand->userFloat("vChi2");
      mumukk_nDof      = mumukk_cand->userFloat("nDof");
      mumukk_charge    = mumukk_cand->charge();

      mumukk_cosAlphaBS = mumukk_cand->userFloat("cosAlphaBS");
      mumukk_ctauPVBS = mumukk_cand->userFloat("ctauPVBS");
      mumukk_ctauErrPVBS = mumukk_cand->userFloat("ctauErrPVBS");

      // highTrackFromPVBS = mumukk_cand->userFloat("highTrackFromPVBS");
      // lowTrackFromPVBS = mumukk_cand->userFloat("lowTrackFromPVBS");

      auto five_two = highMuon_p4 + lowMuon_p4 + highKaon_p4 + lowKaon_p4 + fourthKaon_p4; //3 4 6 + - - (~5)
      auto five_thr = highMuon_p4 + lowMuon_p4 + highKaon_p4 + thirdKaon_p4 + fourthKaon_p4; // 3 5 6 + + -  (~4)
      auto five_fou = highMuon_p4 + lowMuon_p4 + lowKaon_p4  + thirdKaon_p4 + fourthKaon_p4; // 4 5 6 - + - (~3)


      auto four_two = highMuon_p4 + lowMuon_p4 + highKaon_p4 + fourthKaon_p4; //3 6 + -
      auto four_thr = highMuon_p4 + lowMuon_p4 + lowKaon_p4 + thirdKaon_p4; //4 5 - +
      auto four_fou = highMuon_p4 + lowMuon_p4 + fourthKaon_p4 + thirdKaon_p4; //5 6 + -

      auto four_fiv = highMuon_p4 + lowMuon_p4 + highKaon_p4 + thirdKaon_p4; //3 5 + +
      auto four_six = highMuon_p4 + lowMuon_p4 + lowKaon_p4 + fourthKaon_p4; //4 6 - -


      fiveOne_pt       = five_cand->pt();
      fiveOne_eta      = five_cand->eta();
      fiveOne_phi      = five_cand->phi();
      fiveOne_p        = five_cand->p();

      fiveTwo_pt       = five_two.Pt();
      fiveTwo_eta      = five_two.Eta();
      fiveTwo_phi      = five_two.Phi();
      fiveTwo_p        = five_two.P();

      fiveThree_pt       = five_thr.Pt();
      fiveThree_eta      = five_thr.Eta();
      fiveThree_phi      = five_thr.Phi();
      fiveThree_p        = five_thr.P();

      fiveFour_pt       = five_fou.Pt();
      fiveFour_eta      = five_fou.Eta();
      fiveFour_phi      = five_fou.Phi();
      fiveFour_p        = five_fou.P();

      five_vProb    = five_cand->userFloat("vProb");
      five_nDof     = five_cand->userFloat("nDof");
      five_vChi2    = five_cand->userFloat("vChi2");
      five_charge   = five_cand->charge();

      five_m           = five_cand->mass();
      five_m_ref       = five_cand->userFloat("mass_ref_0");
      five_mass_ppk    = five_cand->userFloat("mass_ref_1");
      five_mass_kpp    = five_cand->userFloat("mass_ref_2");
      five_mass_pkp    = five_cand->userFloat("mass_ref_3");
      five_mass_ppp    = five_cand->userFloat("mass_ref_4");

      dca_m1t3 = five_cand->userFloat("dca_m1t3");
      dca_m2t3 = five_cand->userFloat("dca_m2t3");
      dca_t1t3 = five_cand->userFloat("dca_t1t3");
      dca_t2t3 = five_cand->userFloat("dca_t2t3");
      dca_m1t4 = six_cand.userFloat("dca_m1t4");
      dca_m2t4 = six_cand.userFloat("dca_m2t4");
      dca_t1t4 = six_cand.userFloat("dca_t1t4");
      dca_t2t4 = six_cand.userFloat("dca_t2t4");
      dca_t3t4 = six_cand.userFloat("dca_t3t4");

      thirdTrackFromPV        = five_cand->userFloat("thirdTrackFromPV");
      thirdTrackFromPVDZ      = five_cand->userFloat("thirdTrackFromPVDZ");
      // thirdTrackFromPVBS      = five_cand->userFloat("thirdTrackFromPVBS");
      thirdTrackFromPVCA      = five_cand->userFloat("thirdTrackFromPVCA");


      mumukk_cosAlpha = mumukk_cand->userFloat("cosAlpha");
      mumukk_ctauPV = mumukk_cand->userFloat("ctauPV");
      mumukk_ctauErrPV = mumukk_cand->userFloat("ctauErrPV");

      highTrackFromPV = mumukk_cand->userFloat("highTrackFromPV");
      lowTrackFromPV = mumukk_cand->userFloat("lowTrackFromPV");

      mumukk_cosAlphaDZ = mumukk_cand->userFloat("cosAlphaDZ");
      mumukk_ctauPVDZ = mumukk_cand->userFloat("ctauPVDZ");
      mumukk_ctauErrPVDZ = mumukk_cand->userFloat("ctauErrPVDZ");

      highTrackFromPVDZ = mumukk_cand->userFloat("highTrackFromPVDZ");
      lowTrackFromPVDZ = mumukk_cand->userFloat("lowTrackFromPVDZ");

      mumukk_cosAlphaCA = mumukk_cand->userFloat("cosAlphaCA");
      mumukk_ctauPVCA = mumukk_cand->userFloat("ctauPVCA");
      mumukk_ctauErrPVCA = mumukk_cand->userFloat("ctauErrPVCA");

      highTrackFromPVCA = mumukk_cand->userFloat("highTrackFromPVCA");
      lowTrackFromPVCA = mumukk_cand->userFloat("lowTrackFromPVCA");

      dca_m1m2 = mumukk_cand->userFloat("dca_m1m2");
      dca_m1t1 = mumukk_cand->userFloat("dca_m1t1");
      dca_m1t2 = mumukk_cand->userFloat("dca_m1t2");
      dca_m2t1 = mumukk_cand->userFloat("dca_m2t1");
      dca_m2t2 = mumukk_cand->userFloat("dca_m2t2");
      dca_t1t2 = mumukk_cand->userFloat("dca_t1t2");

      dimuon_vProb        = mumu_cand->userFloat("vProb");
      dimuon_vChi2        = mumu_cand->userFloat("vNChi2");
      dimuon_DCA          = mumu_cand->userFloat("DCA");
      dimuon_ctauPV       = mumu_cand->userFloat("ppdlPV");
      dimuon_ctauErrPV    = mumu_cand->userFloat("ppdlErrPV");
      dimuon_cosAlpha     = mumu_cand->userFloat("cosAlpha");

      dimuonDiTrkOne_mmpp   = (highMuon_p4 + lowMuon_p4 + highPion_p4 + lowPion_p4).M(); // 34
      dimuonDiTrkTwo_mmpp   = (highMuon_p4 + lowMuon_p4 + highPion_p4 + fourthPion_p4).M(); // 36
      dimuonDiTrkThree_mmpp = (highMuon_p4 + lowMuon_p4 + lowPion_p4  + thirdPion_p4).M(); //45
      dimuonDiTrkFour_mmpp  = (highMuon_p4 + lowMuon_p4 + fourthPion_p4 + thirdPion_p4).M(); //64

      dimuonDiTrkOne_mmkk   = (highMuon_p4 + lowMuon_p4 + highKaon_p4 + lowKaon_p4).M();
      dimuonDiTrkTwo_mmkk   = (highMuon_p4 + lowMuon_p4 + highKaon_p4 + fourthKaon_p4).M();
      dimuonDiTrkThree_mmkk = (highMuon_p4 + lowMuon_p4 + lowKaon_p4  + thirdKaon_p4).M();
      dimuonDiTrkFour_mmkk  = (highMuon_p4 + lowMuon_p4 + fourthKaon_p4 + thirdKaon_p4).M();

      dimuonDiTrkOne_mmpk   = (highMuon_p4 + lowMuon_p4 + highPion_p4 + lowKaon_p4).M();
      dimuonDiTrkTwo_mmpk   = (highMuon_p4 + lowMuon_p4 + highPion_p4 + fourthPion_p4).M();
      dimuonDiTrkThree_mmpk = (highMuon_p4 + lowMuon_p4 + lowPion_p4  + thirdKaon_p4).M();
      dimuonDiTrkFour_mmpk  = (highMuon_p4 + lowMuon_p4 + fourthPion_p4 + thirdKaon_p4).M();

      dimuonDiTrkOne_mmkp   = (highMuon_p4 + lowMuon_p4 + highKaon_p4 + lowPion_p4).M();
      dimuonDiTrkTwo_mmkp   = (highMuon_p4 + lowMuon_p4 + highKaon_p4 + fourthKaon_p4).M();
      dimuonDiTrkThree_mmkp = (highMuon_p4 + lowMuon_p4 + lowKaon_p4  + thirdPion_p4).M();
      dimuonDiTrkFour_mmkp  = (highMuon_p4 + lowMuon_p4 + fourthKaon_p4 + thirdPion_p4).M();

      dimuonDiTrkSix_mmpp = (highMuon_p4 + lowMuon_p4 + highPion_p4 + thirdPion_p4).M();
      dimuonDiTrkSix_mmkk = (highMuon_p4 + lowMuon_p4 + highKaon_p4  + thirdKaon_p4).M();
      dimuonDiTrkSix_mmpk = (highMuon_p4 + lowMuon_p4 + highPion_p4 + thirdKaon_p4).M();
      dimuonDiTrkSix_mmkp = (highMuon_p4 + lowMuon_p4 + highKaon_p4 + thirdPion_p4).M();

      dimuonDiTrkFive_mmpp = (highMuon_p4 + lowMuon_p4 + lowPion_p4 + fourthPion_p4).M();
      dimuonDiTrkFive_mmkk = (highMuon_p4 + lowMuon_p4 + lowKaon_p4 + fourthKaon_p4).M();
      dimuonDiTrkFive_mmpk = (highMuon_p4 + lowMuon_p4 + lowPion_p4 + fourthKaon_p4).M();
      dimuonDiTrkFive_mmkp = (highMuon_p4 + lowMuon_p4 + lowKaon_p4 + fourthPion_p4).M();

      // fiveTrackOne_cand    = dynamic_cast<const pat::CompositeCandidate*>(six_ref->daughter("fiveTrackOne"));
      // fiveTrackTwo_cand    = dynamic_cast<const pat::CompositeCandidate*>(six_ref->daughter("fiveTrackTwo"));
      // fiveTrackThree_cand  = dynamic_cast<const pat::CompositeCandidate*>(six_ref->daughter("fiveTrackThree"));
      // fiveTrackFour_cand   = dynamic_cast<const pat::CompositeCandidate*>(six_ref->daughter("fiveTrackFour"));

      // triTrack_cand = dynamic_cast<const pat::CompositeCandidate*>(six_ref->daughter("triTrack"));

      mumukk_p4.SetPtEtaPhiM(mumukk_cand->pt(),mumukk_cand->eta(),mumukk_cand->phi(),mumukk_cand->mass());
      mumukk_m = mumukk_cand->mass();
      mumukk_pt = mumukk_cand->pt();
      mumukk_eta = mumukk_cand->eta();
      mumukk_phi = mumukk_cand->phi();
      mumukk_p = mumukk_cand->p();

      dimuon_p4.SetPtEtaPhiM(mumu_cand->pt(),mumu_cand->eta(),mumu_cand->phi(),mumu_cand->mass());
      dimuon_m = mumu_cand->mass();
      dimuon_pt = mumu_cand->pt();
      dimuon_eta = mumu_cand->eta();
      dimuon_phi = mumu_cand->phi();
      dimuon_p = mumu_cand->p();

      ditrack_p4.SetPtEtaPhiM(kk_cand->pt(),kk_cand->eta(),kk_cand->phi(),kk_cand->mass());
      ditrack_m = kk_cand->mass();

      diTrackOne_pt = kk_cand->pt();
      diTrackOne_eta = kk_cand->eta();
      diTrackOne_phi = kk_cand->phi();
      diTrackOne_p = kk_cand->p();

      auto kk_two = highKaon_p4 + thirdKaon_p4; //3 5 + +
      auto kk_thr = highKaon_p4 + fourthKaon_p4; //3 6 + -
      auto kk_fou = lowKaon_p4 + thirdKaon_p4; //4 5 - +
      auto kk_fiv = lowKaon_p4 + fourthKaon_p4; //4 6 - -
      auto kk_six = fourthKaon_p4 + thirdKaon_p4; //5 6 + -

      diTrackOne_kk = ditrack_m;
      diTrackTwo_kk = kk_two.M();
      diTrackThree_kk = kk_thr.M();
      diTrackFour_kk = kk_fou.M();
      diTrackFive_kk = kk_fiv.M();
      diTrackSix_kk = kk_six.M();

      diTrackTwo_pt = kk_two.Pt();
      diTrackTwo_eta = kk_two.Eta();
      diTrackTwo_phi = kk_two.Phi();
      diTrackTwo_p = kk_two.P();

      diTrackThree_pt = kk_thr.Pt();
      diTrackThree_eta = kk_thr.Eta();
      diTrackThree_phi = kk_thr.Phi();
      diTrackThree_p = kk_thr.P();

      diTrackFour_pt = kk_fou.Pt();
      diTrackFour_eta = kk_fou.Eta();
      diTrackFour_phi = kk_fou.Phi();
      diTrackFour_p = kk_fou.P();

      diTrackFive_pt = kk_fiv.Pt();
      diTrackFive_eta = kk_fiv.Eta();
      diTrackFive_phi = kk_fiv.Phi();
      diTrackFive_p = kk_fiv.P();

      diTrackSix_pt = kk_six.Pt();
      diTrackSix_eta = kk_six.Eta();
      diTrackSix_phi = kk_six.Phi();
      diTrackSix_p = kk_six.P();

      diTrackOne_pp = (highPion_p4 + lowPion_p4).M();
      diTrackTwo_pp = (highPion_p4 + thirdPion_p4).M(); //3 5 + +
      diTrackThree_pp = (highPion_p4 + fourthPion_p4).M(); //3 6 + -
      diTrackFour_pp = (lowPion_p4 + thirdPion_p4).M(); //4 5 - +
      diTrackFive_pp = (lowPion_p4 + fourthPion_p4).M(); //4 6 - -
      diTrackSix_pp = (fourthPion_p4 + thirdPion_p4).M(); //5 6 + -

      diTrackOne_pk = (highPion_p4 + lowKaon_p4).M();
      diTrackTwo_pk = (highPion_p4 + thirdKaon_p4).M(); //3 5 + +
      diTrackThree_pk = (highPion_p4 + fourthKaon_p4).M(); //3 6 + -
      diTrackFour_pk = (lowPion_p4 + thirdKaon_p4).M(); //4 5 - +
      diTrackFive_pk = (lowPion_p4 + fourthKaon_p4).M(); //4 6 - -
      diTrackSix_pk = (fourthPion_p4 + thirdKaon_p4).M(); //5 6 + -

      diTrackOne_kp = (highKaon_p4 + lowPion_p4).M();
      diTrackTwo_kp = (highKaon_p4 + thirdPion_p4).M(); //3 5 + +
      diTrackThree_kp = (highKaon_p4 + fourthPion_p4).M(); //3 6 + -
      diTrackFour_kp = (lowKaon_p4 + thirdPion_p4).M(); //4 5 - +
      diTrackFive_kp = (lowKaon_p4 + fourthPion_p4).M(); //4 6 - -
      diTrackSix_kp = (fourthKaon_p4 + thirdPion_p4).M(); //5 6 + -

      //trackOne_cand->SetPtEtaPhiM(mumu_cand->pt(),mumu_cand->eta(),mumu_cand->phi(),mumu_cand->mass());

      highTrackMuonDR = mumukk_cand->userFloat("highKaonMuonDR");
      highTrackMuonDP = mumukk_cand->userFloat("highKaonMuonDP");
      highTrackMuonDPt = mumukk_cand->userFloat("highKaonMuonDPt");

      lowTrackMuonDR   = mumukk_cand->userFloat("lowKaonMuonDR");
      lowTrackMuonDP   = mumukk_cand->userFloat("lowKaonMuonDP");
      lowTrackMuonDPt  = mumukk_cand->userFloat("lowKaonMuonDPt");

      thirdTrackMuonDR   = five_cand->userFloat("thirdTrackMuonDR");
      thirdTrackMuonDP   = five_cand->userFloat("thirdTrackMuonDP");
      thirdTrackMuonDPt  = five_cand->userFloat("thirdTrackMuonDPt");


      mumukk_vx    = mumukk_cand->userFloat("vtxX");
      mumukk_vy    = mumukk_cand->userFloat("vtxY");
      mumukk_vz    = mumukk_cand->userFloat("vtxZ");

      five_vx    = five_cand->userFloat("vtxX");
      five_vy    = five_cand->userFloat("vtxY");
      five_vz    = five_cand->userFloat("vtxZ");

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

      fourthTrack_pt     = trackFour_cand->pt();
      fourthTrack_eta    = trackFour_cand->eta();
      fourthTrack_phi    = trackFour_cand->phi();
      fourthTrack_charge = trackFour_cand->charge();
      fourthTrack_dz     = trackFour_cand->bestTrack()->dz();
      fourthTrack_dxy    = trackFour_cand->bestTrack()->dxy();

      fourthTrack_NPixelHits = trackFour_cand->bestTrack()->hitPattern().numberOfValidPixelHits();
      fourthTrack_NStripHits = trackFour_cand->bestTrack()->hitPattern().numberOfValidStripHits();
      fourthTrack_NTrackhits = trackFour_cand->bestTrack()->hitPattern().numberOfValidTrackerHits();
      fourthTrack_NBPixHits  = trackFour_cand->bestTrack()->hitPattern().numberOfValidStripHits();
      fourthTrack_NPixLayers = trackFour_cand->bestTrack()->hitPattern().pixelLayersWithMeasurement();
      fourthTrack_NTraLayers = trackFour_cand->bestTrack()->hitPattern().trackerLayersWithMeasurement();
      fourthTrack_NStrLayers = trackFour_cand->bestTrack()->hitPattern().stripLayersWithMeasurement();
      fourthTrack_NBPixLayers = trackFour_cand->bestTrack()->hitPattern().pixelBarrelLayersWithMeasurement();


      dimuonDiTrkOne_pt     = mumukk_cand->pt();
      dimuonDiTrkOne_eta    = mumukk_cand->eta();
      dimuonDiTrkOne_phi    = mumukk_cand->phi();
      dimuonDiTrkOne_charge = mumukk_cand->charge();

      dimuonDiTrkTwo_pt     = four_two.Pt();
      dimuonDiTrkTwo_eta    = four_two.Eta();
      dimuonDiTrkTwo_phi    = four_two.Phi();
      dimuonDiTrkTwo_charge = highTrack_charge + thirdTrack_charge;

      dimuonDiTrkThree_pt     = four_thr.Pt();
      dimuonDiTrkThree_eta    = four_thr.Eta();
      dimuonDiTrkThree_phi    = four_thr.Phi();
      dimuonDiTrkThree_charge = highTrack_charge + fourthTrack_charge;

      dimuonDiTrkFour_pt     = four_fou.Pt();
      dimuonDiTrkFour_eta    = four_fou.Eta();
      dimuonDiTrkFour_phi    = four_fou.Phi();
      dimuonDiTrkFour_charge = lowTrack_charge + thirdTrack_charge;

      dimuonDiTrkFive_pt     = four_fiv.Pt();
      dimuonDiTrkFive_eta    = four_fiv.Eta();
      dimuonDiTrkFive_phi    = four_fiv.Phi();
      dimuonDiTrkFive_charge = lowTrack_charge + fourthTrack_charge;

      dimuonDiTrkSix_pt     = four_six.Pt();
      dimuonDiTrkSix_eta    = four_six.Eta();
      dimuonDiTrkSix_phi    = four_six.Phi();
      dimuonDiTrkSix_charge = thirdTrack_charge + fourthTrack_charge;

      // triTrack_m      = triTrack_cand->mass();
      auto triTrack_one = (highKaon_p4 + lowKaon_p4 + thirdKaon_p4);
      auto triTrack_two = (highKaon_p4 + lowKaon_p4 + fourthKaon_p4);
      auto triTrack_thr = (highKaon_p4 + fourthKaon_p4 + thirdKaon_p4);
      auto triTrack_fou = (lowKaon_p4 + fourthKaon_p4 + thirdKaon_p4);

      triTrackOne_pt     = triTrack_one.Pt();
      triTrackOne_eta    = triTrack_one.Eta();
      triTrackOne_phi    = triTrack_one.Phi();
      triTrackOne_charge = lowTrack_charge + thirdTrack_charge + highTrack_charge;

      triTrackTwo_pt     = triTrack_two.Pt();
      triTrackTwo_eta    = triTrack_two.Eta();
      triTrackTwo_phi    = triTrack_two.Phi();
      triTrackTwo_charge = lowTrack_charge + highTrack_charge + fourthTrack_charge;

      triTrackThree_pt     = triTrack_thr.Pt();
      triTrackThree_eta    = triTrack_thr.Eta();
      triTrackThree_phi    = triTrack_thr.Phi();
      triTrackThree_charge = highTrack_charge + thirdTrack_charge + fourthTrack_charge;

      triTrackFour_pt     = triTrack_fou.Pt();
      triTrackFour_eta    = triTrack_fou.Eta();
      triTrackFour_phi    = triTrack_fou.Phi();
      triTrackFour_charge = lowTrack_charge + thirdTrack_charge + fourthTrack_charge;


      triTrackOne_kkk = (highKaon_p4 + lowKaon_p4 + thirdKaon_p4).M();
      triTrackTwo_kkk = (highKaon_p4 + lowKaon_p4 + fourthKaon_p4).M();
      triTrackThree_kkk = (highKaon_p4 + fourthKaon_p4 + thirdKaon_p4).M();
      triTrackFour_kkk = (lowKaon_p4 + fourthKaon_p4 + thirdKaon_p4).M();

      triTrackOne_kkp = (highKaon_p4 + lowKaon_p4 + thirdPion_p4).M();
      triTrackTwo_kkp = (highKaon_p4 + lowKaon_p4 + fourthPion_p4).M();
      triTrackThree_kkp = (highKaon_p4 + fourthKaon_p4 + thirdPion_p4).M();
      triTrackFour_kkp = (lowKaon_p4 + fourthKaon_p4 + thirdPion_p4).M();

      triTrackOne_kpp = (highKaon_p4 + lowPion_p4 + thirdPion_p4).M();
      triTrackTwo_kpp = (highKaon_p4 + lowPion_p4 + fourthPion_p4).M();
      triTrackThree_kpp = (highKaon_p4 + fourthPion_p4 + thirdPion_p4).M();
      triTrackFour_kpp = (lowKaon_p4 + fourthPion_p4 + thirdPion_p4).M();

      triTrackOne_ppp = (highPion_p4 + lowPion_p4 + thirdPion_p4).M();
      triTrackTwo_ppp = (highPion_p4 + lowPion_p4 + fourthPion_p4).M();
      triTrackThree_ppp = (highPion_p4 + fourthPion_p4 + thirdPion_p4).M();
      triTrackFour_ppp = (lowPion_p4 + fourthPion_p4 + thirdPion_p4).M();


      if(IsMC_ || OnlyGen_)
      {

        gen_six_p4.SetPtEtaPhiM(-1.0,0.0,0.0,-0.01);
        gen_five_p4.SetPtEtaPhiM(-1.0,0.0,0.0,-0.01);
        gen_mumukk_p4.SetPtEtaPhiM(-1.0,0.,0.,-0.01);
        gen_jpsi_p4.SetPtEtaPhiM(-1.0,0.0,0.0,-0.01);
        gen_phi_p4.SetPtEtaPhiM(-1.0,0.0,0.0,-0.01);
        gen_highKaon_p4.SetPtEtaPhiM(-1.0,0.0,0.0,-0.01);
        gen_lowMuon_p4.SetPtEtaPhiM(-1.0,0.0,0.0,-0.01);
        gen_highMuon_p4.SetPtEtaPhiM(-1.0,0.0,0.0,-0.01);
        gen_lowKaon_p4.SetPtEtaPhiM(-1.0,0.0,0.0,-0.01);
        gen_thirdtrack_p4.SetPtEtaPhiM(-1.0,0.0,0.0,-0.01);

        gen_six_pdg         = 0.0;
        gen_mumukk_pdg      = 0.0;
        gen_phi_pdg         = 0.0;
        gen_jpsi_pdg        = 0.0;
        gen_five_pdg        = 0.0;

        gen_lowMuon_pdg     = 0.0;
        gen_highMuon_pdg    = 0.0;
        gen_highKaon_pdg    = 0.0;
        gen_lowKaon_pdg     = 0.0;
        gen_thirdtrack_pdg  = 0.0;

        gen_lowMuon_mompdg     = 0.0;
        gen_highMuon_mompdg    = 0.0;
        gen_highKaon_mompdg    = 0.0;
        gen_lowKaon_mompdg     = 0.0;

        gen_lowMuon_status     = 0.0;
        gen_highMuon_status    = 0.0;
        gen_highKaon_status    = 0.0;
        gen_lowKaon_status     = 0.0;

        gen_mumukk_prompt   = 0.0;
        gen_five_prompt     = 0.0;
        gen_phi_prompt      = 0.0;
        gen_jpsi_prompt     = 0.0;

        gen_mumukk_pt       = 0.0;
        gen_phi_pt          = 0.0;
        gen_jpsi_pt         = 0.0;
        gen_five_pt         = 0.0;

        gen_mumukk_p   = 0.0;
        gen_phi_p           = 0.0;
        gen_jpsi_p          = 0.0;
        gen_five_p          = 0.0;

        gen_mumukk_eta = 0.0;
        gen_phi_eta         = 0.0;
        gen_jpsi_eta        = 0.0;

        gen_mumukk_phi = 0.0;
        gen_phi_phi         = 0.0;
        gen_jpsi_phi        = 0.0;

        gen_lowMuon_pt     = 0.0;
        gen_highMuon_pt    = 0.0;
        gen_highKaon_pt    = 0.0;
        gen_lowKaon_pt     = 0.0;
        gen_thirdtrack_pt  = 0.0;

        gen_lowMuon_p     = 0.0;
        gen_highMuon_p    = 0.0;
        gen_highKaon_p    = 0.0;
        gen_lowKaon_p     = 0.0;
        gen_thirdtrack_p  = 0.0;

        gen_lowMuon_eta     = 0.0;
        gen_highMuon_eta    = 0.0;
        gen_highKaon_eta    = 0.0;
        gen_lowKaon_eta     = 0.0;
        gen_thirdtrack_eta  = 0.0;

        gen_lowMuon_phi     = 0.0;
        gen_highMuon_phi    = 0.0;
        gen_highKaon_phi    = 0.0;
        gen_lowKaon_phi     = 0.0;

        reco::GenParticleRef genhighMuon  = highMuon->genParticleRef();
        reco::GenParticleRef genlowMuon   = lowMuon->genParticleRef();

        const reco::GenParticle *genhighKaon,*genlowKaon,*gen_thirdtrack,*genFourthTrack,*swap;
        reco::GenParticleRef phiMomHigh, phiMomLow, jpsiMomHigh, jpsiMomLow;
        reco::GenParticleRef jpsiMom, phiMom, thirdMom;

        Double_t hasHighGen   = mumukk_cand->userFloat("hasHighGen");
        Double_t hasLowGen    = mumukk_cand->userFloat("hasLowGen");
        Double_t hasThirdGen  = six_cand.userFloat("hasThirdGen");
        Double_t hasFourthGen = six_cand.userFloat("hasFourthGen");

        if(hasHighGen>0.0)
          genhighKaon    = dynamic_cast <const reco::GenParticle *>(mumukk_cand->daughter("highKaonGen"));
        if(hasLowGen>0.0)
          genlowKaon     = dynamic_cast <const reco::GenParticle *>(mumukk_cand->daughter("lowKaonGen"));
        if(hasThirdGen>0.0)
          gen_thirdtrack  = dynamic_cast <const reco::GenParticle *>(six_cand.daughter("thirdTrackGen"));
        if(hasFourthGen>0.0)
          genFourthTrack = dynamic_cast <const reco::GenParticle *>(six_cand.daughter("fourthTrackGen"));

        if(hasFourthGen>0.0)
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
          gen_thirdtrack_p4.SetPtEtaPhiM(gen_thirdtrack->pt(),gen_thirdtrack->eta(),gen_thirdtrack->phi(),gen_thirdtrack->mass());

          if(gen_thirdtrack->numberOfMothers()>0)
            thirdMom  = gen_thirdtrack->motherRef();

          gen_thirdtrack_pdg     = (Double_t)gen_thirdtrack->pdgId();

          if(phiMomLow.isNonnull() && gen_thirdtrack->numberOfMothers()>0)
            gen_thirdtrack_mompdg  = thirdMom->pdgId();

          gen_thirdtrack_status  = (Double_t)gen_thirdtrack->status();
          gen_thirdtrack_pt      = (Double_t)gen_thirdtrack->pt();
          gen_thirdtrack_p       = (Double_t)gen_thirdtrack->p();
          gen_thirdtrack_eta     = (Double_t)gen_thirdtrack->eta();
          gen_thirdtrack_phi     = (Double_t)gen_thirdtrack->phi();


        }

        if(hasFourthGen>0.0)
        {
          gen_fourthtrackp4.SetPtEtaPhiM(genFourthTrack->pt(),genFourthTrack->eta(),genFourthTrack->phi(),genFourthTrack->mass());

          if(genFourthTrack->numberOfMothers()>0)
            thirdMom  = genFourthTrack->motherRef();

          gen_fourthtrackpdg     = (Double_t)genFourthTrack->pdgId();

          if(phiMomLow.isNonnull() && genFourthTrack->numberOfMothers()>0)
            gen_fourthtrackmompdg  = thirdMom->pdgId();

          gen_fourthtrackstatus  = (Double_t)genFourthTrack->status();
          gen_fourthtrackpt      = (Double_t)genFourthTrack->pt();
          gen_fourthtrackp       = (Double_t)genFourthTrack->p();
          gen_fourthtracketa     = (Double_t)genFourthTrack->eta();
          gen_fourthtrackphi     = (Double_t)genFourthTrack->phi();


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

            gen_mumukk_p4.SetPtEtaPhiM(jpsiMom->pt(),jpsiMom->eta(),jpsiMom->phi(),jpsiMom->mass());
            gen_mumukk_pdg = (Double_t) jpsiMom->pdgId();
            gen_mumukk_prompt = (Double_t) jpsiMom->isPromptDecayed();
            gen_mumukk_p = (Double_t) jpsiMom->p();
            gen_mumukk_pt = (Double_t) jpsiMom->pt();
            gen_mumukk_eta = (Double_t) jpsiMom->eta();
            gen_mumukk_phi = (Double_t) jpsiMom->phi();

          }

        }

      } //IsMC || onlyGen


      sixtracks_tree->Fill();

        // dimuontt candidates are sorted by vProb
    }

  }

}

// ------------ method called once each job just before starting event loop  ------------
void SixTracksRootupler::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void SixTracksRootupler::endJob() {}

// ------------ method called when starting to processes a run  ------------
void SixTracksRootupler::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void SixTracksRootupler::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void SixTracksRootupler::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void SixTracksRootupler::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SixTracksRootupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SixTracksRootupler);
