/*
Package:    DiMuonDiTrakFiveRootuplerFit
Class:      DiMuonDiTrakFiveRootuplerFit

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

class DiMuonDiTrakFiveRootuplerFit : public edm::EDAnalyzer {
public:
  explicit DiMuonDiTrakFiveRootuplerFit(const edm::ParameterSet&);
  ~DiMuonDiTrakFiveRootuplerFit() override;

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
  edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuonditrk_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> fivetrakpos_,fivetrakneg_,fivetrakneu_;
  edm::EDGetTokenT<reco::VertexCollection> pVertices_;
  edm::EDGetTokenT<edm::TriggerResults> triggers_;
  bool isMC_,onlyBest_,OnlyGen_ ;
  UInt_t motherpdgid_,phipdgid_,jpspdgid_;
  std::vector<std::string>  hlts_;
  std::vector<std::string>  hltFilters_;
  std::string treeName_;

  UInt_t run,event,numPrimaryVertices, trigger;

  TLorentzVector dimuonditrk_p4;
  TLorentzVector dimuon_p4;
  TLorentzVector ditrak_p4;
  TLorentzVector lowmuon_p4;
  TLorentzVector highmuon_p4;
  TLorentzVector highkaon_p4;
  TLorentzVector lowkaon_p4;

  TLorentzVector dimuonditrk_rf_p4;
  TLorentzVector dimuonditrk_rf_const_p4;
  TLorentzVector dimuon_rf_p4, dimuon_not_rf_p4;
  TLorentzVector ditrak_rf_p4, ditrak_not_rf_p4;
  TLorentzVector lowmuon_rf_p4;
  TLorentzVector highmuon_rf_p4;
  TLorentzVector kaonp_rf_p4;
  TLorentzVector kaonn_rf_p4;

  TLorentzVector fivetraks_pos_p4,fivetraks_pion_pos_p4;
  TLorentzVector dimuontrak_pos_p4,dimuontrak_pion_pos_p4;
  TLorentzVector fifthpion_pos_p4,fifthkaon_pos_p4;

  TLorentzVector fivetraks_neg_p4,fivetraks_pion_neg_p4;
  TLorentzVector dimuontrak_neg_p4,dimuontrak_pion_neg_p4;
  TLorentzVector fifthpion_neg_p4,fifthkaon_neg_p4;


  TLorentzVector fivetraks_neu_p4,fivetraks_pion_neu_p4;
  TLorentzVector dimuontrak_neu_p4,dimuontrak_pion_neu_p4;
  TLorentzVector fifthpion_neu_p4,fifthkaon_neu_p4;


  Int_t dimuonditrk_charge;

  UInt_t dimuon_triggerMatch, dimuon_triggerMatch_rf;

  Double_t dimuonditrk_vProb,  dimuonditrk_vChi2;
  Double_t dimuonditrk_rf_vProb, dimuonditrk_rf_vChi2, dimuonditrk_rf_nDof, dimuonditrk_rf_cosAlpha, dimuonditrk_rf_ctauPV, dimuonditrk_rf_ctauErrPV;
  Double_t dimuonditrk_rf_c_vProb, dimuonditrk_rf_c_vChi2, dimuonditrk_rf_c_nDof, dimuonditrk_rf_c_cosAlpha, dimuonditrk_rf_c_ctauPV, dimuonditrk_rf_c_ctauErrPV;
  Double_t dimuonditrk_pt, dimuonditrk_eta, dimuonditrk_phi, dimuonditrk_y, dimuonditrk_vx, dimuonditrk_vy, dimuonditrk_vz;

  Double_t pv_x, pv_y, pv_z;

  Double_t dimuonditrk_cosAlpha, dimuonditrk_ctauPV, dimuonditrk_ctauErrPV, dimuonditrk_countTksOfPV, dimuonditrk_vertexWeight;
  Double_t dimuonditrk_sumPTPV, dimuonditrk_mu1FromPV, dimuonditrk_mu2FromPV, dimuonditrk_tPFromPV, dimuonditrk_tMFromPV;
  Double_t dimuonditrk_mu1W, dimuonditrk_mu2W, dimuonditrk_tPW, dimuonditrk_tMW;

  Double_t dimuonditrk_cosAlphaDZ, dimuonditrk_ctauPVDZ, dimuonditrk_ctauErrPVDZ, dimuonditrk_countTksOfPVDZ, dimuonditrk_vertexWeightDZ;
  Double_t dimuonditrk_sumPTPVDZ, dimuonditrk_mu1FromPVDZ, dimuonditrk_mu2FromPVDZ, dimuonditrk_tPFromPVDZ, dimuonditrk_tMFromPVDZ;
  Double_t dimuonditrk_mu1DZW, dimuonditrk_mu2DZW, dimuonditrk_tPDZW, dimuonditrk_tMDZW;

  Double_t dimuonditrk_cosAlphaBS, dimuonditrk_ctauPVBS, dimuonditrk_ctauErrPVBS, dimuonditrk_countTksOfPVBS, dimuonditrk_vertexWeightBS;
  Double_t dimuonditrk_sumPTPVBS, dimuonditrk_mu1FromPVBS, dimuonditrk_mu2FromPVBS, dimuonditrk_tPFromPVBS, dimuonditrk_tMFromPVBS;
  Double_t dimuonditrk_mu1BSW, dimuonditrk_mu2BSW, dimuonditrk_tPBSW, dimuonditrk_tMBSW;

  Double_t dimuonditrk_dca_m1m2, dimuonditrk_dca_m1t1, dimuonditrk_dca_m1t2, dimuonditrk_dca_m2t1, dimuonditrk_dca_m2t2, dimuonditrk_dca_t1t2;
  Double_t dimuon_vProb, dimuon_vChi2, dimuon_DCA, dimuon_ctauPV, dimuon_ctauErrPV, dimuon_cosAlpha;

  Double_t gen_dimuonditrk_m,dimuonditrk_m,dimuon_m,dimuon_pt,ditrak_m,ditrak_pt;
  Double_t highkaon_pt,lowkaon_pt,highmuon_pt,lowmuon_pt,dimuonditrk_nDof,dimuonditrk_m_rf,dimuonditrk_m_rf_c,dimuonditrk_m_rf_d_c;

  Bool_t lowmuon_isLoose, lowmuon_isSoft, lowmuon_isMedium, lowmuon_isHighPt, lowmuon_isTight;
  Bool_t highmuon_isLoose, highmuon_isSoft, highmuon_isMedium, highmuon_isHighPt, highmuon_isTight;

  Bool_t lowmuon_isTracker, lowmuon_isGlobal, highmuon_isTracker, highmuon_isGlobal;
  UInt_t lowmuon_type, highmuon_type;

  Bool_t lowmuon_rf_isLoose, lowmuon_rf_isSoft, lowmuon_rf_isMedium, lowmuon_rf_isHighPt;
  Bool_t highmuon_rf_isLoose, highmuon_rf_isSoft, highmuon_rf_isMedium, highmuon_rf_isHighPt;

  Bool_t lowmuon_rf_isTracker, lowmuon_rf_isGlobal, highmuon_rf_isTracker, highmuon_rf_isGlobal;
  UInt_t lowmuon_rf_type, highmuon_rf_type;

  UInt_t lowmuon_NPixelHits, lowmuon_NStripHits, lowmuon_NTrackhits, lowmuon_NBPixHits, lowmuon_NPixLayers, lowmuon_NTraLayers, lowmuon_NStrLayers, lowmuon_NBPixLayers;
  UInt_t highmuon_NPixelHits, highmuon_NStripHits, highmuon_NTrackhits, highmuon_NBPixHits, highmuon_NPixLayers, highmuon_NTraLayers, highmuon_NStrLayers, highmuon_NBPixLayers;

  Double_t fivetraks_pos_kaon_m_rf, fivetraks_pos_pion_m_rf, fivetraks_pos_kaon_m, fivetraks_pos_pion_m;
  Double_t fivetraks_pos_vProb, fivetraks_pos_vChi2, fivetraks_pos_nDof, fivetraks_pos_charge;
  Double_t fivetraks_pos_cosAlpha, fivetraks_pos_ctauPV, fivetraks_pos_ctauErrPV;
  Double_t fifthtrak_pos_pt, fifthtrak_pos_phi, fifthtrak_pos_eta,fifthtrak_pos_y;
  Double_t fivetraks_pos_pion_trim,fivetraks_pos_kaon_trim,fivetraks_pos_eta;
  Double_t fivetraks_pos_pt,fivetraks_pos_phi,fivetraks_pos_y,fifthtrak_pos_charge;

  Double_t fivetraks_neu_kaon_m_rf, fivetraks_neu_pion_m_rf, fivetraks_neu_kaon_m, fivetraks_neu_pion_m;
  Double_t fivetraks_neu_vProb, fivetraks_neu_vChi2, fivetraks_neu_nDof, fivetraks_neu_charge;
  Double_t fivetraks_neu_cosAlpha, fivetraks_neu_ctauPV, fivetraks_neu_ctauErrPV;
  Double_t fifthtrak_neu_pt, fifthtrak_neu_phi, fifthtrak_neu_eta,fifthtrak_neu_y;
  Double_t fivetraks_neu_pion_trim,fivetraks_neu_kaon_trim,fivetraks_neu_eta;
  Double_t fivetraks_neu_pt,fivetraks_neu_phi,fivetraks_neu_y,fifthtrak_neu_charge;

  Double_t fivetraks_neg_kaon_m_rf, fivetraks_neg_pion_m_rf, fivetraks_neg_kaon_m, fivetraks_neg_pion_m;
  Double_t fivetraks_neg_vProb, fivetraks_neg_vChi2, fivetraks_neg_nDof, fivetraks_neg_charge;
  Double_t fivetraks_neg_cosAlpha, fivetraks_neg_ctauPV, fivetraks_neg_ctauErrPV;
  Double_t fifthtrak_neg_pt, fifthtrak_neg_phi, fifthtrak_neg_eta,fifthtrak_neg_y;
  Double_t fivetraks_neg_pion_trim,fivetraks_neg_kaon_trim,fivetraks_neg_eta;
  Double_t fivetraks_neg_pt,fivetraks_neg_phi,fivetraks_neg_y,fifthtrak_neg_charge;

  // Double_t highkaon_y, lowkaon_y, highmuon_y, lowmuon_y;

  Double_t dimuonditrk_refPK_mass, dimuonditrk_refKP_mass, dimuonditrk_refPP_mass, dimuonditrk_refPK_vChi2;
  Double_t dimuonditrk_refKP_vChi2, dimuonditrk_refPP_vChi2, dimuonditrk_refPK_nDof, dimuonditrk_refKP_nDof;
  Double_t dimuonditrk_refPP_nDof, dimuonditrk_refPK_vProb, dimuonditrk_refKP_vProb, dimuonditrk_refPP_vProb;

  Double_t highkaon_eta, lowkaon_eta, highmuon_eta, lowmuon_eta, highkaon_phi, lowkaon_phi, highmuon_phi, lowmuon_phi;
  Double_t highkaon_dz, lowkaon_dz, highmuon_dz, lowmuon_dz, highkaon_dxy, lowkaon_dxy, highmuon_dxy, lowmuon_dxy;
  // Double_t highkaon_etaError, lowkaon_etaError, highmuon_etaError, lowmuon_etaError, highkaon_phiError, lowkaon_phiError, highmuon_phiError, lowmuon_phiError;
  Int_t highkaon_NPixelHits, highkaon_NStripHits, highkaon_NTrackhits, highkaon_NBPixHits, highkaon_NPixLayers;
  Int_t highkaon_NTraLayers, highkaon_NStrLayers, highkaon_NBPixLayers, lowkaon_NPixelHits, lowkaon_NStripHits;
  Int_t lowkaon_NTrackhits, lowkaon_NBPixHits, lowkaon_NPixLayers, lowkaon_NTraLayers, lowkaon_NStrLayers, lowkaon_NBPixLayers;

  UInt_t highkaonMatch, lowkaonMatch,lowmuonMatch, highmuonMatch;

  Int_t dimuonditrk_rf_bindx;

  Int_t noXCandidates;

  Bool_t isBestCandidate;

  //MC
  TLorentzVector gen_dimuonditrk_p4, gen_jpsi_p4, gen_phi_p4;
  TLorentzVector gen_lowmuon_p4, gen_highmuon_p4, gen_highkaon_p4, gen_lowkaon_p4;

  Double_t gen_dimuonditrk_pdg, gen_phi_pdg, gen_jpsi_pdg;
  Double_t gen_lowmuon_pdg, gen_highmuon_pdg, gen_highkaon_pdg, gen_lowkaon_pdg;
  Double_t gen_lowmuon_mompdg, gen_highmuon_mompdg, gen_highkaon_mompdg, gen_lowkaon_mompdg;
  Double_t gen_lowmuon_status, gen_highmuon_status, gen_highkaon_status, gen_lowkaon_status;
  Double_t gen_dimuonditrk_prompt, gen_phi_prompt, gen_jpsi_prompt;
  Double_t gen_dimuonditrk_pt, gen_dimuonditrk_p, gen_dimuonditrk_eta;
  Double_t gen_phi_pt, gen_phi_p, gen_phi_eta;
  Double_t gen_jpsi_pt, gen_jpsi_p, gen_jpsi_eta;
  Double_t gen_lowmuon_phi, gen_highmuon_phi, gen_highkaon_phi, gen_lowkaon_phi;
  Double_t gen_dimuonditrk_phi, gen_phi_phi, gen_jpsi_phi;
  Double_t gen_lowmuon_p, gen_highmuon_p, gen_highkaon_p, gen_lowkaon_p;
  Double_t gen_lowmuon_pt, gen_highmuon_pt, gen_highkaon_pt, gen_lowkaon_pt;
  Double_t gen_lowmuon_eta, gen_highmuon_eta, gen_highkaon_eta, gen_lowkaon_eta;

  TTree* dimuonditrk_tree, *dimuonditrk_tree_rf;
  edm::EDGetTokenT< std::vector <reco::GenParticle> > genCands_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;
};

UInt_t DiMuonDiTrakFiveRootuplerFit::isTriggerMatched(pat::CompositeCandidate *diMuon_cand) {
  const pat::Muon* highmuon = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("highmuon"));
  const pat::Muon* lowmuon = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("lowmuon"));
  UInt_t matched = 0;  // if no list is given, is not matched

  // if matched a given trigger, set the bit, in the same order as listed
  for (unsigned int iTr = 0; iTr<hltFilters_.size(); iTr++ ) {
    // std::cout << hltFilters_[iTr] << std::endl;
    const pat::TriggerObjectStandAloneCollection mu1HLTMatches = highmuon->triggerObjectMatchesByFilter(hltFilters_[iTr]);
    const pat::TriggerObjectStandAloneCollection mu2HLTMatches = lowmuon->triggerObjectMatchesByFilter(hltFilters_[iTr]);
    if (!mu1HLTMatches.empty() && !mu2HLTMatches.empty()) matched += (1<<iTr);
    // if (!mu1HLTMatches.empty() && !mu2HLTMatches.empty()) std::cout << std::endl << hltFilters_[iTr] << std::endl;
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
DiMuonDiTrakFiveRootuplerFit::DiMuonDiTrakFiveRootuplerFit(const edm::ParameterSet& iConfig):
dimuonditrk_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("DiMuoDiTrak"))),
fivetrakpos_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("FiveTrakPos"))),
fivetrakneg_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("FiveTrakNeg"))),
fivetrakneu_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("FiveTrakNeu"))),
pVertices_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PrimaryVertices"))),
triggers_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
isMC_(iConfig.getParameter<bool>("isMC")),
onlyBest_(iConfig.getParameter<bool>("OnlyBest")),
OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
motherpdgid_(iConfig.getParameter<uint32_t>("Mother_pdg")),
phipdgid_(iConfig.getParameter<uint32_t>("jpsi_pdg")),
jpspdgid_(iConfig.getParameter<uint32_t>("Phi_pdg")),
hlts_(iConfig.getParameter<std::vector<std::string>>("HLTs")),
hltFilters_(iConfig.getParameter<std::vector<std::string>>("Filters")),
treeName_(iConfig.getParameter<std::string>("TreeName"))
{
  edm::Service<TFileService> fs;
  dimuonditrk_tree = fs->make<TTree>(treeName_.data(),"Tree of DiMuon and DiTrak");


  dimuonditrk_tree->Branch("run",                &run,                "run/I");
  dimuonditrk_tree->Branch("event",              &event,              "event/I");
  dimuonditrk_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");
  dimuonditrk_tree->Branch("trigger",            &trigger,            "trigger/I");

  if(!OnlyGen_)
  {
    dimuonditrk_tree->Branch("noXCandidates",      &noXCandidates,      "noXCandidates/I");

    //p4s
    dimuonditrk_tree->Branch("dimuonditrk_p4",   "TLorentzVector", &dimuonditrk_p4);
    dimuonditrk_tree->Branch("ditrak_p4",     "TLorentzVector", &ditrak_p4);
    dimuonditrk_tree->Branch("dimuon_p4",     "TLorentzVector", &dimuon_p4);
    dimuonditrk_tree->Branch("lowmuon_p4",   "TLorentzVector", &lowmuon_p4);
    dimuonditrk_tree->Branch("highmuon_p4",   "TLorentzVector", &highmuon_p4);
    dimuonditrk_tree->Branch("highkaon_p4",   "TLorentzVector", &highkaon_p4);
    dimuonditrk_tree->Branch("lowkaon_p4",   "TLorentzVector", &lowkaon_p4);
    dimuonditrk_tree->Branch("fivetraks_pos_p4",   "TLorentzVector", &fivetraks_pos_p4);
    dimuonditrk_tree->Branch("fivetraks_pion_pos_p4",   "TLorentzVector", &fivetraks_pion_pos_p4);
    dimuonditrk_tree->Branch("dimuontrak_pos_p4",   "TLorentzVector", &dimuontrak_pos_p4);
    dimuonditrk_tree->Branch("dimuontrak_pion_pos_p4",   "TLorentzVector", &dimuontrak_pion_pos_p4);
    dimuonditrk_tree->Branch("fifthkaon_pos_p4",   "TLorentzVector", &fifthkaon_pos_p4);
    dimuonditrk_tree->Branch("fifthpion_pos_p4",   "TLorentzVector", &fifthpion_pos_p4);
    dimuonditrk_tree->Branch("fivetraks_neu_p4",   "TLorentzVector", &fivetraks_neu_p4);
    dimuonditrk_tree->Branch("fivetraks_pion_neu_p4",   "TLorentzVector", &fivetraks_pion_neu_p4);
    dimuonditrk_tree->Branch("dimuontrak_neu_p4",   "TLorentzVector", &dimuontrak_neu_p4);
    dimuonditrk_tree->Branch("dimuontrak_pion_neu_p4",   "TLorentzVector", &dimuontrak_pion_neu_p4);
    dimuonditrk_tree->Branch("fifthkaon_neu_p4",   "TLorentzVector", &fifthkaon_neu_p4);
    dimuonditrk_tree->Branch("fifthpion_neu_p4",   "TLorentzVector", &fifthpion_neu_p4);
    dimuonditrk_tree->Branch("fivetraks_neg_p4",   "TLorentzVector", &fivetraks_neg_p4);
    dimuonditrk_tree->Branch("fivetraks_pion_neg_p4",   "TLorentzVector", &fivetraks_pion_neg_p4);
    dimuonditrk_tree->Branch("dimuontrak_neg_p4",   "TLorentzVector", &dimuontrak_neg_p4);
    dimuonditrk_tree->Branch("dimuontrak_pion_neg_p4",   "TLorentzVector", &dimuontrak_pion_neg_p4);
    dimuonditrk_tree->Branch("fifthkaon_neg_p4",   "TLorentzVector", &fifthkaon_neg_p4);
    dimuonditrk_tree->Branch("fifthpion_neg_p4",   "TLorentzVector", &fifthpion_neg_p4);

    //refitted p4s
    dimuonditrk_tree->Branch("dimuonditrk_rf_p4",   "TLorentzVector", &dimuonditrk_rf_p4);
    dimuonditrk_tree->Branch("ditrak_rf_p4",     "TLorentzVector", &ditrak_rf_p4);
    dimuonditrk_tree->Branch("dimuon_rf_p4",     "TLorentzVector", &dimuon_rf_p4);
    dimuonditrk_tree->Branch("lowmuon_rf_p4",   "TLorentzVector", &lowmuon_rf_p4);
    dimuonditrk_tree->Branch("highmuon_rf_p4",   "TLorentzVector", &highmuon_rf_p4);
    dimuonditrk_tree->Branch("kaonp_rf_p4",   "TLorentzVector", &kaonp_rf_p4);
    dimuonditrk_tree->Branch("kaonn_rf_p4",   "TLorentzVector", &kaonn_rf_p4);

    //kin
    dimuonditrk_tree->Branch("dimuonditrk_m",       &dimuonditrk_m,        "dimuonditrk_m/D");
    dimuonditrk_tree->Branch("dimuonditrk_m_rf",       &dimuonditrk_m_rf,        "dimuonditrk_m_rf/D");
    dimuonditrk_tree->Branch("dimuonditrk_m_rf_d_c",       &dimuonditrk_m_rf_d_c,        "dimuonditrk_m_rf_d_c/D");
    dimuonditrk_tree->Branch("dimuonditrk_m_rf_c",       &dimuonditrk_m_rf_c,        "dimuonditrk_m_rf_c/D");
    dimuonditrk_tree->Branch("dimuonditrk_pt",          &dimuonditrk_pt,          "dimuonditrk_pt/D");
    dimuonditrk_tree->Branch("dimuonditrk_eta",          &dimuonditrk_eta,          "dimuonditrk_eta/D");
    dimuonditrk_tree->Branch("dimuonditrk_phi",          &dimuonditrk_phi,          "dimuonditrk_phi/D");
    dimuonditrk_tree->Branch("dimuonditrk_y",          &dimuonditrk_y,          "dimuonditrk_y/D");

    dimuonditrk_tree->Branch("dimuonditrk_vx",          &dimuonditrk_vx,          "dimuonditrk_vx/D");
    dimuonditrk_tree->Branch("dimuonditrk_vy",          &dimuonditrk_vy,          "dimuonditrk_vy/D");
    dimuonditrk_tree->Branch("dimuonditrk_vz",          &dimuonditrk_vz,          "dimuonditrk_vz/D");

    dimuonditrk_tree->Branch("pv_x",          &pv_x,          "pv_x/D");
    dimuonditrk_tree->Branch("pv_y",          &pv_y,          "pv_y/D");
    dimuonditrk_tree->Branch("pv_z",          &pv_z,          "pv_z/D");

    dimuonditrk_tree->Branch("dimuon_m",       &dimuon_m,       "dimuon_m/D");
    dimuonditrk_tree->Branch("dimuon_pt",    &dimuon_pt,    "dimuon_pt/D");
    dimuonditrk_tree->Branch("ditrak_m",     &ditrak_m,     "ditrak_m/D");
    dimuonditrk_tree->Branch("ditrak_pt",       &ditrak_pt,        "ditrak_pt/D");
    dimuonditrk_tree->Branch("highkaon_pt",          &highkaon_pt,          "highkaon_pt/D");
    dimuonditrk_tree->Branch("lowkaon_pt",       &lowkaon_pt,       "lowkaon_pt/D");
    dimuonditrk_tree->Branch("highmuon_pt",    &highmuon_pt,    "highmuon_pt/D");
    dimuonditrk_tree->Branch("lowmuon_pt",     &lowmuon_pt,     "lowmuon_pt/D");

    dimuonditrk_tree->Branch("highkaon_eta",        &highkaon_eta,        "highkaon_eta/D");
    dimuonditrk_tree->Branch("lowkaon_eta",        &lowkaon_eta,        "lowkaon_eta/D");
    dimuonditrk_tree->Branch("highmuon_eta",        &highmuon_eta,        "highmuon_eta/D");
    dimuonditrk_tree->Branch("lowmuon_eta",        &lowmuon_eta,        "lowmuon_eta/D");

    dimuonditrk_tree->Branch("highkaon_phi",        &highkaon_phi,        "highkaon_phi/D");
    dimuonditrk_tree->Branch("lowkaon_phi",        &lowkaon_phi,        "lowkaon_phi/D");
    dimuonditrk_tree->Branch("highmuon_phi",        &highmuon_phi,        "highmuon_phi/D");
    dimuonditrk_tree->Branch("lowmuon_phi",        &lowmuon_phi,        "lowmuon_phi/D");

    // dimuonditrk_tree->Branch("highkaon_y",        &highkaon_y,        "highkaon_y/D");
    // dimuonditrk_tree->Branch("lowkaon_y",        &lowkaon_y,        "lowkaon_y/D");
    // dimuonditrk_tree->Branch("highmuon_y",        &highmuon_y,        "highmuon_y/D");
    // dimuonditrk_tree->Branch("lowmuon_y",        &lowmuon_y,        "lowmuon_y/D");

    //FiveTracks - Pos

    dimuonditrk_tree->Branch("fivetraks_pos_kaon_m",     &fivetraks_pos_kaon_m,     "fivetraks_pos_kaon_m/D");
    dimuonditrk_tree->Branch("fivetraks_pos_pion_m",     &fivetraks_pos_pion_m,     "fivetraks_pos_pion_m/D");
    dimuonditrk_tree->Branch("fivetraks_pos_kaon_trim",     &fivetraks_pos_kaon_trim,     "fivetraks_pos_kaon_trim/D");
    dimuonditrk_tree->Branch("fivetraks_pos_pion_trim",     &fivetraks_pos_pion_trim,     "fivetraks_pos_pion_trim/D");
    dimuonditrk_tree->Branch("fivetraks_pos_kaon_m_rf",     &fivetraks_pos_kaon_m_rf,     "fivetraks_pos_kaon_m_rf/D");
    dimuonditrk_tree->Branch("fivetraks_pos_pion_m_rf",     &fivetraks_pos_pion_m_rf,     "fivetraks_pos_pion_m_rf/D");
    dimuonditrk_tree->Branch("fivetraks_pos_vProb",      &fivetraks_pos_vProb,        "fivetraks_pos_vProb/D");
    dimuonditrk_tree->Branch("fivetraks_pos_vChi2",      &fivetraks_pos_vChi2,        "fivetraks_pos_vChi2/D");
    dimuonditrk_tree->Branch("fivetraks_pos_nDof",       &fivetraks_pos_nDof,         "fivetraks_pos_nDof/D");
    dimuonditrk_tree->Branch("fivetraks_pos_charge",     &fivetraks_pos_charge,       "fivetraks_pos_charge/I");
    dimuonditrk_tree->Branch("fivetraks_pos_cosAlpha",      &fivetraks_pos_cosAlpha,        "fivetraks_pos_cosAlpha/D");
    dimuonditrk_tree->Branch("fivetraks_pos_ctauPV",      &fivetraks_pos_ctauPV,        "fivetraks_pos_ctauPV/D");
    dimuonditrk_tree->Branch("fivetraks_pos_ctauErrPV",      &fivetraks_pos_ctauErrPV,        "fivetraks_pos_ctauErrPV/D");
    dimuonditrk_tree->Branch("fivetraks_pos_eta",       &fivetraks_pos_eta,        "fivetraks_pos_eta/D");
    dimuonditrk_tree->Branch("fivetraks_pos_pt",       &fivetraks_pos_pt,        "fivetraks_pos_pt/D");
    dimuonditrk_tree->Branch("fivetraks_pos_phi",       &fivetraks_pos_phi,        "fivetraks_pos_phi/D");
    dimuonditrk_tree->Branch("fivetraks_pos_y",       &fivetraks_pos_y,        "fivetraks_pos_y/D");

    dimuonditrk_tree->Branch("fifthtrak_pos_charge",       &fifthtrak_pos_charge,        "fifthtrak_pos_charge/D");
    dimuonditrk_tree->Branch("fifthtrak_pos_eta",       &fifthtrak_pos_eta,        "fifthtrak_pos_eta/D");
    dimuonditrk_tree->Branch("fifthtrak_pos_pt",       &fifthtrak_pos_pt,        "fifthtrak_pos_pt/D");
    dimuonditrk_tree->Branch("fifthtrak_pos_phi",       &fifthtrak_pos_phi,        "fifthtrak_pos_phi/D");
    dimuonditrk_tree->Branch("fifthtrak_pos_y",       &fifthtrak_pos_y,        "fifthtrak_pos_y/D");

    //FiveTracks - Neg

    dimuonditrk_tree->Branch("fivetraks_neg_kaon_m",     &fivetraks_neg_kaon_m,     "fivetraks_neg_kaon_m/D");
    dimuonditrk_tree->Branch("fivetraks_neg_pion_m",     &fivetraks_neg_pion_m,     "fivetraks_neg_pion_m/D");
    dimuonditrk_tree->Branch("fivetraks_neg_kaon_trim",     &fivetraks_neg_kaon_trim,     "fivetraks_neg_kaon_trim/D");
    dimuonditrk_tree->Branch("fivetraks_neg_pion_trim",     &fivetraks_neg_pion_trim,     "fivetraks_neg_pion_trim/D");
    dimuonditrk_tree->Branch("fivetraks_neg_kaon_m_rf",     &fivetraks_neg_kaon_m_rf,     "fivetraks_neg_kaon_m_rf/D");
    dimuonditrk_tree->Branch("fivetraks_neg_pion_m_rf",     &fivetraks_neg_pion_m_rf,     "fivetraks_neg_pion_m_rf/D");
    dimuonditrk_tree->Branch("fivetraks_neg_vProb",      &fivetraks_neg_vProb,        "fivetraks_neg_vProb/D");
    dimuonditrk_tree->Branch("fivetraks_neg_vChi2",      &fivetraks_neg_vChi2,        "fivetraks_neg_vChi2/D");
    dimuonditrk_tree->Branch("fivetraks_neg_nDof",       &fivetraks_neg_nDof,         "fivetraks_neg_nDof/D");
    dimuonditrk_tree->Branch("fivetraks_neg_charge",     &fivetraks_neg_charge,       "fivetraks_neg_charge/I");
    dimuonditrk_tree->Branch("fivetraks_neg_cosAlpha",      &fivetraks_neg_cosAlpha,        "fivetraks_neg_cosAlpha/D");
    dimuonditrk_tree->Branch("fivetraks_neg_ctauPV",      &fivetraks_neg_ctauPV,        "fivetraks_neg_ctauPV/D");
    dimuonditrk_tree->Branch("fivetraks_neg_ctauErrPV",      &fivetraks_neg_ctauErrPV,        "fivetraks_neg_ctauErrPV/D");
    dimuonditrk_tree->Branch("fivetraks_neg_eta",       &fivetraks_neg_eta,        "fivetraks_neg_eta/D");
    dimuonditrk_tree->Branch("fivetraks_neg_pt",       &fivetraks_neg_pt,        "fivetraks_neg_pt/D");
    dimuonditrk_tree->Branch("fivetraks_neg_phi",       &fivetraks_neg_phi,        "fivetraks_neg_phi/D");
    dimuonditrk_tree->Branch("fivetraks_neg_y",       &fivetraks_neg_y,        "fivetraks_neg_y/D");


    dimuonditrk_tree->Branch("fifthtrak_neg_charge",       &fifthtrak_neg_charge,        "fifthtrak_neg_charge/D");
    dimuonditrk_tree->Branch("fifthtrak_neg_eta",       &fifthtrak_neg_eta,        "fifthtrak_neg_eta/D");
    dimuonditrk_tree->Branch("fifthtrak_neg_pt",       &fifthtrak_neg_pt,        "fifthtrak_neg_pt/D");
    dimuonditrk_tree->Branch("fifthtrak_neg_phi",       &fifthtrak_neg_phi,        "fifthtrak_neg_phi/D");
    dimuonditrk_tree->Branch("fifthtrak_neg_y",       &fifthtrak_neg_y,        "fifthtrak_neg_y/D");

    //FiveTracks - Neu

    dimuonditrk_tree->Branch("fivetraks_neu_kaon_m",     &fivetraks_neu_kaon_m,     "fivetraks_neu_kaon_m/D");
    dimuonditrk_tree->Branch("fivetraks_neu_pion_m",     &fivetraks_neu_pion_m,     "fivetraks_neu_pion_m/D");
    dimuonditrk_tree->Branch("fivetraks_neu_kaon_trim",     &fivetraks_neu_kaon_trim,     "fivetraks_neu_kaon_trim/D");
    dimuonditrk_tree->Branch("fivetraks_neu_pion_trim",     &fivetraks_neu_pion_trim,     "fivetraks_neu_pion_trim/D");
    dimuonditrk_tree->Branch("fivetraks_neu_kaon_m_rf",     &fivetraks_neu_kaon_m_rf,     "fivetraks_neu_kaon_m_rf/D");
    dimuonditrk_tree->Branch("fivetraks_neu_pion_m_rf",     &fivetraks_neu_pion_m_rf,     "fivetraks_neu_pion_m_rf/D");
    dimuonditrk_tree->Branch("fivetraks_neu_vProb",      &fivetraks_neu_vProb,        "fivetraks_neu_vProb/D");
    dimuonditrk_tree->Branch("fivetraks_neu_vChi2",      &fivetraks_neu_vChi2,        "fivetraks_neu_vChi2/D");
    dimuonditrk_tree->Branch("fivetraks_neu_nDof",       &fivetraks_neu_nDof,         "fivetraks_neu_nDof/D");
    dimuonditrk_tree->Branch("fivetraks_neu_charge",     &fivetraks_neu_charge,       "fivetraks_neu_charge/I");
    dimuonditrk_tree->Branch("fivetraks_neu_cosAlpha",      &fivetraks_neu_cosAlpha,        "fivetraks_neu_cosAlpha/D");
    dimuonditrk_tree->Branch("fivetraks_neu_ctauPV",      &fivetraks_neu_ctauPV,        "fivetraks_neu_ctauPV/D");
    dimuonditrk_tree->Branch("fivetraks_neu_ctauErrPV",      &fivetraks_neu_ctauErrPV,        "fivetraks_neu_ctauErrPV/D");
    dimuonditrk_tree->Branch("fivetraks_neu_eta",       &fivetraks_neu_eta,        "fivetraks_neu_eta/D");
    dimuonditrk_tree->Branch("fivetraks_neu_pt",       &fivetraks_neu_pt,        "fivetraks_neu_pt/D");
    dimuonditrk_tree->Branch("fivetraks_neu_phi",       &fivetraks_neu_phi,        "fivetraks_neu_phi/D");
    dimuonditrk_tree->Branch("fivetraks_neu_y",       &fivetraks_neu_y,        "fivetraks_neu_y/D");


    dimuonditrk_tree->Branch("fifthtrak_neu_charge",       &fifthtrak_neu_charge,        "fifthtrak_neu_charge/D");
    dimuonditrk_tree->Branch("fifthtrak_neu_eta",       &fifthtrak_neu_eta,        "fifthtrak_neu_eta/D");
    dimuonditrk_tree->Branch("fifthtrak_neu_pt",       &fifthtrak_neu_pt,        "fifthtrak_neu_pt/D");
    dimuonditrk_tree->Branch("fifthtrak_neu_phi",       &fifthtrak_neu_phi,        "fifthtrak_neu_phi/D");
    dimuonditrk_tree->Branch("fifthtrak_neu_y",       &fifthtrak_neu_y,        "fifthtrak_neu_y/D");

    //Pion refits
    dimuonditrk_tree->Branch("dimuonditrk_refPK_mass", &dimuonditrk_refPK_mass, "dimuonditrk_refPK_mass/D");
    dimuonditrk_tree->Branch("dimuonditrk_refKP_mass", &dimuonditrk_refKP_mass, "dimuonditrk_refKP_mass/D");
    dimuonditrk_tree->Branch("dimuonditrk_refPP_mass", &dimuonditrk_refPP_mass, "dimuonditrk_refPP_mass/D");
    dimuonditrk_tree->Branch("dimuonditrk_refPK_vChi2", &dimuonditrk_refPK_vChi2, "dimuonditrk_refPK_vChi2/D");
    dimuonditrk_tree->Branch("dimuonditrk_refKP_vChi2", &dimuonditrk_refKP_vChi2, "dimuonditrk_refKP_vChi2/D");
    dimuonditrk_tree->Branch("dimuonditrk_refPP_vChi2", &dimuonditrk_refPP_vChi2, "dimuonditrk_refPP_vChi2/D");
    dimuonditrk_tree->Branch("dimuonditrk_refPK_nDof", &dimuonditrk_refPK_nDof, "dimuonditrk_refPK_nDof/D");
    dimuonditrk_tree->Branch("dimuonditrk_refKP_nDof", &dimuonditrk_refKP_nDof, "dimuonditrk_refKP_nDof/D");
    dimuonditrk_tree->Branch("dimuonditrk_refPP_nDof", &dimuonditrk_refPP_nDof, "dimuonditrk_refPP_nDof/D");
    dimuonditrk_tree->Branch("dimuonditrk_refPK_vProb", &dimuonditrk_refPK_vProb, "dimuonditrk_refPK_vProb/D");
    dimuonditrk_tree->Branch("dimuonditrk_refKP_vProb", &dimuonditrk_refKP_vProb, "dimuonditrk_refKP_vProb/D");
    dimuonditrk_tree->Branch("dimuonditrk_refPP_vProb", &dimuonditrk_refPP_vProb, "dimuonditrk_refPP_vProb/D");
    // dimuonditrk_tree->Branch("highkaon_trig_pt",          &highkaon_trig_pt,          "highkaon_trig_pt/D");
    // dimuonditrk_tree->Branch("lowkaon_trig_pt",       &lowkaon_trig_pt,       "lowkaon_trig_pt/D");
    // dimuonditrk_tree->Branch("highmuon_trig_pt",    &highmuon_trig_pt,    "highmuon_trig_pt/D");
    // dimuonditrk_tree->Branch("lowmuon_trig_pt",     &lowmuon_trig_pt,     "lowmuon_trig_pt/D");

    //2mu vertexing
    dimuonditrk_tree->Branch("dimuon_vProb",        &dimuon_vProb,        "dimuon_vProb/D");
    dimuonditrk_tree->Branch("dimuon_vChi2",       &dimuon_vChi2,        "dimuon_vChi2/D");
    dimuonditrk_tree->Branch("dimuon_DCA",          &dimuon_DCA,          "dimuon_DCA/D");
    dimuonditrk_tree->Branch("dimuon_ctauPV",       &dimuon_ctauPV,       "dimuon_ctauPV/D");
    dimuonditrk_tree->Branch("dimuon_ctauErrPV",    &dimuon_ctauErrPV,    "dimuon_ctauErrPV/D");
    dimuonditrk_tree->Branch("dimuon_cosAlpha",     &dimuon_cosAlpha,     "dimuon_cosAlpha/D");
    dimuonditrk_tree->Branch("dimuon_triggerMatch", &dimuon_triggerMatch, "dimuon_triggerMatch/I");

    //2mu+2Trk vertexing
    dimuonditrk_tree->Branch("dimuonditrk_vProb",      &dimuonditrk_vProb,        "dimuonditrk_vProb/D");
    dimuonditrk_tree->Branch("dimuonditrk_vChi2",      &dimuonditrk_vChi2,        "dimuonditrk_vChi2/D");
    dimuonditrk_tree->Branch("dimuonditrk_nDof",       &dimuonditrk_nDof,         "dimuonditrk_nDof/D");
    dimuonditrk_tree->Branch("dimuonditrk_charge",     &dimuonditrk_charge,       "dimuonditrk_charge/I");

    dimuonditrk_tree->Branch("dimuonditrk_cosAlpha",      &dimuonditrk_cosAlpha,        "dimuonditrk_cosAlpha/D");
    dimuonditrk_tree->Branch("dimuonditrk_ctauPV",      &dimuonditrk_ctauPV,        "dimuonditrk_ctauPV/D");
    dimuonditrk_tree->Branch("dimuonditrk_ctauErrPV",      &dimuonditrk_ctauErrPV,        "dimuonditrk_ctauErrPV/D");
    // dimuonditrk_tree->Branch("dimuonditrk_countTksOfPV",      &dimuonditrk_countTksOfPV,        "dimuonditrk_countTksOfPV/D");
    // dimuonditrk_tree->Branch("dimuonditrk_vertexWeight",      &dimuonditrk_vertexWeight,        "dimuonditrk_vertexWeight/D");
    // dimuonditrk_tree->Branch("dimuonditrk_sumPTPV",      &dimuonditrk_sumPTPV,        "dimuonditrk_sumPTPV/D");
    // dimuonditrk_tree->Branch("dimuonditrk_mu1FromPV",      &dimuonditrk_mu1FromPV,        "dimuonditrk_mu1FromPV/D");
    // dimuonditrk_tree->Branch("dimuonditrk_mu2FromPV",      &dimuonditrk_mu2FromPV,        "dimuonditrk_mu2FromPV/D");
    dimuonditrk_tree->Branch("dimuonditrk_tPFromPV",      &dimuonditrk_tPFromPV,        "dimuonditrk_tPFromPV/D");
    dimuonditrk_tree->Branch("dimuonditrk_tMFromPV",      &dimuonditrk_tMFromPV,        "dimuonditrk_tMFromPV/D");
    // dimuonditrk_tree->Branch("dimuonditrk_mu1W",      &dimuonditrk_mu1W,        "dimuonditrk_mu1W/D");
    // dimuonditrk_tree->Branch("dimuonditrk_mu2W",      &dimuonditrk_mu2W,        "dimuonditrk_mu2W/D");
    // dimuonditrk_tree->Branch("dimuonditrk_tPW",      &dimuonditrk_tPW,        "dimuonditrk_tPW/D");
    // dimuonditrk_tree->Branch("dimuonditrk_tMW",      &dimuonditrk_tMW,        "dimuonditrk_tMW/D");

    dimuonditrk_tree->Branch("dimuonditrk_cosAlphaDZ",      &dimuonditrk_cosAlphaDZ,        "dimuonditrk_cosAlphaDZ/D");
    dimuonditrk_tree->Branch("dimuonditrk_ctauPVDZ",      &dimuonditrk_ctauPVDZ,        "dimuonditrk_ctauPVDZ/D");
    dimuonditrk_tree->Branch("dimuonditrk_ctauErrPVDZ",      &dimuonditrk_ctauErrPVDZ,        "dimuonditrk_ctauErrPVDZ/D");
    // dimuonditrk_tree->Branch("dimuonditrk_countTksOfPVDZ",      &dimuonditrk_countTksOfPVDZ,        "dimuonditrk_countTksOfPVDZ/D");
    // dimuonditrk_tree->Branch("dimuonditrk_vertexWeightDZ",      &dimuonditrk_vertexWeightDZ,        "dimuonditrk_vertexWeightDZ/D");
    // dimuonditrk_tree->Branch("dimuonditrk_sumPTPVDZ",      &dimuonditrk_sumPTPVDZ,        "dimuonditrk_sumPTPVDZ/D");
    // dimuonditrk_tree->Branch("dimuonditrk_mu1FromPVDZ",      &dimuonditrk_mu1FromPVDZ,        "dimuonditrk_mu1FromPVDZ/D");
    // dimuonditrk_tree->Branch("dimuonditrk_mu2FromPVDZ",      &dimuonditrk_mu2FromPVDZ,        "dimuonditrk_mu2FromPVDZ/D");
    dimuonditrk_tree->Branch("dimuonditrk_tPFromPVDZ",      &dimuonditrk_tPFromPVDZ,        "dimuonditrk_tPFromPVDZ/D");
    dimuonditrk_tree->Branch("dimuonditrk_tMFromPVDZ",      &dimuonditrk_tMFromPVDZ,        "dimuonditrk_tMFromPVDZ/D");
    // dimuonditrk_tree->Branch("dimuonditrk_mu1DZW",      &dimuonditrk_mu1DZW,        "dimuonditrk_mu1DZW/D");
    // dimuonditrk_tree->Branch("dimuonditrk_mu2DZW",      &dimuonditrk_mu2DZW,        "dimuonditrk_mu2DZW/D");
    // dimuonditrk_tree->Branch("dimuonditrk_tPDZW",      &dimuonditrk_tPDZW,        "dimuonditrk_tPDZW/D");
    // dimuonditrk_tree->Branch("dimuonditrk_tMDZW",      &dimuonditrk_tMDZW,        "dimuonditrk_tMDZW/D");

    dimuonditrk_tree->Branch("dimuonditrk_cosAlphaBS",      &dimuonditrk_cosAlphaBS,        "dimuonditrk_cosAlphaBS/D");
    dimuonditrk_tree->Branch("dimuonditrk_ctauPVBS",      &dimuonditrk_ctauPVBS,        "dimuonditrk_ctauPVBS/D");
    dimuonditrk_tree->Branch("dimuonditrk_ctauErrPVBS",      &dimuonditrk_ctauErrPVBS,        "dimuonditrk_ctauErrPVBS/D");
    // dimuonditrk_tree->Branch("dimuonditrk_countTksOfPVBS",      &dimuonditrk_countTksOfPVBS,        "dimuonditrk_countTksOfPVBS/D");
    // dimuonditrk_tree->Branch("dimuonditrk_vertexWeightBS",      &dimuonditrk_vertexWeightBS,        "dimuonditrk_vertexWeightBS/D");
    // dimuonditrk_tree->Branch("dimuonditrk_sumPTPVBS",      &dimuonditrk_sumPTPVBS,        "dimuonditrk_sumPTPVBS/D");
    // dimuonditrk_tree->Branch("dimuonditrk_mu1FromPVBS",      &dimuonditrk_mu1FromPVBS,        "dimuonditrk_mu1FromPVBS/D");
    // dimuonditrk_tree->Branch("dimuonditrk_mu2FromPVBS",      &dimuonditrk_mu2FromPVBS,        "dimuonditrk_mu2FromPVBS/D");
    dimuonditrk_tree->Branch("dimuonditrk_tPFromPVBS",      &dimuonditrk_tPFromPVBS,        "dimuonditrk_tPFromPVBS/D");
    dimuonditrk_tree->Branch("dimuonditrk_tMFromPVBS",      &dimuonditrk_tMFromPVBS,        "dimuonditrk_tMFromPVBS/D");
    // dimuonditrk_tree->Branch("dimuonditrk_mu1BSW",      &dimuonditrk_mu1BSW,        "dimuonditrk_mu1BSW/D");
    // dimuonditrk_tree->Branch("dimuonditrk_mu2BSW",      &dimuonditrk_mu2BSW,        "dimuonditrk_mu2BSW/D");
    // dimuonditrk_tree->Branch("dimuonditrk_tPBSW",      &dimuonditrk_tPBSW,        "dimuonditrk_tPBSW/D");
    // dimuonditrk_tree->Branch("dimuonditrk_tMBSW",      &dimuonditrk_tMBSW,        "dimuonditrk_tMBSW/D");

    dimuonditrk_tree->Branch("dimuonditrk_dca_m1m2",      &dimuonditrk_vProb,        "dimuonditrk_dca_m1m2/D");
    dimuonditrk_tree->Branch("dimuonditrk_dca_m1t1",      &dimuonditrk_vProb,        "dimuonditrk_dca_m1t1/D");
    dimuonditrk_tree->Branch("dimuonditrk_dca_m1t2",      &dimuonditrk_vProb,        "dimuonditrk_dca_m1t2/D");
    dimuonditrk_tree->Branch("dimuonditrk_dca_m2t1",      &dimuonditrk_vProb,        "dimuonditrk_dca_m2t1/D");
    dimuonditrk_tree->Branch("dimuonditrk_dca_m2t2",      &dimuonditrk_vProb,        "dimuonditrk_dca_m2t2/D");
    dimuonditrk_tree->Branch("dimuonditrk_dca_t1t2",      &dimuonditrk_vProb,        "dimuonditrk_dca_t1t2/D");

    dimuonditrk_tree->Branch("dimuonditrk_rf_vProb",      &dimuonditrk_rf_vProb,        "dimuonditrk_rf_vProb/D");
    dimuonditrk_tree->Branch("dimuonditrk_rf_vChi2",      &dimuonditrk_rf_vChi2,        "dimuonditrk_rf_vChi2/D");
    dimuonditrk_tree->Branch("dimuonditrk_rf_nDof",       &dimuonditrk_rf_nDof,         "dimuonditrk_rf_nDof/D");
    dimuonditrk_tree->Branch("dimuonditrk_rf_cosAlpha",   &dimuonditrk_rf_cosAlpha,     "dimuonditrk_rf_cosAlpha/D");
    dimuonditrk_tree->Branch("dimuonditrk_rf_ctauPV",     &dimuonditrk_rf_ctauPV,       "dimuonditrk_rf_ctauPV/D");
    dimuonditrk_tree->Branch("dimuonditrk_rf_ctauErrPV",  &dimuonditrk_rf_ctauErrPV,    "dimuonditrk_rf_ctauErrPV/D");

    dimuonditrk_tree->Branch("dimuonditrk_rf_c_vProb",      &dimuonditrk_rf_c_vProb,        "dimuonditrk_rf_c_vProb/D");
    dimuonditrk_tree->Branch("dimuonditrk_rf_c_vChi2",      &dimuonditrk_rf_c_vChi2,        "dimuonditrk_rf_c_vChi2/D");
    dimuonditrk_tree->Branch("dimuonditrk_rf_c_nDof",       &dimuonditrk_rf_c_nDof,         "dimuonditrk_rf_c_nDof/D");
    dimuonditrk_tree->Branch("dimuonditrk_rf_c_cosAlpha",   &dimuonditrk_rf_c_cosAlpha,     "dimuonditrk_rf_c_cosAlpha/D");
    dimuonditrk_tree->Branch("dimuonditrk_rf_c_ctauPV",     &dimuonditrk_rf_c_ctauPV,       "dimuonditrk_rf_c_ctauPV/D");
    dimuonditrk_tree->Branch("dimuonditrk_rf_c_ctauErrPV",  &dimuonditrk_rf_c_ctauErrPV,    "dimuonditrk_rf_c_ctauErrPV/D");

    dimuonditrk_tree->Branch("highkaonMatch",     &highkaonMatch,       "highkaonMatch/I");
    dimuonditrk_tree->Branch("lowkaonMatch",     &lowkaonMatch,       "lowkaonMatch/I");
    dimuonditrk_tree->Branch("lowmuonMatch",     &lowmuonMatch,       "lowmuonMatch/I");
    dimuonditrk_tree->Branch("highmuonMatch",     &highmuonMatch,       "highmuonMatch/I");

    //Muon flags
    dimuonditrk_tree->Branch("lowmuon_isTight",        &lowmuon_isTight,        "lowmuon_isTight/O");
    dimuonditrk_tree->Branch("lowmuon_isLoose",        &lowmuon_isLoose,        "lowmuon_isLoose/O");
    dimuonditrk_tree->Branch("lowmuon_isSoft",        &lowmuon_isSoft,        "lowmuon_isSoft/O");
    dimuonditrk_tree->Branch("lowmuon_isMedium",        &lowmuon_isMedium,        "lowmuon_isMedium/O");
    dimuonditrk_tree->Branch("lowmuon_isHighPt",        &lowmuon_isHighPt,        "lowmuon_isHighPt/O");

    dimuonditrk_tree->Branch("lowmuon_isTracker",        &lowmuon_isTracker,        "lowmuon_isTracker/O");
    dimuonditrk_tree->Branch("lowmuon_isGlobal",        &lowmuon_isGlobal,        "lowmuon_isGlobal/O");

    dimuonditrk_tree->Branch("lowmuon_NPixelHits",        &lowmuon_NPixelHits,        "lowmuon_NPixelHits/I");
    dimuonditrk_tree->Branch("lowmuon_NStripHits",        &lowmuon_NStripHits,        "lowmuon_NStripHits/I");
    dimuonditrk_tree->Branch("lowmuon_NTrackhits",        &lowmuon_NTrackhits,        "lowmuon_NTrackhits/I");
    dimuonditrk_tree->Branch("lowmuon_NBPixHits",        &lowmuon_NBPixHits,        "lowmuon_NBPixHits/I");

    dimuonditrk_tree->Branch("lowmuon_NPixLayers",        &lowmuon_NPixLayers,        "lowmuon_NPixLayers/I");
    dimuonditrk_tree->Branch("lowmuon_NTraLayers",        &lowmuon_NTraLayers,        "lowmuon_NTraLayers/I");
    dimuonditrk_tree->Branch("lowmuon_NStrLayers",        &lowmuon_NStrLayers,        "lowmuon_NStrLayers/I");
    dimuonditrk_tree->Branch("lowmuon_NBPixLayers",        &lowmuon_NBPixLayers,        "lowmuon_NBPixLayers/I");

    dimuonditrk_tree->Branch("highmuon_isTight",        &highmuon_isTight,        "highmuon_isTight/O");
    dimuonditrk_tree->Branch("highmuon_isLoose",        &highmuon_isLoose,        "highmuon_isLoose/O");
    dimuonditrk_tree->Branch("highmuon_isSoft",        &highmuon_isSoft,        "highmuon_isSoft/O");
    dimuonditrk_tree->Branch("highmuon_isMedium",        &highmuon_isMedium,        "highmuon_isMedium/O");
    dimuonditrk_tree->Branch("highmuon_isHighPt",        &highmuon_isHighPt,        "highmuon_isHighPt/O");

    dimuonditrk_tree->Branch("highmuon_isTracker",        &highmuon_isTracker,        "highmuon_isTracker/O");
    dimuonditrk_tree->Branch("highmuon_isGlobal",        &highmuon_isGlobal,        "highmuon_isGlobal/O");

    dimuonditrk_tree->Branch("highmuon_NPixelHits",        &highmuon_NPixelHits,        "highmuon_NPixelHits/I");
    dimuonditrk_tree->Branch("highmuon_NStripHits",        &highmuon_NStripHits,        "highmuon_NStripHits/I");
    dimuonditrk_tree->Branch("highmuon_NTrackhits",        &highmuon_NTrackhits,        "highmuon_NTrackhits/I");
    dimuonditrk_tree->Branch("highmuon_NBPixHits",        &highmuon_NBPixHits,        "highmuon_NBPixHits/I");

    dimuonditrk_tree->Branch("highmuon_NPixLayers",        &highmuon_NPixLayers,        "highmuon_NPixLayers/I");
    dimuonditrk_tree->Branch("highmuon_NTraLayers",        &highmuon_NTraLayers,        "highmuon_NTraLayers/I");
    dimuonditrk_tree->Branch("highmuon_NStrLayers",        &highmuon_NStrLayers,        "highmuon_NStrLayers/I");
    dimuonditrk_tree->Branch("highmuon_NBPixLayers",        &highmuon_NBPixLayers,        "highmuon_NBPixLayers/I");

    dimuonditrk_tree->Branch("highkaon_eta",        &highkaon_eta,        "highkaon_eta/D");
    dimuonditrk_tree->Branch("highkaon_phi",        &highkaon_phi,        "highkaon_phi/D");
    dimuonditrk_tree->Branch("highkaon_dz",        &highkaon_dz,        "highkaon_dz/D");
    dimuonditrk_tree->Branch("highkaon_dxy",        &highkaon_dxy,        "highkaon_dxy/D");

    dimuonditrk_tree->Branch("lowkaon_eta",        &lowkaon_eta,        "lowkaon_eta/D");
    dimuonditrk_tree->Branch("lowkaon_phi",        &lowkaon_phi,        "lowkaon_phi/D");
    dimuonditrk_tree->Branch("lowkaon_dz",        &lowkaon_dz,        "lowkaon_dz/D");
    dimuonditrk_tree->Branch("lowkaon_dxy",        &lowkaon_dxy,        "lowkaon_dxy/D");

    dimuonditrk_tree->Branch("lowmuon_type",     &lowmuon_type,       "lowmuon_type/i");
    dimuonditrk_tree->Branch("highmuon_type",     &highmuon_type,       "highmuon_type/i");

    //Tracks Flags

    dimuonditrk_tree->Branch("highkaon_eta",        &highkaon_eta,        "highkaon_eta/D");
    dimuonditrk_tree->Branch("highkaon_phi",        &highkaon_phi,        "highkaon_phi/D");
    dimuonditrk_tree->Branch("highkaon_dz",        &highkaon_dz,        "highkaon_dz/D");
    dimuonditrk_tree->Branch("highkaon_dxy",        &highkaon_dxy,        "highkaon_dxy/D");

    dimuonditrk_tree->Branch("lowkaon_eta",        &lowkaon_eta,        "lowkaon_eta/D");
    dimuonditrk_tree->Branch("lowkaon_phi",        &lowkaon_phi,        "lowkaon_phi/D");
    dimuonditrk_tree->Branch("lowkaon_dz",        &lowkaon_dz,        "lowkaon_dz/D");
    dimuonditrk_tree->Branch("lowkaon_dxy",        &lowkaon_dxy,        "lowkaon_dxy/D");

    dimuonditrk_tree->Branch("highkaon_NPixelHits",        &highkaon_NPixelHits,        "highkaon_NPixelHits/I");
    dimuonditrk_tree->Branch("highkaon_NStripHits",        &highkaon_NStripHits,        "highkaon_NStripHits/I");
    dimuonditrk_tree->Branch("highkaon_NTrackhits",        &highkaon_NTrackhits,        "highkaon_NTrackhits/I");
    dimuonditrk_tree->Branch("highkaon_NBPixHits",        &highkaon_NBPixHits,        "highkaon_NBPixHits/I");

    dimuonditrk_tree->Branch("highkaon_NPixLayers",        &highkaon_NPixLayers,        "highkaon_NPixLayers/I");
    dimuonditrk_tree->Branch("highkaon_NTraLayers",        &highkaon_NTraLayers,        "highkaon_NTraLayers/I");
    dimuonditrk_tree->Branch("highkaon_NStrLayers",        &highkaon_NStrLayers,        "highkaon_NStrLayers/I");
    dimuonditrk_tree->Branch("highkaon_NBPixLayers",        &highkaon_NBPixLayers,        "highkaon_NBPixLayers/I");

    dimuonditrk_tree->Branch("lowkaon_NPixelHits",        &lowkaon_NPixelHits,        "lowkaon_NPixelHits/I");
    dimuonditrk_tree->Branch("lowkaon_NStripHits",        &lowkaon_NStripHits,        "lowkaon_NStripHits/I");
    dimuonditrk_tree->Branch("lowkaon_NTrackhits",        &lowkaon_NTrackhits,        "lowkaon_NTrackhits/I");
    dimuonditrk_tree->Branch("lowkaon_NBPixHits",        &lowkaon_NBPixHits,        "lowkaon_NBPixHits/I");

    dimuonditrk_tree->Branch("lowkaon_NPixLayers",        &lowkaon_NPixLayers,        "lowkaon_NPixLayers/I");
    dimuonditrk_tree->Branch("lowkaon_NTraLayers",        &lowkaon_NTraLayers,        "lowkaon_NTraLayers/I");
    dimuonditrk_tree->Branch("lowkaon_NStrLayers",        &lowkaon_NStrLayers,        "lowkaon_NStrLayers/I");
    dimuonditrk_tree->Branch("lowkaon_NBPixLayers",        &lowkaon_NBPixLayers,        "lowkaon_NBPixLayers/I");

  }
  int pdgid_ = 0;

  if (isMC_ || OnlyGen_) {

    dimuonditrk_tree->Branch("gen_dimuonditrk_p4", "TLorentzVector",  &gen_dimuonditrk_p4);
    dimuonditrk_tree->Branch("gen_jpsi_p4", "TLorentzVector",  &gen_jpsi_p4);
    dimuonditrk_tree->Branch("gen_phi_p4", "TLorentzVector",  &gen_phi_p4);

    dimuonditrk_tree->Branch("gen_highkaon_p4", "TLorentzVector",  &gen_highkaon_p4);
    dimuonditrk_tree->Branch("gen_lowmuon_p4",  "TLorentzVector",  &gen_lowmuon_p4);
    dimuonditrk_tree->Branch("gen_highmuon_p4",  "TLorentzVector",  &gen_highmuon_p4);
    dimuonditrk_tree->Branch("gen_lowkaon_p4",  "TLorentzVector",  &gen_lowkaon_p4);

    dimuonditrk_tree->Branch("gen_dimuonditrk_pdg",&gen_dimuonditrk_pdg,"gen_dimuonditrk_pdg/D");
    dimuonditrk_tree->Branch("gen_phi_pdg",&gen_phi_pdg,"gen_phi_pdg/D");

    dimuonditrk_tree->Branch("gen_lowmuon_pdg",&gen_lowmuon_pdg,"gen_lowmuon_pdg/D");
    dimuonditrk_tree->Branch("gen_highmuon_pdg",&gen_highmuon_pdg,"gen_highmuon_pdg/D");
    dimuonditrk_tree->Branch("gen_highkaon_pdg",&gen_highkaon_pdg,"gen_highkaon_pdg/D");
    dimuonditrk_tree->Branch("gen_lowkaon_pdg",&gen_lowkaon_pdg,"gen_lowkaon_pdg/D");

    dimuonditrk_tree->Branch("gen_lowmuon_mompdg",&gen_lowmuon_mompdg,"gen_lowmuon_mompdg/D");
    dimuonditrk_tree->Branch("gen_highmuon_mompdg",&gen_highmuon_mompdg,"gen_highmuon_mompdg/D");
    dimuonditrk_tree->Branch("gen_highkaon_mompdg",&gen_highkaon_mompdg,"gen_highkaon_mompdg/D");
    dimuonditrk_tree->Branch("gen_lowkaon_mompdg",&gen_lowkaon_mompdg,"gen_lowkaon_mompdg/D");

    dimuonditrk_tree->Branch("gen_lowmuon_status",&gen_lowmuon_status,"gen_lowmuon_status/D");
    dimuonditrk_tree->Branch("gen_highmuon_status",&gen_highmuon_status,"gen_highmuon_status/D");
    dimuonditrk_tree->Branch("gen_highkaon_status",&gen_highkaon_status,"gen_highkaon_status/D");
    dimuonditrk_tree->Branch("gen_lowkaon_status",&gen_lowkaon_status,"gen_lowkaon_status/D");

    dimuonditrk_tree->Branch("gen_lowmuon_p",&gen_lowmuon_p,"gen_lowmuon_p/D");
    dimuonditrk_tree->Branch("gen_highmuon_p",&gen_highmuon_p,"gen_highmuon_p/D");
    dimuonditrk_tree->Branch("gen_highkaon_p",&gen_highkaon_p,"gen_highkaon_p/D");
    dimuonditrk_tree->Branch("gen_lowkaon_p",&gen_lowkaon_p,"gen_lowkaon_p/D");

    dimuonditrk_tree->Branch("gen_lowmuon_pt",&gen_lowmuon_pt,"gen_lowmuon_pt/D");
    dimuonditrk_tree->Branch("gen_highmuon_pt",&gen_highmuon_pt,"gen_highmuon_pt/D");
    dimuonditrk_tree->Branch("gen_highkaon_pt",&gen_highkaon_pt,"gen_highkaon_pt/D");
    dimuonditrk_tree->Branch("gen_lowkaon_pt",&gen_lowkaon_pt,"gen_lowkaon_pt/D");

    dimuonditrk_tree->Branch("gen_lowmuon_eta",&gen_lowmuon_eta,"gen_lowmuon_eta/D");
    dimuonditrk_tree->Branch("gen_highmuon_eta",&gen_highmuon_eta,"gen_highmuon_eta/D");
    dimuonditrk_tree->Branch("gen_highkaon_eta",&gen_highkaon_eta,"gen_highkaon_eta/D");
    dimuonditrk_tree->Branch("gen_lowkaon_eta",&gen_lowkaon_eta,"gen_lowkaon_eta/D");

    dimuonditrk_tree->Branch("gen_lowmuon_phi",&gen_lowmuon_phi,"gen_lowmuon_phi/D");
    dimuonditrk_tree->Branch("gen_highmuon_phi",&gen_highmuon_phi,"gen_highmuon_phi/D");
    dimuonditrk_tree->Branch("gen_highkaon_phi",&gen_highkaon_phi,"gen_highkaon_phi/D");
    dimuonditrk_tree->Branch("gen_lowkaon_phi",&gen_lowkaon_phi,"gen_lowkaon_phi/D");

    dimuonditrk_tree->Branch("gen_dimuonditrk_prompt",&gen_dimuonditrk_prompt,"gen_dimuonditrk_prompt/D");
    dimuonditrk_tree->Branch("gen_phi_prompt",&gen_phi_prompt,"gen_phi_prompt/D");
    dimuonditrk_tree->Branch("gen_jpsi_prompt",&gen_jpsi_prompt,"gen_jpsi_prompt/D");

    dimuonditrk_tree->Branch("gen_dimuonditrk_pt",&gen_dimuonditrk_pt,"gen_dimuonditrk_pt/D");
    dimuonditrk_tree->Branch("gen_phi_pt",&gen_phi_pt,"phigen_phi_pt_pt/D");
    dimuonditrk_tree->Branch("gen_jpsi_pt",&gen_jpsi_pt,"gen_jpsi_pt/D");

    dimuonditrk_tree->Branch("gen_dimuonditrk_p",&gen_dimuonditrk_p,"gen_dimuonditrk_p/D");
    dimuonditrk_tree->Branch("gen_phi_p",&gen_phi_p,"phigen_phi_p_p/D");
    dimuonditrk_tree->Branch("gen_jpsi_p",&gen_jpsi_p,"gen_jpsi_p/D");

    dimuonditrk_tree->Branch("gen_dimuonditrk_eta",&gen_dimuonditrk_eta,"gen_dimuonditrk_eta/D");
    dimuonditrk_tree->Branch("gen_phi_eta",&gen_phi_eta,"gen_phi_eta/D");
    dimuonditrk_tree->Branch("gen_jpsi_eta",&gen_jpsi_eta,"gen_jpsi_eta/D");

    dimuonditrk_tree->Branch("gen_dimuonditrk_phi",&gen_dimuonditrk_phi,"gen_dimuonditrk_phi/D");
    dimuonditrk_tree->Branch("gen_phi_phi",&gen_phi_phi,"gen_phi_phi/D");
    dimuonditrk_tree->Branch("gen_jpsi_phi",&gen_jpsi_phi,"gen_jpsi_phi/D");
  }

  //Track flags


  dimuonditrk_tree->Branch("isBestCandidate",        &isBestCandidate,        "isBestCandidate/O");

  genCands_ = consumes< std::vector <reco::GenParticle> >((edm::InputTag)"prunedGenParticles");
  packCands_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");

}

DiMuonDiTrakFiveRootuplerFit::~DiMuonDiTrakFiveRootuplerFit() {}

//
// member functions
//

bool DiMuonDiTrakFiveRootuplerFit::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
  if (ancestor == particle ) return true;
  for (size_t i=0; i< particle->numberOfMothers(); i++) {
    if (isAncestor(ancestor, particle->mother(i))) return true;
  }
  return false;
}


// ------------ method called for each event  ------------
void DiMuonDiTrakFiveRootuplerFit::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //  using namespace edm;
  using namespace std;

  edm::Handle<std::vector<pat::CompositeCandidate>> dimuonditrk_cand_handle;
  iEvent.getByToken(dimuonditrk_, dimuonditrk_cand_handle);

  edm::Handle<std::vector<pat::CompositeCandidate>> fivepos_cand_handle;
  iEvent.getByToken(fivetrakpos_, fivepos_cand_handle);

  edm::Handle<std::vector<pat::CompositeCandidate>> fiveneg_cand_handle;
  iEvent.getByToken(fivetrakneg_, fiveneg_cand_handle);

  edm::Handle<std::vector<pat::CompositeCandidate>> fiveneu_cand_handle;
  iEvent.getByToken(fivetrakneu_, fiveneu_cand_handle);

  edm::Handle<std::vector<reco::Vertex >> primaryVertices_handle;
  iEvent.getByToken(pVertices_, primaryVertices_handle);

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( triggers_ , triggerResults_handle);

  numPrimaryVertices = primaryVertices_handle->size();
  run = iEvent.id().run();
  event = iEvent.id().event();

  reco::Vertex theBeamSpotV;

  // edm::Handle<reco::BeamSpot> theBeamSpot;
  // iEvent.getByToken(thebeamspot_,theBeamSpot);
  // reco::BeamSpot bs = *theBeamSpot;


  trigger = 0;

  if (triggerResults_handle.isValid()) {
    const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
    unsigned int NTRIGGERS = hlts_.size();

    for (unsigned int i = 0; i < NTRIGGERS; i++) {
      for (int version = 1; version < 20; version++) {
        std::stringstream ss;
        ss << hlts_[i] << "_v" << version;
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
  //   gen_highmuon_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
  //   foundit++;
  // }
  // if ( motherInPrunedCollection != nullptr && (d->pdgId() == -13 ) && isAncestor(aditrkdimu , motherInPrunedCollection) ) {
  //   gen_lowmuon_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
  //   foundit++;
  // }

  gen_dimuonditrk_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_phi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_lowmuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_highmuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_highkaon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_lowkaon_p4.SetPtEtaPhiM(0.,0.,0.,0.);


  std::map <unsigned int,const pat::CompositeCandidate*> fourToFiveMapPos,fourToFiveMapNeu,fourToFiveMapNeg;

  if(!OnlyGen_)
  {
    if (!dimuonditrk_cand_handle.isValid()) std::cout<< "No dimuontt information " << run << "," << event <<std::endl;
    if (!fivepos_cand_handle.isValid()) std::cout<< "No fivetrack pos information " << run << "," << event <<std::endl;
    if (!fiveneu_cand_handle.isValid()) std::cout<< "No fivetrack neu information " << run << "," << event <<std::endl;
    if (!fiveneg_cand_handle.isValid()) std::cout<< "No fivetrack neg information " << run << "," << event <<std::endl;

    if (dimuonditrk_cand_handle.isValid() && fivepos_cand_handle.isValid())
    {
      for (unsigned int i=0; i< fivepos_cand_handle->size(); i++)
      fourToFiveMapPos[(fivepos_cand_handle->at(i).userInt("index"))] = &(fivepos_cand_handle->at(i));
    }

    if (dimuonditrk_cand_handle.isValid() && fiveneg_cand_handle.isValid())
    {
      for (unsigned int i=0; i< fiveneg_cand_handle->size(); i++)
      fourToFiveMapNeg[(fiveneg_cand_handle->at(i).userInt("index"))] = &(fiveneg_cand_handle->at(i));
    }

    if (dimuonditrk_cand_handle.isValid() && fiveneu_cand_handle.isValid())
    {
      for (unsigned int i=0; i< fiveneu_cand_handle->size(); i++)
      fourToFiveMapNeu[(fiveneu_cand_handle->at(i).userInt("index"))] = &(fiveneu_cand_handle->at(i));
    }

  }

  // get rf information. Notice we are just keeping combinations with succesfull vertex fit
  if(!OnlyGen_)
  {
    if (!dimuonditrk_cand_handle.isValid()) std::cout<< "No dimuontt information " << run << "," << event <<std::endl;
    if (!fivepos_cand_handle.isValid()) std::cout<< "No fivetrack pos information " << run << "," << event <<std::endl;
    if (!fiveneu_cand_handle.isValid()) std::cout<< "No fivetrack neu information " << run << "," << event <<std::endl;
    if (!fiveneg_cand_handle.isValid()) std::cout<< "No fivetrack neg information " << run << "," << event <<std::endl;
    if (dimuonditrk_cand_handle.isValid()) {

      fivetraks_pos_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      fivetraks_pion_pos_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      dimuontrak_pos_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      dimuontrak_pion_pos_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      fifthpion_pos_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      fifthkaon_pos_p4.SetPtEtaPhiM(0.,0.,0.,0.);

      fivetraks_neu_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      fivetraks_pion_neu_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      dimuontrak_neu_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      dimuontrak_pion_neu_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      fifthpion_neu_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      fifthkaon_neu_p4.SetPtEtaPhiM(0.,0.,0.,0.);

      fivetraks_neg_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      fivetraks_pion_neg_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      dimuontrak_neg_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      dimuontrak_pion_neg_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      fifthpion_neg_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      fifthkaon_neg_p4.SetPtEtaPhiM(0.,0.,0.,0.);

      dimuonditrk_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      ditrak_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      lowmuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      highmuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      highkaon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      lowkaon_p4.SetPtEtaPhiM(0.,0.,0.,0.);

      dimuonditrk_rf_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      ditrak_rf_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      dimuon_rf_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      lowmuon_rf_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      highmuon_rf_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      // highkaon_rf_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      // lowkaon_rf_p4.SetPtEtaPhiM(0.,0.,0.,0.);

      pat::CompositeCandidate *dimuonditrk_rf_cand, dimuonditrk_cand, *dimuon_cand;
      pat::CompositeCandidate *ditrak_cand, *dimuon_cand_rf, *ditrak_cand_rf;
      const pat::CompositeCandidate *fivetrak_cand;
      noXCandidates = (Int_t)(dimuonditrk_cand_handle->size());

      for (unsigned int i=0; i< dimuonditrk_cand_handle->size(); i++){

        dimuonditrk_cand   = dimuonditrk_cand_handle->at(i);

        const reco::Vertex thePrimaryV = *(dimuonditrk_cand.userData<reco::Vertex>("bestPV"));

        dimuonditrk_vProb     = dimuonditrk_cand.userFloat("vProb");
        dimuonditrk_vChi2     = dimuonditrk_cand.userFloat("vChi2");
        dimuonditrk_nDof      = dimuonditrk_cand.userFloat("nDof");
        dimuonditrk_charge    = dimuonditrk_cand.charge();

        dimuonditrk_cosAlphaBS = dimuonditrk_cand.userFloat("cosAlphaBS");
        dimuonditrk_ctauPVBS = dimuonditrk_cand.userFloat("ctauPVBS");
        dimuonditrk_ctauErrPVBS = dimuonditrk_cand.userFloat("ctauErrPVBS");
        // dimuonditrk_countTksOfPVBS = dimuonditrk_cand.userFloat("countTksOfPVBS");
        // dimuonditrk_vertexWeightBS = dimuonditrk_cand.userFloat("vertexWeightBS");
        // dimuonditrk_sumPTPVBS = dimuonditrk_cand.userFloat("sumPTPVBS");
        // dimuonditrk_mu1FromPVBS = dimuonditrk_cand.userFloat("mu1FromPVBS");
        // dimuonditrk_mu2FromPVBS = dimuonditrk_cand.userFloat("mu2FromPVBS");
        dimuonditrk_tPFromPVBS = dimuonditrk_cand.userFloat("tPFromPVBS");
        dimuonditrk_tMFromPVBS = dimuonditrk_cand.userFloat("tMFromPVBS");
        // dimuonditrk_mu1BSW = dimuonditrk_cand.userFloat("mu1BSW");
        // dimuonditrk_mu2BSW = dimuonditrk_cand.userFloat("mu2BSW");
        // dimuonditrk_tPBSW = dimuonditrk_cand.userFloat("tPBSW");
        // dimuonditrk_tMBSW = dimuonditrk_cand.userFloat("tMBSW");

        dimuonditrk_cosAlpha = dimuonditrk_cand.userFloat("cosAlpha");
        dimuonditrk_ctauPV = dimuonditrk_cand.userFloat("ctauPV");
        dimuonditrk_ctauErrPV = dimuonditrk_cand.userFloat("ctauErrPV");
        // dimuonditrk_countTksOfPV = dimuonditrk_cand.userFloat("countTksOfPV");
        // dimuonditrk_vertexWeight = dimuonditrk_cand.userFloat("vertexWeight");
        // dimuonditrk_sumPTPV = dimuonditrk_cand.userFloat("sumPTPV");
        // dimuonditrk_mu1FromPV = dimuonditrk_cand.userFloat("mu1FromPV");
        // dimuonditrk_mu2FromPV = dimuonditrk_cand.userFloat("mu2FromPV");
        dimuonditrk_tPFromPV = dimuonditrk_cand.userFloat("tPFromPV");
        dimuonditrk_tMFromPV = dimuonditrk_cand.userFloat("tMFromPV");
        // dimuonditrk_mu1W = dimuonditrk_cand.userFloat("mu1W");
        // dimuonditrk_mu2W = dimuonditrk_cand.userFloat("mu2W");
        // dimuonditrk_tPW = dimuonditrk_cand.userFloat("tPW");
        // dimuonditrk_tMW = dimuonditrk_cand.userFloat("tMW");

        dimuonditrk_cosAlphaDZ = dimuonditrk_cand.userFloat("cosAlphaDZ");
        dimuonditrk_ctauPVDZ = dimuonditrk_cand.userFloat("ctauPVDZ");
        dimuonditrk_ctauErrPVDZ = dimuonditrk_cand.userFloat("ctauErrPVDZ");
        // dimuonditrk_countTksOfPVDZ = dimuonditrk_cand.userFloat("countTksOfPVDZ");
        // dimuonditrk_vertexWeightDZ = dimuonditrk_cand.userFloat("vertexWeightDZ");
        // dimuonditrk_sumPTPVDZ = dimuonditrk_cand.userFloat("sumPTPVDZ");
        // dimuonditrk_mu1FromPVDZ = dimuonditrk_cand.userFloat("mu1FromPVDZ");
        // dimuonditrk_mu2FromPVDZ = dimuonditrk_cand.userFloat("mu2FromPVDZ");
        dimuonditrk_tPFromPVDZ = dimuonditrk_cand.userFloat("tPFromPVDZ");
        dimuonditrk_tMFromPVDZ = dimuonditrk_cand.userFloat("tMFromPVDZ");
        // dimuonditrk_mu1DZW = dimuonditrk_cand.userFloat("mu1DZW");
        // dimuonditrk_mu1DZW = dimuonditrk_cand.userFloat("mu2DZW");
        // dimuonditrk_tPDZW = dimuonditrk_cand.userFloat("tPDZW");
        // dimuonditrk_tMDZW = dimuonditrk_cand.userFloat("tMDZW");

        dimuonditrk_dca_m1m2 = dimuonditrk_cand.userFloat("dca_m1m2");
        dimuonditrk_dca_m1t1 = dimuonditrk_cand.userFloat("dca_m1t1");
        dimuonditrk_dca_m1t2 = dimuonditrk_cand.userFloat("dca_m1t2");
        dimuonditrk_dca_m2t1 = dimuonditrk_cand.userFloat("dca_m2t1");
        dimuonditrk_dca_m2t2 = dimuonditrk_cand.userFloat("dca_m2t2");
        dimuonditrk_dca_t1t2 = dimuonditrk_cand.userFloat("dca_t1t2");

        highkaonMatch = dimuonditrk_cand.userInt("highkaonMatch");
        lowkaonMatch = dimuonditrk_cand.userInt("lowkaonMatch");

        //unref corresponding

        dimuon_cand = dynamic_cast <pat::CompositeCandidate *>(dimuonditrk_cand.daughter("dimuon"));
        ditrak_cand = dynamic_cast <pat::CompositeCandidate *>(dimuonditrk_cand.daughter("ditrak"));

        const pat::Muon *lowmuon, *highmuon;

        reco::Candidate::LorentzVector vhighmuon = dimuon_cand->daughter("highmuon")->p4();
        reco::Candidate::LorentzVector vlowmuon = dimuon_cand->daughter("lowmuon")->p4();

        // if (dimuon_cand->daughter("highmuon")->charge() < 0) {
        if (vhighmuon.pt() < vlowmuon.pt()) {
          vhighmuon = dimuon_cand->daughter("lowmuon")->p4();
          vlowmuon = dimuon_cand->daughter("highmuon")->p4();
          highmuon = dynamic_cast<const pat::Muon*>(dimuon_cand->daughter("highmuon"));
          lowmuon = dynamic_cast<const pat::Muon*>(dimuon_cand->daughter("lowmuon"));
        } else
        {
          lowmuon = dynamic_cast<const pat::Muon*>(dimuon_cand->daughter("highmuon"));
          highmuon = dynamic_cast<const pat::Muon*>(dimuon_cand->daughter("lowmuon"));
        }

        lowmuon_isTight    = lowmuon->isTightMuon(thePrimaryV);
        lowmuon_isLoose    = lowmuon->isLooseMuon();
        lowmuon_isSoft     = lowmuon->isSoftMuon(thePrimaryV);
        lowmuon_isMedium   = lowmuon->isMediumMuon();
        lowmuon_isHighPt   = lowmuon->isHighPtMuon(thePrimaryV);
        lowmuon_isTracker  = lowmuon->isTrackerMuon();
        lowmuon_isGlobal   = lowmuon->isGlobalMuon();
        lowmuon_NPixelHits = lowmuon->innerTrack()->hitPattern().numberOfValidPixelHits();
        lowmuon_NStripHits = lowmuon->innerTrack()->hitPattern().numberOfValidStripHits();
        lowmuon_NTrackhits = lowmuon->innerTrack()->hitPattern().numberOfValidTrackerHits();
        lowmuon_NBPixHits  = lowmuon->innerTrack()->hitPattern().numberOfValidStripHits();
        lowmuon_NPixLayers = lowmuon->innerTrack()->hitPattern().pixelLayersWithMeasurement();
        lowmuon_NTraLayers = lowmuon->innerTrack()->hitPattern().trackerLayersWithMeasurement();
        lowmuon_NStrLayers = lowmuon->innerTrack()->hitPattern().stripLayersWithMeasurement();
        lowmuon_NBPixLayers = lowmuon->innerTrack()->hitPattern().pixelBarrelLayersWithMeasurement();

        highmuon_isTight    = highmuon->isTightMuon(thePrimaryV);
        highmuon_isLoose    = highmuon->isLooseMuon();
        highmuon_isSoft     = highmuon->isSoftMuon(thePrimaryV);
        highmuon_isMedium   = highmuon->isMediumMuon();
        highmuon_isHighPt   = highmuon->isHighPtMuon(thePrimaryV);
        highmuon_isTracker  = highmuon->isTrackerMuon();
        highmuon_isGlobal   = highmuon->isGlobalMuon();
        highmuon_NPixelHits = highmuon->innerTrack()->hitPattern().numberOfValidPixelHits();
        highmuon_NStripHits = highmuon->innerTrack()->hitPattern().numberOfValidStripHits();
        highmuon_NTrackhits = highmuon->innerTrack()->hitPattern().numberOfValidTrackerHits();
        highmuon_NBPixHits  = highmuon->innerTrack()->hitPattern().numberOfValidStripHits();
        highmuon_NPixLayers = highmuon->innerTrack()->hitPattern().pixelLayersWithMeasurement();
        highmuon_NTraLayers = highmuon->innerTrack()->hitPattern().trackerLayersWithMeasurement();
        highmuon_NStrLayers = highmuon->innerTrack()->hitPattern().stripLayersWithMeasurement();
        highmuon_NBPixLayers = highmuon->innerTrack()->hitPattern().pixelBarrelLayersWithMeasurement();

        highmuon_eta  = highmuon->innerTrack()->eta();
        highmuon_phi  = highmuon->innerTrack()->phi();
        highmuon_dz  = highmuon->innerTrack()->dz();
        highmuon_dxy  = highmuon->innerTrack()->dxy();
        lowmuon_eta  = lowmuon->innerTrack()->eta();
        lowmuon_phi  = lowmuon->innerTrack()->phi();
        lowmuon_dz  = lowmuon->innerTrack()->dz();
        lowmuon_dxy  = lowmuon->innerTrack()->dxy();

        lowmuon_type       = lowmuon->type();
        highmuon_type       = highmuon->type();

        lowmuon_p4.SetPtEtaPhiM(vhighmuon.pt(), vhighmuon.eta(), vhighmuon.phi(), vhighmuon.mass());
        highmuon_p4.SetPtEtaPhiM(vlowmuon.pt(), vlowmuon.eta(), vlowmuon.phi(), vlowmuon.mass());

        reco::Candidate::LorentzVector kP = ditrak_cand->daughter("highTrak")->p4();
        reco::Candidate::LorentzVector kM = ditrak_cand->daughter("lowTrak")->p4();

        pat::PackedCandidate *highkaon  = dynamic_cast <pat::PackedCandidate *>(ditrak_cand->daughter("highTrak"));
        pat::PackedCandidate *lowkaon  = dynamic_cast <pat::PackedCandidate *>(ditrak_cand->daughter("lowTrak"));

        highkaon_p4.SetPtEtaPhiM(kP.pt(), kP.eta(), kP.phi(), kP.mass());
        lowkaon_p4.SetPtEtaPhiM(kM.pt(), kM.eta(), kM.phi(), kM.mass());

        highkaon_pt     = std::max(kP.pt(),kM.pt());
        lowkaon_pt      = -std::max(-kP.pt(),-kM.pt());
        highmuon_pt     = std::max(vlowmuon.pt(),vhighmuon.pt());
        lowmuon_pt      = -std::max(-vlowmuon.pt(),-vhighmuon.pt());

        lowkaon_NPixelHits = lowkaon->bestTrack()->hitPattern().numberOfValidPixelHits();
        lowkaon_NStripHits = lowkaon->bestTrack()->hitPattern().numberOfValidStripHits();
        lowkaon_NTrackhits = lowkaon->bestTrack()->hitPattern().numberOfValidTrackerHits();
        lowkaon_NBPixHits  = lowkaon->bestTrack()->hitPattern().numberOfValidStripHits();
        lowkaon_NPixLayers = lowkaon->bestTrack()->hitPattern().pixelLayersWithMeasurement();
        lowkaon_NTraLayers = lowkaon->bestTrack()->hitPattern().trackerLayersWithMeasurement();
        lowkaon_NStrLayers = lowkaon->bestTrack()->hitPattern().stripLayersWithMeasurement();
        lowkaon_NBPixLayers = lowkaon->bestTrack()->hitPattern().pixelBarrelLayersWithMeasurement();
        highkaon_NPixelHits = highkaon->bestTrack()->hitPattern().numberOfValidPixelHits();
        highkaon_NStripHits = highkaon->bestTrack()->hitPattern().numberOfValidStripHits();
        highkaon_NTrackhits = highkaon->bestTrack()->hitPattern().numberOfValidTrackerHits();
        highkaon_NBPixHits  = highkaon->bestTrack()->hitPattern().numberOfValidStripHits();
        highkaon_NPixLayers = highkaon->bestTrack()->hitPattern().pixelLayersWithMeasurement();
        highkaon_NTraLayers = highkaon->bestTrack()->hitPattern().trackerLayersWithMeasurement();
        highkaon_NStrLayers = highkaon->bestTrack()->hitPattern().stripLayersWithMeasurement();
        highkaon_NBPixLayers = highkaon->bestTrack()->hitPattern().pixelBarrelLayersWithMeasurement();

        highkaon_eta  = highkaon->bestTrack()->eta();
        highkaon_phi  = highkaon->bestTrack()->phi();
        highkaon_dz   = highkaon->bestTrack()->dz();
        highkaon_dxy  = highkaon->bestTrack()->dxy();
        lowkaon_eta   = lowkaon->bestTrack()->eta();
        lowkaon_phi   = lowkaon->bestTrack()->phi();
        lowkaon_dz    = lowkaon->bestTrack()->dz();
        lowkaon_dxy   = lowkaon->bestTrack()->dxy();

        //double kmass = 0.4936770;
        dimuonditrk_p4.SetPtEtaPhiM(dimuonditrk_cand.pt(),dimuonditrk_cand.eta(),dimuonditrk_cand.phi(),dimuonditrk_cand.mass());
        dimuon_p4.SetPtEtaPhiM(dimuon_cand->pt(),dimuon_cand->eta(),dimuon_cand->phi(),dimuon_cand->mass());
        ditrak_p4.SetPtEtaPhiM(ditrak_cand->pt(), ditrak_cand->eta(), ditrak_cand->phi(), ditrak_cand->mass());

        dimuon_vProb        = dimuon_cand->userFloat("vProb");
        dimuon_vChi2        = dimuon_cand->userFloat("vNChi2");
        dimuon_DCA          = dimuon_cand->userFloat("DCA");
        dimuon_ctauPV       = dimuon_cand->userFloat("ppdlPV");
        dimuon_ctauErrPV    = dimuon_cand->userFloat("ppdlErrPV");
        dimuon_cosAlpha     = dimuon_cand->userFloat("cosAlpha");

        lowmuonMatch   = dimuon_cand->userInt("highmuonTMatch");
        highmuonMatch    = dimuon_cand->userInt("lowmuonTMatch");
        dimuon_triggerMatch = -std::max(-lowmuonMatch,highmuonMatch);//DiMuonDiTrakFiveRootuplerFit::isTriggerMatched(dimuon_cand);

        dimuonditrk_rf_p4.SetPtEtaPhiM(0.0,0.0,0.0,3.9);
        dimuonditrk_rf_const_p4.SetPtEtaPhiM(0.0,0.0,0.0,3.9);
        dimuon_rf_p4.SetPtEtaPhiM(0.0,0.0,0.0,2.4);
        ditrak_rf_p4.SetPtEtaPhiM(0.0,0.0,0.0,0.5);

        lowmuon_rf_p4.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
        highmuon_rf_p4.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
        kaonp_rf_p4.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
        kaonn_rf_p4.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
        dimuonditrk_m_rf_c = 3.9;
        dimuonditrk_m_rf_d_c = 3.9;

        if (dimuonditrk_cand.userFloat("has_ref") >= 0)
        {
          dimuonditrk_rf_cand = dynamic_cast <pat::CompositeCandidate *>(dimuonditrk_cand.daughter("ref_cand"));
          dimuonditrk_rf_p4.SetPtEtaPhiM(dimuonditrk_rf_cand->pt(),dimuonditrk_rf_cand->eta(),dimuonditrk_rf_cand->phi(),dimuonditrk_rf_cand->mass());
          dimuon_rf_p4.SetPtEtaPhiM(dimuonditrk_rf_cand->daughter("dimuon")->pt(),dimuonditrk_rf_cand->daughter("dimuon")->eta(),
          dimuonditrk_rf_cand->daughter("dimuon")->phi(),dimuonditrk_rf_cand->daughter("dimuon")->mass());
          ditrak_rf_p4.SetPtEtaPhiM(dimuonditrk_rf_cand->daughter("ditrak")->pt(),dimuonditrk_rf_cand->daughter("ditrak")->eta(),
          dimuonditrk_rf_cand->daughter("ditrak")->phi(),dimuonditrk_rf_cand->daughter("ditrak")->mass());

          dimuon_cand_rf = dynamic_cast <pat::CompositeCandidate *>(dimuonditrk_rf_cand->daughter("dimuon"));
          ditrak_cand_rf = dynamic_cast <pat::CompositeCandidate *>(dimuonditrk_rf_cand->daughter("ditrak"));

          vhighmuon = dimuon_cand_rf->daughter("highmuon")->p4();
          vlowmuon = dimuon_cand_rf->daughter("lowmuon")->p4();

          lowmuon_rf_p4.SetPtEtaPhiM(vhighmuon.pt(), vhighmuon.eta(), vhighmuon.phi(), vhighmuon.mass());
          highmuon_rf_p4.SetPtEtaPhiM(vlowmuon.pt(), vlowmuon.eta(), vlowmuon.phi(), vlowmuon.mass());

          kP = ditrak_cand_rf->daughter("highTrak")->p4();
          kM = ditrak_cand_rf->daughter("lowTrak")->p4();

          kaonp_rf_p4.SetPtEtaPhiM(kP.pt(), kP.eta(), kP.phi(), kP.mass());
          kaonn_rf_p4.SetPtEtaPhiM(kM.pt(), kM.eta(), kM.phi(), kM.mass());

          dimuonditrk_m_rf_c= dimuonditrk_rf_cand->mass();

          dimuonditrk_rf_vProb = dimuonditrk_cand.userFloat("vProb_ref");
          dimuonditrk_rf_vChi2 = dimuonditrk_cand.userFloat("vChi2_ref");
          dimuonditrk_rf_nDof = dimuonditrk_cand.userFloat("nDof_ref");
          dimuonditrk_rf_cosAlpha = dimuonditrk_cand.userFloat("cosAlpha_ref");
          dimuonditrk_rf_ctauPV = dimuonditrk_cand.userFloat("ctauPV_ref");
          dimuonditrk_rf_ctauErrPV = dimuonditrk_cand.userFloat("ctauErrPV_ref");


        }

        if (dimuonditrk_cand.userFloat("has_const_ref") >= 0)
        {
          dimuonditrk_rf_cand = dynamic_cast <pat::CompositeCandidate *>(dimuonditrk_cand.daughter("ref_const_cand"));
          dimuonditrk_rf_const_p4.SetPtEtaPhiM(dimuonditrk_rf_cand->pt(),dimuonditrk_rf_cand->eta(),dimuonditrk_rf_cand->phi(),dimuonditrk_rf_cand->mass());

          dimuonditrk_m_rf_d_c= dimuonditrk_rf_cand->mass();

          dimuonditrk_rf_vProb = dimuonditrk_cand.userFloat("vProb_const_ref");
          dimuonditrk_rf_vChi2 = dimuonditrk_cand.userFloat("vChi2_const_ref");
          dimuonditrk_rf_nDof = dimuonditrk_cand.userFloat("nDof_const_ref");
          dimuonditrk_rf_cosAlpha = dimuonditrk_cand.userFloat("cosAlpha_const_ref");
          dimuonditrk_rf_ctauPV = dimuonditrk_cand.userFloat("ctauPV_const_ref");
          dimuonditrk_rf_ctauErrPV = dimuonditrk_cand.userFloat("ctauErrPV_const_ref");

        }

        dimuonditrk_m    = dimuonditrk_cand.mass();
        dimuonditrk_m_rf = dimuonditrk_cand.userFloat("mass_rf");
        dimuonditrk_pt   = dimuonditrk_cand.pt();

        dimuonditrk_eta   = dimuonditrk_cand.eta();
        dimuonditrk_phi   = dimuonditrk_cand.phi();
        dimuonditrk_y     = dimuonditrk_cand.y();
        dimuonditrk_vx    = dimuonditrk_cand.userFloat("vtxX");
        dimuonditrk_vy    = dimuonditrk_cand.userFloat("vtxY");
        dimuonditrk_vz    = dimuonditrk_cand.userFloat("vtxZ");

        dimuon_m         = dimuon_cand->mass();
        dimuon_pt        = dimuon_cand->pt();
        ditrak_m         = ditrak_cand->mass();
        ditrak_pt        = ditrak_cand->pt();

        dimuonditrk_refPK_mass  = dimuonditrk_cand.userFloat("massPKRefit");
        dimuonditrk_refKP_mass  = dimuonditrk_cand.userFloat("massKPRefit");
        dimuonditrk_refPP_mass  = dimuonditrk_cand.userFloat("massPPRefit");

        dimuonditrk_refPK_vChi2 = dimuonditrk_cand.userFloat("vChi2PKRefit");
        dimuonditrk_refKP_vChi2 = dimuonditrk_cand.userFloat("vChi2KPRefit");
        dimuonditrk_refPP_vChi2 = dimuonditrk_cand.userFloat("vChi2PPRefit");

        dimuonditrk_refPK_nDof  = dimuonditrk_cand.userFloat("nDofPKRefit");
        dimuonditrk_refKP_nDof  = dimuonditrk_cand.userFloat("nDofKPRefit");
        dimuonditrk_refPP_nDof  = dimuonditrk_cand.userFloat("nDofPPRefit");

        dimuonditrk_refPK_vProb = dimuonditrk_cand.userFloat("vProbPKRefit");
        dimuonditrk_refKP_vProb = dimuonditrk_cand.userFloat("vProbKPRefit");
        dimuonditrk_refPP_vProb = dimuonditrk_cand.userFloat("vProbPPRefit");

        fivetraks_pos_kaon_m    = -1.0;
        fivetraks_pos_pion_m    = -1.0;
        fivetraks_pos_kaon_trim    = -1.0;
        fivetraks_pos_pion_trim    = -1.0;
        fivetraks_pos_kaon_m_rf    = -1.0;
        fivetraks_pos_pion_m_rf    = -1.0;
        fivetraks_pos_vProb    = -1.0;
        fivetraks_pos_vChi2    = -1.0;
        fivetraks_pos_nDof    = -1.0;
        fivetraks_pos_charge    = -1.0;
        fivetraks_pos_cosAlpha    = -1.0;
        fivetraks_pos_ctauPV    = -1.0;
        fivetraks_pos_ctauErrPV    = -1.0;
        fivetraks_pos_eta    = -1.0;
        fivetraks_pos_pt    = -1.0;
        fivetraks_pos_phi    = -1.0;
        fivetraks_pos_y    = -1.0;
        fifthtrak_pos_charge = -1.0;
        fifthtrak_pos_eta = -1.0;
        fifthtrak_pos_pt = -1.0;
        fifthtrak_pos_phi = -1.0;
        fifthtrak_pos_y = -1.0;

        if((fourToFiveMapPos.find(i))!=fourToFiveMapPos.end() && fivepos_cand_handle.isValid() )
        {

          fivetrak_cand = fourToFiveMapPos[(i)];

          // dimuontrak_cand = dynamic_cast <const pat::CompositeCandidate *>(fivetrak_cand->daughter("dimuontrak"));
          const pat::CompositeCandidate* fivetrakpion_cand = dynamic_cast <const pat::CompositeCandidate* >(fivetrak_cand->daughter("withpion"));
          // dimuontrakpion_cand = dynamic_cast <const pat::CompositeCandidate *>(fivetrakpion_cand->daughter("dimuontrak"));

          fivetraks_pos_p4.SetPtEtaPhiM(fivetrak_cand->pt(), fivetrak_cand->eta(), fivetrak_cand->phi(), fivetrak_cand->mass());
          fivetraks_pion_pos_p4.SetPtEtaPhiM(fivetrakpion_cand->pt(), fivetrakpion_cand->eta(), fivetrakpion_cand->phi(), fivetrakpion_cand->mass());
          // dimuontrak_pos_p4.SetPtEtaPhiM(fivetrakpion_cand->pt(), fivetrakpion_cand->eta(), fivetrakpion_cand->phi(), fivetrakpion_cand->mass());
          // dimuontrak_pos_p4.SetPtEtaPhiM(dimuontrakpion_cand->pt(), dimuontrakpion_cand->eta(), dimuontrakpion_cand->phi(), dimuontrakpion_cand->mass());

          reco::Candidate::LorentzVector fifth = fivetrak_cand->daughter("fifth")->p4();
          reco::Candidate::LorentzVector fifthpion = fivetrakpion_cand->daughter("fifth")->p4();
          fifthkaon_pos_p4.SetPtEtaPhiM(fifth.pt(), fifth.eta(), fifth.phi(), fifth.mass());
          fifthpion_pos_p4.SetPtEtaPhiM(fifthpion.pt(), fifthpion.eta(), fifthpion.phi(), fifthpion.mass());
          dimuontrak_pos_p4.SetPtEtaPhiM(fifthkaon_pos_p4.Pt() + dimuon_cand->pt(),fifthkaon_pos_p4.Eta() + dimuon_cand->eta(),fifthkaon_pos_p4.Phi() + dimuon_cand->phi(),fifthkaon_pos_p4.M() + dimuon_cand->mass());
          dimuontrak_pion_pos_p4.SetPtEtaPhiM(fifthpion_pos_p4.Pt() + dimuon_cand->pt(),fifthpion_pos_p4.Eta() + dimuon_cand->eta(),fifthpion_pos_p4.Phi() + dimuon_cand->phi(),fifthpion_pos_p4.M() + dimuon_cand->mass());

          fivetraks_pos_kaon_m    = fivetrak_cand->mass();
          fivetraks_pos_pion_m    = fivetrak_cand->userFloat("mass_pion");
          fivetraks_pos_kaon_trim    = dimuontrak_pos_p4.M();
          fivetraks_pos_pion_trim    = dimuontrak_pion_pos_p4.M();
          fivetraks_pos_kaon_m_rf    = fivetrak_cand->userFloat("mass_kaon_rf");
          fivetraks_pos_pion_m_rf    = fivetrak_cand->userFloat("mass_pion_rf");
          fivetraks_pos_vProb    = fivetrak_cand->userFloat("vProb");
          fivetraks_pos_vChi2    = fivetrak_cand->userFloat("vChi2");
          fivetraks_pos_nDof    = fivetrak_cand->userFloat("nDof");
          fivetraks_pos_charge    = fivetrak_cand->mass();
          fivetraks_pos_cosAlpha    = fivetrak_cand->userFloat("cosAlpha");
          fivetraks_pos_ctauPV    = fivetrak_cand->userFloat("ctauPV");
          fivetraks_pos_ctauErrPV    = fivetrak_cand->userFloat("ctauErrPV");
          fivetraks_pos_eta    = fivetrak_cand->eta();
          fivetraks_pos_pt    = fivetrak_cand->pt();
          fivetraks_pos_phi    = fivetrak_cand->phi();
          fivetraks_pos_y    = fivetrak_cand->y();

          fifthtrak_pos_charge = fivetrak_cand->daughter("fifth")->charge();
          fifthtrak_pos_eta = fifth.Eta();
          fifthtrak_pos_pt = fifth.Pt();
          fifthtrak_pos_phi = fifth.Phi();
          fifthtrak_pos_y = fifth.Rapidity();
        }

        fivetraks_neu_kaon_m    = -1.0;
        fivetraks_neu_pion_m    = -1.0;
        fivetraks_neu_kaon_trim    = -1.0;
        fivetraks_neu_pion_trim    = -1.0;
        fivetraks_neu_kaon_m_rf    = -1.0;
        fivetraks_neu_pion_m_rf    = -1.0;
        fivetraks_neu_vProb    = -1.0;
        fivetraks_neu_vChi2    = -1.0;
        fivetraks_neu_nDof    = -1.0;
        fivetraks_neu_charge    = -1.0;
        fivetraks_neu_cosAlpha    = -1.0;
        fivetraks_neu_ctauPV    = -1.0;
        fivetraks_neu_ctauErrPV    = -1.0;
        fivetraks_neu_eta    = -1.0;
        fivetraks_neu_pt    = -1.0;
        fivetraks_neu_phi    = -1.0;
        fivetraks_neu_y    = -1.0;
        fifthtrak_neu_charge = -1.0;
        fifthtrak_neu_eta = -1.0;
        fifthtrak_neu_pt = -1.0;
        fifthtrak_neu_phi = -1.0;
        fifthtrak_neu_y = -1.0;

        if((fourToFiveMapNeu.find(i))!=fourToFiveMapNeu.end() && fiveneu_cand_handle.isValid() )
        {

          fivetrak_cand = fourToFiveMapNeu[(i)];

          // dimuontrak_cand = dynamic_cast <const pat::CompositeCandidate *>(fivetrak_cand->daughter("dimuontrak"));
          const pat::CompositeCandidate* fivetrakpion_cand = dynamic_cast <const pat::CompositeCandidate *>(fivetrak_cand->daughter("withpion"));
          // dimuontrakpion_cand = dynamic_cast <const pat::CompositeCandidate *>(fivetrakpion_cand->daughter("dimuontrak"));

          fivetraks_neu_p4.SetPtEtaPhiM(fivetrak_cand->pt(), fivetrak_cand->eta(), fivetrak_cand->phi(), fivetrak_cand->mass());
          fivetraks_pion_neu_p4.SetPtEtaPhiM(fivetrakpion_cand->pt(), fivetrakpion_cand->eta(), fivetrakpion_cand->phi(), fivetrakpion_cand->mass());
          // dimuontrak_neu_p4.SetPtEtaPhiM(fivetrakpion_cand->pt(), fivetrakpion_cand->eta(), fivetrakpion_cand->phi(), fivetrakpion_cand->mass());
          // dimuontrak_neu_p4.SetPtEtaPhiM(dimuontrakpion_cand->pt(), dimuontrakpion_cand->eta(), dimuontrakpion_cand->phi(), dimuontrakpion_cand->mass());

          reco::Candidate::LorentzVector fifth = fivetrak_cand->daughter("fifth")->p4();
          reco::Candidate::LorentzVector fifthpion = fivetrakpion_cand->daughter("fifth")->p4();
          fifthkaon_neu_p4.SetPtEtaPhiM(fifth.pt(), fifth.eta(), fifth.phi(), fifth.mass());
          fifthpion_neu_p4.SetPtEtaPhiM(fifthpion.pt(), fifthpion.eta(), fifthpion.phi(), fifthpion.mass());
          dimuontrak_neu_p4.SetPtEtaPhiM(fifthkaon_neu_p4.Pt() + dimuon_cand->pt(),fifthkaon_neu_p4.Eta() + dimuon_cand->eta(),fifthkaon_neu_p4.Phi() + dimuon_cand->phi(),fifthkaon_neu_p4.M() + dimuon_cand->mass());
          dimuontrak_pion_neu_p4.SetPtEtaPhiM(fifthkaon_neu_p4.Pt() + dimuon_cand->pt(),fifthkaon_neu_p4.Eta() + dimuon_cand->eta(),fifthkaon_neu_p4.Phi() + dimuon_cand->phi(),fifthkaon_neu_p4.M() + dimuon_cand->mass());

          fivetraks_neu_kaon_m    = fivetrak_cand->mass();
          fivetraks_neu_pion_m    = fivetrak_cand->userFloat("mass_pion");
          fivetraks_neu_kaon_trim    = dimuontrak_neu_p4.M();
          fivetraks_neu_pion_trim    = dimuontrak_pion_neu_p4.M();
          fivetraks_neu_kaon_m_rf    = fivetrak_cand->userFloat("mass_kaon_rf");
          fivetraks_neu_pion_m_rf    = fivetrak_cand->userFloat("mass_pion_rf");
          fivetraks_neu_vProb    = fivetrak_cand->userFloat("vProb");
          fivetraks_neu_vChi2    = fivetrak_cand->userFloat("vChi2");
          fivetraks_neu_nDof    = fivetrak_cand->userFloat("nDof");
          fivetraks_neu_charge    = fivetrak_cand->mass();
          fivetraks_neu_cosAlpha    = fivetrak_cand->userFloat("cosAlpha");
          fivetraks_neu_ctauPV    = fivetrak_cand->userFloat("ctauPV");
          fivetraks_neu_ctauErrPV    = fivetrak_cand->userFloat("ctauErrPV");
          fivetraks_neu_eta    = fivetrak_cand->eta();
          fivetraks_neu_pt    = fivetrak_cand->pt();
          fivetraks_neu_phi    = fivetrak_cand->phi();
          fivetraks_neu_y    = fivetrak_cand->y();

          fifthtrak_neu_charge = fivetrak_cand->daughter("fifth")->charge();
          fifthtrak_neu_eta = fifth.Eta();
          fifthtrak_neu_pt = fifth.Pt();
          fifthtrak_neu_phi = fifth.Phi();
          fifthtrak_neu_y = fifth.Rapidity();
        }

        fivetraks_neg_kaon_m    = -1.0;
        fivetraks_neg_pion_m    = -1.0;
        fivetraks_neg_kaon_trim    = -1.0;
        fivetraks_neg_pion_trim    = -1.0;
        fivetraks_neg_kaon_m_rf    = -1.0;
        fivetraks_neg_pion_m_rf    = -1.0;
        fivetraks_neg_vProb    = -1.0;
        fivetraks_neg_vChi2    = -1.0;
        fivetraks_neg_nDof    = -1.0;
        fivetraks_neg_charge    = -1.0;
        fivetraks_neg_cosAlpha    = -1.0;
        fivetraks_neg_ctauPV    = -1.0;
        fivetraks_neg_ctauErrPV    = -1.0;
        fivetraks_neg_eta    = -1.0;
        fivetraks_neg_pt    = -1.0;
        fivetraks_neg_phi    = -1.0;
        fivetraks_neg_y    = -1.0;
        fifthtrak_neg_charge = -1.0;
        fifthtrak_neg_eta = -1.0;
        fifthtrak_neg_pt = -1.0;
        fifthtrak_neg_phi = -1.0;
        fifthtrak_neg_y = -1.0;

        if((fourToFiveMapNeg.find(i))!=fourToFiveMapNeg.end() && fiveneg_cand_handle.isValid())
        {

          fivetrak_cand = fourToFiveMapNeg[(i)];

          // dimuontrak_cand = dynamic_cast <const pat::CompositeCandidate *>(fivetrak_cand->daughter("dimuontrak"));
          const pat::CompositeCandidate* fivetrakpion_cand = dynamic_cast <const pat::CompositeCandidate *>(fivetrak_cand->daughter("withpion"));
          // dimuontrakpion_cand = dynamic_cast <const pat::CompositeCandidate *>(fivetrakpion_cand->daughter("dimuontrak"));

          fivetraks_neg_p4.SetPtEtaPhiM(fivetrak_cand->pt(), fivetrak_cand->eta(), fivetrak_cand->phi(), fivetrak_cand->mass());
          fivetraks_pion_neg_p4.SetPtEtaPhiM(fivetrakpion_cand->pt(), fivetrakpion_cand->eta(), fivetrakpion_cand->phi(), fivetrakpion_cand->mass());
          // dimuontrak_neg_p4.SetPtEtaPhiM(fivetrakpion_cand->pt(), fivetrakpion_cand->eta(), fivetrakpion_cand->phi(), fivetrakpion_cand->mass());
          // dimuontrak_neg_p4.SetPtEtaPhiM(dimuontrakpion_cand->pt(), dimuontrakpion_cand->eta(), dimuontrakpion_cand->phi(), dimuontrakpion_cand->mass());

          reco::Candidate::LorentzVector fifth = fivetrak_cand->daughter("fifth")->p4();
          reco::Candidate::LorentzVector fifthpion = fivetrakpion_cand->daughter("fifth")->p4();
          fifthkaon_neg_p4.SetPtEtaPhiM(fifth.pt(), fifth.eta(), fifth.phi(), fifth.mass());
          fifthkaon_neg_p4.SetPtEtaPhiM(fifthpion.pt(), fifthpion.eta(), fifthpion.phi(), fifthpion.mass());
          dimuontrak_neg_p4.SetPtEtaPhiM(fifthkaon_neg_p4.Pt() + dimuon_cand->pt(),fifthkaon_neg_p4.Eta() + dimuon_cand->eta(),fifthkaon_neg_p4.Phi() + dimuon_cand->phi(),fifthkaon_neg_p4.M() + dimuon_cand->mass());
          dimuontrak_pion_neg_p4.SetPtEtaPhiM(fifthpion_neg_p4.Pt() + dimuon_cand->pt(),fifthpion_neg_p4.Eta() + dimuon_cand->eta(),fifthpion_neg_p4.Phi() + dimuon_cand->phi(),fifthpion_neg_p4.M() + dimuon_cand->mass());

          fivetraks_neg_kaon_m    = fivetrak_cand->mass();
          fivetraks_neg_pion_m    = fivetrak_cand->userFloat("mass_pion");
          fivetraks_neg_kaon_trim    = dimuontrak_neg_p4.M();
          fivetraks_neg_pion_trim    = dimuontrak_pion_neg_p4.M();
          fivetraks_neg_kaon_m_rf    = fivetrak_cand->userFloat("mass_kaon_rf");
          fivetraks_neg_pion_m_rf    = fivetrak_cand->userFloat("mass_pion_rf");
          fivetraks_neg_vProb    = fivetrak_cand->userFloat("vProb");
          fivetraks_neg_vChi2    = fivetrak_cand->userFloat("vChi2");
          fivetraks_neg_nDof    = fivetrak_cand->userFloat("nDof");
          fivetraks_neg_charge    = fivetrak_cand->mass();
          fivetraks_neg_cosAlpha    = fivetrak_cand->userFloat("cosAlpha");
          fivetraks_neg_ctauPV    = fivetrak_cand->userFloat("ctauPV");
          fivetraks_neg_ctauErrPV    = fivetrak_cand->userFloat("ctauErrPV");
          fivetraks_neg_eta    = fivetrak_cand->eta();
          fivetraks_neg_pt    = fivetrak_cand->pt();
          fivetraks_neg_phi    = fivetrak_cand->phi();
          fivetraks_neg_y    = fivetrak_cand->y();

          fifthtrak_neg_charge = fivetrak_cand->daughter("fifth")->charge();
          fifthtrak_neg_eta = fifth.Eta();
          fifthtrak_neg_pt = fifth.Pt();
          fifthtrak_neg_phi = fifth.Phi();
          fifthtrak_neg_y = fifth.Rapidity();
        }



        if(isMC_ || OnlyGen_)
        {

          gen_dimuonditrk_p4.SetPtEtaPhiM(-1.0,0.0,0.0,3.9);
          gen_jpsi_p4.SetPtEtaPhiM(-1.0,0.0,0.0,2.5);
          gen_phi_p4.SetPtEtaPhiM(-1.0,0.0,0.0,0.9);
          gen_highkaon_p4.SetPtEtaPhiM(-1.0,0.0,0.0,0.2);
          gen_lowmuon_p4.SetPtEtaPhiM(-1.0,0.0,0.0,0.0);
          gen_highmuon_p4.SetPtEtaPhiM(-1.0,0.0,0.0,0.0);
          gen_lowkaon_p4.SetPtEtaPhiM(-1.0,0.0,0.0,0.2);

          gen_dimuonditrk_pdg = 0.0;
          gen_phi_pdg         = 0.0;

          gen_lowmuon_pdg     = 0.0;
          gen_highmuon_pdg    = 0.0;
          gen_highkaon_pdg    = 0.0;
          gen_lowkaon_pdg     = 0.0;

          gen_lowmuon_mompdg     = 0.0;
          gen_highmuon_mompdg    = 0.0;
          gen_highkaon_mompdg    = 0.0;
          gen_lowkaon_mompdg     = 0.0;

          gen_lowmuon_status     = 0.0;
          gen_highmuon_status    = 0.0;
          gen_highkaon_status    = 0.0;
          gen_lowkaon_status     = 0.0;

          gen_dimuonditrk_prompt = 0.0;
          gen_phi_prompt      = 0.0;
          gen_jpsi_prompt     = 0.0;

          gen_dimuonditrk_pt  = 0.0;
          gen_phi_pt          = 0.0;
          gen_jpsi_pt         = 0.0;

          gen_dimuonditrk_p   = 0.0;
          gen_phi_p           = 0.0;
          gen_jpsi_p          = 0.0;

          gen_dimuonditrk_eta = 0.0;
          gen_phi_eta         = 0.0;
          gen_jpsi_eta        = 0.0;

          gen_dimuonditrk_phi = 0.0;
          gen_phi_phi         = 0.0;
          gen_jpsi_phi        = 0.0;

          gen_lowmuon_pt     = 0.0;
          gen_highmuon_pt    = 0.0;
          gen_highkaon_pt    = 0.0;
          gen_lowkaon_pt     = 0.0;

          gen_lowmuon_p     = 0.0;
          gen_highmuon_p    = 0.0;
          gen_highkaon_p    = 0.0;
          gen_lowkaon_p     = 0.0;

          gen_lowmuon_eta     = 0.0;
          gen_highmuon_eta    = 0.0;
          gen_highkaon_eta    = 0.0;
          gen_lowkaon_eta     = 0.0;

          gen_lowmuon_phi     = 0.0;
          gen_highmuon_phi    = 0.0;
          gen_highkaon_phi    = 0.0;
          gen_lowkaon_phi     = 0.0;

          reco::GenParticleRef genhighmuon  = highmuon->genParticleRef();
          reco::GenParticleRef genlowmuon   = lowmuon->genParticleRef();

          const reco::GenParticle *genhighkaon,*genlowkaon;
          reco::GenParticleRef phiMomHigh, phiMomLow, jpsiMomHigh, jpsiMomLow;
          reco::GenParticleRef jpsiMom, phiMom;

          float hasHighGen = dimuonditrk_cand.userFloat("hasHighGen");
          float hasLowGen = dimuonditrk_cand.userFloat("hasLowGen");

          if(hasHighGen>0.0)
          genhighkaon = dynamic_cast <const reco::GenParticle *>(dimuonditrk_cand.daughter("highkaonGen"));
          if(hasLowGen>0.0)
          genlowkaon = dynamic_cast <const reco::GenParticle *>(dimuonditrk_cand.daughter("lowkaonGen"));

          if(hasHighGen>0.0)
          {
            gen_highkaon_p4.SetPtEtaPhiM(genhighkaon->pt(),genhighkaon->eta(),genhighkaon->phi(),genhighkaon->mass());
            if(genhighkaon->numberOfMothers()>0)
            phiMomLow  = genhighkaon->motherRef();

            gen_highkaon_pdg     = (float)genhighkaon->pdgId();

            if(phiMomLow.isNonnull())
            gen_highkaon_mompdg  = phiMomLow->pdgId();

            gen_highkaon_status  = (float)genhighkaon->status();
            gen_highkaon_pt      = (float)genhighkaon->pt();
            gen_highkaon_p       = (float)genhighkaon->p();
            gen_highkaon_eta     = (float)genhighkaon->eta();
            gen_highkaon_phi     = (float)genhighkaon->phi();

          }

          if(hasLowGen>0.0)
          {
            gen_lowkaon_p4.SetPtEtaPhiM(genhighkaon->pt(),genhighkaon->eta(),genhighkaon->phi(),genhighkaon->mass());
            if(genhighkaon->numberOfMothers()>0)
            phiMomLow  = genhighkaon->motherRef();

            gen_lowkaon_pdg     = (float)genhighkaon->pdgId();

            if(phiMomLow.isNonnull())
            gen_lowkaon_mompdg  = phiMomLow->pdgId();

            gen_lowkaon_status  = (float)genhighkaon->status();
            gen_lowkaon_pt      = (float)genhighkaon->pt();
            gen_lowkaon_p       = (float)genhighkaon->p();
            gen_lowkaon_eta     = (float)genhighkaon->eta();
            gen_lowkaon_phi     = (float)genhighkaon->phi();
          }

          if(genhighmuon.isNonnull())
          {
            gen_highmuon_p4.SetPtEtaPhiM(genhighkaon->pt(),genhighkaon->eta(),genhighkaon->phi(),genhighkaon->mass());
            if(genhighkaon->numberOfMothers()>0)
            jpsiMomHigh  = genhighkaon->motherRef();

            gen_highmuon_pdg     = (float)genhighkaon->pdgId();

            if(jpsiMomHigh.isNonnull())
            gen_highmuon_mompdg  = jpsiMomHigh->pdgId();

            gen_highmuon_status  = (float)genhighkaon->status();
            gen_highmuon_pt      = (float)genhighkaon->pt();
            gen_highmuon_p       = (float)genhighkaon->p();
            gen_highmuon_eta     = (float)genhighkaon->eta();
            gen_highmuon_phi     = (float)genhighkaon->phi();
          }

          if(genlowmuon.isNonnull())
          {
            gen_lowmuon_p4.SetPtEtaPhiM(genlowkaon->pt(),genlowkaon->eta(),genlowkaon->phi(),genlowkaon->mass());
            if(genlowkaon->numberOfMothers()>0)
            jpsiMomLow  = genlowkaon->motherRef();

            gen_lowmuon_pdg     = (float)genlowkaon->pdgId();

            if(jpsiMomLow.isNonnull())
            gen_lowmuon_mompdg  = jpsiMomLow->pdgId();

            gen_lowmuon_status  = (float)genlowkaon->status();
            gen_lowmuon_pt      = (float)genlowkaon->pt();
            gen_lowmuon_p       = (float)genlowkaon->p();
            gen_lowmuon_eta     = (float)genlowkaon->eta();
            gen_lowmuon_phi     = (float)genlowkaon->phi();
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
              gen_jpsi_p4.SetPtEtaPhiM(jpsiMomHigh->pt(),jpsiMomHigh->eta(),jpsiMomHigh->jpsi(),jpsiMomHigh->mass());
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

            if(jpsiMom==phiMom && jpsiMom.isNonnull() && phiMom.isNonnull())
            {
              gen_dimuonditrk_p4.SetPtEtaPhiM(jpsiMom->pt(),jpsiMom->eta(),jpsiMom->phi(),jpsiMom->mass());
              gen_dimuonditrk_pdg = (float) jpsiMom->pdgId();
              gen_dimuonditrk_prompt = (float) jpsiMom->isPromptDecayed();
              gen_dimuonditrk_p = (float) jpsiMom->p();
              gen_dimuonditrk_pt = (float) jpsiMom->pt();
              gen_dimuonditrk_eta = (float) jpsiMom->eta();
              gen_dimuonditrk_phi = (float) jpsiMom->phi();
            }

          }

        }

        dimuonditrk_tree->Fill();

        if (onlyBest_) break;
        else
        isBestCandidate = false;

        // dimuontt candidates are sorted by vProb
      }

    }
  }

}

// ------------ method called once each job just before starting event loop  ------------
void DiMuonDiTrakFiveRootuplerFit::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void DiMuonDiTrakFiveRootuplerFit::endJob() {}

// ------------ method called when starting to processes a run  ------------
void DiMuonDiTrakFiveRootuplerFit::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void DiMuonDiTrakFiveRootuplerFit::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void DiMuonDiTrakFiveRootuplerFit::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void DiMuonDiTrakFiveRootuplerFit::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DiMuonDiTrakFiveRootuplerFit::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonDiTrakFiveRootuplerFit);
