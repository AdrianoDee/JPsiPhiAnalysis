/*
Package:    DiMuonDiTrakRootupler
Class:      DiMuonDiTrakRootupler

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

class DiMuonDiTrakRootupler : public edm::EDAnalyzer {
public:
  explicit DiMuonDiTrakRootupler(const edm::ParameterSet&);
  ~DiMuonDiTrakRootupler() override;

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
  edm::EDGetTokenT<reco::VertexCollection> pVertices_;
  edm::EDGetTokenT<edm::TriggerResults> triggers_;
  bool isMC_,onlyBest_,OnlyGen_ ;
  std::vector<std::string>  hlts_;
  std::vector<std::string>  hltFilters_;
  std::string treeName_;

  UInt_t run, event, lumi, numPrimaryVertices, trigger;

  TLorentzVector dimuonditrk_p4;
  TLorentzVector dimuon_p4;
  TLorentzVector ditrak_p4;
  TLorentzVector lowMuon_p4;
  TLorentzVector highMuon_p4;
  TLorentzVector highKaon_p4;
  TLorentzVector lowKaon_p4;

  TLorentzVector dimuonditrk_rf_p4;
  TLorentzVector dimuonditrk_rf_const_p4;
  TLorentzVector dimuon_rf_p4, dimuon_not_rf_p4;
  TLorentzVector ditrak_rf_p4, ditrak_not_rf_p4;
  TLorentzVector lowMuon_rf_p4;
  TLorentzVector highMuon_rf_p4;
  TLorentzVector kaonp_rf_p4;
  TLorentzVector kaonn_rf_p4;

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
  Double_t highKaon_pt,lowKaon_pt,highMuon_pt,lowMuon_pt,dimuonditrk_nDof,dimuonditrk_m_rf,dimuonditrk_m_rf_c,dimuonditrk_m_rf_d_c;

  Bool_t lowMuon_isLoose, lowMuon_isSoft, lowMuon_isMedium, lowMuon_isHighPt, lowMuon_isTight;
  Bool_t highMuon_isLoose, highMuon_isSoft, highMuon_isMedium, highMuon_isHighPt, highMuon_isTight;

  Bool_t lowMuon_isTracker, lowMuon_isGlobal, highMuon_isTracker, highMuon_isGlobal;
  UInt_t lowMuon_type, highMuon_type;

  Bool_t lowMuon_rf_isLoose, lowMuon_rf_isSoft, lowMuon_rf_isMedium, lowMuon_rf_isHighPt;
  Bool_t highMuon_rf_isLoose, highMuon_rf_isSoft, highMuon_rf_isMedium, highMuon_rf_isHighPt;

  Bool_t lowMuon_rf_isTracker, lowMuon_rf_isGlobal, highMuon_rf_isTracker, highMuon_rf_isGlobal;
  UInt_t lowMuon_rf_type, highMuon_rf_type;

  UInt_t lowMuon_NPixelHits, lowMuon_NStripHits, lowMuon_NTrackhits, lowMuon_NBPixHits, lowMuon_NPixLayers, lowMuon_NTraLayers, lowMuon_NStrLayers, lowMuon_NBPixLayers;
  UInt_t highMuon_NPixelHits, highMuon_NStripHits, highMuon_NTrackhits, highMuon_NBPixHits, highMuon_NPixLayers, highMuon_NTraLayers, highMuon_NStrLayers, highMuon_NBPixLayers;

  // Double_t highKaon_m, lowKaon_m, highMuon_y, lowMuon_y;

  Double_t dimuonditrk_refPK_mass, dimuonditrk_refKP_mass, dimuonditrk_refPP_mass, dimuonditrk_refPK_vChi2;
  Double_t dimuonditrk_refKP_vChi2, dimuonditrk_refPP_vChi2, dimuonditrk_refPK_nDof, dimuonditrk_refKP_nDof;
  Double_t dimuonditrk_refPP_nDof, dimuonditrk_refPK_vProb, dimuonditrk_refKP_vProb, dimuonditrk_refPP_vProb;

  Double_t highKaon_eta, lowKaon_eta, highMuon_eta, lowMuon_eta, highKaon_phi, lowKaon_phi, highMuon_phi, lowMuon_phi;
  Double_t highKaon_dz, lowKaon_dz, highMuon_dz, lowMuon_dz, highKaon_dxy, lowKaon_dxy, highMuon_dxy, lowMuon_dxy;
  // Double_t highKaon_etaError, lowKaon_etaError, highMuon_etaError, lowMuon_etaError, highKaon_phiError, lowKaon_phiError, highMuon_phiError, lowMuon_phiError;
  Int_t highKaon_NPixelHits, highKaon_NStripHits, highKaon_NTrackhits, highKaon_NBPixHits, highKaon_NPixLayers;
  Int_t highKaon_NTraLayers, highKaon_NStrLayers, highKaon_NBPixLayers, lowKaon_NPixelHits, lowKaon_NStripHits;
  Int_t lowKaon_NTrackhits, lowKaon_NBPixHits, lowKaon_NPixLayers, lowKaon_NTraLayers, lowKaon_NStrLayers, lowKaon_NBPixLayers;

  UInt_t highKaonMatch, lowKaonMatch,lowMuonMatch, highMuonMatch;

  Int_t dimuonditrk_rf_bindx;

  Int_t noXCandidates;

  Bool_t isBestCandidate;

  char *hltword;
  //MC
  TLorentzVector gen_dimuonditrk_p4, gen_jpsi_p4, gen_phi_p4;
  TLorentzVector gen_lowMuon_p4, gen_highMuon_p4, gen_highKaon_p4, gen_lowKaon_p4;

  Double_t gen_dimuonditrk_pdg, gen_phi_pdg, gen_jpsi_pdg;
  Double_t gen_lowMuon_pdg, gen_highMuon_pdg, gen_highKaon_pdg, gen_lowKaon_pdg;
  Double_t gen_lowMuon_mompdg, gen_highMuon_mompdg, gen_highKaon_mompdg, gen_lowKaon_mompdg;
  Double_t gen_lowMuon_status, gen_highMuon_status, gen_highKaon_status, gen_lowKaon_status;
  Double_t gen_dimuonditrk_prompt, gen_phi_prompt, gen_jpsi_prompt;
  Double_t gen_dimuonditrk_pt, gen_dimuonditrk_p, gen_dimuonditrk_eta;
  Double_t gen_phi_pt, gen_phi_p, gen_phi_eta;
  Double_t gen_jpsi_pt, gen_jpsi_p, gen_jpsi_eta;
  Double_t gen_lowMuon_phi, gen_highMuon_phi, gen_highKaon_phi, gen_lowKaon_phi;
  Double_t gen_dimuonditrk_phi, gen_phi_phi, gen_jpsi_phi;
  Double_t gen_lowMuon_p, gen_highMuon_p, gen_highKaon_p, gen_lowKaon_p;
  Double_t gen_lowMuon_pt, gen_highMuon_pt, gen_highKaon_pt, gen_lowKaon_pt;
  Double_t gen_lowMuon_eta, gen_highMuon_eta, gen_highKaon_eta, gen_lowKaon_eta;

  TTree* dimuonditrk_tree, *dimuonditrk_tree_rf;
  edm::EDGetTokenT< std::vector <reco::GenParticle> > genCands_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;
};

UInt_t DiMuonDiTrakRootupler::isTriggerMatched(pat::CompositeCandidate *diMuon_cand) {
  const pat::Muon* highMuon = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("highMuon"));
  const pat::Muon* lowMuon = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("lowMuon"));
  UInt_t matched = 0;  // if no list is given, is not matched

  // if matched a given trigger, set the bit, in the same order as listed
  for (unsigned int iTr = 0; iTr<hltFilters_.size(); iTr++ ) {
    // std::cout << hltFilters_[iTr] << std::endl;
    const pat::TriggerObjectStandAloneCollection mu1HLTMatches = highMuon->triggerObjectMatchesByFilter(hltFilters_[iTr]);
    const pat::TriggerObjectStandAloneCollection mu2HLTMatches = lowMuon->triggerObjectMatchesByFilter(hltFilters_[iTr]);
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
DiMuonDiTrakRootupler::DiMuonDiTrakRootupler(const edm::ParameterSet& iConfig):
dimuonditrk_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("DiMuoDiTrak"))),
pVertices_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PrimaryVertices"))),
triggers_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
isMC_(iConfig.getParameter<bool>("isMC")),
onlyBest_(iConfig.getParameter<bool>("OnlyBest")),
OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
hlts_(iConfig.getParameter<std::vector<std::string>>("HLTs")),
hltFilters_(iConfig.getParameter<std::vector<std::string>>("Filters")),
treeName_(iConfig.getParameter<std::string>("TreeName"))
{
  edm::Service<TFileService> fs;
  dimuonditrk_tree = fs->make<TTree>(treeName_.data(),"Tree of DiMuon and DiTrak");


  dimuonditrk_tree->Branch("run",                &run,                "run/I");
  dimuonditrk_tree->Branch("event",              &event,              "event/I");
  dimuonditrk_tree->Branch("lumi",              &lumi,              "lumi/I");
  dimuonditrk_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");
  dimuonditrk_tree->Branch("trigger",            &trigger,            "trigger/I");

  // dimuonditrk_tree->Branch("hltword",            hltword,            "hltword/C");

  if(!OnlyGen_)
  {
    dimuonditrk_tree->Branch("noXCandidates",      &noXCandidates,      "noXCandidates/I");

    //p4s
    dimuonditrk_tree->Branch("dimuonditrk_p4",   "TLorentzVector", &dimuonditrk_p4);
    dimuonditrk_tree->Branch("ditrak_p4",     "TLorentzVector", &ditrak_p4);
    dimuonditrk_tree->Branch("dimuon_p4",     "TLorentzVector", &dimuon_p4);
    dimuonditrk_tree->Branch("lowMuon_p4",   "TLorentzVector", &lowMuon_p4);
    dimuonditrk_tree->Branch("highMuon_p4",   "TLorentzVector", &highMuon_p4);
    dimuonditrk_tree->Branch("highKaon_p4",   "TLorentzVector", &highKaon_p4);
    dimuonditrk_tree->Branch("lowKaon_p4",   "TLorentzVector", &lowKaon_p4);


    //refitted p4s
    dimuonditrk_tree->Branch("dimuonditrk_rf_p4",   "TLorentzVector", &dimuonditrk_rf_p4);
    dimuonditrk_tree->Branch("ditrak_rf_p4",     "TLorentzVector", &ditrak_rf_p4);
    dimuonditrk_tree->Branch("dimuon_rf_p4",     "TLorentzVector", &dimuon_rf_p4);
    dimuonditrk_tree->Branch("lowMuon_rf_p4",   "TLorentzVector", &lowMuon_rf_p4);
    dimuonditrk_tree->Branch("highMuon_rf_p4",   "TLorentzVector", &highMuon_rf_p4);
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

    dimuonditrk_tree->Branch("highKaon_pt",          &highKaon_pt,          "highKaon_pt/D");
    dimuonditrk_tree->Branch("lowKaon_pt",       &lowKaon_pt,       "lowKaon_pt/D");
    dimuonditrk_tree->Branch("highMuon_pt",    &highMuon_pt,    "highMuon_pt/D");
    dimuonditrk_tree->Branch("lowMuon_pt",     &lowMuon_pt,     "lowMuon_pt/D");

    dimuonditrk_tree->Branch("highKaon_eta",        &highKaon_eta,        "highKaon_eta/D");
    dimuonditrk_tree->Branch("lowKaon_eta",        &lowKaon_eta,        "lowKaon_eta/D");
    dimuonditrk_tree->Branch("highMuon_eta",        &highMuon_eta,        "highMuon_eta/D");
    dimuonditrk_tree->Branch("lowMuon_eta",        &lowMuon_eta,        "lowMuon_eta/D");

    dimuonditrk_tree->Branch("highKaon_phi",        &highKaon_phi,        "highKaon_phi/D");
    dimuonditrk_tree->Branch("lowKaon_phi",        &lowKaon_phi,        "lowKaon_phi/D");
    dimuonditrk_tree->Branch("highMuon_phi",        &highMuon_phi,        "highMuon_phi/D");
    dimuonditrk_tree->Branch("lowMuon_phi",        &lowMuon_phi,        "lowMuon_phi/D");

    dimuonditrk_tree->Branch("highKaon_dxy",        &highKaon_dxy,        "highKaon_dxy/D");
    dimuonditrk_tree->Branch("lowKaon_dxy",        &lowKaon_dxy,        "lowKaon_dxy/D");
    dimuonditrk_tree->Branch("highMuon_dxy",        &highMuon_dxy,        "highMuon_dxy/D");
    dimuonditrk_tree->Branch("lowMuon_dxy",        &lowMuon_dxy,        "lowMuon_dxy/D");

    dimuonditrk_tree->Branch("highKaon_dz",        &highKaon_dz,        "highKaon_dz/D");
    dimuonditrk_tree->Branch("lowKaon_dz",        &lowKaon_dz,        "lowKaon_dz/D");
    dimuonditrk_tree->Branch("highMuon_dz",        &highMuon_dz,        "highMuon_dz/D");
    dimuonditrk_tree->Branch("lowMuon_dz",        &lowMuon_dz,        "lowMuon_dz/D");


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
    // dimuonditrk_tree->Branch("highKaon_trig_pt",          &highKaon_trig_pt,          "highKaon_trig_pt/D");
    // dimuonditrk_tree->Branch("lowKaon_trig_pt",       &lowKaon_trig_pt,       "lowKaon_trig_pt/D");
    // dimuonditrk_tree->Branch("highMuon_trig_pt",    &highMuon_trig_pt,    "highMuon_trig_pt/D");
    // dimuonditrk_tree->Branch("lowMuon_trig_pt",     &lowMuon_trig_pt,     "lowMuon_trig_pt/D");

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

    dimuonditrk_tree->Branch("highKaonMatch",     &highKaonMatch,       "highKaonMatch/I");
    dimuonditrk_tree->Branch("lowKaonMatch",     &lowKaonMatch,       "lowKaonMatch/I");
    dimuonditrk_tree->Branch("lowMuonMatch",     &lowMuonMatch,       "lowMuonMatch/I");
    dimuonditrk_tree->Branch("highMuonMatch",     &highMuonMatch,       "highMuonMatch/I");

    //Muon flags
    dimuonditrk_tree->Branch("lowMuon_isTight",        &lowMuon_isTight,        "lowMuon_isTight/O");
    dimuonditrk_tree->Branch("lowMuon_isLoose",        &lowMuon_isLoose,        "lowMuon_isLoose/O");
    dimuonditrk_tree->Branch("lowMuon_isSoft",        &lowMuon_isSoft,        "lowMuon_isSoft/O");
    dimuonditrk_tree->Branch("lowMuon_isMedium",        &lowMuon_isMedium,        "lowMuon_isMedium/O");
    dimuonditrk_tree->Branch("lowMuon_isHighPt",        &lowMuon_isHighPt,        "lowMuon_isHighPt/O");

    dimuonditrk_tree->Branch("lowMuon_isTracker",        &lowMuon_isTracker,        "lowMuon_isTracker/O");
    dimuonditrk_tree->Branch("lowMuon_isGlobal",        &lowMuon_isGlobal,        "lowMuon_isGlobal/O");

    dimuonditrk_tree->Branch("lowMuon_NPixelHits",        &lowMuon_NPixelHits,        "lowMuon_NPixelHits/I");
    dimuonditrk_tree->Branch("lowMuon_NStripHits",        &lowMuon_NStripHits,        "lowMuon_NStripHits/I");
    dimuonditrk_tree->Branch("lowMuon_NTrackhits",        &lowMuon_NTrackhits,        "lowMuon_NTrackhits/I");
    dimuonditrk_tree->Branch("lowMuon_NBPixHits",        &lowMuon_NBPixHits,        "lowMuon_NBPixHits/I");

    dimuonditrk_tree->Branch("lowMuon_NPixLayers",        &lowMuon_NPixLayers,        "lowMuon_NPixLayers/I");
    dimuonditrk_tree->Branch("lowMuon_NTraLayers",        &lowMuon_NTraLayers,        "lowMuon_NTraLayers/I");
    dimuonditrk_tree->Branch("lowMuon_NStrLayers",        &lowMuon_NStrLayers,        "lowMuon_NStrLayers/I");
    dimuonditrk_tree->Branch("lowMuon_NBPixLayers",        &lowMuon_NBPixLayers,        "lowMuon_NBPixLayers/I");

    dimuonditrk_tree->Branch("highMuon_isTight",        &highMuon_isTight,        "highMuon_isTight/O");
    dimuonditrk_tree->Branch("highMuon_isLoose",        &highMuon_isLoose,        "highMuon_isLoose/O");
    dimuonditrk_tree->Branch("highMuon_isSoft",        &highMuon_isSoft,        "highMuon_isSoft/O");
    dimuonditrk_tree->Branch("highMuon_isMedium",        &highMuon_isMedium,        "highMuon_isMedium/O");
    dimuonditrk_tree->Branch("highMuon_isHighPt",        &highMuon_isHighPt,        "highMuon_isHighPt/O");

    dimuonditrk_tree->Branch("highMuon_isTracker",        &highMuon_isTracker,        "highMuon_isTracker/O");
    dimuonditrk_tree->Branch("highMuon_isGlobal",        &highMuon_isGlobal,        "highMuon_isGlobal/O");

    dimuonditrk_tree->Branch("highMuon_NPixelHits",        &highMuon_NPixelHits,        "highMuon_NPixelHits/I");
    dimuonditrk_tree->Branch("highMuon_NStripHits",        &highMuon_NStripHits,        "highMuon_NStripHits/I");
    dimuonditrk_tree->Branch("highMuon_NTrackhits",        &highMuon_NTrackhits,        "highMuon_NTrackhits/I");
    dimuonditrk_tree->Branch("highMuon_NBPixHits",        &highMuon_NBPixHits,        "highMuon_NBPixHits/I");

    dimuonditrk_tree->Branch("highMuon_NPixLayers",        &highMuon_NPixLayers,        "highMuon_NPixLayers/I");
    dimuonditrk_tree->Branch("highMuon_NTraLayers",        &highMuon_NTraLayers,        "highMuon_NTraLayers/I");
    dimuonditrk_tree->Branch("highMuon_NStrLayers",        &highMuon_NStrLayers,        "highMuon_NStrLayers/I");
    dimuonditrk_tree->Branch("highMuon_NBPixLayers",        &highMuon_NBPixLayers,        "highMuon_NBPixLayers/I");

    dimuonditrk_tree->Branch("lowMuon_type",     &lowMuon_type,       "lowMuon_type/i");
    dimuonditrk_tree->Branch("highMuon_type",     &highMuon_type,       "highMuon_type/i");

    //Tracks Flags

    dimuonditrk_tree->Branch("highKaon_NPixelHits",        &highKaon_NPixelHits,        "highKaon_NPixelHits/I");
    dimuonditrk_tree->Branch("highKaon_NStripHits",        &highKaon_NStripHits,        "highKaon_NStripHits/I");
    dimuonditrk_tree->Branch("highKaon_NTrackhits",        &highKaon_NTrackhits,        "highKaon_NTrackhits/I");
    dimuonditrk_tree->Branch("highKaon_NBPixHits",        &highKaon_NBPixHits,        "highKaon_NBPixHits/I");

    dimuonditrk_tree->Branch("highKaon_NPixLayers",        &highKaon_NPixLayers,        "highKaon_NPixLayers/I");
    dimuonditrk_tree->Branch("highKaon_NTraLayers",        &highKaon_NTraLayers,        "highKaon_NTraLayers/I");
    dimuonditrk_tree->Branch("highKaon_NStrLayers",        &highKaon_NStrLayers,        "highKaon_NStrLayers/I");
    dimuonditrk_tree->Branch("highKaon_NBPixLayers",        &highKaon_NBPixLayers,        "highKaon_NBPixLayers/I");

    dimuonditrk_tree->Branch("lowKaon_NPixelHits",        &lowKaon_NPixelHits,        "lowKaon_NPixelHits/I");
    dimuonditrk_tree->Branch("lowKaon_NStripHits",        &lowKaon_NStripHits,        "lowKaon_NStripHits/I");
    dimuonditrk_tree->Branch("lowKaon_NTrackhits",        &lowKaon_NTrackhits,        "lowKaon_NTrackhits/I");
    dimuonditrk_tree->Branch("lowKaon_NBPixHits",        &lowKaon_NBPixHits,        "lowKaon_NBPixHits/I");

    dimuonditrk_tree->Branch("lowKaon_NPixLayers",        &lowKaon_NPixLayers,        "lowKaon_NPixLayers/I");
    dimuonditrk_tree->Branch("lowKaon_NTraLayers",        &lowKaon_NTraLayers,        "lowKaon_NTraLayers/I");
    dimuonditrk_tree->Branch("lowKaon_NStrLayers",        &lowKaon_NStrLayers,        "lowKaon_NStrLayers/I");
    dimuonditrk_tree->Branch("lowKaon_NBPixLayers",        &lowKaon_NBPixLayers,        "lowKaon_NBPixLayers/I");

  }

  if (isMC_ || OnlyGen_) {

    dimuonditrk_tree->Branch("gen_dimuonditrk_p4", "TLorentzVector",  &gen_dimuonditrk_p4);
    dimuonditrk_tree->Branch("gen_jpsi_p4", "TLorentzVector",  &gen_jpsi_p4);
    dimuonditrk_tree->Branch("gen_phi_p4", "TLorentzVector",  &gen_phi_p4);

    dimuonditrk_tree->Branch("gen_highKaon_p4", "TLorentzVector",  &gen_highKaon_p4);
    dimuonditrk_tree->Branch("gen_lowMuon_p4",  "TLorentzVector",  &gen_lowMuon_p4);
    dimuonditrk_tree->Branch("gen_highMuon_p4",  "TLorentzVector",  &gen_highMuon_p4);
    dimuonditrk_tree->Branch("gen_lowKaon_p4",  "TLorentzVector",  &gen_lowKaon_p4);

    dimuonditrk_tree->Branch("gen_dimuonditrk_pdg",&gen_dimuonditrk_pdg,"gen_dimuonditrk_pdg/D");
    dimuonditrk_tree->Branch("gen_phi_pdg",&gen_phi_pdg,"gen_phi_pdg/D");
    dimuonditrk_tree->Branch("gen_jpsi_pdg",&gen_jpsi_pdg,"gen_jpsi_pdg/D");

    dimuonditrk_tree->Branch("gen_lowMuon_pdg",&gen_lowMuon_pdg,"gen_lowMuon_pdg/D");
    dimuonditrk_tree->Branch("gen_highMuon_pdg",&gen_highMuon_pdg,"gen_highMuon_pdg/D");
    dimuonditrk_tree->Branch("gen_highKaon_pdg",&gen_highKaon_pdg,"gen_highKaon_pdg/D");
    dimuonditrk_tree->Branch("gen_lowKaon_pdg",&gen_lowKaon_pdg,"gen_lowKaon_pdg/D");

    dimuonditrk_tree->Branch("gen_lowMuon_mompdg",&gen_lowMuon_mompdg,"gen_lowMuon_mompdg/D");
    dimuonditrk_tree->Branch("gen_highMuon_mompdg",&gen_highMuon_mompdg,"gen_highMuon_mompdg/D");
    dimuonditrk_tree->Branch("gen_highKaon_mompdg",&gen_highKaon_mompdg,"gen_highKaon_mompdg/D");
    dimuonditrk_tree->Branch("gen_lowKaon_mompdg",&gen_lowKaon_mompdg,"gen_lowKaon_mompdg/D");

    dimuonditrk_tree->Branch("gen_lowMuon_status",&gen_lowMuon_status,"gen_lowMuon_status/D");
    dimuonditrk_tree->Branch("gen_highMuon_status",&gen_highMuon_status,"gen_highMuon_status/D");
    dimuonditrk_tree->Branch("gen_highKaon_status",&gen_highKaon_status,"gen_highKaon_status/D");
    dimuonditrk_tree->Branch("gen_lowKaon_status",&gen_lowKaon_status,"gen_lowKaon_status/D");

    dimuonditrk_tree->Branch("gen_lowMuon_p",&gen_lowMuon_p,"gen_lowMuon_p/D");
    dimuonditrk_tree->Branch("gen_highMuon_p",&gen_highMuon_p,"gen_highMuon_p/D");
    dimuonditrk_tree->Branch("gen_highKaon_p",&gen_highKaon_p,"gen_highKaon_p/D");
    dimuonditrk_tree->Branch("gen_lowKaon_p",&gen_lowKaon_p,"gen_lowKaon_p/D");

    dimuonditrk_tree->Branch("gen_lowMuon_pt",&gen_lowMuon_pt,"gen_lowMuon_pt/D");
    dimuonditrk_tree->Branch("gen_highMuon_pt",&gen_highMuon_pt,"gen_highMuon_pt/D");
    dimuonditrk_tree->Branch("gen_highKaon_pt",&gen_highKaon_pt,"gen_highKaon_pt/D");
    dimuonditrk_tree->Branch("gen_lowKaon_pt",&gen_lowKaon_pt,"gen_lowKaon_pt/D");

    dimuonditrk_tree->Branch("gen_lowMuon_eta",&gen_lowMuon_eta,"gen_lowMuon_eta/D");
    dimuonditrk_tree->Branch("gen_highMuon_eta",&gen_highMuon_eta,"gen_highMuon_eta/D");
    dimuonditrk_tree->Branch("gen_highKaon_eta",&gen_highKaon_eta,"gen_highKaon_eta/D");
    dimuonditrk_tree->Branch("gen_lowKaon_eta",&gen_lowKaon_eta,"gen_lowKaon_eta/D");

    dimuonditrk_tree->Branch("gen_lowMuon_phi",&gen_lowMuon_phi,"gen_lowMuon_phi/D");
    dimuonditrk_tree->Branch("gen_highMuon_phi",&gen_highMuon_phi,"gen_highMuon_phi/D");
    dimuonditrk_tree->Branch("gen_highKaon_phi",&gen_highKaon_phi,"gen_highKaon_phi/D");
    dimuonditrk_tree->Branch("gen_lowKaon_phi",&gen_lowKaon_phi,"gen_lowKaon_phi/D");

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

DiMuonDiTrakRootupler::~DiMuonDiTrakRootupler() {}

//
// member functions
//

bool DiMuonDiTrakRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
  if (ancestor == particle ) return true;
  for (size_t i=0; i< particle->numberOfMothers(); i++) {
    if (isAncestor(ancestor, particle->mother(i))) return true;
  }
  return false;
}


// ------------ method called for each event  ------------
void DiMuonDiTrakRootupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //  using namespace edm;
  using namespace std;

  edm::Handle<std::vector<pat::CompositeCandidate>> dimuonditrk_cand_handle;
  iEvent.getByToken(dimuonditrk_, dimuonditrk_cand_handle);

  edm::Handle<std::vector<reco::Vertex >> primaryVertices_handle;
  iEvent.getByToken(pVertices_, primaryVertices_handle);

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( triggers_ , triggerResults_handle);

  numPrimaryVertices = primaryVertices_handle->size();
  run = iEvent.id().run();
  event = iEvent.id().event();
  lumi = iEvent.id().luminosityBlock();

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

    // std::string hltstring = "";
    //
    // for (unsigned int i = 0; i < NTRIGGERS; i++)
    // {
    //   hltstring = hltstring + hlts_[i];
    // }
    // strcpy(hltword, hltstring.c_str());

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
  //   gen_highMuon_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
  //   foundit++;
  // }
  // if ( motherInPrunedCollection != nullptr && (d->pdgId() == -13 ) && isAncestor(aditrkdimu , motherInPrunedCollection) ) {
  //   gen_lowMuon_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
  //   foundit++;
  // }

  gen_dimuonditrk_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_phi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_lowMuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_highMuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_highKaon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_lowKaon_p4.SetPtEtaPhiM(0.,0.,0.,0.);

  // get rf information. Notice we are just keeping combinations with succesfull vertex fit
  if(!OnlyGen_)
  {
    if (!dimuonditrk_cand_handle.isValid()) std::cout<< "No dimuontt information " << run << "," << event <<std::endl;

      dimuonditrk_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      ditrak_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      lowMuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      highMuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      highKaon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      lowKaon_p4.SetPtEtaPhiM(0.,0.,0.,0.);

      dimuonditrk_rf_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      ditrak_rf_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      dimuon_rf_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      lowMuon_rf_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      highMuon_rf_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      // highKaon_rf_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      // lowKaon_rf_p4.SetPtEtaPhiM(0.,0.,0.,0.);

      pat::CompositeCandidate *dimuonditrk_rf_cand, dimuonditrk_cand, *dimuon_cand;
      pat::CompositeCandidate *ditrak_cand, *dimuon_cand_rf, *ditrak_cand_rf;

      noXCandidates = (Int_t)(dimuonditrk_cand_handle->size());

      for (unsigned int i=0; i< dimuonditrk_cand_handle->size(); i++){

        dimuonditrk_cand   = dimuonditrk_cand_handle->at(i);

        const reco::Vertex thePrimaryV = *(dimuonditrk_cand.userData<reco::Vertex>("bestPV"));

        pv_x = thePrimaryV.position().x();
        pv_y = thePrimaryV.position().y();
        pv_z = thePrimaryV.position().z();

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

        highKaonMatch = dimuonditrk_cand.userInt("highKaonMatch");
        lowKaonMatch = dimuonditrk_cand.userInt("lowKaonMatch");

        //unref corresponding

        dimuon_cand = dynamic_cast <pat::CompositeCandidate *>(dimuonditrk_cand.daughter("dimuon"));
        ditrak_cand = dynamic_cast <pat::CompositeCandidate *>(dimuonditrk_cand.daughter("ditrak"));

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

        lowMuon_isTight    = lowMuon->isTightMuon(thePrimaryV);
        lowMuon_isLoose    = lowMuon->isLooseMuon();
        lowMuon_isSoft     = lowMuon->isSoftMuon(thePrimaryV);
        lowMuon_isMedium   = lowMuon->isMediumMuon();
        lowMuon_isHighPt   = lowMuon->isHighPtMuon(thePrimaryV);
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

        highMuon_isTight    = highMuon->isTightMuon(thePrimaryV);
        highMuon_isLoose    = highMuon->isLooseMuon();
        highMuon_isSoft     = highMuon->isSoftMuon(thePrimaryV);
        highMuon_isMedium   = highMuon->isMediumMuon();
        highMuon_isHighPt   = highMuon->isHighPtMuon(thePrimaryV);
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

        highMuon_eta  = highMuon->innerTrack()->eta();
        highMuon_phi  = highMuon->innerTrack()->phi();
        highMuon_dz  = highMuon->innerTrack()->dz();
        highMuon_dxy  = highMuon->innerTrack()->dxy();
        lowMuon_eta  = lowMuon->innerTrack()->eta();
        lowMuon_phi  = lowMuon->innerTrack()->phi();
        lowMuon_dz  = lowMuon->innerTrack()->dz();
        lowMuon_dxy  = lowMuon->innerTrack()->dxy();

        lowMuon_type       = lowMuon->type();
        highMuon_type       = highMuon->type();

        lowMuon_p4.SetPtEtaPhiM(vhighMuon.pt(), vhighMuon.eta(), vhighMuon.phi(), vhighMuon.mass());
        highMuon_p4.SetPtEtaPhiM(vlowMuon.pt(), vlowMuon.eta(), vlowMuon.phi(), vlowMuon.mass());

        reco::Candidate::LorentzVector kP = ditrak_cand->daughter("highTrak")->p4();
        reco::Candidate::LorentzVector kM = ditrak_cand->daughter("lowTrak")->p4();

        pat::PackedCandidate *highKaon  = dynamic_cast <pat::PackedCandidate *>(ditrak_cand->daughter("highTrak"));
        pat::PackedCandidate *lowKaon  = dynamic_cast <pat::PackedCandidate *>(ditrak_cand->daughter("lowTrak"));

        highKaon_p4.SetPtEtaPhiM(kP.pt(), kP.eta(), kP.phi(), kP.mass());
        lowKaon_p4.SetPtEtaPhiM(kM.pt(), kM.eta(), kM.phi(), kM.mass());

        highKaon_pt     = std::max(kP.pt(),kM.pt());
        lowKaon_pt      = -std::max(-kP.pt(),-kM.pt());
        highMuon_pt     = std::max(vlowMuon.pt(),vhighMuon.pt());
        lowMuon_pt      = -std::max(-vlowMuon.pt(),-vhighMuon.pt());

        lowKaon_NPixelHits = lowKaon->bestTrack()->hitPattern().numberOfValidPixelHits();
        lowKaon_NStripHits = lowKaon->bestTrack()->hitPattern().numberOfValidStripHits();
        lowKaon_NTrackhits = lowKaon->bestTrack()->hitPattern().numberOfValidTrackerHits();
        lowKaon_NBPixHits  = lowKaon->bestTrack()->hitPattern().numberOfValidStripHits();
        lowKaon_NPixLayers = lowKaon->bestTrack()->hitPattern().pixelLayersWithMeasurement();
        lowKaon_NTraLayers = lowKaon->bestTrack()->hitPattern().trackerLayersWithMeasurement();
        lowKaon_NStrLayers = lowKaon->bestTrack()->hitPattern().stripLayersWithMeasurement();
        lowKaon_NBPixLayers = lowKaon->bestTrack()->hitPattern().pixelBarrelLayersWithMeasurement();
        highKaon_NPixelHits = highKaon->bestTrack()->hitPattern().numberOfValidPixelHits();
        highKaon_NStripHits = highKaon->bestTrack()->hitPattern().numberOfValidStripHits();
        highKaon_NTrackhits = highKaon->bestTrack()->hitPattern().numberOfValidTrackerHits();
        highKaon_NBPixHits  = highKaon->bestTrack()->hitPattern().numberOfValidStripHits();
        highKaon_NPixLayers = highKaon->bestTrack()->hitPattern().pixelLayersWithMeasurement();
        highKaon_NTraLayers = highKaon->bestTrack()->hitPattern().trackerLayersWithMeasurement();
        highKaon_NStrLayers = highKaon->bestTrack()->hitPattern().stripLayersWithMeasurement();
        highKaon_NBPixLayers = highKaon->bestTrack()->hitPattern().pixelBarrelLayersWithMeasurement();

        highKaon_eta  = highKaon->bestTrack()->eta();
        highKaon_phi  = highKaon->bestTrack()->phi();
        highKaon_dz   = highKaon->bestTrack()->dz();
        highKaon_dxy  = highKaon->bestTrack()->dxy();
        lowKaon_eta   = lowKaon->bestTrack()->eta();
        lowKaon_phi   = lowKaon->bestTrack()->phi();
        lowKaon_dz    = lowKaon->bestTrack()->dz();
        lowKaon_dxy   = lowKaon->bestTrack()->dxy();

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

        lowMuonMatch   = dimuon_cand->userInt("highMuonTMatch");
        highMuonMatch    = dimuon_cand->userInt("lowMuonTMatch");
        dimuon_triggerMatch = -std::max(-lowMuonMatch,highMuonMatch);

        dimuonditrk_rf_p4.SetPtEtaPhiM(0.0,0.0,0.0,3.9);
        dimuonditrk_rf_const_p4.SetPtEtaPhiM(0.0,0.0,0.0,3.9);
        dimuon_rf_p4.SetPtEtaPhiM(0.0,0.0,0.0,2.4);
        ditrak_rf_p4.SetPtEtaPhiM(0.0,0.0,0.0,0.5);

        lowMuon_rf_p4.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
        highMuon_rf_p4.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
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

          vhighMuon = dimuon_cand_rf->daughter("highMuon")->p4();
          vlowMuon = dimuon_cand_rf->daughter("lowMuon")->p4();

          lowMuon_rf_p4.SetPtEtaPhiM(vhighMuon.pt(), vhighMuon.eta(), vhighMuon.phi(), vhighMuon.mass());
          highMuon_rf_p4.SetPtEtaPhiM(vlowMuon.pt(), vlowMuon.eta(), vlowMuon.phi(), vlowMuon.mass());

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

        if(isMC_ || OnlyGen_)
        {

          gen_dimuonditrk_p4.SetPtEtaPhiM(-1.0,0.0,0.0,3.9);
          gen_jpsi_p4.SetPtEtaPhiM(-1.0,0.0,0.0,2.5);
          gen_phi_p4.SetPtEtaPhiM(-1.0,0.0,0.0,0.9);
          gen_highKaon_p4.SetPtEtaPhiM(-1.0,0.0,0.0,0.2);
          gen_lowMuon_p4.SetPtEtaPhiM(-1.0,0.0,0.0,0.0);
          gen_highMuon_p4.SetPtEtaPhiM(-1.0,0.0,0.0,0.0);
          gen_lowKaon_p4.SetPtEtaPhiM(-1.0,0.0,0.0,0.2);

          gen_dimuonditrk_pdg = 0.0;
          gen_phi_pdg         = 0.0;
          gen_jpsi_pdg        = 0.0;

          gen_lowMuon_pdg     = 0.0;
          gen_highMuon_pdg    = 0.0;
          gen_highKaon_pdg    = 0.0;
          gen_lowKaon_pdg     = 0.0;

          gen_lowMuon_mompdg     = 0.0;
          gen_highMuon_mompdg    = 0.0;
          gen_highKaon_mompdg    = 0.0;
          gen_lowKaon_mompdg     = 0.0;

          gen_lowMuon_status     = 0.0;
          gen_highMuon_status    = 0.0;
          gen_highKaon_status    = 0.0;
          gen_lowKaon_status     = 0.0;

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

          gen_lowMuon_pt     = 0.0;
          gen_highMuon_pt    = 0.0;
          gen_highKaon_pt    = 0.0;
          gen_lowKaon_pt     = 0.0;

          gen_lowMuon_p     = 0.0;
          gen_highMuon_p    = 0.0;
          gen_highKaon_p    = 0.0;
          gen_lowKaon_p     = 0.0;

          gen_lowMuon_eta     = 0.0;
          gen_highMuon_eta    = 0.0;
          gen_highKaon_eta    = 0.0;
          gen_lowKaon_eta     = 0.0;

          gen_lowMuon_phi     = 0.0;
          gen_highMuon_phi    = 0.0;
          gen_highKaon_phi    = 0.0;
          gen_lowKaon_phi     = 0.0;

          reco::GenParticleRef genhighMuon  = highMuon->genParticleRef();
          reco::GenParticleRef genlowMuon   = lowMuon->genParticleRef();

          const reco::GenParticle *genhighKaon,*genlowKaon;
          reco::GenParticleRef phiMomHigh, phiMomLow, jpsiMomHigh, jpsiMomLow;
          reco::GenParticleRef jpsiMom, phiMom;

          Double_t hasHighGen = dimuonditrk_cand.userFloat("hasHighGen");
          Double_t hasLowGen = dimuonditrk_cand.userFloat("hasLowGen");

          if(hasHighGen>0.0)
            genhighKaon = dynamic_cast <const reco::GenParticle *>(dimuonditrk_cand.daughter("highKaonGen"));
          if(hasLowGen>0.0)
            genlowKaon = dynamic_cast <const reco::GenParticle *>(dimuonditrk_cand.daughter("lowKaonGen"));

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

            if(jpsiMom==phiMom && jpsiMom.isNonnull() && phiMom.isNonnull())
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

        }

        dimuonditrk_tree->Fill();

        if (onlyBest_) break;
        else
        isBestCandidate = false;

        // dimuontt candidates are sorted by vProb
      }

    }
  }


// ------------ method called once each job just before starting event loop  ------------
void DiMuonDiTrakRootupler::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void DiMuonDiTrakRootupler::endJob() {}

// ------------ method called when starting to processes a run  ------------
void DiMuonDiTrakRootupler::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void DiMuonDiTrakRootupler::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void DiMuonDiTrakRootupler::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void DiMuonDiTrakRootupler::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DiMuonDiTrakRootupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonDiTrakRootupler);
