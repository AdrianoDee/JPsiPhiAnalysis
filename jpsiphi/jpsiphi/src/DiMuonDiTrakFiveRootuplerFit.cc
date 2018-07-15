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
  edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuonditrk_cand_Label;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> fivetracks_cand_Label;
  edm::EDGetTokenT<reco::BeamSpot> thebeamspot_;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  bool isMC_,OnlyBest_,OnlyGen_ ;
  UInt_t motherpdgid_,phipdgid_,jpspdgid_;
  std::vector<std::string>  HLTs_;
  std::vector<std::string>  HLTFilters_;
  std::string treeName_;

  UInt_t run,event,numPrimaryVertices, trigger;

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

  TLorentzVector fivetraks_p4;
  TLorentzVector fivetrakspion_p4;
  TLorentzVector fifthkaon_p4;
  TLorentzVector dimuontrak_p4;
  TLorentzVector dimuontrakpion_p4;

  Int_t dimuonditrk_charge;

  UInt_t dimuon_triggerMatch, dimuon_triggerMatch_rf;

  Double_t dimuonditrk_vProb,  dimuonditrk_vChi2;
  Double_t dimuonditrk_rf_vProb, dimuonditrk_rf_vChi2, dimuonditrk_rf_nDof, dimuonditrk_rf_cosAlpha, dimuonditrk_rf_ctauPV, dimuonditrk_rf_ctauErrPV;
  Double_t dimuonditrk_rf_c_vProb, dimuonditrk_rf_c_vChi2, dimuonditrk_rf_c_nDof, dimuonditrk_rf_c_cosAlpha, dimuonditrk_rf_c_ctauPV, dimuonditrk_rf_c_ctauErrPV;

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

  Double_t gen_dimuonditrk_m,dimuonditrk_m,dimuonditrk_pt,dimuon_m,dimuon_pt,ditrak_m,ditrak_pt;
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

  Double_t fivetraks_kaon_m, fivetraks_kaon_vProb, fivetraks_kaon_vChi2, fivetraks_kaon_nDof, fivetraks_kaon_charge;
  Double_t fivetraks_kaon_cosAlpha, fivetraks_kaon_ctauPV, fivetraks_kaon_ctauErrPV, fifthtrak_kaon_m, fifthtrak_kaon_eta;
  Double_t fifthtrak_kaon_pt, fifthtrak_kaon_phi;

  Double_t highKaon_y, lowKaon_y, highMuon_y, lowMuon_y;

  Double_t fivetraks_pion_m, fivetraks_pion_vProb, fivetraks_pion_vChi2, fivetraks_pion_nDof, fivetraks_pion_charge;
  Double_t fivetraks_pion_cosAlpha, fivetraks_pion_ctauPV, fivetraks_pion_ctauErrPV, fifthtrak_pion_m, fifthtrak_pion_eta;
  Double_t fifthtrak_pion_pt, fifthtrak_pion_phi;

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

  UInt_t dimuonditrk_pdgid,  dimuonditrk_jpsipdg;
  Double_t dimuonditrk_isprompt,dimuonditrk_jpsippdl;

  Int_t          gen_dimuonditrk_pdgId;
  TLorentzVector gen_dimuonditrk_p4;
  TLorentzVector gen_b_p4;
  TLorentzVector gen_dimuon_p4;
  TLorentzVector gen_ditrak_p4;
  TLorentzVector gen_lowMuon_p4;
  TLorentzVector gen_highMuon_p4;
  TLorentzVector gen_highKaon_p4;
  TLorentzVector gen_lowKaon_p4;

  TLorentzVector gen_b4_p4;
  TLorentzVector gen_d1_p4;
  TLorentzVector gen_d2_p4;
  TLorentzVector gen_gd1_p4;
  TLorentzVector gen_gd2_p4;
  TLorentzVector gen_gd3_p4;
  TLorentzVector gen_gd4_p4;
  TLorentzVector gen_gd5_p4;
  TLorentzVector gen_gd6_p4;

  TTree* dimuonditrk_tree, *dimuonditrk_tree_rf;
  edm::EDGetTokenT< std::vector <reco::GenParticle> > genCands_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;
};

UInt_t DiMuonDiTrakFiveRootuplerFit::isTriggerMatched(pat::CompositeCandidate *diMuon_cand) {
  const pat::Muon* highMuon = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("highMuon"));
  const pat::Muon* lowMuon = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("lowMuon"));
  UInt_t matched = 0;  // if no list is given, is not matched

  // if matched a given trigger, set the bit, in the same order as listed
  for (unsigned int iTr = 0; iTr<HLTFilters_.size(); iTr++ ) {
    // std::cout << HLTFilters_[iTr] << std::endl;
    const pat::TriggerObjectStandAloneCollection mu1HLTMatches = highMuon->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
    const pat::TriggerObjectStandAloneCollection mu2HLTMatches = lowMuon->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
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
DiMuonDiTrakFiveRootuplerFit::DiMuonDiTrakFiveRootuplerFit(const edm::ParameterSet& iConfig):
        dimuonditrk_cand_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("dimuonditrk_cand"))),
        fivetracks_cand_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("fivetracks_cand"))),
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
          dimuonditrk_tree->Branch("lowMuon_p4",   "TLorentzVector", &lowMuon_p4);
          dimuonditrk_tree->Branch("highMuon_p4",   "TLorentzVector", &highMuon_p4);
          dimuonditrk_tree->Branch("highKaon_p4",   "TLorentzVector", &highKaon_p4);
          dimuonditrk_tree->Branch("lowKaon_p4",   "TLorentzVector", &lowKaon_p4);
          dimuonditrk_tree->Branch("fivetraks_p4",   "TLorentzVector", &fivetraks_p4);
          dimuonditrk_tree->Branch("fivetrakspion_p4",   "TLorentzVector", &fivetrakspion_p4);
          dimuonditrk_tree->Branch("fifthkaon_p4",   "TLorentzVector", &fifthkaon_p4);
          dimuonditrk_tree->Branch("dimuontrak_p4",   "TLorentzVector", &dimuontrak_p4);
          dimuonditrk_tree->Branch("dimuontrakpion_p4",   "TLorentzVector", &dimuontrakpion_p4);

          //refitted p4s
          dimuonditrk_tree->Branch("dimuonditrk_rf_p4",   "TLorentzVector", &dimuonditrk_rf_p4);
          dimuonditrk_tree->Branch("ditrak_rf_p4",     "TLorentzVector", &ditrak_rf_p4);
          dimuonditrk_tree->Branch("dimuon_rf_p4",     "TLorentzVector", &dimuon_rf_p4);
          dimuonditrk_tree->Branch("lowMuon_rf_p4",   "TLorentzVector", &lowMuon_rf_p4);
          dimuonditrk_tree->Branch("highMuon_rf_p4",   "TLorentzVector", &highMuon_rf_p4);
          dimuonditrk_tree->Branch("kaonp_rf_p4",   "TLorentzVector", &kaonp_rf_p4);
          dimuonditrk_tree->Branch("kaonn_rf_p4",   "TLorentzVector", &kaonn_rf_p4);

          //kin
          dimuonditrk_tree->Branch("gen_dimuonditrk_m",        &gen_dimuonditrk_m,        "gen_dimuonditrk_m/D");
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

          dimuonditrk_tree->Branch("highKaon_y",        &highKaon_y,        "highKaon_y/D");
          dimuonditrk_tree->Branch("lowKaon_y",        &lowKaon_y,        "lowKaon_y/D");
          dimuonditrk_tree->Branch("highMuon_y",        &highMuon_y,        "highMuon_y/D");
          dimuonditrk_tree->Branch("lowMuon_y",        &lowMuon_y,        "lowMuon_y/D");

          //FiveTracks - Kaon

          dimuonditrk_tree->Branch("fivetraks_kaon_m",     &fivetraks_kaon_m,     "fivetraks_kaon_m/D");
          dimuonditrk_tree->Branch("fivetraks_pion_m",     &fivetraks_pion_m,     "fivetraks_pion_m/D");
          dimuonditrk_tree->Branch("fivetraks_kaon_trim",     &fivetraks_kaon_trim,     "fivetraks_kaon_trim/D");
          dimuonditrk_tree->Branch("fivetraks_pion_trim",     &fivetraks_pion_trim,     "fivetraks_pion_trim/D");
          dimuonditrk_tree->Branch("fivetraks_kaon_m_rf",     &fivetraks_kaon_m_rf,     "fivetraks_kaon_m_rf/D");
          dimuonditrk_tree->Branch("fivetraks_pion_m_rf",     &fivetraks_pion_m_rf,     "fivetraks_pion_m_rf/D");
          dimuonditrk_tree->Branch("fivetraks_vProb",      &fivetraks_vProb,        "fivetraks_vProb/D");
          dimuonditrk_tree->Branch("fivetraks_vChi2",      &fivetraks_vChi2,        "fivetraks_vChi2/D");
          dimuonditrk_tree->Branch("fivetraks_nDof",       &fivetraks_nDof,         "fivetraks_nDof/D");
          dimuonditrk_tree->Branch("fivetraks_charge",     &fivetraks_charge,       "fivetraks_charge/I");
          dimuonditrk_tree->Branch("fivetraks_cosAlpha",      &fivetraks_cosAlpha,        "fivetraks_cosAlpha/D");
          dimuonditrk_tree->Branch("fivetraks_ctauPV",      &fivetraks_ctauPV,        "fivetraks_ctauPV/D");
          dimuonditrk_tree->Branch("fivetraks_ctauErrPV",      &fivetraks_ctauErrPV,        "fivetraks_ctauErrPV/D");
          dimuonditrk_tree->Branch("fivetraks_eta",       &fivetraks_eta,        "fivetraks_eta/D");
          dimuonditrk_tree->Branch("fivetraks_pt",       &fivetraks_pt,        "fivetraks_pt/D");
          dimuonditrk_tree->Branch("fivetraks_phi",       &fivetraks_phi,        "fivetraks_phi/D");
          dimuonditrk_tree->Branch("fivetraks_y",       &fivetraks_y,        "fivetraks_y/D");


          dimuonditrk_tree->Branch("fifthtrak_charge",       &fifthtrak_charge,        "fifthtrak_charge/D");
          dimuonditrk_tree->Branch("fifthtrak_eta",       &fifthtrak_eta,        "fifthtrak_eta/D");
          dimuonditrk_tree->Branch("fifthtrak_pt",       &fifthtrak_pt,        "fifthtrak_pt/D");
          dimuonditrk_tree->Branch("fifthtrak_phi",       &fifthtrak_phi,        "fifthtrak_phi/D");
          dimuonditrk_tree->Branch("fifthtrak_y",       &fifthtrak_y,        "fifthtrak_y/D");


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
          dimuonditrk_tree->Branch("dimuonditrk_countTksOfPV",      &dimuonditrk_countTksOfPV,        "dimuonditrk_countTksOfPV/D");
          dimuonditrk_tree->Branch("dimuonditrk_vertexWeight",      &dimuonditrk_vertexWeight,        "dimuonditrk_vertexWeight/D");
          dimuonditrk_tree->Branch("dimuonditrk_sumPTPV",      &dimuonditrk_sumPTPV,        "dimuonditrk_sumPTPV/D");
          dimuonditrk_tree->Branch("dimuonditrk_mu1FromPV",      &dimuonditrk_mu1FromPV,        "dimuonditrk_mu1FromPV/D");
          dimuonditrk_tree->Branch("dimuonditrk_mu2FromPV",      &dimuonditrk_mu2FromPV,        "dimuonditrk_mu2FromPV/D");
          dimuonditrk_tree->Branch("dimuonditrk_tPFromPV",      &dimuonditrk_tPFromPV,        "dimuonditrk_tPFromPV/D");
          dimuonditrk_tree->Branch("dimuonditrk_tMFromPV",      &dimuonditrk_tMFromPV,        "dimuonditrk_tMFromPV/D");
          dimuonditrk_tree->Branch("dimuonditrk_mu1W",      &dimuonditrk_mu1W,        "dimuonditrk_mu1W/D");
          dimuonditrk_tree->Branch("dimuonditrk_mu2W",      &dimuonditrk_mu2W,        "dimuonditrk_mu2W/D");
          dimuonditrk_tree->Branch("dimuonditrk_tPW",      &dimuonditrk_tPW,        "dimuonditrk_tPW/D");
          dimuonditrk_tree->Branch("dimuonditrk_tMW",      &dimuonditrk_tMW,        "dimuonditrk_tMW/D");

          dimuonditrk_tree->Branch("dimuonditrk_cosAlphaDZ",      &dimuonditrk_cosAlphaDZ,        "dimuonditrk_cosAlphaDZ/D");
          dimuonditrk_tree->Branch("dimuonditrk_ctauPVDZ",      &dimuonditrk_ctauPVDZ,        "dimuonditrk_ctauPVDZ/D");
          dimuonditrk_tree->Branch("dimuonditrk_ctauErrPVDZ",      &dimuonditrk_ctauErrPVDZ,        "dimuonditrk_ctauErrPVDZ/D");
          dimuonditrk_tree->Branch("dimuonditrk_countTksOfPVDZ",      &dimuonditrk_countTksOfPVDZ,        "dimuonditrk_countTksOfPVDZ/D");
          dimuonditrk_tree->Branch("dimuonditrk_vertexWeightDZ",      &dimuonditrk_vertexWeightDZ,        "dimuonditrk_vertexWeightDZ/D");
          dimuonditrk_tree->Branch("dimuonditrk_sumPTPVDZ",      &dimuonditrk_sumPTPVDZ,        "dimuonditrk_sumPTPVDZ/D");
          dimuonditrk_tree->Branch("dimuonditrk_mu1FromPVDZ",      &dimuonditrk_mu1FromPVDZ,        "dimuonditrk_mu1FromPVDZ/D");
          dimuonditrk_tree->Branch("dimuonditrk_mu2FromPVDZ",      &dimuonditrk_mu2FromPVDZ,        "dimuonditrk_mu2FromPVDZ/D");
          dimuonditrk_tree->Branch("dimuonditrk_tPFromPVDZ",      &dimuonditrk_tPFromPVDZ,        "dimuonditrk_tPFromPVDZ/D");
          dimuonditrk_tree->Branch("dimuonditrk_tMFromPVDZ",      &dimuonditrk_tMFromPVDZ,        "dimuonditrk_tMFromPVDZ/D");
          dimuonditrk_tree->Branch("dimuonditrk_mu1DZW",      &dimuonditrk_mu1DZW,        "dimuonditrk_mu1DZW/D");
          dimuonditrk_tree->Branch("dimuonditrk_mu2DZW",      &dimuonditrk_mu2DZW,        "dimuonditrk_mu2DZW/D");
          dimuonditrk_tree->Branch("dimuonditrk_tPDZW",      &dimuonditrk_tPDZW,        "dimuonditrk_tPDZW/D");
          dimuonditrk_tree->Branch("dimuonditrk_tMDZW",      &dimuonditrk_tMDZW,        "dimuonditrk_tMDZW/D");

          dimuonditrk_tree->Branch("dimuonditrk_cosAlphaBS",      &dimuonditrk_cosAlphaBS,        "dimuonditrk_cosAlphaBS/D");
          dimuonditrk_tree->Branch("dimuonditrk_ctauPVBS",      &dimuonditrk_ctauPVBS,        "dimuonditrk_ctauPVBS/D");
          dimuonditrk_tree->Branch("dimuonditrk_ctauErrPVBS",      &dimuonditrk_ctauErrPVBS,        "dimuonditrk_ctauErrPVBS/D");
          dimuonditrk_tree->Branch("dimuonditrk_countTksOfPVBS",      &dimuonditrk_countTksOfPVBS,        "dimuonditrk_countTksOfPVBS/D");
          dimuonditrk_tree->Branch("dimuonditrk_vertexWeightBS",      &dimuonditrk_vertexWeightBS,        "dimuonditrk_vertexWeightBS/D");
          dimuonditrk_tree->Branch("dimuonditrk_sumPTPVBS",      &dimuonditrk_sumPTPVBS,        "dimuonditrk_sumPTPVBS/D");
          dimuonditrk_tree->Branch("dimuonditrk_mu1FromPVBS",      &dimuonditrk_mu1FromPVBS,        "dimuonditrk_mu1FromPVBS/D");
          dimuonditrk_tree->Branch("dimuonditrk_mu2FromPVBS",      &dimuonditrk_mu2FromPVBS,        "dimuonditrk_mu2FromPVBS/D");
          dimuonditrk_tree->Branch("dimuonditrk_tPFromPVBS",      &dimuonditrk_tPFromPVBS,        "dimuonditrk_tPFromPVBS/D");
          dimuonditrk_tree->Branch("dimuonditrk_tMFromPVBS",      &dimuonditrk_tMFromPVBS,        "dimuonditrk_tMFromPVBS/D");
          dimuonditrk_tree->Branch("dimuonditrk_mu1BSW",      &dimuonditrk_mu1BSW,        "dimuonditrk_mu1BSW/D");
          dimuonditrk_tree->Branch("dimuonditrk_mu2BSW",      &dimuonditrk_mu2BSW,        "dimuonditrk_mu2BSW/D");
          dimuonditrk_tree->Branch("dimuonditrk_tPBSW",      &dimuonditrk_tPBSW,        "dimuonditrk_tPBSW/D");
          dimuonditrk_tree->Branch("dimuonditrk_tMBSW",      &dimuonditrk_tMBSW,        "dimuonditrk_tMBSW/D");

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

          dimuonditrk_tree->Branch("highKaon_eta",        &highKaon_eta,        "highKaon_eta/D");
          dimuonditrk_tree->Branch("highKaon_phi",        &highKaon_phi,        "highKaon_phi/D");
          dimuonditrk_tree->Branch("highKaon_dz",        &highKaon_dz,        "highKaon_dz/D");
          dimuonditrk_tree->Branch("highKaon_dxy",        &highKaon_dxy,        "highKaon_dxy/D");

          dimuonditrk_tree->Branch("lowKaon_eta",        &lowKaon_eta,        "lowKaon_eta/D");
          dimuonditrk_tree->Branch("lowKaon_phi",        &lowKaon_phi,        "lowKaon_phi/D");
          dimuonditrk_tree->Branch("lowKaon_dz",        &lowKaon_dz,        "lowKaon_dz/D");
          dimuonditrk_tree->Branch("lowKaon_dxy",        &lowKaon_dxy,        "lowKaon_dxy/D");

          dimuonditrk_tree->Branch("lowMuon_type",     &lowMuon_type,       "lowMuon_type/i");
          dimuonditrk_tree->Branch("highMuon_type",     &highMuon_type,       "highMuon_type/i");

          //Tracks Flags

          dimuonditrk_tree->Branch("highKaon_eta",        &highKaon_eta,        "highKaon_eta/D");
          dimuonditrk_tree->Branch("highKaon_phi",        &highKaon_phi,        "highKaon_phi/D");
          dimuonditrk_tree->Branch("highKaon_dz",        &highKaon_dz,        "highKaon_dz/D");
          dimuonditrk_tree->Branch("highKaon_dxy",        &highKaon_dxy,        "highKaon_dxy/D");

          dimuonditrk_tree->Branch("lowKaon_eta",        &lowKaon_eta,        "lowKaon_eta/D");
          dimuonditrk_tree->Branch("lowKaon_phi",        &lowKaon_phi,        "lowKaon_phi/D");
          dimuonditrk_tree->Branch("lowKaon_dz",        &lowKaon_dz,        "lowKaon_dz/D");
          dimuonditrk_tree->Branch("lowKaon_dxy",        &lowKaon_dxy,        "lowKaon_dxy/D");

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
        int pdgid_ = 0;

        if (isMC_ || OnlyGen_) {
           std::cout << "DiMuonRootupler::DiMuonRootupler: Dimuon id " << pdgid_ << std::endl;

           dimuonditrk_tree->Branch("gen_dimuonditrk_pdgId",  &gen_dimuonditrk_pdgId,     "gen_dimuonditrk_pdgId/I");
           dimuonditrk_tree->Branch("gen_dimuonditrk_p4", "TLorentzVector",  &gen_dimuonditrk_p4);
           dimuonditrk_tree->Branch("gen_b_p4", "TLorentzVector",  &gen_b_p4);
           dimuonditrk_tree->Branch("gen_dimuon_p4", "TLorentzVector",  &gen_dimuon_p4);
           dimuonditrk_tree->Branch("gen_ditrak_p4", "TLorentzVector",  &gen_ditrak_p4);
           dimuonditrk_tree->Branch("gen_lowMuon_p4",  "TLorentzVector",  &gen_lowMuon_p4);
           dimuonditrk_tree->Branch("gen_highMuon_p4",  "TLorentzVector",  &gen_highMuon_p4);
           dimuonditrk_tree->Branch("gen_highKaon_p4",  "TLorentzVector",  &gen_highKaon_p4);
           dimuonditrk_tree->Branch("gen_lowKaon_p4",  "TLorentzVector",  &gen_lowKaon_p4);

           dimuonditrk_tree->Branch("gen_dimuonditrk_m",  &gen_dimuonditrk_m,    "gen_dimuonditrk_m/D");

           dimuonditrk_tree->Branch("dimuonditrk_pdgid",  &dimuonditrk_pdgid,    "dimuonditrk_pdgid/I");
           dimuonditrk_tree->Branch("dimuonditrk_jpsipdg",  &dimuonditrk_jpsipdg,    "dimuonditrk_jpsipdg/I");
           dimuonditrk_tree->Branch("dimuonditrk_isprompt",  &dimuonditrk_isprompt,    "dimuonditrk_isprompt/D");
           dimuonditrk_tree->Branch("dimuonditrk_jpsippdl",  &dimuonditrk_jpsippdl,    "dimuonditrk_jpsippdl/D");

           dimuonditrk_tree->Branch("gen_b4_p4", "TLorentzVector",  &gen_b4_p4);
           dimuonditrk_tree->Branch("gen_d1_p4",  "TLorentzVector",  &gen_d1_p4);
           dimuonditrk_tree->Branch("gen_d2_p4",  "TLorentzVector",  &gen_d2_p4);
           dimuonditrk_tree->Branch("gen_gd1_p4",  "TLorentzVector",  &gen_gd1_p4);
           dimuonditrk_tree->Branch("gen_gd2_p4",  "TLorentzVector",  &gen_gd2_p4);
           dimuonditrk_tree->Branch("gen_gd3_p4", "TLorentzVector",  &gen_gd3_p4);
           dimuonditrk_tree->Branch("gen_gd4_p4",  "TLorentzVector",  &gen_gd4_p4);
           dimuonditrk_tree->Branch("gen_gd5_p4",  "TLorentzVector",  &gen_gd5_p4);
           dimuonditrk_tree->Branch("gen_gd6_p4",  "TLorentzVector",  &gen_gd6_p4);

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
  iEvent.getByToken(dimuonditrk_cand_Label, dimuonditrk_cand_handle);

  edm::Handle<std::vector<pat::CompositeCandidate>> fivetracks_cand_hanlde;
  iEvent.getByToken(fivetracks_cand_Label, fivetracks_cand_hanlde);

  edm::Handle<std::vector<reco::Vertex >> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( triggerResults_Label , triggerResults_handle);

  numPrimaryVertices = primaryVertices_handle->size();
  run = iEvent.id().run();
  event = iEvent.id().event();

  reco::Vertex theBeamSpotV;

  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(thebeamspot_,theBeamSpot);
  reco::BeamSpot bs = *theBeamSpot;

  const reco::Vertex thePrimaryV = *dimuonditrakCand.userData<reco::Vertex>("bestPV");

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
//   gen_highMuon_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
//   foundit++;
// }
// if ( motherInPrunedCollection != nullptr && (d->pdgId() == -13 ) && isAncestor(aditrkdimu , motherInPrunedCollection) ) {
//   gen_lowMuon_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
//   foundit++;
// }

gen_dimuonditrk_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_ditrak_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_lowMuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_highMuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_highKaon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
gen_lowKaon_p4.SetPtEtaPhiM(0.,0.,0.,0.);

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

fivetraks_p4.SetPtEtaPhiM(0.,0.,0.,0.);
fivetrakspion_p4.SetPtEtaPhiM(0.,0.,0.,0.);
fifthkaon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
dimuontrak_p4.SetPtEtaPhiM(0.,0.,0.,0.);
dimuontrakpion_p4.SetPtEtaPhiM(0.,0.,0.,0.);

dimuonditrk_p4.SetPtEtaPhiM(0.,0.,0.,0.);
ditrak_p4.SetPtEtaPhiM(0.,0.,0.,0.);
dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
lowMuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
highMuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
highKaon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
lowKaon_p4.SetPtEtaPhiM(0.,0.,0.,0.);

dimuonditrk_pdgid      = 0;
dimuonditrk_isprompt   = -99.0;
dimuonditrk_jpsipdg    = 0;
dimuonditrk_jpsippdl   = -99.0;

gen_dimuonditrk_pdgId = 0;

if ( (isMC_ || OnlyGen_) && packed.isValid() && pruned.isValid() ) {
  for (size_t i=0; i<pruned->size(); i++) {
    // std::cout << "Valid"<<std::endl;
    const reco::Candidate *aditrkdimu = &(*pruned)[i];

    if ( (abs(aditrkdimu->pdgId()) == motherpdgid_) && (aditrkdimu->status() == 2))
      gen_b_p4.SetPtEtaPhiM(aditrkdimu->pt(),aditrkdimu->eta(),aditrkdimu->phi(),aditrkdimu->mass());

    if ( (abs(aditrkdimu->pdgId()) == motherpdgid_) && (aditrkdimu->status() == 2) && (aditrkdimu->numberOfDaughters() > 1) && (aditrkdimu->numberOfDaughters() < 7) ) {
      //asking for decay (status==2) && two daughters
      bool goToJPsi = false;
      bool goToPhi = false;

      bool muP = false, muN = false, kP = false, kN = false;

      int noDaughters = 0, noGDaughters = 0;
      int theJPsi = 0, thePhi = 0, theMuP = 0, theMuN = 0, theKP = 0, theKN = 0;

      std::vector<const reco::Candidate *> daughters,gdaughters;
      for(size_t j = 0; j < aditrkdimu->numberOfDaughters(); ++j)
      {
        const reco::Candidate * daughter = aditrkdimu->daughter(j);
        if(daughter->mother(daughter->numberOfMothers()-1) != aditrkdimu) continue;
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
              theMuP = j;
              muP=true;
            }
            if(gdaughter->pdgId()==13)
            {
              theMuN = j;
              muN=true;
            }
            if(gdaughter->pdgId()==321)
            {
              theKP = j;
              kP=true;
            }
            if(gdaughter->pdgId()==-321)
            {
              theKN = j;
              kN=true;
            }
          }

          ++noGDaughters;
        }
      }

      if(noDaughters == 2 && noGDaughters > 3 && noGDaughters < 7 && goToJPsi && goToPhi)
      {

        // for (size_t j = 0; j < daughters.size(); j++)
        //   std::cout << "Daughter no. " << j << " - id : " << daughters[j]->pdgId() << std::endl;
        //
        // for (size_t j = 0; j < gdaughters.size(); j++)
        //   std::cout << "GrandDaughter no. " << j << " - id : " << gdaughters[j]->pdgId() << std::endl;

        gen_b4_p4.SetPtEtaPhiM(aditrkdimu->pt(),aditrkdimu->eta(),aditrkdimu->phi(),aditrkdimu->mass());
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

      if(muP && muN && kP && kN)
      {

        gen_dimuonditrk_p4.SetPtEtaPhiM(aditrkdimu->pt(),aditrkdimu->eta(),aditrkdimu->phi(),aditrkdimu->mass());
        gen_dimuon_p4.SetPtEtaPhiM(daughters[theJPsi]->pt(),daughters[theJPsi]->eta(),daughters[theJPsi]->phi(),daughters[theJPsi]->mass());
        gen_ditrak_p4.SetPtEtaPhiM(daughters[thePhi]->pt(),daughters[thePhi]->eta(),daughters[thePhi]->phi(),daughters[thePhi]->mass());
        gen_highMuon_p4.SetPtEtaPhiM(gdaughters[theMuN]->pt(),gdaughters[theMuN]->eta(),gdaughters[theMuN]->phi(),gdaughters[theMuN]->mass());
        gen_lowMuon_p4.SetPtEtaPhiM(gdaughters[theMuP]->pt(),gdaughters[theMuP]->eta(),gdaughters[theMuP]->phi(),gdaughters[theMuP]->mass());
        gen_lowKaon_p4.SetPtEtaPhiM(gdaughters[theKN]->pt(),gdaughters[theKN]->eta(),gdaughters[theKN]->phi(),gdaughters[theKN]->mass());
        gen_highKaon_p4.SetPtEtaPhiM(gdaughters[theKP]->pt(),gdaughters[theKP]->eta(),gdaughters[theKP]->phi(),gdaughters[theKP]->mass());
        gen_dimuonditrk_pdgId = aditrkdimu->pdgId();
    }
  } // for (size
}
}  // end if isMC

if(OnlyGen_) dimuonditrk_tree->Fill();

std::map <int,const pat::CompositeCandidate*> fourToFiveMap;

if(!OnlyGen_)
{
  if (!dimuonditrk_cand_handle.isValid()) std::cout<< "No dimuontt information " << run << "," << event <<std::endl;
  if (!fivetracks_cand_hanlde.isValid()) std::cout<< "No fivetrack information " << run << "," << event <<std::endl;
  if (dimuonditrk_cand_handle.isValid() && fivetracks_cand_hanlde.isValid())
  {
    for (unsigned int i=0; i< fivetracks_cand_hanlde->size(); i++)
      fourToFiveMap[fivetracks_cand_hanlde->at(i).userInt("index")] = &(fivetracks_cand_hanlde->at(i));
  }
}

// get rf information. Notice we are just keeping combinations with succesfull vertex fit
  if(!OnlyGen_)
  {
  if (!dimuonditrk_cand_handle.isValid()) std::cout<< "No dimuontt information " << run << "," << event <<std::endl;
  if (!fivetracks_cand_hanlde.isValid()) std::cout<< "No fivetrack information " << run << "," << event <<std::endl;
  if (dimuonditrk_cand_handle.isValid() && fivetracks_cand_hanlde.isValid()) {

    pat::CompositeCandidate *dimuonditrk_rf_cand, dimuonditrk_cand, *dimuon_cand;
    pat::CompositeCandidate *ditrak_cand, *dimuon_cand_rf, *ditrak_cand_rf;
    pat::CompositeCandidate *fivetrak_cand, *dimuontrak_cand,*fivetrakpion_cand, *dimuontrakpion_cand;

    noXCandidates = (Int_t)(dimuonditrk_cand_handle->size());
    //UnRefitted Handle
    for (unsigned int i=0; i< dimuonditrk_cand_handle->size(); i++){

      dimuonditrk_cand   = dimuonditrk_cand_handle->at(i);

      fivetraks_kaon_m    = -1.0;
      fivetraks_pion_m    = -1.0;
      fivetraks_kaon_trim    = -1.0;
      fivetraks_pion_trim    = -1.0;
      fivetraks_kaon_m_rf    = -1.0;
      fivetraks_pion_m_rf    = -1.0;
      fivetraks_vProb    = -1.0;
      fivetraks_vChi2    = -1.0;
      fivetraks_nDof    = -1.0;
      fivetraks_charge    = -1.0;
      fivetraks_cosAlpha    = -1.0;
      fivetraks_ctauPV    = -1.0;
      fivetraks_ctauErrPV    = -1.0;
      fivetraks_eta    = -1.0;
      fivetraks_pt    = -1.0;
      fivetraks_phi    = -1.0;
      fivetraks_y    = -1.0;
      fifthtrak_charge = -1.0;
      fifthtrak_eta = -1.0;
      fifthtrak_pt = -1.0;
      fifthtrak_phi = -1.0;
      fifthtrak_y = -1.0;

      if(std::find(fourToFiveMap.begin(),fourToFiveMap.end(),i)!=fourToFiveMap.end())
      {
        fivetrak_cand = fourToFiveMap[i];
        dimuontrak_cand = dynamic_cast <pat::CompositeCandidate *>(dimuonditrk_cand.daughter("dimuontrak"));
        fivetrakpion_cand = dynamic_cast <pat::CompositeCandidate *>(dimuonditrk_cand.daughter("withpion"));
        dimuontrakpion_cand = dynamic_cast <pat::CompositeCandidate *>(dimuonditrk_cand.daughter("dimuontrakpion"));

        fivetraks_p4.SetPtEtaPhiM(fivetrak_cand.pt(), fivetrak_cand.eta(), fivetrak_cand.phi(), fivetrak_cand.mass());
        fivetrakspion_p4.SetPtEtaPhiM(dimuontrak_cand.pt(), dimuontrak_cand.eta(), dimuontrak_cand.phi(), dimuontrak_cand.mass());
        dimuontrak_p4.SetPtEtaPhiM(fivetrakpion_cand.pt(), fivetrakpion_cand.eta(), fivetrakpion_cand.phi(), fivetrakpion_cand.mass());
        dimuontrakpion_p4.SetPtEtaPhiM(dimuontrakpion_cand.pt(), dimuontrakpion_cand.eta(), dimuontrakpion_cand.phi(), dimuontrakpion_cand.mass());

        reco::Candidate::LorentzVector fifth = fivetrak_cand->daughter("fifth")->p4();
        fifthkaon_p4.SetPtEtaPhiM(fifth.pt(), fifth.eta(), fifth.phi(), fifth.mass());

        fivetraks_kaon_m    = fivetrak_cand.mass();
        fivetraks_pion_m    = fivetrak_cand.userFloat("mass_pion");
        fivetraks_kaon_trim    = fivetrak_cand.userFloat("trim");
        fivetraks_pion_trim    = fivetrak_cand.userFloat("trim_pion");
        fivetraks_kaon_m_rf    = fivetrak_cand.userFloat("mass_kaon_rf");
        fivetraks_pion_m_rf    = fivetrak_cand.userFloat("mass_pion_rf");
        fivetraks_vProb    = fivetrak_cand.userFloat("vProb");
        fivetraks_vChi2    = fivetrak_cand.userFloat("vChi2");
        fivetraks_nDof    = fivetrak_cand.userFloat("nDof");
        fivetraks_charge    = fivetrak_cand.mass();
        fivetraks_cosAlpha    = fivetrak_cand.userFloat("cosAlpha");
        fivetraks_ctauPV    = fivetrak_cand.userFloat("ctauPV");
        fivetraks_ctauErrPV    = fivetrak_cand.userFloat("ctauErrPV");
        fivetraks_eta    = fivetrak_cand.eta();
        fivetraks_pt    = fivetrak_cand.pt();
        fivetraks_phi    = fivetrak_cand.phi();
        fivetraks_y    = fivetrak_cand.y();

        fifthtrak_charge = fivetrak_cand->daughter("fifth")->charge();
        fifthtrak_eta = fifth.Eta();
        fifthtrak_pt = fifth.Pt();
        fifthtrak_phi = fifth.Phi();
        fifthtrak_y = fifth.Rapidity();
      }

      dimuonditrk_vProb     = dimuonditrk_cand.userFloat("vProb");
      dimuonditrk_vChi2     = dimuonditrk_cand.userFloat("vChi2");
      dimuonditrk_nDof      = dimuonditrk_cand.userFloat("nDof");
      dimuonditrk_charge    = dimuonditrk_cand.charge();

      dimuonditrk_cosAlphaBS = dimuonditrk_cand.userFloat("cosAlphaBS");
      dimuonditrk_ctauPVBS = dimuonditrk_cand.userFloat("ctauPVBS");
      dimuonditrk_ctauErrPVBS = dimuonditrk_cand.userFloat("ctauErrPVBS");
      dimuonditrk_countTksOfPVBS = dimuonditrk_cand.userFloat("countTksOfPVBS");
      dimuonditrk_vertexWeightBS = dimuonditrk_cand.userFloat("vertexWeightBS");
      dimuonditrk_sumPTPVBS = dimuonditrk_cand.userFloat("sumPTPVBS");
      dimuonditrk_mu1FromPVBS = dimuonditrk_cand.userFloat("mu1FromPVBS");
      dimuonditrk_mu2FromPVBS = dimuonditrk_cand.userFloat("mu2FromPVBS");
      dimuonditrk_tPFromPVBS = dimuonditrk_cand.userFloat("tPFromPVBS");
      dimuonditrk_tMFromPVBS = dimuonditrk_cand.userFloat("tMFromPVBS");
      dimuonditrk_mu1BSW = dimuonditrk_cand.userFloat("mu1BSW");
      dimuonditrk_mu2BSW = dimuonditrk_cand.userFloat("mu2BSW");
      dimuonditrk_tPBSW = dimuonditrk_cand.userFloat("tPBSW");
      dimuonditrk_tMBSW = dimuonditrk_cand.userFloat("tMBSW");

      dimuonditrk_cosAlpha = dimuonditrk_cand.userFloat("cosAlpha");
      dimuonditrk_ctauPV = dimuonditrk_cand.userFloat("ctauPV");
      dimuonditrk_ctauErrPV = dimuonditrk_cand.userFloat("ctauErrPV");
      dimuonditrk_countTksOfPV = dimuonditrk_cand.userFloat("countTksOfPV");
      dimuonditrk_vertexWeight = dimuonditrk_cand.userFloat("vertexWeight");
      dimuonditrk_sumPTPV = dimuonditrk_cand.userFloat("sumPTPV");
      dimuonditrk_mu1FromPV = dimuonditrk_cand.userFloat("mu1FromPV");
      dimuonditrk_mu2FromPV = dimuonditrk_cand.userFloat("mu2FromPV");
      dimuonditrk_tPFromPV = dimuonditrk_cand.userFloat("tPFromPV");
      dimuonditrk_tMFromPV = dimuonditrk_cand.userFloat("tMFromPV");
      dimuonditrk_mu1W = dimuonditrk_cand.userFloat("mu1W");
      dimuonditrk_mu2W = dimuonditrk_cand.userFloat("mu2W");
      dimuonditrk_tPW = dimuonditrk_cand.userFloat("tPW");
      dimuonditrk_tMW = dimuonditrk_cand.userFloat("tMW");

      dimuonditrk_cosAlphaDZ = dimuonditrk_cand.userFloat("cosAlphaDZ");
      dimuonditrk_ctauPVDZ = dimuonditrk_cand.userFloat("ctauPVDZ");
      dimuonditrk_ctauErrPVDZ = dimuonditrk_cand.userFloat("ctauErrPVDZ");
      dimuonditrk_countTksOfPVDZ = dimuonditrk_cand.userFloat("countTksOfPVDZ");
      dimuonditrk_vertexWeightDZ = dimuonditrk_cand.userFloat("vertexWeightDZ");
      dimuonditrk_sumPTPVDZ = dimuonditrk_cand.userFloat("sumPTPVDZ");
      dimuonditrk_mu1FromPVDZ = dimuonditrk_cand.userFloat("mu1FromPVDZ");
      dimuonditrk_mu2FromPVDZ = dimuonditrk_cand.userFloat("mu2FromPVDZ");
      dimuonditrk_tPFromPVDZ = dimuonditrk_cand.userFloat("tPFromPVDZ");
      dimuonditrk_tMFromPVDZ = dimuonditrk_cand.userFloat("tMFromPVDZ");
      dimuonditrk_mu1DZW = dimuonditrk_cand.userFloat("mu1DZW");
      dimuonditrk_mu1DZW = dimuonditrk_cand.userFloat("mu2DZW");
      dimuonditrk_tPDZW = dimuonditrk_cand.userFloat("tPDZW");
      dimuonditrk_tMDZW = dimuonditrk_cand.userFloat("tMDZW");

      dimuonditrk_dca_m1m2 = dimuonditrk_cand.userFloat("dca_m1m2");
      dimuonditrk_dca_m1t1 = dimuonditrk_cand.userFloat("dca_m1t1");
      dimuonditrk_dca_m1t2 = dimuonditrk_cand.userFloat("dca_m1t2");
      dimuonditrk_dca_m2t1 = dimuonditrk_cand.userFloat("dca_m2t1");
      dimuonditrk_dca_m2t2 = dimuonditrk_cand.userFloat("dca_m2t2");
      dimuonditrk_dca_t1t2 = dimuonditrk_cand.userFloat("dca_t1t2");


      if(isMC_ || OnlyGen_)
      {
        dimuonditrk_jpsipdg    = dimuonditrk_cand.userInt("jPsiGenPdgId");
        dimuonditrk_jpsippdl   = dimuonditrk_cand.userFloat("jPsiPpdlTrue");
        dimuonditrk_pdgid      = dimuonditrk_cand.userInt("xGenPdgId");
        dimuonditrk_isprompt   = dimuonditrk_cand.userFloat("xGenIsPrompt");
      }

      highKaonMatch = dimuonditrk_cand.userInt("highKaonMatch");
      lowKaonMatch = dimuonditrk_cand.userInt("lowKaonMatch");

      //unref corresponding

      dimuon_cand = dynamic_cast <pat::CompositeCandidate *>(dimuonditrk_cand.daughter("dimuon"));
      ditrak_cand = dynamic_cast <pat::CompositeCandidate *>(dimuonditrk_cand.daughter("ditrak"));

      const pat::Muon *lowMuon, *highMuon;

      reco::Candidate::LorentzVector vHighMuon = dimuon_cand->daughter("highMuon")->p4();
      reco::Candidate::LorentzVector vLowMuon = dimuon_cand->daughter("lowMuon")->p4();

      // if (dimuon_cand->daughter("highMuon")->charge() < 0) {
      if (vHighMuon.pt() < vLowMuon.pt()) {
         vHighMuon = dimuon_cand->daughter("lowMuon")->p4();
         vLowMuon = dimuon_cand->daughter("highMuon")->p4();
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

      lowMuon_p4.SetPtEtaPhiM(vHighMuon.pt(), vHighMuon.eta(), vHighMuon.phi(), vHighMuon.mass());
      highMuon_p4.SetPtEtaPhiM(vLowMuon.pt(), vLowMuon.eta(), vLowMuon.phi(), vLowMuon.mass());

      reco::Candidate::LorentzVector kP = ditrak_cand->daughter("highTrak")->p4();
      reco::Candidate::LorentzVector kM = ditrak_cand->daughter("lowTrak")->p4();

      pat::CompositeCandidate *highTrack = dynamic_cast <pat::CompositeCandidate *>(ditrak_cand->daughter("highTrak"));
      pat::CompositeCandidate *lowTrack  = dynamic_cast <pat::CompositeCandidate *>(ditrak_cand->daughter("lowTrak"));

      highKaon_p4.SetPtEtaPhiM(kP.pt(), kP.eta(), kP.phi(), kP.mass());
      lowKaon_p4.SetPtEtaPhiM(kM.pt(), kM.eta(), kM.phi(), kM.mass());

      highKaon_pt     = std::max(kP.pt(),kM.pt());
      lowKaon_pt      = -std::max(-kP.pt(),-kM.pt());
      highMuon_pt     = std::max(vLowMuon.pt(),vHighMuon.pt());
      lowMuon_pt      = -std::max(-vLowMuon.pt(),-vHighMuon.pt());

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
      dimuon_triggerMatch = DiMuonDiTrakFiveRootuplerFit::isTriggerMatched(dimuon_cand);
      lowMuonMatch   = dimuon_cand->userInt("highMuonTMatch");
      highMuonMatch    = dimuon_cand->userInt("lowMuonTMatch");

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

        vHighMuon = dimuon_cand_rf->daughter("highMuon")->p4();
        vLowMuon = dimuon_cand_rf->daughter("lowMuon")->p4();

        lowMuon_rf_p4.SetPtEtaPhiM(vHighMuon.pt(), vHighMuon.eta(), vHighMuon.phi(), vHighMuon.mass());
        highMuon_rf_p4.SetPtEtaPhiM(vLowMuon.pt(), vLowMuon.eta(), vLowMuon.phi(), vLowMuon.mass());

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

      dimuonditrk_tree->Fill();

      if (OnlyBest_) break;
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