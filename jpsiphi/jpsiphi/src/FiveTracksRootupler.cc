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


  std::vector < float > fiveTracksMass, fiveTracksVProb, fiveTracksCTau, fiveTracksCTauErr, fiveTracksCosAlpha, PsiPrimeOne, PsiPrimeTwo;

  UInt_t run,event,numPrimaryVertices, trigger, dimuonditrak_id;

  TLorentzVector dimuonditrk_p4;
  TLorentzVector dimuon_p4;
  TLorentzVector ditrak_p4;
  TLorentzVector muonp_p4;
  TLorentzVector muonn_p4;
  TLorentzVector kaonp_p4;
  TLorentzVector kaonn_p4;

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

  Double_t dimuonditrk_vProb,  dimuonditrk_vChi2, dimuonditrk_cosAlpha, dimuonditrk_ctauPV, dimuonditrk_ctauErrPV;

  Double_t dimuon_vProb, dimuon_vChi2, dimuon_DCA, dimuon_ctauPV, dimuon_ctauErrPV, dimuon_cosAlpha;

  Double_t gen_dimuonditrk_m,dimuonditrk_m,dimuonditrk_pt,dimuon_m,dimuon_pt,ditrak_m,ditrak_pt;
  Double_t highKaon_pt,lowKaon_pt,highMuon_pt,lowMuon_pt,dimuonditrk_nDof,dimuonditrk_m_rf;

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

  UInt_t tPMatch, tNMatch;

  Int_t track_KN_nvsh, track_KN_nvph;

  Int_t dimuonditrk_rf_bindx;

  Int_t noXCandidates;

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

  TTree* fivetracks_tree, *fivetracks_tree_rf;
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
        fivetracks_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");
        fivetracks_tree->Branch("trigger",            &trigger,            "trigger/I");
        fivetracks_tree->Branch("noXCandidates",      &noXCandidates,      "noXCandidates/I");

        fivetracks_tree->Branch("dimuonditrak_id",      &dimuonditrak_id,      "dimuonditrak_id/I");
        fivetracks_tree->Branch("dimuonditrak_p4",     "TLorentzVector", &dimuonditrak_p4);
        fivetracks_tree->Branch("dimuonditrk_m",       &dimuonditrk_m,        "dimuonditrk_m/D");
        fivetracks_tree->Branch("dimuonditrk_pt",          &dimuonditrk_pt,          "dimuonditrk_pt/D");

        fivetracks_tree->Branch("dimuon_m",       &dimuon_m,       "dimuon_m/D");
        fivetracks_tree->Branch("dimuon_pt",    &dimuon_pt,    "dimuon_pt/D");

        //The kinematic doesn't change, only mass
        fivetracks_tree->Branch("highTrack_pt",          &highTrack_pt,          "highTrack_pt/D");
        fivetracks_tree->Branch("highTrack_eta",        &highTrack_eta,        "highTrack_eta/D");
        fivetracks_tree->Branch("highTrack_phi",        &highTrack_phi,        "highTrack_phi/D");

        fivetracks_tree->Branch("lowTrack_pt",          &lowTrack_pt,          "lowTrack_pt/D");
        fivetracks_tree->Branch("lowTrack_eta",        &lowTrack_eta,        "lowTrack_eta/D");
        fivetracks_tree->Branch("lowTrack_phi",        &lowTrack_phi,        "lowTrack_phi/D");

        fivetracks_tree->Branch("thirdTrack_pt",          &thirdTrack_pt,          "thirdTrack_pt/D");
        fivetracks_tree->Branch("thirdTrack_eta",        &thirdTrack_eta,        "thirdTrack_eta/D");
        fivetracks_tree->Branch("thirdTrack_phi",        &thirdTrack_phi,        "thirdTrack_phi/D");

        numMasses = 5;

        for(size_t i = 0; i<numMasses;i++)
        {
          fiveTracksMass.push_back(2.5);
          fiveTracksVNDof.push_back(-1.0);
          fiveTracksVChi2.push_back(-1.0);
          fiveTracksVProb.push_back(-0.1);
          fiveTracksCTau.push_back(-1000.0);
          fiveTracksCTauErr.push_back(-1000.0);
          fiveTracksCosAlpha.push_back(-1.1);
        }

        std::vector < std::string > refNames;
        refNames.push_back("kkk"); refNames.push_back("ppk"); refNames.push_back("kpp");
        refNames.push_back("pkp"); refNames.push_back("ppp");


        for(size_t i = 0; i<numMasses;i++)
        {
          std::string name = "mass_" + refNames[i];
          std::string var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&fiveTracksMass[i],var.c_str());
          name = "vProb_" + refNames[i]; var = name + "/D";

          fivetracks_tree->Branch(name.c_str(),&fiveTracksVProb[i],var.c_str());
          name = "ctau_" + refNames[i]; var = name + "/D";

          fivetracks_tree->Branch(name.c_str(),&fiveTracksCTau[i],var.c_str());
          name = "nDof_" + refNames[i]; var = name + "/D";

          fivetracks_tree->Branch(name.c_str(),&fiveTracksCTau[i],var.c_str());
          name = "vChi2_" + refNames[i]; var = name + "/D";

          fivetracks_tree->Branch(name.c_str(),&fiveTracksCTau[i],var.c_str());
          name = "ctauErr_" + refNames[i]; var = name + "/D";

          fivetracks_tree->Branch(name.c_str(),&fiveTracksCTauErr[i],var.c_str());
          name = "cosAlpha_" + refNames[i]; var = name + "/D";

          fivetracks_tree->Branch(name.c_str(),&fiveTracksCosAlpha[i],var.c_str());
          name = "onePsiPrime_" + refNames[i]; var = name + "/D";

          fivetracks_tree->Branch(name.c_str(),&PsiPrimeOne[i],var.c_str());
          name = "twoPsiPrime_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&PsiPrimeTwo[i],var.c_str());

          fivetracks_tree->Branch("five_p4",   "TLorentzVector", &five_p4);

          fivetracks_tree->Branch("psiPrimeOne_p4",     "TLorentzVector", &psiPrimeOne_p4);
          fivetracks_tree->Branch("psiPrimeTwo_p4",     "TLorentzVector", &psiPrimeTwo_p4);


        }

        // //p4s
        // fivetracks_tree->Branch("dimuonditrk_p4",   "TLorentzVector", &dimuonditrk_p4);
        // fivetracks_tree->Branch("ditrak_p4",     "TLorentzVector", &ditrak_p4);
        // fivetracks_tree->Branch("dimuon_p4",     "TLorentzVector", &dimuon_p4);
        // fivetracks_tree->Branch("muonp_p4",   "TLorentzVector", &muonp_p4);
        // fivetracks_tree->Branch("muonn_p4",   "TLorentzVector", &muonn_p4);
        // fivetracks_tree->Branch("kaonp_p4",   "TLorentzVector", &kaonp_p4);
        // fivetracks_tree->Branch("kaonn_p4",   "TLorentzVector", &kaonn_p4);
        //
        // //refitted p4s
        // fivetracks_tree->Branch("dimuonditrk_rf_p4",   "TLorentzVector", &dimuonditrk_rf_p4);
        // fivetracks_tree->Branch("ditrak_rf_p4",     "TLorentzVector", &ditrak_rf_p4);
        // fivetracks_tree->Branch("dimuon_rf_p4",     "TLorentzVector", &dimuon_rf_p4);
        // fivetracks_tree->Branch("muonp_rf_p4",   "TLorentzVector", &muonp_rf_p4);
        // fivetracks_tree->Branch("muonn_rf_p4",   "TLorentzVector", &muonn_rf_p4);
        // fivetracks_tree->Branch("kaonp_rf_p4",   "TLorentzVector", &kaonp_rf_p4);
        // fivetracks_tree->Branch("kaonn_rf_p4",   "TLorentzVector", &kaonn_rf_p4);
        //
        // //kin
        // fivetracks_tree->Branch("gen_dimuonditrk_m",        &gen_dimuonditrk_m,        "gen_dimuonditrk_m/D");
        // fivetracks_tree->Branch("dimuonditrk_m",       &dimuonditrk_m,        "dimuonditrk_m/D");
        // fivetracks_tree->Branch("dimuonditrk_m_rf",       &dimuonditrk_m_rf,        "dimuonditrk_m_rf/D");
        // fivetracks_tree->Branch("dimuonditrk_pt",          &dimuonditrk_pt,          "dimuonditrk_pt/D");
        // fivetracks_tree->Branch("dimuon_m",       &dimuon_m,       "dimuon_m/D");
        // fivetracks_tree->Branch("dimuon_pt",    &dimuon_pt,    "dimuon_pt/D");
        // fivetracks_tree->Branch("ditrak_m",     &ditrak_m,     "ditrak_m/D");
        // fivetracks_tree->Branch("ditrak_pt",       &ditrak_pt,        "ditrak_pt/D");
        // fivetracks_tree->Branch("highKaon_pt",          &highKaon_pt,          "highKaon_pt/D");
        // fivetracks_tree->Branch("lowKaon_pt",       &lowKaon_pt,       "lowKaon_pt/D");
        // fivetracks_tree->Branch("highMuon_pt",    &highMuon_pt,    "highMuon_pt/D");
        // fivetracks_tree->Branch("lowMuon_pt",     &lowMuon_pt,     "lowMuon_pt/D");
        //
        // //2mu vertexing
        // fivetracks_tree->Branch("dimuon_vProb",        &dimuon_vProb,        "dimuon_vProb/D");
        // fivetracks_tree->Branch("dimuon_vNChi2",       &dimuon_vChi2,        "dimuon_vNChi2/D");
        // fivetracks_tree->Branch("dimuon_DCA",          &dimuon_DCA,          "dimuon_DCA/D");
        // fivetracks_tree->Branch("dimuon_ctauPV",       &dimuon_ctauPV,       "dimuon_ctauPV/D");
        // fivetracks_tree->Branch("dimuon_ctauErrPV",    &dimuon_ctauErrPV,    "dimuon_ctauErrPV/D");
        // fivetracks_tree->Branch("dimuon_cosAlpha",     &dimuon_cosAlpha,     "dimuon_cosAlpha/D");
        // fivetracks_tree->Branch("dimuon_triggerMatch", &dimuon_triggerMatch, "dimuon_triggerMatch/I");
        //
        // //2mu+2Trk vertexing
        // fivetracks_tree->Branch("dimuonditrk_vProb",      &dimuonditrk_vProb,        "dimuonditrk_vProb/D");
        // fivetracks_tree->Branch("dimuonditrk_vChi2",      &dimuonditrk_vChi2,        "dimuonditrk_vChi2/D");
        // fivetracks_tree->Branch("dimuonditrk_nDof",       &dimuonditrk_nDof,         "dimuonditrk_nDof/D");
        // fivetracks_tree->Branch("dimuonditrk_cosAlpha",   &dimuonditrk_cosAlpha,     "dimuonditrk_cosAlpha/D");
        // fivetracks_tree->Branch("dimuonditrk_ctauPV",     &dimuonditrk_ctauPV,       "dimuonditrk_ctauPV/D");
        // fivetracks_tree->Branch("dimuonditrk_ctauErrPV",  &dimuonditrk_ctauErrPV,    "dimuonditrk_ctauErrPV/D");
        // fivetracks_tree->Branch("dimuonditrk_charge",     &dimuonditrk_charge,       "dimuonditrk_charge/I");
        //
        // fivetracks_tree->Branch("tPMatch",     &tPMatch,       "tPMatch/I");
        // fivetracks_tree->Branch("tNMatch",     &tNMatch,       "tNMatch/I");
        // //Muon flags
        // fivetracks_tree->Branch("muonP_isLoose",        &muonP_isLoose,        "muonP_isLoose/O");
        // fivetracks_tree->Branch("muonP_isSoft",        &muonP_isSoft,        "muonP_isSoft/O");
        // fivetracks_tree->Branch("muonP_isMedium",        &muonP_isMedium,        "muonP_isMedium/O");
        // fivetracks_tree->Branch("muonP_isHighPt",        &muonP_isHighPt,        "muonP_isHighPt/O");
        //
        // fivetracks_tree->Branch("muonP_isTracker",        &muonP_isTracker,        "muonP_isTracker/O");
        // fivetracks_tree->Branch("muonP_isGlobal",        &muonP_isGlobal,        "muonP_isGlobal/O");
        //
        // fivetracks_tree->Branch("muonN_isLoose",        &muonN_isLoose,        "muonN_isLoose/O");
        // fivetracks_tree->Branch("muonN_isSoft",        &muonN_isSoft,        "muonN_isSoft/O");
        // fivetracks_tree->Branch("muonN_isMedium",        &muonN_isMedium,        "muonN_isMedium/O");
        // fivetracks_tree->Branch("muonN_isHighPt",        &muonN_isHighPt,        "muonN_isHighPt/O");
        //
        // fivetracks_tree->Branch("muonN_isTracker",        &muonN_isTracker,        "muonN_isTracker/O");
        // fivetracks_tree->Branch("muonN_isGlobal",        &muonN_isGlobal,        "muonN_isGlobal/O");
        //
        // fivetracks_tree->Branch("muonP_type",     &muonP_type,       "muonP_type/i");
        // fivetracks_tree->Branch("muonN_type",     &muonN_type,       "muonN_type/i");
        //
        // int pdgid_ = 0;
        //
        // if (isMC_ ) {
        //    std::cout << "DiMuonRootupler::DiMuonRootupler: Dimuon id " << pdgid_ << std::endl;
        //    fivetracks_tree->Branch("gen_dimuonditrk_pdgId",  &gen_dimuonditrk_pdgId,     "gen_dimuonditrk_pdgId/I");
        //    fivetracks_tree->Branch("gen_dimuonditrk_p4", "TLorentzVector",  &gen_dimuonditrk_p4);
        //    fivetracks_tree->Branch("gen_b_p4", "TLorentzVector",  &gen_b_p4);
        //    fivetracks_tree->Branch("gen_dimuon_p4", "TLorentzVector",  &gen_dimuon_p4);
        //    fivetracks_tree->Branch("gen_ditrak_p4", "TLorentzVector",  &gen_ditrak_p4);
        //    fivetracks_tree->Branch("gen_muonp_p4",  "TLorentzVector",  &gen_muonp_p4);
        //    fivetracks_tree->Branch("gen_muonn_p4",  "TLorentzVector",  &gen_muonn_p4);
        //    fivetracks_tree->Branch("gen_kaonp_p4",  "TLorentzVector",  &gen_kaonp_p4);
        //    fivetracks_tree->Branch("gen_kaonn_p4",  "TLorentzVector",  &gen_kaonn_p4);
        //
        //    fivetracks_tree->Branch("gen_dimuonditrk_m",  &gen_dimuonditrk_m,    "gen_dimuonditrk_m/D");
        //
        //    fivetracks_tree->Branch("gen_b4_p4", "TLorentzVector",  &gen_b4_p4);
        //    fivetracks_tree->Branch("gen_d1_p4",  "TLorentzVector",  &gen_d1_p4);
        //    fivetracks_tree->Branch("gen_d2_p4",  "TLorentzVector",  &gen_d2_p4);
        //    fivetracks_tree->Branch("gen_gd1_p4",  "TLorentzVector",  &gen_gd1_p4);
        //    fivetracks_tree->Branch("gen_gd2_p4",  "TLorentzVector",  &gen_gd2_p4);
        //    fivetracks_tree->Branch("gen_gd3_p4", "TLorentzVector",  &gen_gd3_p4);
        //    fivetracks_tree->Branch("gen_gd4_p4",  "TLorentzVector",  &gen_gd4_p4);
        //    fivetracks_tree->Branch("gen_gd5_p4",  "TLorentzVector",  &gen_gd5_p4);
        //    fivetracks_tree->Branch("gen_gd6_p4",  "TLorentzVector",  &gen_gd6_p4);
        //
        // }
        //
        // //Track flags
        //
        //
        // fivetracks_tree->Branch("isBestCandidate",        &isBestCandidate,        "isBestCandidate/O");

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

if(!OnlyGen_)
  if (!fivetracks_cand_handle.isValid()) std::cout<< "No five tracks information " << run << "," << event <<std::endl;

  if (fivetracks_cand_handle.isValid()) {

    pat::CompositeCandidate dimuonditrk_rf_cand, fivetracks_cand, *dimuon_cand, *ditrak_cand, *dimuon_cand_rf, *ditrak_cand_rf;
    for (unsigned int i=0; i< fivetracks_cand_handle->size(); i++)
    {
      fivetracks_cand  = fivetracks_cand_handle->at(i);


      fiveCandKaon.addUserInt("dimuontt_index",d);

      std::string name;
      for(size_t i = 0; i<numMasses;i++)
      {

        fiveCands[i].addUserFloat("vProb",fiveTracksVProb[i]);
        fiveCands[i].addUserFloat("vChi2",fiveTracksCTau[i]);
        fiveCands[i].addUserFloat("nDof",fiveTracksCTauErr[i]);
        fiveCands[i].addUserFloat("cosAlpha",fiveTracksCosAlpha[i]);
        fiveCands[i].addUserFloat("ctauPV",fiveTracksVNDof[i]);
        fiveCands[i].addUserFloat("ctauErrPV",fiveTracksVChi2[i]);

        name = "fiveCand_" + std::to_string(i);

        fiveCandKaon.addDaughter(fiveCands[i],name);
      }

      for(size_t i = 0; i<numMasses;i++)
      {
        std::string name = "mass_" + std::to_string(i);
        fiveTracksMass[i] = fivetracks_cand.userFloat(name);
        name = "vProb_" + std::to_string(i);
        fiveTracksVProb[i] = fivetracks_cand.userFloat(name);
        name = "ctau_" + std::to_string(i);
        fiveTracksCTau[i] = fivetracks_cand.userFloat(name);
        name = "ctauErr_" + std::to_string(i);
        fiveTracksCTauErr[i] = fivetracks_cand.userFloat(name);
        name = "cosAlpha_" + std::to_string(i);
        fiveTracksCosAlpha[i] = fivetracks_cand.userFloat(name);
        name = "onePsiPrime_" + std::to_string(i);
        PsiPrimeOne[i] = fivetracks_cand.userFloat(name);
        name = "twoPsiPrime_" + std::to_string(i);
        PsiPrimeTwo[i] = fivetracks_cand.userFloat(name);
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
