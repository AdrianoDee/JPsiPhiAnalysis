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

  std::vector < float > fiveTracksMass, fiveTracksMassRef,
  std::vector < float > fiveTracksVProb, fiveTracksCTau, fiveTracksCTauErr;
  std::vector < float > fiveTracksCosAlpha, psiPrimeOne, psiPrimeTwo, psiPrimeThree;
  std::vector < float > trackOneMass, trackTwoMass, trackThreeMass;

  std::vector < TLorentzVector > five_p4, five_ref_p4, psiPrimeOne_p4, psiPrimeTwo_p4;
  std::vector < TLorentzVector >  dimuonDiTrakOne_p4, dimuonDiTrakTwo_p4, dimuonDiTrakThree_p4;
  UInt_t run, event, lumi, numPrimaryVertices, trigger, dimuonditrak_id;

  TLorentzVector , dimuonditrk_p4, dimuon_p4, ditrak_p4;
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
        fivetracks_tree->Branch("lumi",              &lumi,              "lumi/I");
        fivetracks_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");
        fivetracks_tree->Branch("trigger",            &trigger,            "trigger/I");

        fivetracks_tree->Branch("noFiveCandidates",      &noFiveCandidates,      "noFiveCandidates/I");

        fivetracks_tree->Branch("dimuonditrak_id",      &dimuonditrak_id,      "dimuonditrak_id/I");
        fivetracks_tree->Branch("dimuonditrak_p4",     "TLorentzVector", &dimuonditrak_p4);

        fivetracks_tree->Branch("dimuonditrk_m",       &dimuonditrk_m,        "dimuonditrk_m/D");
        fivetracks_tree->Branch("dimuonditrk_pt",      &dimuonditrk_pt,       "dimuonditrk_pt/D");
        fivetracks_tree->Branch("dimuonditrk_eta",     &dimuonditrk_eta,      "dimuonditrk_eta/D");
        fivetracks_tree->Branch("dimuonditrk_phi",     &dimuonditrk_phi,      "dimuonditrk_phi/D");

        fivetracks_tree->Branch("dimuon_m",      &dimuon_m,     "dimuon_m/D");
        fivetracks_tree->Branch("dimuon_pt",     &dimuon_pt,    "dimuon_pt/D");
        fivetracks_tree->Branch("dimuon_eta",    &dimuon_pt,    "dimuon_eta/D");
        fivetracks_tree->Branch("dimuon_phi",    &dimuon_phi,   "dimuon_phi/D");

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

        //J/Psi TrTr system
        fivetracks_tree->Branch("dimuonDiTrakOne_pt",          &dimuonDiTrakOne_pt,          "dimuonDiTrakOne_pt/D");
        fivetracks_tree->Branch("dimuonDiTrakOne_eta",        &dimuonDiTrakOne_eta,        "dimuonDiTrakOne_eta/D");
        fivetracks_tree->Branch("dimuonDiTrakOne_phi",        &dimuonDiTrakOne_phi,        "dimuonDiTrakOne_phi/D");
        fivetracks_tree->Branch("dimuonDiTrakOne_charge",        &dimuonDiTrakOne_charge,        "dimuonDiTrakOne_charge/D");

        fivetracks_tree->Branch("dimuonDiTrakTwo_pt",          &dimuonDiTrakTwo_pt,          "dimuonDiTrakTwo_pt/D");
        fivetracks_tree->Branch("dimuonDiTrakTwo_eta",        &dimuonDiTrakTwo_eta,        "dimuonDiTrakTwo_eta/D");
        fivetracks_tree->Branch("dimuonDiTrakTwo_phi",        &dimuonDiTrakTwo_phi,        "dimuonDiTrakTwo_phi/D");
        fivetracks_tree->Branch("dimuonDiTrakTwo_charge",        &dimuonDiTrakTwo_charge,        "dimuonDiTrakTwo_charge/D");

        fivetracks_tree->Branch("dimuonDiTrakThree_pt",          &dimuonDiTrakThree_pt,          "dimuonDiTrakThree_pt/D");
        fivetracks_tree->Branch("dimuonDiTrakThree_eta",        &dimuonDiTrakThree_eta,        "dimuonDiTrakThree_eta/D");
        fivetracks_tree->Branch("dimuonDiTrakThree_phi",        &dimuonDiTrakThree_phi,        "dimuonDiTrakThree_phi/D");
        fivetracks_tree->Branch("dimuonDiTrakThree_charge",        &dimuonDiTrakThree_charge,        "dimuonDiTrakThree_charge/D");

        //PsiPrime
        //One
        fivetracks_tree->Branch("psiPrimeOne_pt",         &psiPrimeOne_pt,        "psiPrimeOne_pt/D");
        fivetracks_tree->Branch("psiPrimeOne_eta",        &psiPrimeOne_eta,       "psiPrimeOne_eta/D");
        fivetracks_tree->Branch("psiPrimeOne_phi",        &psiPrimeOne_phi,       "psiPrimeOne_phi/D");

        fivetracks_tree->Branch("psiPrimeOne_p_pt",       &psiPrimeOne_p_pt,       "psiPrimeOne_p_pt/D");
        fivetracks_tree->Branch("psiPrimeOne_p_eta",      &psiPrimeOne_p_eta,      "psiPrimeOne_p_eta/D");
        fivetracks_tree->Branch("psiPrimeOne_p_phi",      &psiPrimeOne_p_phi,      "psiPrimeOne_p_phi/D");
        fivetracks_tree->Branch("psiPrimeOne_m_n",        &psiPrimeOne_m_n,        "psiPrimeOne_m_n/D");

        fivetracks_tree->Branch("psiPrimeOne_m_pt",       &psiPrimeOne_m_pt,       "psiPrimeOne_m_pt/D");
        fivetracks_tree->Branch("psiPrimeOne_m_eta",      &psiPrimeOne_m_eta,      "psiPrimeOne_m_eta/D");
        fivetracks_tree->Branch("psiPrimeOne_m_phi",      &psiPrimeOne_m_phi,      "psiPrimeOne_m_phi/D");
        fivetracks_tree->Branch("psiPrimeOne_m_n",        &psiPrimeOne_m_n,        "psiPrimeOne_m_n/D");

        //Two
        fivetracks_tree->Branch("psiPrimeTwo_pt",         &psiPrimeTwo_pt,        "psiPrimeTwo_pt/D");
        fivetracks_tree->Branch("psiPrimeTwo_eta",        &psiPrimeTwo_eta,       "psiPrimeTwo_eta/D");
        fivetracks_tree->Branch("psiPrimeTwo_phi",        &psiPrimeTwo_phi,       "psiPrimeTwo_phi/D");

        fivetracks_tree->Branch("psiPrimeTwo_p_pt",       &psiPrimeTwo_p_pt,       "psiPrimeTwo_p_pt/D");
        fivetracks_tree->Branch("psiPrimeTwo_p_eta",      &psiPrimeTwo_p_eta,      "psiPrimeTwo_p_eta/D");
        fivetracks_tree->Branch("psiPrimeTwo_p_phi",      &psiPrimeTwo_p_phi,      "psiPrimeTwo_p_phi/D");
        fivetracks_tree->Branch("psiPrimeTwo_m_n",        &psiPrimeTwo_m_n,        "psiPrimeTwo_m_n/D");

        fivetracks_tree->Branch("psiPrimeTwo_m_pt",       &psiPrimeTwo_m_pt,       "psiPrimeTwo_m_pt/D");
        fivetracks_tree->Branch("psiPrimeTwo_m_eta",      &psiPrimeTwo_m_eta,      "psiPrimeTwo_m_eta/D");
        fivetracks_tree->Branch("psiPrimeTwo_m_phi",      &psiPrimeTwo_m_phi,      "psiPrimeTwo_m_phi/D");
        fivetracks_tree->Branch("psiPrimeTwo_m_n",        &psiPrimeTwo_m_n,        "psiPrimeTwo_m_n/D");

        //TriTrak system
        fivetracks_tree->Branch("triTrak_pt",         &triTrak_pt,         "triTrak_pt/D");
        fivetracks_tree->Branch("triTrak_eta",        &triTrak_eta,        "triTrak_eta/D");
        fivetracks_tree->Branch("triTrak_phi",        &triTrak_phi,        "triTrak_phi/D");

        numMasses = 5;

        TLorentzVector zero;

        for(size_t i = 0; i<numMasses;i++)
        {
          fiveTracksMass.push_back(2.0);
          fiveTracksMassRef.push_back(2.0);

          fiveTracksVNDof.push_back(-1.0);
          fiveTracksVChi2.push_back(-1.0);
          fiveTracksVProb.push_back(-0.1);

          fiveTracksCTau.push_back(-1000.0);
          fiveTracksCTauErr.push_back(-1000.0);
          fiveTracksCosAlpha.push_back(-1.1);

          trackOneMass.push_back(0.0);
          trackTwoMass.push_back(0.0);
          trackThreeMass.push_back(0.0);

          five_p4.push_back(zero);
          five_ref_p4.push_back(zero);
          dimuonDiTrakOne_p4.push_back(zero);
          dimuonDiTrakTwo_p4.push_back(zero);
          dimuonDiTrakThree_p4.push_back(zero);
          psiPrimeOne_p4.push_back(zero);
          psiPrimeTwo_p4.push_back(zero);

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


          name = "onePsiPrime_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&psiPrimeOne[i],var.c_str());
          name = "twoPsiPrime_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&psiPrimeTwo[i],var.c_str());

          name = "onePsiPrime_p_m" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&psiPrimeOne_p_Mass[i],var.c_str());
          name = "onePsiPrime_n_m" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&psiPrimeOne_n_Mass[i],var.c_str());
          name = "twoPsiPrime_p_m" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&psiPrimeTwo_p_Mass[i],var.c_str());
          name = "twoPsiPrime_n_m" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&psiPrimeTwo_n_Mass[i],var.c_str());


          name = "onedimuonDiTrak_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&dimuonDiTrakOne[i],var.c_str());
          name = "twodimuonDiTrak_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&dimuonDiTrakTwo[i],var.c_str());
          name = "threedimuonDiTrak_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&dimuonDiTrakThree[i],var.c_str());


          name = "trackOne_m_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&trackOneMass[i],var.c_str());
          name = "trackTwo_m_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&trackTwoMass[i],var.c_str());
          name = "trackThree_m_" + refNames[i]; var = name + "/D";
          fivetracks_tree->Branch(name.c_str(),&trackThreeMass[i],var.c_str());

          name = "five_p4_" + refNames[i];
          fivetracks_tree->Branch(name.c_str(),"TLorentzVector",&five_p4[i]);
          name = "five_ref_p4_" + refNames[i];
          fivetracks_tree->Branch(name.c_str(),"TLorentzVector",&five_ref_p4[i]);

          name = "dimuonDiTrakOne_p4_" + refNames[i];
          fivetracks_tree->Branch(name.c_str(),"TLorentzVector",&dimuonDiTrakOne_p4[i]);
          name = "dimuonDiTrakTwo_p4_" + refNames[i];
          fivetracks_tree->Branch(name.c_str(),"TLorentzVector",&dimuonDiTrakTwo_p4[i]);
          name = "dimuonDiTrakThree_p4_" + refNames[i];
          fivetracks_tree->Branch(name.c_str(),"TLorentzVector",&dimuonDiTrakThree_p4[i]);

          name = "psiPrimeOne_p4_" + refNames[i];
          fivetracks_tree->Branch(name.c_str(),"TLorentzVector",&psiPrimeOne_p4[i]);
          name = "psiPrimeTwo_p4_" + refNames[i];
          fivetracks_tree->Branch(name.c_str(),"TLorentzVector",&psiPrimeTwo_p4[i]);


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

if(!OnlyGen_)
  if (!fivetracks_cand_handle.isValid()) std::cout<< "No five tracks information " << run << "," << event <<std::endl;

  if (fivetracks_cand_handle.isValid()) {

    pat::CompositeCandidate five_cand,
    const pat::PackedCandidate *trakOne_cand, *trakTwo_cand, *trakThree_cand;
    const pat::CompositeCandidate *dimuonDiTrakOne_cand, *dimuonDiTrakTwo_cand, *dimuonDiTrakThree_cand, *dimuonditrak_cand, *dimuon_cand;

    noFiveCandidates = (Int_t)(fivetracks_cand_handle.size());

    for (unsigned int i=0; i< fivetracks_cand_handle->size(); i++)
    {
      five_cand  = fivetracks_cand_handle->at(i);
      dimuonditrak_id = userInt("dimuontt_index");

      dimuonditrak_cand = dynamic_cast<const pat::CompositeCandidate*>(five_cand.daughter("dimuonditrak"));
      dimuon_cand = dynamic_cast<const pat::CompositeCandidate*>(dimuonditrak_cand->daughter("dimuon"));

      trakOne_cand = dynamic_cast<const pat::PackedCandidate*>(five_cand.daughter("trakOne"));
      trakTwo_cand = dynamic_cast<const pat::PackedCandidate*>(five_cand.daughter("trakTwo"));
      trakThree_cand = dynamic_cast<const pat::PackedCandidate*>(five_cand.daughter("trakThree"));

      dimuonDiTrakOne_cand = dynamic_cast<const pat::CompositeCandidate*>(five_cand.daughter("dimuonDiTrakOne"));
      dimuonDiTrakTwo_cand = dynamic_cast<const pat::CompositeCandidate*>(five_cand.daughter("dimuonDiTrakTwo"));
      dimuonDiTrakThree_cand = dynamic_cast<const pat::CompositeCandidate*>(five_cand.daughter("dimuonDiTrakThree"));

      dimuonditrak_p4.SetPtEtaPhiM(dimuonditrak_cand->pt(),dimuonditrak_cand->eta(),dimuonditrak_cand->phi(),dimuonditrak_cand->mass());
      dimuonditrk_m = dimuonditrak_cand->mass();
      dimuonditrk_pt = dimuonditrak_cand->pt();
      dimuonditrk_eta = dimuonditrak_cand->eta();
      dimuonditrk_phi = dimuonditrak_cand->phi();

      dimuon_p4.SetPtEtaPhiM(dimuon_cand->pt(),dimuon_cand->eta(),dimuon_cand->phi(),dimuon_cand->mass());
      dimuon_m = dimuon_cand->mass();
      dimuon_pt = dimuon_cand->pt();
      dimuon_eta = dimuon_cand->eta();
      dimuon_phi = dimuon_cand->phi();

      //trakOne_cand.SetPtEtaPhiM(dimuon_cand->pt(),dimuon_cand->eta(),dimuon_cand->phi(),dimuon_cand->mass());
      highTrack_pt = trakTwo_cand->pt();
      highTrack_eta = trakTwo_cand->eta();
      highTrack_phi = trakTwo_cand->phi();
      highTrack_charge = trakTwo_cand->charge();

      lowTrack_pt = trakOne_cand->pt();
      lowTrack_eta = trakOne_cand->eta();
      lowTrack_phi = trakOne_cand->phi();
      lowTrack_charge = trakOne_cand->charge();

      thirdTrack_pt = trakThree_cand->pt();
      thirdTrack_eta = trakThree_cand->eta();
      thirdTrack_phi = trakThree_cand->phi();
      thirdTrack_charge = trakThree_cand->charge();

      dimuonDiTrakOne_pt     = dimuonDiTrakOne_cand->pt();
      dimuonDiTrakOne_eta    = dimuonDiTrakOne_cand->eta();
      dimuonDiTrakOne_phi    = dimuonDiTrakOne_cand->phi();
      dimuonDiTrakOne_charge = dimuonDiTrakOne_cand->charge();

      dimuonDiTrakTwo_pt     = dimuonDiTrakTwo_cand->pt();
      dimuonDiTrakTwo_eta    = dimuonDiTrakTwo_cand->eta();
      dimuonDiTrakTwo_phi    = dimuonDiTrakTwo_cand->phi();
      dimuonDiTrakTwo_charge = dimuonDiTrakTwo_cand->charge();

      dimuonDiTrakThree_pt     = dimuonDiTrakThree_cand->pt();
      dimuonDiTrakThree_eta    = dimuonDiTrakThree_cand->eta();
      dimuonDiTrakThree_phi    = dimuonDiTrakThree_cand->phi();
      dimuonDiTrakThree_charge = dimuonDiTrakThree_cand->charge();


      for(size_t i = 0; i<numMasses;i++)
      {

        const pat::CompositeCandidate* five_cand_ref, five_cand_ref_ref;
        const pat::PackedCandidate *trakOne_cand_ref, *trakTwo_cand_ref, *trakThree_cand_ref;
        const pat::CompositeCandidate *dimuonDiTrakOne_cand_ref, *dimuonDiTrakTwo_cand_ref, *dimuonDiTrakThree_cand_ref;

        std::string name = "fiveCand_" + std::to_string(i);

        five_cand_ref     = dynamic_cast<const pat::CompositeCandidate*>(five_cand.daughter(name));
        five_cand_ref_ref = dynamic_cast<const pat::CompositeCandidate*>(five_cand_ref->daughter("fiveRef"));

        fiveTracksMass[i] = five_cand_ref->mass();

        fiveTracksMassRef[i] = five_cand.userFloat("mass_ref");

        fiveTracksVProb[i] = five_cand.userFloat("vProb");
        fiveTracksVNDof[i] = five_cand.userFloat("nDof");
        fiveTracksVChi2[i] = five_cand.userFloat("vChi2");

        fiveTracksCTau[i] = five_cand.userFloat("ctauPV");
        fiveTracksCTauErr[i] = five_cand.userFloat("ctauErrPV");
        fiveTracksCosAlpha[i] = five_cand.userFloat("cosAlpha");

        trakOne_cand_ref = dynamic_cast<const pat::PackedCandidate*>(five_cand_ref->daughter("trakOne"));
        trakTwo_cand_ref = dynamic_cast<const pat::PackedCandidate*>(five_cand_ref->daughter("trakTwo"));
        trakThree_cand_ref = dynamic_cast<const pat::PackedCandidate*>(five_cand_ref->daughter("trakThree"));

        dimuonDiTrakOne_cand_ref = dynamic_cast<const pat::CompositeCandidate*>(five_cand_ref->daughter("dimuonDiTrakOne"));
        dimuonDiTrakTwo_cand_ref = dynamic_cast<const pat::CompositeCandidate*>(five_cand_ref->daughter("dimuonDiTrakTwo"));
        dimuonDiTrakThree_cand_ref = dynamic_cast<const pat::CompositeCandidate*>(five_cand_ref->daughter("dimuonDiTrakThree"));

        dimuonDiTrakOne[i]  = dimuonDiTrakOne_cand_ref->mass();
        dimuonDiTrakTwo[i]  = dimuonDiTrakTwo_cand_ref->mass();
        dimuonDiTrakThree[i]  = dimuonDiTrakThree_cand_ref->mass();

        trackOneMass[i]  = trakOne_cand_ref.mass();
        trackTwoMass[i]  = trakTwo_cand_ref.mass();
        trackThreeMass[i]  = trakThree_cand_ref.mass();

        five_p4[i].SetPtEtaPhiM(five_cand_ref->pt(),five_cand_ref->eta(),five_cand_ref->phi(),five_cand_ref->mass());
        five_ref_p4[i].SetPtEtaPhiM(five_cand_ref_ref->pt(),five_cand_ref_ref->eta(),five_cand_ref_ref->phi(),five_cand_ref_ref->mass());

        dimuonDiTrakOne_p4[i].SetPtEtaPhiM(dimuonDiTrakOne_cand_ref->pt(),dimuonDiTrakOne_cand_ref->eta(),dimuonDiTrakOne_cand_ref->phi(),dimuonDiTrakOne_cand_ref->mass());
        dimuonDiTrakTwo_p4[i].SetPtEtaPhiM(dimuonDiTrakTwo_cand_ref->pt(),dimuonDiTrakTwo_cand_ref->eta(),dimuonDiTrakTwo_cand_ref->phi(),dimuonDiTrakTwo_cand_ref->mass());
        dimuonDiTrakThree_p4[i].SetPtEtaPhiM(dimuonDiTrakThree_cand_ref->pt(),dimuonDiTrakThree_cand_ref->eta(),dimuonDiTrakThree_cand_ref->phi(),dimuonDiTrakThree_cand_ref->mass());


        if(dimuonDiTrakOne_cand_ref->charge() != 0)
        {
          psiPrimeOne[i] = dimuonDiTrakTwo_cand_ref->mass();
          psiPrimeTwo[i] = dimuonDiTrakThree_cand_ref->mass();

          psiPrimeOne_p4[i].SetPtEtaPhiM(dimuonDiTrakTwo_cand_ref->pt(),dimuonDiTrakTwo_cand_ref->eta(),dimuonDiTrakTwo_cand_ref->phi(),dimuonDiTrakTwo_cand_ref->mass());
          psiPrimeTwo_p4[i].SetPtEtaPhiM(dimuonDiTrakThree_cand_ref->pt(),dimuonDiTrakThree_cand_ref->eta(),dimuonDiTrakThree_cand_ref->phi(),dimuonDiTrakThree_cand_ref->mass());

          if(dimuonDiTrakTwo_cand_ref->daughter("trak1")->charge()>0)
          {
            psiPrimeOne_p_Mass[i] = dimuonDiTrakTwo_cand_ref->daughter("trak1")->mass();
            psiPrimeOne_n_Mass[i] = dimuonDiTrakTwo_cand_ref->daughter("trak2")->mass();
          }
          else
          {
            psiPrimeOne_p_Mass[i] = dimuonDiTrakTwo_cand_ref->daughter("trak2")->mass();
            psiPrimeOne_n_Mass[i] = dimuonDiTrakTwo_cand_ref->daughter("trak1")->mass();
          }

          if(dimuonDiTrakThree_cand_ref->daughter("trak1")->charge()>0)
          {
            psiPrimeTwo_p_Mass[i] = dimuonDiTrakThree_cand_ref->daughter("trak1")->mass();
            psiPrimeTwo_n_Mass[i] = dimuonDiTrakThree_cand_ref->daughter("trak2")->mass();
          }
          else
          {
            psiPrimeTwo_p_Mass[i] = dimuonDiTrakThree_cand_ref->daughter("trak2")->mass();
            psiPrimeTwo_n_Mass[i] = dimuonDiTrakThree_cand_ref->daughter("trak1")->mass();
          }

        }

        if(dimuonDiTrakTwo_cand_ref->charge() != 0)
        {
          psiPrimeOne[i] = dimuonDiTrakOne_cand_ref->mass();
          psiPrimeTwo[i] = dimuonDiTrakThree_cand_ref->mass();

          psiPrimeOne_p4[i].SetPtEtaPhiM(dimuonDiTrakOne_cand_ref->pt(),dimuonDiTrakOne_cand_ref->eta(),dimuonDiTrakOne_cand_ref->phi(),dimuonDiTrakOne_cand_ref->mass());
          psiPrimeTwo_p4[i].SetPtEtaPhiM(dimuonDiTrakThree_cand_ref->pt(),dimuonDiTrakThree_cand_ref->eta(),dimuonDiTrakThree_cand_ref->phi(),dimuonDiTrakThree_cand_ref->mass());


          if(dimuonDiTrakOne_cand_ref->daughter("trak1")->charge()>0)
          {
            psiPrimeOne_p_Mass[i] = dimuonDiTrakOne_cand_ref->daughter("trak1")->mass();
            psiPrimeOne_n_Mass[i] = dimuonDiTrakOne_cand_ref->daughter("trak2")->mass();
          }
          else
          {
            psiPrimeOne_p_Mass[i] = dimuonDiTrakOne_cand_ref->daughter("trak2")->mass();
            psiPrimeOne_n_Mass[i] = dimuonDiTrakOne_cand_ref->daughter("trak1")->mass();
          }

          if(dimuonDiTrakThree_cand_ref->daughter("trak1")->charge()>0)
          {
            psiPrimeTwo_p_Mass[i] = dimuonDiTrakThree_cand_ref->daughter("trak1")->mass();
            psiPrimeTwo_n_Mass[i] = dimuonDiTrakThree_cand_ref->daughter("trak2")->mass();
          }
          else
          {
            psiPrimeTwo_p_Mass[i] = dimuonDiTrakThree_cand_ref->daughter("trak2")->mass();
            psiPrimeTwo_n_Mass[i] = dimuonDiTrakThree_cand_ref->daughter("trak1")->mass();
          }
        }
        if(dimuonDiTrakThree_cand_ref->charge() != 0)
        {
          psiPrimeOne[i] = dimuonDiTrakOne_cand_ref->mass();
          psiPrimeTwo[i] = dimuonDiTrakTwo_cand_ref->mass();

          psiPrimeOne_p4[i].SetPtEtaPhiM(dimuonDiTrakOne_cand_ref->pt(),dimuonDiTrakOne_cand_ref->eta(),dimuonDiTrakOne_cand_ref->phi(),dimuonDiTrakOne_cand_ref->mass());
          psiPrimeTwo_p4[i].SetPtEtaPhiM(dimuonDiTrakTwo_cand_ref->pt(),dimuonDiTrakTwo_cand_ref->eta(),dimuonDiTrakTwo_cand_ref->phi(),dimuonDiTrakTwo_cand_ref->mass());

          if(dimuonDiTrakOne_cand_ref->daughter("trak1")->charge()>0)
          {
            psiPrimeOne_p_Mass[i] = dimuonDiTrakOne_cand_ref->daughter("trak1")->mass();
            psiPrimeOne_n_Mass[i] = dimuonDiTrakOne_cand_ref->daughter("trak2")->mass();
          }
          else
          {
            psiPrimeOne_p_Mass[i] = dimuonDiTrakOne_cand_ref->daughter("trak2")->mass();
            psiPrimeOne_n_Mass[i] = dimuonDiTrakOne_cand_ref->daughter("trak1")->mass();
          }

          if(dimuonDiTrakTwo_cand_ref->daughter("trak1")->charge()>0)
          {
            psiPrimeTwo_p_Mass[i] = dimuonDiTrakTwo_cand_ref->daughter("trak1")->mass();
            psiPrimeTwo_n_Mass[i] = dimuonDiTrakTwo_cand_ref->daughter("trak2")->mass();
          }
          else
          {
            psiPrimeTwo_p_Mass[i] = dimuonDiTrakTwo_cand_ref->daughter("trak2")->mass();
            psiPrimeTwo_n_Mass[i] = dimuonDiTrakTwo_cand_ref->daughter("trak1")->mass();
          }
        }


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
