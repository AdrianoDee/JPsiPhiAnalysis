/*
   Package:    DoubleDiMuonRootupler
   Class:      DoubleDiMuonRootupler

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


class DoubleDiMuonRootupler : public edm::EDAnalyzer {
   public:
      explicit DoubleDiMuonRootupler(const edm::ParameterSet&);
      ~DoubleDiMuonRootupler() override;

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
  edm::EDGetTokenT<pat::CompositeCandidateCollection> doubledimuons_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> doubledimuon_rf_cand_Label;
  edm::EDGetTokenT<reco::BeamSpot> thebeamspot_;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  int  doubledimuon_pdgid_, jpsi_pdgid_, phi_pdgid_;
  bool isMC_,OnlyBest_,OnlyGen_,AddMC_;
  UInt_t MomPdgId_, JPsiPdgId_, PhiPdgId_;
  std::vector<std::string>                            HLTs_;
  std::vector<std::string>                            HLTFilters_;

  UInt_t run,event,numPrimaryVertices, trigger;

  TLorentzVector doubledimuon_p4, jpsi_p4, phi_p4, mHighPhi_p4;
  TLorentzVector mLowJPsi_p4, mHighJPsi_p4, mLowPhi_p4;

  TLorentzVector doubledimuon_rf_p4,doubledimuon_rf_const_p4, jpsi_rf_p4, phi_rf_p4;
  TLorentzVector mHighPhi_rf_p4, mLowJPsi_rf_p4, mHighJPsi_rf_p4, mLowPhi_rf_p4;

  Double_t jpsi_triggerMatch, phi_triggerMatch;
  Double_t mHighJPsiMatch, mLowJPsiMatch, mHighPhiMatch, mLowPhiMatch;

  Double_t noXCandidates;

  //////////////////////////
  // Four muons variables

  Double_t doubledimuon_charge;

  //Kin variables
  Double_t doubledimuon_m,doubledimuon_m_rf,doubledimuon_m_rf_c,doubledimuon_m_rf_d_c;
  Double_t doubledimuon_p, doubledimuon_pt, doubledimuon_eta, doubledimuon_theta, doubledimuon_y, doubledimuon_phi;
  Double_t doubledimuon_e, doubledimuon_dxy,doubledimuon_dxyErr, doubledimuon_dz,doubledimuon_dzErr;

  //Vertexing variables
  Double_t doubledimuon_vProb,  doubledimuon_vChi2,doubledimuon_nDof;
  Double_t doubledimuon_cosAlpha, doubledimuon_ctauPV, doubledimuon_ctauErrPV, doubledimuon_cosAlpha3D;
  Double_t doubledimuon_lxy, doubledimuon_lxyErr, doubledimuon_lxyz, doubledimuon_lxyzErr;

  Double_t doubledimuon_rf_lxy,doubledimuon_rf_lxyErr,doubledimuon_rf_lxyz,doubledimuon_rf_lxyzErr;
  Double_t doubledimuon_rf_vProb,  doubledimuon_rf_vChi2, doubledimuon_rf_nDof;
  Double_t doubledimuon_rf_cosAlpha, doubledimuon_rf_ctauPV, doubledimuon_rf_ctauErrPV;

  Double_t doubledimuon_rf_c_vProb, doubledimuon_rf_c_vChi2, doubledimuon_rf_c_nDof;
  Double_t doubledimuon_rf_c_cosAlpha, doubledimuon_rf_c_ctauPV, doubledimuon_rf_c_ctauErrPV;
  Double_t doubledimuon_rf_c_lxy, doubledimuon_rf_c_lxyErr, doubledimuon_rf_c_lxyz, doubledimuon_rf_c_lxyzErr;

  Double_t doubledimuon_vx, doubledimuon_vy, doubledimuon_vz;
  Double_t pv_x, pv_y, pv_z;

  Double_t doubledimuon_dca_m1m2, doubledimuon_dca_m1t1, doubledimuon_dca_m1t2;
  Double_t doubledimuon_dca_t1t2,doubledimuon_dca_m2t1, doubledimuon_dca_m2t2;

  Double_t doubledimuon_cosAlphaDZ, doubledimuon_cosAlpha3DDZ, doubledimuon_ctauPVDZ, doubledimuon_ctauErrPVDZ;
  Double_t doubledimuon_cosAlphaBS, doubledimuon_cosAlpha3DBS, doubledimuon_ctauPVBS, doubledimuon_ctauErrPVBS;
  Double_t doubledimuon_lxyDZ, doubledimuon_lxyErrDZ, doubledimuon_lxyzDZ, doubledimuon_lxyzErrDZ;
  Double_t doubledimuon_lxyBS, doubledimuon_lxyErrBS, doubledimuon_lxyzBS, doubledimuon_lxyzErrBS;

  //DCAs

  Double_t doubledimuon_dca_mp1mp2, doubledimuon_dca_mp1mj1, doubledimuon_dca_mp1mj2;
  Double_t doubledimuon_dca_mp2mj1, doubledimuon_dca_mp2mj2, doubledimuon_dca_mj1mj2;
  //////////////////////////
  // Single muons variables

  //Kin variables
  Double_t mHighPhi_p,mLowPhi_p,mHighJPsi_p,mLowJPsi_p;
  Double_t mHighPhi_pt,mLowPhi_pt,mHighJPsi_pt,mLowJPsi_pt;
  Double_t mHighPhi_ptErr,mLowPhi_ptErr,mHighJPsi_ptErr,mLowJPsi_ptErr;

  //Angular variables - Eta;Phi,Theata,Lambda
  Double_t mHighJPsi_eta, mLowJPsi_eta, mHighPhi_eta, mLowPhi_eta;
  Double_t mHighJPsi_etaErr, mLowJPsi_etaErr, mHighPhi_etaErr, mLowPhi_etaErr;

  Double_t mHighJPsi_phi, mLowJPsi_phi, mHighPhi_phi, mLowPhi_phi;
  Double_t mHighJPsi_phiErr, mLowJPsi_phiErr, mHighPhi_phiErr, mLowPhi_phiErr;

  Double_t mHighJPsi_theta, mLowJPsi_theta, mHighPhi_theta, mLowPhi_theta;
  Double_t mHighJPsi_thetaErr, mLowJPsi_thetaErr, mHighPhi_thetaErr, mLowPhi_thetaErr;

  Double_t mHighJPsi_lambda, mLowJPsi_lambda, mHighPhi_lambda, mLowPhi_lambda;
  Double_t mHighJPsi_lambdaErr, mLowJPsi_lambdaErr, mHighPhi_lambdaErr, mLowPhi_lambdaErr;

  //Impact variables

  Double_t mLowPhi_dxy, mLowPhi_dxyErr, mLowPhi_dz, mLowPhi_dzErr;
  Double_t mHighPhi_dxy, mHighPhi_dxyErr, mHighPhi_dz, mHighPhi_dzErr;

  Double_t mHighJPsi_dxy, mHighJPsi_dxyErr, mHighJPsi_dz, mHighJPsi_dzErr;
  Double_t mLowJPsi_dxy, mLowJPsi_dxyErr, mLowJPsi_dz, mLowJPsi_dzErr;

  //Tracker hits, layers,
  Double_t mHighJPsi_NPixelHits, mHighJPsi_NStripHits, mHighJPsi_NTrackhits, mHighJPsi_NBPixHits;
  Double_t mHighJPsi_NPixLayers, mHighJPsi_NTraLayers, mHighJPsi_NStrLayers, mHighJPsi_NBPixLayers;

  Double_t mLowJPsi_NPixelHits, mLowJPsi_NStripHits, mLowJPsi_NTrackhits, mLowJPsi_NBPixHits;
  Double_t mLowJPsi_NPixLayers, mLowJPsi_NTraLayers, mLowJPsi_NStrLayers, mLowJPsi_NBPixLayers;

  Double_t mHighPhi_NPixelHits, mHighPhi_NStripHits, mHighPhi_NTrackhits, mHighPhi_NBPixHits;
  Double_t mHighPhi_NPixLayers, mHighPhi_NTraLayers, mHighPhi_NStrLayers, mHighPhi_NBPixLayers;

  Double_t mLowPhi_NPixelHits, mLowPhi_NStripHits, mLowPhi_NTrackhits, mLowPhi_NBPixHits;
  Double_t mLowPhi_NPixLayers, mLowPhi_NTraLayers, mLowPhi_NStrLayers, mLowPhi_NBPixLayers;

  //Muon Flags
  Double_t mHighJPsi_isLoose, mHighJPsi_isSoft, mHighJPsi_isMedium, mHighJPsi_isHighPt;
  Double_t mLowJPsi_isLoose, mLowJPsi_isSoft, mLowJPsi_isMedium, mLowJPsi_isHighPt;

  Double_t mHighPhi_isLoose, mHighPhi_isSoft, mHighPhi_isMedium, mHighPhi_isHighPt;
  Double_t mLowPhi_isLoose, mLowPhi_isSoft, mLowPhi_isMedium, mLowPhi_isHighPt;

  Double_t mHighJPsi_isTracker, mHighJPsi_isGlobal, mLowJPsi_isTracker, mLowJPsi_isGlobal;
  Double_t mHighPhi_isTracker, mHighPhi_isGlobal, mLowPhi_isTracker, mLowPhi_isGlobal;

  Double_t mHighJPsi_type, mLowJPsi_type, mHighPhi_type, mLowPhi_type;

  ///////////////////
  //DiMuons Variables

  Double_t jpsi_triggerMatch_rf, phi_triggerMatch_rf;

  //Kin variables
  Double_t jpsi_m,jpsi_m_rf,jpsi_m_rf_c,jpsi_m_rf_d_c;
  Double_t jpsi_p, jpsi_theta, jpsi_eta, jpsi_pt, jpsi_y;
  Double_t jpsi_e, jpsi_dxy,jpsi_dxyErr, jpsi_dz,jpsi_dzErr;
  Double_t jpsi_vProb, jpsi_vChi2, jpsi_DCA, jpsi_ctauPV, jpsi_ctauErrPV, jpsi_cosAlpha;
  Double_t jpsi_lxy,jpsi_lxyz,jpsi_lxyErr,jpsi_lxyzErr,jpsi_cosAlpha3D;

  Double_t phi_m,phi_m_rf,phi_m_rf_c,phi_m_rf_d_c;
  Double_t phi_p, phi_pt, phi_eta, phi_theta, phi_y;
  Double_t phi_e, phi_dxy,phi_dxyErr, phi_dz,phi_dzErr;
  Double_t phi_vProb, phi_vChi2, phi_DCA, phi_ctauPV, phi_ctauErrPV, phi_cosAlpha;
  Double_t phi_lxy,phi_lxyz,phi_lxyErr,phi_lxyzErr,phi_cosAlpha3D;


  //////////////
  //MC variables

  Double_t gen_doubledimuon_pdg, gen_phi_pdg, gen_jpsi_pdg;
  Double_t gen_doubledimuon_prompt, gen_phi_prompt, gen_jpsi_prompt;
  Double_t gen_doubledimuon_pt, gen_doubledimuon_p, gen_doubledimuon_eta;
  Double_t gen_phi_pt, gen_phi_p, gen_phi_eta;
  Double_t gen_jpsi_pt, gen_jpsi_p, gen_jpsi_eta;
  Double_t isBestCandidate;

  TLorentzVector gen_doubledimuon_p4, gen_jpsi_p4, gen_phi_p4;
  TLorentzVector gen_mHighPhi_p4, gen_mLowJPsi_p4, gen_mHighJPsi_p4, gen_mLowPhi_p4;

  TTree* fourmuon_tree, *fourmuon_tree_rf;
  edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;
};

//
// constants, enums and typedefs
//

UInt_t DoubleDiMuonRootupler::isTriggerMatched(pat::CompositeCandidate *diMuon_cand) {
  const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("highMuon"));
  const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("lowMuon"));
  UInt_t matched = 0;  // if no list is given, is not matched

  // if matched a given trigger, set the bit, in the same order as listed
  for (unsigned int iTr = 0; iTr<HLTFilters_.size(); iTr++ ) {
    // //std::cout << HLTFilters_[iTr] << std::endl;
    const pat::TriggerObjectStandAloneCollection mu1HLTMatches = muon1->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
    const pat::TriggerObjectStandAloneCollection mu2HLTMatches = muon2->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
    if (!mu1HLTMatches.empty() && !mu2HLTMatches.empty()) matched += (1<<iTr);
    // if (!mu1HLTMatches.empty() && !mu2HLTMatches.empty()) //std::cout << std::endl << HLTFilters_[iTr] << std::endl;
  }

  return matched;
}


// constructors and destructor
//
DoubleDiMuonRootupler::DoubleDiMuonRootupler(const edm::ParameterSet& iConfig):
        doubledimuons_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("FourMuons"))),
        // thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpot"))),
        primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PrimaryVertices"))),
        triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
	      //isMC_(iConfig.getParameter<bool>("isMC")),
        OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
        OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
        AddMC_(iConfig.getParameter<bool>("AddMC")),
        HLTs_(iConfig.getParameter<std::vector<std::string>>("HLTs")),
        HLTFilters_(iConfig.getParameter<std::vector<std::string>>("Filters"))
{
	      edm::Service<TFileService> fs;
        fourmuon_tree = fs->make<TTree>("FourMuonTree","Tree of JPsi and Phi in 4 Muons");

        fourmuon_tree->Branch("run",                &run,                "run/D");
        fourmuon_tree->Branch("event",              &event,              "event/D");
        fourmuon_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/D");
        fourmuon_tree->Branch("trigger",            &trigger,            "trigger/D");

        fourmuon_tree->Branch("noXCandidates",      &noXCandidates,      "noXCandidates/D");

        fourmuon_tree->Branch("doubledimuon_p4",  "TLorentzVector", &doubledimuon_p4);
        fourmuon_tree->Branch("phi_p4",           "TLorentzVector", &phi_p4);
        fourmuon_tree->Branch("jpsi_p4",          "TLorentzVector", &jpsi_p4);
        fourmuon_tree->Branch("mLowPhi_p4",       "TLorentzVector", &mLowPhi_p4);
        fourmuon_tree->Branch("mLowJPsi_p4",      "TLorentzVector", &mLowJPsi_p4);
        fourmuon_tree->Branch("mHighJPsi_p4",     "TLorentzVector", &mHighJPsi_p4);
        fourmuon_tree->Branch("mLowPhi_p4",       "TLorentzVector", &mLowPhi_p4);

        fourmuon_tree->Branch("doubledimuon_rf_p4", "TLorentzVector", &doubledimuon_rf_p4);
        fourmuon_tree->Branch("phi_rf_p4",          "TLorentzVector", &phi_rf_p4);
        fourmuon_tree->Branch("jpsi_rf_p4",         "TLorentzVector", &jpsi_rf_p4);
        fourmuon_tree->Branch("mLowPhi_rf_p4",      "TLorentzVector", &mLowPhi_rf_p4);
        fourmuon_tree->Branch("mLowJPsi_rf_p4",     "TLorentzVector", &mLowJPsi_rf_p4);
        fourmuon_tree->Branch("mHighJPsi_rf_p4",    "TLorentzVector", &mHighJPsi_rf_p4);
        fourmuon_tree->Branch("mLowPhi_rf_p4",      "TLorentzVector", &mLowPhi_rf_p4);

        //Trigger Matching
        fourmuon_tree->Branch("jpsi_triggerMatch",  &jpsi_triggerMatch, "jpsi_triggerMatch/D");
        fourmuon_tree->Branch("phi_triggerMatch",   &phi_triggerMatch,  "phi_triggerMatch/D");
        fourmuon_tree->Branch("mHighJPsiMatch",     &mHighJPsiMatch,    "mHighJPsiMatch/D");
        fourmuon_tree->Branch("mLowJPsiMatch",      &mLowJPsiMatch,     "mLowJPsiMatch/D");
        fourmuon_tree->Branch("mHighPhiMatch",      &mHighPhiMatch,     "mHighPhiMatch/D");
        fourmuon_tree->Branch("mLowPhiMatch",       &mLowPhiMatch,      "mLowPhiMatch/D");

        //Four Muons Variables
        fourmuon_tree->Branch("doubledimuon_charge",&doubledimuon_charge,"doubledimuon_charge/D");
        //Kin
        fourmuon_tree->Branch("doubledimuon_m",&doubledimuon_m,"doubledimuon_m/D");
        fourmuon_tree->Branch("doubledimuon_m_rf",&doubledimuon_m_rf,"doubledimuon_m_rf/D");
        fourmuon_tree->Branch("doubledimuon_m_rf_c",&doubledimuon_m_rf_c,"doubledimuon_m_rf_c/D");
        fourmuon_tree->Branch("doubledimuon_m_rf_d_c",&doubledimuon_m_rf_d_c,"doubledimuon_m_rf_d_c/D");

        fourmuon_tree->Branch("doubledimuon_p",&doubledimuon_p,"doubledimuon_p/D");
        fourmuon_tree->Branch("doubledimuon_pt",&doubledimuon_pt,"doubledimuon_pt/D");
        fourmuon_tree->Branch("doubledimuon_e",&doubledimuon_e,"doubledimuon_e/D");

        //Angular
        fourmuon_tree->Branch("doubledimuon_eta",&doubledimuon_eta,"doubledimuon_eta/D");
        fourmuon_tree->Branch("doubledimuon_theta",&doubledimuon_theta,"doubledimuon_theta/D");
        fourmuon_tree->Branch("doubledimuon_y",&doubledimuon_y,"doubledimuon_y/D");

        //Impact
        fourmuon_tree->Branch("doubledimuon_dxy",&doubledimuon_dxy,"doubledimuon_dxy/D");
        fourmuon_tree->Branch("doubledimuon_dxyErr",&doubledimuon_dxyErr,"doubledimuon_dxyErr/D");
        fourmuon_tree->Branch("doubledimuon_dz",&doubledimuon_dz,"doubledimuon_dz/D");
        fourmuon_tree->Branch("doubledimuon_dzErr",&doubledimuon_dzErr,"doubledimuon_dzErr/D");

        //Vertex
        fourmuon_tree->Branch("doubledimuon_vProb",&doubledimuon_vProb,"doubledimuon_vProb/D");
        fourmuon_tree->Branch("doubledimuon_vChi2",&doubledimuon_vChi2,"doubledimuon_vChi2/D");
        fourmuon_tree->Branch("doubledimuon_nDof",&doubledimuon_nDof,"doubledimuon_nDof/D");

        fourmuon_tree->Branch("doubledimuon_rf_vProb",&doubledimuon_rf_vProb,"doubledimuon_rf_vProb/D");
        fourmuon_tree->Branch(" doubledimuon_rf_vChi2",&doubledimuon_rf_vChi2," doubledimuon_rf_vChi2/D");
        fourmuon_tree->Branch("doubledimuon_rf_nDof",&doubledimuon_rf_nDof,"doubledimuon_rf_nDof/D");

        fourmuon_tree->Branch("doubledimuon_rf_c_vProb",&doubledimuon_rf_c_vProb,"doubledimuon_rf_c_vProb/D");
        fourmuon_tree->Branch("doubledimuon_rf_c_vChi2",&doubledimuon_rf_c_vChi2,"doubledimuon_rf_c_vChi2/D");
        fourmuon_tree->Branch("doubledimuon_rf_c_nDof",&doubledimuon_rf_c_nDof,"doubledimuon_rf_c_nDof/D");

        fourmuon_tree->Branch("doubledimuon_vx",&doubledimuon_vx,"doubledimuon_vx/D");
        fourmuon_tree->Branch("doubledimuon_vy",&doubledimuon_vy,"doubledimuon_vy/D");
        fourmuon_tree->Branch("doubledimuon_vz",&doubledimuon_vz,"doubledimuon_vz/D");
        fourmuon_tree->Branch("pv_x",&pv_x,"pv_x/D");
        fourmuon_tree->Branch("pv_y",&pv_y,"pv_y/D");
        fourmuon_tree->Branch("pv_z",&pv_z,"pv_z/D");

        //Pointing Angle & Flight
        fourmuon_tree->Branch("doubledimuon_cosAlpha",&doubledimuon_cosAlpha,"doubledimuon_cosAlpha/D");
        fourmuon_tree->Branch("doubledimuon_cosAlpha3D",&doubledimuon_cosAlpha3D,"doubledimuon_cosAlpha3D/D");
        fourmuon_tree->Branch("doubledimuon_ctauPV",&doubledimuon_ctauPV,"doubledimuon_ctauPV/D");
        fourmuon_tree->Branch("doubledimuon_ctauErrPV",&doubledimuon_ctauErrPV,"doubledimuon_ctauErrPV/D");
        fourmuon_tree->Branch("doubledimuon_lxy",&doubledimuon_lxy,"doubledimuon_lxy/D");
        fourmuon_tree->Branch("doubledimuon_lxyErr",&doubledimuon_lxyErr,"doubledimuon_lxyErr/D");
        fourmuon_tree->Branch("doubledimuon_lxyz",&doubledimuon_lxyz,"doubledimuon_lxyz/D");
        fourmuon_tree->Branch("doubledimuon_lxyzErr",&doubledimuon_lxyzErr,"doubledimuon_lxyzErr/D");

        fourmuon_tree->Branch("doubledimuon_rf_cosAlpha",&doubledimuon_rf_cosAlpha,"doubledimuon_rf_cosAlpha/D");
        fourmuon_tree->Branch("doubledimuon_rf_ctauPV",&doubledimuon_rf_ctauPV,"doubledimuon_rf_ctauPV/D");
        fourmuon_tree->Branch("doubledimuon_rf_ctauErrPV",&doubledimuon_rf_ctauErrPV,"doubledimuon_rf_ctauErrPV/D");
        fourmuon_tree->Branch("doubledimuon_rf_lxy",&doubledimuon_rf_lxy,"doubledimuon_rf_lxy/D");
        fourmuon_tree->Branch("doubledimuon_rf_lxyErr",&doubledimuon_rf_lxyErr,"doubledimuon_rf_lxyErr/D");
        fourmuon_tree->Branch("doubledimuon_rf_lxyz",&doubledimuon_rf_lxyz,"doubledimuon_rf_lxyz/D");
        fourmuon_tree->Branch("doubledimuon_rf_lxyzErr",&doubledimuon_rf_lxyzErr,"doubledimuon_rf_lxyzErr/D");

        fourmuon_tree->Branch("doubledimuon_rf_c_cosAlpha",&doubledimuon_rf_c_cosAlpha,"doubledimuon_rf_c_cosAlpha/D");
        fourmuon_tree->Branch("doubledimuon_rf_c_ctauPV",&doubledimuon_rf_c_ctauPV,"doubledimuon_rf_c_ctauPV/D");
        fourmuon_tree->Branch("doubledimuon_rf_c_ctauErrPV",&doubledimuon_rf_c_ctauErrPV,"doubledimuon_rf_c_ctauErrPV/D");
        fourmuon_tree->Branch("doubledimuon_rf_c_lxy",&doubledimuon_rf_c_lxy,"doubledimuon_rf_c_lxy/D");
        fourmuon_tree->Branch("doubledimuon_rf_c_lxyErr",&doubledimuon_rf_c_lxyErr,"doubledimuon_rf_c_lxyErr/D");
        fourmuon_tree->Branch("doubledimuon_rf_c_lxyz",&doubledimuon_rf_c_lxyz,"doubledimuon_rf_c_lxyz/D");
        fourmuon_tree->Branch("doubledimuon_rf_c_lxyzErr",&doubledimuon_rf_c_lxyzErr,"doubledimuon_rf_c_lxyzErr/D");

        //DCA
        fourmuon_tree->Branch("doubledimuon_dca_mp1mp2",&doubledimuon_dca_mp1mp2,"doubledimuon_dca_mp1mp2/D");
        fourmuon_tree->Branch("doubledimuon_dca_mp1mj1",&doubledimuon_dca_mp1mj1,"doubledimuon_dca_mp1mj1/D");
        fourmuon_tree->Branch("doubledimuon_dca_mp1mj2",&doubledimuon_dca_mp1mj2,"doubledimuon_dca_mp1mj2/D");
        fourmuon_tree->Branch("doubledimuon_dca_mp2mj1",&doubledimuon_dca_mp2mj1,"doubledimuon_dca_mp2mj1/D");
        fourmuon_tree->Branch("doubledimuon_dca_mp2mj2",&doubledimuon_dca_mp2mj2,"doubledimuon_dca_mp2mj2/D");
        fourmuon_tree->Branch("doubledimuon_dca_mj1mj2",&doubledimuon_dca_mj1mj2,"doubledimuon_dca_mj1mj2/D");

        //More Vertex
        fourmuon_tree->Branch("doubledimuon_cosAlphaDZ",&doubledimuon_cosAlphaDZ,"doubledimuon_cosAlphaDZ/D");
        fourmuon_tree->Branch("doubledimuon_cosAlpha3DDZ",&doubledimuon_cosAlpha3DDZ,"doubledimuon_cosAlpha3DDZ/D");
        fourmuon_tree->Branch("doubledimuon_ctauPVDZ",&doubledimuon_ctauPVDZ,"doubledimuon_ctauPVDZ/D");
        fourmuon_tree->Branch("doubledimuon_ctauErrPVDZ",&doubledimuon_ctauErrPVDZ,"doubledimuon_ctauErrPVDZ/D");
        fourmuon_tree->Branch("doubledimuon_lxyDZ",&doubledimuon_lxyDZ,"doubledimuon_lxyDZ/D");
        fourmuon_tree->Branch("doubledimuon_lxyErrDZ",&doubledimuon_lxyErrDZ,"doubledimuon_lxyErrDZ/D");
        fourmuon_tree->Branch("doubledimuon_lxyzDZ",&doubledimuon_lxyzDZ,"doubledimuon_lxyzDZ/D");
        fourmuon_tree->Branch("doubledimuon_lxyzErrDZ",&doubledimuon_lxyzErrDZ,"doubledimuon_lxyzErrDZ/D");

        fourmuon_tree->Branch("doubledimuon_cosAlphaBS",&doubledimuon_cosAlphaBS,"doubledimuon_cosAlphaBS/D");
        fourmuon_tree->Branch("doubledimuon_cosAlpha3DBS",&doubledimuon_cosAlpha3DBS,"doubledimuon_cosAlpha3DBS/D");
        fourmuon_tree->Branch("doubledimuon_ctauPVBS",&doubledimuon_ctauPVBS,"doubledimuon_ctauPVBS/D");
        fourmuon_tree->Branch("doubledimuon_ctauErrPVBS",&doubledimuon_ctauErrPVBS,"doubledimuon_ctauErrPVBS/D");
        fourmuon_tree->Branch("doubledimuon_lxyBS",&doubledimuon_lxyBS,"doubledimuon_lxyBS/D");
        fourmuon_tree->Branch("doubledimuon_lxyErrBS",&doubledimuon_lxyErrBS,"doubledimuon_lxyErrBS/D");
        fourmuon_tree->Branch("doubledimuon_lxyzBS",&doubledimuon_lxyzBS,"doubledimuon_lxyzBS/D");
        fourmuon_tree->Branch("doubledimuon_lxyzErrBS",&doubledimuon_lxyzErrBS,"doubledimuon_lxyzErrBS/D");

        //Single muons variable
        //Kin
        fourmuon_tree->Branch("mHighPhi_p",&mHighPhi_p,"mHighPhi_p/D");
        fourmuon_tree->Branch("mLowPhi_p",&mLowPhi_p,"mLowPhi_p/D");
        fourmuon_tree->Branch("mHighJPsi_p",&mHighJPsi_p,"mHighJPsi_p/D");
        fourmuon_tree->Branch("mLowJPsi_p",&mLowJPsi_p,"mLowJPsi_p/D");
        fourmuon_tree->Branch("mHighPhi_pt",&mHighPhi_pt,"mHighPhi_pt/D");
        fourmuon_tree->Branch("mLowPhi_pt",&mLowPhi_pt,"mLowPhi_pt/D");
        fourmuon_tree->Branch("mHighJPsi_pt",&mHighJPsi_pt,"mHighJPsi_pt/D");
        fourmuon_tree->Branch("mLowJPsi_pt",&mLowJPsi_pt,"mLowJPsi_pt/D");
        fourmuon_tree->Branch("mHighPhi_ptErr",&mHighPhi_ptErr,"mHighPhi_ptErr/D");
        fourmuon_tree->Branch("mLowPhi_ptErr",&mLowPhi_ptErr,"mLowPhi_ptErr/D");
        fourmuon_tree->Branch("mHighJPsi_ptErr",&mHighJPsi_ptErr,"mHighJPsi_ptErr/D");
        fourmuon_tree->Branch("mLowJPsi_ptErr",&mLowJPsi_ptErr,"mLowJPsi_ptErr/D");

        //Angular
        fourmuon_tree->Branch("mHighJPsi_eta",&mHighJPsi_eta,"mHighJPsi_eta/D");
        fourmuon_tree->Branch("mLowJPsi_eta",&mLowJPsi_eta,"mLowJPsi_eta/D");
        fourmuon_tree->Branch("mHighPhi_eta",&mHighPhi_eta,"mHighPhi_eta/D");
        fourmuon_tree->Branch("mLowPhi_eta",&mLowPhi_eta,"mLowPhi_eta/D");

        fourmuon_tree->Branch("mHighJPsi_etaErr",&mHighJPsi_etaErr,"mHighJPsi_etaErr/D");
        fourmuon_tree->Branch("mLowJPsi_etaErr",&mLowJPsi_etaErr,"mLowJPsi_etaErr/D");
        fourmuon_tree->Branch("mHighPhi_etaErr",&mHighPhi_etaErr,"mHighPhi_etaErr/D");
        fourmuon_tree->Branch("mLowPhi_etaErr",&mLowPhi_etaErr,"mLowPhi_etaErr/D");

        fourmuon_tree->Branch("mHighJPsi_phi",&mHighJPsi_phi,"mHighJPsi_phi/D");
        fourmuon_tree->Branch("mLowJPsi_phi",&mLowJPsi_phi,"mLowJPsi_phi/D");
        fourmuon_tree->Branch("mHighPhi_phi",&mHighPhi_phi,"mHighPhi_phi/D");
        fourmuon_tree->Branch("mLowPhi_phi",&mLowPhi_phi,"mLowPhi_phi/D");

        fourmuon_tree->Branch("mHighJPsi_phiErr",&mHighJPsi_phiErr,"mHighJPsi_phiErr/D");
        fourmuon_tree->Branch("mLowJPsi_phiErr",&mLowJPsi_phiErr,"mLowJPsi_phiErr/D");
        fourmuon_tree->Branch("mHighPhi_phiErr",&mHighPhi_phiErr,"mHighPhi_phiErr/D");
        fourmuon_tree->Branch("mLowPhi_phiErr",&mLowPhi_phiErr,"mLowPhi_phiErr/D");

        fourmuon_tree->Branch("mHighJPsi_theta",&mHighJPsi_theta,"mHighJPsi_theta/D");
        fourmuon_tree->Branch("mLowJPsi_theta",&mLowJPsi_theta,"mLowJPsi_theta/D");
        fourmuon_tree->Branch("mHighPhi_theta",&mHighPhi_theta,"mHighPhi_theta/D");
        fourmuon_tree->Branch("mLowPhi_theta",&mLowPhi_theta,"mLowPhi_theta/D");

        fourmuon_tree->Branch("mHighJPsi_thetaErr",&mHighJPsi_thetaErr,"mHighJPsi_thetaErr/D");
        fourmuon_tree->Branch("mLowJPsi_thetaErr",&mLowJPsi_thetaErr,"mLowJPsi_thetaErr/D");
        fourmuon_tree->Branch("mHighPhi_thetaErr",&mHighPhi_thetaErr,"mHighPhi_thetaErr/D");
        fourmuon_tree->Branch("mLowPhi_thetaErr",&mLowPhi_thetaErr,"mLowPhi_thetaErr/D");

        fourmuon_tree->Branch("mHighJPsi_lambda",&mHighJPsi_lambda,"mHighJPsi_lambda/D");
        fourmuon_tree->Branch("mLowJPsi_lambda",&mLowJPsi_lambda,"mLowJPsi_lambda/D");
        fourmuon_tree->Branch("mHighPhi_lambda",&mHighPhi_lambda,"mHighPhi_lambda/D");
        fourmuon_tree->Branch("mLowPhi_lambda",&mLowPhi_lambda,"mLowPhi_lambda/D");

        fourmuon_tree->Branch("mHighJPsi_lambdaErr",&mHighJPsi_lambdaErr,"mHighJPsi_lambdaErr/D");
        fourmuon_tree->Branch("mLowJPsi_lambdaErr",&mLowJPsi_lambdaErr,"mLowJPsi_lambdaErr/D");
        fourmuon_tree->Branch("mHighPhi_lambdaErr",&mHighPhi_lambdaErr,"mHighPhi_lambdaErr/D");
        fourmuon_tree->Branch("mLowPhi_lambdaErr",&mLowPhi_lambdaErr,"mLowPhi_lambdaErr/D");

        //Impact
        fourmuon_tree->Branch("mLowPhi_dxy",&mLowPhi_dxy,"mLowPhi_dxy/D");
        fourmuon_tree->Branch("mLowPhi_dxyErr",&mLowPhi_dxyErr,"mLowPhi_dxyErr/D");
        fourmuon_tree->Branch("mLowPhi_dz",&mLowPhi_dz,"mLowPhi_dz/D");
        fourmuon_tree->Branch("mLowPhi_dzErr",&mLowPhi_dzErr,"mLowPhi_dzErr/D");

        fourmuon_tree->Branch("mHighPhi_dxy",&mHighPhi_dxy,"mHighPhi_dxy/D");
        fourmuon_tree->Branch("mHighPhi_dxyErr",&mHighPhi_dxyErr,"mHighPhi_dxyErr/D");
        fourmuon_tree->Branch("mHighPhi_dz",&mHighPhi_dz,"mHighPhi_dz/D");
        fourmuon_tree->Branch("mHighPhi_dzErr",&mHighPhi_dzErr,"mHighPhi_dzErr/D");

        fourmuon_tree->Branch("mHighJPsi_dxy",&mHighJPsi_dxy,"mHighJPsi_dxy/D");
        fourmuon_tree->Branch("mHighJPsi_dxyErr",&mHighJPsi_dxyErr,"mHighJPsi_dxyErr/D");
        fourmuon_tree->Branch("mHighJPsi_dz",&mHighJPsi_dz,"mHighJPsi_dz/D");
        fourmuon_tree->Branch("mHighJPsi_dzErr",&mHighJPsi_dzErr,"mHighJPsi_dzErr/D");

        fourmuon_tree->Branch("mLowJPsi_dxy",&mLowJPsi_dxy,"mLowJPsi_dxy/D");
        fourmuon_tree->Branch("mLowJPsi_dxyErr",&mLowJPsi_dxyErr,"mLowJPsi_dxyErr/D");
        fourmuon_tree->Branch("mLowJPsi_dz",&mLowJPsi_dz,"mLowJPsi_dz/D");
        fourmuon_tree->Branch("mLowJPsi_dzErr",&mLowJPsi_dzErr,"mLowJPsi_dzErr/D");

        //Tracker
        fourmuon_tree->Branch("mHighJPsi_NPixelHits",&mHighJPsi_NPixelHits,"mHighJPsi_NPixelHits/D");
        fourmuon_tree->Branch("mHighJPsi_NStripHits",&mHighJPsi_NStripHits,"mHighJPsi_NStripHits/D");
        fourmuon_tree->Branch("mHighJPsi_NTrackhits",&mHighJPsi_NTrackhits,"mHighJPsi_NTrackhits/D");
        fourmuon_tree->Branch("mHighJPsi_NBPixHits",&mHighJPsi_NBPixHits,"mHighJPsi_NBPixHits/D");

        fourmuon_tree->Branch("mHighJPsi_NPixLayers",&mHighJPsi_NPixLayers,"mHighJPsi_NPixLayers/D");
        fourmuon_tree->Branch("mHighJPsi_NTraLayers",&mHighJPsi_NTraLayers,"mHighJPsi_NTraLayers/D");
        fourmuon_tree->Branch("mHighJPsi_NStrLayers",&mHighJPsi_NStrLayers,"mHighJPsi_NStrLayers/D");
        fourmuon_tree->Branch("mHighJPsi_NBPixLayers",&mHighJPsi_NBPixLayers,"mHighJPsi_NBPixLayers/D");

        fourmuon_tree->Branch("mLowJPsi_NPixelHits",&mLowJPsi_NPixelHits,"mLowJPsi_NPixelHits/D");
        fourmuon_tree->Branch("mLowJPsi_NStripHits",&mLowJPsi_NStripHits,"mLowJPsi_NStripHits/D");
        fourmuon_tree->Branch("mLowJPsi_NTrackhits",&mLowJPsi_NTrackhits,"mLowJPsi_NTrackhits/D");
        fourmuon_tree->Branch("mLowJPsi_NBPixHits",&mLowJPsi_NBPixHits,"mLowJPsi_NBPixHits/D");

        fourmuon_tree->Branch("mLowJPsi_NPixLayers",&mLowJPsi_NPixLayers,"mLowJPsi_NPixLayers/D");
        fourmuon_tree->Branch("mLowJPsi_NTraLayers",&mLowJPsi_NTraLayers,"mLowJPsi_NTraLayers/D");
        fourmuon_tree->Branch("mLowJPsi_NStrLayers",&mLowJPsi_NStrLayers,"mLowJPsi_NStrLayers/D");
        fourmuon_tree->Branch("mLowJPsi_NBPixLayers",&mLowJPsi_NBPixLayers,"mLowJPsi_NBPixLayers/D");

        fourmuon_tree->Branch("mHighPhi_NPixelHits",&mHighPhi_NPixelHits,"mHighPhi_NPixelHits/D");
        fourmuon_tree->Branch("mHighPhi_NStripHits",&mHighPhi_NStripHits,"mHighPhi_NStripHits/D");
        fourmuon_tree->Branch("mHighPhi_NTrackhits",&mHighPhi_NTrackhits,"mHighPhi_NTrackhits/D");
        fourmuon_tree->Branch("mHighPhi_NBPixHits",&mHighPhi_NBPixHits,"mHighPhi_NBPixHits/D");

        fourmuon_tree->Branch("mHighPhi_NPixLayers",&mHighPhi_NPixLayers,"mHighPhi_NPixLayers/D");
        fourmuon_tree->Branch("mHighPhi_NTraLayers",&mHighPhi_NTraLayers,"mHighPhi_NTraLayers/D");
        fourmuon_tree->Branch("mHighPhi_NStrLayers",&mHighPhi_NStrLayers,"mHighPhi_NStrLayers/D");
        fourmuon_tree->Branch("mHighPhi_NBPixLayers",&mHighPhi_NBPixLayers,"mHighPhi_NBPixLayers/D");

        fourmuon_tree->Branch("mLowPhi_NPixelHits",&mLowPhi_NPixelHits,"mLowPhi_NPixelHits/D");
        fourmuon_tree->Branch("mLowPhi_NStripHits",&mLowPhi_NStripHits,"mLowPhi_NStripHits/D");
        fourmuon_tree->Branch("mLowPhi_NTrackhits",&mLowPhi_NTrackhits,"mLowPhi_NTrackhits/D");
        fourmuon_tree->Branch("mLowPhi_NBPixHits",&mLowPhi_NBPixHits,"mLowPhi_NBPixHits/D");

        fourmuon_tree->Branch("mLowPhi_NPixLayers",&mLowPhi_NPixLayers,"mLowPhi_NPixLayers/D");
        fourmuon_tree->Branch("mLowPhi_NTraLayers",&mLowPhi_NTraLayers,"mLowPhi_NTraLayers/D");
        fourmuon_tree->Branch("mLowPhi_NStrLayers",&mLowPhi_NStrLayers,"mLowPhi_NStrLayers/D");
        fourmuon_tree->Branch("mLowPhi_NBPixLayers",&mLowPhi_NBPixLayers,"mLowPhi_NBPixLayers/D");

        //Quality flags
        fourmuon_tree->Branch("mHighJPsi_isLoose",&mHighJPsi_isLoose,"mHighJPsi_isLoose/D");
        fourmuon_tree->Branch("mHighJPsi_isSoft",&mHighJPsi_isSoft,"mHighJPsi_isSoft/D");
        fourmuon_tree->Branch("mHighJPsi_isMedium",&mHighJPsi_isMedium,"mHighJPsi_isMedium/D");
        fourmuon_tree->Branch("mHighJPsi_isHighPt",&mHighJPsi_isHighPt,"mHighJPsi_isHighPt/D");

        fourmuon_tree->Branch("mLowJPsi_isLoose",&mLowJPsi_isLoose,"mLowJPsi_isLoose/D");
        fourmuon_tree->Branch("mLowJPsi_isSoft",&mLowJPsi_isSoft,"mLowJPsi_isSoft/D");
        fourmuon_tree->Branch("mLowJPsi_isMedium",&mLowJPsi_isMedium,"mLowJPsi_isMedium/D");
        fourmuon_tree->Branch("mLowJPsi_isHighPt",&mLowJPsi_isHighPt,"mLowJPsi_isHighPt/D");

        fourmuon_tree->Branch("mHighPhi_isLoose",&mHighPhi_isLoose,"mHighPhi_isLoose/D");
        fourmuon_tree->Branch("mHighPhi_isSoft",&mHighPhi_isSoft,"mHighPhi_isSoft/D");
        fourmuon_tree->Branch("mHighPhi_isMedium",&mHighPhi_isMedium,"mHighPhi_isMedium/D");
        fourmuon_tree->Branch("mHighPhi_isHighPt",&mHighPhi_isHighPt,"mHighPhi_isHighPt/D");

        fourmuon_tree->Branch("mLowPhi_isLoose",&mLowPhi_isLoose,"mLowPhi_isLoose/D");
        fourmuon_tree->Branch("mLowPhi_isSoft",&mLowPhi_isSoft,"mLowPhi_isSoft/D");
        fourmuon_tree->Branch("mLowPhi_isMedium",&mLowPhi_isMedium,"mLowPhi_isMedium/D");
        fourmuon_tree->Branch("mLowPhi_isHighPt",&mLowPhi_isHighPt,"mLowPhi_isHighPt/D");

        fourmuon_tree->Branch("mHighJPsi_isTracker",&mHighJPsi_isTracker,"mHighJPsi_isTracker/D");
        fourmuon_tree->Branch("mHighJPsi_isGlobal",&mHighJPsi_isGlobal,"mHighJPsi_isGlobal/D");
        fourmuon_tree->Branch("mLowJPsi_isTracker",&mLowJPsi_isTracker,"mLowJPsi_isTracker/D");
        fourmuon_tree->Branch("mLowJPsi_isGlobal",&mLowJPsi_isGlobal,"mLowJPsi_isGlobal/D");

        fourmuon_tree->Branch("mHighPhi_isTracker",&mHighPhi_isTracker,"mHighPhi_isTracker/D");
        fourmuon_tree->Branch("mHighPhi_isGlobal",&mHighPhi_isGlobal,"mHighPhi_isGlobal/D");
        fourmuon_tree->Branch("mLowPhi_isTracker",&mLowPhi_isTracker,"mLowPhi_isTracker/D");
        fourmuon_tree->Branch("mLowPhi_isGlobal",&mLowPhi_isGlobal,"mLowPhi_isGlobal/D");

        fourmuon_tree->Branch("mHighJPsi_type",&mHighJPsi_type,"mHighJPsi_type/D");
        fourmuon_tree->Branch("mLowJPsi_type",&mLowJPsi_type,"mLowJPsi_type/D");
        fourmuon_tree->Branch("mHighPhi_type",&mHighPhi_type,"mHighPhi_type/D");
        fourmuon_tree->Branch("mLowPhi_type",&mLowPhi_type,"mLowPhi_type/D");

        // //Dimuons
        //
        // fourmuon_tree->Branch("jpsi_triggerMatch_rf",&jpsi_triggerMatch_rf,"jpsi_triggerMatch_rf/D");
        // fourmuon_tree->Branch("phi_triggerMatch_rf",&phi_triggerMatch_rf,"phi_triggerMatch_rf/D");

        //JPsi
        fourmuon_tree->Branch("jpsi_m",&jpsi_m,"jpsi_m/D");
        fourmuon_tree->Branch("jpsi_m_rf",&jpsi_m_rf,"jpsi_m_rf/D");
        fourmuon_tree->Branch("jpsi_m_rf_c",&jpsi_m_rf_c,"jpsi_m_rf_c/D");
        fourmuon_tree->Branch("jpsi_m_rf_d_c",&jpsi_m_rf_d_c,"jpsi_m_rf_d_c/D");

        fourmuon_tree->Branch("jpsi_p",&jpsi_p,"jpsi_p/D");
        fourmuon_tree->Branch("jpsi_pt",&jpsi_pt,"jpsi_pt/D");
        fourmuon_tree->Branch("jpsi_eta",&jpsi_eta,"jpsi_eta/D");
        fourmuon_tree->Branch("jpsi_theta",&jpsi_theta,"jpsi_theta/D");
        fourmuon_tree->Branch("jpsi_y",&jpsi_y,"jpsi_y/D");
        fourmuon_tree->Branch("jpsi_e",&jpsi_e,"jpsi_e/D");

        fourmuon_tree->Branch("jpsi_dxy",&jpsi_dxy,"jpsi_dxy/D");
        fourmuon_tree->Branch("jpsi_dxyErr",&jpsi_dxyErr,"jpsi_dxyErr/D");
        fourmuon_tree->Branch("jpsi_dz",&jpsi_dz,"jpsi_dz/D");
        fourmuon_tree->Branch("jpsi_dzErr",&jpsi_dzErr,"jpsi_dzErr/D");

        fourmuon_tree->Branch("jpsi_vProb",&jpsi_vProb,"jpsi_vProb/D");
        fourmuon_tree->Branch("jpsi_vChi2",&jpsi_vChi2,"jpsi_vChi2/D");
        fourmuon_tree->Branch("jpsi_DCA",&jpsi_DCA,"jpsi_DCA/D");
        fourmuon_tree->Branch("jpsi_ctauPV",&jpsi_ctauPV,"jpsi_ctauPV/D");
        fourmuon_tree->Branch("jpsi_ctauErrPV",&jpsi_ctauErrPV,"jpsi_ctauErrPV/D");
        fourmuon_tree->Branch("jpsi_cosAlpha",&jpsi_cosAlpha,"jpsi_cosAlpha/D");
        fourmuon_tree->Branch("jpsi_lxy",&jpsi_lxy,"jpsi_lxy/D");
        fourmuon_tree->Branch("jpsi_lxyz",&jpsi_lxyz,"jpsi_lxyz/D");
        fourmuon_tree->Branch("jpsi_lxyErr",&jpsi_lxyErr,"jpsi_lxyErr/D");
        fourmuon_tree->Branch("jpsi_lxyzErr",&jpsi_lxyzErr,"jpsi_lxyzErr/D");
        fourmuon_tree->Branch("jpsi_cosAlpha3D",&jpsi_cosAlpha3D,"jpsi_cosAlpha3D/D");
        //Phi
        fourmuon_tree->Branch("phi_m",&phi_m,"phi_m/D");
        fourmuon_tree->Branch("phi_m_rf",&phi_m_rf,"phi_m_rf/D");
        fourmuon_tree->Branch("phi_m_rf_c",&phi_m_rf_c,"phi_m_rf_c/D");
        fourmuon_tree->Branch("phi_m_rf_d_c",&phi_m_rf_d_c,"phi_m_rf_d_c/D");

        fourmuon_tree->Branch("phi_p",&phi_p,"phi_p/D");
        fourmuon_tree->Branch("phi_pt",&phi_pt,"phi_pt/D");
        fourmuon_tree->Branch("phi_eta",&phi_eta,"phi_eta/D");
        fourmuon_tree->Branch("phi_theta",&phi_theta,"phi_theta/D");
        fourmuon_tree->Branch("phi_y",&phi_y,"phi_y/D");
        fourmuon_tree->Branch("phi_e",&phi_e,"phi_e/D");

        fourmuon_tree->Branch("phi_dxy",&phi_dxy,"phi_dxy/D");
        fourmuon_tree->Branch("phi_dxyErr",&phi_dxyErr,"phi_dxyErr/D");
        fourmuon_tree->Branch("phi_dz",&phi_dz,"phi_dz/D");
        fourmuon_tree->Branch("phi_dzErr",&phi_dzErr,"phi_dzErr/D");

        fourmuon_tree->Branch("phi_vProb",&phi_vProb,"phi_vProb/D");
        fourmuon_tree->Branch("phi_vChi2",&phi_vChi2,"phi_vChi2/D");
        fourmuon_tree->Branch("phi_DCA",&phi_DCA,"phi_DCA/D");
        fourmuon_tree->Branch("phi_ctauPV",&phi_ctauPV,"phi_ctauPV/D");
        fourmuon_tree->Branch("phi_ctauErrPV",&phi_ctauErrPV,"phi_ctauErrPV/D");
        fourmuon_tree->Branch("phi_cosAlpha",&phi_cosAlpha,"phi_cosAlpha/D");
        fourmuon_tree->Branch("phi_lxy",&phi_lxy,"phi_lxy/D");
        fourmuon_tree->Branch("phi_lxyz",&phi_lxyz,"phi_lxyz/D");
        fourmuon_tree->Branch("phi_lxyErr",&phi_lxyErr,"phi_lxyErr/D");
        fourmuon_tree->Branch("phi_lxyzErr",&phi_lxyzErr,"phi_lxyzErr/D");
        fourmuon_tree->Branch("phi_cosAlpha3D",&phi_cosAlpha3D,"phi_cosAlpha3D/D");


        fourmuon_tree->Branch("isBestCandidate",&isBestCandidate,"isBestCandidate/D");

        if(AddMC_)
          {

            fourmuon_tree->Branch("gen_doubledimuon_p4", "TLorentzVector",  &gen_doubledimuon_p4);
            fourmuon_tree->Branch("gen_jpsi_p4", "TLorentzVector",  &gen_jpsi_p4);
            fourmuon_tree->Branch("gen_phi_p4", "TLorentzVector",  &gen_phi_p4);

            fourmuon_tree->Branch("gen_mHighPhi_p4", "TLorentzVector",  &gen_mHighPhi_p4);
            fourmuon_tree->Branch("gen_mLowJPsi_p4",  "TLorentzVector",  &gen_mLowJPsi_p4);
            fourmuon_tree->Branch("gen_mHighJPsi_p4",  "TLorentzVector",  &gen_mHighJPsi_p4);
            fourmuon_tree->Branch("gen_mLowPhi_p4",  "TLorentzVector",  &gen_mLowPhi_p4);

            fourmuon_tree->Branch("gen_doubledimuon_pdg",&gen_doubledimuon_pdg,"gen_doubledimuon_pdg/D");
            fourmuon_tree->Branch("gen_phi_pdg",&gen_phi_pdg,"gen_phi_pdg/D");
            fourmuon_tree->Branch("gen_jpsi_pdg",&gen_jpsi_pdg,"gen_jpsi_pdg/D");

            fourmuon_tree->Branch("gen_doubledimuon_prompt",&gen_doubledimuon_prompt,"gen_doubledimuon_prompt/D");
            fourmuon_tree->Branch("gen_phi_prompt",&gen_phi_prompt,"gen_phi_prompt/D");
            fourmuon_tree->Branch("gen_jpsi_prompt",&gen_jpsi_prompt,"gen_jpsi_prompt/D");

            fourmuon_tree->Branch("gen_doubledimuon_pt",&gen_doubledimuon_pt,"gen_doubledimuon_pt/D");
            fourmuon_tree->Branch("gen_phi_pt",&gen_phi_pt,"phigen_phi_pt_pt/D");
            fourmuon_tree->Branch("gen_jpsi_pt",&gen_jpsi_pt,"gen_jpsi_pt/D");

            fourmuon_tree->Branch("gen_doubledimuon_p",&gen_doubledimuon_p,"gen_doubledimuon_p/D");
            fourmuon_tree->Branch("gen_phi_p",&gen_phi_p,"phigen_phi_p_p/D");
            fourmuon_tree->Branch("gen_jpsi_p",&gen_jpsi_p,"gen_jpsi_p/D");

            fourmuon_tree->Branch("gen_doubledimuon_eta",&gen_doubledimuon_eta,"gen_doubledimuon_eta/D");
            fourmuon_tree->Branch("gen_phi_eta",&gen_phi_eta,"gen_phi_eta/D");
            fourmuon_tree->Branch("gen_jpsi_eta",&gen_jpsi_eta,"gen_jpsi_eta/D");

          }

    genCands_ = consumes< std::vector <reco::GenParticle> >((edm::InputTag)"prunedGenParticles");
    packCands_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");
}

DoubleDiMuonRootupler::~DoubleDiMuonRootupler() {}

//
// member functions
//

bool DoubleDiMuonRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   for (size_t i=0; i< particle->numberOfMothers(); i++) {
      if (isAncestor(ancestor, particle->mother(i))) return true;
   }
   return false;
}


// ------------ method called for each event  ------------
void DoubleDiMuonRootupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
//  using namespace edm;
  using namespace std;

  edm::Handle<std::vector<pat::CompositeCandidate>> doubledimuon_cand_handle;
  iEvent.getByToken(doubledimuons_, doubledimuon_cand_handle);

  edm::Handle<std::vector<reco::Vertex >> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( triggerResults_Label , triggerResults_handle);

  numPrimaryVertices = primaryVertices_handle->size();
  run = iEvent.id().run();
  event = iEvent.id().event();

  // reco::Vertex thePrimaryV;
  // reco::Vertex theBeamSpotV;
  //
  // edm::Handle<reco::BeamSpot> theBeamSpot;
  // iEvent.getByToken(thebeamspot_,theBeamSpot);
  // reco::BeamSpot bs = *theBeamSpot;
  //
  // if ( primaryVertices_handle->begin() != primaryVertices_handle->end() ) {
  //   thePrimaryV = reco::Vertex(*(primaryVertices_handle->begin()));
  // }
  // else {
  //   thePrimaryV = reco::Vertex(bs.position(), bs.covariance3D());
  // }
  //

	// grab Trigger information
	// save it in variable trigger, trigger is an int between 0 and 7, in binary it is:
	// (pass 10)(pass 8)(pass 0)
	// ex. 7 = pass 0, 8 and 10
	// ex. 6 = pass 8, 10
        // ex. 1 = pass 0

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
   } else //std::cout << "*** NO triggerResults found " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;

   gen_doubledimuon_pdg = 0;

   edm::Handle< std::vector <reco::GenParticle> > pruned;
   iEvent.getByToken(genCands_, pruned);

   // Packed particles are all the status 1. The navigation to pruned is possible (the other direction should be made by hand)
   edm::Handle<pat::PackedGenParticleCollection> packed;
   iEvent.getByToken(packCands_,  packed);

// grabbing doubledimuon information

  //std::cout << "Debug  1" << std::endl;

  if (!doubledimuon_cand_handle.isValid()) std::cout<< "No doubledimuon information " << run << "," << event <<std::endl;
// get rf information. Notice we are just keeping combinations with succesfull vertex fit
  if (!OnlyGen_ && doubledimuon_cand_handle.isValid()) {

    doubledimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    phi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    mLowPhi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    mLowJPsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    mHighJPsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    mLowPhi_p4.SetPtEtaPhiM(0.,0.,0.,0.);

    doubledimuon_rf_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    phi_rf_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    jpsi_rf_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    mLowPhi_rf_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    mLowJPsi_rf_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    mHighJPsi_rf_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    mLowPhi_rf_p4.SetPtEtaPhiM(0.,0.,0.,0.);

    gen_doubledimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_phi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_mHighJPsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_mLowJPsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_mHighPhi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_mLowPhi_p4.SetPtEtaPhiM(0.,0.,0.,0.);

    noXCandidates = (Double_t)(doubledimuon_cand_handle->size());

    //std::cout << "Debug  2" << std::endl;

    pat::CompositeCandidate *doubledimuon_rf_cand, doubledimuon_cand, *jpsi_cand, *phi_cand;//  , *jpsi_cand, *jpsi_cand;
    pat::Muon *phiMuHigh, *phiMuLow, *jPsiMuLow, *jPsiMuHigh;

    for (unsigned int i=0; i< doubledimuon_cand_handle->size(); i++){

      doubledimuon_cand      = doubledimuon_cand_handle->at(i);

      const reco::Vertex thePrimaryV = *(doubledimuon_cand.userData<reco::Vertex>("bestPV"));

      jpsi_cand = dynamic_cast <pat::CompositeCandidate *>(doubledimuon_cand.daughter("jpsi"));
      phi_cand = dynamic_cast <pat::CompositeCandidate *>(doubledimuon_cand.daughter("phi"));

      phiMuHigh  = dynamic_cast <pat::Muon *>(phi_cand->daughter("highMuon"));
      phiMuLow   = dynamic_cast <pat::Muon *>(phi_cand->daughter("lowMuon"));
      jPsiMuLow  = dynamic_cast <pat::Muon *>(jpsi_cand->daughter("lowMuon"));
      jPsiMuHigh = dynamic_cast <pat::Muon *>(jpsi_cand->daughter("highMuon"));

      // if (doubledimuon_rf_bindx<0 || doubledimuon_rf_bindx>(int) doubledimuon_cand_handle->size()) {
      //   //std::cout << "Incorrect index for oniatt combination " << run << "," << event <<"," << doubledimuon_rf_bindx << std::endl;
      //   continue;
      // }

      doubledimuon_charge    = doubledimuon_cand.charge();
      doubledimuon_vProb     = doubledimuon_cand.userFloat("vProb");
      doubledimuon_vChi2     = doubledimuon_cand.userFloat("vChi2");
      doubledimuon_nDof      = doubledimuon_cand.userFloat("nDof");
      doubledimuon_cosAlpha  = doubledimuon_cand.userFloat("cosAlpha");
      doubledimuon_ctauPV    = doubledimuon_cand.userFloat("ctauPV");
      doubledimuon_ctauErrPV = doubledimuon_cand.userFloat("ctauErrPV");

      doubledimuon_dca_mp1mp2 = doubledimuon_cand.userFloat("dca_mp1mp2");
      doubledimuon_dca_mp1mj1 = doubledimuon_cand.userFloat("dca_mp1mj1");
      doubledimuon_dca_mp1mj2 = doubledimuon_cand.userFloat("dca_mp1mj2");
      doubledimuon_dca_mp2mj1 = doubledimuon_cand.userFloat("dca_mp2mj1");
      doubledimuon_dca_mp2mj2 = doubledimuon_cand.userFloat("dca_mp2mj2");
      doubledimuon_dca_mj1mj2 = doubledimuon_cand.userFloat("dca_mj1mj2");

      doubledimuon_cosAlphaBS = doubledimuon_cand.userFloat("cosAlphaBS");
      doubledimuon_ctauPVBS = doubledimuon_cand.userFloat("ctauPVBS");
      doubledimuon_ctauErrPVBS = doubledimuon_cand.userFloat("ctauErrPVBS");
      doubledimuon_cosAlpha = doubledimuon_cand.userFloat("cosAlpha");
      doubledimuon_ctauPV = doubledimuon_cand.userFloat("ctauPV");
      doubledimuon_ctauErrPV = doubledimuon_cand.userFloat("ctauErrPV");
      doubledimuon_cosAlphaDZ = doubledimuon_cand.userFloat("cosAlphaDZ");
      doubledimuon_ctauPVDZ = doubledimuon_cand.userFloat("ctauPVDZ");
      doubledimuon_ctauErrPVDZ = doubledimuon_cand.userFloat("ctauErrPVDZ");

      //std::cout << "Debug  3" << std::endl;
      doubledimuon_m     = doubledimuon_cand.mass();
      doubledimuon_m_rf  = doubledimuon_cand.userFloat("mass_rf");
      doubledimuon_pt    = doubledimuon_cand.pt();
      doubledimuon_eta   = doubledimuon_cand.eta();
      doubledimuon_phi   = doubledimuon_cand.phi();
      doubledimuon_y     = doubledimuon_cand.y();
      doubledimuon_vx    = doubledimuon_cand.userFloat("vtxX");
      doubledimuon_vy    = doubledimuon_cand.userFloat("vtxY");
      doubledimuon_vz    = doubledimuon_cand.userFloat("vtxZ");


      // doubledimuon_pdgid    = doubledimuon_rf_cand->userInt("phiGenPdgId");
      // doubledimuon_phipdg   = doubledimuon_rf_cand->userFloat("phiPpdlTrue");
      // doubledimuon_isprompt = doubledimuon_rf_cand->userInt("xGenPdgId");
      // doubledimuon_phippdl  = doubledimuon_rf_cand->userFloat("xGenIsPrompt");


      if (AddMC_)
        {

        gen_doubledimuon_pdg  = -1.0;
        gen_doubledimuon_prompt  = -1.0;
        gen_doubledimuon_p  = -1.0;
        gen_doubledimuon_pt  = -1.0;
        gen_doubledimuon_eta  = -1.0;

        gen_jpsi_pdg    = -1.0;
        gen_jpsi_prompt   = -1.0;
        gen_jpsi_p    = -1.0;
        gen_jpsi_pt   = -1.0;
        gen_jpsi_eta   = -1.0;
        gen_phi_pdg   = -1.0;
        gen_phi_prompt = -1.0;
        gen_phi_p =  -1.0;
        gen_phi_pt = -1.0;
        gen_phi_eta = -1.0;
        //std::cout << "Debug  4" << std::endl;
        reco::GenParticleRef genPhiMuHigh  = phiMuHigh->genParticleRef();
        reco::GenParticleRef genPhiMuLow   = phiMuLow->genParticleRef();
        reco::GenParticleRef genJPsiMuHigh = jPsiMuLow->genParticleRef();
        reco::GenParticleRef genJPsiMuLow  = jPsiMuHigh->genParticleRef();

        if(genJPsiMuLow.isNonnull())
          gen_mHighJPsi_p4.SetPtEtaPhiM(genJPsiMuLow->pt(),genJPsiMuLow->eta(),genJPsiMuLow->phi(),genJPsiMuLow->mass());
        if(genJPsiMuHigh.isNonnull())
          gen_mLowJPsi_p4.SetPtEtaPhiM(genJPsiMuHigh->pt(),genJPsiMuHigh->eta(),genJPsiMuHigh->phi(),genJPsiMuHigh->mass());
        if(genPhiMuHigh.isNonnull())
          gen_mHighPhi_p4.SetPtEtaPhiM(genPhiMuHigh->pt(),genPhiMuHigh->eta(),genPhiMuHigh->phi(),genPhiMuHigh->mass());
        if(genPhiMuLow.isNonnull())
          gen_mLowPhi_p4.SetPtEtaPhiM(genPhiMuLow->pt(),genPhiMuLow->eta(),genPhiMuLow->phi(),genPhiMuLow->mass());

        //std::cout << "Debug  5" << std::endl;
        if (genJPsiMuLow.isNonnull() && genJPsiMuHigh.isNonnull() && genPhiMuHigh.isNonnull() && genPhiMuLow.isNonnull())
        {
          if (genPhiMuHigh->numberOfMothers()>0 && genPhiMuLow->numberOfMothers()>0)
          {
          if (genJPsiMuHigh->numberOfMothers()>0 && genJPsiMuLow->numberOfMothers()>0)
          {
            reco::GenParticleRef phiMomHigh  = genPhiMuHigh->motherRef();
            reco::GenParticleRef phiMomLow   = genPhiMuLow->motherRef();
            reco::GenParticleRef jPsiMomHigh = genJPsiMuHigh->motherRef();
            reco::GenParticleRef jPsiMomLow  = genJPsiMuLow->motherRef();

            if( phiMomHigh.isNonnull() &&  phiMomLow.isNonnull() &&  jPsiMomHigh.isNonnull() &&  jPsiMomLow.isNonnull() )
            {
            bool sameJPsiMom = (jPsiMomHigh == jPsiMomLow);
            bool samePhiMom = (phiMomHigh == phiMomLow);

            if(sameJPsiMom)
            {
              gen_jpsi_p4.SetPtEtaPhiM(jPsiMomHigh->pt(),jPsiMomHigh->eta(),jPsiMomHigh->phi(),jPsiMomHigh->mass());
              gen_jpsi_pdg   = jPsiMomHigh->pdgId();
              gen_jpsi_prompt   = jPsiMomHigh->isPromptDecayed();
              gen_jpsi_p   = jPsiMomHigh->p();
              gen_jpsi_pt   = jPsiMomHigh->pt();
              gen_jpsi_eta   = jPsiMomHigh->eta();
            }

            if(samePhiMom)
            {
              gen_phi_p4.SetPtEtaPhiM(phiMomHigh->pt(),phiMomHigh->eta(),phiMomHigh->phi(),phiMomHigh->mass());
              gen_phi_pdg   = phiMomHigh->pdgId();
              gen_phi_prompt   = phiMomHigh->isPromptDecayed();
              gen_phi_p   = phiMomHigh->p();
              gen_phi_pt   = phiMomHigh->pt();
              gen_phi_eta   = phiMomHigh->eta();
            }

            //std::cout << "Debug  6" << std::endl;
            if(sameJPsiMom && samePhiMom && jPsiMomHigh->numberOfMothers()>0 && phiMomLow->numberOfMothers()>0)
            {
              reco::GenParticleRef jspiMom  = jPsiMomHigh->motherRef();
              reco::GenParticleRef phiMom   = phiMomHigh->motherRef();

              if(jspiMom==phiMom && jspiMom.isNonnull() && phiMom.isNonnull())
              {
                gen_doubledimuon_p4.SetPtEtaPhiM(jspiMom->pt(),jspiMom->eta(),jspiMom->phi(),jspiMom->mass());
                gen_doubledimuon_pdg = (Double_t) jspiMom->pdgId();
                gen_doubledimuon_prompt = (Double_t) jspiMom->isPromptDecayed();
                gen_doubledimuon_p = (Double_t) jspiMom->p();
                gen_doubledimuon_pt = (Double_t) jspiMom->pt();
                gen_doubledimuon_eta = (Double_t) jspiMom->eta();
              }

            }
          }
          }
       }
     }
   }
      reco::Candidate::LorentzVector vJPsiHigh = jpsi_cand->daughter("highMuon")->p4();
      reco::Candidate::LorentzVector vJPsiLow  = jpsi_cand->daughter("lowMuon")->p4();

      if (jpsi_cand->daughter("highMuon")->charge() < 0) {
         vJPsiHigh = jpsi_cand->daughter("lowMuon")->p4();
         vJPsiLow = jpsi_cand->daughter("highMuon")->p4();
      }
      //std::cout << "Debug  7" << std::endl;
      mLowPhi_rf_p4.SetPtEtaPhiM(vJPsiHigh.pt(), vJPsiHigh.eta(), vJPsiHigh.phi(), vJPsiHigh.mass());
      mLowJPsi_rf_p4.SetPtEtaPhiM(vJPsiLow.pt(), vJPsiLow.eta(), vJPsiLow.phi(), vJPsiLow.mass());

      reco::Candidate::LorentzVector vPhiHigh = jpsi_cand->daughter("highMuon")->p4();
      reco::Candidate::LorentzVector vPhiLow  = jpsi_cand->daughter("lowMuon")->p4();

      if (jpsi_cand->daughter("highMuon")->charge() < 0) {
         vPhiHigh = jpsi_cand->daughter("lowMuon")->p4();
         vPhiLow = jpsi_cand->daughter("highMuon")->p4();
      }

      jpsi_vProb        = (Double_t) jpsi_cand->userFloat("vProb");
      jpsi_vChi2        = (Double_t) jpsi_cand->userFloat("vNChi2");
      jpsi_DCA          = (Double_t) jpsi_cand->userFloat("DCA");
      jpsi_ctauPV       = (Double_t) jpsi_cand->userFloat("ppdlPV");
      jpsi_ctauErrPV    = (Double_t) jpsi_cand->userFloat("ppdlErrPV");
      jpsi_cosAlpha     = (Double_t) jpsi_cand->userFloat("cosAlpha");
      jpsi_lxy          = (Double_t) jpsi_cand->userFloat("l_xy");
      jpsi_lxyz         = (Double_t) jpsi_cand->userFloat("lErr_xy");
      jpsi_lxyErr       = (Double_t) jpsi_cand->userFloat("l_xyz");
      jpsi_lxyzErr      = (Double_t) jpsi_cand->userFloat("lErr_xyz");
      jpsi_cosAlpha3D   = (Double_t) jpsi_cand->userFloat("cosAlpha3D");

      mHighJPsiMatch    = jpsi_cand->userInt("highMuonTMatch");
      mLowJPsiMatch     = jpsi_cand->userInt("lowMuonTMatch");
      jpsi_triggerMatch = -std::max(-mLowJPsiMatch,mHighJPsiMatch);

      //jpsi_triggerMatch = DoubleDiMuonRootupler::isTriggerMatched(jpsi_cand);

      phi_vProb        = (Double_t) phi_cand->userFloat("vProb");
      phi_vChi2        = (Double_t) phi_cand->userFloat("vNChi2");
      phi_DCA          = (Double_t) phi_cand->userFloat("DCA");
      phi_ctauPV       = (Double_t) phi_cand->userFloat("ppdlPV");
      phi_ctauErrPV    = (Double_t) phi_cand->userFloat("ppdlErrPV");
      phi_cosAlpha     = (Double_t) phi_cand->userFloat("cosAlpha");
      phi_lxy          = (Double_t) phi_cand->userFloat("l_xy");
      phi_lxyz         = (Double_t) phi_cand->userFloat("lErr_xy");
      phi_lxyErr       = (Double_t) phi_cand->userFloat("l_xyz");
      phi_lxyzErr      = (Double_t) phi_cand->userFloat("lErr_xyz");
      phi_cosAlpha3D   = (Double_t) phi_cand->userFloat("cosAlpha3D");

      mHighPhiMatch    = (Double_t) phi_cand->userInt("highMuonTMatch");
      mLowPhiMatch     = (Double_t) phi_cand->userInt("lowMuonTMatch");
      phi_triggerMatch = -std::max(-mLowPhiMatch,mHighPhiMatch);

      const pat::Muon *jpsiHighMuon, *jspiLowMuon, *phiHighMuon, *phiLowMuon;

      jpsiHighMuon = dynamic_cast<const pat::Muon*>(jpsi_cand->daughter("highMuon"));
      jspiLowMuon = dynamic_cast<const pat::Muon*>(jpsi_cand->daughter("lowMuon"));

      mHighJPsi_p4.SetPtEtaPhiM(vJPsiHigh.pt(), vJPsiHigh.eta(), vJPsiHigh.phi(), vJPsiHigh.mass());
      mLowJPsi_p4.SetPtEtaPhiM(vJPsiLow.pt(), vJPsiLow.eta(), vJPsiLow.phi(), vJPsiLow.mass());
      //std::cout << "Debug  8" << std::endl;
      mHighJPsi_isLoose   = (Double_t) jpsiHighMuon->isLooseMuon();
      mHighJPsi_isSoft    = (Double_t) jpsiHighMuon->isSoftMuon(thePrimaryV);
      mHighJPsi_isMedium  = (Double_t) jpsiHighMuon->isMediumMuon();
      mHighJPsi_isHighPt  = (Double_t) jpsiHighMuon->isHighPtMuon(thePrimaryV);
      mHighJPsi_isTracker = (Double_t) jpsiHighMuon->isTrackerMuon();
      mHighJPsi_isGlobal  = (Double_t) jpsiHighMuon->isGlobalMuon();
      mLowJPsi_isLoose    = (Double_t) jspiLowMuon->isLooseMuon();
      mLowJPsi_isSoft     = (Double_t) jspiLowMuon->isSoftMuon(thePrimaryV);
      mLowJPsi_isMedium   = (Double_t) jspiLowMuon->isMediumMuon();
      mLowJPsi_isHighPt   = (Double_t) jspiLowMuon->isHighPtMuon(thePrimaryV);
      mLowJPsi_isTracker  = (Double_t) jspiLowMuon->isTrackerMuon();
      mLowJPsi_isGlobal   = (Double_t) jspiLowMuon->isGlobalMuon();
      mHighJPsi_type      = (Double_t) jpsiHighMuon->type();
      mLowJPsi_type       = (Double_t) jspiLowMuon->type();

      //double kmass = 0.4936770;
      doubledimuon_p4.SetPtEtaPhiM(doubledimuon_cand.pt(),doubledimuon_cand.eta(),doubledimuon_cand.phi(),doubledimuon_cand.mass());
      jpsi_p4.SetPtEtaPhiM(jpsi_cand->pt(),jpsi_cand->eta(),jpsi_cand->phi(),jpsi_cand->mass());
      phi_p4.SetPtEtaPhiM(phi_cand->pt(), phi_cand->eta(), phi_cand->phi(), phi_cand->mass());

      vPhiHigh = phi_cand->daughter("highMuon")->p4();
      vPhiLow = phi_cand->daughter("lowMuon")->p4();

      phiLowMuon = dynamic_cast<const pat::Muon*>(phi_cand->daughter("highMuon"));
      phiHighMuon = dynamic_cast<const pat::Muon*>(phi_cand->daughter("lowMuon"));

      // highDiM_m        = (Double_t) jpsi_cand->mass();
      // highDiM_pt       = (Double_t) jpsi_cand->pt();
      // lowDiM_m         = (Double_t) phi_cand->mass();
      // lowDiM_pt        = (Double_t) phi_cand->pt();

      jpsi_m      =   (Double_t) jpsi_cand->mass() ;
      jpsi_p      =   (Double_t) jpsi_cand->p() ;
      jpsi_theta  =   (Double_t) jpsi_cand->theta() ;
      jpsi_eta    =   (Double_t) jpsi_cand->eta() ;
      jpsi_y      =   (Double_t) jpsi_cand->y() ;
      jpsi_e      =   (Double_t) jpsi_cand->energy();
      // jpsi_dxy      =   jpsi_cand->dxy();

      phi_m      =   (Double_t) phi_cand->mass() ;
      phi_p      =   (Double_t) phi_cand->p() ;
      phi_theta  =   (Double_t) phi_cand->theta() ;
      phi_eta    =   (Double_t) phi_cand->eta() ;
      phi_y      =   (Double_t) phi_cand->y() ;
      phi_e      =   (Double_t) phi_cand->energy();
      // phi_dxy      =   phi_cand->dxy();

      //jpsi_dxyEr      =   jpsi_cand. ;
      //jpsi_d      =   jpsi_cand. ;
      //jpsi_dzErr      =   jpsi_cand. ;
      mHighJPsi_pt     = vJPsiHigh.pt();
      mLowJPsi_pt      = vJPsiLow.pt();
      mHighPhi_pt      = vPhiHigh.pt();
      mLowPhi_pt       = vPhiLow.pt();

      mHighPhi_isLoose    = (Double_t)  phiHighMuon->isLooseMuon();
      mHighPhi_isSoft     = (Double_t)  phiHighMuon->isSoftMuon(thePrimaryV);
      mHighPhi_isMedium   = (Double_t) phiHighMuon->isMediumMuon();
      mHighPhi_isHighPt   = (Double_t) phiHighMuon->isHighPtMuon(thePrimaryV);
      mHighPhi_isTracker  = (Double_t) phiHighMuon->isTrackerMuon();
      mHighPhi_isGlobal   = (Double_t) phiHighMuon->isGlobalMuon();
      mLowPhi_isLoose     = (Double_t) phiLowMuon->isLooseMuon();
      mLowPhi_isSoft      = (Double_t) phiLowMuon->isSoftMuon(thePrimaryV);
      mLowPhi_isMedium    = (Double_t) phiLowMuon->isMediumMuon();
      mLowPhi_isHighPt    = (Double_t) phiLowMuon->isHighPtMuon(thePrimaryV);
      mLowPhi_isTracker   = (Double_t) phiLowMuon->isTrackerMuon();
      mLowPhi_isGlobal    = (Double_t) phiLowMuon->isGlobalMuon();
      mHighPhi_type       = (Double_t)  phiHighMuon->type();
      mLowPhi_type        = (Double_t) phiLowMuon->type();

      mHighPhi_p4.SetPtEtaPhiM(vPhiHigh.pt(), vPhiHigh.eta(), vPhiHigh.phi(), vPhiHigh.mass());
      mLowPhi_p4.SetPtEtaPhiM(vPhiLow.pt(), vPhiLow.eta(), vPhiLow.phi(), vPhiLow.mass());

      if (doubledimuon_cand.userFloat("has_ref") >= 0)
      {
        doubledimuon_rf_cand = dynamic_cast <pat::CompositeCandidate *>(doubledimuon_cand.daughter("ref_cand"));
        doubledimuon_rf_p4.SetPtEtaPhiM(doubledimuon_rf_cand->pt(),doubledimuon_rf_cand->eta(),doubledimuon_rf_cand->phi(),doubledimuon_rf_cand->mass());
        jpsi_rf_p4.SetPtEtaPhiM(doubledimuon_rf_cand->daughter("jpsi")->pt(),doubledimuon_rf_cand->daughter("jpsi")->eta(),
                                doubledimuon_rf_cand->daughter("jpsi")->phi(),doubledimuon_rf_cand->daughter("jpsi")->mass());
        jpsi_rf_p4.SetPtEtaPhiM(doubledimuon_rf_cand->daughter("phi")->pt(),doubledimuon_rf_cand->daughter("phi")->eta(),
                                doubledimuon_rf_cand->daughter("phi")->phi(),doubledimuon_rf_cand->daughter("phi")->mass());

        jpsi_cand = dynamic_cast <pat::CompositeCandidate *>(doubledimuon_rf_cand->daughter("jpsi"));
        jpsi_cand = dynamic_cast <pat::CompositeCandidate *>(doubledimuon_rf_cand->daughter("phi"));

        vJPsiHigh = jpsi_cand->daughter("highMuon")->p4();
        vJPsiLow = jpsi_cand->daughter("lowMuon")->p4();

        mLowJPsi_rf_p4.SetPtEtaPhiM(vJPsiHigh.pt(), vJPsiHigh.eta(), vJPsiHigh.phi(), vJPsiHigh.mass());
        mLowJPsi_rf_p4.SetPtEtaPhiM(vJPsiLow.pt(), vJPsiLow.eta(), vJPsiLow.phi(), vJPsiLow.mass());

        vPhiHigh = jpsi_cand->daughter("highMuon")->p4();
        vPhiLow = jpsi_cand->daughter("lowMuon")->p4();

        mHighPhi_rf_p4.SetPtEtaPhiM(vPhiHigh.pt(), vPhiHigh.eta(), vPhiHigh.phi(), vPhiHigh.mass());
        mLowPhi_rf_p4.SetPtEtaPhiM(vPhiLow.pt(), vPhiLow.eta(), vPhiLow.phi(), vPhiLow.mass());

        doubledimuon_m_rf_c= doubledimuon_rf_cand->mass();

        doubledimuon_rf_vProb = doubledimuon_cand.userFloat("vProb_ref");
        doubledimuon_rf_vChi2 = doubledimuon_cand.userFloat("vChi2_ref");
        doubledimuon_rf_nDof = doubledimuon_cand.userFloat("nDof_ref");
        doubledimuon_rf_cosAlpha = doubledimuon_cand.userFloat("cosAlpha_ref");
        doubledimuon_rf_ctauPV = doubledimuon_cand.userFloat("ctauPV_ref");
        doubledimuon_rf_ctauErrPV = doubledimuon_cand.userFloat("ctauErrPV_ref");
        //std::cout << "Debug  9" << std::endl;
      }

      if (doubledimuon_cand.userFloat("has_const_ref") >= 0)
      {
        doubledimuon_rf_cand = dynamic_cast <pat::CompositeCandidate *>(doubledimuon_cand.daughter("ref_const_cand"));
        doubledimuon_rf_const_p4.SetPtEtaPhiM(doubledimuon_rf_cand->pt(),doubledimuon_rf_cand->eta(),doubledimuon_rf_cand->phi(),doubledimuon_rf_cand->mass());

        doubledimuon_m_rf_d_c= doubledimuon_rf_cand->mass();

        doubledimuon_rf_vProb = doubledimuon_cand.userFloat("vProb_const_ref");
        doubledimuon_rf_vChi2 = doubledimuon_cand.userFloat("vChi2_const_ref");
        doubledimuon_rf_nDof = doubledimuon_cand.userFloat("nDof_const_ref");
        doubledimuon_rf_cosAlpha = doubledimuon_cand.userFloat("cosAlpha_const_ref");
        doubledimuon_rf_ctauPV = doubledimuon_cand.userFloat("ctauPV_const_ref");
        doubledimuon_rf_ctauErrPV = doubledimuon_cand.userFloat("ctauErrPV_const_ref");

      }
      //std::cout << "Debug  10" << std::endl;
      fourmuon_tree->Fill();

      if (OnlyBest_) break;
      else if(i==0)
      isBestCandidate = false;

      }

  }

}

// ------------ method called once each job just before starting event loop  ------------
void DoubleDiMuonRootupler::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void DoubleDiMuonRootupler::endJob() {}

// ------------ method called when starting to processes a run  ------------
void DoubleDiMuonRootupler::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void DoubleDiMuonRootupler::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void DoubleDiMuonRootupler::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void DoubleDiMuonRootupler::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DoubleDiMuonRootupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DoubleDiMuonRootupler);
