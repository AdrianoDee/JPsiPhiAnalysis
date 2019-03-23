#define FiveTracks_cxx

// The class definition in FiveTracks.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("FiveTracks.C")
// root> T->Process("FiveTracks.C","some options")
// root> T->Process("FiveTracks.C+")
//


#include "FiveTracks.h"
#include <TH2.h>
#include <TStyle.h>


Double_t out;

//TTree* outTuple = new TTree("2mu2kSkimmedTree","2mu2kSkimmedTree");
//TBranch* b = outTuple->Branch("out",       &out,        "out/D");

void FiveTracks::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
}

void FiveTracks::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  std::string outputString = "2mu3trk_five_tree.root";
  OutFile = new TProofOutputFile( outputString.data() );
  fOut = OutFile->OpenFile("RECREATE");
  if (!(fOut=OutFile->OpenFile("RECREATE")))
  {
    Warning("SlaveBegin","Problems opening file: %s%s", OutFile->GetDir(), OutFile->GetFileName() );
  }

  ////////////////// Histograms //////////////////
  JPsi_mass = 3.096916; /// pdg mass
  Phi_mass = 1.019455; /// pdg mass
  Phi_mean = 1.019723;
  Phi_sigma = 2.35607e-03;//2.28400e-03;

  outTree = new TTree("2mu2kSkimmedTree","2mu2kSkimmedTree");

  outTree->Branch("event",        &out_event,     "event");
  outTree->Branch("lumi",         &out_lumi,      "lumi");
  outTree->Branch("numPrimaryVertices",   &out_numPrimaryVertices,        "numPrimaryVertices");
  outTree->Branch("trigger",      &out_trigger,   "trigger");
  outTree->Branch("noFiveCandidates",     &out_noFiveCandidates,  "noFiveCandidates");
  outTree->Branch("dimuonditrk_id",       &out_dimuonditrk_id,    "dimuonditrk_id");



  outTree->Branch("dimuonditrk_m",        &out_dimuonditrk_m,     "dimuonditrk_m");
  outTree->Branch("dimuonditrk_pt",       &out_dimuonditrk_pt,    "dimuonditrk_pt");
  outTree->Branch("dimuonditrk_eta",      &out_dimuonditrk_eta,   "dimuonditrk_eta");
  outTree->Branch("dimuonditrk_phi",      &out_dimuonditrk_phi,   "dimuonditrk_phi");
  outTree->Branch("dimuonditrk_p",        &out_dimuonditrk_p,     "dimuonditrk_p");
  outTree->Branch("dimuon_m",     &out_dimuon_m,  "dimuon_m");
  outTree->Branch("dimuon_pt",    &out_dimuon_pt,         "dimuon_pt");
  outTree->Branch("dimuon_eta",   &out_dimuon_eta,        "dimuon_eta");
  outTree->Branch("dimuon_phi",   &out_dimuon_phi,        "dimuon_phi");
  outTree->Branch("dimuon_p",     &out_dimuon_p,  "dimuon_p");
  outTree->Branch("highKaonMatch",        &out_highKaonMatch,     "highKaonMatch");
  outTree->Branch("lowKaonMatch",         &out_lowKaonMatch,      "lowKaonMatch");
  outTree->Branch("lowMuonMatch",         &out_lowMuonMatch,      "lowMuonMatch");
  outTree->Branch("highMuonMatch",        &out_highMuonMatch,     "highMuonMatch");
  outTree->Branch("ditrak_m",     &out_ditrak_m,  "ditrak_m");
  outTree->Branch("ditrakOne_pt",         &out_ditrakOne_pt,      "ditrakOne_pt");
  outTree->Branch("ditrakOne_eta",        &out_ditrakOne_eta,     "ditrakOne_eta");
  outTree->Branch("ditrakOne_phi",        &out_ditrakOne_phi,     "ditrakOne_phi");
  outTree->Branch("ditrakOne_p",  &out_ditrakOne_p,       "ditrakOne_p");
  outTree->Branch("ditrakTwo_pt",         &out_ditrakTwo_pt,      "ditrakTwo_pt");
  outTree->Branch("ditrakTwo_eta",        &out_ditrakTwo_eta,     "ditrakTwo_eta");
  outTree->Branch("ditrakTwo_phi",        &out_ditrakTwo_phi,     "ditrakTwo_phi");
  outTree->Branch("ditrakTwo_p",  &out_ditrakTwo_p,       "ditrakTwo_p");
  outTree->Branch("ditrakThree_pt",       &out_ditrakThree_pt,    "ditrakThree_pt");
  outTree->Branch("ditrakThree_eta",      &out_ditrakThree_eta,   "ditrakThree_eta");
  outTree->Branch("ditrakThree_phi",      &out_ditrakThree_phi,   "ditrakThree_phi");
  outTree->Branch("ditrakThree_p",        &out_ditrakThree_p,     "ditrakThree_p");
  outTree->Branch("highTrack_pt",         &out_highTrack_pt,      "highTrack_pt");
  outTree->Branch("highTrack_eta",        &out_highTrack_eta,     "highTrack_eta");
  outTree->Branch("highTrack_phi",        &out_highTrack_phi,     "highTrack_phi");
  outTree->Branch("highTrack_charge",     &out_highTrack_charge,  "highTrack_charge");
  outTree->Branch("lowTrack_pt",  &out_lowTrack_pt,       "lowTrack_pt");
  outTree->Branch("lowTrack_eta",         &out_lowTrack_eta,      "lowTrack_eta");
  outTree->Branch("lowTrack_phi",         &out_lowTrack_phi,      "lowTrack_phi");
  outTree->Branch("lowTrack_charge",      &out_lowTrack_charge,   "lowTrack_charge");
  outTree->Branch("thirdTrack_pt",        &out_thirdTrack_pt,     "thirdTrack_pt");
  outTree->Branch("thirdTrack_eta",       &out_thirdTrack_eta,    "thirdTrack_eta");
  outTree->Branch("thirdTrack_phi",       &out_thirdTrack_phi,    "thirdTrack_phi");
  outTree->Branch("thirdTrack_charge",    &out_thirdTrack_charge,         "thirdTrack_charge");
  outTree->Branch("thirdTrack_dz",        &out_thirdTrack_dz,     "thirdTrack_dz");
  outTree->Branch("thirdTrack_dxy",       &out_thirdTrack_dxy,    "thirdTrack_dxy");
  outTree->Branch("dimuonDiTrkOne_pt",    &out_dimuonDiTrkOne_pt,         "dimuonDiTrkOne_pt");
  outTree->Branch("dimuonDiTrkOne_eta",   &out_dimuonDiTrkOne_eta,        "dimuonDiTrkOne_eta");
  outTree->Branch("dimuonDiTrkOne_phi",   &out_dimuonDiTrkOne_phi,        "dimuonDiTrkOne_phi");
  outTree->Branch("dimuonDiTrkOne_charge",        &out_dimuonDiTrkOne_charge,     "dimuonDiTrkOne_charge");
  outTree->Branch("dimuonDiTrkTwo_pt",    &out_dimuonDiTrkTwo_pt,         "dimuonDiTrkTwo_pt");
  outTree->Branch("dimuonDiTrkTwo_eta",   &out_dimuonDiTrkTwo_eta,        "dimuonDiTrkTwo_eta");
  outTree->Branch("dimuonDiTrkTwo_phi",   &out_dimuonDiTrkTwo_phi,        "dimuonDiTrkTwo_phi");
  outTree->Branch("dimuonDiTrkTwo_charge",        &out_dimuonDiTrkTwo_charge,     "dimuonDiTrkTwo_charge");
  outTree->Branch("dimuonDiTrkThree_pt",  &out_dimuonDiTrkThree_pt,       "dimuonDiTrkThree_pt");
  outTree->Branch("dimuonDiTrkThree_eta",         &out_dimuonDiTrkThree_eta,      "dimuonDiTrkThree_eta");
  outTree->Branch("dimuonDiTrkThree_phi",         &out_dimuonDiTrkThree_phi,      "dimuonDiTrkThree_phi");
  outTree->Branch("dimuonDiTrkThree_charge",      &out_dimuonDiTrkThree_charge,   "dimuonDiTrkThree_charge");
  outTree->Branch("psiPrimeSame_pt",      &out_psiPrimeSame_pt,   "psiPrimeSame_pt");
  outTree->Branch("psiPrimeSame_eta",     &out_psiPrimeSame_eta,  "psiPrimeSame_eta");
  outTree->Branch("psiPrimeSame_phi",     &out_psiPrimeSame_phi,  "psiPrimeSame_phi");
  outTree->Branch("psiPrimeSame_n",       &out_psiPrimeSame_n,    "psiPrimeSame_n");
  outTree->Branch("psiPrimeSame_p_pt",    &out_psiPrimeSame_p_pt,         "psiPrimeSame_p_pt");
  outTree->Branch("psiPrimeSame_p_eta",   &out_psiPrimeSame_p_eta,        "psiPrimeSame_p_eta");
  outTree->Branch("psiPrimeSame_p_phi",   &out_psiPrimeSame_p_phi,        "psiPrimeSame_p_phi");
  outTree->Branch("psiPrimeSame_p_n",     &out_psiPrimeSame_p_n,  "psiPrimeSame_p_n");
  outTree->Branch("psiPrimeSame_m_pt",    &out_psiPrimeSame_m_pt,         "psiPrimeSame_m_pt");
  outTree->Branch("psiPrimeSame_m_eta",   &out_psiPrimeSame_m_eta,        "psiPrimeSame_m_eta");
  outTree->Branch("psiPrimeSame_m_phi",   &out_psiPrimeSame_m_phi,        "psiPrimeSame_m_phi");
  outTree->Branch("psiPrimeSame_m_n",     &out_psiPrimeSame_m_n,  "psiPrimeSame_m_n");
  outTree->Branch("psiPrimeMixed_pt",     &out_psiPrimeMixed_pt,  "psiPrimeMixed_pt");
  outTree->Branch("psiPrimeMixed_eta",    &out_psiPrimeMixed_eta,         "psiPrimeMixed_eta");
  outTree->Branch("psiPrimeMixed_phi",    &out_psiPrimeMixed_phi,         "psiPrimeMixed_phi");
  outTree->Branch("psiPrimeMixed_n",      &out_psiPrimeMixed_n,   "psiPrimeMixed_n");
  outTree->Branch("psiPrimeMixed_p_pt",   &out_psiPrimeMixed_p_pt,        "psiPrimeMixed_p_pt");
  outTree->Branch("psiPrimeMixed_p_eta",  &out_psiPrimeMixed_p_eta,       "psiPrimeMixed_p_eta");
  outTree->Branch("psiPrimeMixed_p_phi",  &out_psiPrimeMixed_p_phi,       "psiPrimeMixed_p_phi");
  outTree->Branch("psiPrimeMixed_p_n",    &out_psiPrimeMixed_p_n,         "psiPrimeMixed_p_n");
  outTree->Branch("psiPrimeMixed_m_pt",   &out_psiPrimeMixed_m_pt,        "psiPrimeMixed_m_pt");
  outTree->Branch("psiPrimeMixed_m_eta",  &out_psiPrimeMixed_m_eta,       "psiPrimeMixed_m_eta");
  outTree->Branch("psiPrimeMixed_m_phi",  &out_psiPrimeMixed_m_phi,       "psiPrimeMixed_m_phi");
  outTree->Branch("psiPrimeMixed_m_n",    &out_psiPrimeMixed_m_n,         "psiPrimeMixed_m_n");
  outTree->Branch("psiPrimeSame_ditrak_pt",       &out_psiPrimeSame_ditrak_pt,    "psiPrimeSame_ditrak_pt");
  outTree->Branch("psiPrimeSame_ditrak_eta",      &out_psiPrimeSame_ditrak_eta,   "psiPrimeSame_ditrak_eta");
  outTree->Branch("psiPrimeSame_ditrak_phi",      &out_psiPrimeSame_ditrak_phi,   "psiPrimeSame_ditrak_phi");
  outTree->Branch("psiPrimeSame_ditrak_n",        &out_psiPrimeSame_ditrak_n,     "psiPrimeSame_ditrak_n");
  outTree->Branch("psiPrimeMixed_ditrak_pt",      &out_psiPrimeMixed_ditrak_pt,   "psiPrimeMixed_ditrak_pt");
  outTree->Branch("psiPrimeMixed_ditrak_eta",     &out_psiPrimeMixed_ditrak_eta,  "psiPrimeMixed_ditrak_eta");
  outTree->Branch("psiPrimeMixed_ditrak_phi",     &out_psiPrimeMixed_ditrak_phi,  "psiPrimeMixed_ditrak_phi");
  outTree->Branch("psiPrimeMixed_ditrak_n",       &out_psiPrimeMixed_ditrak_n,    "psiPrimeMixed_ditrak_n");
  outTree->Branch("triTrak_pt",   &out_triTrak_pt,        "triTrak_pt");
  outTree->Branch("triTrak_eta",  &out_triTrak_eta,       "triTrak_eta");
  outTree->Branch("triTrak_phi",  &out_triTrak_phi,       "triTrak_phi");
  outTree->Branch("triTrak_charge",       &out_triTrak_charge,    "triTrak_charge");
  outTree->Branch("mass_kkk",     &out_mass_kkk,  "mass_kkk");
  outTree->Branch("mass_ref_kkk",         &out_mass_ref_kkk,      "mass_ref_kkk");
  outTree->Branch("vProb_kkk",    &out_vProb_kkk,         "vProb_kkk");
  outTree->Branch("nDof_kkk",     &out_nDof_kkk,  "nDof_kkk");
  outTree->Branch("vChi2_kkk",    &out_vChi2_kkk,         "vChi2_kkk");
  outTree->Branch("ctau_kkk",     &out_ctau_kkk,  "ctau_kkk");
  outTree->Branch("ctauErr_kkk",  &out_ctauErr_kkk,       "ctauErr_kkk");
  outTree->Branch("cosAlpha_kkk",         &out_cosAlpha_kkk,      "cosAlpha_kkk");
  outTree->Branch("onePsiPrime_m_kkk",    &out_onePsiPrime_m_kkk,         "onePsiPrime_m_kkk");
  outTree->Branch("twoPsiPrime_m_kkk",    &out_twoPsiPrime_m_kkk,         "twoPsiPrime_m_kkk");
  outTree->Branch("onePsiPrime_p_mkkk",   &out_onePsiPrime_p_mkkk,        "onePsiPrime_p_mkkk");
  outTree->Branch("onePsiPrime_m_mkkk",   &out_onePsiPrime_m_mkkk,        "onePsiPrime_m_mkkk");
  outTree->Branch("twoPsiPrime_p_mkkk",   &out_twoPsiPrime_p_mkkk,        "twoPsiPrime_p_mkkk");
  outTree->Branch("twoPsiPrime_m_mkkk",   &out_twoPsiPrime_m_mkkk,        "twoPsiPrime_m_mkkk");
  outTree->Branch("dimuonDiTrkOne_m_kkk",         &out_dimuonDiTrkOne_m_kkk,      "dimuonDiTrkOne_m_kkk");
  outTree->Branch("dimuonDiTrkTwo_m_kkk",         &out_dimuonDiTrkTwo_m_kkk,      "dimuonDiTrkTwo_m_kkk");
  outTree->Branch("dimuonDiTrkThree_m_kkk",       &out_dimuonDiTrkThree_m_kkk,    "dimuonDiTrkThree_m_kkk");
  outTree->Branch("ditrakOne_m_kkk",      &out_ditrakOne_m_kkk,   "ditrakOne_m_kkk");
  outTree->Branch("ditrakTwo_m_kkk",      &out_ditrakTwo_m_kkk,   "ditrakTwo_m_kkk");
  outTree->Branch("ditrakThree_m_kkk",    &out_ditrakThree_m_kkk,         "ditrakThree_m_kkk");
  outTree->Branch("trackOne_m_kkk",       &out_trackOne_m_kkk,    "trackOne_m_kkk");
  outTree->Branch("trackTwo_m_kkk",       &out_trackTwo_m_kkk,    "trackTwo_m_kkk");
  outTree->Branch("trackThree_m_kkk",     &out_trackThree_m_kkk,  "trackThree_m_kkk");

  outTree->Branch("triTrak_m_kkk",        &out_triTrak_m_kkk,     "triTrak_m_kkk");
  outTree->Branch("mass_ppk",     &out_mass_ppk,  "mass_ppk");
  outTree->Branch("mass_ref_ppk",         &out_mass_ref_ppk,      "mass_ref_ppk");
  outTree->Branch("vProb_ppk",    &out_vProb_ppk,         "vProb_ppk");
  outTree->Branch("nDof_ppk",     &out_nDof_ppk,  "nDof_ppk");
  outTree->Branch("vChi2_ppk",    &out_vChi2_ppk,         "vChi2_ppk");
  outTree->Branch("ctau_ppk",     &out_ctau_ppk,  "ctau_ppk");
  outTree->Branch("ctauErr_ppk",  &out_ctauErr_ppk,       "ctauErr_ppk");
  outTree->Branch("cosAlpha_ppk",         &out_cosAlpha_ppk,      "cosAlpha_ppk");
  outTree->Branch("onePsiPrime_m_ppk",    &out_onePsiPrime_m_ppk,         "onePsiPrime_m_ppk");
  outTree->Branch("twoPsiPrime_m_ppk",    &out_twoPsiPrime_m_ppk,         "twoPsiPrime_m_ppk");
  outTree->Branch("onePsiPrime_p_mppk",   &out_onePsiPrime_p_mppk,        "onePsiPrime_p_mppk");
  outTree->Branch("onePsiPrime_m_mppk",   &out_onePsiPrime_m_mppk,        "onePsiPrime_m_mppk");
  outTree->Branch("twoPsiPrime_p_mppk",   &out_twoPsiPrime_p_mppk,        "twoPsiPrime_p_mppk");
  outTree->Branch("twoPsiPrime_m_mppk",   &out_twoPsiPrime_m_mppk,        "twoPsiPrime_m_mppk");
  outTree->Branch("dimuonDiTrkOne_m_ppk",         &out_dimuonDiTrkOne_m_ppk,      "dimuonDiTrkOne_m_ppk");
  outTree->Branch("dimuonDiTrkTwo_m_ppk",         &out_dimuonDiTrkTwo_m_ppk,      "dimuonDiTrkTwo_m_ppk");
  outTree->Branch("dimuonDiTrkThree_m_ppk",       &out_dimuonDiTrkThree_m_ppk,    "dimuonDiTrkThree_m_ppk");
  outTree->Branch("ditrakOne_m_ppk",      &out_ditrakOne_m_ppk,   "ditrakOne_m_ppk");
  outTree->Branch("ditrakTwo_m_ppk",      &out_ditrakTwo_m_ppk,   "ditrakTwo_m_ppk");
  outTree->Branch("ditrakThree_m_ppk",    &out_ditrakThree_m_ppk,         "ditrakThree_m_ppk");
  outTree->Branch("trackOne_m_ppk",       &out_trackOne_m_ppk,    "trackOne_m_ppk");
  outTree->Branch("trackTwo_m_ppk",       &out_trackTwo_m_ppk,    "trackTwo_m_ppk");
  outTree->Branch("trackThree_m_ppk",     &out_trackThree_m_ppk,  "trackThree_m_ppk");

  outTree->Branch("triTrak_m_ppk",        &out_triTrak_m_ppk,     "triTrak_m_ppk");
  outTree->Branch("mass_kpp",     &out_mass_kpp,  "mass_kpp");
  outTree->Branch("mass_ref_kpp",         &out_mass_ref_kpp,      "mass_ref_kpp");
  outTree->Branch("vProb_kpp",    &out_vProb_kpp,         "vProb_kpp");
  outTree->Branch("nDof_kpp",     &out_nDof_kpp,  "nDof_kpp");
  outTree->Branch("vChi2_kpp",    &out_vChi2_kpp,         "vChi2_kpp");
  outTree->Branch("ctau_kpp",     &out_ctau_kpp,  "ctau_kpp");
  outTree->Branch("ctauErr_kpp",  &out_ctauErr_kpp,       "ctauErr_kpp");
  outTree->Branch("cosAlpha_kpp",         &out_cosAlpha_kpp,      "cosAlpha_kpp");
  outTree->Branch("onePsiPrime_m_kpp",    &out_onePsiPrime_m_kpp,         "onePsiPrime_m_kpp");
  outTree->Branch("twoPsiPrime_m_kpp",    &out_twoPsiPrime_m_kpp,         "twoPsiPrime_m_kpp");
  outTree->Branch("onePsiPrime_p_mkpp",   &out_onePsiPrime_p_mkpp,        "onePsiPrime_p_mkpp");
  outTree->Branch("onePsiPrime_m_mkpp",   &out_onePsiPrime_m_mkpp,        "onePsiPrime_m_mkpp");
  outTree->Branch("twoPsiPrime_p_mkpp",   &out_twoPsiPrime_p_mkpp,        "twoPsiPrime_p_mkpp");
  outTree->Branch("twoPsiPrime_m_mkpp",   &out_twoPsiPrime_m_mkpp,        "twoPsiPrime_m_mkpp");
  outTree->Branch("dimuonDiTrkOne_m_kpp",         &out_dimuonDiTrkOne_m_kpp,      "dimuonDiTrkOne_m_kpp");
  outTree->Branch("dimuonDiTrkTwo_m_kpp",         &out_dimuonDiTrkTwo_m_kpp,      "dimuonDiTrkTwo_m_kpp");
  outTree->Branch("dimuonDiTrkThree_m_kpp",       &out_dimuonDiTrkThree_m_kpp,    "dimuonDiTrkThree_m_kpp");
  outTree->Branch("ditrakOne_m_kpp",      &out_ditrakOne_m_kpp,   "ditrakOne_m_kpp");
  outTree->Branch("ditrakTwo_m_kpp",      &out_ditrakTwo_m_kpp,   "ditrakTwo_m_kpp");
  outTree->Branch("ditrakThree_m_kpp",    &out_ditrakThree_m_kpp,         "ditrakThree_m_kpp");
  outTree->Branch("trackOne_m_kpp",       &out_trackOne_m_kpp,    "trackOne_m_kpp");
  outTree->Branch("trackTwo_m_kpp",       &out_trackTwo_m_kpp,    "trackTwo_m_kpp");
  outTree->Branch("trackThree_m_kpp",     &out_trackThree_m_kpp,  "trackThree_m_kpp");

  outTree->Branch("triTrak_m_kpp",        &out_triTrak_m_kpp,     "triTrak_m_kpp");
  outTree->Branch("mass_pkp",     &out_mass_pkp,  "mass_pkp");
  outTree->Branch("mass_ref_pkp",         &out_mass_ref_pkp,      "mass_ref_pkp");
  outTree->Branch("vProb_pkp",    &out_vProb_pkp,         "vProb_pkp");
  outTree->Branch("nDof_pkp",     &out_nDof_pkp,  "nDof_pkp");
  outTree->Branch("vChi2_pkp",    &out_vChi2_pkp,         "vChi2_pkp");
  outTree->Branch("ctau_pkp",     &out_ctau_pkp,  "ctau_pkp");
  outTree->Branch("ctauErr_pkp",  &out_ctauErr_pkp,       "ctauErr_pkp");
  outTree->Branch("cosAlpha_pkp",         &out_cosAlpha_pkp,      "cosAlpha_pkp");
  outTree->Branch("onePsiPrime_m_pkp",    &out_onePsiPrime_m_pkp,         "onePsiPrime_m_pkp");
  outTree->Branch("twoPsiPrime_m_pkp",    &out_twoPsiPrime_m_pkp,         "twoPsiPrime_m_pkp");
  outTree->Branch("onePsiPrime_p_mpkp",   &out_onePsiPrime_p_mpkp,        "onePsiPrime_p_mpkp");
  outTree->Branch("onePsiPrime_m_mpkp",   &out_onePsiPrime_m_mpkp,        "onePsiPrime_m_mpkp");
  outTree->Branch("twoPsiPrime_p_mpkp",   &out_twoPsiPrime_p_mpkp,        "twoPsiPrime_p_mpkp");
  outTree->Branch("twoPsiPrime_m_mpkp",   &out_twoPsiPrime_m_mpkp,        "twoPsiPrime_m_mpkp");
  outTree->Branch("dimuonDiTrkOne_m_pkp",         &out_dimuonDiTrkOne_m_pkp,      "dimuonDiTrkOne_m_pkp");
  outTree->Branch("dimuonDiTrkTwo_m_pkp",         &out_dimuonDiTrkTwo_m_pkp,      "dimuonDiTrkTwo_m_pkp");
  outTree->Branch("dimuonDiTrkThree_m_pkp",       &out_dimuonDiTrkThree_m_pkp,    "dimuonDiTrkThree_m_pkp");
  outTree->Branch("ditrakOne_m_pkp",      &out_ditrakOne_m_pkp,   "ditrakOne_m_pkp");
  outTree->Branch("ditrakTwo_m_pkp",      &out_ditrakTwo_m_pkp,   "ditrakTwo_m_pkp");
  outTree->Branch("ditrakThree_m_pkp",    &out_ditrakThree_m_pkp,         "ditrakThree_m_pkp");
  outTree->Branch("trackOne_m_pkp",       &out_trackOne_m_pkp,    "trackOne_m_pkp");
  outTree->Branch("trackTwo_m_pkp",       &out_trackTwo_m_pkp,    "trackTwo_m_pkp");
  outTree->Branch("trackThree_m_pkp",     &out_trackThree_m_pkp,  "trackThree_m_pkp");

  outTree->Branch("triTrak_m_pkp",        &out_triTrak_m_pkp,     "triTrak_m_pkp");
  outTree->Branch("mass_ppp",     &out_mass_ppp,  "mass_ppp");
  outTree->Branch("mass_ref_ppp",         &out_mass_ref_ppp,      "mass_ref_ppp");
  outTree->Branch("vProb_ppp",    &out_vProb_ppp,         "vProb_ppp");
  outTree->Branch("nDof_ppp",     &out_nDof_ppp,  "nDof_ppp");
  outTree->Branch("vChi2_ppp",    &out_vChi2_ppp,         "vChi2_ppp");
  outTree->Branch("ctau_ppp",     &out_ctau_ppp,  "ctau_ppp");
  outTree->Branch("ctauErr_ppp",  &out_ctauErr_ppp,       "ctauErr_ppp");
  outTree->Branch("cosAlpha_ppp",         &out_cosAlpha_ppp,      "cosAlpha_ppp");
  outTree->Branch("onePsiPrime_m_ppp",    &out_onePsiPrime_m_ppp,         "onePsiPrime_m_ppp");
  outTree->Branch("twoPsiPrime_m_ppp",    &out_twoPsiPrime_m_ppp,         "twoPsiPrime_m_ppp");
  outTree->Branch("onePsiPrime_p_mppp",   &out_onePsiPrime_p_mppp,        "onePsiPrime_p_mppp");
  outTree->Branch("onePsiPrime_m_mppp",   &out_onePsiPrime_m_mppp,        "onePsiPrime_m_mppp");
  outTree->Branch("twoPsiPrime_p_mppp",   &out_twoPsiPrime_p_mppp,        "twoPsiPrime_p_mppp");
  outTree->Branch("twoPsiPrime_m_mppp",   &out_twoPsiPrime_m_mppp,        "twoPsiPrime_m_mppp");
  outTree->Branch("dimuonDiTrkOne_m_ppp",         &out_dimuonDiTrkOne_m_ppp,      "dimuonDiTrkOne_m_ppp");
  outTree->Branch("dimuonDiTrkTwo_m_ppp",         &out_dimuonDiTrkTwo_m_ppp,      "dimuonDiTrkTwo_m_ppp");
  outTree->Branch("dimuonDiTrkThree_m_ppp",       &out_dimuonDiTrkThree_m_ppp,    "dimuonDiTrkThree_m_ppp");
  outTree->Branch("ditrakOne_m_ppp",      &out_ditrakOne_m_ppp,   "ditrakOne_m_ppp");
  outTree->Branch("ditrakTwo_m_ppp",      &out_ditrakTwo_m_ppp,   "ditrakTwo_m_ppp");
  outTree->Branch("ditrakThree_m_ppp",    &out_ditrakThree_m_ppp,         "ditrakThree_m_ppp");
  outTree->Branch("trackOne_m_ppp",       &out_trackOne_m_ppp,    "trackOne_m_ppp");
  outTree->Branch("trackTwo_m_ppp",       &out_trackTwo_m_ppp,    "trackTwo_m_ppp");
  outTree->Branch("trackThree_m_ppp",     &out_trackThree_m_ppp,  "trackThree_m_ppp");


}

Bool_t FiveTracks::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  fReader.SetEntry(entry);

  bool test = true;

  test = test && (*lowMuonMatch>0.0) && (*highMuonMatch>0.0);

  test = test && (*vProb_ppk> 0.01 || *vProb_kkk > 0.01 || *vProb_ppp > 0.01 || *vProb_kpp > 0.01 || *vProb_pkp > 0.01);
  test = test && (*dimuon_pt> 3.0); //&& (*dimuonDiTrkOne_pt>5.0) ; 
  //int a = (int) (*trigger);
  //std::cout << (*trigger);
  //

  if(test)
  {
    out_run =       (Float_t)(*run);
    out_event =     (Float_t)(*event);
    out_lumi =      (Float_t)(*lumi);
    out_numPrimaryVertices =        (Float_t)(*numPrimaryVertices);
    out_trigger =   (Float_t)(*trigger);
    out_noFiveCandidates =  (Float_t)(*noFiveCandidates);
    out_dimuonditrk_id =    (Float_t)(*dimuonditrk_id);
    out_dimuonditrk_m =     (Float_t)(*dimuonditrk_m);
    out_dimuonditrk_pt =    (Float_t)(*dimuonditrk_pt);
    out_dimuonditrk_eta =   (Float_t)(*dimuonditrk_eta);
    out_dimuonditrk_phi =   (Float_t)(*dimuonditrk_phi);
    out_dimuonditrk_p =     (Float_t)(*dimuonditrk_p);
    out_dimuon_m =  (Float_t)(*dimuon_m);
    out_dimuon_pt =         (Float_t)(*dimuon_pt);
    out_dimuon_eta =        (Float_t)(*dimuon_eta);
    out_dimuon_phi =        (Float_t)(*dimuon_phi);
    out_dimuon_p =  (Float_t)(*dimuon_p);
    out_highKaonMatch =     (Float_t)(*highKaonMatch);
    out_lowKaonMatch =      (Float_t)(*lowKaonMatch);
    out_lowMuonMatch =      (Float_t)(*lowMuonMatch);
    out_highMuonMatch =     (Float_t)(*highMuonMatch);
    out_ditrak_m =  (Float_t)(*ditrak_m);
    out_ditrakOne_pt =      (Float_t)(*ditrakOne_pt);
    out_ditrakOne_eta =     (Float_t)(*ditrakOne_eta);
    out_ditrakOne_phi =     (Float_t)(*ditrakOne_phi);
    out_ditrakOne_p =       (Float_t)(*ditrakOne_p);
    out_ditrakTwo_pt =      (Float_t)(*ditrakTwo_pt);
    out_ditrakTwo_eta =     (Float_t)(*ditrakTwo_eta);
    out_ditrakTwo_phi =     (Float_t)(*ditrakTwo_phi);
    out_ditrakTwo_p =       (Float_t)(*ditrakTwo_p);
    out_ditrakThree_pt =    (Float_t)(*ditrakThree_pt);
    out_ditrakThree_eta =   (Float_t)(*ditrakThree_eta);
    out_ditrakThree_phi =   (Float_t)(*ditrakThree_phi);
    out_ditrakThree_p =     (Float_t)(*ditrakThree_p);
    out_highTrack_pt =      (Float_t)(*highTrack_pt);
    out_highTrack_eta =     (Float_t)(*highTrack_eta);
    out_highTrack_phi =     (Float_t)(*highTrack_phi);
    out_highTrack_charge =  (Float_t)(*highTrack_charge);
    out_lowTrack_pt =       (Float_t)(*lowTrack_pt);
    out_lowTrack_eta =      (Float_t)(*lowTrack_eta);
    out_lowTrack_phi =      (Float_t)(*lowTrack_phi);
    out_lowTrack_charge =   (Float_t)(*lowTrack_charge);
    out_thirdTrack_pt =     (Float_t)(*thirdTrack_pt);
    out_thirdTrack_eta =    (Float_t)(*thirdTrack_eta);
    out_thirdTrack_phi =    (Float_t)(*thirdTrack_phi);
    out_thirdTrack_charge =         (Float_t)(*thirdTrack_charge);
    out_thirdTrack_dz =     (Float_t)(*thirdTrack_dz);
    out_thirdTrack_dxy =    (Float_t)(*thirdTrack_dxy);
    out_dimuonDiTrkOne_pt =         (Float_t)(*dimuonDiTrkOne_pt);
    out_dimuonDiTrkOne_eta =        (Float_t)(*dimuonDiTrkOne_eta);
    out_dimuonDiTrkOne_phi =        (Float_t)(*dimuonDiTrkOne_phi);
    out_dimuonDiTrkOne_charge =     (Float_t)(*dimuonDiTrkOne_charge);
    out_dimuonDiTrkTwo_pt =         (Float_t)(*dimuonDiTrkTwo_pt);
    out_dimuonDiTrkTwo_eta =        (Float_t)(*dimuonDiTrkTwo_eta);
    out_dimuonDiTrkTwo_phi =        (Float_t)(*dimuonDiTrkTwo_phi);
    out_dimuonDiTrkTwo_charge =     (Float_t)(*dimuonDiTrkTwo_charge);
    out_dimuonDiTrkThree_pt =       (Float_t)(*dimuonDiTrkThree_pt);
    out_dimuonDiTrkThree_eta =      (Float_t)(*dimuonDiTrkThree_eta);
    out_dimuonDiTrkThree_phi =      (Float_t)(*dimuonDiTrkThree_phi);
    out_dimuonDiTrkThree_charge =   (Float_t)(*dimuonDiTrkThree_charge);
    out_psiPrimeSame_pt =   (Float_t)(*psiPrimeSame_pt);
    out_psiPrimeSame_eta =  (Float_t)(*psiPrimeSame_eta);
    out_psiPrimeSame_phi =  (Float_t)(*psiPrimeSame_phi);
    out_psiPrimeSame_n =    (Float_t)(*psiPrimeSame_n);
    out_psiPrimeSame_p_pt =         (Float_t)(*psiPrimeSame_p_pt);
    out_psiPrimeSame_p_eta =        (Float_t)(*psiPrimeSame_p_eta);
    out_psiPrimeSame_p_phi =        (Float_t)(*psiPrimeSame_p_phi);
    out_psiPrimeSame_p_n =  (Float_t)(*psiPrimeSame_p_n);
    out_psiPrimeSame_m_pt =         (Float_t)(*psiPrimeSame_m_pt);
    out_psiPrimeSame_m_eta =        (Float_t)(*psiPrimeSame_m_eta);
    out_psiPrimeSame_m_phi =        (Float_t)(*psiPrimeSame_m_phi);
    out_psiPrimeSame_m_n =  (Float_t)(*psiPrimeSame_m_n);
    out_psiPrimeMixed_pt =  (Float_t)(*psiPrimeMixed_pt);
    out_psiPrimeMixed_eta =         (Float_t)(*psiPrimeMixed_eta);
    out_psiPrimeMixed_phi =         (Float_t)(*psiPrimeMixed_phi);
    out_psiPrimeMixed_n =   (Float_t)(*psiPrimeMixed_n);
    out_psiPrimeMixed_p_pt =        (Float_t)(*psiPrimeMixed_p_pt);
    out_psiPrimeMixed_p_eta =       (Float_t)(*psiPrimeMixed_p_eta);
    out_psiPrimeMixed_p_phi =       (Float_t)(*psiPrimeMixed_p_phi);
    out_psiPrimeMixed_p_n =         (Float_t)(*psiPrimeMixed_p_n);
    out_psiPrimeMixed_m_pt =        (Float_t)(*psiPrimeMixed_m_pt);
    out_psiPrimeMixed_m_eta =       (Float_t)(*psiPrimeMixed_m_eta);
    out_psiPrimeMixed_m_phi =       (Float_t)(*psiPrimeMixed_m_phi);
    out_psiPrimeMixed_m_n =         (Float_t)(*psiPrimeMixed_m_n);
    out_psiPrimeSame_ditrak_pt =    (Float_t)(*psiPrimeSame_ditrak_pt);
    out_psiPrimeSame_ditrak_eta =   (Float_t)(*psiPrimeSame_ditrak_eta);
    out_psiPrimeSame_ditrak_phi =   (Float_t)(*psiPrimeSame_ditrak_phi);
    out_psiPrimeSame_ditrak_n =     (Float_t)(*psiPrimeSame_ditrak_n);
    out_psiPrimeMixed_ditrak_pt =   (Float_t)(*psiPrimeMixed_ditrak_pt);
    out_psiPrimeMixed_ditrak_eta =  (Float_t)(*psiPrimeMixed_ditrak_eta);
    out_psiPrimeMixed_ditrak_phi =  (Float_t)(*psiPrimeMixed_ditrak_phi);
    out_psiPrimeMixed_ditrak_n =    (Float_t)(*psiPrimeMixed_ditrak_n);
    out_triTrak_pt =        (Float_t)(*triTrak_pt);
    out_triTrak_eta =       (Float_t)(*triTrak_eta);
    out_triTrak_phi =       (Float_t)(*triTrak_phi);
    out_triTrak_charge =    (Float_t)(*triTrak_charge);
    out_mass_kkk =  (Float_t)(*mass_kkk);
    out_mass_ref_kkk =      (Float_t)(*mass_ref_kkk);
    out_vProb_kkk =         (Float_t)(*vProb_kkk);
    out_nDof_kkk =  (Float_t)(*nDof_kkk);
    out_vChi2_kkk =         (Float_t)(*vChi2_kkk);
    out_ctau_kkk =  (Float_t)(*ctau_kkk);
    out_ctauErr_kkk =       (Float_t)(*ctauErr_kkk);
    out_cosAlpha_kkk =      (Float_t)(*cosAlpha_kkk);
    out_onePsiPrime_m_kkk =         (Float_t)(*onePsiPrime_m_kkk);
    out_twoPsiPrime_m_kkk =         (Float_t)(*twoPsiPrime_m_kkk);
    out_onePsiPrime_p_mkkk =        (Float_t)(*onePsiPrime_p_mkkk);
    out_onePsiPrime_m_mkkk =        (Float_t)(*onePsiPrime_m_mkkk);
    out_twoPsiPrime_p_mkkk =        (Float_t)(*twoPsiPrime_p_mkkk);
    out_twoPsiPrime_m_mkkk =        (Float_t)(*twoPsiPrime_m_mkkk);
    out_dimuonDiTrkOne_m_kkk =      (Float_t)(*dimuonDiTrkOne_m_kkk);
    out_dimuonDiTrkTwo_m_kkk =      (Float_t)(*dimuonDiTrkTwo_m_kkk);
    out_dimuonDiTrkThree_m_kkk =    (Float_t)(*dimuonDiTrkThree_m_kkk);
    out_ditrakOne_m_kkk =   (Float_t)(*ditrakOne_m_kkk);
    out_ditrakTwo_m_kkk =   (Float_t)(*ditrakTwo_m_kkk);
    out_ditrakThree_m_kkk =         (Float_t)(*ditrakThree_m_kkk);
    out_trackOne_m_kkk =    (Float_t)(*trackOne_m_kkk);
    out_trackTwo_m_kkk =    (Float_t)(*trackTwo_m_kkk);
    out_trackThree_m_kkk =  (Float_t)(*trackThree_m_kkk);

    out_triTrak_m_kkk =     (Float_t)(*triTrak_m_kkk);
    out_mass_ppk =  (Float_t)(*mass_ppk);
    out_mass_ref_ppk =      (Float_t)(*mass_ref_ppk);
    out_vProb_ppk =         (Float_t)(*vProb_ppk);
    out_nDof_ppk =  (Float_t)(*nDof_ppk);
    out_vChi2_ppk =         (Float_t)(*vChi2_ppk);
    out_ctau_ppk =  (Float_t)(*ctau_ppk);
    out_ctauErr_ppk =       (Float_t)(*ctauErr_ppk);
    out_cosAlpha_ppk =      (Float_t)(*cosAlpha_ppk);
    out_onePsiPrime_m_ppk =         (Float_t)(*onePsiPrime_m_ppk);
    out_twoPsiPrime_m_ppk =         (Float_t)(*twoPsiPrime_m_ppk);
    out_onePsiPrime_p_mppk =        (Float_t)(*onePsiPrime_p_mppk);
    out_onePsiPrime_m_mppk =        (Float_t)(*onePsiPrime_m_mppk);
    out_twoPsiPrime_p_mppk =        (Float_t)(*twoPsiPrime_p_mppk);
    out_twoPsiPrime_m_mppk =        (Float_t)(*twoPsiPrime_m_mppk);
    out_dimuonDiTrkOne_m_ppk =      (Float_t)(*dimuonDiTrkOne_m_ppk);
    out_dimuonDiTrkTwo_m_ppk =      (Float_t)(*dimuonDiTrkTwo_m_ppk);
    out_dimuonDiTrkThree_m_ppk =    (Float_t)(*dimuonDiTrkThree_m_ppk);
    out_ditrakOne_m_ppk =   (Float_t)(*ditrakOne_m_ppk);
    out_ditrakTwo_m_ppk =   (Float_t)(*ditrakTwo_m_ppk);
    out_ditrakThree_m_ppk =         (Float_t)(*ditrakThree_m_ppk);
    out_trackOne_m_ppk =    (Float_t)(*trackOne_m_ppk);
    out_trackTwo_m_ppk =    (Float_t)(*trackTwo_m_ppk);
    out_trackThree_m_ppk =  (Float_t)(*trackThree_m_ppk);

    out_triTrak_m_ppk =     (Float_t)(*triTrak_m_ppk);
    out_mass_kpp =  (Float_t)(*mass_kpp);
    out_mass_ref_kpp =      (Float_t)(*mass_ref_kpp);
    out_vProb_kpp =         (Float_t)(*vProb_kpp);
    out_nDof_kpp =  (Float_t)(*nDof_kpp);
    out_vChi2_kpp =         (Float_t)(*vChi2_kpp);
    out_ctau_kpp =  (Float_t)(*ctau_kpp);
    out_ctauErr_kpp =       (Float_t)(*ctauErr_kpp);
    out_cosAlpha_kpp =      (Float_t)(*cosAlpha_kpp);
    out_onePsiPrime_m_kpp =         (Float_t)(*onePsiPrime_m_kpp);
    out_twoPsiPrime_m_kpp =         (Float_t)(*twoPsiPrime_m_kpp);
    out_onePsiPrime_p_mkpp =        (Float_t)(*onePsiPrime_p_mkpp);
    out_onePsiPrime_m_mkpp =        (Float_t)(*onePsiPrime_m_mkpp);
    out_twoPsiPrime_p_mkpp =        (Float_t)(*twoPsiPrime_p_mkpp);
    out_twoPsiPrime_m_mkpp =        (Float_t)(*twoPsiPrime_m_mkpp);
    out_dimuonDiTrkOne_m_kpp =      (Float_t)(*dimuonDiTrkOne_m_kpp);
    out_dimuonDiTrkTwo_m_kpp =      (Float_t)(*dimuonDiTrkTwo_m_kpp);
    out_dimuonDiTrkThree_m_kpp =    (Float_t)(*dimuonDiTrkThree_m_kpp);
    out_ditrakOne_m_kpp =   (Float_t)(*ditrakOne_m_kpp);
    out_ditrakTwo_m_kpp =   (Float_t)(*ditrakTwo_m_kpp);
    out_ditrakThree_m_kpp =         (Float_t)(*ditrakThree_m_kpp);
    out_trackOne_m_kpp =    (Float_t)(*trackOne_m_kpp);
    out_trackTwo_m_kpp =    (Float_t)(*trackTwo_m_kpp);
    out_trackThree_m_kpp =  (Float_t)(*trackThree_m_kpp);

    out_triTrak_m_kpp =     (Float_t)(*triTrak_m_kpp);
    out_mass_pkp =  (Float_t)(*mass_pkp);
    out_mass_ref_pkp =      (Float_t)(*mass_ref_pkp);
    out_vProb_pkp =         (Float_t)(*vProb_pkp);
    out_nDof_pkp =  (Float_t)(*nDof_pkp);
    out_vChi2_pkp =         (Float_t)(*vChi2_pkp);
    out_ctau_pkp =  (Float_t)(*ctau_pkp);
    out_ctauErr_pkp =       (Float_t)(*ctauErr_pkp);
    out_cosAlpha_pkp =      (Float_t)(*cosAlpha_pkp);
    out_onePsiPrime_m_pkp =         (Float_t)(*onePsiPrime_m_pkp);
    out_twoPsiPrime_m_pkp =         (Float_t)(*twoPsiPrime_m_pkp);
    out_onePsiPrime_p_mpkp =        (Float_t)(*onePsiPrime_p_mpkp);
    out_onePsiPrime_m_mpkp =        (Float_t)(*onePsiPrime_m_mpkp);
    out_twoPsiPrime_p_mpkp =        (Float_t)(*twoPsiPrime_p_mpkp);
    out_twoPsiPrime_m_mpkp =        (Float_t)(*twoPsiPrime_m_mpkp);
    out_dimuonDiTrkOne_m_pkp =      (Float_t)(*dimuonDiTrkOne_m_pkp);
    out_dimuonDiTrkTwo_m_pkp =      (Float_t)(*dimuonDiTrkTwo_m_pkp);
    out_dimuonDiTrkThree_m_pkp =    (Float_t)(*dimuonDiTrkThree_m_pkp);
    out_ditrakOne_m_pkp =   (Float_t)(*ditrakOne_m_pkp);
    out_ditrakTwo_m_pkp =   (Float_t)(*ditrakTwo_m_pkp);
    out_ditrakThree_m_pkp =         (Float_t)(*ditrakThree_m_pkp);
    out_trackOne_m_pkp =    (Float_t)(*trackOne_m_pkp);
    out_trackTwo_m_pkp =    (Float_t)(*trackTwo_m_pkp);
    out_trackThree_m_pkp =  (Float_t)(*trackThree_m_pkp);

    out_triTrak_m_pkp =     (Float_t)(*triTrak_m_pkp);
    out_mass_ppp =  (Float_t)(*mass_ppp);
    out_mass_ref_ppp =      (Float_t)(*mass_ref_ppp);
    out_vProb_ppp =         (Float_t)(*vProb_ppp);
    out_nDof_ppp =  (Float_t)(*nDof_ppp);
    out_vChi2_ppp =         (Float_t)(*vChi2_ppp);
    out_ctau_ppp =  (Float_t)(*ctau_ppp);
    out_ctauErr_ppp =       (Float_t)(*ctauErr_ppp);
    out_cosAlpha_ppp =      (Float_t)(*cosAlpha_ppp);
    out_onePsiPrime_m_ppp =         (Float_t)(*onePsiPrime_m_ppp);
    out_twoPsiPrime_m_ppp =         (Float_t)(*twoPsiPrime_m_ppp);
    out_onePsiPrime_p_mppp =        (Float_t)(*onePsiPrime_p_mppp);
    out_onePsiPrime_m_mppp =        (Float_t)(*onePsiPrime_m_mppp);
    out_twoPsiPrime_p_mppp =        (Float_t)(*twoPsiPrime_p_mppp);
    out_twoPsiPrime_m_mppp =        (Float_t)(*twoPsiPrime_m_mppp);
    out_dimuonDiTrkOne_m_ppp =      (Float_t)(*dimuonDiTrkOne_m_ppp);
    out_dimuonDiTrkTwo_m_ppp =      (Float_t)(*dimuonDiTrkTwo_m_ppp);
    out_dimuonDiTrkThree_m_ppp =    (Float_t)(*dimuonDiTrkThree_m_ppp);
    out_ditrakOne_m_ppp =   (Float_t)(*ditrakOne_m_ppp);
    out_ditrakTwo_m_ppp =   (Float_t)(*ditrakTwo_m_ppp);
    out_ditrakThree_m_ppp =         (Float_t)(*ditrakThree_m_ppp);
    out_trackOne_m_ppp =    (Float_t)(*trackOne_m_ppp);
    out_trackTwo_m_ppp =    (Float_t)(*trackTwo_m_ppp);
    out_trackThree_m_ppp =  (Float_t)(*trackThree_m_ppp);


    outTree->Fill();
  }

  return kTRUE;
}

void FiveTracks::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  TDirectory *savedir = gDirectory;
  if (fOut)
  {
    fOut->cd();
    gStyle->SetOptStat(111111) ;


    outTree->Write();
    OutFile->Print();
    fOutput->Add(OutFile);
    gDirectory = savedir;
    fOut->Close();

  }

}

void FiveTracks::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

}
