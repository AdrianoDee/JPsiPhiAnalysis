#define FourMuons_cxx
// The class definition in FourMuons.h has been generated automatically
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
// root> T->Process("FourMuons.C")
// root> T->Process("FourMuons.C","some options")
// root> T->Process("FourMuons.C+")
//


#include "FourMuons.h"
#include <TH2.h>
#include <TStyle.h>

void FourMuons::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void FourMuons::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   std::string outputString = "2mu2mu_tree.root";
   OutFile = new TProofOutputFile( outputString.data() );
   fOut = OutFile->OpenFile("RECREATE/F");
   if (!(fOut=OutFile->OpenFile("RECREATE")))
   {
     Warning("SlaveBegin","Problems opening file: %s%s", OutFile->GetDir(), OutFile->GetFileName() );
   }

   ////////////////// Histograms //////////////////
   JPsi_mass = 3.096916; /// pdg mass
   Phi_mass = 1.019455; /// pdg mass
   Phi_mean = 1.019723;
   Phi_sigma = 2.35607e-03;//2.28400e-03;

   outTree = new TTree("2mu2muSkimmedTree","2mu2muSkimmedTree/F");

    outTree->Branch("run",  &out_run,       "run/F");
    outTree->Branch("event",        &out_event,     "event/F");
    outTree->Branch("numPrimaryVertices",   &out_numPrimaryVertices,        "numPrimaryVertices/F");
    outTree->Branch("trigger",      &out_trigger,   "trigger/F");
    outTree->Branch("noXCandidates",        &out_noXCandidates,     "noXCandidates/F");

    outTree->Branch("jpsi_triggerMatch",    &out_jpsi_triggerMatch,         "jpsi_triggerMatch/F");
    outTree->Branch("phi_triggerMatch",     &out_phi_triggerMatch,  "phi_triggerMatch/F");
    outTree->Branch("mHighJPsiMatch",       &out_mHighJPsiMatch,    "mHighJPsiMatch/F");
    outTree->Branch("mLowJPsiMatch",        &out_mLowJPsiMatch,     "mLowJPsiMatch/F");
    outTree->Branch("mHighPhiMatch",        &out_mHighPhiMatch,     "mHighPhiMatch/F");
    outTree->Branch("mLowPhiMatch",         &out_mLowPhiMatch,      "mLowPhiMatch/F");
    outTree->Branch("doubledimuon_charge",  &out_doubledimuon_charge,       "doubledimuon_charge/F");
    outTree->Branch("doubledimuon_m",       &out_doubledimuon_m,    "doubledimuon_m/F");
    outTree->Branch("doubledimuon_m_rf",    &out_doubledimuon_m_rf,         "doubledimuon_m_rf/F");
    outTree->Branch("doubledimuon_m_rf_c",  &out_doubledimuon_m_rf_c,       "doubledimuon_m_rf_c/F");
    outTree->Branch("doubledimuon_m_rf_d_c",        &out_doubledimuon_m_rf_d_c,     "doubledimuon_m_rf_d_c/F");
    outTree->Branch("doubledimuon_p",       &out_doubledimuon_p,    "doubledimuon_p/F");
    outTree->Branch("doubledimuon_pt",      &out_doubledimuon_pt,   "doubledimuon_pt/F");
    outTree->Branch("doubledimuon_e",       &out_doubledimuon_e,    "doubledimuon_e/F");
    outTree->Branch("doubledimuon_eta",     &out_doubledimuon_eta,  "doubledimuon_eta/F");
    outTree->Branch("doubledimuon_theta",   &out_doubledimuon_theta,        "doubledimuon_theta/F");
    outTree->Branch("doubledimuon_y",       &out_doubledimuon_y,    "doubledimuon_y/F");
    outTree->Branch("doubledimuon_dxy",     &out_doubledimuon_dxy,  "doubledimuon_dxy/F");
    outTree->Branch("doubledimuon_dxyErr",  &out_doubledimuon_dxyErr,       "doubledimuon_dxyErr/F");
    outTree->Branch("doubledimuon_dz",      &out_doubledimuon_dz,   "doubledimuon_dz/F");
    outTree->Branch("doubledimuon_dzErr",   &out_doubledimuon_dzErr,        "doubledimuon_dzErr/F");
    outTree->Branch("doubledimuon_vProb",   &out_doubledimuon_vProb,        "doubledimuon_vProb/F");
    outTree->Branch("doubledimuon_vChi2",   &out_doubledimuon_vChi2,        "doubledimuon_vChi2/F");
    outTree->Branch("doubledimuon_nDof",    &out_doubledimuon_nDof,         "doubledimuon_nDof/F");
    outTree->Branch("doubledimuon_rf_vProb",        &out_doubledimuon_rf_vProb,     "doubledimuon_rf_vProb/F");
    outTree->Branch("doubledimuon_rf_nDof",         &out_doubledimuon_rf_nDof,      "doubledimuon_rf_nDof/F");
    outTree->Branch("doubledimuon_rf_c_vProb",      &out_doubledimuon_rf_c_vProb,   "doubledimuon_rf_c_vProb/F");
    outTree->Branch("doubledimuon_rf_c_vChi2",      &out_doubledimuon_rf_c_vChi2,   "doubledimuon_rf_c_vChi2/F");
    outTree->Branch("doubledimuon_rf_c_nDof",       &out_doubledimuon_rf_c_nDof,    "doubledimuon_rf_c_nDof/F");
    outTree->Branch("doubledimuon_vx",      &out_doubledimuon_vx,   "doubledimuon_vx/F");
    outTree->Branch("doubledimuon_vy",      &out_doubledimuon_vy,   "doubledimuon_vy/F");
    outTree->Branch("doubledimuon_vz",      &out_doubledimuon_vz,   "doubledimuon_vz/F");
    outTree->Branch("pv_x",         &out_pv_x,      "pv_x/F");
    outTree->Branch("pv_y",         &out_pv_y,      "pv_y/F");
    outTree->Branch("pv_z",         &out_pv_z,      "pv_z/F");
    outTree->Branch("doubledimuon_cosAlpha",        &out_doubledimuon_cosAlpha,     "doubledimuon_cosAlpha/F");
    outTree->Branch("doubledimuon_cosAlpha3D",      &out_doubledimuon_cosAlpha3D,   "doubledimuon_cosAlpha3D/F");
    outTree->Branch("doubledimuon_ctauPV",  &out_doubledimuon_ctauPV,       "doubledimuon_ctauPV/F");
    outTree->Branch("doubledimuon_ctauErrPV",       &out_doubledimuon_ctauErrPV,    "doubledimuon_ctauErrPV/F");
    outTree->Branch("doubledimuon_lxy",     &out_doubledimuon_lxy,  "doubledimuon_lxy/F");
    outTree->Branch("doubledimuon_lxyErr",  &out_doubledimuon_lxyErr,       "doubledimuon_lxyErr/F");
    outTree->Branch("doubledimuon_lxyz",    &out_doubledimuon_lxyz,         "doubledimuon_lxyz/F");
    outTree->Branch("doubledimuon_lxyzErr",         &out_doubledimuon_lxyzErr,      "doubledimuon_lxyzErr/F");
    outTree->Branch("doubledimuon_rf_cosAlpha",     &out_doubledimuon_rf_cosAlpha,  "doubledimuon_rf_cosAlpha/F");
    outTree->Branch("doubledimuon_rf_ctauPV",       &out_doubledimuon_rf_ctauPV,    "doubledimuon_rf_ctauPV/F");
    outTree->Branch("doubledimuon_rf_ctauErrPV",    &out_doubledimuon_rf_ctauErrPV,         "doubledimuon_rf_ctauErrPV/F");
    outTree->Branch("doubledimuon_rf_lxy",  &out_doubledimuon_rf_lxy,       "doubledimuon_rf_lxy/F");
    outTree->Branch("doubledimuon_rf_lxyErr",       &out_doubledimuon_rf_lxyErr,    "doubledimuon_rf_lxyErr/F");
    outTree->Branch("doubledimuon_rf_lxyz",         &out_doubledimuon_rf_lxyz,      "doubledimuon_rf_lxyz/F");
    outTree->Branch("doubledimuon_rf_lxyzErr",      &out_doubledimuon_rf_lxyzErr,   "doubledimuon_rf_lxyzErr/F");
    outTree->Branch("doubledimuon_rf_c_cosAlpha",   &out_doubledimuon_rf_c_cosAlpha,        "doubledimuon_rf_c_cosAlpha/F");
    outTree->Branch("doubledimuon_rf_c_ctauPV",     &out_doubledimuon_rf_c_ctauPV,  "doubledimuon_rf_c_ctauPV/F");
    outTree->Branch("doubledimuon_rf_c_ctauErrPV",  &out_doubledimuon_rf_c_ctauErrPV,       "doubledimuon_rf_c_ctauErrPV/F");
    outTree->Branch("doubledimuon_rf_c_lxy",        &out_doubledimuon_rf_c_lxy,     "doubledimuon_rf_c_lxy/F");
    outTree->Branch("doubledimuon_rf_c_lxyErr",     &out_doubledimuon_rf_c_lxyErr,  "doubledimuon_rf_c_lxyErr/F");
    outTree->Branch("doubledimuon_rf_c_lxyz",       &out_doubledimuon_rf_c_lxyz,    "doubledimuon_rf_c_lxyz/F");
    outTree->Branch("doubledimuon_rf_c_lxyzErr",    &out_doubledimuon_rf_c_lxyzErr,         "doubledimuon_rf_c_lxyzErr/F");
    outTree->Branch("doubledimuon_dca_mp1mp2",      &out_doubledimuon_dca_mp1mp2,   "doubledimuon_dca_mp1mp2/F");
    outTree->Branch("doubledimuon_dca_mp1mj1",      &out_doubledimuon_dca_mp1mj1,   "doubledimuon_dca_mp1mj1/F");
    outTree->Branch("doubledimuon_dca_mp1mj2",      &out_doubledimuon_dca_mp1mj2,   "doubledimuon_dca_mp1mj2/F");
    outTree->Branch("doubledimuon_dca_mp2mj1",      &out_doubledimuon_dca_mp2mj1,   "doubledimuon_dca_mp2mj1/F");
    outTree->Branch("doubledimuon_dca_mp2mj2",      &out_doubledimuon_dca_mp2mj2,   "doubledimuon_dca_mp2mj2/F");
    outTree->Branch("doubledimuon_dca_mj1mj2",      &out_doubledimuon_dca_mj1mj2,   "doubledimuon_dca_mj1mj2/F");
    outTree->Branch("doubledimuon_cosAlphaDZ",      &out_doubledimuon_cosAlphaDZ,   "doubledimuon_cosAlphaDZ/F");
    outTree->Branch("doubledimuon_cosAlpha3DDZ",    &out_doubledimuon_cosAlpha3DDZ,         "doubledimuon_cosAlpha3DDZ/F");
    outTree->Branch("doubledimuon_ctauPVDZ",        &out_doubledimuon_ctauPVDZ,     "doubledimuon_ctauPVDZ/F");
    outTree->Branch("doubledimuon_ctauErrPVDZ",     &out_doubledimuon_ctauErrPVDZ,  "doubledimuon_ctauErrPVDZ/F");
    outTree->Branch("doubledimuon_lxyDZ",   &out_doubledimuon_lxyDZ,        "doubledimuon_lxyDZ/F");
    outTree->Branch("doubledimuon_lxyErrDZ",        &out_doubledimuon_lxyErrDZ,     "doubledimuon_lxyErrDZ/F");
    outTree->Branch("doubledimuon_lxyzDZ",  &out_doubledimuon_lxyzDZ,       "doubledimuon_lxyzDZ/F");
    outTree->Branch("doubledimuon_lxyzErrDZ",       &out_doubledimuon_lxyzErrDZ,    "doubledimuon_lxyzErrDZ/F");
    outTree->Branch("doubledimuon_cosAlphaBS",      &out_doubledimuon_cosAlphaBS,   "doubledimuon_cosAlphaBS/F");
    outTree->Branch("doubledimuon_cosAlpha3DBS",    &out_doubledimuon_cosAlpha3DBS,         "doubledimuon_cosAlpha3DBS/F");
    outTree->Branch("doubledimuon_ctauPVBS",        &out_doubledimuon_ctauPVBS,     "doubledimuon_ctauPVBS/F");
    outTree->Branch("doubledimuon_ctauErrPVBS",     &out_doubledimuon_ctauErrPVBS,  "doubledimuon_ctauErrPVBS/F");
    outTree->Branch("doubledimuon_lxyBS",   &out_doubledimuon_lxyBS,        "doubledimuon_lxyBS/F");
    outTree->Branch("doubledimuon_lxyErrBS",        &out_doubledimuon_lxyErrBS,     "doubledimuon_lxyErrBS/F");
    outTree->Branch("doubledimuon_lxyzBS",  &out_doubledimuon_lxyzBS,       "doubledimuon_lxyzBS/F");
    outTree->Branch("doubledimuon_lxyzErrBS",       &out_doubledimuon_lxyzErrBS,    "doubledimuon_lxyzErrBS/F");
    outTree->Branch("mHighPhi_p",   &out_mHighPhi_p,        "mHighPhi_p/F");
    outTree->Branch("mLowPhi_p",    &out_mLowPhi_p,         "mLowPhi_p/F");
    outTree->Branch("mHighJPsi_p",  &out_mHighJPsi_p,       "mHighJPsi_p/F");
    outTree->Branch("mLowJPsi_p",   &out_mLowJPsi_p,        "mLowJPsi_p/F");
    outTree->Branch("mHighPhi_pt",  &out_mHighPhi_pt,       "mHighPhi_pt/F");
    outTree->Branch("mLowPhi_pt",   &out_mLowPhi_pt,        "mLowPhi_pt/F");
    outTree->Branch("mHighJPsi_pt",         &out_mHighJPsi_pt,      "mHighJPsi_pt/F");
    outTree->Branch("mLowJPsi_pt",  &out_mLowJPsi_pt,       "mLowJPsi_pt/F");
    outTree->Branch("mHighPhi_ptErr",       &out_mHighPhi_ptErr,    "mHighPhi_ptErr/F");
    outTree->Branch("mLowPhi_ptErr",        &out_mLowPhi_ptErr,     "mLowPhi_ptErr/F");
    outTree->Branch("mHighJPsi_ptErr",      &out_mHighJPsi_ptErr,   "mHighJPsi_ptErr/F");
    outTree->Branch("mLowJPsi_ptErr",       &out_mLowJPsi_ptErr,    "mLowJPsi_ptErr/F");
    outTree->Branch("mHighJPsi_eta",        &out_mHighJPsi_eta,     "mHighJPsi_eta/F");
    outTree->Branch("mLowJPsi_eta",         &out_mLowJPsi_eta,      "mLowJPsi_eta/F");
    outTree->Branch("mHighPhi_eta",         &out_mHighPhi_eta,      "mHighPhi_eta/F");
    outTree->Branch("mLowPhi_eta",  &out_mLowPhi_eta,       "mLowPhi_eta/F");
    outTree->Branch("mHighJPsi_etaErr",     &out_mHighJPsi_etaErr,  "mHighJPsi_etaErr/F");
    outTree->Branch("mLowJPsi_etaErr",      &out_mLowJPsi_etaErr,   "mLowJPsi_etaErr/F");
    outTree->Branch("mHighPhi_etaErr",      &out_mHighPhi_etaErr,   "mHighPhi_etaErr/F");
    outTree->Branch("mLowPhi_etaErr",       &out_mLowPhi_etaErr,    "mLowPhi_etaErr/F");
    outTree->Branch("mHighJPsi_phi",        &out_mHighJPsi_phi,     "mHighJPsi_phi/F");
    outTree->Branch("mLowJPsi_phi",         &out_mLowJPsi_phi,      "mLowJPsi_phi/F");
    outTree->Branch("mHighPhi_phi",         &out_mHighPhi_phi,      "mHighPhi_phi/F");
    outTree->Branch("mLowPhi_phi",  &out_mLowPhi_phi,       "mLowPhi_phi/F");
    outTree->Branch("mHighJPsi_phiErr",     &out_mHighJPsi_phiErr,  "mHighJPsi_phiErr/F");
    outTree->Branch("mLowJPsi_phiErr",      &out_mLowJPsi_phiErr,   "mLowJPsi_phiErr/F");
    outTree->Branch("mHighPhi_phiErr",      &out_mHighPhi_phiErr,   "mHighPhi_phiErr/F");
    outTree->Branch("mLowPhi_phiErr",       &out_mLowPhi_phiErr,    "mLowPhi_phiErr/F");
    outTree->Branch("mHighJPsi_theta",      &out_mHighJPsi_theta,   "mHighJPsi_theta/F");
    outTree->Branch("mLowJPsi_theta",       &out_mLowJPsi_theta,    "mLowJPsi_theta/F");
    outTree->Branch("mHighPhi_theta",       &out_mHighPhi_theta,    "mHighPhi_theta/F");
    outTree->Branch("mLowPhi_theta",        &out_mLowPhi_theta,     "mLowPhi_theta/F");
    outTree->Branch("mHighJPsi_thetaErr",   &out_mHighJPsi_thetaErr,        "mHighJPsi_thetaErr/F");
    outTree->Branch("mLowJPsi_thetaErr",    &out_mLowJPsi_thetaErr,         "mLowJPsi_thetaErr/F");
    outTree->Branch("mHighPhi_thetaErr",    &out_mHighPhi_thetaErr,         "mHighPhi_thetaErr/F");
    outTree->Branch("mLowPhi_thetaErr",     &out_mLowPhi_thetaErr,  "mLowPhi_thetaErr/F");
    outTree->Branch("mHighJPsi_lambda",     &out_mHighJPsi_lambda,  "mHighJPsi_lambda/F");
    outTree->Branch("mLowJPsi_lambda",      &out_mLowJPsi_lambda,   "mLowJPsi_lambda/F");
    outTree->Branch("mHighPhi_lambda",      &out_mHighPhi_lambda,   "mHighPhi_lambda/F");
    outTree->Branch("mLowPhi_lambda",       &out_mLowPhi_lambda,    "mLowPhi_lambda/F");
    outTree->Branch("mHighJPsi_lambdaErr",  &out_mHighJPsi_lambdaErr,       "mHighJPsi_lambdaErr/F");
    outTree->Branch("mLowJPsi_lambdaErr",   &out_mLowJPsi_lambdaErr,        "mLowJPsi_lambdaErr/F");
    outTree->Branch("mHighPhi_lambdaErr",   &out_mHighPhi_lambdaErr,        "mHighPhi_lambdaErr/F");
    outTree->Branch("mLowPhi_lambdaErr",    &out_mLowPhi_lambdaErr,         "mLowPhi_lambdaErr/F");
    outTree->Branch("mLowPhi_dxy",  &out_mLowPhi_dxy,       "mLowPhi_dxy/F");
    outTree->Branch("mLowPhi_dxyErr",       &out_mLowPhi_dxyErr,    "mLowPhi_dxyErr/F");
    outTree->Branch("mLowPhi_dz",   &out_mLowPhi_dz,        "mLowPhi_dz/F");
    outTree->Branch("mLowPhi_dzErr",        &out_mLowPhi_dzErr,     "mLowPhi_dzErr/F");
    outTree->Branch("mHighPhi_dxy",         &out_mHighPhi_dxy,      "mHighPhi_dxy/F");
    outTree->Branch("mHighPhi_dxyErr",      &out_mHighPhi_dxyErr,   "mHighPhi_dxyErr/F");
    outTree->Branch("mHighPhi_dz",  &out_mHighPhi_dz,       "mHighPhi_dz/F");
    outTree->Branch("mHighPhi_dzErr",       &out_mHighPhi_dzErr,    "mHighPhi_dzErr/F");
    outTree->Branch("mHighJPsi_dxy",        &out_mHighJPsi_dxy,     "mHighJPsi_dxy/F");
    outTree->Branch("mHighJPsi_dxyErr",     &out_mHighJPsi_dxyErr,  "mHighJPsi_dxyErr/F");
    outTree->Branch("mHighJPsi_dz",         &out_mHighJPsi_dz,      "mHighJPsi_dz/F");
    outTree->Branch("mHighJPsi_dzErr",      &out_mHighJPsi_dzErr,   "mHighJPsi_dzErr/F");
    outTree->Branch("mLowJPsi_dxy",         &out_mLowJPsi_dxy,      "mLowJPsi_dxy/F");
    outTree->Branch("mLowJPsi_dxyErr",      &out_mLowJPsi_dxyErr,   "mLowJPsi_dxyErr/F");
    outTree->Branch("mLowJPsi_dz",  &out_mLowJPsi_dz,       "mLowJPsi_dz/F");
    outTree->Branch("mLowJPsi_dzErr",       &out_mLowJPsi_dzErr,    "mLowJPsi_dzErr/F");
    outTree->Branch("mHighJPsi_NPixelHits",         &out_mHighJPsi_NPixelHits,      "mHighJPsi_NPixelHits/F");
    outTree->Branch("mHighJPsi_NStripHits",         &out_mHighJPsi_NStripHits,      "mHighJPsi_NStripHits/F");
    outTree->Branch("mHighJPsi_NTrackhits",         &out_mHighJPsi_NTrackhits,      "mHighJPsi_NTrackhits/F");
    outTree->Branch("mHighJPsi_NBPixHits",  &out_mHighJPsi_NBPixHits,       "mHighJPsi_NBPixHits/F");
    outTree->Branch("mHighJPsi_NPixLayers",         &out_mHighJPsi_NPixLayers,      "mHighJPsi_NPixLayers/F");
    outTree->Branch("mHighJPsi_NTraLayers",         &out_mHighJPsi_NTraLayers,      "mHighJPsi_NTraLayers/F");
    outTree->Branch("mHighJPsi_NStrLayers",         &out_mHighJPsi_NStrLayers,      "mHighJPsi_NStrLayers/F");
    outTree->Branch("mHighJPsi_NBPixLayers",        &out_mHighJPsi_NBPixLayers,     "mHighJPsi_NBPixLayers/F");
    outTree->Branch("mLowJPsi_NPixelHits",  &out_mLowJPsi_NPixelHits,       "mLowJPsi_NPixelHits/F");
    outTree->Branch("mLowJPsi_NStripHits",  &out_mLowJPsi_NStripHits,       "mLowJPsi_NStripHits/F");
    outTree->Branch("mLowJPsi_NTrackhits",  &out_mLowJPsi_NTrackhits,       "mLowJPsi_NTrackhits/F");
    outTree->Branch("mLowJPsi_NBPixHits",   &out_mLowJPsi_NBPixHits,        "mLowJPsi_NBPixHits/F");
    outTree->Branch("mLowJPsi_NPixLayers",  &out_mLowJPsi_NPixLayers,       "mLowJPsi_NPixLayers/F");
    outTree->Branch("mLowJPsi_NTraLayers",  &out_mLowJPsi_NTraLayers,       "mLowJPsi_NTraLayers/F");
    outTree->Branch("mLowJPsi_NStrLayers",  &out_mLowJPsi_NStrLayers,       "mLowJPsi_NStrLayers/F");
    outTree->Branch("mLowJPsi_NBPixLayers",         &out_mLowJPsi_NBPixLayers,      "mLowJPsi_NBPixLayers/F");
    outTree->Branch("mHighPhi_NPixelHits",  &out_mHighPhi_NPixelHits,       "mHighPhi_NPixelHits/F");
    outTree->Branch("mHighPhi_NStripHits",  &out_mHighPhi_NStripHits,       "mHighPhi_NStripHits/F");
    outTree->Branch("mHighPhi_NTrackhits",  &out_mHighPhi_NTrackhits,       "mHighPhi_NTrackhits/F");
    outTree->Branch("mHighPhi_NBPixHits",   &out_mHighPhi_NBPixHits,        "mHighPhi_NBPixHits/F");
    outTree->Branch("mHighPhi_NPixLayers",  &out_mHighPhi_NPixLayers,       "mHighPhi_NPixLayers/F");
    outTree->Branch("mHighPhi_NTraLayers",  &out_mHighPhi_NTraLayers,       "mHighPhi_NTraLayers/F");
    outTree->Branch("mHighPhi_NStrLayers",  &out_mHighPhi_NStrLayers,       "mHighPhi_NStrLayers/F");
    outTree->Branch("mHighPhi_NBPixLayers",         &out_mHighPhi_NBPixLayers,      "mHighPhi_NBPixLayers/F");
    outTree->Branch("mLowPhi_NPixelHits",   &out_mLowPhi_NPixelHits,        "mLowPhi_NPixelHits/F");
    outTree->Branch("mLowPhi_NStripHits",   &out_mLowPhi_NStripHits,        "mLowPhi_NStripHits/F");
    outTree->Branch("mLowPhi_NTrackhits",   &out_mLowPhi_NTrackhits,        "mLowPhi_NTrackhits/F");
    outTree->Branch("mLowPhi_NBPixHits",    &out_mLowPhi_NBPixHits,         "mLowPhi_NBPixHits/F");
    outTree->Branch("mLowPhi_NPixLayers",   &out_mLowPhi_NPixLayers,        "mLowPhi_NPixLayers/F");
    outTree->Branch("mLowPhi_NTraLayers",   &out_mLowPhi_NTraLayers,        "mLowPhi_NTraLayers/F");
    outTree->Branch("mLowPhi_NStrLayers",   &out_mLowPhi_NStrLayers,        "mLowPhi_NStrLayers/F");
    outTree->Branch("mLowPhi_NBPixLayers",  &out_mLowPhi_NBPixLayers,       "mLowPhi_NBPixLayers/F");
    outTree->Branch("mHighJPsi_isLoose",    &out_mHighJPsi_isLoose,         "mHighJPsi_isLoose/F");
    outTree->Branch("mHighJPsi_isSoft",     &out_mHighJPsi_isSoft,  "mHighJPsi_isSoft/F");
    outTree->Branch("mHighJPsi_isMedium",   &out_mHighJPsi_isMedium,        "mHighJPsi_isMedium/F");
    outTree->Branch("mHighJPsi_isHighPt",   &out_mHighJPsi_isHighPt,        "mHighJPsi_isHighPt/F");
    outTree->Branch("mLowJPsi_isLoose",     &out_mLowJPsi_isLoose,  "mLowJPsi_isLoose/F");
    outTree->Branch("mLowJPsi_isSoft",      &out_mLowJPsi_isSoft,   "mLowJPsi_isSoft/F");
    outTree->Branch("mLowJPsi_isMedium",    &out_mLowJPsi_isMedium,         "mLowJPsi_isMedium/F");
    outTree->Branch("mLowJPsi_isHighPt",    &out_mLowJPsi_isHighPt,         "mLowJPsi_isHighPt/F");
    outTree->Branch("mHighPhi_isLoose",     &out_mHighPhi_isLoose,  "mHighPhi_isLoose/F");
    outTree->Branch("mHighPhi_isSoft",      &out_mHighPhi_isSoft,   "mHighPhi_isSoft/F");
    outTree->Branch("mHighPhi_isMedium",    &out_mHighPhi_isMedium,         "mHighPhi_isMedium/F");
    outTree->Branch("mHighPhi_isHighPt",    &out_mHighPhi_isHighPt,         "mHighPhi_isHighPt/F");
    outTree->Branch("mLowPhi_isLoose",      &out_mLowPhi_isLoose,   "mLowPhi_isLoose/F");
    outTree->Branch("mLowPhi_isSoft",       &out_mLowPhi_isSoft,    "mLowPhi_isSoft/F");
    outTree->Branch("mLowPhi_isMedium",     &out_mLowPhi_isMedium,  "mLowPhi_isMedium/F");
    outTree->Branch("mLowPhi_isHighPt",     &out_mLowPhi_isHighPt,  "mLowPhi_isHighPt/F");
    outTree->Branch("mHighJPsi_isTracker",  &out_mHighJPsi_isTracker,       "mHighJPsi_isTracker/F");
    outTree->Branch("mHighJPsi_isGlobal",   &out_mHighJPsi_isGlobal,        "mHighJPsi_isGlobal/F");
    outTree->Branch("mLowJPsi_isTracker",   &out_mLowJPsi_isTracker,        "mLowJPsi_isTracker/F");
    outTree->Branch("mLowJPsi_isGlobal",    &out_mLowJPsi_isGlobal,         "mLowJPsi_isGlobal/F");
    outTree->Branch("mHighPhi_isTracker",   &out_mHighPhi_isTracker,        "mHighPhi_isTracker/F");
    outTree->Branch("mHighPhi_isGlobal",    &out_mHighPhi_isGlobal,         "mHighPhi_isGlobal/F");
    outTree->Branch("mLowPhi_isTracker",    &out_mLowPhi_isTracker,         "mLowPhi_isTracker/F");
    outTree->Branch("mLowPhi_isGlobal",     &out_mLowPhi_isGlobal,  "mLowPhi_isGlobal/F");
    outTree->Branch("mHighJPsi_type",       &out_mHighJPsi_type,    "mHighJPsi_type/F");
    outTree->Branch("mLowJPsi_type",        &out_mLowJPsi_type,     "mLowJPsi_type/F");
    outTree->Branch("mHighPhi_type",        &out_mHighPhi_type,     "mHighPhi_type/F");
    outTree->Branch("mLowPhi_type",         &out_mLowPhi_type,      "mLowPhi_type/F");
    outTree->Branch("jpsi_m",       &out_jpsi_m,    "jpsi_m/F");
    outTree->Branch("jpsi_m_rf",    &out_jpsi_m_rf,         "jpsi_m_rf/F");
    outTree->Branch("jpsi_m_rf_c",  &out_jpsi_m_rf_c,       "jpsi_m_rf_c/F");
    outTree->Branch("jpsi_m_rf_d_c",        &out_jpsi_m_rf_d_c,     "jpsi_m_rf_d_c/F");
    outTree->Branch("jpsi_p",       &out_jpsi_p,    "jpsi_p/F");
    outTree->Branch("jpsi_pt",      &out_jpsi_pt,   "jpsi_pt/F");
    outTree->Branch("jpsi_eta",     &out_jpsi_eta,  "jpsi_eta/F");
    outTree->Branch("jpsi_theta",   &out_jpsi_theta,        "jpsi_theta/F");
    outTree->Branch("jpsi_y",       &out_jpsi_y,    "jpsi_y/F");
    outTree->Branch("jpsi_e",       &out_jpsi_e,    "jpsi_e/F");
    outTree->Branch("jpsi_dxy",     &out_jpsi_dxy,  "jpsi_dxy/F");
    outTree->Branch("jpsi_dxyErr",  &out_jpsi_dxyErr,       "jpsi_dxyErr/F");
    outTree->Branch("jpsi_dz",      &out_jpsi_dz,   "jpsi_dz/F");
    outTree->Branch("jpsi_dzErr",   &out_jpsi_dzErr,        "jpsi_dzErr/F");
    outTree->Branch("jpsi_vProb",   &out_jpsi_vProb,        "jpsi_vProb/F");
    outTree->Branch("jpsi_vChi2",   &out_jpsi_vChi2,        "jpsi_vChi2/F");
    outTree->Branch("jpsi_DCA",     &out_jpsi_DCA,  "jpsi_DCA/F");
    outTree->Branch("jpsi_ctauPV",  &out_jpsi_ctauPV,       "jpsi_ctauPV/F");
    outTree->Branch("jpsi_ctauErrPV",       &out_jpsi_ctauErrPV,    "jpsi_ctauErrPV/F");
    outTree->Branch("jpsi_cosAlpha",        &out_jpsi_cosAlpha,     "jpsi_cosAlpha/F");
    outTree->Branch("jpsi_lxy",     &out_jpsi_lxy,  "jpsi_lxy/F");
    outTree->Branch("jpsi_lxyz",    &out_jpsi_lxyz,         "jpsi_lxyz/F");
    outTree->Branch("jpsi_lxyErr",  &out_jpsi_lxyErr,       "jpsi_lxyErr/F");
    outTree->Branch("jpsi_lxyzErr",         &out_jpsi_lxyzErr,      "jpsi_lxyzErr/F");
    outTree->Branch("jpsi_cosAlpha3D",      &out_jpsi_cosAlpha3D,   "jpsi_cosAlpha3D/F");
    outTree->Branch("phi_m",        &out_phi_m,     "phi_m/F");
    outTree->Branch("phi_m_rf",     &out_phi_m_rf,  "phi_m_rf/F");
    outTree->Branch("phi_m_rf_c",   &out_phi_m_rf_c,        "phi_m_rf_c/F");
    outTree->Branch("phi_m_rf_d_c",         &out_phi_m_rf_d_c,      "phi_m_rf_d_c/F");
    outTree->Branch("phi_p",        &out_phi_p,     "phi_p/F");
    outTree->Branch("phi_pt",       &out_phi_pt,    "phi_pt/F");
    outTree->Branch("phi_eta",      &out_phi_eta,   "phi_eta/F");
    outTree->Branch("phi_theta",    &out_phi_theta,         "phi_theta/F");
    outTree->Branch("phi_y",        &out_phi_y,     "phi_y/F");
    outTree->Branch("phi_e",        &out_phi_e,     "phi_e/F");
    outTree->Branch("phi_dxy",      &out_phi_dxy,   "phi_dxy/F");
    outTree->Branch("phi_dxyErr",   &out_phi_dxyErr,        "phi_dxyErr/F");
    outTree->Branch("phi_dz",       &out_phi_dz,    "phi_dz/F");
    outTree->Branch("phi_dzErr",    &out_phi_dzErr,         "phi_dzErr/F");
    outTree->Branch("phi_vProb",    &out_phi_vProb,         "phi_vProb/F");
    outTree->Branch("phi_vChi2",    &out_phi_vChi2,         "phi_vChi2/F");
    outTree->Branch("phi_DCA",      &out_phi_DCA,   "phi_DCA/F");
    outTree->Branch("phi_ctauPV",   &out_phi_ctauPV,        "phi_ctauPV/F");
    outTree->Branch("phi_ctauErrPV",        &out_phi_ctauErrPV,     "phi_ctauErrPV/F");
    outTree->Branch("phi_cosAlpha",         &out_phi_cosAlpha,      "phi_cosAlpha/F");
    outTree->Branch("phi_lxy",      &out_phi_lxy,   "phi_lxy/F");
    outTree->Branch("phi_lxyz",     &out_phi_lxyz,  "phi_lxyz/F");
    outTree->Branch("phi_lxyErr",   &out_phi_lxyErr,        "phi_lxyErr/F");
    outTree->Branch("phi_lxyzErr",  &out_phi_lxyzErr,       "phi_lxyzErr/F");
    outTree->Branch("phi_cosAlpha3D",       &out_phi_cosAlpha3D,    "phi_cosAlpha3D/F");
    outTree->Branch("isBestCandidate",      &out_isBestCandidate,   "isBestCandidate/F");


}

Bool_t FourMuons::Process(Long64_t entry)
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



    out_run =       (Float_t)(*run);
if(false)
{
    out_event =     (Float_t)(*event);
    out_numPrimaryVertices =        (Float_t)(*numPrimaryVertices);
    out_trigger =   (Float_t)(*trigger);
    out_noXCandidates =     (Float_t)(*noXCandidates);

    out_jpsi_triggerMatch =         (Float_t)(*jpsi_triggerMatch);
    out_phi_triggerMatch =  (Float_t)(*phi_triggerMatch);
    out_mHighJPsiMatch =    (Float_t)(*mHighJPsiMatch);
    out_mLowJPsiMatch =     (Float_t)(*mLowJPsiMatch);
    out_mHighPhiMatch =     (Float_t)(*mHighPhiMatch);
    out_mLowPhiMatch =      (Float_t)(*mLowPhiMatch);
    out_doubledimuon_charge =       (Float_t)(*doubledimuon_charge);
    out_doubledimuon_m =    (Float_t)(*doubledimuon_m);
    out_doubledimuon_m_rf =         (Float_t)(*doubledimuon_m_rf);
    out_doubledimuon_m_rf_c =       (Float_t)(*doubledimuon_m_rf_c);
    out_doubledimuon_m_rf_d_c =     (Float_t)(*doubledimuon_m_rf_d_c);
    out_doubledimuon_p =    (Float_t)(*doubledimuon_p);
    out_doubledimuon_pt =   (Float_t)(*doubledimuon_pt);
    out_doubledimuon_e =    (Float_t)(*doubledimuon_e);
    out_doubledimuon_eta =  (Float_t)(*doubledimuon_eta);
    out_doubledimuon_theta =        (Float_t)(*doubledimuon_theta);
    out_doubledimuon_y =    (Float_t)(*doubledimuon_y);
    out_doubledimuon_dxy =  (Float_t)(*doubledimuon_dxy);
    out_doubledimuon_dxyErr =       (Float_t)(*doubledimuon_dxyErr);
    out_doubledimuon_dz =   (Float_t)(*doubledimuon_dz);
    out_doubledimuon_dzErr =        (Float_t)(*doubledimuon_dzErr);
    out_doubledimuon_vProb =        (Float_t)(*doubledimuon_vProb);
    out_doubledimuon_vChi2 =        (Float_t)(*doubledimuon_vChi2);
    out_doubledimuon_nDof =         (Float_t)(*doubledimuon_nDof);
    out_doubledimuon_rf_vProb =     (Float_t)(*doubledimuon_rf_vProb);
    out_doubledimuon_rf_nDof =      (Float_t)(*doubledimuon_rf_nDof);
    out_doubledimuon_rf_c_vProb =   (Float_t)(*doubledimuon_rf_c_vProb);
    out_doubledimuon_rf_c_vChi2 =   (Float_t)(*doubledimuon_rf_c_vChi2);
    out_doubledimuon_rf_c_nDof =    (Float_t)(*doubledimuon_rf_c_nDof);
    out_doubledimuon_vx =   (Float_t)(*doubledimuon_vx);
    out_doubledimuon_vy =   (Float_t)(*doubledimuon_vy);
    out_doubledimuon_vz =   (Float_t)(*doubledimuon_vz);
    out_pv_x =      (Float_t)(*pv_x);
    out_pv_y =      (Float_t)(*pv_y);
    out_pv_z =      (Float_t)(*pv_z);
    out_doubledimuon_cosAlpha =     (Float_t)(*doubledimuon_cosAlpha);
    out_doubledimuon_cosAlpha3D =   (Float_t)(*doubledimuon_cosAlpha3D);
    out_doubledimuon_ctauPV =       (Float_t)(*doubledimuon_ctauPV);
    out_doubledimuon_ctauErrPV =    (Float_t)(*doubledimuon_ctauErrPV);
    out_doubledimuon_lxy =  (Float_t)(*doubledimuon_lxy);
    out_doubledimuon_lxyErr =       (Float_t)(*doubledimuon_lxyErr);
    out_doubledimuon_lxyz =         (Float_t)(*doubledimuon_lxyz);
    out_doubledimuon_lxyzErr =      (Float_t)(*doubledimuon_lxyzErr);
    out_doubledimuon_rf_cosAlpha =  (Float_t)(*doubledimuon_rf_cosAlpha);
    out_doubledimuon_rf_ctauPV =    (Float_t)(*doubledimuon_rf_ctauPV);
    out_doubledimuon_rf_ctauErrPV =         (Float_t)(*doubledimuon_rf_ctauErrPV);
    out_doubledimuon_rf_lxy =       (Float_t)(*doubledimuon_rf_lxy);
    out_doubledimuon_rf_lxyErr =    (Float_t)(*doubledimuon_rf_lxyErr);
    out_doubledimuon_rf_lxyz =      (Float_t)(*doubledimuon_rf_lxyz);
    out_doubledimuon_rf_lxyzErr =   (Float_t)(*doubledimuon_rf_lxyzErr);
    out_doubledimuon_rf_c_cosAlpha =        (Float_t)(*doubledimuon_rf_c_cosAlpha);
    out_doubledimuon_rf_c_ctauPV =  (Float_t)(*doubledimuon_rf_c_ctauPV);
    out_doubledimuon_rf_c_ctauErrPV =       (Float_t)(*doubledimuon_rf_c_ctauErrPV);
    out_doubledimuon_rf_c_lxy =     (Float_t)(*doubledimuon_rf_c_lxy);
    out_doubledimuon_rf_c_lxyErr =  (Float_t)(*doubledimuon_rf_c_lxyErr);
    out_doubledimuon_rf_c_lxyz =    (Float_t)(*doubledimuon_rf_c_lxyz);
    out_doubledimuon_rf_c_lxyzErr =         (Float_t)(*doubledimuon_rf_c_lxyzErr);
    out_doubledimuon_dca_mp1mp2 =   (Float_t)(*doubledimuon_dca_mp1mp2);
    out_doubledimuon_dca_mp1mj1 =   (Float_t)(*doubledimuon_dca_mp1mj1);
    out_doubledimuon_dca_mp1mj2 =   (Float_t)(*doubledimuon_dca_mp1mj2);
    out_doubledimuon_dca_mp2mj1 =   (Float_t)(*doubledimuon_dca_mp2mj1);
    out_doubledimuon_dca_mp2mj2 =   (Float_t)(*doubledimuon_dca_mp2mj2);
    out_doubledimuon_dca_mj1mj2 =   (Float_t)(*doubledimuon_dca_mj1mj2);
    out_doubledimuon_cosAlphaDZ =   (Float_t)(*doubledimuon_cosAlphaDZ);
    out_doubledimuon_cosAlpha3DDZ =         (Float_t)(*doubledimuon_cosAlpha3DDZ);
    out_doubledimuon_ctauPVDZ =     (Float_t)(*doubledimuon_ctauPVDZ);
    out_doubledimuon_ctauErrPVDZ =  (Float_t)(*doubledimuon_ctauErrPVDZ);
    out_doubledimuon_lxyDZ =        (Float_t)(*doubledimuon_lxyDZ);
    out_doubledimuon_lxyErrDZ =     (Float_t)(*doubledimuon_lxyErrDZ);
    out_doubledimuon_lxyzDZ =       (Float_t)(*doubledimuon_lxyzDZ);
    out_doubledimuon_lxyzErrDZ =    (Float_t)(*doubledimuon_lxyzErrDZ);
    out_doubledimuon_cosAlphaBS =   (Float_t)(*doubledimuon_cosAlphaBS);
    out_doubledimuon_cosAlpha3DBS =         (Float_t)(*doubledimuon_cosAlpha3DBS);
    out_doubledimuon_ctauPVBS =     (Float_t)(*doubledimuon_ctauPVBS);
    out_doubledimuon_ctauErrPVBS =  (Float_t)(*doubledimuon_ctauErrPVBS);
    out_doubledimuon_lxyBS =        (Float_t)(*doubledimuon_lxyBS);
    out_doubledimuon_lxyErrBS =     (Float_t)(*doubledimuon_lxyErrBS);
    out_doubledimuon_lxyzBS =       (Float_t)(*doubledimuon_lxyzBS);
    out_doubledimuon_lxyzErrBS =    (Float_t)(*doubledimuon_lxyzErrBS);
    out_mHighPhi_p =        (Float_t)(*mHighPhi_p);
    out_mLowPhi_p =         (Float_t)(*mLowPhi_p);
    out_mHighJPsi_p =       (Float_t)(*mHighJPsi_p);
    out_mLowJPsi_p =        (Float_t)(*mLowJPsi_p);
    out_mHighPhi_pt =       (Float_t)(*mHighPhi_pt);
    out_mLowPhi_pt =        (Float_t)(*mLowPhi_pt);
    out_mHighJPsi_pt =      (Float_t)(*mHighJPsi_pt);
    out_mLowJPsi_pt =       (Float_t)(*mLowJPsi_pt);
    out_mHighPhi_ptErr =    (Float_t)(*mHighPhi_ptErr);
    out_mLowPhi_ptErr =     (Float_t)(*mLowPhi_ptErr);
    out_mHighJPsi_ptErr =   (Float_t)(*mHighJPsi_ptErr);
    out_mLowJPsi_ptErr =    (Float_t)(*mLowJPsi_ptErr);
    out_mHighJPsi_eta =     (Float_t)(*mHighJPsi_eta);
    out_mLowJPsi_eta =      (Float_t)(*mLowJPsi_eta);
    out_mHighPhi_eta =      (Float_t)(*mHighPhi_eta);
    out_mLowPhi_eta =       (Float_t)(*mLowPhi_eta);
    out_mHighJPsi_etaErr =  (Float_t)(*mHighJPsi_etaErr);
    out_mLowJPsi_etaErr =   (Float_t)(*mLowJPsi_etaErr);
    out_mHighPhi_etaErr =   (Float_t)(*mHighPhi_etaErr);
    out_mLowPhi_etaErr =    (Float_t)(*mLowPhi_etaErr);
    out_mHighJPsi_phi =     (Float_t)(*mHighJPsi_phi);
    out_mLowJPsi_phi =      (Float_t)(*mLowJPsi_phi);
    out_mHighPhi_phi =      (Float_t)(*mHighPhi_phi);
    out_mLowPhi_phi =       (Float_t)(*mLowPhi_phi);
    out_mHighJPsi_phiErr =  (Float_t)(*mHighJPsi_phiErr);
    out_mLowJPsi_phiErr =   (Float_t)(*mLowJPsi_phiErr);
    out_mHighPhi_phiErr =   (Float_t)(*mHighPhi_phiErr);
    out_mLowPhi_phiErr =    (Float_t)(*mLowPhi_phiErr);
    out_mHighJPsi_theta =   (Float_t)(*mHighJPsi_theta);
    out_mLowJPsi_theta =    (Float_t)(*mLowJPsi_theta);
    out_mHighPhi_theta =    (Float_t)(*mHighPhi_theta);
    out_mLowPhi_theta =     (Float_t)(*mLowPhi_theta);
    out_mHighJPsi_thetaErr =        (Float_t)(*mHighJPsi_thetaErr);
    out_mLowJPsi_thetaErr =         (Float_t)(*mLowJPsi_thetaErr);
    out_mHighPhi_thetaErr =         (Float_t)(*mHighPhi_thetaErr);
    out_mLowPhi_thetaErr =  (Float_t)(*mLowPhi_thetaErr);
    out_mHighJPsi_lambda =  (Float_t)(*mHighJPsi_lambda);
    out_mLowJPsi_lambda =   (Float_t)(*mLowJPsi_lambda);
    out_mHighPhi_lambda =   (Float_t)(*mHighPhi_lambda);
    out_mLowPhi_lambda =    (Float_t)(*mLowPhi_lambda);
    out_mHighJPsi_lambdaErr =       (Float_t)(*mHighJPsi_lambdaErr);
    out_mLowJPsi_lambdaErr =        (Float_t)(*mLowJPsi_lambdaErr);
    out_mHighPhi_lambdaErr =        (Float_t)(*mHighPhi_lambdaErr);
    out_mLowPhi_lambdaErr =         (Float_t)(*mLowPhi_lambdaErr);
    out_mLowPhi_dxy =       (Float_t)(*mLowPhi_dxy);
    out_mLowPhi_dxyErr =    (Float_t)(*mLowPhi_dxyErr);
    out_mLowPhi_dz =        (Float_t)(*mLowPhi_dz);
    out_mLowPhi_dzErr =     (Float_t)(*mLowPhi_dzErr);
    out_mHighPhi_dxy =      (Float_t)(*mHighPhi_dxy);
    out_mHighPhi_dxyErr =   (Float_t)(*mHighPhi_dxyErr);
    out_mHighPhi_dz =       (Float_t)(*mHighPhi_dz);
    out_mHighPhi_dzErr =    (Float_t)(*mHighPhi_dzErr);
    out_mHighJPsi_dxy =     (Float_t)(*mHighJPsi_dxy);
    out_mHighJPsi_dxyErr =  (Float_t)(*mHighJPsi_dxyErr);
    out_mHighJPsi_dz =      (Float_t)(*mHighJPsi_dz);
    out_mHighJPsi_dzErr =   (Float_t)(*mHighJPsi_dzErr);
    out_mLowJPsi_dxy =      (Float_t)(*mLowJPsi_dxy);
    out_mLowJPsi_dxyErr =   (Float_t)(*mLowJPsi_dxyErr);
    out_mLowJPsi_dz =       (Float_t)(*mLowJPsi_dz);
    out_mLowJPsi_dzErr =    (Float_t)(*mLowJPsi_dzErr);
    out_mHighJPsi_NPixelHits =      (Float_t)(*mHighJPsi_NPixelHits);
    out_mHighJPsi_NStripHits =      (Float_t)(*mHighJPsi_NStripHits);
    out_mHighJPsi_NTrackhits =      (Float_t)(*mHighJPsi_NTrackhits);
    out_mHighJPsi_NBPixHits =       (Float_t)(*mHighJPsi_NBPixHits);
    out_mHighJPsi_NPixLayers =      (Float_t)(*mHighJPsi_NPixLayers);
    out_mHighJPsi_NTraLayers =      (Float_t)(*mHighJPsi_NTraLayers);
    out_mHighJPsi_NStrLayers =      (Float_t)(*mHighJPsi_NStrLayers);
    out_mHighJPsi_NBPixLayers =     (Float_t)(*mHighJPsi_NBPixLayers);
    out_mLowJPsi_NPixelHits =       (Float_t)(*mLowJPsi_NPixelHits);
    out_mLowJPsi_NStripHits =       (Float_t)(*mLowJPsi_NStripHits);
    out_mLowJPsi_NTrackhits =       (Float_t)(*mLowJPsi_NTrackhits);
    out_mLowJPsi_NBPixHits =        (Float_t)(*mLowJPsi_NBPixHits);
    out_mLowJPsi_NPixLayers =       (Float_t)(*mLowJPsi_NPixLayers);
    out_mLowJPsi_NTraLayers =       (Float_t)(*mLowJPsi_NTraLayers);
    out_mLowJPsi_NStrLayers =       (Float_t)(*mLowJPsi_NStrLayers);
    out_mLowJPsi_NBPixLayers =      (Float_t)(*mLowJPsi_NBPixLayers);
    out_mHighPhi_NPixelHits =       (Float_t)(*mHighPhi_NPixelHits);
    out_mHighPhi_NStripHits =       (Float_t)(*mHighPhi_NStripHits);
    out_mHighPhi_NTrackhits =       (Float_t)(*mHighPhi_NTrackhits);
    out_mHighPhi_NBPixHits =        (Float_t)(*mHighPhi_NBPixHits);
    out_mHighPhi_NPixLayers =       (Float_t)(*mHighPhi_NPixLayers);
    out_mHighPhi_NTraLayers =       (Float_t)(*mHighPhi_NTraLayers);
    out_mHighPhi_NStrLayers =       (Float_t)(*mHighPhi_NStrLayers);
    out_mHighPhi_NBPixLayers =      (Float_t)(*mHighPhi_NBPixLayers);
    out_mLowPhi_NPixelHits =        (Float_t)(*mLowPhi_NPixelHits);
    out_mLowPhi_NStripHits =        (Float_t)(*mLowPhi_NStripHits);
    out_mLowPhi_NTrackhits =        (Float_t)(*mLowPhi_NTrackhits);
    out_mLowPhi_NBPixHits =         (Float_t)(*mLowPhi_NBPixHits);
    out_mLowPhi_NPixLayers =        (Float_t)(*mLowPhi_NPixLayers);
    out_mLowPhi_NTraLayers =        (Float_t)(*mLowPhi_NTraLayers);
    out_mLowPhi_NStrLayers =        (Float_t)(*mLowPhi_NStrLayers);
    out_mLowPhi_NBPixLayers =       (Float_t)(*mLowPhi_NBPixLayers);
    out_mHighJPsi_isLoose =         (Float_t)(*mHighJPsi_isLoose);
    out_mHighJPsi_isSoft =  (Float_t)(*mHighJPsi_isSoft);
    out_mHighJPsi_isMedium =        (Float_t)(*mHighJPsi_isMedium);
    out_mHighJPsi_isHighPt =        (Float_t)(*mHighJPsi_isHighPt);
    out_mLowJPsi_isLoose =  (Float_t)(*mLowJPsi_isLoose);
    out_mLowJPsi_isSoft =   (Float_t)(*mLowJPsi_isSoft);
    out_mLowJPsi_isMedium =         (Float_t)(*mLowJPsi_isMedium);
    out_mLowJPsi_isHighPt =         (Float_t)(*mLowJPsi_isHighPt);
    out_mHighPhi_isLoose =  (Float_t)(*mHighPhi_isLoose);
    out_mHighPhi_isSoft =   (Float_t)(*mHighPhi_isSoft);
    out_mHighPhi_isMedium =         (Float_t)(*mHighPhi_isMedium);
    out_mHighPhi_isHighPt =         (Float_t)(*mHighPhi_isHighPt);
    out_mLowPhi_isLoose =   (Float_t)(*mLowPhi_isLoose);
    out_mLowPhi_isSoft =    (Float_t)(*mLowPhi_isSoft);
    out_mLowPhi_isMedium =  (Float_t)(*mLowPhi_isMedium);
    out_mLowPhi_isHighPt =  (Float_t)(*mLowPhi_isHighPt);
    out_mHighJPsi_isTracker =       (Float_t)(*mHighJPsi_isTracker);
    out_mHighJPsi_isGlobal =        (Float_t)(*mHighJPsi_isGlobal);
    out_mLowJPsi_isTracker =        (Float_t)(*mLowJPsi_isTracker);
    out_mLowJPsi_isGlobal =         (Float_t)(*mLowJPsi_isGlobal);
    out_mHighPhi_isTracker =        (Float_t)(*mHighPhi_isTracker);
    out_mHighPhi_isGlobal =         (Float_t)(*mHighPhi_isGlobal);
    out_mLowPhi_isTracker =         (Float_t)(*mLowPhi_isTracker);
    out_mLowPhi_isGlobal =  (Float_t)(*mLowPhi_isGlobal);
    out_mHighJPsi_type =    (Float_t)(*mHighJPsi_type);
    out_mLowJPsi_type =     (Float_t)(*mLowJPsi_type);
    out_mHighPhi_type =     (Float_t)(*mHighPhi_type);
    out_mLowPhi_type =      (Float_t)(*mLowPhi_type);
    out_jpsi_m =    (Float_t)(*jpsi_m);
    out_jpsi_m_rf =         (Float_t)(*jpsi_m_rf);
    out_jpsi_m_rf_c =       (Float_t)(*jpsi_m_rf_c);
    out_jpsi_m_rf_d_c =     (Float_t)(*jpsi_m_rf_d_c);
    out_jpsi_p =    (Float_t)(*jpsi_p);
    out_jpsi_pt =   (Float_t)(*jpsi_pt);
    out_jpsi_eta =  (Float_t)(*jpsi_eta);
    out_jpsi_theta =        (Float_t)(*jpsi_theta);
    out_jpsi_y =    (Float_t)(*jpsi_y);
    out_jpsi_e =    (Float_t)(*jpsi_e);
    out_jpsi_dxy =  (Float_t)(*jpsi_dxy);
    out_jpsi_dxyErr =       (Float_t)(*jpsi_dxyErr);
    out_jpsi_dz =   (Float_t)(*jpsi_dz);
    out_jpsi_dzErr =        (Float_t)(*jpsi_dzErr);
    out_jpsi_vProb =        (Float_t)(*jpsi_vProb);
    out_jpsi_vChi2 =        (Float_t)(*jpsi_vChi2);
    out_jpsi_DCA =  (Float_t)(*jpsi_DCA);
    out_jpsi_ctauPV =       (Float_t)(*jpsi_ctauPV);
    out_jpsi_ctauErrPV =    (Float_t)(*jpsi_ctauErrPV);
    out_jpsi_cosAlpha =     (Float_t)(*jpsi_cosAlpha);
    out_jpsi_lxy =  (Float_t)(*jpsi_lxy);
    out_jpsi_lxyz =         (Float_t)(*jpsi_lxyz);
    out_jpsi_lxyErr =       (Float_t)(*jpsi_lxyErr);
    out_jpsi_lxyzErr =      (Float_t)(*jpsi_lxyzErr);
    out_jpsi_cosAlpha3D =   (Float_t)(*jpsi_cosAlpha3D);
    out_phi_m =     (Float_t)(*phi_m);
    out_phi_m_rf =  (Float_t)(*phi_m_rf);
    out_phi_m_rf_c =        (Float_t)(*phi_m_rf_c);
    out_phi_m_rf_d_c =      (Float_t)(*phi_m_rf_d_c);
    out_phi_p =     (Float_t)(*phi_p);
    out_phi_pt =    (Float_t)(*phi_pt);
    out_phi_eta =   (Float_t)(*phi_eta);
    out_phi_theta =         (Float_t)(*phi_theta);
    out_phi_y =     (Float_t)(*phi_y);
    out_phi_e =     (Float_t)(*phi_e);
    out_phi_dxy =   (Float_t)(*phi_dxy);
    out_phi_dxyErr =        (Float_t)(*phi_dxyErr);
    out_phi_dz =    (Float_t)(*phi_dz);
    out_phi_dzErr =         (Float_t)(*phi_dzErr);
    out_phi_vProb =         (Float_t)(*phi_vProb);
    out_phi_vChi2 =         (Float_t)(*phi_vChi2);
    out_phi_DCA =   (Float_t)(*phi_DCA);
    out_phi_ctauPV =        (Float_t)(*phi_ctauPV);
    out_phi_ctauErrPV =     (Float_t)(*phi_ctauErrPV);
    out_phi_cosAlpha =      (Float_t)(*phi_cosAlpha);
    out_phi_lxy =   (Float_t)(*phi_lxy);
    out_phi_lxyz =  (Float_t)(*phi_lxyz);
    out_phi_lxyErr =        (Float_t)(*phi_lxyErr);
    out_phi_lxyzErr =       (Float_t)(*phi_lxyzErr);
    out_phi_cosAlpha3D =    (Float_t)(*phi_cosAlpha3D);
    out_isBestCandidate =   (Float_t)(*isBestCandidate);
  }

  outTree->Fill();

   return kTRUE;
}

void FourMuons::SlaveTerminate()
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

void FourMuons::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}
