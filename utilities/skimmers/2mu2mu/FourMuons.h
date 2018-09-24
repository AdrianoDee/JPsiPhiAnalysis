//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 24 18:54:25 2018 by ROOT version 6.12/07
// from TTree FourMuonTree/Tree of JPsi and Phi in 4 Muons
// found on file: /lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20180919_163853/180919_143902/0009/sum.root
//////////////////////////////////////////////////////////

#ifndef FourMuons_h
#define FourMuons_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include "TLorentzVector.h"



class FourMuons : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   Float_t JPsi_mass, Phi_mass, Phi_mean, Phi_sigma;
   TTree *outTree;

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Double_t> run = {fReader, "run"};
   TTreeReaderValue<Double_t> event = {fReader, "event"};
   TTreeReaderValue<Double_t> numPrimaryVertices = {fReader, "numPrimaryVertices"};
   TTreeReaderValue<Double_t> trigger = {fReader, "trigger"};
   TTreeReaderValue<Double_t> noXCandidates = {fReader, "noXCandidates"};
   TTreeReaderValue<TLorentzVector> doubledimuon_p4 = {fReader, "doubledimuon_p4"};
   TTreeReaderValue<TLorentzVector> phi_p4 = {fReader, "phi_p4"};
   TTreeReaderValue<TLorentzVector> jpsi_p4 = {fReader, "jpsi_p4"};
   TTreeReaderValue<TLorentzVector> mLowPhi_p4 = {fReader, "mLowPhi_p4"};
   TTreeReaderValue<TLorentzVector> mLowJPsi_p4 = {fReader, "mLowJPsi_p4"};
   TTreeReaderValue<TLorentzVector> mHighJPsi_p4 = {fReader, "mHighJPsi_p4"};
   TTreeReaderValue<TLorentzVector> doubledimuon_rf_p4 = {fReader, "doubledimuon_rf_p4"};
   TTreeReaderValue<TLorentzVector> phi_rf_p4 = {fReader, "phi_rf_p4"};
   TTreeReaderValue<TLorentzVector> jpsi_rf_p4 = {fReader, "jpsi_rf_p4"};
   TTreeReaderValue<TLorentzVector> mLowPhi_rf_p4 = {fReader, "mLowPhi_rf_p4"};
   TTreeReaderValue<TLorentzVector> mLowJPsi_rf_p4 = {fReader, "mLowJPsi_rf_p4"};
   TTreeReaderValue<TLorentzVector> mHighJPsi_rf_p4 = {fReader, "mHighJPsi_rf_p4"};
   TTreeReaderValue<Double_t> jpsi_triggerMatch = {fReader, "jpsi_triggerMatch"};
   TTreeReaderValue<Double_t> phi_triggerMatch = {fReader, "phi_triggerMatch"};
   TTreeReaderValue<Double_t> mHighJPsiMatch = {fReader, "mHighJPsiMatch"};
   TTreeReaderValue<Double_t> mLowJPsiMatch = {fReader, "mLowJPsiMatch"};
   TTreeReaderValue<Double_t> mHighPhiMatch = {fReader, "mHighPhiMatch"};
   TTreeReaderValue<Double_t> mLowPhiMatch = {fReader, "mLowPhiMatch"};
   TTreeReaderValue<Double_t> doubledimuon_charge = {fReader, "doubledimuon_charge"};
   TTreeReaderValue<Double_t> doubledimuon_m = {fReader, "doubledimuon_m"};
   TTreeReaderValue<Double_t> doubledimuon_m_rf = {fReader, "doubledimuon_m_rf"};
   TTreeReaderValue<Double_t> doubledimuon_m_rf_c = {fReader, "doubledimuon_m_rf_c"};
   TTreeReaderValue<Double_t> doubledimuon_m_rf_d_c = {fReader, "doubledimuon_m_rf_d_c"};
   TTreeReaderValue<Double_t> doubledimuon_p = {fReader, "doubledimuon_p"};
   TTreeReaderValue<Double_t> doubledimuon_pt = {fReader, "doubledimuon_pt"};
   TTreeReaderValue<Double_t> doubledimuon_e = {fReader, "doubledimuon_e"};
   TTreeReaderValue<Double_t> doubledimuon_eta = {fReader, "doubledimuon_eta"};
   TTreeReaderValue<Double_t> doubledimuon_theta = {fReader, "doubledimuon_theta"};
   TTreeReaderValue<Double_t> doubledimuon_y = {fReader, "doubledimuon_y"};
   TTreeReaderValue<Double_t> doubledimuon_dxy = {fReader, "doubledimuon_dxy"};
   TTreeReaderValue<Double_t> doubledimuon_dxyErr = {fReader, "doubledimuon_dxyErr"};
   TTreeReaderValue<Double_t> doubledimuon_dz = {fReader, "doubledimuon_dz"};
   TTreeReaderValue<Double_t> doubledimuon_dzErr = {fReader, "doubledimuon_dzErr"};
   TTreeReaderValue<Double_t> doubledimuon_vProb = {fReader, "doubledimuon_vProb"};
   TTreeReaderValue<Double_t> doubledimuon_vChi2 = {fReader, "doubledimuon_vChi2"};
   TTreeReaderValue<Double_t> doubledimuon_nDof = {fReader, "doubledimuon_nDof"};
   TTreeReaderValue<Double_t> doubledimuon_rf_vProb = {fReader, "doubledimuon_rf_vProb"};
   TTreeReaderValue<Double_t>  doubledimuon_rf_vChi2 = {fReader, " doubledimuon_rf_vChi2"};
   TTreeReaderValue<Double_t> doubledimuon_rf_nDof = {fReader, "doubledimuon_rf_nDof"};
   TTreeReaderValue<Double_t> doubledimuon_rf_c_vProb = {fReader, "doubledimuon_rf_c_vProb"};
   TTreeReaderValue<Double_t> doubledimuon_rf_c_vChi2 = {fReader, "doubledimuon_rf_c_vChi2"};
   TTreeReaderValue<Double_t> doubledimuon_rf_c_nDof = {fReader, "doubledimuon_rf_c_nDof"};
   TTreeReaderValue<Double_t> doubledimuon_vx = {fReader, "doubledimuon_vx"};
   TTreeReaderValue<Double_t> doubledimuon_vy = {fReader, "doubledimuon_vy"};
   TTreeReaderValue<Double_t> doubledimuon_vz = {fReader, "doubledimuon_vz"};
   TTreeReaderValue<Double_t> pv_x = {fReader, "pv_x"};
   TTreeReaderValue<Double_t> pv_y = {fReader, "pv_y"};
   TTreeReaderValue<Double_t> pv_z = {fReader, "pv_z"};
   TTreeReaderValue<Double_t> doubledimuon_cosAlpha = {fReader, "doubledimuon_cosAlpha"};
   TTreeReaderValue<Double_t> doubledimuon_cosAlpha3D = {fReader, "doubledimuon_cosAlpha3D"};
   TTreeReaderValue<Double_t> doubledimuon_ctauPV = {fReader, "doubledimuon_ctauPV"};
   TTreeReaderValue<Double_t> doubledimuon_ctauErrPV = {fReader, "doubledimuon_ctauErrPV"};
   TTreeReaderValue<Double_t> doubledimuon_lxy = {fReader, "doubledimuon_lxy"};
   TTreeReaderValue<Double_t> doubledimuon_lxyErr = {fReader, "doubledimuon_lxyErr"};
   TTreeReaderValue<Double_t> doubledimuon_lxyz = {fReader, "doubledimuon_lxyz"};
   TTreeReaderValue<Double_t> doubledimuon_lxyzErr = {fReader, "doubledimuon_lxyzErr"};
   TTreeReaderValue<Double_t> doubledimuon_rf_cosAlpha = {fReader, "doubledimuon_rf_cosAlpha"};
   TTreeReaderValue<Double_t> doubledimuon_rf_ctauPV = {fReader, "doubledimuon_rf_ctauPV"};
   TTreeReaderValue<Double_t> doubledimuon_rf_ctauErrPV = {fReader, "doubledimuon_rf_ctauErrPV"};
   TTreeReaderValue<Double_t> doubledimuon_rf_lxy = {fReader, "doubledimuon_rf_lxy"};
   TTreeReaderValue<Double_t> doubledimuon_rf_lxyErr = {fReader, "doubledimuon_rf_lxyErr"};
   TTreeReaderValue<Double_t> doubledimuon_rf_lxyz = {fReader, "doubledimuon_rf_lxyz"};
   TTreeReaderValue<Double_t> doubledimuon_rf_lxyzErr = {fReader, "doubledimuon_rf_lxyzErr"};
   TTreeReaderValue<Double_t> doubledimuon_rf_c_cosAlpha = {fReader, "doubledimuon_rf_c_cosAlpha"};
   TTreeReaderValue<Double_t> doubledimuon_rf_c_ctauPV = {fReader, "doubledimuon_rf_c_ctauPV"};
   TTreeReaderValue<Double_t> doubledimuon_rf_c_ctauErrPV = {fReader, "doubledimuon_rf_c_ctauErrPV"};
   TTreeReaderValue<Double_t> doubledimuon_rf_c_lxy = {fReader, "doubledimuon_rf_c_lxy"};
   TTreeReaderValue<Double_t> doubledimuon_rf_c_lxyErr = {fReader, "doubledimuon_rf_c_lxyErr"};
   TTreeReaderValue<Double_t> doubledimuon_rf_c_lxyz = {fReader, "doubledimuon_rf_c_lxyz"};
   TTreeReaderValue<Double_t> doubledimuon_rf_c_lxyzErr = {fReader, "doubledimuon_rf_c_lxyzErr"};
   TTreeReaderValue<Double_t> doubledimuon_dca_mp1mp2 = {fReader, "doubledimuon_dca_mp1mp2"};
   TTreeReaderValue<Double_t> doubledimuon_dca_mp1mj1 = {fReader, "doubledimuon_dca_mp1mj1"};
   TTreeReaderValue<Double_t> doubledimuon_dca_mp1mj2 = {fReader, "doubledimuon_dca_mp1mj2"};
   TTreeReaderValue<Double_t> doubledimuon_dca_mp2mj1 = {fReader, "doubledimuon_dca_mp2mj1"};
   TTreeReaderValue<Double_t> doubledimuon_dca_mp2mj2 = {fReader, "doubledimuon_dca_mp2mj2"};
   TTreeReaderValue<Double_t> doubledimuon_dca_mj1mj2 = {fReader, "doubledimuon_dca_mj1mj2"};
   TTreeReaderValue<Double_t> doubledimuon_cosAlphaDZ = {fReader, "doubledimuon_cosAlphaDZ"};
   TTreeReaderValue<Double_t> doubledimuon_cosAlpha3DDZ = {fReader, "doubledimuon_cosAlpha3DDZ"};
   TTreeReaderValue<Double_t> doubledimuon_ctauPVDZ = {fReader, "doubledimuon_ctauPVDZ"};
   TTreeReaderValue<Double_t> doubledimuon_ctauErrPVDZ = {fReader, "doubledimuon_ctauErrPVDZ"};
   TTreeReaderValue<Double_t> doubledimuon_lxyDZ = {fReader, "doubledimuon_lxyDZ"};
   TTreeReaderValue<Double_t> doubledimuon_lxyErrDZ = {fReader, "doubledimuon_lxyErrDZ"};
   TTreeReaderValue<Double_t> doubledimuon_lxyzDZ = {fReader, "doubledimuon_lxyzDZ"};
   TTreeReaderValue<Double_t> doubledimuon_lxyzErrDZ = {fReader, "doubledimuon_lxyzErrDZ"};
   TTreeReaderValue<Double_t> doubledimuon_cosAlphaBS = {fReader, "doubledimuon_cosAlphaBS"};
   TTreeReaderValue<Double_t> doubledimuon_cosAlpha3DBS = {fReader, "doubledimuon_cosAlpha3DBS"};
   TTreeReaderValue<Double_t> doubledimuon_ctauPVBS = {fReader, "doubledimuon_ctauPVBS"};
   TTreeReaderValue<Double_t> doubledimuon_ctauErrPVBS = {fReader, "doubledimuon_ctauErrPVBS"};
   TTreeReaderValue<Double_t> doubledimuon_lxyBS = {fReader, "doubledimuon_lxyBS"};
   TTreeReaderValue<Double_t> doubledimuon_lxyErrBS = {fReader, "doubledimuon_lxyErrBS"};
   TTreeReaderValue<Double_t> doubledimuon_lxyzBS = {fReader, "doubledimuon_lxyzBS"};
   TTreeReaderValue<Double_t> doubledimuon_lxyzErrBS = {fReader, "doubledimuon_lxyzErrBS"};
   TTreeReaderValue<Double_t> mHighPhi_p = {fReader, "mHighPhi_p"};
   TTreeReaderValue<Double_t> mLowPhi_p = {fReader, "mLowPhi_p"};
   TTreeReaderValue<Double_t> mHighJPsi_p = {fReader, "mHighJPsi_p"};
   TTreeReaderValue<Double_t> mLowJPsi_p = {fReader, "mLowJPsi_p"};
   TTreeReaderValue<Double_t> mHighPhi_pt = {fReader, "mHighPhi_pt"};
   TTreeReaderValue<Double_t> mLowPhi_pt = {fReader, "mLowPhi_pt"};
   TTreeReaderValue<Double_t> mHighJPsi_pt = {fReader, "mHighJPsi_pt"};
   TTreeReaderValue<Double_t> mLowJPsi_pt = {fReader, "mLowJPsi_pt"};
   TTreeReaderValue<Double_t> mHighPhi_ptErr = {fReader, "mHighPhi_ptErr"};
   TTreeReaderValue<Double_t> mLowPhi_ptErr = {fReader, "mLowPhi_ptErr"};
   TTreeReaderValue<Double_t> mHighJPsi_ptErr = {fReader, "mHighJPsi_ptErr"};
   TTreeReaderValue<Double_t> mLowJPsi_ptErr = {fReader, "mLowJPsi_ptErr"};
   TTreeReaderValue<Double_t> mHighJPsi_eta = {fReader, "mHighJPsi_eta"};
   TTreeReaderValue<Double_t> mLowJPsi_eta = {fReader, "mLowJPsi_eta"};
   TTreeReaderValue<Double_t> mHighPhi_eta = {fReader, "mHighPhi_eta"};
   TTreeReaderValue<Double_t> mLowPhi_eta = {fReader, "mLowPhi_eta"};
   TTreeReaderValue<Double_t> mHighJPsi_etaErr = {fReader, "mHighJPsi_etaErr"};
   TTreeReaderValue<Double_t> mLowJPsi_etaErr = {fReader, "mLowJPsi_etaErr"};
   TTreeReaderValue<Double_t> mHighPhi_etaErr = {fReader, "mHighPhi_etaErr"};
   TTreeReaderValue<Double_t> mLowPhi_etaErr = {fReader, "mLowPhi_etaErr"};
   TTreeReaderValue<Double_t> mHighJPsi_phi = {fReader, "mHighJPsi_phi"};
   TTreeReaderValue<Double_t> mLowJPsi_phi = {fReader, "mLowJPsi_phi"};
   TTreeReaderValue<Double_t> mHighPhi_phi = {fReader, "mHighPhi_phi"};
   TTreeReaderValue<Double_t> mLowPhi_phi = {fReader, "mLowPhi_phi"};
   TTreeReaderValue<Double_t> mHighJPsi_phiErr = {fReader, "mHighJPsi_phiErr"};
   TTreeReaderValue<Double_t> mLowJPsi_phiErr = {fReader, "mLowJPsi_phiErr"};
   TTreeReaderValue<Double_t> mHighPhi_phiErr = {fReader, "mHighPhi_phiErr"};
   TTreeReaderValue<Double_t> mLowPhi_phiErr = {fReader, "mLowPhi_phiErr"};
   TTreeReaderValue<Double_t> mHighJPsi_theta = {fReader, "mHighJPsi_theta"};
   TTreeReaderValue<Double_t> mLowJPsi_theta = {fReader, "mLowJPsi_theta"};
   TTreeReaderValue<Double_t> mHighPhi_theta = {fReader, "mHighPhi_theta"};
   TTreeReaderValue<Double_t> mLowPhi_theta = {fReader, "mLowPhi_theta"};
   TTreeReaderValue<Double_t> mHighJPsi_thetaErr = {fReader, "mHighJPsi_thetaErr"};
   TTreeReaderValue<Double_t> mLowJPsi_thetaErr = {fReader, "mLowJPsi_thetaErr"};
   TTreeReaderValue<Double_t> mHighPhi_thetaErr = {fReader, "mHighPhi_thetaErr"};
   TTreeReaderValue<Double_t> mLowPhi_thetaErr = {fReader, "mLowPhi_thetaErr"};
   TTreeReaderValue<Double_t> mHighJPsi_lambda = {fReader, "mHighJPsi_lambda"};
   TTreeReaderValue<Double_t> mLowJPsi_lambda = {fReader, "mLowJPsi_lambda"};
   TTreeReaderValue<Double_t> mHighPhi_lambda = {fReader, "mHighPhi_lambda"};
   TTreeReaderValue<Double_t> mLowPhi_lambda = {fReader, "mLowPhi_lambda"};
   TTreeReaderValue<Double_t> mHighJPsi_lambdaErr = {fReader, "mHighJPsi_lambdaErr"};
   TTreeReaderValue<Double_t> mLowJPsi_lambdaErr = {fReader, "mLowJPsi_lambdaErr"};
   TTreeReaderValue<Double_t> mHighPhi_lambdaErr = {fReader, "mHighPhi_lambdaErr"};
   TTreeReaderValue<Double_t> mLowPhi_lambdaErr = {fReader, "mLowPhi_lambdaErr"};
   TTreeReaderValue<Double_t> mLowPhi_dxy = {fReader, "mLowPhi_dxy"};
   TTreeReaderValue<Double_t> mLowPhi_dxyErr = {fReader, "mLowPhi_dxyErr"};
   TTreeReaderValue<Double_t> mLowPhi_dz = {fReader, "mLowPhi_dz"};
   TTreeReaderValue<Double_t> mLowPhi_dzErr = {fReader, "mLowPhi_dzErr"};
   TTreeReaderValue<Double_t> mHighPhi_dxy = {fReader, "mHighPhi_dxy"};
   TTreeReaderValue<Double_t> mHighPhi_dxyErr = {fReader, "mHighPhi_dxyErr"};
   TTreeReaderValue<Double_t> mHighPhi_dz = {fReader, "mHighPhi_dz"};
   TTreeReaderValue<Double_t> mHighPhi_dzErr = {fReader, "mHighPhi_dzErr"};
   TTreeReaderValue<Double_t> mHighJPsi_dxy = {fReader, "mHighJPsi_dxy"};
   TTreeReaderValue<Double_t> mHighJPsi_dxyErr = {fReader, "mHighJPsi_dxyErr"};
   TTreeReaderValue<Double_t> mHighJPsi_dz = {fReader, "mHighJPsi_dz"};
   TTreeReaderValue<Double_t> mHighJPsi_dzErr = {fReader, "mHighJPsi_dzErr"};
   TTreeReaderValue<Double_t> mLowJPsi_dxy = {fReader, "mLowJPsi_dxy"};
   TTreeReaderValue<Double_t> mLowJPsi_dxyErr = {fReader, "mLowJPsi_dxyErr"};
   TTreeReaderValue<Double_t> mLowJPsi_dz = {fReader, "mLowJPsi_dz"};
   TTreeReaderValue<Double_t> mLowJPsi_dzErr = {fReader, "mLowJPsi_dzErr"};
   TTreeReaderValue<Double_t> mHighJPsi_NPixelHits = {fReader, "mHighJPsi_NPixelHits"};
   TTreeReaderValue<Double_t> mHighJPsi_NStripHits = {fReader, "mHighJPsi_NStripHits"};
   TTreeReaderValue<Double_t> mHighJPsi_NTrackhits = {fReader, "mHighJPsi_NTrackhits"};
   TTreeReaderValue<Double_t> mHighJPsi_NBPixHits = {fReader, "mHighJPsi_NBPixHits"};
   TTreeReaderValue<Double_t> mHighJPsi_NPixLayers = {fReader, "mHighJPsi_NPixLayers"};
   TTreeReaderValue<Double_t> mHighJPsi_NTraLayers = {fReader, "mHighJPsi_NTraLayers"};
   TTreeReaderValue<Double_t> mHighJPsi_NStrLayers = {fReader, "mHighJPsi_NStrLayers"};
   TTreeReaderValue<Double_t> mHighJPsi_NBPixLayers = {fReader, "mHighJPsi_NBPixLayers"};
   TTreeReaderValue<Double_t> mLowJPsi_NPixelHits = {fReader, "mLowJPsi_NPixelHits"};
   TTreeReaderValue<Double_t> mLowJPsi_NStripHits = {fReader, "mLowJPsi_NStripHits"};
   TTreeReaderValue<Double_t> mLowJPsi_NTrackhits = {fReader, "mLowJPsi_NTrackhits"};
   TTreeReaderValue<Double_t> mLowJPsi_NBPixHits = {fReader, "mLowJPsi_NBPixHits"};
   TTreeReaderValue<Double_t> mLowJPsi_NPixLayers = {fReader, "mLowJPsi_NPixLayers"};
   TTreeReaderValue<Double_t> mLowJPsi_NTraLayers = {fReader, "mLowJPsi_NTraLayers"};
   TTreeReaderValue<Double_t> mLowJPsi_NStrLayers = {fReader, "mLowJPsi_NStrLayers"};
   TTreeReaderValue<Double_t> mLowJPsi_NBPixLayers = {fReader, "mLowJPsi_NBPixLayers"};
   TTreeReaderValue<Double_t> mHighPhi_NPixelHits = {fReader, "mHighPhi_NPixelHits"};
   TTreeReaderValue<Double_t> mHighPhi_NStripHits = {fReader, "mHighPhi_NStripHits"};
   TTreeReaderValue<Double_t> mHighPhi_NTrackhits = {fReader, "mHighPhi_NTrackhits"};
   TTreeReaderValue<Double_t> mHighPhi_NBPixHits = {fReader, "mHighPhi_NBPixHits"};
   TTreeReaderValue<Double_t> mHighPhi_NPixLayers = {fReader, "mHighPhi_NPixLayers"};
   TTreeReaderValue<Double_t> mHighPhi_NTraLayers = {fReader, "mHighPhi_NTraLayers"};
   TTreeReaderValue<Double_t> mHighPhi_NStrLayers = {fReader, "mHighPhi_NStrLayers"};
   TTreeReaderValue<Double_t> mHighPhi_NBPixLayers = {fReader, "mHighPhi_NBPixLayers"};
   TTreeReaderValue<Double_t> mLowPhi_NPixelHits = {fReader, "mLowPhi_NPixelHits"};
   TTreeReaderValue<Double_t> mLowPhi_NStripHits = {fReader, "mLowPhi_NStripHits"};
   TTreeReaderValue<Double_t> mLowPhi_NTrackhits = {fReader, "mLowPhi_NTrackhits"};
   TTreeReaderValue<Double_t> mLowPhi_NBPixHits = {fReader, "mLowPhi_NBPixHits"};
   TTreeReaderValue<Double_t> mLowPhi_NPixLayers = {fReader, "mLowPhi_NPixLayers"};
   TTreeReaderValue<Double_t> mLowPhi_NTraLayers = {fReader, "mLowPhi_NTraLayers"};
   TTreeReaderValue<Double_t> mLowPhi_NStrLayers = {fReader, "mLowPhi_NStrLayers"};
   TTreeReaderValue<Double_t> mLowPhi_NBPixLayers = {fReader, "mLowPhi_NBPixLayers"};
   TTreeReaderValue<Double_t> mHighJPsi_isLoose = {fReader, "mHighJPsi_isLoose"};
   TTreeReaderValue<Double_t> mHighJPsi_isSoft = {fReader, "mHighJPsi_isSoft"};
   TTreeReaderValue<Double_t> mHighJPsi_isMedium = {fReader, "mHighJPsi_isMedium"};
   TTreeReaderValue<Double_t> mHighJPsi_isHighPt = {fReader, "mHighJPsi_isHighPt"};
   TTreeReaderValue<Double_t> mLowJPsi_isLoose = {fReader, "mLowJPsi_isLoose"};
   TTreeReaderValue<Double_t> mLowJPsi_isSoft = {fReader, "mLowJPsi_isSoft"};
   TTreeReaderValue<Double_t> mLowJPsi_isMedium = {fReader, "mLowJPsi_isMedium"};
   TTreeReaderValue<Double_t> mLowJPsi_isHighPt = {fReader, "mLowJPsi_isHighPt"};
   TTreeReaderValue<Double_t> mHighPhi_isLoose = {fReader, "mHighPhi_isLoose"};
   TTreeReaderValue<Double_t> mHighPhi_isSoft = {fReader, "mHighPhi_isSoft"};
   TTreeReaderValue<Double_t> mHighPhi_isMedium = {fReader, "mHighPhi_isMedium"};
   TTreeReaderValue<Double_t> mHighPhi_isHighPt = {fReader, "mHighPhi_isHighPt"};
   TTreeReaderValue<Double_t> mLowPhi_isLoose = {fReader, "mLowPhi_isLoose"};
   TTreeReaderValue<Double_t> mLowPhi_isSoft = {fReader, "mLowPhi_isSoft"};
   TTreeReaderValue<Double_t> mLowPhi_isMedium = {fReader, "mLowPhi_isMedium"};
   TTreeReaderValue<Double_t> mLowPhi_isHighPt = {fReader, "mLowPhi_isHighPt"};
   TTreeReaderValue<Double_t> mHighJPsi_isTracker = {fReader, "mHighJPsi_isTracker"};
   TTreeReaderValue<Double_t> mHighJPsi_isGlobal = {fReader, "mHighJPsi_isGlobal"};
   TTreeReaderValue<Double_t> mLowJPsi_isTracker = {fReader, "mLowJPsi_isTracker"};
   TTreeReaderValue<Double_t> mLowJPsi_isGlobal = {fReader, "mLowJPsi_isGlobal"};
   TTreeReaderValue<Double_t> mHighPhi_isTracker = {fReader, "mHighPhi_isTracker"};
   TTreeReaderValue<Double_t> mHighPhi_isGlobal = {fReader, "mHighPhi_isGlobal"};
   TTreeReaderValue<Double_t> mLowPhi_isTracker = {fReader, "mLowPhi_isTracker"};
   TTreeReaderValue<Double_t> mLowPhi_isGlobal = {fReader, "mLowPhi_isGlobal"};
   TTreeReaderValue<Double_t> mHighJPsi_type = {fReader, "mHighJPsi_type"};
   TTreeReaderValue<Double_t> mLowJPsi_type = {fReader, "mLowJPsi_type"};
   TTreeReaderValue<Double_t> mHighPhi_type = {fReader, "mHighPhi_type"};
   TTreeReaderValue<Double_t> mLowPhi_type = {fReader, "mLowPhi_type"};
   TTreeReaderValue<Double_t> jpsi_m = {fReader, "jpsi_m"};
   TTreeReaderValue<Double_t> jpsi_m_rf = {fReader, "jpsi_m_rf"};
   TTreeReaderValue<Double_t> jpsi_m_rf_c = {fReader, "jpsi_m_rf_c"};
   TTreeReaderValue<Double_t> jpsi_m_rf_d_c = {fReader, "jpsi_m_rf_d_c"};
   TTreeReaderValue<Double_t> jpsi_p = {fReader, "jpsi_p"};
   TTreeReaderValue<Double_t> jpsi_pt = {fReader, "jpsi_pt"};
   TTreeReaderValue<Double_t> jpsi_eta = {fReader, "jpsi_eta"};
   TTreeReaderValue<Double_t> jpsi_theta = {fReader, "jpsi_theta"};
   TTreeReaderValue<Double_t> jpsi_y = {fReader, "jpsi_y"};
   TTreeReaderValue<Double_t> jpsi_e = {fReader, "jpsi_e"};
   TTreeReaderValue<Double_t> jpsi_dxy = {fReader, "jpsi_dxy"};
   TTreeReaderValue<Double_t> jpsi_dxyErr = {fReader, "jpsi_dxyErr"};
   TTreeReaderValue<Double_t> jpsi_dz = {fReader, "jpsi_dz"};
   TTreeReaderValue<Double_t> jpsi_dzErr = {fReader, "jpsi_dzErr"};
   TTreeReaderValue<Double_t> jpsi_vProb = {fReader, "jpsi_vProb"};
   TTreeReaderValue<Double_t> jpsi_vChi2 = {fReader, "jpsi_vChi2"};
   TTreeReaderValue<Double_t> jpsi_DCA = {fReader, "jpsi_DCA"};
   TTreeReaderValue<Double_t> jpsi_ctauPV = {fReader, "jpsi_ctauPV"};
   TTreeReaderValue<Double_t> jpsi_ctauErrPV = {fReader, "jpsi_ctauErrPV"};
   TTreeReaderValue<Double_t> jpsi_cosAlpha = {fReader, "jpsi_cosAlpha"};
   TTreeReaderValue<Double_t> jpsi_lxy = {fReader, "jpsi_lxy"};
   TTreeReaderValue<Double_t> jpsi_lxyz = {fReader, "jpsi_lxyz"};
   TTreeReaderValue<Double_t> jpsi_lxyErr = {fReader, "jpsi_lxyErr"};
   TTreeReaderValue<Double_t> jpsi_lxyzErr = {fReader, "jpsi_lxyzErr"};
   TTreeReaderValue<Double_t> jpsi_cosAlpha3D = {fReader, "jpsi_cosAlpha3D"};
   TTreeReaderValue<Double_t> phi_m = {fReader, "phi_m"};
   TTreeReaderValue<Double_t> phi_m_rf = {fReader, "phi_m_rf"};
   TTreeReaderValue<Double_t> phi_m_rf_c = {fReader, "phi_m_rf_c"};
   TTreeReaderValue<Double_t> phi_m_rf_d_c = {fReader, "phi_m_rf_d_c"};
   TTreeReaderValue<Double_t> phi_p = {fReader, "phi_p"};
   TTreeReaderValue<Double_t> phi_pt = {fReader, "phi_pt"};
   TTreeReaderValue<Double_t> phi_eta = {fReader, "phi_eta"};
   TTreeReaderValue<Double_t> phi_theta = {fReader, "phi_theta"};
   TTreeReaderValue<Double_t> phi_y = {fReader, "phi_y"};
   TTreeReaderValue<Double_t> phi_e = {fReader, "phi_e"};
   TTreeReaderValue<Double_t> phi_dxy = {fReader, "phi_dxy"};
   TTreeReaderValue<Double_t> phi_dxyErr = {fReader, "phi_dxyErr"};
   TTreeReaderValue<Double_t> phi_dz = {fReader, "phi_dz"};
   TTreeReaderValue<Double_t> phi_dzErr = {fReader, "phi_dzErr"};
   TTreeReaderValue<Double_t> phi_vProb = {fReader, "phi_vProb"};
   TTreeReaderValue<Double_t> phi_vChi2 = {fReader, "phi_vChi2"};
   TTreeReaderValue<Double_t> phi_DCA = {fReader, "phi_DCA"};
   TTreeReaderValue<Double_t> phi_ctauPV = {fReader, "phi_ctauPV"};
   TTreeReaderValue<Double_t> phi_ctauErrPV = {fReader, "phi_ctauErrPV"};
   TTreeReaderValue<Double_t> phi_cosAlpha = {fReader, "phi_cosAlpha"};
   TTreeReaderValue<Double_t> phi_lxy = {fReader, "phi_lxy"};
   TTreeReaderValue<Double_t> phi_lxyz = {fReader, "phi_lxyz"};
   TTreeReaderValue<Double_t> phi_lxyErr = {fReader, "phi_lxyErr"};
   TTreeReaderValue<Double_t> phi_lxyzErr = {fReader, "phi_lxyzErr"};
   TTreeReaderValue<Double_t> phi_cosAlpha3D = {fReader, "phi_cosAlpha3D"};
   TTreeReaderValue<Double_t> isBestCandidate = {fReader, "isBestCandidate"};


   FourMuons(TTree * /*tree*/ =0) { }
   virtual ~FourMuons() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(FourMuons,0);

};

#endif

#ifdef FourMuons_cxx
void FourMuons::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t FourMuons::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}
