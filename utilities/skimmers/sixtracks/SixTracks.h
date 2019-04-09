//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Mar 23 16:45:36 2019 by ROOT version 6.12/07
// from TTree SixTracksTree/Tree of DiMuon and DiTrack
// found on file: rootuple-2018-dimuondiTrack_fivedataD2018_2_0_1.root
//////////////////////////////////////////////////////////

#ifndef SixTracks_h
#define SixTracks_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <TSystem.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TBranch.h>
//#include <TCint.h>
#include <TRandom.h>
#include <TMath.h>
#include <TDirectory.h>
#include "TEnv.h"
#include <TString.h>
#include <TSelector.h>
#include <TProof.h>
#include <TProofOutputFile.h>

#include "TPoint.h"
#include <TH1.h>
#include <TH2.h>
#include <TH2F.h>
#include <TF1.h>
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <map>


// Headers needed by this particular selector
#include "TLorentzVector.h"



class SixTracks : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   Float_t JPsi_mass, Phi_mass, Phi_mean, Phi_sigma;
   TTree *outTree;

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> run = {fReader, "run"};
   TTreeReaderValue<Int_t> event = {fReader, "event"};
   TTreeReaderValue<Int_t> lumi = {fReader, "lumi"};
   TTreeReaderValue<Int_t> numPrimaryVertices = {fReader, "numPrimaryVertices"};
   TTreeReaderValue<Int_t> trigger = {fReader, "trigger"};
   TTreeReaderValue<Int_t> noSixCandidates = {fReader, "noSixCandidates"};
   TTreeReaderValue<Int_t> five_id = {fReader, "five_id"};
   TTreeReaderValue<Int_t> dimuon_id = {fReader, "dimuon_id"};
   TTreeReaderValue<Int_t> p_id = {fReader, "p_id"};
   TTreeReaderValue<Int_t> m_id = {fReader, "m_id"};
   TTreeReaderValue<Int_t> t_id = {fReader, "t_id"};
   TTreeReaderValue<Int_t> f_id = {fReader, "f_id"};
   TTreeReaderValue<TLorentzVector> six_p4 = {fReader, "six_p4"};
   TTreeReaderValue<TLorentzVector> five_p4 = {fReader, "five_p4"};
   TTreeReaderValue<TLorentzVector> dimuonditrk_p4 = {fReader, "dimuonditrk_p4"};
   TTreeReaderValue<TLorentzVector> ditrack_p4 = {fReader, "ditrack_p4"};
   TTreeReaderValue<TLorentzVector> dimuon_p4 = {fReader, "dimuon_p4"};
   TTreeReaderValue<TLorentzVector> lowMuon_p4 = {fReader, "lowMuon_p4"};
   TTreeReaderValue<TLorentzVector> highMuon_p4 = {fReader, "highMuon_p4"};
   TTreeReaderValue<TLorentzVector> highKaon_p4 = {fReader, "highKaon_p4"};
   TTreeReaderValue<TLorentzVector> lowKaon_p4 = {fReader, "lowKaon_p4"};
   TTreeReaderValue<TLorentzVector> thirdKaon_p4 = {fReader, "thirdKaon_p4"};
   TTreeReaderValue<TLorentzVector> fourthKaon_p4 = {fReader, "fourthKaon_p4"};
   TTreeReaderValue<TLorentzVector> highPion_p4 = {fReader, "highPion_p4"};
   TTreeReaderValue<TLorentzVector> lowPion_p4 = {fReader, "lowPion_p4"};
   TTreeReaderValue<TLorentzVector> thirdPion_p4 = {fReader, "thirdPion_p4"};
   TTreeReaderValue<TLorentzVector> fourthPion_p4 = {fReader, "fourthPion_p4"};
   TTreeReaderValue<TLorentzVector> highProton_p4 = {fReader, "highProton_p4"};
   TTreeReaderValue<TLorentzVector> lowProton_p4 = {fReader, "lowProton_p4"};
   TTreeReaderValue<TLorentzVector> thirdProton_p4 = {fReader, "thirdProton_p4"};
   TTreeReaderValue<TLorentzVector> fourthProton_p4 = {fReader, "fourthProton_p4"};
   TTreeReaderValue<Double_t> dimuonditrk_m = {fReader, "dimuonditrk_m"};
   TTreeReaderValue<Double_t> dimuonditrk_pt = {fReader, "dimuonditrk_pt"};
   TTreeReaderValue<Double_t> dimuonditrk_eta = {fReader, "dimuonditrk_eta"};
   TTreeReaderValue<Double_t> dimuonditrk_phi = {fReader, "dimuonditrk_phi"};
   TTreeReaderValue<Double_t> dimuonditrk_p = {fReader, "dimuonditrk_p"};
   TTreeReaderValue<Double_t> dimuon_m = {fReader, "dimuon_m"};
   TTreeReaderValue<Double_t> dimuon_pt = {fReader, "dimuon_pt"};
   TTreeReaderValue<Double_t> dimuon_eta = {fReader, "dimuon_eta"};
   TTreeReaderValue<Double_t> dimuon_phi = {fReader, "dimuon_phi"};
   TTreeReaderValue<Double_t> dimuon_p = {fReader, "dimuon_p"};
   TTreeReaderValue<Double_t> highTrackMatch = {fReader, "highTrackMatch"};
   TTreeReaderValue<Double_t> lowTrackMatch = {fReader, "lowTrackMatch"};
   TTreeReaderValue<Double_t> lowMuonMatch = {fReader, "lowMuonMatch"};
   TTreeReaderValue<Double_t> highMuonMatch = {fReader, "highMuonMatch"};
   TTreeReaderValue<Double_t> thirdTrackMatch = {fReader, "thirdTrackMatch"};
   TTreeReaderValue<Double_t> fourthTrackMatch = {fReader, "fourthTrackMatch"};
   TTreeReaderValue<Double_t> ditrack_m = {fReader, "ditrack_m"};
   TTreeReaderValue<Double_t> diTrackOne_pt = {fReader, "diTrackOne_pt"};
   TTreeReaderValue<Double_t> diTrackOne_eta = {fReader, "diTrackOne_eta"};
   TTreeReaderValue<Double_t> diTrackOne_phi = {fReader, "diTrackOne_phi"};
   TTreeReaderValue<Double_t> diTrackOne_p = {fReader, "diTrackOne_p"};
   TTreeReaderValue<Double_t> diTrackTwo_pt = {fReader, "diTrackTwo_pt"};
   TTreeReaderValue<Double_t> diTrackTwo_eta = {fReader, "diTrackTwo_eta"};
   TTreeReaderValue<Double_t> diTrackTwo_phi = {fReader, "diTrackTwo_phi"};
   TTreeReaderValue<Double_t> diTrackTwo_p = {fReader, "diTrackTwo_p"};
   TTreeReaderValue<Double_t> diTrackThree_pt = {fReader, "diTrackThree_pt"};
   TTreeReaderValue<Double_t> diTrackThree_eta = {fReader, "diTrackThree_eta"};
   TTreeReaderValue<Double_t> diTrackThree_phi = {fReader, "diTrackThree_phi"};
   TTreeReaderValue<Double_t> diTrackThree_p = {fReader, "diTrackThree_p"};
   TTreeReaderValue<Double_t> diTrackFour_pt = {fReader, "diTrackFour_pt"};
   TTreeReaderValue<Double_t> diTrackFour_eta = {fReader, "diTrackFour_eta"};
   TTreeReaderValue<Double_t> diTrackFour_phi = {fReader, "diTrackFour_phi"};
   TTreeReaderValue<Double_t> diTrackFour_p = {fReader, "diTrackFour_p"};
   TTreeReaderValue<Double_t> diTrackFive_pt = {fReader, "diTrackFive_pt"};
   TTreeReaderValue<Double_t> diTrackFive_eta = {fReader, "diTrackFive_eta"};
   TTreeReaderValue<Double_t> diTrackFive_phi = {fReader, "diTrackFive_phi"};
   TTreeReaderValue<Double_t> diTrackFive_p = {fReader, "diTrackFive_p"};
   TTreeReaderValue<Double_t> diTrackSix_pt = {fReader, "diTrackSix_pt"};
   TTreeReaderValue<Double_t> diTrackSix_eta = {fReader, "diTrackSix_eta"};
   TTreeReaderValue<Double_t> diTrackSix_phi = {fReader, "diTrackSix_phi"};
   TTreeReaderValue<Double_t> diTrackSix_p = {fReader, "diTrackSix_p"};
   TTreeReaderValue<Double_t> dimuonDiTrkOne_mmpp = {fReader, "dimuonDiTrkOne_mmpp"};
   TTreeReaderValue<Double_t> dimuonDiTrkTwo_mmpp = {fReader, "dimuonDiTrkTwo_mmpp"};
   TTreeReaderValue<Double_t> dimuonDiTrkThree_mmpp = {fReader, "dimuonDiTrkThree_mmpp"};
   TTreeReaderValue<Double_t> dimuonDiTrkFour_mmpp = {fReader, "dimuonDiTrkFour_mmpp"};

   // TTreeReaderValue<Double_t> dimuonDiTrkOne_mmkk = {fReader, "dimuonDiTrkOne_mmkk"};
   // TTreeReaderValue<Double_t> dimuonDiTrkTwo_mmkk = {fReader, "dimuonDiTrkTwo_mmkk"};
   // TTreeReaderValue<Double_t> dimuonDiTrkThree_mmkk = {fReader, "dimuonDiTrkThree_mmkk"};
   // TTreeReaderValue<Double_t> dimuonDiTrkFour_mmkk = {fReader, "dimuonDiTrkFour_mmkk"};

   TTreeReaderValue<Double_t> highMuon_pt = {fReader, "highMuon_pt"};
   TTreeReaderValue<Double_t> highMuon_eta = {fReader, "highMuon_eta"};
   TTreeReaderValue<Double_t> highMuon_phi = {fReader, "highMuon_phi"};
   TTreeReaderValue<Double_t> highMuon_charge = {fReader, "highMuon_charge"};
   TTreeReaderValue<Double_t> highMuon_dz = {fReader, "highMuon_dz"};
   TTreeReaderValue<Double_t> highMuon_dxy = {fReader, "highMuon_dxy"};
   TTreeReaderValue<Double_t> lowMuon_pt = {fReader, "lowMuon_pt"};
   TTreeReaderValue<Double_t> lowMuon_eta = {fReader, "lowMuon_eta"};
   TTreeReaderValue<Double_t> lowMuon_phi = {fReader, "lowMuon_phi"};
   TTreeReaderValue<Double_t> lowMuon_charge = {fReader, "lowMuon_charge"};
   TTreeReaderValue<Double_t> lowMuon_dz = {fReader, "lowMuon_dz"};
   TTreeReaderValue<Double_t> lowMuon_dxy = {fReader, "lowMuon_dxy"};
   TTreeReaderValue<Double_t> highTrack_pt = {fReader, "highTrack_pt"};
   TTreeReaderValue<Double_t> highTrack_eta = {fReader, "highTrack_eta"};
   TTreeReaderValue<Double_t> highTrack_phi = {fReader, "highTrack_phi"};
   TTreeReaderValue<Double_t> highTrack_charge = {fReader, "highTrack_charge"};
   TTreeReaderValue<Double_t> highTrack_dz = {fReader, "highTrack_dz"};
   TTreeReaderValue<Double_t> highTrack_dxy = {fReader, "highTrack_dxy"};
   TTreeReaderValue<Double_t> lowTrack_pt = {fReader, "lowTrack_pt"};
   TTreeReaderValue<Double_t> lowTrack_eta = {fReader, "lowTrack_eta"};
   TTreeReaderValue<Double_t> lowTrack_phi = {fReader, "lowTrack_phi"};
   TTreeReaderValue<Double_t> lowTrack_charge = {fReader, "lowTrack_charge"};
   TTreeReaderValue<Double_t> lowTrack_dz = {fReader, "lowTrack_dz"};
   TTreeReaderValue<Double_t> lowTrack_dxy = {fReader, "lowTrack_dxy"};
   TTreeReaderValue<Double_t> thirdTrack_pt = {fReader, "thirdTrack_pt"};
   TTreeReaderValue<Double_t> thirdTrack_eta = {fReader, "thirdTrack_eta"};
   TTreeReaderValue<Double_t> thirdTrack_phi = {fReader, "thirdTrack_phi"};
   TTreeReaderValue<Double_t> thirdTrack_charge = {fReader, "thirdTrack_charge"};
   TTreeReaderValue<Double_t> thirdTrack_dz = {fReader, "thirdTrack_dz"};
   TTreeReaderValue<Double_t> thirdTrack_dxy = {fReader, "thirdTrack_dxy"};

   TTreeReaderValue<Double_t> dimuonDiTrkOne_pt = {fReader, "dimuonDiTrkOne_pt"};
   TTreeReaderValue<Double_t> dimuonDiTrkOne_eta = {fReader, "dimuonDiTrkOne_eta"};
   TTreeReaderValue<Double_t> dimuonDiTrkOne_phi = {fReader, "dimuonDiTrkOne_phi"};
   TTreeReaderValue<Double_t> dimuonDiTrkOne_charge = {fReader, "dimuonDiTrkOne_charge"};
   TTreeReaderValue<Double_t> dimuonDiTrkOne_p = {fReader, "dimuonDiTrkOne_p"};
   TTreeReaderValue<Double_t> dimuonDiTrkTwo_pt = {fReader, "dimuonDiTrkTwo_pt"};
   TTreeReaderValue<Double_t> dimuonDiTrkTwo_eta = {fReader, "dimuonDiTrkTwo_eta"};
   TTreeReaderValue<Double_t> dimuonDiTrkTwo_phi = {fReader, "dimuonDiTrkTwo_phi"};
   TTreeReaderValue<Double_t> dimuonDiTrkTwo_charge = {fReader, "dimuonDiTrkTwo_charge"};
   TTreeReaderValue<Double_t> dimuonDiTrkTwo_p = {fReader, "dimuonDiTrkTwo_p"};
   TTreeReaderValue<Double_t> dimuonDiTrkThree_pt = {fReader, "dimuonDiTrkThree_pt"};
   TTreeReaderValue<Double_t> dimuonDiTrkThree_eta = {fReader, "dimuonDiTrkThree_eta"};
   TTreeReaderValue<Double_t> dimuonDiTrkThree_phi = {fReader, "dimuonDiTrkThree_phi"};
   TTreeReaderValue<Double_t> dimuonDiTrkThree_charge = {fReader, "dimuonDiTrkThree_charge"};
   TTreeReaderValue<Double_t> dimuonDiTrkThree_p = {fReader, "dimuonDiTrkThree_p"};

   TTreeReaderValue<Double_t> dimuonDiTrkFour_pt = {fReader, "dimuonDiTrkFour_pt"};
   TTreeReaderValue<Double_t> dimuonDiTrkFour_eta = {fReader, "dimuonDiTrkFour_eta"};
   TTreeReaderValue<Double_t> dimuonDiTrkFour_phi = {fReader, "dimuonDiTrkFour_phi"};
   TTreeReaderValue<Double_t> dimuonDiTrkFour_charge = {fReader, "dimuonDiTrkFour_charge"};
   TTreeReaderValue<Double_t> dimuonDiTrkFour_p = {fReader, "dimuonDiTrkFour_p"};
   TTreeReaderValue<Double_t> dimuonDiTrkFive_pt = {fReader, "dimuonDiTrkFive_pt"};
   TTreeReaderValue<Double_t> dimuonDiTrkFive_eta = {fReader, "dimuonDiTrkFive_eta"};
   TTreeReaderValue<Double_t> dimuonDiTrkFive_phi = {fReader, "dimuonDiTrkFive_phi"};
   TTreeReaderValue<Double_t> dimuonDiTrkFive_charge = {fReader, "dimuonDiTrkFive_charge"};
   TTreeReaderValue<Double_t> dimuonDiTrkFive_p = {fReader, "dimuonDiTrkFive_p"};

   TTreeReaderValue<Double_t> dimuonDiTrkSix_pt = {fReader, "dimuonDiTrkSix_pt"};
   TTreeReaderValue<Double_t> dimuonDiTrkSix_eta = {fReader, "dimuonDiTrkSix_eta"};
   TTreeReaderValue<Double_t> dimuonDiTrkSix_phi = {fReader, "dimuonDiTrkSix_phi"};
   TTreeReaderValue<Double_t> dimuonDiTrkSix_charge = {fReader, "dimuonDiTrkSix_charge"};
   TTreeReaderValue<Double_t> dimuonDiTrkSix_p = {fReader, "dimuonDiTrkSix_p"};

   //
   TTreeReaderValue<Double_t> dimuon_vProb = {fReader, "dimuon_vProb"};
   TTreeReaderValue<Double_t> dimuon_vChi2 = {fReader, "dimuon_vChi2"};
   TTreeReaderValue<Double_t> dimuon_DCA = {fReader, "dimuon_DCA"};
   TTreeReaderValue<Double_t> dimuon_ctauPV = {fReader, "dimuon_ctauPV"};
   TTreeReaderValue<Double_t> dimuon_ctauErrPV = {fReader, "dimuon_ctauErrPV"};
   TTreeReaderValue<Double_t> dimuon_cosAlpha = {fReader, "dimuon_cosAlpha"};
   TTreeReaderValue<Double_t> triTrack_m = {fReader, "triTrack_m"};
   TTreeReaderValue<Double_t> triTrack_pt = {fReader, "triTrack_pt"};
   TTreeReaderValue<Double_t> triTrack_eta = {fReader, "triTrack_eta"};
   TTreeReaderValue<Double_t> triTrack_phi = {fReader, "triTrack_phi"};
   TTreeReaderValue<Double_t> triTrack_charge = {fReader, "triTrack_charge"};
   TTreeReaderValue<Double_t> dimuonditrk_vProb = {fReader, "dimuonditrk_vProb"};
   TTreeReaderValue<Double_t> dimuonditrk_vChi2 = {fReader, "dimuonditrk_vChi2"};
   TTreeReaderValue<Double_t> dimuonditrk_nDof = {fReader, "dimuonditrk_nDof"};
   TTreeReaderValue<Int_t> dimuonditrk_charge = {fReader, "dimuonditrk_charge"};
   TTreeReaderValue<Double_t> dimuonditrk_cosAlpha = {fReader, "dimuonditrk_cosAlpha"};
   TTreeReaderValue<Double_t> dimuonditrk_ctauPV = {fReader, "dimuonditrk_ctauPV"};
   TTreeReaderValue<Double_t> dimuonditrk_ctauErrPV = {fReader, "dimuonditrk_ctauErrPV"};
   TTreeReaderValue<Double_t> dimuonditrk_cosAlphaCA = {fReader, "dimuonditrk_cosAlphaCA"};
   TTreeReaderValue<Double_t> dimuonditrk_ctauPVCA = {fReader, "dimuonditrk_ctauPVCA"};
   TTreeReaderValue<Double_t> dimuonditrk_ctauErrPVCA = {fReader, "dimuonditrk_ctauErrPVCA"};
   TTreeReaderValue<Double_t> dimuonditrk_cosAlphaDZ = {fReader, "dimuonditrk_cosAlphaDZ"};
   TTreeReaderValue<Double_t> dimuonditrk_ctauPVDZ = {fReader, "dimuonditrk_ctauPVDZ"};
   TTreeReaderValue<Double_t> dimuonditrk_ctauErrPVDZ = {fReader, "dimuonditrk_ctauErrPVDZ"};
   TTreeReaderValue<Double_t> dimuonditrk_cosAlphaBS = {fReader, "dimuonditrk_cosAlphaBS"};
   TTreeReaderValue<Double_t> dimuonditrk_ctauPVBS = {fReader, "dimuonditrk_ctauPVBS"};
   TTreeReaderValue<Double_t> dimuonditrk_ctauErrPVBS = {fReader, "dimuonditrk_ctauErrPVBS"};
   TTreeReaderValue<Double_t> dimuonditrk_vx = {fReader, "dimuonditrk_vx"};
   TTreeReaderValue<Double_t> dimuonditrk_vy = {fReader, "dimuonditrk_vy"};
   TTreeReaderValue<Double_t> dimuonditrk_vz = {fReader, "dimuonditrk_vz"};
   // TTreeReaderValue<Double_t> dca_m1m2 = {fReader, "dca_m1m2"};
   // TTreeReaderValue<Double_t> dca_m1t1 = {fReader, "dca_m1t1"};
   // TTreeReaderValue<Double_t> dca_m1t2 = {fReader, "dca_m1t2"};
   // TTreeReaderValue<Double_t> dca_m2t1 = {fReader, "dca_m2t1"};
   // TTreeReaderValue<Double_t> dca_m2t2 = {fReader, "dca_m2t2"};
   // TTreeReaderValue<Double_t> dca_t1t2 = {fReader, "dca_t1t2"};
   // TTreeReaderValue<Double_t> dca_m1t3 = {fReader, "dca_m1t3"};
   // TTreeReaderValue<Double_t> dca_m2t3 = {fReader, "dca_m2t3"};
   // TTreeReaderValue<Double_t> dca_t1t3 = {fReader, "dca_t1t3"};
   // TTreeReaderValue<Double_t> dca_t2t3 = {fReader, "dca_t2t3"};
   // TTreeReaderValue<Double_t> dca_m1t4 = {fReader, "dca_m1t4"};
   // TTreeReaderValue<Double_t> dca_m2t4 = {fReader, "dca_m2t4"};
   // TTreeReaderValue<Double_t> dca_t1t4 = {fReader, "dca_t1t4"};
   // TTreeReaderValue<Double_t> dca_t2t4 = {fReader, "dca_t2t4"};
   // TTreeReaderValue<Double_t> dca_t3t4 = {fReader, "dca_t3t4"};
   // TTreeReaderValue<Double_t> highTrackMuonDR = {fReader, "highTrackMuonDR"};
   // TTreeReaderValue<Double_t> highTrackMuonDP = {fReader, "highTrackMuonDP"};
   // TTreeReaderValue<Double_t> highTrackMuonDPt = {fReader, "highTrackMuonDPt"};
   // TTreeReaderValue<Double_t> lowTrackMuonDR = {fReader, "lowTrackMuonDR"};
   // TTreeReaderValue<Double_t> lowTrackMuonDP = {fReader, "lowTrackMuonDP"};
   // TTreeReaderValue<Double_t> lowTrackMuonDPt = {fReader, "lowTrackMuonDPt"};
   // TTreeReaderValue<Double_t> thirdTrackMuonDR = {fReader, "thirdTrackMuonDR"};
   // TTreeReaderValue<Double_t> thirdTrackMuonDP = {fReader, "thirdTrackMuonDP"};
   // TTreeReaderValue<Double_t> thirdTrackMuonDPt = {fReader, "thirdTrackMuonDPt"};
   // TTreeReaderValue<Double_t> fourthTrackMuonDR = {fReader, "fourthTrackMuonDR"};
   // TTreeReaderValue<Double_t> fourthTrackMuonDP = {fReader, "fourthTrackMuonDP"};
   // TTreeReaderValue<Double_t> fourthTrackMuonDPt = {fReader, "fourthTrackMuonDPt"};
   // TTreeReaderValue<Double_t> tPFromPV = {fReader, "tPFromPV"};
   // TTreeReaderValue<Double_t> tMFromPV = {fReader, "tMFromPV"};
   // TTreeReaderValue<Double_t> tTFromPV = {fReader, "tTFromPV"};
   // TTreeReaderValue<Double_t> tFFromPV = {fReader, "tFFromPV"};
   // TTreeReaderValue<Double_t> tPFromPVCA = {fReader, "tPFromPVCA"};
   // TTreeReaderValue<Double_t> tMFromPVCA = {fReader, "tMFromPVCA"};
   // TTreeReaderValue<Double_t> tTFromPVCA = {fReader, "tTFromPVCA"};
   // TTreeReaderValue<Double_t> tFFromPVCA = {fReader, "tFFromPVCA"};
   // TTreeReaderValue<Double_t> tPFromPVDZ = {fReader, "tPFromPVDZ"};
   // TTreeReaderValue<Double_t> tMFromPVDZ = {fReader, "tMFromPVDZ"};
   // TTreeReaderValue<Double_t> tTFromPVDZ = {fReader, "tTFromPVDZ"};
   // TTreeReaderValue<Double_t> tFFromPVDZ = {fReader, "tFFromPVDZ"};

   // TTreeReaderValue<Double_t> five_m = {fReader, "five_m"};
   // TTreeReaderValue<Double_t> five_m_ref = {fReader, "five_m_ref"};
   // TTreeReaderValue<Double_t> five_mass_ppk = {fReader, "five_mass_ppk"};
   // TTreeReaderValue<Double_t> five_mass_kpp = {fReader, "five_mass_kpp"};
   // TTreeReaderValue<Double_t> five_mass_pkp = {fReader, "five_mass_pkp"};
   // TTreeReaderValue<Double_t> five_mass_ppp = {fReader, "five_mass_ppp"};
   // TTreeReaderValue<Double_t> fiveOne_pt = {fReader, "fiveOne_pt"};
   // TTreeReaderValue<Double_t> fiveOne_eta = {fReader, "fiveOne_eta"};
   // TTreeReaderValue<Double_t> fiveOne_phi = {fReader, "fiveOne_phi"};
   // TTreeReaderValue<Double_t> fiveOne_p = {fReader, "fiveOne_p"};
   // TTreeReaderValue<Double_t> fiveTwo_pt = {fReader, "fiveTwo_pt"};
   // TTreeReaderValue<Double_t> fiveTwo_eta = {fReader, "fiveTwo_eta"};
   // TTreeReaderValue<Double_t> fiveTwo_phi = {fReader, "fiveTwo_phi"};
   // TTreeReaderValue<Double_t> fiveTwo_p = {fReader, "fiveTwo_p"};
   // TTreeReaderValue<Double_t> fiveThree_pt = {fReader, "fiveThree_pt"};
   // TTreeReaderValue<Double_t> fiveThree_eta = {fReader, "fiveThree_eta"};
   // TTreeReaderValue<Double_t> fiveThree_phi = {fReader, "fiveThree_phi"};
   // TTreeReaderValue<Double_t> fiveThree_p = {fReader, "fiveThree_p"};
   // TTreeReaderValue<Double_t> fiveFour_pt = {fReader, "fiveFour_pt"};
   // TTreeReaderValue<Double_t> fiveFour_eta = {fReader, "fiveFour_eta"};
   // TTreeReaderValue<Double_t> fiveFour_phi = {fReader, "fiveFour_phi"};
   // TTreeReaderValue<Double_t> fiveFour_p = {fReader, "fiveFour_p"};
   // TTreeReaderValue<Double_t> fiveFive_pt = {fReader, "fiveFive_pt"};
   // TTreeReaderValue<Double_t> fiveFive_eta = {fReader, "fiveFive_eta"};
   // TTreeReaderValue<Double_t> fiveFive_phi = {fReader, "fiveFive_phi"};
   // TTreeReaderValue<Double_t> fiveFive_p = {fReader, "fiveFive_p"};
   // TTreeReaderValue<Double_t> five_cosAlpha = {fReader, "five_cosAlpha"};
   // TTreeReaderValue<Double_t> five_ctauPV = {fReader, "five_ctauPV"};
   // TTreeReaderValue<Double_t> five_ctauErrPV = {fReader, "five_ctauErrPV"};
   // TTreeReaderValue<Double_t> five_cosAlphaCA = {fReader, "five_cosAlphaCA"};
   // TTreeReaderValue<Double_t> five_ctauPVCA = {fReader, "five_ctauPVCA"};
   // TTreeReaderValue<Double_t> five_ctauErrPVCA = {fReader, "five_ctauErrPVCA"};
   // TTreeReaderValue<Double_t> five_cosAlphaDZ = {fReader, "five_cosAlphaDZ"};
   // TTreeReaderValue<Double_t> five_ctauPVDZ = {fReader, "five_ctauPVDZ"};
   // TTreeReaderValue<Double_t> five_ctauErrPVDZ = {fReader, "five_ctauErrPVDZ"};
   // TTreeReaderValue<Double_t> five_cosAlphaBS = {fReader, "five_cosAlphaBS"};
   // TTreeReaderValue<Double_t> five_ctauPVBS = {fReader, "five_ctauPVBS"};
   // TTreeReaderValue<Double_t> five_ctauErrPVBS = {fReader, "five_ctauErrPVBS"};
   // TTreeReaderValue<Double_t> five_vProb = {fReader, "five_vProb"};
   // TTreeReaderValue<Double_t> five_nDof = {fReader, "five_nDof"};
   // TTreeReaderValue<Double_t> five_vChi2 = {fReader, "five_vChi2"};
   // TTreeReaderValue<Double_t> five_vx = {fReader, "five_vx"};
   // TTreeReaderValue<Double_t> five_vy = {fReader, "five_vy"};
   // TTreeReaderValue<Double_t> five_vz = {fReader, "five_vz"};
   // TTreeReaderValue<Int_t> five_charge = {fReader, "five_charge"};
   // TTreeReaderValue<Double_t> bestPV_X = {fReader, "bestPV_X"};
   // TTreeReaderValue<Double_t> bestPV_Y = {fReader, "bestPV_Y"};
   // TTreeReaderValue<Double_t> bestPV_Z = {fReader, "bestPV_Z"};
   // TTreeReaderValue<Double_t> cosAlphaPV_X = {fReader, "cosAlphaPV_X"};
   // TTreeReaderValue<Double_t> cosAlphaPV_Y = {fReader, "cosAlphaPV_Y"};
   // TTreeReaderValue<Double_t> cosAlphaPV_Z = {fReader, "cosAlphaPV_Z"};
   // TTreeReaderValue<Double_t> bS_X = {fReader, "bS_X"};
   // TTreeReaderValue<Double_t> bS_Y = {fReader, "bS_Y"};
   // TTreeReaderValue<Double_t> bS_Z = {fReader, "bS_Z"};
   // TTreeReaderValue<Double_t> zPV_X = {fReader, "zPV_X"};
   // TTreeReaderValue<Double_t> zPV_Y = {fReader, "zPV_Y"};
   // TTreeReaderValue<Double_t> zPV_Z = {fReader, "zPV_Z"};

   // TTreeReaderValue<Bool_t> lowMuon_isTight = {fReader, "lowMuon_isTight"};
   // TTreeReaderValue<Bool_t> lowMuon_isLoose = {fReader, "lowMuon_isLoose"};
   // TTreeReaderValue<Bool_t> lowMuon_isSoft = {fReader, "lowMuon_isSoft"};
   // TTreeReaderValue<Bool_t> lowMuon_isMedium = {fReader, "lowMuon_isMedium"};
   // TTreeReaderValue<Bool_t> lowMuon_isHighPt = {fReader, "lowMuon_isHighPt"};
   // TTreeReaderValue<Bool_t> lowMuon_isTracker = {fReader, "lowMuon_isTracker"};
   // TTreeReaderValue<Bool_t> lowMuon_isGlobal = {fReader, "lowMuon_isGlobal"};
   // TTreeReaderValue<Int_t> lowMuon_NPixelHits = {fReader, "lowMuon_NPixelHits"};
   // TTreeReaderValue<Int_t> lowMuon_NStripHits = {fReader, "lowMuon_NStripHits"};
   // TTreeReaderValue<Int_t> lowMuon_NTrackhits = {fReader, "lowMuon_NTrackhits"};
   // TTreeReaderValue<Int_t> lowMuon_NBPixHits = {fReader, "lowMuon_NBPixHits"};
   // TTreeReaderValue<Int_t> lowMuon_NPixLayers = {fReader, "lowMuon_NPixLayers"};
   // TTreeReaderValue<Int_t> lowMuon_NTraLayers = {fReader, "lowMuon_NTraLayers"};
   // TTreeReaderValue<Int_t> lowMuon_NStrLayers = {fReader, "lowMuon_NStrLayers"};
   // TTreeReaderValue<Int_t> lowMuon_NBPixLayers = {fReader, "lowMuon_NBPixLayers"};
   // TTreeReaderValue<Bool_t> highMuon_isTight = {fReader, "highMuon_isTight"};
   // TTreeReaderValue<Bool_t> highMuon_isLoose = {fReader, "highMuon_isLoose"};
   // TTreeReaderValue<Bool_t> highMuon_isSoft = {fReader, "highMuon_isSoft"};
   // TTreeReaderValue<Bool_t> highMuon_isMedium = {fReader, "highMuon_isMedium"};
   // TTreeReaderValue<Bool_t> highMuon_isHighPt = {fReader, "highMuon_isHighPt"};
   // TTreeReaderValue<Bool_t> highMuon_isTracker = {fReader, "highMuon_isTracker"};
   // TTreeReaderValue<Bool_t> highMuon_isGlobal = {fReader, "highMuon_isGlobal"};
   // TTreeReaderValue<Int_t> highMuon_NPixelHits = {fReader, "highMuon_NPixelHits"};
   // TTreeReaderValue<Int_t> highMuon_NStripHits = {fReader, "highMuon_NStripHits"};
   // TTreeReaderValue<Int_t> highMuon_NTrackhits = {fReader, "highMuon_NTrackhits"};
   // TTreeReaderValue<Int_t> highMuon_NBPixHits = {fReader, "highMuon_NBPixHits"};
   // TTreeReaderValue<Int_t> highMuon_NPixLayers = {fReader, "highMuon_NPixLayers"};
   // TTreeReaderValue<Int_t> highMuon_NTraLayers = {fReader, "highMuon_NTraLayers"};
   // TTreeReaderValue<Int_t> highMuon_NStrLayers = {fReader, "highMuon_NStrLayers"};
   // TTreeReaderValue<Int_t> highMuon_NBPixLayers = {fReader, "highMuon_NBPixLayers"};
   // TTreeReaderValue<UInt_t> lowMuon_type = {fReader, "lowMuon_type"};
   // TTreeReaderValue<UInt_t> highMuon_type = {fReader, "highMuon_type"};
   // TTreeReaderValue<Int_t> highTrack_NPixelHits = {fReader, "highTrack_NPixelHits"};
   // TTreeReaderValue<Int_t> highTrack_NStripHits = {fReader, "highTrack_NStripHits"};
   // TTreeReaderValue<Int_t> highTrack_NTrackhits = {fReader, "highTrack_NTrackhits"};
   // TTreeReaderValue<Int_t> highTrack_NBPixHits = {fReader, "highTrack_NBPixHits"};
   // TTreeReaderValue<Int_t> highTrack_NPixLayers = {fReader, "highTrack_NPixLayers"};
   // TTreeReaderValue<Int_t> highTrack_NTraLayers = {fReader, "highTrack_NTraLayers"};
   // TTreeReaderValue<Int_t> highTrack_NStrLayers = {fReader, "highTrack_NStrLayers"};
   // TTreeReaderValue<Int_t> highTrack_NBPixLayers = {fReader, "highTrack_NBPixLayers"};
   // TTreeReaderValue<Int_t> lowTrack_NPixelHits = {fReader, "lowTrack_NPixelHits"};
   // TTreeReaderValue<Int_t> lowTrack_NStripHits = {fReader, "lowTrack_NStripHits"};
   // TTreeReaderValue<Int_t> lowTrack_NTrackhits = {fReader, "lowTrack_NTrackhits"};
   // TTreeReaderValue<Int_t> lowTrack_NBPixHits = {fReader, "lowTrack_NBPixHits"};
   // TTreeReaderValue<Int_t> lowTrack_NPixLayers = {fReader, "lowTrack_NPixLayers"};
   // TTreeReaderValue<Int_t> lowTrack_NTraLayers = {fReader, "lowTrack_NTraLayers"};
   // TTreeReaderValue<Int_t> lowTrack_NStrLayers = {fReader, "lowTrack_NStrLayers"};
   // TTreeReaderValue<Int_t> lowTrack_NBPixLayers = {fReader, "lowTrack_NBPixLayers"};
   // TTreeReaderValue<Int_t> thirdTrack_NPixelHits = {fReader, "thirdTrack_NPixelHits"};
   // TTreeReaderValue<Int_t> thirdTrack_NStripHits = {fReader, "thirdTrack_NStripHits"};
   // TTreeReaderValue<Int_t> thirdTrack_NTrackhits = {fReader, "thirdTrack_NTrackhits"};
   // TTreeReaderValue<Int_t> thirdTrack_NBPixHits = {fReader, "thirdTrack_NBPixHits"};
   // TTreeReaderValue<Int_t> thirdTrack_NPixLayers = {fReader, "thirdTrack_NPixLayers"};
   // TTreeReaderValue<Int_t> thirdTrack_NTraLayers = {fReader, "thirdTrack_NTraLayers"};
   // TTreeReaderValue<Int_t> thirdTrack_NStrLayers = {fReader, "thirdTrack_NStrLayers"};
   // TTreeReaderValue<Int_t> thirdTrack_NBPixLayers = {fReader, "thirdTrack_NBPixLayers"};
   // TTreeReaderValue<Int_t> fourthTrack_NPixLayers = {fReader, "fourthTrack_NPixLayers"};
   // TTreeReaderValue<Int_t> fourthTrack_NTraLayers = {fReader, "fourthTrack_NTraLayers"};
   // TTreeReaderValue<Int_t> fourthTrack_NStrLayers = {fReader, "fourthTrack_NStrLayers"};
   // TTreeReaderValue<Int_t> fourthTrack_NBPixLayers = {fReader, "fourthTrack_NBPixLayers"};
   // TTreeReaderValue<Int_t> fourthTrack_NPixelHits = {fReader, "fourthTrack_NPixelHits"};
   // TTreeReaderValue<Int_t> fourthTrack_NStripHits = {fReader, "fourthTrack_NStripHits"};
   // TTreeReaderValue<Int_t> fourthTrack_NTrackhits = {fReader, "fourthTrack_NTrackhits"};
   // TTreeReaderValue<Int_t> fourthTrack_NBPixHits = {fReader, "fourthTrack_NBPixHits"};
   // TTreeReaderValue<Double_t> six_m = {fReader, "six_m"};
   // TTreeReaderValue<Double_t> six_m_ref = {fReader, "six_m_ref"};
   // TTreeReaderValue<Double_t> six_mass_ppkk = {fReader, "six_mass_ppkk"};
   // TTreeReaderValue<Double_t> six_mass_pkpk = {fReader, "six_mass_pkpk"};
   // TTreeReaderValue<Double_t> six_mass_pkkk = {fReader, "six_mass_pkkk"};
   // TTreeReaderValue<Double_t> six_mass_kpkp = {fReader, "six_mass_kpkp"};
   // TTreeReaderValue<Double_t> six_mass_kppk = {fReader, "six_mass_kppk"};
   // TTreeReaderValue<Double_t> six_mass_kkkk = {fReader, "six_mass_kkkk"};
   // TTreeReaderValue<Double_t> six_pt = {fReader, "six_pt"};
   // TTreeReaderValue<Double_t> six_eta = {fReader, "six_eta"};
   // TTreeReaderValue<Double_t> six_phi = {fReader, "six_phi"};
   // TTreeReaderValue<Double_t> six_p = {fReader, "six_p"};
   // TTreeReaderValue<Double_t> six_cosAlpha = {fReader, "six_cosAlpha"};
   // TTreeReaderValue<Double_t> six_ctauPV = {fReader, "six_ctauPV"};
   // TTreeReaderValue<Double_t> six_ctauErrPV = {fReader, "six_ctauErrPV"};
   // TTreeReaderValue<Double_t> six_cosAlphaCA = {fReader, "six_cosAlphaCA"};
   // TTreeReaderValue<Double_t> six_ctauPVCA = {fReader, "six_ctauPVCA"};
   // TTreeReaderValue<Double_t> six_ctauErrPVCA = {fReader, "six_ctauErrPVCA"};
   // TTreeReaderValue<Double_t> six_cosAlphaDZ = {fReader, "six_cosAlphaDZ"};
   // TTreeReaderValue<Double_t> six_ctauPVDZ = {fReader, "six_ctauPVDZ"};
   // TTreeReaderValue<Double_t> six_ctauErrPVDZ = {fReader, "six_ctauErrPVDZ"};
   // TTreeReaderValue<Double_t> six_cosAlphaBS = {fReader, "six_cosAlphaBS"};
   // TTreeReaderValue<Double_t> six_ctauPVBS = {fReader, "six_ctauPVBS"};
   // TTreeReaderValue<Double_t> six_ctauErrPVBS = {fReader, "six_ctauErrPVBS"};
   // TTreeReaderValue<Double_t> six_vProb = {fReader, "six_vProb"};
   // TTreeReaderValue<Double_t> six_nDof = {fReader, "six_nDof"};
   // TTreeReaderValue<Double_t> six_vChi2 = {fReader, "six_vChi2"};
   // TTreeReaderValue<Double_t> six_vx = {fReader, "six_vx"};
   // TTreeReaderValue<Double_t> six_vy = {fReader, "six_vy"};
   // TTreeReaderValue<Double_t> six_vz = {fReader, "six_vz"};
   // TTreeReaderValue<Int_t> six_charge = {fReader, "six_charge"};

   Float_t out_run, out_event, out_lumi, out_numPrimaryVertices, out_trigger;
   Float_t out_noSixCandidates, out_five_id, out_dimuon_id, out_p_id, out_m_id;
   Float_t out_t_id, out_f_id;
   TLorentzVector out_six_p4, out_five_p4, out_dimuonditrk_p4;
   TLorentzVector out_ditrack_p4, out_dimuon_p4, out_lowMuon_p4, out_highMuon_p4, out_highKaon_p4;
   TLorentzVector out_lowKaon_p4, out_thirdKaon_p4, out_fourthKaon_p4, out_highPion_p4, out_lowPion_p4;
   TLorentzVector out_thirdPion_p4, out_fourthPion_p4, out_highProton_p4, out_lowProton_p4, out_thirdProton_p4;
   TLorentzVector out_fourthProton_p4;
   Float_t out_dimuonditrk_m, out_dimuonditrk_pt, out_dimuonditrk_eta, out_dimuonditrk_phi;
   Float_t out_dimuonditrk_p, out_dimuon_m, out_dimuon_pt, out_dimuon_eta, out_dimuon_phi;
   Float_t out_dimuon_p, out_highTrackMatch, out_lowTrackMatch, out_lowMuonMatch, out_highMuonMatch;
   Float_t out_thirdTrackMatch, out_fourthTrackMatch, out_ditrack_m, out_diTrackOne_pt, out_diTrackOne_eta;
   Float_t out_diTrackOne_phi, out_diTrackOne_p, out_diTrackTwo_pt, out_diTrackTwo_eta, out_diTrackTwo_phi;
   Float_t out_diTrackTwo_p, out_diTrackThree_pt, out_diTrackThree_eta, out_diTrackThree_phi, out_diTrackThree_p;
   Float_t out_diTrackFour_pt, out_diTrackFour_eta, out_diTrackFour_phi, out_diTrackFour_p, out_diTrackFive_pt;
   Float_t out_diTrackFive_eta, out_diTrackFive_phi, out_diTrackFive_p, out_diTrackSix_pt, out_diTrackSix_eta;
   Float_t out_diTrackSix_phi, out_diTrackSix_p, out_dimuonDiTrkOne_mmpp, out_dimuonDiTrkTwo_mmpp, out_dimuonDiTrkThree_mmpp;
   Float_t out_dimuonDiTrkFour_mmpp, out_dimuonDiTrkOne_mmkk, out_dimuonDiTrkTwo_mmkk, out_dimuonDiTrkThree_mmkk, out_dimuonDiTrkFour_mmkk;
   Float_t out_highMuon_pt, out_highMuon_eta, out_highMuon_phi, out_highMuon_charge, out_highMuon_dz;
   Float_t out_highMuon_dxy, out_lowMuon_pt, out_lowMuon_eta, out_lowMuon_phi, out_lowMuon_charge;
   Float_t out_lowMuon_dz, out_lowMuon_dxy, out_highTrack_pt, out_highTrack_eta, out_highTrack_phi;
   Float_t out_highTrack_charge, out_highTrack_dz, out_highTrack_dxy, out_lowTrack_pt, out_lowTrack_eta;
   Float_t out_lowTrack_phi, out_lowTrack_charge, out_lowTrack_dz, out_lowTrack_dxy, out_thirdTrack_pt;
   Float_t out_thirdTrack_eta, out_thirdTrack_phi, out_thirdTrack_charge, out_thirdTrack_dz, out_thirdTrack_dxy;
   Float_t out_dimuonDiTrkOne_pt, out_dimuonDiTrkOne_eta, out_dimuonDiTrkOne_phi, out_dimuonDiTrkOne_charge, out_dimuonDiTrkOne_p;
   Float_t out_dimuonDiTrkTwo_pt, out_dimuonDiTrkTwo_eta, out_dimuonDiTrkTwo_phi, out_dimuonDiTrkTwo_charge, out_dimuonDiTrkTwo_p;
   Float_t out_dimuonDiTrkThree_pt, out_dimuonDiTrkThree_eta, out_dimuonDiTrkThree_phi, out_dimuonDiTrkThree_charge, out_dimuonDiTrkThree_p;
   Float_t out_dimuonDiTrkFour_pt, out_dimuonDiTrkFour_eta, out_dimuonDiTrkFour_phi, out_dimuonDiTrkFour_charge, out_dimuonDiTrkFour_p;
   Float_t out_dimuonDiTrkFive_pt, out_dimuonDiTrkFive_eta, out_dimuonDiTrkFive_phi, out_dimuonDiTrkFive_charge, out_dimuonDiTrkFive_p;
   Float_t out_dimuonDiTrkSix_pt, out_dimuonDiTrkSix_eta, out_dimuonDiTrkSix_phi, out_dimuonDiTrkSix_charge, out_dimuonDiTrkSix_p;
   Float_t out_dimuon_vProb, out_dimuon_vChi2, out_dimuon_DCA, out_dimuon_ctauPV, out_dimuon_ctauErrPV;
   Float_t out_dimuon_cosAlpha, out_triTrack_m, out_triTrack_pt, out_triTrack_eta, out_triTrack_phi;
   Float_t out_triTrack_charge, out_dimuonditrk_vProb, out_dimuonditrk_vChi2, out_dimuonditrk_nDof, out_dimuonditrk_charge;
   Float_t out_dimuonditrk_cosAlpha, out_dimuonditrk_ctauPV, out_dimuonditrk_ctauErrPV, out_dimuonditrk_cosAlphaCA, out_dimuonditrk_ctauPVCA;
   Float_t out_dimuonditrk_ctauErrPVCA, out_dimuonditrk_cosAlphaDZ, out_dimuonditrk_ctauPVDZ, out_dimuonditrk_ctauErrPVDZ, out_dimuonditrk_cosAlphaBS;
   Float_t out_dimuonditrk_ctauPVBS, out_dimuonditrk_ctauErrPVBS, out_dimuonditrk_vx, out_dimuonditrk_vy, out_dimuonditrk_vz;
   Float_t out_dca_m1m2, out_dca_m1t1, out_dca_m1t2, out_dca_m2t1, out_dca_m2t2;
   Float_t out_dca_t1t2, out_dca_m1t3, out_dca_m2t3, out_dca_t1t3, out_dca_t2t3;
   Float_t out_dca_m1t4, out_dca_m2t4, out_dca_t1t4, out_dca_t2t4, out_dca_t3t4;
   Float_t out_highTrackMuonDR, out_highTrackMuonDP, out_highTrackMuonDPt, out_lowTrackMuonDR, out_lowTrackMuonDP;
   Float_t out_lowTrackMuonDPt, out_thirdTrackMuonDR, out_thirdTrackMuonDP, out_thirdTrackMuonDPt, out_fourthTrackMuonDR;
   Float_t out_fourthTrackMuonDP, out_fourthTrackMuonDPt, out_tPFromPV, out_tMFromPV, out_tTFromPV;
   Float_t out_tFFromPV, out_tPFromPVCA, out_tMFromPVCA, out_tTFromPVCA, out_tFFromPVCA;
   Float_t out_tPFromPVDZ, out_tMFromPVDZ, out_tTFromPVDZ, out_tFFromPVDZ, out_five_m;
   Float_t out_five_m_ref, out_five_mass_ppk, out_five_mass_kpp, out_five_mass_pkp, out_five_mass_ppp;
   Float_t out_fiveOne_pt, out_fiveOne_eta, out_fiveOne_phi, out_fiveOne_p, out_fiveTwo_pt;
   Float_t out_fiveTwo_eta, out_fiveTwo_phi, out_fiveTwo_p, out_fiveThree_pt, out_fiveThree_eta;
   Float_t out_fiveThree_phi, out_fiveThree_p, out_fiveFour_pt, out_fiveFour_eta, out_fiveFour_phi;
   Float_t out_fiveFour_p, out_fiveFive_pt, out_fiveFive_eta, out_fiveFive_phi, out_fiveFive_p;
   Float_t out_five_cosAlpha, out_five_ctauPV, out_five_ctauErrPV, out_five_cosAlphaCA, out_five_ctauPVCA;
   Float_t out_five_ctauErrPVCA, out_five_cosAlphaDZ, out_five_ctauPVDZ, out_five_ctauErrPVDZ, out_five_cosAlphaBS;
   Float_t out_five_ctauPVBS, out_five_ctauErrPVBS, out_five_vProb, out_five_nDof, out_five_vChi2;
   Float_t out_five_vx, out_five_vy, out_five_vz, out_five_charge, out_bestPV_X;
   Float_t out_bestPV_Y, out_bestPV_Z, out_cosAlphaPV_X, out_cosAlphaPV_Y, out_cosAlphaPV_Z;
   Float_t out_bS_X, out_bS_Y, out_bS_Z, out_zPV_X, out_zPV_Y;
   Float_t out_zPV_Z, out_lowMuon_isTight, out_lowMuon_isLoose, out_lowMuon_isSoft, out_lowMuon_isMedium;
   Float_t out_lowMuon_isHighPt, out_lowMuon_isTracker, out_lowMuon_isGlobal, out_lowMuon_NPixelHits, out_lowMuon_NStripHits;
   Float_t out_lowMuon_NTrackhits, out_lowMuon_NBPixHits, out_lowMuon_NPixLayers, out_lowMuon_NTraLayers, out_lowMuon_NStrLayers;
   Float_t out_lowMuon_NBPixLayers, out_highMuon_isTight, out_highMuon_isLoose, out_highMuon_isSoft, out_highMuon_isMedium;
   Float_t out_highMuon_isHighPt, out_highMuon_isTracker, out_highMuon_isGlobal, out_highMuon_NPixelHits, out_highMuon_NStripHits;
   Float_t out_highMuon_NTrackhits, out_highMuon_NBPixHits, out_highMuon_NPixLayers, out_highMuon_NTraLayers, out_highMuon_NStrLayers;
   Float_t out_highMuon_NBPixLayers, out_lowMuon_type, out_highMuon_type, out_highTrack_NPixelHits, out_highTrack_NStripHits;
   Float_t out_highTrack_NTrackhits, out_highTrack_NBPixHits, out_highTrack_NPixLayers, out_highTrack_NTraLayers, out_highTrack_NStrLayers;
   Float_t out_highTrack_NBPixLayers, out_lowTrack_NPixelHits, out_lowTrack_NStripHits, out_lowTrack_NTrackhits, out_lowTrack_NBPixHits;
   Float_t out_lowTrack_NPixLayers, out_lowTrack_NTraLayers, out_lowTrack_NStrLayers, out_lowTrack_NBPixLayers, out_thirdTrack_NPixelHits;
   Float_t out_thirdTrack_NStripHits, out_thirdTrack_NTrackhits, out_thirdTrack_NBPixHits, out_thirdTrack_NPixLayers, out_thirdTrack_NTraLayers;
   Float_t out_thirdTrack_NStrLayers, out_thirdTrack_NBPixLayers, out_fourthTrack_NPixLayers, out_fourthTrack_NTraLayers, out_fourthTrack_NStrLayers;
   Float_t out_fourthTrack_NBPixLayers, out_fourthTrack_NPixelHits, out_fourthTrack_NStripHits, out_fourthTrack_NTrackhits, out_fourthTrack_NBPixHits;
   Float_t out_six_m, out_six_m_ref, out_six_mass_ppkk, out_six_mass_pkpk, out_six_mass_pkkk;
   Float_t out_six_mass_kpkp, out_six_mass_kppk, out_six_mass_kkkk, out_six_pt, out_six_eta;
   Float_t out_six_phi, out_six_p, out_six_cosAlpha, out_six_ctauPV, out_six_ctauErrPV;
   Float_t out_six_cosAlphaCA, out_six_ctauPVCA, out_six_ctauErrPVCA, out_six_cosAlphaDZ, out_six_ctauPVDZ;
   Float_t out_six_ctauErrPVDZ, out_six_cosAlphaBS, out_six_ctauPVBS, out_six_ctauErrPVBS, out_six_vProb;
   Float_t out_six_nDof, out_six_vChi2, out_six_vx, out_six_vy, out_six_vz;
   Float_t out_six_charge;


   SixTracks(TTree * /*tree*/ =0) { }
   virtual ~SixTracks() { }
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

   TProofOutputFile *OutFile;
   TFile            *fOut;

   ClassDef(SixTracks,0);

};

#endif

#ifdef SixTracks_cxx
void SixTracks::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t SixTracks::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef SixTracks_cxx
