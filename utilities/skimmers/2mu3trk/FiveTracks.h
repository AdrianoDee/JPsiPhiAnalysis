//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep 13 02:14:52 2018 by ROOT version 6.12/07
// from TTree JPsiPhiTree/Tree of DiMuon and DiTrak
// found on file: ../../../../../CMSSW_10_2_1/src/jpsiphi/jpsiphi/test/rootuple-2018-dimuonditrak_bbbar_hard_0.root
//////////////////////////////////////////////////////////

#ifndef FiveTracks_h
#define FiveTracks_h

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



class FiveTracks : public TSelector {
  public :
  TTreeReader     fReader;  //!the tree reader
  TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

  Float_t JPsi_mass, Phi_mass, Phi_mean, Phi_sigma;
  TTree *outTree;

  //Double_t out;

  // Readers to access the data (delete the ones you do not need).

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> run = {fReader, "run"};
   TTreeReaderValue<Int_t> event = {fReader, "event"};
   TTreeReaderValue<Int_t> lumi = {fReader, "lumi"};
   TTreeReaderValue<Int_t> numPrimaryVertices = {fReader, "numPrimaryVertices"};
   TTreeReaderValue<Int_t> trigger = {fReader, "trigger"};
   TTreeReaderValue<Int_t> noFiveCandidates = {fReader, "noFiveCandidates"};
   TTreeReaderValue<Int_t> dimuonditrk_id = {fReader, "dimuonditrk_id"};
   TTreeReaderValue<TLorentzVector> dimuonditrk_p4 = {fReader, "dimuonditrk_p4"};
   TTreeReaderValue<TLorentzVector> ditrak_p4 = {fReader, "ditrak_p4"};
   TTreeReaderValue<TLorentzVector> dimuon_p4 = {fReader, "dimuon_p4"};
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
   TTreeReaderValue<Double_t> highKaonMatch = {fReader, "highKaonMatch"};
   TTreeReaderValue<Double_t> lowKaonMatch = {fReader, "lowKaonMatch"};
   TTreeReaderValue<Double_t> lowMuonMatch = {fReader, "lowMuonMatch"};
   TTreeReaderValue<Double_t> highMuonMatch = {fReader, "highMuonMatch"};
   TTreeReaderValue<Double_t> ditrak_m = {fReader, "ditrak_m"};
   TTreeReaderValue<Double_t> ditrakOne_pt = {fReader, "ditrakOne_pt"};
   TTreeReaderValue<Double_t> ditrakOne_eta = {fReader, "ditrakOne_eta"};
   TTreeReaderValue<Double_t> ditrakOne_phi = {fReader, "ditrakOne_phi"};
   TTreeReaderValue<Double_t> ditrakOne_p = {fReader, "ditrakOne_p"};
   TTreeReaderValue<Double_t> ditrakTwo_pt = {fReader, "ditrakTwo_pt"};
   TTreeReaderValue<Double_t> ditrakTwo_eta = {fReader, "ditrakTwo_eta"};
   TTreeReaderValue<Double_t> ditrakTwo_phi = {fReader, "ditrakTwo_phi"};
   TTreeReaderValue<Double_t> ditrakTwo_p = {fReader, "ditrakTwo_p"};
   TTreeReaderValue<Double_t> ditrakThree_pt = {fReader, "ditrakThree_pt"};
   TTreeReaderValue<Double_t> ditrakThree_eta = {fReader, "ditrakThree_eta"};
   TTreeReaderValue<Double_t> ditrakThree_phi = {fReader, "ditrakThree_phi"};
   TTreeReaderValue<Double_t> ditrakThree_p = {fReader, "ditrakThree_p"};
   TTreeReaderValue<Double_t> highTrack_pt = {fReader, "highTrack_pt"};
   TTreeReaderValue<Double_t> highTrack_eta = {fReader, "highTrack_eta"};
   TTreeReaderValue<Double_t> highTrack_phi = {fReader, "highTrack_phi"};
   TTreeReaderValue<Double_t> highTrack_charge = {fReader, "highTrack_charge"};
   TTreeReaderValue<Double_t> lowTrack_pt = {fReader, "lowTrack_pt"};
   TTreeReaderValue<Double_t> lowTrack_eta = {fReader, "lowTrack_eta"};
   TTreeReaderValue<Double_t> lowTrack_phi = {fReader, "lowTrack_phi"};
   TTreeReaderValue<Double_t> lowTrack_charge = {fReader, "lowTrack_charge"};
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
   TTreeReaderValue<Double_t> dimuonDiTrkTwo_pt = {fReader, "dimuonDiTrkTwo_pt"};
   TTreeReaderValue<Double_t> dimuonDiTrkTwo_eta = {fReader, "dimuonDiTrkTwo_eta"};
   TTreeReaderValue<Double_t> dimuonDiTrkTwo_phi = {fReader, "dimuonDiTrkTwo_phi"};
   TTreeReaderValue<Double_t> dimuonDiTrkTwo_charge = {fReader, "dimuonDiTrkTwo_charge"};
   TTreeReaderValue<Double_t> dimuonDiTrkThree_pt = {fReader, "dimuonDiTrkThree_pt"};
   TTreeReaderValue<Double_t> dimuonDiTrkThree_eta = {fReader, "dimuonDiTrkThree_eta"};
   TTreeReaderValue<Double_t> dimuonDiTrkThree_phi = {fReader, "dimuonDiTrkThree_phi"};
   TTreeReaderValue<Double_t> dimuonDiTrkThree_charge = {fReader, "dimuonDiTrkThree_charge"};
   TTreeReaderValue<Double_t> psiPrimeSame_pt = {fReader, "psiPrimeSame_pt"};
   TTreeReaderValue<Double_t> psiPrimeSame_eta = {fReader, "psiPrimeSame_eta"};
   TTreeReaderValue<Double_t> psiPrimeSame_phi = {fReader, "psiPrimeSame_phi"};
   TTreeReaderValue<Double_t> psiPrimeSame_n = {fReader, "psiPrimeSame_n"};
   TTreeReaderValue<Double_t> psiPrimeSame_p_pt = {fReader, "psiPrimeSame_p_pt"};
   TTreeReaderValue<Double_t> psiPrimeSame_p_eta = {fReader, "psiPrimeSame_p_eta"};
   TTreeReaderValue<Double_t> psiPrimeSame_p_phi = {fReader, "psiPrimeSame_p_phi"};
   TTreeReaderValue<Double_t> psiPrimeSame_p_n = {fReader, "psiPrimeSame_p_n"};
   TTreeReaderValue<Double_t> psiPrimeSame_m_pt = {fReader, "psiPrimeSame_m_pt"};
   TTreeReaderValue<Double_t> psiPrimeSame_m_eta = {fReader, "psiPrimeSame_m_eta"};
   TTreeReaderValue<Double_t> psiPrimeSame_m_phi = {fReader, "psiPrimeSame_m_phi"};
   TTreeReaderValue<Double_t> psiPrimeSame_m_n = {fReader, "psiPrimeSame_m_n"};
   TTreeReaderValue<Double_t> psiPrimeMixed_pt = {fReader, "psiPrimeMixed_pt"};
   TTreeReaderValue<Double_t> psiPrimeMixed_eta = {fReader, "psiPrimeMixed_eta"};
   TTreeReaderValue<Double_t> psiPrimeMixed_phi = {fReader, "psiPrimeMixed_phi"};
   TTreeReaderValue<Double_t> psiPrimeMixed_n = {fReader, "psiPrimeMixed_n"};
   TTreeReaderValue<Double_t> psiPrimeMixed_p_pt = {fReader, "psiPrimeMixed_p_pt"};
   TTreeReaderValue<Double_t> psiPrimeMixed_p_eta = {fReader, "psiPrimeMixed_p_eta"};
   TTreeReaderValue<Double_t> psiPrimeMixed_p_phi = {fReader, "psiPrimeMixed_p_phi"};
   TTreeReaderValue<Double_t> psiPrimeMixed_p_n = {fReader, "psiPrimeMixed_p_n"};
   TTreeReaderValue<Double_t> psiPrimeMixed_m_pt = {fReader, "psiPrimeMixed_m_pt"};
   TTreeReaderValue<Double_t> psiPrimeMixed_m_eta = {fReader, "psiPrimeMixed_m_eta"};
   TTreeReaderValue<Double_t> psiPrimeMixed_m_phi = {fReader, "psiPrimeMixed_m_phi"};
   TTreeReaderValue<Double_t> psiPrimeMixed_m_n = {fReader, "psiPrimeMixed_m_n"};
   TTreeReaderValue<Double_t> psiPrimeSame_ditrak_pt = {fReader, "psiPrimeSame_ditrak_pt"};
   TTreeReaderValue<Double_t> psiPrimeSame_ditrak_eta = {fReader, "psiPrimeSame_ditrak_eta"};
   TTreeReaderValue<Double_t> psiPrimeSame_ditrak_phi = {fReader, "psiPrimeSame_ditrak_phi"};
   TTreeReaderValue<Double_t> psiPrimeSame_ditrak_n = {fReader, "psiPrimeSame_ditrak_n"};
   TTreeReaderValue<Double_t> psiPrimeMixed_ditrak_pt = {fReader, "psiPrimeMixed_ditrak_pt"};
   TTreeReaderValue<Double_t> psiPrimeMixed_ditrak_eta = {fReader, "psiPrimeMixed_ditrak_eta"};
   TTreeReaderValue<Double_t> psiPrimeMixed_ditrak_phi = {fReader, "psiPrimeMixed_ditrak_phi"};
   TTreeReaderValue<Double_t> psiPrimeMixed_ditrak_n = {fReader, "psiPrimeMixed_ditrak_n"};
   TTreeReaderValue<Double_t> triTrak_pt = {fReader, "triTrak_pt"};
   TTreeReaderValue<Double_t> triTrak_eta = {fReader, "triTrak_eta"};
   TTreeReaderValue<Double_t> triTrak_phi = {fReader, "triTrak_phi"};
   TTreeReaderValue<Double_t> triTrak_charge = {fReader, "triTrak_charge"};
   TTreeReaderValue<Double_t> mass_kkk = {fReader, "mass_kkk"};
   TTreeReaderValue<Double_t> mass_ref_kkk = {fReader, "mass_ref_kkk"};
   TTreeReaderValue<Double_t> vProb_kkk = {fReader, "vProb_kkk"};
   TTreeReaderValue<Double_t> nDof_kkk = {fReader, "nDof_kkk"};
   TTreeReaderValue<Double_t> vChi2_kkk = {fReader, "vChi2_kkk"};
   TTreeReaderValue<Double_t> ctau_kkk = {fReader, "ctau_kkk"};
   TTreeReaderValue<Double_t> ctauErr_kkk = {fReader, "ctauErr_kkk"};
   TTreeReaderValue<Double_t> cosAlpha_kkk = {fReader, "cosAlpha_kkk"};
   TTreeReaderValue<Double_t> onePsiPrime_m_kkk = {fReader, "onePsiPrime_m_kkk"};
   TTreeReaderValue<Double_t> twoPsiPrime_m_kkk = {fReader, "twoPsiPrime_m_kkk"};
   TTreeReaderValue<Double_t> onePsiPrime_p_mkkk = {fReader, "onePsiPrime_p_mkkk"};
   TTreeReaderValue<Double_t> onePsiPrime_m_mkkk = {fReader, "onePsiPrime_m_mkkk"};
   TTreeReaderValue<Double_t> twoPsiPrime_p_mkkk = {fReader, "twoPsiPrime_p_mkkk"};
   TTreeReaderValue<Double_t> twoPsiPrime_m_mkkk = {fReader, "twoPsiPrime_m_mkkk"};
   TTreeReaderValue<Double_t> dimuonDiTrkOne_m_kkk = {fReader, "dimuonDiTrkOne_m_kkk"};
   TTreeReaderValue<Double_t> dimuonDiTrkTwo_m_kkk = {fReader, "dimuonDiTrkTwo_m_kkk"};
   TTreeReaderValue<Double_t> dimuonDiTrkThree_m_kkk = {fReader, "dimuonDiTrkThree_m_kkk"};
   TTreeReaderValue<Double_t> ditrakOne_m_kkk = {fReader, "ditrakOne_m_kkk"};
   TTreeReaderValue<Double_t> ditrakTwo_m_kkk = {fReader, "ditrakTwo_m_kkk"};
   TTreeReaderValue<Double_t> ditrakThree_m_kkk = {fReader, "ditrakThree_m_kkk"};
   TTreeReaderValue<Double_t> trackOne_m_kkk = {fReader, "trackOne_m_kkk"};
   TTreeReaderValue<Double_t> trackTwo_m_kkk = {fReader, "trackTwo_m_kkk"};
   TTreeReaderValue<Double_t> trackThree_m_kkk = {fReader, "trackThree_m_kkk"};
   TTreeReaderValue<TLorentzVector> five_p4_kkk = {fReader, "five_p4_kkk"};
   TTreeReaderValue<TLorentzVector> dimuonDiTrkOne_p4_kkk = {fReader, "dimuonDiTrkOne_p4_kkk"};
   TTreeReaderValue<TLorentzVector> dimuonDiTrkTwo_p4_kkk = {fReader, "dimuonDiTrkTwo_p4_kkk"};
   TTreeReaderValue<TLorentzVector> dimuonDiTrkThree_p4_kkk = {fReader, "dimuonDiTrkThree_p4_kkk"};
   TTreeReaderValue<TLorentzVector> ditrakOne_p4_kkk = {fReader, "ditrakOne_p4_kkk"};
   TTreeReaderValue<TLorentzVector> ditrakTwo_p4_kkk = {fReader, "ditrakTwo_p4_kkk"};
   TTreeReaderValue<TLorentzVector> ditrakThree_p4_kkk = {fReader, "ditrakThree_p4_kkk"};
   TTreeReaderValue<TLorentzVector> psiPrimeSame_p4_kkk = {fReader, "psiPrimeSame_p4_kkk"};
   TTreeReaderValue<TLorentzVector> psiPrimeMixed_p4_kkk = {fReader, "psiPrimeMixed_p4_kkk"};
   TTreeReaderValue<TLorentzVector> psiPrimeSame_ditrak_p4_kkk = {fReader, "psiPrimeSame_ditrak_p4_kkk"};
   TTreeReaderValue<TLorentzVector> psiPrimeMixed_ditrak_p4_kkk = {fReader, "psiPrimeMixed_ditrak_p4_kkk"};
   TTreeReaderValue<Double_t> triTrak_m_kkk = {fReader, "triTrak_m_kkk"};
   TTreeReaderValue<Double_t> mass_ppk = {fReader, "mass_ppk"};
   TTreeReaderValue<Double_t> mass_ref_ppk = {fReader, "mass_ref_ppk"};
   TTreeReaderValue<Double_t> vProb_ppk = {fReader, "vProb_ppk"};
   TTreeReaderValue<Double_t> nDof_ppk = {fReader, "nDof_ppk"};
   TTreeReaderValue<Double_t> vChi2_ppk = {fReader, "vChi2_ppk"};
   TTreeReaderValue<Double_t> ctau_ppk = {fReader, "ctau_ppk"};
   TTreeReaderValue<Double_t> ctauErr_ppk = {fReader, "ctauErr_ppk"};
   TTreeReaderValue<Double_t> cosAlpha_ppk = {fReader, "cosAlpha_ppk"};
   TTreeReaderValue<Double_t> onePsiPrime_m_ppk = {fReader, "onePsiPrime_m_ppk"};
   TTreeReaderValue<Double_t> twoPsiPrime_m_ppk = {fReader, "twoPsiPrime_m_ppk"};
   TTreeReaderValue<Double_t> onePsiPrime_p_mppk = {fReader, "onePsiPrime_p_mppk"};
   TTreeReaderValue<Double_t> onePsiPrime_m_mppk = {fReader, "onePsiPrime_m_mppk"};
   TTreeReaderValue<Double_t> twoPsiPrime_p_mppk = {fReader, "twoPsiPrime_p_mppk"};
   TTreeReaderValue<Double_t> twoPsiPrime_m_mppk = {fReader, "twoPsiPrime_m_mppk"};
   TTreeReaderValue<Double_t> dimuonDiTrkOne_m_ppk = {fReader, "dimuonDiTrkOne_m_ppk"};
   TTreeReaderValue<Double_t> dimuonDiTrkTwo_m_ppk = {fReader, "dimuonDiTrkTwo_m_ppk"};
   TTreeReaderValue<Double_t> dimuonDiTrkThree_m_ppk = {fReader, "dimuonDiTrkThree_m_ppk"};
   TTreeReaderValue<Double_t> ditrakOne_m_ppk = {fReader, "ditrakOne_m_ppk"};
   TTreeReaderValue<Double_t> ditrakTwo_m_ppk = {fReader, "ditrakTwo_m_ppk"};
   TTreeReaderValue<Double_t> ditrakThree_m_ppk = {fReader, "ditrakThree_m_ppk"};
   TTreeReaderValue<Double_t> trackOne_m_ppk = {fReader, "trackOne_m_ppk"};
   TTreeReaderValue<Double_t> trackTwo_m_ppk = {fReader, "trackTwo_m_ppk"};
   TTreeReaderValue<Double_t> trackThree_m_ppk = {fReader, "trackThree_m_ppk"};
   TTreeReaderValue<TLorentzVector> five_p4_ppk = {fReader, "five_p4_ppk"};
   TTreeReaderValue<TLorentzVector> dimuonDiTrkOne_p4_ppk = {fReader, "dimuonDiTrkOne_p4_ppk"};
   TTreeReaderValue<TLorentzVector> dimuonDiTrkTwo_p4_ppk = {fReader, "dimuonDiTrkTwo_p4_ppk"};
   TTreeReaderValue<TLorentzVector> dimuonDiTrkThree_p4_ppk = {fReader, "dimuonDiTrkThree_p4_ppk"};
   TTreeReaderValue<TLorentzVector> ditrakOne_p4_ppk = {fReader, "ditrakOne_p4_ppk"};
   TTreeReaderValue<TLorentzVector> ditrakTwo_p4_ppk = {fReader, "ditrakTwo_p4_ppk"};
   TTreeReaderValue<TLorentzVector> ditrakThree_p4_ppk = {fReader, "ditrakThree_p4_ppk"};
   TTreeReaderValue<TLorentzVector> psiPrimeSame_p4_ppk = {fReader, "psiPrimeSame_p4_ppk"};
   TTreeReaderValue<TLorentzVector> psiPrimeMixed_p4_ppk = {fReader, "psiPrimeMixed_p4_ppk"};
   TTreeReaderValue<TLorentzVector> psiPrimeSame_ditrak_p4_ppk = {fReader, "psiPrimeSame_ditrak_p4_ppk"};
   TTreeReaderValue<TLorentzVector> psiPrimeMixed_ditrak_p4_ppk = {fReader, "psiPrimeMixed_ditrak_p4_ppk"};
   TTreeReaderValue<Double_t> triTrak_m_ppk = {fReader, "triTrak_m_ppk"};
   TTreeReaderValue<Double_t> mass_kpp = {fReader, "mass_kpp"};
   TTreeReaderValue<Double_t> mass_ref_kpp = {fReader, "mass_ref_kpp"};
   TTreeReaderValue<Double_t> vProb_kpp = {fReader, "vProb_kpp"};
   TTreeReaderValue<Double_t> nDof_kpp = {fReader, "nDof_kpp"};
   TTreeReaderValue<Double_t> vChi2_kpp = {fReader, "vChi2_kpp"};
   TTreeReaderValue<Double_t> ctau_kpp = {fReader, "ctau_kpp"};
   TTreeReaderValue<Double_t> ctauErr_kpp = {fReader, "ctauErr_kpp"};
   TTreeReaderValue<Double_t> cosAlpha_kpp = {fReader, "cosAlpha_kpp"};
   TTreeReaderValue<Double_t> onePsiPrime_m_kpp = {fReader, "onePsiPrime_m_kpp"};
   TTreeReaderValue<Double_t> twoPsiPrime_m_kpp = {fReader, "twoPsiPrime_m_kpp"};
   TTreeReaderValue<Double_t> onePsiPrime_p_mkpp = {fReader, "onePsiPrime_p_mkpp"};
   TTreeReaderValue<Double_t> onePsiPrime_m_mkpp = {fReader, "onePsiPrime_m_mkpp"};
   TTreeReaderValue<Double_t> twoPsiPrime_p_mkpp = {fReader, "twoPsiPrime_p_mkpp"};
   TTreeReaderValue<Double_t> twoPsiPrime_m_mkpp = {fReader, "twoPsiPrime_m_mkpp"};
   TTreeReaderValue<Double_t> dimuonDiTrkOne_m_kpp = {fReader, "dimuonDiTrkOne_m_kpp"};
   TTreeReaderValue<Double_t> dimuonDiTrkTwo_m_kpp = {fReader, "dimuonDiTrkTwo_m_kpp"};
   TTreeReaderValue<Double_t> dimuonDiTrkThree_m_kpp = {fReader, "dimuonDiTrkThree_m_kpp"};
   TTreeReaderValue<Double_t> ditrakOne_m_kpp = {fReader, "ditrakOne_m_kpp"};
   TTreeReaderValue<Double_t> ditrakTwo_m_kpp = {fReader, "ditrakTwo_m_kpp"};
   TTreeReaderValue<Double_t> ditrakThree_m_kpp = {fReader, "ditrakThree_m_kpp"};
   TTreeReaderValue<Double_t> trackOne_m_kpp = {fReader, "trackOne_m_kpp"};
   TTreeReaderValue<Double_t> trackTwo_m_kpp = {fReader, "trackTwo_m_kpp"};
   TTreeReaderValue<Double_t> trackThree_m_kpp = {fReader, "trackThree_m_kpp"};
   TTreeReaderValue<TLorentzVector> five_p4_kpp = {fReader, "five_p4_kpp"};
   TTreeReaderValue<TLorentzVector> dimuonDiTrkOne_p4_kpp = {fReader, "dimuonDiTrkOne_p4_kpp"};
   TTreeReaderValue<TLorentzVector> dimuonDiTrkTwo_p4_kpp = {fReader, "dimuonDiTrkTwo_p4_kpp"};
   TTreeReaderValue<TLorentzVector> dimuonDiTrkThree_p4_kpp = {fReader, "dimuonDiTrkThree_p4_kpp"};
   TTreeReaderValue<TLorentzVector> ditrakOne_p4_kpp = {fReader, "ditrakOne_p4_kpp"};
   TTreeReaderValue<TLorentzVector> ditrakTwo_p4_kpp = {fReader, "ditrakTwo_p4_kpp"};
   TTreeReaderValue<TLorentzVector> ditrakThree_p4_kpp = {fReader, "ditrakThree_p4_kpp"};
   TTreeReaderValue<TLorentzVector> psiPrimeSame_p4_kpp = {fReader, "psiPrimeSame_p4_kpp"};
   TTreeReaderValue<TLorentzVector> psiPrimeMixed_p4_kpp = {fReader, "psiPrimeMixed_p4_kpp"};
   TTreeReaderValue<TLorentzVector> psiPrimeSame_ditrak_p4_kpp = {fReader, "psiPrimeSame_ditrak_p4_kpp"};
   TTreeReaderValue<TLorentzVector> psiPrimeMixed_ditrak_p4_kpp = {fReader, "psiPrimeMixed_ditrak_p4_kpp"};
   TTreeReaderValue<Double_t> triTrak_m_kpp = {fReader, "triTrak_m_kpp"};
   TTreeReaderValue<Double_t> mass_pkp = {fReader, "mass_pkp"};
   TTreeReaderValue<Double_t> mass_ref_pkp = {fReader, "mass_ref_pkp"};
   TTreeReaderValue<Double_t> vProb_pkp = {fReader, "vProb_pkp"};
   TTreeReaderValue<Double_t> nDof_pkp = {fReader, "nDof_pkp"};
   TTreeReaderValue<Double_t> vChi2_pkp = {fReader, "vChi2_pkp"};
   TTreeReaderValue<Double_t> ctau_pkp = {fReader, "ctau_pkp"};
   TTreeReaderValue<Double_t> ctauErr_pkp = {fReader, "ctauErr_pkp"};
   TTreeReaderValue<Double_t> cosAlpha_pkp = {fReader, "cosAlpha_pkp"};
   TTreeReaderValue<Double_t> onePsiPrime_m_pkp = {fReader, "onePsiPrime_m_pkp"};
   TTreeReaderValue<Double_t> twoPsiPrime_m_pkp = {fReader, "twoPsiPrime_m_pkp"};
   TTreeReaderValue<Double_t> onePsiPrime_p_mpkp = {fReader, "onePsiPrime_p_mpkp"};
   TTreeReaderValue<Double_t> onePsiPrime_m_mpkp = {fReader, "onePsiPrime_m_mpkp"};
   TTreeReaderValue<Double_t> twoPsiPrime_p_mpkp = {fReader, "twoPsiPrime_p_mpkp"};
   TTreeReaderValue<Double_t> twoPsiPrime_m_mpkp = {fReader, "twoPsiPrime_m_mpkp"};
   TTreeReaderValue<Double_t> dimuonDiTrkOne_m_pkp = {fReader, "dimuonDiTrkOne_m_pkp"};
   TTreeReaderValue<Double_t> dimuonDiTrkTwo_m_pkp = {fReader, "dimuonDiTrkTwo_m_pkp"};
   TTreeReaderValue<Double_t> dimuonDiTrkThree_m_pkp = {fReader, "dimuonDiTrkThree_m_pkp"};
   TTreeReaderValue<Double_t> ditrakOne_m_pkp = {fReader, "ditrakOne_m_pkp"};
   TTreeReaderValue<Double_t> ditrakTwo_m_pkp = {fReader, "ditrakTwo_m_pkp"};
   TTreeReaderValue<Double_t> ditrakThree_m_pkp = {fReader, "ditrakThree_m_pkp"};
   TTreeReaderValue<Double_t> trackOne_m_pkp = {fReader, "trackOne_m_pkp"};
   TTreeReaderValue<Double_t> trackTwo_m_pkp = {fReader, "trackTwo_m_pkp"};
   TTreeReaderValue<Double_t> trackThree_m_pkp = {fReader, "trackThree_m_pkp"};
   TTreeReaderValue<TLorentzVector> five_p4_pkp = {fReader, "five_p4_pkp"};
   TTreeReaderValue<TLorentzVector> dimuonDiTrkOne_p4_pkp = {fReader, "dimuonDiTrkOne_p4_pkp"};
   TTreeReaderValue<TLorentzVector> dimuonDiTrkTwo_p4_pkp = {fReader, "dimuonDiTrkTwo_p4_pkp"};
   TTreeReaderValue<TLorentzVector> dimuonDiTrkThree_p4_pkp = {fReader, "dimuonDiTrkThree_p4_pkp"};
   TTreeReaderValue<TLorentzVector> ditrakOne_p4_pkp = {fReader, "ditrakOne_p4_pkp"};
   TTreeReaderValue<TLorentzVector> ditrakTwo_p4_pkp = {fReader, "ditrakTwo_p4_pkp"};
   TTreeReaderValue<TLorentzVector> ditrakThree_p4_pkp = {fReader, "ditrakThree_p4_pkp"};
   TTreeReaderValue<TLorentzVector> psiPrimeSame_p4_pkp = {fReader, "psiPrimeSame_p4_pkp"};
   TTreeReaderValue<TLorentzVector> psiPrimeMixed_p4_pkp = {fReader, "psiPrimeMixed_p4_pkp"};
   TTreeReaderValue<TLorentzVector> psiPrimeSame_ditrak_p4_pkp = {fReader, "psiPrimeSame_ditrak_p4_pkp"};
   TTreeReaderValue<TLorentzVector> psiPrimeMixed_ditrak_p4_pkp = {fReader, "psiPrimeMixed_ditrak_p4_pkp"};
   TTreeReaderValue<Double_t> triTrak_m_pkp = {fReader, "triTrak_m_pkp"};
   TTreeReaderValue<Double_t> mass_ppp = {fReader, "mass_ppp"};
   TTreeReaderValue<Double_t> mass_ref_ppp = {fReader, "mass_ref_ppp"};
   TTreeReaderValue<Double_t> vProb_ppp = {fReader, "vProb_ppp"};
   TTreeReaderValue<Double_t> nDof_ppp = {fReader, "nDof_ppp"};
   TTreeReaderValue<Double_t> vChi2_ppp = {fReader, "vChi2_ppp"};
   TTreeReaderValue<Double_t> ctau_ppp = {fReader, "ctau_ppp"};
   TTreeReaderValue<Double_t> ctauErr_ppp = {fReader, "ctauErr_ppp"};
   TTreeReaderValue<Double_t> cosAlpha_ppp = {fReader, "cosAlpha_ppp"};
   TTreeReaderValue<Double_t> onePsiPrime_m_ppp = {fReader, "onePsiPrime_m_ppp"};
   TTreeReaderValue<Double_t> twoPsiPrime_m_ppp = {fReader, "twoPsiPrime_m_ppp"};
   TTreeReaderValue<Double_t> onePsiPrime_p_mppp = {fReader, "onePsiPrime_p_mppp"};
   TTreeReaderValue<Double_t> onePsiPrime_m_mppp = {fReader, "onePsiPrime_m_mppp"};
   TTreeReaderValue<Double_t> twoPsiPrime_p_mppp = {fReader, "twoPsiPrime_p_mppp"};
   TTreeReaderValue<Double_t> twoPsiPrime_m_mppp = {fReader, "twoPsiPrime_m_mppp"};
   TTreeReaderValue<Double_t> dimuonDiTrkOne_m_ppp = {fReader, "dimuonDiTrkOne_m_ppp"};
   TTreeReaderValue<Double_t> dimuonDiTrkTwo_m_ppp = {fReader, "dimuonDiTrkTwo_m_ppp"};
   TTreeReaderValue<Double_t> dimuonDiTrkThree_m_ppp = {fReader, "dimuonDiTrkThree_m_ppp"};
   TTreeReaderValue<Double_t> ditrakOne_m_ppp = {fReader, "ditrakOne_m_ppp"};
   TTreeReaderValue<Double_t> ditrakTwo_m_ppp = {fReader, "ditrakTwo_m_ppp"};
   TTreeReaderValue<Double_t> ditrakThree_m_ppp = {fReader, "ditrakThree_m_ppp"};
   TTreeReaderValue<Double_t> trackOne_m_ppp = {fReader, "trackOne_m_ppp"};
   TTreeReaderValue<Double_t> trackTwo_m_ppp = {fReader, "trackTwo_m_ppp"};
   TTreeReaderValue<Double_t> trackThree_m_ppp = {fReader, "trackThree_m_ppp"};
   TTreeReaderValue<TLorentzVector> five_p4_ppp = {fReader, "five_p4_ppp"};
   TTreeReaderValue<TLorentzVector> dimuonDiTrkOne_p4_ppp = {fReader, "dimuonDiTrkOne_p4_ppp"};
   TTreeReaderValue<TLorentzVector> dimuonDiTrkTwo_p4_ppp = {fReader, "dimuonDiTrkTwo_p4_ppp"};
   TTreeReaderValue<TLorentzVector> dimuonDiTrkThree_p4_ppp = {fReader, "dimuonDiTrkThree_p4_ppp"};
   TTreeReaderValue<TLorentzVector> ditrakOne_p4_ppp = {fReader, "ditrakOne_p4_ppp"};
   TTreeReaderValue<TLorentzVector> ditrakTwo_p4_ppp = {fReader, "ditrakTwo_p4_ppp"};
   TTreeReaderValue<TLorentzVector> ditrakThree_p4_ppp = {fReader, "ditrakThree_p4_ppp"};
   TTreeReaderValue<TLorentzVector> psiPrimeSame_p4_ppp = {fReader, "psiPrimeSame_p4_ppp"};
   TTreeReaderValue<TLorentzVector> psiPrimeMixed_p4_ppp = {fReader, "psiPrimeMixed_p4_ppp"};
   TTreeReaderValue<TLorentzVector> psiPrimeSame_ditrak_p4_ppp = {fReader, "psiPrimeSame_ditrak_p4_ppp"};
   TTreeReaderValue<TLorentzVector> psiPrimeMixed_ditrak_p4_ppp = {fReader, "psiPrimeMixed_ditrak_p4_ppp"};
   TTreeReaderValue<Double_t> triTrak_m_ppp = {fReader, "triTrak_m_ppp"};


  // In data the gen variables are empty and if
  // uncommented rise conflicts in PROOF
  // (it seems it douesn't like empty leaves)

  /*
  TTreeReaderValue<TLorentzVector> gen_dimuonditrk_p4 = {fReader, "gen_dimuonditrk_p4"};
  TTreeReaderValue<TLorentzVector> gen_jpsi_p4 = {fReader, "gen_jpsi_p4"};
  TTreeReaderValue<TLorentzVector> gen_phi_p4 = {fReader, "gen_phi_p4"};
  TTreeReaderValue<TLorentzVector> gen_highKaon_p4 = {fReader, "gen_highKaon_p4"};
  TTreeReaderValue<TLorentzVector> gen_lowMuon_p4 = {fReader, "gen_lowMuon_p4"};
  TTreeReaderValue<TLorentzVector> gen_highMuon_p4 = {fReader, "gen_highMuon_p4"};
  TTreeReaderValue<TLorentzVector> gen_lowKaon_p4 = {fReader, "gen_lowKaon_p4"};

  TTreeReaderValue<Double_t> gen_dimuonditrk_pdg = {fReader, "gen_dimuonditrk_pdg"};
  TTreeReaderValue<Double_t> gen_phi_pdg = {fReader, "gen_phi_pdg"};
  TTreeReaderValue<Double_t> gen_jpsi_pdg = {fReader, "gen_jpsi_pdg"};
  TTreeReaderValue<Double_t> gen_lowMuon_pdg = {fReader, "gen_lowMuon_pdg"};
  TTreeReaderValue<Double_t> gen_highMuon_pdg = {fReader, "gen_highMuon_pdg"};
  TTreeReaderValue<Double_t> gen_highKaon_pdg = {fReader, "gen_highKaon_pdg"};
  TTreeReaderValue<Double_t> gen_lowKaon_pdg = {fReader, "gen_lowKaon_pdg"};
  TTreeReaderValue<Double_t> gen_lowMuon_mompdg = {fReader, "gen_lowMuon_mompdg"};
  TTreeReaderValue<Double_t> gen_highMuon_mompdg = {fReader, "gen_highMuon_mompdg"};
  TTreeReaderValue<Double_t> gen_highKaon_mompdg = {fReader, "gen_highKaon_mompdg"};
  TTreeReaderValue<Double_t> gen_lowKaon_mompdg = {fReader, "gen_lowKaon_mompdg"};
  TTreeReaderValue<Double_t> gen_lowMuon_status = {fReader, "gen_lowMuon_status"};
  TTreeReaderValue<Double_t> gen_highMuon_status = {fReader, "gen_highMuon_status"};
  TTreeReaderValue<Double_t> gen_highKaon_status = {fReader, "gen_highKaon_status"};
  TTreeReaderValue<Double_t> gen_lowKaon_status = {fReader, "gen_lowKaon_status"};
  TTreeReaderValue<Double_t> gen_lowMuon_p = {fReader, "gen_lowMuon_p"};
  TTreeReaderValue<Double_t> gen_highMuon_p = {fReader, "gen_highMuon_p"};
  TTreeReaderValue<Double_t> gen_highKaon_p = {fReader, "gen_highKaon_p"};
  TTreeReaderValue<Double_t> gen_lowKaon_p = {fReader, "gen_lowKaon_p"};
  TTreeReaderValue<Double_t> gen_lowMuon_pt = {fReader, "gen_lowMuon_pt"};
  TTreeReaderValue<Double_t> gen_highMuon_pt = {fReader, "gen_highMuon_pt"};
  TTreeReaderValue<Double_t> gen_highKaon_pt = {fReader, "gen_highKaon_pt"};
  TTreeReaderValue<Double_t> gen_lowKaon_pt = {fReader, "gen_lowKaon_pt"};
  TTreeReaderValue<Double_t> gen_lowMuon_eta = {fReader, "gen_lowMuon_eta"};
  TTreeReaderValue<Double_t> gen_highMuon_eta = {fReader, "gen_highMuon_eta"};
  TTreeReaderValue<Double_t> gen_highKaon_eta = {fReader, "gen_highKaon_eta"};
  TTreeReaderValue<Double_t> gen_lowKaon_eta = {fReader, "gen_lowKaon_eta"};
  TTreeReaderValue<Double_t> gen_lowMuon_phi = {fReader, "gen_lowMuon_phi"};
  TTreeReaderValue<Double_t> gen_highMuon_phi = {fReader, "gen_highMuon_phi"};
  TTreeReaderValue<Double_t> gen_highKaon_phi = {fReader, "gen_highKaon_phi"};
  TTreeReaderValue<Double_t> gen_lowKaon_phi = {fReader, "gen_lowKaon_phi"};
  TTreeReaderValue<Double_t> gen_dimuonditrk_prompt = {fReader, "gen_dimuonditrk_prompt"};
  TTreeReaderValue<Double_t> gen_phi_prompt = {fReader, "gen_phi_prompt"};
  TTreeReaderValue<Double_t> gen_jpsi_prompt = {fReader, "gen_jpsi_prompt"};
  TTreeReaderValue<Double_t> gen_dimuonditrk_pt = {fReader, "gen_dimuonditrk_pt"};
  TTreeReaderValue<Double_t> phigen_phi_pt_pt = {fReader, "gen_phi_pt"};
  TTreeReaderValue<Double_t> gen_jpsi_pt = {fReader, "gen_jpsi_pt"};
  TTreeReaderValue<Double_t> gen_dimuonditrk_p = {fReader, "gen_dimuonditrk_p"};
  TTreeReaderValue<Double_t> phigen_phi_p_p = {fReader, "gen_phi_p"};
  TTreeReaderValue<Double_t> gen_jpsi_p = {fReader, "gen_jpsi_p"};
  TTreeReaderValue<Double_t> gen_dimuonditrk_eta = {fReader, "gen_dimuonditrk_eta"};
  TTreeReaderValue<Double_t> gen_phi_eta = {fReader, "gen_phi_eta"};
  TTreeReaderValue<Double_t> gen_jpsi_eta = {fReader, "gen_jpsi_eta"};
  TTreeReaderValue<Double_t> gen_dimuonditrk_phi = {fReader, "gen_dimuonditrk_phi"};
  TTreeReaderValue<Double_t> gen_phi_phi = {fReader, "gen_phi_phi"};
  TTreeReaderValue<Double_t> gen_jpsi_phi = {fReader, "gen_jpsi_phi"};
  */


  //Output variables
  Float_t out_run, out_event, out_lumi, out_numPrimaryVertices, out_trigger;
Float_t out_noFiveCandidates, out_dimuonditrk_id, out_dimuonditrk_, out_ditrak_, out_dimuon_;
Float_t out_dimuonditrk_m, out_dimuonditrk_pt, out_dimuonditrk_eta, out_dimuonditrk_phi, out_dimuonditrk_p;
Float_t out_dimuon_m, out_dimuon_pt, out_dimuon_eta, out_dimuon_phi, out_dimuon_p;
Float_t out_highKaonMatch, out_lowKaonMatch, out_lowMuonMatch, out_highMuonMatch, out_ditrak_m;
Float_t out_ditrakOne_pt, out_ditrakOne_eta, out_ditrakOne_phi, out_ditrakOne_p, out_ditrakTwo_pt;
Float_t out_ditrakTwo_eta, out_ditrakTwo_phi, out_ditrakTwo_p, out_ditrakThree_pt, out_ditrakThree_eta;
Float_t out_ditrakThree_phi, out_ditrakThree_p, out_highTrack_pt, out_highTrack_eta, out_highTrack_phi;
Float_t out_highTrack_charge, out_lowTrack_pt, out_lowTrack_eta, out_lowTrack_phi, out_lowTrack_charge;
Float_t out_thirdTrack_pt, out_thirdTrack_eta, out_thirdTrack_phi, out_thirdTrack_charge, out_thirdTrack_dz;
Float_t out_thirdTrack_dxy, out_dimuonDiTrkOne_pt, out_dimuonDiTrkOne_eta, out_dimuonDiTrkOne_phi, out_dimuonDiTrkOne_charge;
Float_t out_dimuonDiTrkTwo_pt, out_dimuonDiTrkTwo_eta, out_dimuonDiTrkTwo_phi, out_dimuonDiTrkTwo_charge, out_dimuonDiTrkThree_pt;
Float_t out_dimuonDiTrkThree_eta, out_dimuonDiTrkThree_phi, out_dimuonDiTrkThree_charge, out_psiPrimeSame_pt, out_psiPrimeSame_eta;
Float_t out_psiPrimeSame_phi, out_psiPrimeSame_n, out_psiPrimeSame_p_pt, out_psiPrimeSame_p_eta, out_psiPrimeSame_p_phi;
Float_t out_psiPrimeSame_p_n, out_psiPrimeSame_m_pt, out_psiPrimeSame_m_eta, out_psiPrimeSame_m_phi, out_psiPrimeSame_m_n;
Float_t out_psiPrimeMixed_pt, out_psiPrimeMixed_eta, out_psiPrimeMixed_phi, out_psiPrimeMixed_n, out_psiPrimeMixed_p_pt;
Float_t out_psiPrimeMixed_p_eta, out_psiPrimeMixed_p_phi, out_psiPrimeMixed_p_n, out_psiPrimeMixed_m_pt, out_psiPrimeMixed_m_eta;
Float_t out_psiPrimeMixed_m_phi, out_psiPrimeMixed_m_n, out_psiPrimeSame_ditrak_pt, out_psiPrimeSame_ditrak_eta, out_psiPrimeSame_ditrak_phi;
Float_t out_psiPrimeSame_ditrak_n, out_psiPrimeMixed_ditrak_pt, out_psiPrimeMixed_ditrak_eta, out_psiPrimeMixed_ditrak_phi, out_psiPrimeMixed_ditrak_n;
Float_t out_triTrak_pt, out_triTrak_eta, out_triTrak_phi, out_triTrak_charge, out_mass_kkk;
Float_t out_mass_ref_kkk, out_vProb_kkk, out_nDof_kkk, out_vChi2_kkk, out_ctau_kkk;
Float_t out_ctauErr_kkk, out_cosAlpha_kkk, out_onePsiPrime_m_kkk, out_twoPsiPrime_m_kkk, out_onePsiPrime_p_mkkk;
Float_t out_onePsiPrime_m_mkkk, out_twoPsiPrime_p_mkkk, out_twoPsiPrime_m_mkkk, out_dimuonDiTrkOne_m_kkk, out_dimuonDiTrkTwo_m_kkk;
Float_t out_dimuonDiTrkThree_m_kkk, out_ditrakOne_m_kkk, out_ditrakTwo_m_kkk, out_ditrakThree_m_kkk, out_trackOne_m_kkk;
Float_t out_trackTwo_m_kkk, out_trackThree_m_kkk, out_five__kkk, out_dimuonDiTrkOne__kkk, out_dimuonDiTrkTwo__kkk;
Float_t out_dimuonDiTrkThree__kkk, out_ditrakOne__kkk, out_ditrakTwo__kkk, out_ditrakThree__kkk, out_psiPrimeSame__kkk;
Float_t out_psiPrimeMixed__kkk, out_psiPrimeSame_ditrak__kkk, out_psiPrimeMixed_ditrak__kkk, out_triTrak_m_kkk, out_mass_ppk;
Float_t out_mass_ref_ppk, out_vProb_ppk, out_nDof_ppk, out_vChi2_ppk, out_ctau_ppk;
Float_t out_ctauErr_ppk, out_cosAlpha_ppk, out_onePsiPrime_m_ppk, out_twoPsiPrime_m_ppk, out_onePsiPrime_p_mppk;
Float_t out_onePsiPrime_m_mppk, out_twoPsiPrime_p_mppk, out_twoPsiPrime_m_mppk, out_dimuonDiTrkOne_m_ppk, out_dimuonDiTrkTwo_m_ppk;
Float_t out_dimuonDiTrkThree_m_ppk, out_ditrakOne_m_ppk, out_ditrakTwo_m_ppk, out_ditrakThree_m_ppk, out_trackOne_m_ppk;
Float_t out_trackTwo_m_ppk, out_trackThree_m_ppk, out_five__ppk, out_dimuonDiTrkOne__ppk, out_dimuonDiTrkTwo__ppk;
Float_t out_dimuonDiTrkThree__ppk, out_ditrakOne__ppk, out_ditrakTwo__ppk, out_ditrakThree__ppk, out_psiPrimeSame__ppk;
Float_t out_psiPrimeMixed__ppk, out_psiPrimeSame_ditrak__ppk, out_psiPrimeMixed_ditrak__ppk, out_triTrak_m_ppk, out_mass_kpp;
Float_t out_mass_ref_kpp, out_vProb_kpp, out_nDof_kpp, out_vChi2_kpp, out_ctau_kpp;
Float_t out_ctauErr_kpp, out_cosAlpha_kpp, out_onePsiPrime_m_kpp, out_twoPsiPrime_m_kpp, out_onePsiPrime_p_mkpp;
Float_t out_onePsiPrime_m_mkpp, out_twoPsiPrime_p_mkpp, out_twoPsiPrime_m_mkpp, out_dimuonDiTrkOne_m_kpp, out_dimuonDiTrkTwo_m_kpp;
Float_t out_dimuonDiTrkThree_m_kpp, out_ditrakOne_m_kpp, out_ditrakTwo_m_kpp, out_ditrakThree_m_kpp, out_trackOne_m_kpp;
Float_t out_trackTwo_m_kpp, out_trackThree_m_kpp, out_five__kpp, out_dimuonDiTrkOne__kpp, out_dimuonDiTrkTwo__kpp;
Float_t out_dimuonDiTrkThree__kpp, out_ditrakOne__kpp, out_ditrakTwo__kpp, out_ditrakThree__kpp, out_psiPrimeSame__kpp;
Float_t out_psiPrimeMixed__kpp, out_psiPrimeSame_ditrak__kpp, out_psiPrimeMixed_ditrak__kpp, out_triTrak_m_kpp, out_mass_pkp;
Float_t out_mass_ref_pkp, out_vProb_pkp, out_nDof_pkp, out_vChi2_pkp, out_ctau_pkp;
Float_t out_ctauErr_pkp, out_cosAlpha_pkp, out_onePsiPrime_m_pkp, out_twoPsiPrime_m_pkp, out_onePsiPrime_p_mpkp;
Float_t out_onePsiPrime_m_mpkp, out_twoPsiPrime_p_mpkp, out_twoPsiPrime_m_mpkp, out_dimuonDiTrkOne_m_pkp, out_dimuonDiTrkTwo_m_pkp;
Float_t out_dimuonDiTrkThree_m_pkp, out_ditrakOne_m_pkp, out_ditrakTwo_m_pkp, out_ditrakThree_m_pkp, out_trackOne_m_pkp;
Float_t out_trackTwo_m_pkp, out_trackThree_m_pkp, out_five__pkp, out_dimuonDiTrkOne__pkp, out_dimuonDiTrkTwo__pkp;
Float_t out_dimuonDiTrkThree__pkp, out_ditrakOne__pkp, out_ditrakTwo__pkp, out_ditrakThree__pkp, out_psiPrimeSame__pkp;
Float_t out_psiPrimeMixed__pkp, out_psiPrimeSame_ditrak__pkp, out_psiPrimeMixed_ditrak__pkp, out_triTrak_m_pkp, out_mass_ppp;
Float_t out_mass_ref_ppp, out_vProb_ppp, out_nDof_ppp, out_vChi2_ppp, out_ctau_ppp;
Float_t out_ctauErr_ppp, out_cosAlpha_ppp, out_onePsiPrime_m_ppp, out_twoPsiPrime_m_ppp, out_onePsiPrime_p_mppp;
Float_t out_onePsiPrime_m_mppp, out_twoPsiPrime_p_mppp, out_twoPsiPrime_m_mppp, out_dimuonDiTrkOne_m_ppp, out_dimuonDiTrkTwo_m_ppp;
Float_t out_dimuonDiTrkThree_m_ppp, out_ditrakOne_m_ppp, out_ditrakTwo_m_ppp, out_ditrakThree_m_ppp, out_trackOne_m_ppp;
Float_t out_trackTwo_m_ppp, out_trackThree_m_ppp, out_five__ppp, out_dimuonDiTrkOne__ppp, out_dimuonDiTrkTwo__ppp;
Float_t out_dimuonDiTrkThree__ppp, out_ditrakOne__ppp, out_ditrakTwo__ppp, out_ditrakThree__ppp, out_psiPrimeSame__ppp;
Float_t out_psiPrimeMixed__ppp, out_psiPrimeSame_ditrak__ppp, out_psiPrimeMixed_ditrak__ppp, out_triTrak_m_ppp;

  FiveTracks(TTree * /*tree*/ =0) { }
  virtual ~FiveTracks() { }
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

  ClassDef(FiveTracks,0);

};

#endif

#ifdef FiveTracks_cxx
void FiveTracks::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the reader is initialized.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  fReader.SetTree(tree);
}

Bool_t FiveTracks::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}


#endif // #ifdef FiveTracks_cxx
