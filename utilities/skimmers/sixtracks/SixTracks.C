#define SixTracks_cxx
// The class definition in SixTracks.h has been generated automatically
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
// root> T->Process("SixTracks.C")
// root> T->Process("SixTracks.C","some options")
// root> T->Process("SixTracks.C+")
//


#include "SixTracks.h"
#include <TH2.h>
#include <TStyle.h>

void SixTracks::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void SixTracks::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t SixTracks::Process(Long64_t entry)
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

   out_run = 	(Float_t)(*run);
   out_event = 	(Float_t)(*event);
   out_lumi = 	(Float_t)(*lumi);
   out_numPrimaryVertices = 	(Float_t)(*numPrimaryVertices);
   out_trigger = 	(Float_t)(*trigger);
   out_noSixCandidates = 	(Float_t)(*noSixCandidates);
   out_five_id = 	(Float_t)(*five_id);
   out_dimuon_id = 	(Float_t)(*dimuon_id);
   out_p_id = 	(Float_t)(*p_id);
   out_m_id = 	(Float_t)(*m_id);
   out_t_id = 	(Float_t)(*t_id);
   out_f_id = 	(Float_t)(*f_id);
   out_six_p4 = 	(Float_t)(*six_p4);
   out_five_p4 = 	(Float_t)(*five_p4);
   out_dimuonditrk_p4 = 	(Float_t)(*dimuonditrk_p4);
   out_ditrack_p4 = 	(Float_t)(*ditrack_p4);
   out_dimuon_p4 = 	(Float_t)(*dimuon_p4);
   out_lowMuon_p4 = 	(Float_t)(*lowMuon_p4);
   out_highMuon_p4 = 	(Float_t)(*highMuon_p4);
   out_highKaon_p4 = 	(Float_t)(*highKaon_p4);
   out_lowKaon_p4 = 	(Float_t)(*lowKaon_p4);
   out_thirdKaon_p4 = 	(Float_t)(*thirdKaon_p4);
   out_fourthKaon_p4 = 	(Float_t)(*fourthKaon_p4);
   out_highPion_p4 = 	(Float_t)(*highPion_p4);
   out_lowPion_p4 = 	(Float_t)(*lowPion_p4);
   out_thirdPion_p4 = 	(Float_t)(*thirdPion_p4);
   out_fourthPion_p4 = 	(Float_t)(*fourthPion_p4);
   out_highProton_p4 = 	(Float_t)(*highProton_p4);
   out_lowProton_p4 = 	(Float_t)(*lowProton_p4);
   out_thirdProton_p4 = 	(Float_t)(*thirdProton_p4);
   out_fourthProton_p4 = 	(Float_t)(*fourthProton_p4);
   out_dimuonditrk_m = 	(Float_t)(*dimuonditrk_m);
   out_dimuonditrk_pt = 	(Float_t)(*dimuonditrk_pt);
   out_dimuonditrk_eta = 	(Float_t)(*dimuonditrk_eta);
   out_dimuonditrk_phi = 	(Float_t)(*dimuonditrk_phi);
   out_dimuonditrk_p = 	(Float_t)(*dimuonditrk_p);
   out_dimuon_m = 	(Float_t)(*dimuon_m);
   out_dimuon_pt = 	(Float_t)(*dimuon_pt);
   out_dimuon_eta = 	(Float_t)(*dimuon_eta);
   out_dimuon_phi = 	(Float_t)(*dimuon_phi);
   out_dimuon_p = 	(Float_t)(*dimuon_p);
   out_highTrackMatch = 	(Float_t)(*highTrackMatch);
   out_lowTrackMatch = 	(Float_t)(*lowTrackMatch);
   out_lowMuonMatch = 	(Float_t)(*lowMuonMatch);
   out_highMuonMatch = 	(Float_t)(*highMuonMatch);
   out_thirdTrackMatch = 	(Float_t)(*thirdTrackMatch);
   out_fourthTrackMatch = 	(Float_t)(*fourthTrackMatch);
   out_ditrack_m = 	(Float_t)(*ditrack_m);
   out_diTrackOne_pt = 	(Float_t)(*diTrackOne_pt);
   out_diTrackOne_eta = 	(Float_t)(*diTrackOne_eta);
   out_diTrackOne_phi = 	(Float_t)(*diTrackOne_phi);
   out_diTrackOne_p = 	(Float_t)(*diTrackOne_p);
   out_diTrackTwo_pt = 	(Float_t)(*diTrackTwo_pt);
   out_diTrackTwo_eta = 	(Float_t)(*diTrackTwo_eta);
   out_diTrackTwo_phi = 	(Float_t)(*diTrackTwo_phi);
   out_diTrackTwo_p = 	(Float_t)(*diTrackTwo_p);
   out_diTrackThree_pt = 	(Float_t)(*diTrackThree_pt);
   out_diTrackThree_eta = 	(Float_t)(*diTrackThree_eta);
   out_diTrackThree_phi = 	(Float_t)(*diTrackThree_phi);
   out_diTrackThree_p = 	(Float_t)(*diTrackThree_p);
   out_diTrackFour_pt = 	(Float_t)(*diTrackFour_pt);
   out_diTrackFour_eta = 	(Float_t)(*diTrackFour_eta);
   out_diTrackFour_phi = 	(Float_t)(*diTrackFour_phi);
   out_diTrackFour_p = 	(Float_t)(*diTrackFour_p);
   out_diTrackFive_pt = 	(Float_t)(*diTrackFive_pt);
   out_diTrackFive_eta = 	(Float_t)(*diTrackFive_eta);
   out_diTrackFive_phi = 	(Float_t)(*diTrackFive_phi);
   out_diTrackFive_p = 	(Float_t)(*diTrackFive_p);
   out_diTrackSix_pt = 	(Float_t)(*diTrackSix_pt);
   out_diTrackSix_eta = 	(Float_t)(*diTrackSix_eta);
   out_diTrackSix_phi = 	(Float_t)(*diTrackSix_phi);
   out_diTrackSix_p = 	(Float_t)(*diTrackSix_p);
   out_dimuonDiTrkOne_mmpp = 	(Float_t)(*dimuonDiTrkOne_mmpp);
   out_dimuonDiTrkTwo_mmpp = 	(Float_t)(*dimuonDiTrkTwo_mmpp);
   out_dimuonDiTrkThree_mmpp = 	(Float_t)(*dimuonDiTrkThree_mmpp);
   out_dimuonDiTrkFour_mmpp = 	(Float_t)(*dimuonDiTrkFour_mmpp);
   out_dimuonDiTrkOne_mm = 	(Float_t)(*dimuonDiTrkOne_mm);
   out_dimuonDiTrkTwo_mmkk = 	(Float_t)(*dimuonDiTrkTwo_mmkk);
   out_dimuonDiTrkThree_mmkk = 	(Float_t)(*dimuonDiTrkThree_mmkk);
   out_dimuonDiTrkFour_mmkk = 	(Float_t)(*dimuonDiTrkFour_mmkk);
   out_highMuon_pt = 	(Float_t)(*highMuon_pt);
   out_highMuon_eta = 	(Float_t)(*highMuon_eta);
   out_highMuon_phi = 	(Float_t)(*highMuon_phi);
   out_highMuon_charge = 	(Float_t)(*highMuon_charge);
   out_highMuon_dz = 	(Float_t)(*highMuon_dz);
   out_highMuon_dxy = 	(Float_t)(*highMuon_dxy);
   out_lowMuon_pt = 	(Float_t)(*lowMuon_pt);
   out_lowMuon_eta = 	(Float_t)(*lowMuon_eta);
   out_lowMuon_phi = 	(Float_t)(*lowMuon_phi);
   out_lowMuon_charge = 	(Float_t)(*lowMuon_charge);
   out_lowMuon_dz = 	(Float_t)(*lowMuon_dz);
   out_lowMuon_dxy = 	(Float_t)(*lowMuon_dxy);
   out_highTrack_pt = 	(Float_t)(*highTrack_pt);
   out_highTrack_eta = 	(Float_t)(*highTrack_eta);
   out_highTrack_phi = 	(Float_t)(*highTrack_phi);
   out_highTrack_charge = 	(Float_t)(*highTrack_charge);
   out_highTrack_dz = 	(Float_t)(*highTrack_dz);
   out_highTrack_dxy = 	(Float_t)(*highTrack_dxy);
   out_lowTrack_pt = 	(Float_t)(*lowTrack_pt);
   out_lowTrack_eta = 	(Float_t)(*lowTrack_eta);
   out_lowTrack_phi = 	(Float_t)(*lowTrack_phi);
   out_lowTrack_charge = 	(Float_t)(*lowTrack_charge);
   out_lowTrack_dz = 	(Float_t)(*lowTrack_dz);
   out_lowTrack_dxy = 	(Float_t)(*lowTrack_dxy);
   out_thirdTrack_pt = 	(Float_t)(*thirdTrack_pt);
   out_thirdTrack_eta = 	(Float_t)(*thirdTrack_eta);
   out_thirdTrack_phi = 	(Float_t)(*thirdTrack_phi);
   out_thirdTrack_charge = 	(Float_t)(*thirdTrack_charge);
   out_thirdTrack_dz = 	(Float_t)(*thirdTrack_dz);
   out_thirdTrack_dxy = 	(Float_t)(*thirdTrack_dxy);
   out_dimuonDiTrkOne_pt = 	(Float_t)(*dimuonDiTrkOne_pt);
   out_dimuonDiTrkOne_eta = 	(Float_t)(*dimuonDiTrkOne_eta);
   out_dimuonDiTrkOne_phi = 	(Float_t)(*dimuonDiTrkOne_phi);
   out_dimuonDiTrkOne_charge = 	(Float_t)(*dimuonDiTrkOne_charge);
   out_dimuonDiTrkOne_p = 	(Float_t)(*dimuonDiTrkOne_p);
   out_dimuonDiTrkTwo_pt = 	(Float_t)(*dimuonDiTrkTwo_pt);
   out_dimuonDiTrkTwo_eta = 	(Float_t)(*dimuonDiTrkTwo_eta);
   out_dimuonDiTrkTwo_phi = 	(Float_t)(*dimuonDiTrkTwo_phi);
   out_dimuonDiTrkTwo_charge = 	(Float_t)(*dimuonDiTrkTwo_charge);
   out_dimuonDiTrkTwo_p = 	(Float_t)(*dimuonDiTrkTwo_p);
   out_dimuonDiTrkThree_pt = 	(Float_t)(*dimuonDiTrkThree_pt);
   out_dimuonDiTrkThree_eta = 	(Float_t)(*dimuonDiTrkThree_eta);
   out_dimuonDiTrkThree_phi = 	(Float_t)(*dimuonDiTrkThree_phi);
   out_dimuonDiTrkThree_charge = 	(Float_t)(*dimuonDiTrkThree_charge);
   out_dimuonDiTrkThree_p = 	(Float_t)(*dimuonDiTrkThree_p);
   out_dimuonDiTrkFour_pt = 	(Float_t)(*dimuonDiTrkFour_pt);
   out_dimuonDiTrkFour_eta = 	(Float_t)(*dimuonDiTrkFour_eta);
   out_dimuonDiTrkFour_phi = 	(Float_t)(*dimuonDiTrkFour_phi);
   out_dimuonDiTrkFour_charge = 	(Float_t)(*dimuonDiTrkFour_charge);
   out_dimuonDiTrkFour_p = 	(Float_t)(*dimuonDiTrkFour_p);
   out_dimuonDiTrkFive_pt = 	(Float_t)(*dimuonDiTrkFive_pt);
   out_dimuonDiTrkFive_eta = 	(Float_t)(*dimuonDiTrkFive_eta);
   out_dimuonDiTrkFive_phi = 	(Float_t)(*dimuonDiTrkFive_phi);
   out_dimuonDiTrkFive_charge = 	(Float_t)(*dimuonDiTrkFive_charge);
   out_dimuonDiTrkFive_p = 	(Float_t)(*dimuonDiTrkFive_p);
   out_dimuonDiTrkSix_pt = 	(Float_t)(*dimuonDiTrkSix_pt);
   out_dimuonDiTrkSix_eta = 	(Float_t)(*dimuonDiTrkSix_eta);
   out_dimuonDiTrkSix_phi = 	(Float_t)(*dimuonDiTrkSix_phi);
   out_dimuonDiTrkSix_charge = 	(Float_t)(*dimuonDiTrkSix_charge);
   out_dimuonDiTrkSix_p = 	(Float_t)(*dimuonDiTrkSix_p);
   out_dimuon_vProb = 	(Float_t)(*dimuon_vProb);
   out_dimuon_vChi2 = 	(Float_t)(*dimuon_vChi2);
   out_dimuon_DCA = 	(Float_t)(*dimuon_DCA);
   out_dimuon_ctauPV = 	(Float_t)(*dimuon_ctauPV);
   out_dimuon_ctauErrPV = 	(Float_t)(*dimuon_ctauErrPV);
   out_dimuon_cosAlpha = 	(Float_t)(*dimuon_cosAlpha);
   out_triTrack_m = 	(Float_t)(*triTrack_m);
   out_triTrack_pt = 	(Float_t)(*triTrack_pt);
   out_triTrack_eta = 	(Float_t)(*triTrack_eta);
   out_triTrack_phi = 	(Float_t)(*triTrack_phi);
   out_triTrack_charge = 	(Float_t)(*triTrack_charge);
   out_dimuonditrk_vProb = 	(Float_t)(*dimuonditrk_vProb);
   out_dimuonditrk_vChi2 = 	(Float_t)(*dimuonditrk_vChi2);
   out_dimuonditrk_nDof = 	(Float_t)(*dimuonditrk_nDof);
   out_dimuonditrk_charge = 	(Float_t)(*dimuonditrk_charge);
   out_dimuonditrk_cosAlpha = 	(Float_t)(*dimuonditrk_cosAlpha);
   out_dimuonditrk_ctauPV = 	(Float_t)(*dimuonditrk_ctauPV);
   out_dimuonditrk_ctauErrPV = 	(Float_t)(*dimuonditrk_ctauErrPV);
   out_dimuonditrk_cosAlphaCA = 	(Float_t)(*dimuonditrk_cosAlphaCA);
   out_dimuonditrk_ctauPVCA = 	(Float_t)(*dimuonditrk_ctauPVCA);
   out_dimuonditrk_ctauErrPVCA = 	(Float_t)(*dimuonditrk_ctauErrPVCA);
   out_dimuonditrk_cosAlphaDZ = 	(Float_t)(*dimuonditrk_cosAlphaDZ);
   out_dimuonditrk_ctauPVDZ = 	(Float_t)(*dimuonditrk_ctauPVDZ);
   out_dimuonditrk_ctauErrPVDZ = 	(Float_t)(*dimuonditrk_ctauErrPVDZ);
   out_dimuonditrk_cosAlphaBS = 	(Float_t)(*dimuonditrk_cosAlphaBS);
   out_dimuonditrk_ctauPVBS = 	(Float_t)(*dimuonditrk_ctauPVBS);
   out_dimuonditrk_ctauErrPVBS = 	(Float_t)(*dimuonditrk_ctauErrPVBS);
   out_dimuonditrk_vx = 	(Float_t)(*dimuonditrk_vx);
   out_dimuonditrk_vy = 	(Float_t)(*dimuonditrk_vy);
   out_dimuonditrk_vz = 	(Float_t)(*dimuonditrk_vz);
   out_dca_m1m2 = 	(Float_t)(*dca_m1m2);
   out_dca_m1t1 = 	(Float_t)(*dca_m1t1);
   out_dca_m1t2 = 	(Float_t)(*dca_m1t2);
   out_dca_m2t1 = 	(Float_t)(*dca_m2t1);
   out_dca_m2t2 = 	(Float_t)(*dca_m2t2);
   out_dca_t1t2 = 	(Float_t)(*dca_t1t2);
   out_dca_m1t3 = 	(Float_t)(*dca_m1t3);
   out_dca_m2t3 = 	(Float_t)(*dca_m2t3);
   out_dca_t1t3 = 	(Float_t)(*dca_t1t3);
   out_dca_t2t3 = 	(Float_t)(*dca_t2t3);
   out_dca_m1t4 = 	(Float_t)(*dca_m1t4);
   out_dca_m2t4 = 	(Float_t)(*dca_m2t4);
   out_dca_t1t4 = 	(Float_t)(*dca_t1t4);
   out_dca_t2t4 = 	(Float_t)(*dca_t2t4);
   out_dca_t3t4 = 	(Float_t)(*dca_t3t4);
   out_highTrackMuonDR = 	(Float_t)(*highTrackMuonDR);
   out_highTrackMuonDP = 	(Float_t)(*highTrackMuonDP);
   out_highTrackMuonDPt = 	(Float_t)(*highTrackMuonDPt);
   out_lowTrackMuonDR = 	(Float_t)(*lowTrackMuonDR);
   out_lowTrackMuonDP = 	(Float_t)(*lowTrackMuonDP);
   out_lowTrackMuonDPt = 	(Float_t)(*lowTrackMuonDPt);
   out_thirdTrackMuonDR = 	(Float_t)(*thirdTrackMuonDR);
   out_thirdTrackMuonDP = 	(Float_t)(*thirdTrackMuonDP);
   out_thirdTrackMuonDPt = 	(Float_t)(*thirdTrackMuonDPt);
   out_fourthTrackMuonDR = 	(Float_t)(*fourthTrackMuonDR);
   out_fourthTrackMuonDP = 	(Float_t)(*fourthTrackMuonDP);
   out_fourthTrackMuonDPt = 	(Float_t)(*fourthTrackMuonDPt);
   out_tPFromPV = 	(Float_t)(*tPFromPV);
   out_tMFromPV = 	(Float_t)(*tMFromPV);
   out_tTFromPV = 	(Float_t)(*tTFromPV);
   out_tFFromPV = 	(Float_t)(*tFFromPV);
   out_tPFromPVCA = 	(Float_t)(*tPFromPVCA);
   out_tMFromPVCA = 	(Float_t)(*tMFromPVCA);
   out_tTFromPVCA = 	(Float_t)(*tTFromPVCA);
   out_tFFromPVCA = 	(Float_t)(*tFFromPVCA);
   out_tPFromPVDZ = 	(Float_t)(*tPFromPVDZ);
   out_tMFromPVDZ = 	(Float_t)(*tMFromPVDZ);
   out_tTFromPVDZ = 	(Float_t)(*tTFromPVDZ);
   out_tFFromPVDZ = 	(Float_t)(*tFFromPVDZ);
   out_five_m = 	(Float_t)(*five_m);
   out_five_m_ref = 	(Float_t)(*five_m_ref);
   out_five_mass_ppk = 	(Float_t)(*five_mass_ppk);
   out_five_mass_kpp = 	(Float_t)(*five_mass_kpp);
   out_five_mass_pkp = 	(Float_t)(*five_mass_pkp);
   out_five_mass_ppp = 	(Float_t)(*five_mass_ppp);
   out_fiveOne_pt = 	(Float_t)(*fiveOne_pt);
   out_fiveOne_eta = 	(Float_t)(*fiveOne_eta);
   out_fiveOne_phi = 	(Float_t)(*fiveOne_phi);
   out_fiveOne_p = 	(Float_t)(*fiveOne_p);
   out_fiveTwo_pt = 	(Float_t)(*fiveTwo_pt);
   out_fiveTwo_eta = 	(Float_t)(*fiveTwo_eta);
   out_fiveTwo_phi = 	(Float_t)(*fiveTwo_phi);
   out_fiveTwo_p = 	(Float_t)(*fiveTwo_p);
   out_fiveThree_pt = 	(Float_t)(*fiveThree_pt);
   out_fiveThree_eta = 	(Float_t)(*fiveThree_eta);
   out_fiveThree_phi = 	(Float_t)(*fiveThree_phi);
   out_fiveThree_p = 	(Float_t)(*fiveThree_p);
   out_fiveFour_pt = 	(Float_t)(*fiveFour_pt);
   out_fiveFour_eta = 	(Float_t)(*fiveFour_eta);
   out_fiveFour_phi = 	(Float_t)(*fiveFour_phi);
   out_fiveFour_p = 	(Float_t)(*fiveFour_p);
   out_fiveFive_pt = 	(Float_t)(*fiveFive_pt);
   out_fiveFive_eta = 	(Float_t)(*fiveFive_eta);
   out_fiveFive_phi = 	(Float_t)(*fiveFive_phi);
   out_fiveFive_p = 	(Float_t)(*fiveFive_p);
   out_five_cosAlpha = 	(Float_t)(*five_cosAlpha);
   out_five_ctauPV = 	(Float_t)(*five_ctauPV);
   out_five_ctauErrPV = 	(Float_t)(*five_ctauErrPV);
   out_five_cosAlphaCA = 	(Float_t)(*five_cosAlphaCA);
   out_five_ctauPVCA = 	(Float_t)(*five_ctauPVCA);
   out_five_ctauErrPVCA = 	(Float_t)(*five_ctauErrPVCA);
   out_five_cosAlphaDZ = 	(Float_t)(*five_cosAlphaDZ);
   out_five_ctauPVDZ = 	(Float_t)(*five_ctauPVDZ);
   out_five_ctauErrPVDZ = 	(Float_t)(*five_ctauErrPVDZ);
   out_five_cosAlphaBS = 	(Float_t)(*five_cosAlphaBS);
   out_five_ctauPVBS = 	(Float_t)(*five_ctauPVBS);
   out_five_ctauErrPVBS = 	(Float_t)(*five_ctauErrPVBS);
   out_five_vProb = 	(Float_t)(*five_vProb);
   out_five_nDof = 	(Float_t)(*five_nDof);
   out_five_vChi2 = 	(Float_t)(*five_vChi2);
   out_five_vx = 	(Float_t)(*five_vx);
   out_five_vy = 	(Float_t)(*five_vy);
   out_five_vz = 	(Float_t)(*five_vz);
   out_five_charge = 	(Float_t)(*five_charge);
   out_bestPV_X = 	(Float_t)(*bestPV_X);
   out_bestPV_Y = 	(Float_t)(*bestPV_Y);
   out_bestPV_Z = 	(Float_t)(*bestPV_Z);
   out_cosAlphaPV_X = 	(Float_t)(*cosAlphaPV_X);
   out_cosAlphaPV_Y = 	(Float_t)(*cosAlphaPV_Y);
   out_cosAlphaPV_Z = 	(Float_t)(*cosAlphaPV_Z);
   out_bS_X = 	(Float_t)(*bS_X);
   out_bS_Y = 	(Float_t)(*bS_Y);
   out_bS_Z = 	(Float_t)(*bS_Z);
   out_zPV_X = 	(Float_t)(*zPV_X);
   out_zPV_Y = 	(Float_t)(*zPV_Y);
   out_zPV_Z = 	(Float_t)(*zPV_Z);
   out_lowMuon_isTight = 	(Float_t)(*lowMuon_isTight);
   out_lowMuon_isLoose = 	(Float_t)(*lowMuon_isLoose);
   out_lowMuon_isSoft = 	(Float_t)(*lowMuon_isSoft);
   out_lowMuon_isMedium = 	(Float_t)(*lowMuon_isMedium);
   out_lowMuon_isHighPt = 	(Float_t)(*lowMuon_isHighPt);
   out_lowMuon_isTracker = 	(Float_t)(*lowMuon_isTracker);
   out_lowMuon_isGlobal = 	(Float_t)(*lowMuon_isGlobal);
   out_lowMuon_NPixelHits = 	(Float_t)(*lowMuon_NPixelHits);
   out_lowMuon_NStripHits = 	(Float_t)(*lowMuon_NStripHits);
   out_lowMuon_NTrackhits = 	(Float_t)(*lowMuon_NTrackhits);
   out_lowMuon_NBPixHits = 	(Float_t)(*lowMuon_NBPixHits);
   out_lowMuon_NPixLayers = 	(Float_t)(*lowMuon_NPixLayers);
   out_lowMuon_NTraLayers = 	(Float_t)(*lowMuon_NTraLayers);
   out_lowMuon_NStrLayers = 	(Float_t)(*lowMuon_NStrLayers);
   out_lowMuon_NBPixLayers = 	(Float_t)(*lowMuon_NBPixLayers);
   out_highMuon_isTight = 	(Float_t)(*highMuon_isTight);
   out_highMuon_isLoose = 	(Float_t)(*highMuon_isLoose);
   out_highMuon_isSoft = 	(Float_t)(*highMuon_isSoft);
   out_highMuon_isMedium = 	(Float_t)(*highMuon_isMedium);
   out_highMuon_isHighPt = 	(Float_t)(*highMuon_isHighPt);
   out_highMuon_isTracker = 	(Float_t)(*highMuon_isTracker);
   out_highMuon_isGlobal = 	(Float_t)(*highMuon_isGlobal);
   out_highMuon_NPixelHits = 	(Float_t)(*highMuon_NPixelHits);
   out_highMuon_NStripHits = 	(Float_t)(*highMuon_NStripHits);
   out_highMuon_NTrackhits = 	(Float_t)(*highMuon_NTrackhits);
   out_highMuon_NBPixHits = 	(Float_t)(*highMuon_NBPixHits);
   out_highMuon_NPixLayers = 	(Float_t)(*highMuon_NPixLayers);
   out_highMuon_NTraLayers = 	(Float_t)(*highMuon_NTraLayers);
   out_highMuon_NStrLayers = 	(Float_t)(*highMuon_NStrLayers);
   out_highMuon_NBPixLayers = 	(Float_t)(*highMuon_NBPixLayers);
   out_lowMuon_type = 	(Float_t)(*lowMuon_type);
   out_highMuon_type = 	(Float_t)(*highMuon_type);
   out_highTrack_NPixelHits = 	(Float_t)(*highTrack_NPixelHits);
   out_highTrack_NStripHits = 	(Float_t)(*highTrack_NStripHits);
   out_highTrack_NTrackhits = 	(Float_t)(*highTrack_NTrackhits);
   out_highTrack_NBPixHits = 	(Float_t)(*highTrack_NBPixHits);
   out_highTrack_NPixLayers = 	(Float_t)(*highTrack_NPixLayers);
   out_highTrack_NTraLayers = 	(Float_t)(*highTrack_NTraLayers);
   out_highTrack_NStrLayers = 	(Float_t)(*highTrack_NStrLayers);
   out_highTrack_NBPixLayers = 	(Float_t)(*highTrack_NBPixLayers);
   out_lowTrack_NPixelHits = 	(Float_t)(*lowTrack_NPixelHits);
   out_lowTrack_NStripHits = 	(Float_t)(*lowTrack_NStripHits);
   out_lowTrack_NTrackhits = 	(Float_t)(*lowTrack_NTrackhits);
   out_lowTrack_NBPixHits = 	(Float_t)(*lowTrack_NBPixHits);
   out_lowTrack_NPixLayers = 	(Float_t)(*lowTrack_NPixLayers);
   out_lowTrack_NTraLayers = 	(Float_t)(*lowTrack_NTraLayers);
   out_lowTrack_NStrLayers = 	(Float_t)(*lowTrack_NStrLayers);
   out_lowTrack_NBPixLayers = 	(Float_t)(*lowTrack_NBPixLayers);
   out_thirdTrack_NPixelHits = 	(Float_t)(*thirdTrack_NPixelHits);
   out_thirdTrack_NStripHits = 	(Float_t)(*thirdTrack_NStripHits);
   out_thirdTrack_NTrackhits = 	(Float_t)(*thirdTrack_NTrackhits);
   out_thirdTrack_NBPixHits = 	(Float_t)(*thirdTrack_NBPixHits);
   out_thirdTrack_NPixLayers = 	(Float_t)(*thirdTrack_NPixLayers);
   out_thirdTrack_NTraLayers = 	(Float_t)(*thirdTrack_NTraLayers);
   out_thirdTrack_NStrLayers = 	(Float_t)(*thirdTrack_NStrLayers);
   out_thirdTrack_NBPixLayers = 	(Float_t)(*thirdTrack_NBPixLayers);
   out_fourthTrack_NPixLayers = 	(Float_t)(*fourthTrack_NPixLayers);
   out_fourthTrack_NTraLayers = 	(Float_t)(*fourthTrack_NTraLayers);
   out_fourthTrack_NStrLayers = 	(Float_t)(*fourthTrack_NStrLayers);
   out_fourthTrack_NBPixLayers = 	(Float_t)(*fourthTrack_NBPixLayers);
   out_fourthTrack_NPixelHits = 	(Float_t)(*fourthTrack_NPixelHits);
   out_fourthTrack_NStripHits = 	(Float_t)(*fourthTrack_NStripHits);
   out_fourthTrack_NTrackhits = 	(Float_t)(*fourthTrack_NTrackhits);
   out_fourthTrack_NBPixHits = 	(Float_t)(*fourthTrack_NBPixHits);
   out_six_m = 	(Float_t)(*six_m);
   out_six_m_ref = 	(Float_t)(*six_m_ref);
   out_six_mass_ppkk = 	(Float_t)(*six_mass_ppkk);
   out_six_mass_pkpk = 	(Float_t)(*six_mass_pkpk);
   out_six_mass_pkkk = 	(Float_t)(*six_mass_pkkk);
   out_six_mass_kpkp = 	(Float_t)(*six_mass_kpkp);
   out_six_mass_kppk = 	(Float_t)(*six_mass_kppk);
   out_six_mass_kkkk = 	(Float_t)(*six_mass_kkkk);
   out_six_pt = 	(Float_t)(*six_pt);
   out_six_eta = 	(Float_t)(*six_eta);
   out_six_phi = 	(Float_t)(*six_phi);
   out_six_p = 	(Float_t)(*six_p);
   out_six_cosAlpha = 	(Float_t)(*six_cosAlpha);
   out_six_ctauPV = 	(Float_t)(*six_ctauPV);
   out_six_ctauErrPV = 	(Float_t)(*six_ctauErrPV);
   out_six_cosAlphaCA = 	(Float_t)(*six_cosAlphaCA);
   out_six_ctauPVCA = 	(Float_t)(*six_ctauPVCA);
   out_six_ctauErrPVCA = 	(Float_t)(*six_ctauErrPVCA);
   out_six_cosAlphaDZ = 	(Float_t)(*six_cosAlphaDZ);
   out_six_ctauPVDZ = 	(Float_t)(*six_ctauPVDZ);
   out_six_ctauErrPVDZ = 	(Float_t)(*six_ctauErrPVDZ);
   out_six_cosAlphaBS = 	(Float_t)(*six_cosAlphaBS);
   out_six_ctauPVBS = 	(Float_t)(*six_ctauPVBS);
   out_six_ctauErrPVBS = 	(Float_t)(*six_ctauErrPVBS);
   out_six_vProb = 	(Float_t)(*six_vProb);
   out_six_nDof = 	(Float_t)(*six_nDof);
   out_six_vChi2 = 	(Float_t)(*six_vChi2);
   out_six_vx = 	(Float_t)(*six_vx);
   out_six_vy = 	(Float_t)(*six_vy);
   out_six_vz = 	(Float_t)(*six_vz);
   out_six_charge = 	(Float_t)(*six_charge);


   return kTRUE;
}

void SixTracks::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void SixTracks::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}
