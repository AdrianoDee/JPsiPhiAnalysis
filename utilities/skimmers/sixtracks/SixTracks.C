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

  std::string outputString = "2mu4k_six_tree.root";
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

  outTree = new TTree("SixTrackSkimmedTree","SixTrackSkimmedTree");

  outTree->Branch("run", 	&out_run, 	"run/F");
  outTree->Branch("event", 	&out_event, 	"event/F");
  outTree->Branch("lumi", 	&out_lumi, 	"lumi/F");
  outTree->Branch("numPrimaryVertices", 	&out_numPrimaryVertices, 	"numPrimaryVertices/F");
  outTree->Branch("trigger", 	&out_trigger, 	"trigger/F");
  outTree->Branch("noSixCandidates", 	&out_noSixCandidates, 	"noSixCandidates/F");
  outTree->Branch("five_id", 	&out_five_id, 	"five_id/F");
  outTree->Branch("dimuon_id", 	&out_dimuon_id, 	"dimuon_id/F");
  outTree->Branch("p_id", 	&out_p_id, 	"p_id/F");
  outTree->Branch("m_id", 	&out_m_id, 	"m_id/F");
  outTree->Branch("t_id", 	&out_t_id, 	"t_id/F");
  outTree->Branch("f_id", 	&out_f_id, 	"f_id/F");
  outTree->Branch("six_p4", "TLorentzVector",     &six_p4);
  outTree->Branch("five_p4", "TLorentzVector",    &five_p4);
  outTree->Branch("dimuonditrk_p4", "TLorentzVector",     &dimuonditrk_p4);
  outTree->Branch("ditrack_p4", "TLorentzVector",         &ditrack_p4);
  outTree->Branch("dimuon_p4", "TLorentzVector",  &dimuon_p4);
  outTree->Branch("lowMuon_p4", "TLorentzVector",         &lowMuon_p4);
  outTree->Branch("highMuon_p4", "TLorentzVector",        &highMuon_p4);
  outTree->Branch("highKaon_p4", "TLorentzVector",        &highKaon_p4);
  outTree->Branch("lowKaon_p4", "TLorentzVector",         &lowKaon_p4);
  outTree->Branch("thirdKaon_p4", "TLorentzVector",       &thirdKaon_p4);
  outTree->Branch("fourthKaon_p4", "TLorentzVector",      &fourthKaon_p4);
  outTree->Branch("highPion_p4", "TLorentzVector",        &highPion_p4);
  outTree->Branch("lowPion_p4", "TLorentzVector",         &lowPion_p4);
  outTree->Branch("thirdPion_p4", "TLorentzVector",       &thirdPion_p4);
  outTree->Branch("fourthPion_p4", "TLorentzVector",      &fourthPion_p4);
  outTree->Branch("highProton_p4", "TLorentzVector",      &highProton_p4);
  outTree->Branch("lowProton_p4", "TLorentzVector",       &lowProton_p4);
  outTree->Branch("thirdProton_p4", "TLorentzVector",     &thirdProton_p4);
  outTree->Branch("fourthProton_p4", "TLorentzVector",    &fourthProton_p4);
  outTree->Branch("dimuonditrk_m", 	&out_dimuonditrk_m, 	"dimuonditrk_m/F");
  outTree->Branch("dimuonditrk_pt", 	&out_dimuonditrk_pt, 	"dimuonditrk_pt/F");
  outTree->Branch("dimuonditrk_eta", 	&out_dimuonditrk_eta, 	"dimuonditrk_eta/F");
  outTree->Branch("dimuonditrk_phi", 	&out_dimuonditrk_phi, 	"dimuonditrk_phi/F");
  outTree->Branch("dimuonditrk_p", 	&out_dimuonditrk_p, 	"dimuonditrk_p/F");
  outTree->Branch("dimuon_m", 	&out_dimuon_m, 	"dimuon_m/F");
  outTree->Branch("dimuon_pt", 	&out_dimuon_pt, 	"dimuon_pt/F");
  outTree->Branch("dimuon_eta", 	&out_dimuon_eta, 	"dimuon_eta/F");
  outTree->Branch("dimuon_phi", 	&out_dimuon_phi, 	"dimuon_phi/F");
  outTree->Branch("dimuon_p", 	&out_dimuon_p, 	"dimuon_p/F");
  outTree->Branch("highTrackMatch", 	&out_highTrackMatch, 	"highTrackMatch/F");
  outTree->Branch("lowTrackMatch", 	&out_lowTrackMatch, 	"lowTrackMatch/F");
  outTree->Branch("lowMuonMatch", 	&out_lowMuonMatch, 	"lowMuonMatch/F");
  outTree->Branch("highMuonMatch", 	&out_highMuonMatch, 	"highMuonMatch/F");
  outTree->Branch("thirdTrackMatch", 	&out_thirdTrackMatch, 	"thirdTrackMatch/F");
  outTree->Branch("fourthTrackMatch", 	&out_fourthTrackMatch, 	"fourthTrackMatch/F");
  outTree->Branch("ditrack_m", 	&out_ditrack_m, 	"ditrack_m/F");
  outTree->Branch("diTrackOne_pt", 	&out_diTrackOne_pt, 	"diTrackOne_pt/F");
  outTree->Branch("diTrackOne_eta", 	&out_diTrackOne_eta, 	"diTrackOne_eta/F");
  outTree->Branch("diTrackOne_phi", 	&out_diTrackOne_phi, 	"diTrackOne_phi/F");
  outTree->Branch("diTrackOne_p", 	&out_diTrackOne_p, 	"diTrackOne_p/F");
  outTree->Branch("diTrackTwo_pt", 	&out_diTrackTwo_pt, 	"diTrackTwo_pt/F");
  outTree->Branch("diTrackTwo_eta", 	&out_diTrackTwo_eta, 	"diTrackTwo_eta/F");
  outTree->Branch("diTrackTwo_phi", 	&out_diTrackTwo_phi, 	"diTrackTwo_phi/F");
  outTree->Branch("diTrackTwo_p", 	&out_diTrackTwo_p, 	"diTrackTwo_p/F");
  outTree->Branch("diTrackThree_pt", 	&out_diTrackThree_pt, 	"diTrackThree_pt/F");
  outTree->Branch("diTrackThree_eta", 	&out_diTrackThree_eta, 	"diTrackThree_eta/F");
  outTree->Branch("diTrackThree_phi", 	&out_diTrackThree_phi, 	"diTrackThree_phi/F");
  outTree->Branch("diTrackThree_p", 	&out_diTrackThree_p, 	"diTrackThree_p/F");
  outTree->Branch("diTrackFour_pt", 	&out_diTrackFour_pt, 	"diTrackFour_pt/F");
  outTree->Branch("diTrackFour_eta", 	&out_diTrackFour_eta, 	"diTrackFour_eta/F");
  outTree->Branch("diTrackFour_phi", 	&out_diTrackFour_phi, 	"diTrackFour_phi/F");
  outTree->Branch("diTrackFour_p", 	&out_diTrackFour_p, 	"diTrackFour_p/F");
  outTree->Branch("diTrackFive_pt", 	&out_diTrackFive_pt, 	"diTrackFive_pt/F");
  outTree->Branch("diTrackFive_eta", 	&out_diTrackFive_eta, 	"diTrackFive_eta/F");
  outTree->Branch("diTrackFive_phi", 	&out_diTrackFive_phi, 	"diTrackFive_phi/F");
  outTree->Branch("diTrackFive_p", 	&out_diTrackFive_p, 	"diTrackFive_p/F");
  outTree->Branch("diTrackSix_pt", 	&out_diTrackSix_pt, 	"diTrackSix_pt/F");
  outTree->Branch("diTrackSix_eta", 	&out_diTrackSix_eta, 	"diTrackSix_eta/F");
  outTree->Branch("diTrackSix_phi", 	&out_diTrackSix_phi, 	"diTrackSix_phi/F");
  outTree->Branch("diTrackSix_p", 	&out_diTrackSix_p, 	"diTrackSix_p/F");
  outTree->Branch("dimuonDiTrkOne_mmpp", 	&out_dimuonDiTrkOne_mmpp, 	"dimuonDiTrkOne_mmpp/F");
  outTree->Branch("dimuonDiTrkTwo_mmpp", 	&out_dimuonDiTrkTwo_mmpp, 	"dimuonDiTrkTwo_mmpp/F");
  outTree->Branch("dimuonDiTrkThree_mmpp", 	&out_dimuonDiTrkThree_mmpp, 	"dimuonDiTrkThree_mmpp/F");
  outTree->Branch("dimuonDiTrkFour_mmpp", 	&out_dimuonDiTrkFour_mmpp, 	"dimuonDiTrkFour_mmpp/F");
  outTree->Branch("dimuonDiTrkOne_mm", 	&out_dimuonDiTrkOne_mmkk, 	"dimuonDiTrkOne_mmkk/F");
  outTree->Branch("dimuonDiTrkTwo_mmkk", 	&out_dimuonDiTrkTwo_mmkk, 	"dimuonDiTrkTwo_mmkk/F");
  outTree->Branch("dimuonDiTrkThree_mmkk", 	&out_dimuonDiTrkThree_mmkk, 	"dimuonDiTrkThree_mmkk/F");
  outTree->Branch("dimuonDiTrkFour_mmkk", 	&out_dimuonDiTrkFour_mmkk, 	"dimuonDiTrkFour_mmkk/F");
  outTree->Branch("highMuon_pt", 	&out_highMuon_pt, 	"highMuon_pt/F");
  outTree->Branch("highMuon_eta", 	&out_highMuon_eta, 	"highMuon_eta/F");
  outTree->Branch("highMuon_phi", 	&out_highMuon_phi, 	"highMuon_phi/F");
  outTree->Branch("highMuon_charge", 	&out_highMuon_charge, 	"highMuon_charge/F");
  outTree->Branch("highMuon_dz", 	&out_highMuon_dz, 	"highMuon_dz/F");
  outTree->Branch("highMuon_dxy", 	&out_highMuon_dxy, 	"highMuon_dxy/F");
  outTree->Branch("lowMuon_pt", 	&out_lowMuon_pt, 	"lowMuon_pt/F");
  outTree->Branch("lowMuon_eta", 	&out_lowMuon_eta, 	"lowMuon_eta/F");
  outTree->Branch("lowMuon_phi", 	&out_lowMuon_phi, 	"lowMuon_phi/F");
  outTree->Branch("lowMuon_charge", 	&out_lowMuon_charge, 	"lowMuon_charge/F");
  outTree->Branch("lowMuon_dz", 	&out_lowMuon_dz, 	"lowMuon_dz/F");
  outTree->Branch("lowMuon_dxy", 	&out_lowMuon_dxy, 	"lowMuon_dxy/F");
  outTree->Branch("highTrack_pt", 	&out_highTrack_pt, 	"highTrack_pt/F");
  outTree->Branch("highTrack_eta", 	&out_highTrack_eta, 	"highTrack_eta/F");
  outTree->Branch("highTrack_phi", 	&out_highTrack_phi, 	"highTrack_phi/F");
  outTree->Branch("highTrack_charge", 	&out_highTrack_charge, 	"highTrack_charge/F");
  outTree->Branch("highTrack_dz", 	&out_highTrack_dz, 	"highTrack_dz/F");
  outTree->Branch("highTrack_dxy", 	&out_highTrack_dxy, 	"highTrack_dxy/F");
  outTree->Branch("lowTrack_pt", 	&out_lowTrack_pt, 	"lowTrack_pt/F");
  outTree->Branch("lowTrack_eta", 	&out_lowTrack_eta, 	"lowTrack_eta/F");
  outTree->Branch("lowTrack_phi", 	&out_lowTrack_phi, 	"lowTrack_phi/F");
  outTree->Branch("lowTrack_charge", 	&out_lowTrack_charge, 	"lowTrack_charge/F");
  outTree->Branch("lowTrack_dz", 	&out_lowTrack_dz, 	"lowTrack_dz/F");
  outTree->Branch("lowTrack_dxy", 	&out_lowTrack_dxy, 	"lowTrack_dxy/F");
  outTree->Branch("thirdTrack_pt", 	&out_thirdTrack_pt, 	"thirdTrack_pt/F");
  outTree->Branch("thirdTrack_eta", 	&out_thirdTrack_eta, 	"thirdTrack_eta/F");
  outTree->Branch("thirdTrack_phi", 	&out_thirdTrack_phi, 	"thirdTrack_phi/F");
  outTree->Branch("thirdTrack_charge", 	&out_thirdTrack_charge, 	"thirdTrack_charge/F");
  outTree->Branch("thirdTrack_dz", 	&out_thirdTrack_dz, 	"thirdTrack_dz/F");
  outTree->Branch("thirdTrack_dxy", 	&out_thirdTrack_dxy, 	"thirdTrack_dxy/F");
  outTree->Branch("dimuonDiTrkOne_pt", 	&out_dimuonDiTrkOne_pt, 	"dimuonDiTrkOne_pt/F");
  outTree->Branch("dimuonDiTrkOne_eta", 	&out_dimuonDiTrkOne_eta, 	"dimuonDiTrkOne_eta/F");
  outTree->Branch("dimuonDiTrkOne_phi", 	&out_dimuonDiTrkOne_phi, 	"dimuonDiTrkOne_phi/F");
  outTree->Branch("dimuonDiTrkOne_charge", 	&out_dimuonDiTrkOne_charge, 	"dimuonDiTrkOne_charge/F");
  outTree->Branch("dimuonDiTrkOne_p", 	&out_dimuonDiTrkOne_p, 	"dimuonDiTrkOne_p/F");
  outTree->Branch("dimuonDiTrkTwo_pt", 	&out_dimuonDiTrkTwo_pt, 	"dimuonDiTrkTwo_pt/F");
  outTree->Branch("dimuonDiTrkTwo_eta", 	&out_dimuonDiTrkTwo_eta, 	"dimuonDiTrkTwo_eta/F");
  outTree->Branch("dimuonDiTrkTwo_phi", 	&out_dimuonDiTrkTwo_phi, 	"dimuonDiTrkTwo_phi/F");
  outTree->Branch("dimuonDiTrkTwo_charge", 	&out_dimuonDiTrkTwo_charge, 	"dimuonDiTrkTwo_charge/F");
  outTree->Branch("dimuonDiTrkTwo_p", 	&out_dimuonDiTrkTwo_p, 	"dimuonDiTrkTwo_p/F");
  outTree->Branch("dimuonDiTrkThree_pt", 	&out_dimuonDiTrkThree_pt, 	"dimuonDiTrkThree_pt/F");
  outTree->Branch("dimuonDiTrkThree_eta", 	&out_dimuonDiTrkThree_eta, 	"dimuonDiTrkThree_eta/F");
  outTree->Branch("dimuonDiTrkThree_phi", 	&out_dimuonDiTrkThree_phi, 	"dimuonDiTrkThree_phi/F");
  outTree->Branch("dimuonDiTrkThree_charge", 	&out_dimuonDiTrkThree_charge, 	"dimuonDiTrkThree_charge/F");
  outTree->Branch("dimuonDiTrkThree_p", 	&out_dimuonDiTrkThree_p, 	"dimuonDiTrkThree_p/F");
  outTree->Branch("dimuonDiTrkFour_pt", 	&out_dimuonDiTrkFour_pt, 	"dimuonDiTrkFour_pt/F");
  outTree->Branch("dimuonDiTrkFour_eta", 	&out_dimuonDiTrkFour_eta, 	"dimuonDiTrkFour_eta/F");
  outTree->Branch("dimuonDiTrkFour_phi", 	&out_dimuonDiTrkFour_phi, 	"dimuonDiTrkFour_phi/F");
  outTree->Branch("dimuonDiTrkFour_charge", 	&out_dimuonDiTrkFour_charge, 	"dimuonDiTrkFour_charge/F");
  outTree->Branch("dimuonDiTrkFour_p", 	&out_dimuonDiTrkFour_p, 	"dimuonDiTrkFour_p/F");
  outTree->Branch("dimuonDiTrkFive_pt", 	&out_dimuonDiTrkFive_pt, 	"dimuonDiTrkFive_pt/F");
  outTree->Branch("dimuonDiTrkFive_eta", 	&out_dimuonDiTrkFive_eta, 	"dimuonDiTrkFive_eta/F");
  outTree->Branch("dimuonDiTrkFive_phi", 	&out_dimuonDiTrkFive_phi, 	"dimuonDiTrkFive_phi/F");
  outTree->Branch("dimuonDiTrkFive_charge", 	&out_dimuonDiTrkFive_charge, 	"dimuonDiTrkFive_charge/F");
  outTree->Branch("dimuonDiTrkFive_p", 	&out_dimuonDiTrkFive_p, 	"dimuonDiTrkFive_p/F");
  outTree->Branch("dimuonDiTrkSix_pt", 	&out_dimuonDiTrkSix_pt, 	"dimuonDiTrkSix_pt/F");
  outTree->Branch("dimuonDiTrkSix_eta", 	&out_dimuonDiTrkSix_eta, 	"dimuonDiTrkSix_eta/F");
  outTree->Branch("dimuonDiTrkSix_phi", 	&out_dimuonDiTrkSix_phi, 	"dimuonDiTrkSix_phi/F");
  outTree->Branch("dimuonDiTrkSix_charge", 	&out_dimuonDiTrkSix_charge, 	"dimuonDiTrkSix_charge/F");
  outTree->Branch("dimuonDiTrkSix_p", 	&out_dimuonDiTrkSix_p, 	"dimuonDiTrkSix_p/F");
  outTree->Branch("dimuon_vProb", 	&out_dimuon_vProb, 	"dimuon_vProb/F");
  outTree->Branch("dimuon_vChi2", 	&out_dimuon_vChi2, 	"dimuon_vChi2/F");
  outTree->Branch("dimuon_DCA", 	&out_dimuon_DCA, 	"dimuon_DCA/F");
  outTree->Branch("dimuon_ctauPV", 	&out_dimuon_ctauPV, 	"dimuon_ctauPV/F");
  outTree->Branch("dimuon_ctauErrPV", 	&out_dimuon_ctauErrPV, 	"dimuon_ctauErrPV/F");
  outTree->Branch("dimuon_cosAlpha", 	&out_dimuon_cosAlpha, 	"dimuon_cosAlpha/F");
  outTree->Branch("triTrack_m", 	&out_triTrack_m, 	"triTrack_m/F");
  outTree->Branch("triTrack_pt", 	&out_triTrack_pt, 	"triTrack_pt/F");
  outTree->Branch("triTrack_eta", 	&out_triTrack_eta, 	"triTrack_eta/F");
  outTree->Branch("triTrack_phi", 	&out_triTrack_phi, 	"triTrack_phi/F");
  outTree->Branch("triTrack_charge", 	&out_triTrack_charge, 	"triTrack_charge/F");
  outTree->Branch("dimuonditrk_vProb", 	&out_dimuonditrk_vProb, 	"dimuonditrk_vProb/F");
  outTree->Branch("dimuonditrk_vChi2", 	&out_dimuonditrk_vChi2, 	"dimuonditrk_vChi2/F");
  outTree->Branch("dimuonditrk_nDof", 	&out_dimuonditrk_nDof, 	"dimuonditrk_nDof/F");
  outTree->Branch("dimuonditrk_charge", 	&out_dimuonditrk_charge, 	"dimuonditrk_charge/F");
  outTree->Branch("dimuonditrk_cosAlpha", 	&out_dimuonditrk_cosAlpha, 	"dimuonditrk_cosAlpha/F");
  outTree->Branch("dimuonditrk_ctauPV", 	&out_dimuonditrk_ctauPV, 	"dimuonditrk_ctauPV/F");
  outTree->Branch("dimuonditrk_ctauErrPV", 	&out_dimuonditrk_ctauErrPV, 	"dimuonditrk_ctauErrPV/F");
  outTree->Branch("dimuonditrk_cosAlphaCA", 	&out_dimuonditrk_cosAlphaCA, 	"dimuonditrk_cosAlphaCA/F");
  outTree->Branch("dimuonditrk_ctauPVCA", 	&out_dimuonditrk_ctauPVCA, 	"dimuonditrk_ctauPVCA/F");
  outTree->Branch("dimuonditrk_ctauErrPVCA", 	&out_dimuonditrk_ctauErrPVCA, 	"dimuonditrk_ctauErrPVCA/F");
  outTree->Branch("dimuonditrk_cosAlphaDZ", 	&out_dimuonditrk_cosAlphaDZ, 	"dimuonditrk_cosAlphaDZ/F");
  outTree->Branch("dimuonditrk_ctauPVDZ", 	&out_dimuonditrk_ctauPVDZ, 	"dimuonditrk_ctauPVDZ/F");
  outTree->Branch("dimuonditrk_ctauErrPVDZ", 	&out_dimuonditrk_ctauErrPVDZ, 	"dimuonditrk_ctauErrPVDZ/F");
  outTree->Branch("dimuonditrk_cosAlphaBS", 	&out_dimuonditrk_cosAlphaBS, 	"dimuonditrk_cosAlphaBS/F");
  outTree->Branch("dimuonditrk_ctauPVBS", 	&out_dimuonditrk_ctauPVBS, 	"dimuonditrk_ctauPVBS/F");
  outTree->Branch("dimuonditrk_ctauErrPVBS", 	&out_dimuonditrk_ctauErrPVBS, 	"dimuonditrk_ctauErrPVBS/F");
  outTree->Branch("dimuonditrk_vx", 	&out_dimuonditrk_vx, 	"dimuonditrk_vx/F");
  outTree->Branch("dimuonditrk_vy", 	&out_dimuonditrk_vy, 	"dimuonditrk_vy/F");
  outTree->Branch("dimuonditrk_vz", 	&out_dimuonditrk_vz, 	"dimuonditrk_vz/F");
  outTree->Branch("dca_m1m2", 	&out_dca_m1m2, 	"dca_m1m2/F");
  outTree->Branch("dca_m1t1", 	&out_dca_m1t1, 	"dca_m1t1/F");
  outTree->Branch("dca_m1t2", 	&out_dca_m1t2, 	"dca_m1t2/F");
  outTree->Branch("dca_m2t1", 	&out_dca_m2t1, 	"dca_m2t1/F");
  outTree->Branch("dca_m2t2", 	&out_dca_m2t2, 	"dca_m2t2/F");
  outTree->Branch("dca_t1t2", 	&out_dca_t1t2, 	"dca_t1t2/F");
  outTree->Branch("dca_m1t3", 	&out_dca_m1t3, 	"dca_m1t3/F");
  outTree->Branch("dca_m2t3", 	&out_dca_m2t3, 	"dca_m2t3/F");
  outTree->Branch("dca_t1t3", 	&out_dca_t1t3, 	"dca_t1t3/F");
  outTree->Branch("dca_t2t3", 	&out_dca_t2t3, 	"dca_t2t3/F");
  outTree->Branch("dca_m1t4", 	&out_dca_m1t4, 	"dca_m1t4/F");
  outTree->Branch("dca_m2t4", 	&out_dca_m2t4, 	"dca_m2t4/F");
  outTree->Branch("dca_t1t4", 	&out_dca_t1t4, 	"dca_t1t4/F");
  outTree->Branch("dca_t2t4", 	&out_dca_t2t4, 	"dca_t2t4/F");
  outTree->Branch("dca_t3t4", 	&out_dca_t3t4, 	"dca_t3t4/F");
  outTree->Branch("highTrackMuonDR", 	&out_highTrackMuonDR, 	"highTrackMuonDR/F");
  outTree->Branch("highTrackMuonDP", 	&out_highTrackMuonDP, 	"highTrackMuonDP/F");
  outTree->Branch("highTrackMuonDPt", 	&out_highTrackMuonDPt, 	"highTrackMuonDPt/F");
  outTree->Branch("lowTrackMuonDR", 	&out_lowTrackMuonDR, 	"lowTrackMuonDR/F");
  outTree->Branch("lowTrackMuonDP", 	&out_lowTrackMuonDP, 	"lowTrackMuonDP/F");
  outTree->Branch("lowTrackMuonDPt", 	&out_lowTrackMuonDPt, 	"lowTrackMuonDPt/F");
  outTree->Branch("thirdTrackMuonDR", 	&out_thirdTrackMuonDR, 	"thirdTrackMuonDR/F");
  outTree->Branch("thirdTrackMuonDP", 	&out_thirdTrackMuonDP, 	"thirdTrackMuonDP/F");
  outTree->Branch("thirdTrackMuonDPt", 	&out_thirdTrackMuonDPt, 	"thirdTrackMuonDPt/F");
  outTree->Branch("fourthTrackMuonDR", 	&out_fourthTrackMuonDR, 	"fourthTrackMuonDR/F");
  outTree->Branch("fourthTrackMuonDP", 	&out_fourthTrackMuonDP, 	"fourthTrackMuonDP/F");
  outTree->Branch("fourthTrackMuonDPt", 	&out_fourthTrackMuonDPt, 	"fourthTrackMuonDPt/F");
  outTree->Branch("tPFromPV", 	&out_tPFromPV, 	"tPFromPV/F");
  outTree->Branch("tMFromPV", 	&out_tMFromPV, 	"tMFromPV/F");
  outTree->Branch("tTFromPV", 	&out_tTFromPV, 	"tTFromPV/F");
  outTree->Branch("tFFromPV", 	&out_tFFromPV, 	"tFFromPV/F");
  outTree->Branch("tPFromPVCA", 	&out_tPFromPVCA, 	"tPFromPVCA/F");
  outTree->Branch("tMFromPVCA", 	&out_tMFromPVCA, 	"tMFromPVCA/F");
  outTree->Branch("tTFromPVCA", 	&out_tTFromPVCA, 	"tTFromPVCA/F");
  outTree->Branch("tFFromPVCA", 	&out_tFFromPVCA, 	"tFFromPVCA/F");
  outTree->Branch("tPFromPVDZ", 	&out_tPFromPVDZ, 	"tPFromPVDZ/F");
  outTree->Branch("tMFromPVDZ", 	&out_tMFromPVDZ, 	"tMFromPVDZ/F");
  outTree->Branch("tTFromPVDZ", 	&out_tTFromPVDZ, 	"tTFromPVDZ/F");
  outTree->Branch("tFFromPVDZ", 	&out_tFFromPVDZ, 	"tFFromPVDZ/F");
  outTree->Branch("five_m", 	&out_five_m, 	"five_m/F");
  outTree->Branch("five_m_ref", 	&out_five_m_ref, 	"five_m_ref/F");
  outTree->Branch("five_mass_ppk", 	&out_five_mass_ppk, 	"five_mass_ppk/F");
  outTree->Branch("five_mass_kpp", 	&out_five_mass_kpp, 	"five_mass_kpp/F");
  outTree->Branch("five_mass_pkp", 	&out_five_mass_pkp, 	"five_mass_pkp/F");
  outTree->Branch("five_mass_ppp", 	&out_five_mass_ppp, 	"five_mass_ppp/F");
  outTree->Branch("fiveOne_pt", 	&out_fiveOne_pt, 	"fiveOne_pt/F");
  outTree->Branch("fiveOne_eta", 	&out_fiveOne_eta, 	"fiveOne_eta/F");
  outTree->Branch("fiveOne_phi", 	&out_fiveOne_phi, 	"fiveOne_phi/F");
  outTree->Branch("fiveOne_p", 	&out_fiveOne_p, 	"fiveOne_p/F");
  outTree->Branch("fiveTwo_pt", 	&out_fiveTwo_pt, 	"fiveTwo_pt/F");
  outTree->Branch("fiveTwo_eta", 	&out_fiveTwo_eta, 	"fiveTwo_eta/F");
  outTree->Branch("fiveTwo_phi", 	&out_fiveTwo_phi, 	"fiveTwo_phi/F");
  outTree->Branch("fiveTwo_p", 	&out_fiveTwo_p, 	"fiveTwo_p/F");
  outTree->Branch("fiveThree_pt", 	&out_fiveThree_pt, 	"fiveThree_pt/F");
  outTree->Branch("fiveThree_eta", 	&out_fiveThree_eta, 	"fiveThree_eta/F");
  outTree->Branch("fiveThree_phi", 	&out_fiveThree_phi, 	"fiveThree_phi/F");
  outTree->Branch("fiveThree_p", 	&out_fiveThree_p, 	"fiveThree_p/F");
  outTree->Branch("fiveFour_pt", 	&out_fiveFour_pt, 	"fiveFour_pt/F");
  outTree->Branch("fiveFour_eta", 	&out_fiveFour_eta, 	"fiveFour_eta/F");
  outTree->Branch("fiveFour_phi", 	&out_fiveFour_phi, 	"fiveFour_phi/F");
  outTree->Branch("fiveFour_p", 	&out_fiveFour_p, 	"fiveFour_p/F");
  outTree->Branch("fiveFive_pt", 	&out_fiveFive_pt, 	"fiveFive_pt/F");
  outTree->Branch("fiveFive_eta", 	&out_fiveFive_eta, 	"fiveFive_eta/F");
  outTree->Branch("fiveFive_phi", 	&out_fiveFive_phi, 	"fiveFive_phi/F");
  outTree->Branch("fiveFive_p", 	&out_fiveFive_p, 	"fiveFive_p/F");
  outTree->Branch("five_cosAlpha", 	&out_five_cosAlpha, 	"five_cosAlpha/F");
  outTree->Branch("five_ctauPV", 	&out_five_ctauPV, 	"five_ctauPV/F");
  outTree->Branch("five_ctauErrPV", 	&out_five_ctauErrPV, 	"five_ctauErrPV/F");
  outTree->Branch("five_cosAlphaCA", 	&out_five_cosAlphaCA, 	"five_cosAlphaCA/F");
  outTree->Branch("five_ctauPVCA", 	&out_five_ctauPVCA, 	"five_ctauPVCA/F");
  outTree->Branch("five_ctauErrPVCA", 	&out_five_ctauErrPVCA, 	"five_ctauErrPVCA/F");
  outTree->Branch("five_cosAlphaDZ", 	&out_five_cosAlphaDZ, 	"five_cosAlphaDZ/F");
  outTree->Branch("five_ctauPVDZ", 	&out_five_ctauPVDZ, 	"five_ctauPVDZ/F");
  outTree->Branch("five_ctauErrPVDZ", 	&out_five_ctauErrPVDZ, 	"five_ctauErrPVDZ/F");
  outTree->Branch("five_cosAlphaBS", 	&out_five_cosAlphaBS, 	"five_cosAlphaBS/F");
  outTree->Branch("five_ctauPVBS", 	&out_five_ctauPVBS, 	"five_ctauPVBS/F");
  outTree->Branch("five_ctauErrPVBS", 	&out_five_ctauErrPVBS, 	"five_ctauErrPVBS/F");
  outTree->Branch("five_vProb", 	&out_five_vProb, 	"five_vProb/F");
  outTree->Branch("five_nDof", 	&out_five_nDof, 	"five_nDof/F");
  outTree->Branch("five_vChi2", 	&out_five_vChi2, 	"five_vChi2/F");
  outTree->Branch("five_vx", 	&out_five_vx, 	"five_vx/F");
  outTree->Branch("five_vy", 	&out_five_vy, 	"five_vy/F");
  outTree->Branch("five_vz", 	&out_five_vz, 	"five_vz/F");
  outTree->Branch("five_charge", 	&out_five_charge, 	"five_charge/F");
  outTree->Branch("bestPV_X", 	&out_bestPV_X, 	"bestPV_X/F");
  outTree->Branch("bestPV_Y", 	&out_bestPV_Y, 	"bestPV_Y/F");
  outTree->Branch("bestPV_Z", 	&out_bestPV_Z, 	"bestPV_Z/F");
  outTree->Branch("cosAlphaPV_X", 	&out_cosAlphaPV_X, 	"cosAlphaPV_X/F");
  outTree->Branch("cosAlphaPV_Y", 	&out_cosAlphaPV_Y, 	"cosAlphaPV_Y/F");
  outTree->Branch("cosAlphaPV_Z", 	&out_cosAlphaPV_Z, 	"cosAlphaPV_Z/F");
  outTree->Branch("bS_X", 	&out_bS_X, 	"bS_X/F");
  outTree->Branch("bS_Y", 	&out_bS_Y, 	"bS_Y/F");
  outTree->Branch("bS_Z", 	&out_bS_Z, 	"bS_Z/F");
  outTree->Branch("zPV_X", 	&out_zPV_X, 	"zPV_X/F");
  outTree->Branch("zPV_Y", 	&out_zPV_Y, 	"zPV_Y/F");
  outTree->Branch("zPV_Z", 	&out_zPV_Z, 	"zPV_Z/F");
  outTree->Branch("lowMuon_isTight", 	&out_lowMuon_isTight, 	"lowMuon_isTight/F");
  outTree->Branch("lowMuon_isLoose", 	&out_lowMuon_isLoose, 	"lowMuon_isLoose/F");
  outTree->Branch("lowMuon_isSoft", 	&out_lowMuon_isSoft, 	"lowMuon_isSoft/F");
  outTree->Branch("lowMuon_isMedium", 	&out_lowMuon_isMedium, 	"lowMuon_isMedium/F");
  outTree->Branch("lowMuon_isHighPt", 	&out_lowMuon_isHighPt, 	"lowMuon_isHighPt/F");
  outTree->Branch("lowMuon_isTracker", 	&out_lowMuon_isTracker, 	"lowMuon_isTracker/F");
  outTree->Branch("lowMuon_isGlobal", 	&out_lowMuon_isGlobal, 	"lowMuon_isGlobal/F");
  outTree->Branch("lowMuon_NPixelHits", 	&out_lowMuon_NPixelHits, 	"lowMuon_NPixelHits/F");
  outTree->Branch("lowMuon_NStripHits", 	&out_lowMuon_NStripHits, 	"lowMuon_NStripHits/F");
  outTree->Branch("lowMuon_NTrackhits", 	&out_lowMuon_NTrackhits, 	"lowMuon_NTrackhits/F");
  outTree->Branch("lowMuon_NBPixHits", 	&out_lowMuon_NBPixHits, 	"lowMuon_NBPixHits/F");
  outTree->Branch("lowMuon_NPixLayers", 	&out_lowMuon_NPixLayers, 	"lowMuon_NPixLayers/F");
  outTree->Branch("lowMuon_NTraLayers", 	&out_lowMuon_NTraLayers, 	"lowMuon_NTraLayers/F");
  outTree->Branch("lowMuon_NStrLayers", 	&out_lowMuon_NStrLayers, 	"lowMuon_NStrLayers/F");
  outTree->Branch("lowMuon_NBPixLayers", 	&out_lowMuon_NBPixLayers, 	"lowMuon_NBPixLayers/F");
  outTree->Branch("highMuon_isTight", 	&out_highMuon_isTight, 	"highMuon_isTight/F");
  outTree->Branch("highMuon_isLoose", 	&out_highMuon_isLoose, 	"highMuon_isLoose/F");
  outTree->Branch("highMuon_isSoft", 	&out_highMuon_isSoft, 	"highMuon_isSoft/F");
  outTree->Branch("highMuon_isMedium", 	&out_highMuon_isMedium, 	"highMuon_isMedium/F");
  outTree->Branch("highMuon_isHighPt", 	&out_highMuon_isHighPt, 	"highMuon_isHighPt/F");
  outTree->Branch("highMuon_isTracker", 	&out_highMuon_isTracker, 	"highMuon_isTracker/F");
  outTree->Branch("highMuon_isGlobal", 	&out_highMuon_isGlobal, 	"highMuon_isGlobal/F");
  outTree->Branch("highMuon_NPixelHits", 	&out_highMuon_NPixelHits, 	"highMuon_NPixelHits/F");
  outTree->Branch("highMuon_NStripHits", 	&out_highMuon_NStripHits, 	"highMuon_NStripHits/F");
  outTree->Branch("highMuon_NTrackhits", 	&out_highMuon_NTrackhits, 	"highMuon_NTrackhits/F");
  outTree->Branch("highMuon_NBPixHits", 	&out_highMuon_NBPixHits, 	"highMuon_NBPixHits/F");
  outTree->Branch("highMuon_NPixLayers", 	&out_highMuon_NPixLayers, 	"highMuon_NPixLayers/F");
  outTree->Branch("highMuon_NTraLayers", 	&out_highMuon_NTraLayers, 	"highMuon_NTraLayers/F");
  outTree->Branch("highMuon_NStrLayers", 	&out_highMuon_NStrLayers, 	"highMuon_NStrLayers/F");
  outTree->Branch("highMuon_NBPixLayers", 	&out_highMuon_NBPixLayers, 	"highMuon_NBPixLayers/F");
  outTree->Branch("lowMuon_type", 	&out_lowMuon_type, 	"lowMuon_type/F");
  outTree->Branch("highMuon_type", 	&out_highMuon_type, 	"highMuon_type/F");
  outTree->Branch("highTrack_NPixelHits", 	&out_highTrack_NPixelHits, 	"highTrack_NPixelHits/F");
  outTree->Branch("highTrack_NStripHits", 	&out_highTrack_NStripHits, 	"highTrack_NStripHits/F");
  outTree->Branch("highTrack_NTrackhits", 	&out_highTrack_NTrackhits, 	"highTrack_NTrackhits/F");
  outTree->Branch("highTrack_NBPixHits", 	&out_highTrack_NBPixHits, 	"highTrack_NBPixHits/F");
  outTree->Branch("highTrack_NPixLayers", 	&out_highTrack_NPixLayers, 	"highTrack_NPixLayers/F");
  outTree->Branch("highTrack_NTraLayers", 	&out_highTrack_NTraLayers, 	"highTrack_NTraLayers/F");
  outTree->Branch("highTrack_NStrLayers", 	&out_highTrack_NStrLayers, 	"highTrack_NStrLayers/F");
  outTree->Branch("highTrack_NBPixLayers", 	&out_highTrack_NBPixLayers, 	"highTrack_NBPixLayers/F");
  outTree->Branch("lowTrack_NPixelHits", 	&out_lowTrack_NPixelHits, 	"lowTrack_NPixelHits/F");
  outTree->Branch("lowTrack_NStripHits", 	&out_lowTrack_NStripHits, 	"lowTrack_NStripHits/F");
  outTree->Branch("lowTrack_NTrackhits", 	&out_lowTrack_NTrackhits, 	"lowTrack_NTrackhits/F");
  outTree->Branch("lowTrack_NBPixHits", 	&out_lowTrack_NBPixHits, 	"lowTrack_NBPixHits/F");
  outTree->Branch("lowTrack_NPixLayers", 	&out_lowTrack_NPixLayers, 	"lowTrack_NPixLayers/F");
  outTree->Branch("lowTrack_NTraLayers", 	&out_lowTrack_NTraLayers, 	"lowTrack_NTraLayers/F");
  outTree->Branch("lowTrack_NStrLayers", 	&out_lowTrack_NStrLayers, 	"lowTrack_NStrLayers/F");
  outTree->Branch("lowTrack_NBPixLayers", 	&out_lowTrack_NBPixLayers, 	"lowTrack_NBPixLayers/F");
  outTree->Branch("thirdTrack_NPixelHits", 	&out_thirdTrack_NPixelHits, 	"thirdTrack_NPixelHits/F");
  outTree->Branch("thirdTrack_NStripHits", 	&out_thirdTrack_NStripHits, 	"thirdTrack_NStripHits/F");
  outTree->Branch("thirdTrack_NTrackhits", 	&out_thirdTrack_NTrackhits, 	"thirdTrack_NTrackhits/F");
  outTree->Branch("thirdTrack_NBPixHits", 	&out_thirdTrack_NBPixHits, 	"thirdTrack_NBPixHits/F");
  outTree->Branch("thirdTrack_NPixLayers", 	&out_thirdTrack_NPixLayers, 	"thirdTrack_NPixLayers/F");
  outTree->Branch("thirdTrack_NTraLayers", 	&out_thirdTrack_NTraLayers, 	"thirdTrack_NTraLayers/F");
  outTree->Branch("thirdTrack_NStrLayers", 	&out_thirdTrack_NStrLayers, 	"thirdTrack_NStrLayers/F");
  outTree->Branch("thirdTrack_NBPixLayers", 	&out_thirdTrack_NBPixLayers, 	"thirdTrack_NBPixLayers/F");
  outTree->Branch("fourthTrack_NPixLayers", 	&out_fourthTrack_NPixLayers, 	"fourthTrack_NPixLayers/F");
  outTree->Branch("fourthTrack_NTraLayers", 	&out_fourthTrack_NTraLayers, 	"fourthTrack_NTraLayers/F");
  outTree->Branch("fourthTrack_NStrLayers", 	&out_fourthTrack_NStrLayers, 	"fourthTrack_NStrLayers/F");
  outTree->Branch("fourthTrack_NBPixLayers", 	&out_fourthTrack_NBPixLayers, 	"fourthTrack_NBPixLayers/F");
  outTree->Branch("fourthTrack_NPixelHits", 	&out_fourthTrack_NPixelHits, 	"fourthTrack_NPixelHits/F");
  outTree->Branch("fourthTrack_NStripHits", 	&out_fourthTrack_NStripHits, 	"fourthTrack_NStripHits/F");
  outTree->Branch("fourthTrack_NTrackhits", 	&out_fourthTrack_NTrackhits, 	"fourthTrack_NTrackhits/F");
  outTree->Branch("fourthTrack_NBPixHits", 	&out_fourthTrack_NBPixHits, 	"fourthTrack_NBPixHits/F");
  outTree->Branch("six_m", 	&out_six_m, 	"six_m/F");
  outTree->Branch("six_m_ref", 	&out_six_m_ref, 	"six_m_ref/F");
  outTree->Branch("six_mass_ppkk", 	&out_six_mass_ppkk, 	"six_mass_ppkk/F");
  outTree->Branch("six_mass_pkpk", 	&out_six_mass_pkpk, 	"six_mass_pkpk/F");
  outTree->Branch("six_mass_pkkk", 	&out_six_mass_pkkk, 	"six_mass_pkkk/F");
  outTree->Branch("six_mass_kpkp", 	&out_six_mass_kpkp, 	"six_mass_kpkp/F");
  outTree->Branch("six_mass_kppk", 	&out_six_mass_kppk, 	"six_mass_kppk/F");
  outTree->Branch("six_mass_kkkk", 	&out_six_mass_kkkk, 	"six_mass_kkkk/F");
  outTree->Branch("six_pt", 	&out_six_pt, 	"six_pt/F");
  outTree->Branch("six_eta", 	&out_six_eta, 	"six_eta/F");
  outTree->Branch("six_phi", 	&out_six_phi, 	"six_phi/F");
  outTree->Branch("six_p", 	&out_six_p, 	"six_p/F");
  outTree->Branch("six_cosAlpha", 	&out_six_cosAlpha, 	"six_cosAlpha/F");
  outTree->Branch("six_ctauPV", 	&out_six_ctauPV, 	"six_ctauPV/F");
  outTree->Branch("six_ctauErrPV", 	&out_six_ctauErrPV, 	"six_ctauErrPV/F");
  outTree->Branch("six_cosAlphaCA", 	&out_six_cosAlphaCA, 	"six_cosAlphaCA/F");
  outTree->Branch("six_ctauPVCA", 	&out_six_ctauPVCA, 	"six_ctauPVCA/F");
  outTree->Branch("six_ctauErrPVCA", 	&out_six_ctauErrPVCA, 	"six_ctauErrPVCA/F");
  outTree->Branch("six_cosAlphaDZ", 	&out_six_cosAlphaDZ, 	"six_cosAlphaDZ/F");
  outTree->Branch("six_ctauPVDZ", 	&out_six_ctauPVDZ, 	"six_ctauPVDZ/F");
  outTree->Branch("six_ctauErrPVDZ", 	&out_six_ctauErrPVDZ, 	"six_ctauErrPVDZ/F");
  outTree->Branch("six_cosAlphaBS", 	&out_six_cosAlphaBS, 	"six_cosAlphaBS/F");
  outTree->Branch("six_ctauPVBS", 	&out_six_ctauPVBS, 	"six_ctauPVBS/F");
  outTree->Branch("six_ctauErrPVBS", 	&out_six_ctauErrPVBS, 	"six_ctauErrPVBS/F");
  outTree->Branch("six_vProb", 	&out_six_vProb, 	"six_vProb/F");
  outTree->Branch("six_nDof", 	&out_six_nDof, 	"six_nDof/F");
  outTree->Branch("six_vChi2", 	&out_six_vChi2, 	"six_vChi2/F");
  outTree->Branch("six_vx", 	&out_six_vx, 	"six_vx/F");
  outTree->Branch("six_vy", 	&out_six_vy, 	"six_vy/F");
  outTree->Branch("six_vz", 	&out_six_vz, 	"six_vz/F");
  outTree->Branch("six_charge", 	&out_six_charge, 	"six_charge/F");

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

  bool test = true;

  // test = test && ((*highMuon_pt) >= 1.0) && ((*highMuon_pt) >= 1.0);
  //
  // test = test && (*lowMuonMatch>0.0) && (*highMuonMatch>0.0);
  //
  // test = test && (*dimuonditrk_vProb> 0.01);

  //int a = (int) (*trigger);
  //std::cout << (*trigger);

  if(test)
  {

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
    out_six_p4 =    (*six_p4);
    out_five_p4 =   (*five_p4);
    out_dimuonditrk_p4 =    (*dimuonditrk_p4);
    out_ditrack_p4 =        (*ditrack_p4);
    out_dimuon_p4 =         (*dimuon_p4);
    out_lowMuon_p4 =        (*lowMuon_p4);
    out_highMuon_p4 =       (*highMuon_p4);
    out_highKaon_p4 =       (*highKaon_p4);
    out_lowKaon_p4 =        (*lowKaon_p4);
    out_thirdKaon_p4 =      (*thirdKaon_p4);
    out_fourthKaon_p4 =     (*fourthKaon_p4);
    out_highPion_p4 =       (*highPion_p4);
    out_lowPion_p4 =        (*lowPion_p4);
    out_thirdPion_p4 =      (*thirdPion_p4);
    out_fourthPion_p4 =     (*fourthPion_p4);
    out_highProton_p4 =     (*highProton_p4);
    out_lowProton_p4 =      (*lowProton_p4);
    out_thirdProton_p4 =    (*thirdProton_p4);
    out_fourthProton_p4 =   (*fourthProton_p4);
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

    // out_dimuonDiTrkOne_mmkk = 	(Float_t)(*dimuonDiTrkOne_mmkk);
    // out_dimuonDiTrkTwo_mmkk = 	(Float_t)(*dimuonDiTrkTwo_mmkk);
    // out_dimuonDiTrkThree_mmkk = 	(Float_t)(*dimuonDiTrkThree_mmkk);
    // out_dimuonDiTrkFour_mmkk = 	(Float_t)(*dimuonDiTrkFour_mmkk);

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
    //
    // out_tPFromPV = 	(Float_t)(*tPFromPV);
    // out_tMFromPV = 	(Float_t)(*tMFromPV);
    // out_tTFromPV = 	(Float_t)(*tTFromPV);
    // out_tFFromPV = 	(Float_t)(*tFFromPV);
    // out_tPFromPVCA = 	(Float_t)(*tPFromPVCA);
    // out_tMFromPVCA = 	(Float_t)(*tMFromPVCA);
    // out_tTFromPVCA = 	(Float_t)(*tTFromPVCA);
    // out_tFFromPVCA = 	(Float_t)(*tFFromPVCA);
    // out_tPFromPVDZ = 	(Float_t)(*tPFromPVDZ);
    // out_tMFromPVDZ = 	(Float_t)(*tMFromPVDZ);
    // out_tTFromPVDZ = 	(Float_t)(*tTFromPVDZ);
    // out_tFFromPVDZ = 	(Float_t)(*tFFromPVDZ);
    //
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

    outTree->Fill();

  }

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
