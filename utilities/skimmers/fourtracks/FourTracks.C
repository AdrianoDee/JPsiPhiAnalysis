#define FourTracks_cxx
// The class definition in FourTracks.h has been generated automatically
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
// root> T->Process("FourTracks.C")
// root> T->Process("FourTracks.C","some options")
// root> T->Process("FourTracks.C+")
//


#include "FourTracks.h"
#include <TH2.h>
#include <TStyle.h>

void FourTracks::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void FourTracks::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   std::string outputString = "2mu4k_four_tree.root";
   OutFile = new TProofOutputFile( outputString.data() );
   fOut = OutFile->OpenFile("RECREATE");
   if (!(fOut=OutFile->OpenFile("RECREATE")))
   {
     Warning("SlaveBegin","Problems opening file: %s%s", OutFile->GetDir(), OutFile->GetFileName() );
   }

   outTree = new TTree("SixTrackSkimmedTree","SixTrackSkimmedTree");

   outTree->Branch("run", 	&out_run, 	"run/F");
   outTree->Branch("event", 	&out_event, 	"event/F");
   outTree->Branch("lumi", 	&out_lumi, 	"lumi/F");
   outTree->Branch("numPrimaryVertices", 	&out_numPrimaryVertices, 	"numPrimaryVertices/F");
   outTree->Branch("trigger", 	&out_trigger, 	"trigger/F");
   outTree->Branch("noXCandidates", 	&out_noXCandidates, 	"noXCandidates/F");
   outTree->Branch("dimuon_id", 	&out_dimuon_id, 	"dimuon_id/F");
   outTree->Branch("p_id", 	&out_p_id, 	"p_id/F");
   outTree->Branch("m_id", 	&out_m_id, 	"m_id/F");
   // outTree->Branch("dimuonditrk_p4", 	&out_dimuonditrk_p4, 	"dimuonditrk_p4/F");
   // outTree->Branch("ditrack_p4", 	&out_ditrack_p4, 	"ditrack_p4/F");
   // outTree->Branch("dimuon_p4", 	&out_dimuon_p4, 	"dimuon_p4/F");
   // outTree->Branch("lowMuon_p4", 	&out_lowMuon_p4, 	"lowMuon_p4/F");
   // outTree->Branch("highMuon_p4", 	&out_highMuon_p4, 	"highMuon_p4/F");
   // outTree->Branch("highKaon_p4", 	&out_highKaon_p4, 	"highKaon_p4/F");
   // outTree->Branch("lowKaon_p4", 	&out_lowKaon_p4, 	"lowKaon_p4/F");
   // outTree->Branch("dimuonditrk_rf_p4", 	&out_dimuonditrk_rf_p4, 	"dimuonditrk_rf_p4/F");
   // outTree->Branch("ditrack_rf_p4", 	&out_ditrack_rf_p4, 	"ditrack_rf_p4/F");
   // outTree->Branch("dimuon_rf_p4", 	&out_dimuon_rf_p4, 	"dimuon_rf_p4/F");
   // outTree->Branch("lowMuon_rf_p4", 	&out_lowMuon_rf_p4, 	"lowMuon_rf_p4/F");
   // outTree->Branch("highMuon_rf_p4", 	&out_highMuon_rf_p4, 	"highMuon_rf_p4/F");
   // outTree->Branch("kaonp_rf_p4", 	&out_kaonp_rf_p4, 	"kaonp_rf_p4/F");
   // outTree->Branch("kaonn_rf_p4", 	&out_kaonn_rf_p4, 	"kaonn_rf_p4/F");
   outTree->Branch("dimuonditrk_m", 	&out_dimuonditrk_m, 	"dimuonditrk_m/F");
   outTree->Branch("dimuonditrk_m_rf", 	&out_dimuonditrk_m_rf, 	"dimuonditrk_m_rf/F");
   outTree->Branch("dimuonditrk_m_rf_d_c", 	&out_dimuonditrk_m_rf_d_c, 	"dimuonditrk_m_rf_d_c/F");
   outTree->Branch("dimuonditrk_m_rf_c", 	&out_dimuonditrk_m_rf_c, 	"dimuonditrk_m_rf_c/F");
   outTree->Branch("dimuonditrk_pt", 	&out_dimuonditrk_pt, 	"dimuonditrk_pt/F");
   outTree->Branch("dimuonditrk_eta", 	&out_dimuonditrk_eta, 	"dimuonditrk_eta/F");
   outTree->Branch("dimuonditrk_phi", 	&out_dimuonditrk_phi, 	"dimuonditrk_phi/F");
   outTree->Branch("dimuonditrk_y", 	&out_dimuonditrk_y, 	"dimuonditrk_y/F");
   outTree->Branch("dimuonditrk_vx", 	&out_dimuonditrk_vx, 	"dimuonditrk_vx/F");
   outTree->Branch("dimuonditrk_vy", 	&out_dimuonditrk_vy, 	"dimuonditrk_vy/F");
   outTree->Branch("dimuonditrk_vz", 	&out_dimuonditrk_vz, 	"dimuonditrk_vz/F");
   outTree->Branch("pv_x", 	&out_pv_x, 	"pv_x/F");
   outTree->Branch("pv_y", 	&out_pv_y, 	"pv_y/F");
   outTree->Branch("pv_z", 	&out_pv_z, 	"pv_z/F");
   outTree->Branch("dimuon_m", 	&out_dimuon_m, 	"dimuon_m/F");
   outTree->Branch("dimuon_pt", 	&out_dimuon_pt, 	"dimuon_pt/F");
   outTree->Branch("ditrack_m", 	&out_ditrack_m, 	"ditrack_m/F");
   outTree->Branch("ditrack_pt", 	&out_ditrack_pt, 	"ditrack_pt/F");
   outTree->Branch("highKaon_pt", 	&out_highKaon_pt, 	"highKaon_pt/F");
   outTree->Branch("lowKaon_pt", 	&out_lowKaon_pt, 	"lowKaon_pt/F");
   outTree->Branch("highMuon_pt", 	&out_highMuon_pt, 	"highMuon_pt/F");
   outTree->Branch("lowMuon_pt", 	&out_lowMuon_pt, 	"lowMuon_pt/F");
   outTree->Branch("highKaon_eta", 	&out_highKaon_eta, 	"highKaon_eta/F");
   outTree->Branch("lowKaon_eta", 	&out_lowKaon_eta, 	"lowKaon_eta/F");
   outTree->Branch("highMuon_eta", 	&out_highMuon_eta, 	"highMuon_eta/F");
   outTree->Branch("lowMuon_eta", 	&out_lowMuon_eta, 	"lowMuon_eta/F");
   outTree->Branch("highKaon_phi", 	&out_highKaon_phi, 	"highKaon_phi/F");
   outTree->Branch("lowKaon_phi", 	&out_lowKaon_phi, 	"lowKaon_phi/F");
   outTree->Branch("highMuon_phi", 	&out_highMuon_phi, 	"highMuon_phi/F");
   outTree->Branch("lowMuon_phi", 	&out_lowMuon_phi, 	"lowMuon_phi/F");
   outTree->Branch("highKaon_dxy", 	&out_highKaon_dxy, 	"highKaon_dxy/F");
   outTree->Branch("lowKaon_dxy", 	&out_lowKaon_dxy, 	"lowKaon_dxy/F");
   outTree->Branch("highMuon_dxy", 	&out_highMuon_dxy, 	"highMuon_dxy/F");
   outTree->Branch("lowMuon_dxy", 	&out_lowMuon_dxy, 	"lowMuon_dxy/F");
   outTree->Branch("highKaon_dz", 	&out_highKaon_dz, 	"highKaon_dz/F");
   outTree->Branch("lowKaon_dz", 	&out_lowKaon_dz, 	"lowKaon_dz/F");
   outTree->Branch("highMuon_dz", 	&out_highMuon_dz, 	"highMuon_dz/F");
   outTree->Branch("lowMuon_dz", 	&out_lowMuon_dz, 	"lowMuon_dz/F");
   outTree->Branch("highKaonMuonDR", 	&out_highKaonMuonDR, 	"highKaonMuonDR/F");
   outTree->Branch("highKaonMuonDP", 	&out_highKaonMuonDP, 	"highKaonMuonDP/F");
   outTree->Branch("highKaonMuonDPt", 	&out_highKaonMuonDPt, 	"highKaonMuonDPt/F");
   outTree->Branch("lowKaonMuonDR", 	&out_lowKaonMuonDR, 	"lowKaonMuonDR/F");
   outTree->Branch("lowKaonMuonDP", 	&out_lowKaonMuonDP, 	"lowKaonMuonDP/F");
   outTree->Branch("lowKaonMuonDPt", 	&out_lowKaonMuonDPt, 	"lowKaonMuonDPt/F");
   outTree->Branch("dimuonditrk_refPK_mass", 	&out_dimuonditrk_refPK_mass, 	"dimuonditrk_refPK_mass/F");
   outTree->Branch("dimuonditrk_refKP_mass", 	&out_dimuonditrk_refKP_mass, 	"dimuonditrk_refKP_mass/F");
   outTree->Branch("dimuonditrk_refPP_mass", 	&out_dimuonditrk_refPP_mass, 	"dimuonditrk_refPP_mass/F");
   outTree->Branch("dimuonditrk_refPK_vChi2", 	&out_dimuonditrk_refPK_vChi2, 	"dimuonditrk_refPK_vChi2/F");
   outTree->Branch("dimuonditrk_refKP_vChi2", 	&out_dimuonditrk_refKP_vChi2, 	"dimuonditrk_refKP_vChi2/F");
   outTree->Branch("dimuonditrk_refPP_vChi2", 	&out_dimuonditrk_refPP_vChi2, 	"dimuonditrk_refPP_vChi2/F");
   outTree->Branch("dimuonditrk_refPK_nDof", 	&out_dimuonditrk_refPK_nDof, 	"dimuonditrk_refPK_nDof/F");
   outTree->Branch("dimuonditrk_refKP_nDof", 	&out_dimuonditrk_refKP_nDof, 	"dimuonditrk_refKP_nDof/F");
   outTree->Branch("dimuonditrk_refPP_nDof", 	&out_dimuonditrk_refPP_nDof, 	"dimuonditrk_refPP_nDof/F");
   outTree->Branch("dimuonditrk_refPK_vProb", 	&out_dimuonditrk_refPK_vProb, 	"dimuonditrk_refPK_vProb/F");
   outTree->Branch("dimuonditrk_refKP_vProb", 	&out_dimuonditrk_refKP_vProb, 	"dimuonditrk_refKP_vProb/F");
   outTree->Branch("dimuonditrk_refPP_vProb", 	&out_dimuonditrk_refPP_vProb, 	"dimuonditrk_refPP_vProb/F");
   outTree->Branch("dimuon_vProb", 	&out_dimuon_vProb, 	"dimuon_vProb/F");
   outTree->Branch("dimuon_vChi2", 	&out_dimuon_vChi2, 	"dimuon_vChi2/F");
   outTree->Branch("dimuon_DCA", 	&out_dimuon_DCA, 	"dimuon_DCA/F");
   outTree->Branch("dimuon_ctauPV", 	&out_dimuon_ctauPV, 	"dimuon_ctauPV/F");
   outTree->Branch("dimuon_ctauErrPV", 	&out_dimuon_ctauErrPV, 	"dimuon_ctauErrPV/F");
   outTree->Branch("dimuon_cosAlpha", 	&out_dimuon_cosAlpha, 	"dimuon_cosAlpha/F");
   outTree->Branch("dimuon_triggerMatch", 	&out_dimuon_triggerMatch, 	"dimuon_triggerMatch/F");
   outTree->Branch("dimuonditrk_vProb", 	&out_dimuonditrk_vProb, 	"dimuonditrk_vProb/F");
   outTree->Branch("dimuonditrk_vChi2", 	&out_dimuonditrk_vChi2, 	"dimuonditrk_vChi2/F");
   outTree->Branch("dimuonditrk_nDof", 	&out_dimuonditrk_nDof, 	"dimuonditrk_nDof/F");
   outTree->Branch("dimuonditrk_charge", 	&out_dimuonditrk_charge, 	"dimuonditrk_charge/F");
   outTree->Branch("dimuonditrk_cosAlpha", 	&out_dimuonditrk_cosAlpha, 	"dimuonditrk_cosAlpha/F");
   outTree->Branch("dimuonditrk_ctauPV", 	&out_dimuonditrk_ctauPV, 	"dimuonditrk_ctauPV/F");
   outTree->Branch("dimuonditrk_ctauErrPV", 	&out_dimuonditrk_ctauErrPV, 	"dimuonditrk_ctauErrPV/F");
   outTree->Branch("dimuonditrk_tPFromPV", 	&out_dimuonditrk_tPFromPV, 	"dimuonditrk_tPFromPV/F");
   outTree->Branch("dimuonditrk_tMFromPV", 	&out_dimuonditrk_tMFromPV, 	"dimuonditrk_tMFromPV/F");
   outTree->Branch("dimuonditrk_cosAlphaCA", 	&out_dimuonditrk_cosAlphaCA, 	"dimuonditrk_cosAlphaCA/F");
   outTree->Branch("dimuonditrk_ctauPVCA", 	&out_dimuonditrk_ctauPVCA, 	"dimuonditrk_ctauPVCA/F");
   outTree->Branch("dimuonditrk_ctauErrPVCA", 	&out_dimuonditrk_ctauErrPVCA, 	"dimuonditrk_ctauErrPVCA/F");
   outTree->Branch("dimuonditrk_tPFromPVCA", 	&out_dimuonditrk_tPFromPVCA, 	"dimuonditrk_tPFromPVCA/F");
   outTree->Branch("dimuonditrk_tMFromPVCA", 	&out_dimuonditrk_tMFromPVCA, 	"dimuonditrk_tMFromPVCA/F");
   outTree->Branch("dimuonditrk_cosAlphaDZ", 	&out_dimuonditrk_cosAlphaDZ, 	"dimuonditrk_cosAlphaDZ/F");
   outTree->Branch("dimuonditrk_ctauPVDZ", 	&out_dimuonditrk_ctauPVDZ, 	"dimuonditrk_ctauPVDZ/F");
   outTree->Branch("dimuonditrk_ctauErrPVDZ", 	&out_dimuonditrk_ctauErrPVDZ, 	"dimuonditrk_ctauErrPVDZ/F");
   outTree->Branch("dimuonditrk_tPFromPVDZ", 	&out_dimuonditrk_tPFromPVDZ, 	"dimuonditrk_tPFromPVDZ/F");
   outTree->Branch("dimuonditrk_tMFromPVDZ", 	&out_dimuonditrk_tMFromPVDZ, 	"dimuonditrk_tMFromPVDZ/F");
   outTree->Branch("dimuonditrk_cosAlphaBS", 	&out_dimuonditrk_cosAlphaBS, 	"dimuonditrk_cosAlphaBS/F");
   outTree->Branch("dimuonditrk_ctauPVBS", 	&out_dimuonditrk_ctauPVBS, 	"dimuonditrk_ctauPVBS/F");
   outTree->Branch("dimuonditrk_ctauErrPVBS", 	&out_dimuonditrk_ctauErrPVBS, 	"dimuonditrk_ctauErrPVBS/F");
   outTree->Branch("dimuonditrk_dca_m1m2", 	&out_dimuonditrk_dca_m1m2, 	"dimuonditrk_dca_m1m2/F");
   outTree->Branch("dimuonditrk_dca_m1t1", 	&out_dimuonditrk_dca_m1t1, 	"dimuonditrk_dca_m1t1/F");
   outTree->Branch("dimuonditrk_dca_m1t2", 	&out_dimuonditrk_dca_m1t2, 	"dimuonditrk_dca_m1t2/F");
   outTree->Branch("dimuonditrk_dca_m2t1", 	&out_dimuonditrk_dca_m2t1, 	"dimuonditrk_dca_m2t1/F");
   outTree->Branch("dimuonditrk_dca_m2t2", 	&out_dimuonditrk_dca_m2t2, 	"dimuonditrk_dca_m2t2/F");
   outTree->Branch("dimuonditrk_dca_t1t2", 	&out_dimuonditrk_dca_t1t2, 	"dimuonditrk_dca_t1t2/F");
   outTree->Branch("dimuonditrk_rf_vProb", 	&out_dimuonditrk_rf_vProb, 	"dimuonditrk_rf_vProb/F");
   outTree->Branch("dimuonditrk_rf_vChi2", 	&out_dimuonditrk_rf_vChi2, 	"dimuonditrk_rf_vChi2/F");
   outTree->Branch("dimuonditrk_rf_nDof", 	&out_dimuonditrk_rf_nDof, 	"dimuonditrk_rf_nDof/F");
   outTree->Branch("dimuonditrk_rf_cosAlpha", 	&out_dimuonditrk_rf_cosAlpha, 	"dimuonditrk_rf_cosAlpha/F");
   outTree->Branch("dimuonditrk_rf_ctauPV", 	&out_dimuonditrk_rf_ctauPV, 	"dimuonditrk_rf_ctauPV/F");
   outTree->Branch("dimuonditrk_rf_ctauErrPV", 	&out_dimuonditrk_rf_ctauErrPV, 	"dimuonditrk_rf_ctauErrPV/F");
   outTree->Branch("dimuonditrk_rf_c_vProb", 	&out_dimuonditrk_rf_c_vProb, 	"dimuonditrk_rf_c_vProb/F");
   outTree->Branch("dimuonditrk_rf_c_vChi2", 	&out_dimuonditrk_rf_c_vChi2, 	"dimuonditrk_rf_c_vChi2/F");
   outTree->Branch("dimuonditrk_rf_c_nDof", 	&out_dimuonditrk_rf_c_nDof, 	"dimuonditrk_rf_c_nDof/F");
   outTree->Branch("dimuonditrk_rf_c_cosAlpha", 	&out_dimuonditrk_rf_c_cosAlpha, 	"dimuonditrk_rf_c_cosAlpha/F");
   outTree->Branch("dimuonditrk_rf_c_ctauPV", 	&out_dimuonditrk_rf_c_ctauPV, 	"dimuonditrk_rf_c_ctauPV/F");
   outTree->Branch("dimuonditrk_rf_c_ctauErrPV", 	&out_dimuonditrk_rf_c_ctauErrPV, 	"dimuonditrk_rf_c_ctauErrPV/F");
   outTree->Branch("highKaonMatch", 	&out_highKaonMatch, 	"highKaonMatch/F");
   outTree->Branch("lowKaonMatch", 	&out_lowKaonMatch, 	"lowKaonMatch/F");
   outTree->Branch("lowMuonMatch", 	&out_lowMuonMatch, 	"lowMuonMatch/F");
   outTree->Branch("highMuonMatch", 	&out_highMuonMatch, 	"highMuonMatch/F");
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
   outTree->Branch("highKaon_NPixelHits", 	&out_highKaon_NPixelHits, 	"highKaon_NPixelHits/F");
   outTree->Branch("highKaon_NStripHits", 	&out_highKaon_NStripHits, 	"highKaon_NStripHits/F");
   outTree->Branch("highKaon_NTrackhits", 	&out_highKaon_NTrackhits, 	"highKaon_NTrackhits/F");
   outTree->Branch("highKaon_NBPixHits", 	&out_highKaon_NBPixHits, 	"highKaon_NBPixHits/F");
   outTree->Branch("highKaon_NPixLayers", 	&out_highKaon_NPixLayers, 	"highKaon_NPixLayers/F");
   outTree->Branch("highKaon_NTraLayers", 	&out_highKaon_NTraLayers, 	"highKaon_NTraLayers/F");
   outTree->Branch("highKaon_NStrLayers", 	&out_highKaon_NStrLayers, 	"highKaon_NStrLayers/F");
   outTree->Branch("highKaon_NBPixLayers", 	&out_highKaon_NBPixLayers, 	"highKaon_NBPixLayers/F");
   outTree->Branch("lowKaon_NPixelHits", 	&out_lowKaon_NPixelHits, 	"lowKaon_NPixelHits/F");
   outTree->Branch("lowKaon_NStripHits", 	&out_lowKaon_NStripHits, 	"lowKaon_NStripHits/F");
   outTree->Branch("lowKaon_NTrackhits", 	&out_lowKaon_NTrackhits, 	"lowKaon_NTrackhits/F");
   outTree->Branch("lowKaon_NBPixHits", 	&out_lowKaon_NBPixHits, 	"lowKaon_NBPixHits/F");
   outTree->Branch("lowKaon_NPixLayers", 	&out_lowKaon_NPixLayers, 	"lowKaon_NPixLayers/F");
   outTree->Branch("lowKaon_NTraLayers", 	&out_lowKaon_NTraLayers, 	"lowKaon_NTraLayers/F");
   outTree->Branch("lowKaon_NStrLayers", 	&out_lowKaon_NStrLayers, 	"lowKaon_NStrLayers/F");
   outTree->Branch("lowKaon_NBPixLayers", 	&out_lowKaon_NBPixLayers, 	"lowKaon_NBPixLayers/F");

}

Bool_t FourTracks::Process(Long64_t entry)
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

   if(true)
   {

     out_run = 	(Float_t)(*run);
     out_event = 	(Float_t)(*event);
     out_lumi = 	(Float_t)(*lumi);
     out_numPrimaryVertices = 	(Float_t)(*numPrimaryVertices);
     out_trigger = 	(Float_t)(*trigger);
     out_noXCandidates = 	(Float_t)(*noXCandidates);
     out_dimuon_id = 	(Float_t)(*dimuon_id);
     out_p_id = 	(Float_t)(*p_id);
     out_m_id = 	(Float_t)(*m_id);
     // out_dimuonditrk_p4 = 	(Float_t)(*dimuonditrk_p4);
     // out_ditrack_p4 = 	(Float_t)(*ditrack_p4);
     // out_dimuon_p4 = 	(Float_t)(*dimuon_p4);
     // out_lowMuon_p4 = 	(Float_t)(*lowMuon_p4);
     // out_highMuon_p4 = 	(Float_t)(*highMuon_p4);
     // out_highKaon_p4 = 	(Float_t)(*highKaon_p4);
     // out_lowKaon_p4 = 	(Float_t)(*lowKaon_p4);
     // out_dimuonditrk_rf_p4 = 	(Float_t)(*dimuonditrk_rf_p4);
     // out_ditrack_rf_p4 = 	(Float_t)(*ditrack_rf_p4);
     // out_dimuon_rf_p4 = 	(Float_t)(*dimuon_rf_p4);
     // out_lowMuon_rf_p4 = 	(Float_t)(*lowMuon_rf_p4);
     // out_highMuon_rf_p4 = 	(Float_t)(*highMuon_rf_p4);
     // out_kaonp_rf_p4 = 	(Float_t)(*kaonp_rf_p4);
     // out_kaonn_rf_p4 = 	(Float_t)(*kaonn_rf_p4);
     out_dimuonditrk_m = 	(Float_t)(*dimuonditrk_m);
     out_dimuonditrk_m_rf = 	(Float_t)(*dimuonditrk_m_rf);
     out_dimuonditrk_m_rf_d_c = 	(Float_t)(*dimuonditrk_m_rf_d_c);
     out_dimuonditrk_m_rf_c = 	(Float_t)(*dimuonditrk_m_rf_c);
     out_dimuonditrk_pt = 	(Float_t)(*dimuonditrk_pt);
     out_dimuonditrk_eta = 	(Float_t)(*dimuonditrk_eta);
     out_dimuonditrk_phi = 	(Float_t)(*dimuonditrk_phi);
     out_dimuonditrk_y = 	(Float_t)(*dimuonditrk_y);
     out_dimuonditrk_vx = 	(Float_t)(*dimuonditrk_vx);
     out_dimuonditrk_vy = 	(Float_t)(*dimuonditrk_vy);
     out_dimuonditrk_vz = 	(Float_t)(*dimuonditrk_vz);
     out_pv_x = 	(Float_t)(*pv_x);
     out_pv_y = 	(Float_t)(*pv_y);
     out_pv_z = 	(Float_t)(*pv_z);
     out_dimuon_m = 	(Float_t)(*dimuon_m);
     out_dimuon_pt = 	(Float_t)(*dimuon_pt);
     out_ditrack_m = 	(Float_t)(*ditrack_m);
     out_ditrack_pt = 	(Float_t)(*ditrack_pt);
     out_highKaon_pt = 	(Float_t)(*highKaon_pt);
     out_lowKaon_pt = 	(Float_t)(*lowKaon_pt);
     out_highMuon_pt = 	(Float_t)(*highMuon_pt);
     out_lowMuon_pt = 	(Float_t)(*lowMuon_pt);
     out_highKaon_eta = 	(Float_t)(*highKaon_eta);
     out_lowKaon_eta = 	(Float_t)(*lowKaon_eta);
     out_highMuon_eta = 	(Float_t)(*highMuon_eta);
     out_lowMuon_eta = 	(Float_t)(*lowMuon_eta);
     out_highKaon_phi = 	(Float_t)(*highKaon_phi);
     out_lowKaon_phi = 	(Float_t)(*lowKaon_phi);
     out_highMuon_phi = 	(Float_t)(*highMuon_phi);
     out_lowMuon_phi = 	(Float_t)(*lowMuon_phi);
     out_highKaon_dxy = 	(Float_t)(*highKaon_dxy);
     out_lowKaon_dxy = 	(Float_t)(*lowKaon_dxy);
     out_highMuon_dxy = 	(Float_t)(*highMuon_dxy);
     out_lowMuon_dxy = 	(Float_t)(*lowMuon_dxy);
     out_highKaon_dz = 	(Float_t)(*highKaon_dz);
     out_lowKaon_dz = 	(Float_t)(*lowKaon_dz);
     out_highMuon_dz = 	(Float_t)(*highMuon_dz);
     out_lowMuon_dz = 	(Float_t)(*lowMuon_dz);
     out_highKaonMuonDR = 	(Float_t)(*highKaonMuonDR);
     out_highKaonMuonDP = 	(Float_t)(*highKaonMuonDP);
     out_highKaonMuonDPt = 	(Float_t)(*highKaonMuonDPt);
     out_lowKaonMuonDR = 	(Float_t)(*lowKaonMuonDR);
     out_lowKaonMuonDP = 	(Float_t)(*lowKaonMuonDP);
     out_lowKaonMuonDPt = 	(Float_t)(*lowKaonMuonDPt);
     out_dimuonditrk_refPK_mass = 	(Float_t)(*dimuonditrk_refPK_mass);
     out_dimuonditrk_refKP_mass = 	(Float_t)(*dimuonditrk_refKP_mass);
     out_dimuonditrk_refPP_mass = 	(Float_t)(*dimuonditrk_refPP_mass);
     out_dimuonditrk_refPK_vChi2 = 	(Float_t)(*dimuonditrk_refPK_vChi2);
     out_dimuonditrk_refKP_vChi2 = 	(Float_t)(*dimuonditrk_refKP_vChi2);
     out_dimuonditrk_refPP_vChi2 = 	(Float_t)(*dimuonditrk_refPP_vChi2);
     out_dimuonditrk_refPK_nDof = 	(Float_t)(*dimuonditrk_refPK_nDof);
     out_dimuonditrk_refKP_nDof = 	(Float_t)(*dimuonditrk_refKP_nDof);
     out_dimuonditrk_refPP_nDof = 	(Float_t)(*dimuonditrk_refPP_nDof);
     out_dimuonditrk_refPK_vProb = 	(Float_t)(*dimuonditrk_refPK_vProb);
     out_dimuonditrk_refKP_vProb = 	(Float_t)(*dimuonditrk_refKP_vProb);
     out_dimuonditrk_refPP_vProb = 	(Float_t)(*dimuonditrk_refPP_vProb);
     out_dimuon_vProb = 	(Float_t)(*dimuon_vProb);
     out_dimuon_vChi2 = 	(Float_t)(*dimuon_vChi2);
     out_dimuon_DCA = 	(Float_t)(*dimuon_DCA);
     out_dimuon_ctauPV = 	(Float_t)(*dimuon_ctauPV);
     out_dimuon_ctauErrPV = 	(Float_t)(*dimuon_ctauErrPV);
     out_dimuon_cosAlpha = 	(Float_t)(*dimuon_cosAlpha);
     out_dimuon_triggerMatch = 	(Float_t)(*dimuon_triggerMatch);
     out_dimuonditrk_vProb = 	(Float_t)(*dimuonditrk_vProb);
     out_dimuonditrk_vChi2 = 	(Float_t)(*dimuonditrk_vChi2);
     out_dimuonditrk_nDof = 	(Float_t)(*dimuonditrk_nDof);
     out_dimuonditrk_charge = 	(Float_t)(*dimuonditrk_charge);
     out_dimuonditrk_cosAlpha = 	(Float_t)(*dimuonditrk_cosAlpha);
     out_dimuonditrk_ctauPV = 	(Float_t)(*dimuonditrk_ctauPV);
     out_dimuonditrk_ctauErrPV = 	(Float_t)(*dimuonditrk_ctauErrPV);
     out_dimuonditrk_tPFromPV = 	(Float_t)(*dimuonditrk_tPFromPV);
     out_dimuonditrk_tMFromPV = 	(Float_t)(*dimuonditrk_tMFromPV);
     out_dimuonditrk_cosAlphaCA = 	(Float_t)(*dimuonditrk_cosAlphaCA);
     out_dimuonditrk_ctauPVCA = 	(Float_t)(*dimuonditrk_ctauPVCA);
     out_dimuonditrk_ctauErrPVCA = 	(Float_t)(*dimuonditrk_ctauErrPVCA);
     out_dimuonditrk_tPFromPVCA = 	(Float_t)(*dimuonditrk_tPFromPVCA);
     out_dimuonditrk_tMFromPVCA = 	(Float_t)(*dimuonditrk_tMFromPVCA);
     out_dimuonditrk_cosAlphaDZ = 	(Float_t)(*dimuonditrk_cosAlphaDZ);
     out_dimuonditrk_ctauPVDZ = 	(Float_t)(*dimuonditrk_ctauPVDZ);
     out_dimuonditrk_ctauErrPVDZ = 	(Float_t)(*dimuonditrk_ctauErrPVDZ);
     out_dimuonditrk_tPFromPVDZ = 	(Float_t)(*dimuonditrk_tPFromPVDZ);
     out_dimuonditrk_tMFromPVDZ = 	(Float_t)(*dimuonditrk_tMFromPVDZ);
     out_dimuonditrk_cosAlphaBS = 	(Float_t)(*dimuonditrk_cosAlphaBS);
     out_dimuonditrk_ctauPVBS = 	(Float_t)(*dimuonditrk_ctauPVBS);
     out_dimuonditrk_ctauErrPVBS = 	(Float_t)(*dimuonditrk_ctauErrPVBS);
     out_dimuonditrk_dca_m1m2 = 	(Float_t)(*dimuonditrk_dca_m1m2);
     out_dimuonditrk_dca_m1t1 = 	(Float_t)(*dimuonditrk_dca_m1t1);
     out_dimuonditrk_dca_m1t2 = 	(Float_t)(*dimuonditrk_dca_m1t2);
     out_dimuonditrk_dca_m2t1 = 	(Float_t)(*dimuonditrk_dca_m2t1);
     out_dimuonditrk_dca_m2t2 = 	(Float_t)(*dimuonditrk_dca_m2t2);
     out_dimuonditrk_dca_t1t2 = 	(Float_t)(*dimuonditrk_dca_t1t2);
     out_dimuonditrk_rf_vProb = 	(Float_t)(*dimuonditrk_rf_vProb);
     out_dimuonditrk_rf_vChi2 = 	(Float_t)(*dimuonditrk_rf_vChi2);
     out_dimuonditrk_rf_nDof = 	(Float_t)(*dimuonditrk_rf_nDof);
     out_dimuonditrk_rf_cosAlpha = 	(Float_t)(*dimuonditrk_rf_cosAlpha);
     out_dimuonditrk_rf_ctauPV = 	(Float_t)(*dimuonditrk_rf_ctauPV);
     out_dimuonditrk_rf_ctauErrPV = 	(Float_t)(*dimuonditrk_rf_ctauErrPV);
     out_dimuonditrk_rf_c_vProb = 	(Float_t)(*dimuonditrk_rf_c_vProb);
     out_dimuonditrk_rf_c_vChi2 = 	(Float_t)(*dimuonditrk_rf_c_vChi2);
     out_dimuonditrk_rf_c_nDof = 	(Float_t)(*dimuonditrk_rf_c_nDof);
     out_dimuonditrk_rf_c_cosAlpha = 	(Float_t)(*dimuonditrk_rf_c_cosAlpha);
     out_dimuonditrk_rf_c_ctauPV = 	(Float_t)(*dimuonditrk_rf_c_ctauPV);
     out_dimuonditrk_rf_c_ctauErrPV = 	(Float_t)(*dimuonditrk_rf_c_ctauErrPV);
     out_highKaonMatch = 	(Float_t)(*highKaonMatch);
     out_lowKaonMatch = 	(Float_t)(*lowKaonMatch);
     out_lowMuonMatch = 	(Float_t)(*lowMuonMatch);
     out_highMuonMatch = 	(Float_t)(*highMuonMatch);
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
     out_highKaon_NPixelHits = 	(Float_t)(*highKaon_NPixelHits);
     out_highKaon_NStripHits = 	(Float_t)(*highKaon_NStripHits);
     out_highKaon_NTrackhits = 	(Float_t)(*highKaon_NTrackhits);
     out_highKaon_NBPixHits = 	(Float_t)(*highKaon_NBPixHits);
     out_highKaon_NPixLayers = 	(Float_t)(*highKaon_NPixLayers);
     out_highKaon_NTraLayers = 	(Float_t)(*highKaon_NTraLayers);
     out_highKaon_NStrLayers = 	(Float_t)(*highKaon_NStrLayers);
     out_highKaon_NBPixLayers = 	(Float_t)(*highKaon_NBPixLayers);
     out_lowKaon_NPixelHits = 	(Float_t)(*lowKaon_NPixelHits);
     out_lowKaon_NStripHits = 	(Float_t)(*lowKaon_NStripHits);
     out_lowKaon_NTrackhits = 	(Float_t)(*lowKaon_NTrackhits);
     out_lowKaon_NBPixHits = 	(Float_t)(*lowKaon_NBPixHits);
     out_lowKaon_NPixLayers = 	(Float_t)(*lowKaon_NPixLayers);
     out_lowKaon_NTraLayers = 	(Float_t)(*lowKaon_NTraLayers);
     out_lowKaon_NStrLayers = 	(Float_t)(*lowKaon_NStrLayers);
     out_lowKaon_NBPixLayers = 	(Float_t)(*lowKaon_NBPixLayers);

     outTree->Fill();
   }

   return kTRUE;
}

void FourTracks::SlaveTerminate()
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

void FourTracks::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}
