#define TwoMuTwoK_cxx

// The class definition in TwoMuTwoK.h has been generated automatically
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
// root> T->Process("TwoMuTwoK.C")
// root> T->Process("TwoMuTwoK.C","some options")
// root> T->Process("TwoMuTwoK.C+")
//


#include "TwoMuTwoK.h"
#include <TH2.h>
#include <TStyle.h>


Double_t out;

//TTree* outTuple = new TTree("2mu2kSkimmedTree","2mu2kSkimmedTree");
//TBranch* b = outTuple->Branch("out",       &out,        "out/D");

void TwoMuTwoK::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
}

void TwoMuTwoK::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  std::string outputString = "2mu2k_five_tree.root";
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

  outTree->Branch("run",                &out_run,                "run/F");
  outTree->Branch("event",              &out_event,              "event/F");
  outTree->Branch("lumi",              &out_lumi,              "lumi/F");
  outTree->Branch("numPrimaryVertices", &out_numPrimaryVertices, "numPrimaryVertices/F");
  outTree->Branch("trigger",            &out_trigger,            "trigger/F");

  outTree->Branch("noXCandidates",      &out_noXCandidates,      "noXCandidates/F");

  //kin
  outTree->Branch("dimuonditrk_m",       &out_dimuonditrk_m,        "dimuonditrk_m/F");
  outTree->Branch("dimuonditrk_m_rf",       &out_dimuonditrk_m_rf,        "dimuonditrk_m_rf/F");
  outTree->Branch("dimuonditrk_m_rf_d_c",       &out_dimuonditrk_m_rf_d_c,        "dimuonditrk_m_rf_d_c/F");
  outTree->Branch("dimuonditrk_m_rf_c",       &out_dimuonditrk_m_rf_c,        "dimuonditrk_m_rf_c/F");
  outTree->Branch("dimuonditrk_pt",          &out_dimuonditrk_pt,          "dimuonditrk_pt/F");
  outTree->Branch("dimuonditrk_eta",          &out_dimuonditrk_eta,          "dimuonditrk_eta/F");
  outTree->Branch("dimuonditrk_phi",          &out_dimuonditrk_phi,          "dimuonditrk_phi/F");
  outTree->Branch("dimuonditrk_y",          &out_dimuonditrk_y,          "dimuonditrk_y/F");

  outTree->Branch("dimuonditrk_vx",          &out_dimuonditrk_vx,          "dimuonditrk_vx/F");
  outTree->Branch("dimuonditrk_vy",          &out_dimuonditrk_vy,          "dimuonditrk_vy/F");
  outTree->Branch("dimuonditrk_vz",          &out_dimuonditrk_vz,          "dimuonditrk_vz/F");

  outTree->Branch("pv_x",          &out_pv_x,          "pv_x/F");
  outTree->Branch("pv_y",          &out_pv_y,          "pv_y/F");
  outTree->Branch("pv_z",          &out_pv_z,          "pv_z/F");

  outTree->Branch("dimuon_m",       &out_dimuon_m,       "dimuon_m/F");
  outTree->Branch("dimuon_pt",    &out_dimuon_pt,    "dimuon_pt/F");
  outTree->Branch("ditrak_m",     &out_ditrak_m,     "ditrak_m/F");
  outTree->Branch("ditrak_pt",       &out_ditrak_pt,        "ditrak_pt/F");

  outTree->Branch("highKaon_pt",          &out_highKaon_pt,          "highKaon_pt/F");
  outTree->Branch("lowKaon_pt",       &out_lowKaon_pt,       "lowKaon_pt/F");
  outTree->Branch("highMuon_pt",    &out_highMuon_pt,    "highMuon_pt/F");
  outTree->Branch("lowMuon_pt",     &out_lowMuon_pt,     "lowMuon_pt/F");

  outTree->Branch("highKaon_eta",        &out_highKaon_eta,        "highKaon_eta/F");
  outTree->Branch("lowKaon_eta",        &out_lowKaon_eta,        "lowKaon_eta/F");
  outTree->Branch("highMuon_eta",        &out_highMuon_eta,        "highMuon_eta/F");
  outTree->Branch("lowMuon_eta",        &out_lowMuon_eta,        "lowMuon_eta/F");

  outTree->Branch("highKaon_phi",        &out_highKaon_phi,        "highKaon_phi/F");
  outTree->Branch("lowKaon_phi",        &out_lowKaon_phi,        "lowKaon_phi/F");
  outTree->Branch("highMuon_phi",        &out_highMuon_phi,        "highMuon_phi/F");
  outTree->Branch("lowMuon_phi",        &out_lowMuon_phi,        "lowMuon_phi/F");

  outTree->Branch("highKaon_dxy",        &out_highKaon_dxy,        "highKaon_dxy/F");
  outTree->Branch("lowKaon_dxy",        &out_lowKaon_dxy,        "lowKaon_dxy/F");
  outTree->Branch("highMuon_dxy",        &out_highMuon_dxy,        "highMuon_dxy/F");
  outTree->Branch("lowMuon_dxy",        &out_lowMuon_dxy,        "lowMuon_dxy/F");

  outTree->Branch("highKaon_dz",        &out_highKaon_dz,        "highKaon_dz/F");
  outTree->Branch("lowKaon_dz",        &out_lowKaon_dz,        "lowKaon_dz/F");
  outTree->Branch("highMuon_dz",        &out_highMuon_dz,        "highMuon_dz/F");
  outTree->Branch("lowMuon_dz",        &out_lowMuon_dz,        "lowMuon_dz/F");


  //Pion refits
  outTree->Branch("dimuonditrk_refPK_mass", &out_dimuonditrk_refPK_mass, "dimuonditrk_refPK_mass/F");
  outTree->Branch("dimuonditrk_refKP_mass", &out_dimuonditrk_refKP_mass, "dimuonditrk_refKP_mass/F");
  outTree->Branch("dimuonditrk_refPP_mass", &out_dimuonditrk_refPP_mass, "dimuonditrk_refPP_mass/F");
  outTree->Branch("dimuonditrk_refPK_vChi2", &out_dimuonditrk_refPK_vChi2, "dimuonditrk_refPK_vChi2/F");
  outTree->Branch("dimuonditrk_refKP_vChi2", &out_dimuonditrk_refKP_vChi2, "dimuonditrk_refKP_vChi2/F");
  outTree->Branch("dimuonditrk_refPP_vChi2", &out_dimuonditrk_refPP_vChi2, "dimuonditrk_refPP_vChi2/F");
  outTree->Branch("dimuonditrk_refPK_nDof", &out_dimuonditrk_refPK_nDof, "dimuonditrk_refPK_nDof/F");
  outTree->Branch("dimuonditrk_refKP_nDof", &out_dimuonditrk_refKP_nDof, "dimuonditrk_refKP_nDof/F");
  outTree->Branch("dimuonditrk_refPP_nDof", &out_dimuonditrk_refPP_nDof, "dimuonditrk_refPP_nDof/F");
  outTree->Branch("dimuonditrk_refPK_vProb", &out_dimuonditrk_refPK_vProb, "dimuonditrk_refPK_vProb/F");
  outTree->Branch("dimuonditrk_refKP_vProb", &out_dimuonditrk_refKP_vProb, "dimuonditrk_refKP_vProb/F");
  outTree->Branch("dimuonditrk_refPP_vProb", &out_dimuonditrk_refPP_vProb, "dimuonditrk_refPP_vProb/F");

  //2mu vertexing
  outTree->Branch("dimuon_vProb",        &out_dimuon_vProb,        "dimuon_vProb/F");
  outTree->Branch("dimuon_vChi2",       &out_dimuon_vChi2,        "dimuon_vChi2/F");
  outTree->Branch("dimuon_DCA",          &out_dimuon_DCA,          "dimuon_DCA/F");
  outTree->Branch("dimuon_ctauPV",       &out_dimuon_ctauPV,       "dimuon_ctauPV/F");
  outTree->Branch("dimuon_ctauErrPV",    &out_dimuon_ctauErrPV,    "dimuon_ctauErrPV/F");
  outTree->Branch("dimuon_cosAlpha",     &out_dimuon_cosAlpha,     "dimuon_cosAlpha/F");
  outTree->Branch("dimuon_triggerMatch", &out_dimuon_triggerMatch, "dimuon_triggerMatch/F");

  //2mu+2Trk vertexing
  outTree->Branch("dimuonditrk_vProb",      &out_dimuonditrk_vProb,        "dimuonditrk_vProb/F");
  outTree->Branch("dimuonditrk_vChi2",      &out_dimuonditrk_vChi2,        "dimuonditrk_vChi2/F");
  outTree->Branch("dimuonditrk_nDof",       &out_dimuonditrk_nDof,         "dimuonditrk_nDof/F");
  outTree->Branch("dimuonditrk_charge",     &out_dimuonditrk_charge,       "dimuonditrk_charge/F");

  outTree->Branch("dimuonditrk_cosAlpha",      &out_dimuonditrk_cosAlpha,        "dimuonditrk_cosAlpha/F");
  outTree->Branch("dimuonditrk_ctauPV",      &out_dimuonditrk_ctauPV,        "dimuonditrk_ctauPV/F");
  outTree->Branch("dimuonditrk_ctauErrPV",      &out_dimuonditrk_ctauErrPV,        "dimuonditrk_ctauErrPV/F");


  outTree->Branch("dimuonditrk_tPFromPV",      &out_dimuonditrk_tPFromPV,        "dimuonditrk_tPFromPV/F");
  outTree->Branch("dimuonditrk_tMFromPV",      &out_dimuonditrk_tMFromPV,        "dimuonditrk_tMFromPV/F");


  outTree->Branch("dimuonditrk_cosAlphaDZ",      &out_dimuonditrk_cosAlphaDZ,        "dimuonditrk_cosAlphaDZ/F");
  outTree->Branch("dimuonditrk_ctauPVDZ",      &out_dimuonditrk_ctauPVDZ,        "dimuonditrk_ctauPVDZ/F");
  outTree->Branch("dimuonditrk_ctauErrPVDZ",      &out_dimuonditrk_ctauErrPVDZ,        "dimuonditrk_ctauErrPVDZ/F");

  outTree->Branch("dimuonditrk_tPFromPVDZ",      &out_dimuonditrk_tPFromPVDZ,        "dimuonditrk_tPFromPVDZ/F");
  outTree->Branch("dimuonditrk_tMFromPVDZ",      &out_dimuonditrk_tMFromPVDZ,        "dimuonditrk_tMFromPVDZ/F");


  outTree->Branch("dimuonditrk_cosAlphaBS",      &out_dimuonditrk_cosAlphaBS,        "dimuonditrk_cosAlphaBS/F");
  outTree->Branch("dimuonditrk_ctauPVBS",      &out_dimuonditrk_ctauPVBS,        "dimuonditrk_ctauPVBS/F");
  outTree->Branch("dimuonditrk_ctauErrPVBS",      &out_dimuonditrk_ctauErrPVBS,        "dimuonditrk_ctauErrPVBS/F");

  outTree->Branch("dimuonditrk_tPFromPVBS",      &out_dimuonditrk_tPFromPVBS,        "dimuonditrk_tPFromPVBS/F");
  outTree->Branch("dimuonditrk_tMFromPVBS",      &out_dimuonditrk_tMFromPVBS,        "dimuonditrk_tMFromPVBS/F");

  outTree->Branch("dimuonditrk_dca_m1m2",      &out_dimuonditrk_vProb,        "dimuonditrk_dca_m1m2/F");
  outTree->Branch("dimuonditrk_dca_m1t1",      &out_dimuonditrk_vProb,        "dimuonditrk_dca_m1t1/F");
  outTree->Branch("dimuonditrk_dca_m1t2",      &out_dimuonditrk_vProb,        "dimuonditrk_dca_m1t2/F");
  outTree->Branch("dimuonditrk_dca_m2t1",      &out_dimuonditrk_vProb,        "dimuonditrk_dca_m2t1/F");
  outTree->Branch("dimuonditrk_dca_m2t2",      &out_dimuonditrk_vProb,        "dimuonditrk_dca_m2t2/F");
  outTree->Branch("dimuonditrk_dca_t1t2",      &out_dimuonditrk_vProb,        "dimuonditrk_dca_t1t2/F");

  outTree->Branch("dimuonditrk_rf_vProb",      &out_dimuonditrk_rf_vProb,        "dimuonditrk_rf_vProb/F");
  outTree->Branch("dimuonditrk_rf_vChi2",      &out_dimuonditrk_rf_vChi2,        "dimuonditrk_rf_vChi2/F");
  outTree->Branch("dimuonditrk_rf_nDof",       &out_dimuonditrk_rf_nDof,         "dimuonditrk_rf_nDof/F");
  outTree->Branch("dimuonditrk_rf_cosAlpha",   &out_dimuonditrk_rf_cosAlpha,     "dimuonditrk_rf_cosAlpha/F");
  outTree->Branch("dimuonditrk_rf_ctauPV",     &out_dimuonditrk_rf_ctauPV,       "dimuonditrk_rf_ctauPV/F");
  outTree->Branch("dimuonditrk_rf_ctauErrPV",  &out_dimuonditrk_rf_ctauErrPV,    "dimuonditrk_rf_ctauErrPV/F");

  outTree->Branch("dimuonditrk_rf_c_vProb",      &out_dimuonditrk_rf_c_vProb,        "dimuonditrk_rf_c_vProb/F");
  outTree->Branch("dimuonditrk_rf_c_vChi2",      &out_dimuonditrk_rf_c_vChi2,        "dimuonditrk_rf_c_vChi2/F");
  outTree->Branch("dimuonditrk_rf_c_nDof",       &out_dimuonditrk_rf_c_nDof,         "dimuonditrk_rf_c_nDof/F");
  outTree->Branch("dimuonditrk_rf_c_cosAlpha",   &out_dimuonditrk_rf_c_cosAlpha,     "dimuonditrk_rf_c_cosAlpha/F");
  outTree->Branch("dimuonditrk_rf_c_ctauPV",     &out_dimuonditrk_rf_c_ctauPV,       "dimuonditrk_rf_c_ctauPV/F");
  outTree->Branch("dimuonditrk_rf_c_ctauErrPV",  &out_dimuonditrk_rf_c_ctauErrPV,    "dimuonditrk_rf_c_ctauErrPV/F");

  outTree->Branch("highKaonMatch",     &out_highKaonMatch,       "highKaonMatch/F");
  outTree->Branch("lowKaonMatch",     &out_lowKaonMatch,       "lowKaonMatch/F");
  outTree->Branch("lowMuonMatch",     &out_lowMuonMatch,       "lowMuonMatch/F");
  outTree->Branch("highMuonMatch",     &out_highMuonMatch,       "highMuonMatch/F");

  //Muon flags
  outTree->Branch("lowMuon_isTight",        &out_lowMuon_isTight,        "lowMuon_isTight/F");
  outTree->Branch("lowMuon_isLoose",        &out_lowMuon_isLoose,        "lowMuon_isLoose/F");
  outTree->Branch("lowMuon_isSoft",        &out_lowMuon_isSoft,        "lowMuon_isSoft/F");
  outTree->Branch("lowMuon_isMedium",        &out_lowMuon_isMedium,        "lowMuon_isMedium/F");
  outTree->Branch("lowMuon_isHighPt",        &out_lowMuon_isHighPt,        "lowMuon_isHighPt/F");

  outTree->Branch("lowMuon_isTracker",        &out_lowMuon_isTracker,        "lowMuon_isTracker/F");
  outTree->Branch("lowMuon_isGlobal",        &out_lowMuon_isGlobal,        "lowMuon_isGlobal/F");

  outTree->Branch("lowMuon_NPixelHits",        &out_lowMuon_NPixelHits,        "lowMuon_NPixelHits/F");
  outTree->Branch("lowMuon_NStripHits",        &out_lowMuon_NStripHits,        "lowMuon_NStripHits/F");
  outTree->Branch("lowMuon_NTrackhits",        &out_lowMuon_NTrackhits,        "lowMuon_NTrackhits/F");
  outTree->Branch("lowMuon_NBPixHits",        &out_lowMuon_NBPixHits,        "lowMuon_NBPixHits/F");

  outTree->Branch("lowMuon_NPixLayers",        &out_lowMuon_NPixLayers,        "lowMuon_NPixLayers/F");
  outTree->Branch("lowMuon_NTraLayers",        &out_lowMuon_NTraLayers,        "lowMuon_NTraLayers/F");
  outTree->Branch("lowMuon_NStrLayers",        &out_lowMuon_NStrLayers,        "lowMuon_NStrLayers/F");
  outTree->Branch("lowMuon_NBPixLayers",        &out_lowMuon_NBPixLayers,        "lowMuon_NBPixLayers/F");

  outTree->Branch("highMuon_isTight",        &out_highMuon_isTight,        "highMuon_isTight/F");
  outTree->Branch("highMuon_isLoose",        &out_highMuon_isLoose,        "highMuon_isLoose/F");
  outTree->Branch("highMuon_isSoft",        &out_highMuon_isSoft,        "highMuon_isSoft/F");
  outTree->Branch("highMuon_isMedium",        &out_highMuon_isMedium,        "highMuon_isMedium/F");
  outTree->Branch("highMuon_isHighPt",        &out_highMuon_isHighPt,        "highMuon_isHighPt/F");

  outTree->Branch("highMuon_isTracker",        &out_highMuon_isTracker,        "highMuon_isTracker/F");
  outTree->Branch("highMuon_isGlobal",        &out_highMuon_isGlobal,        "highMuon_isGlobal/F");

  outTree->Branch("highMuon_NPixelHits",        &out_highMuon_NPixelHits,        "highMuon_NPixelHits/F");
  outTree->Branch("highMuon_NStripHits",        &out_highMuon_NStripHits,        "highMuon_NStripHits/F");
  outTree->Branch("highMuon_NTrackhits",        &out_highMuon_NTrackhits,        "highMuon_NTrackhits/F");
  outTree->Branch("highMuon_NBPixHits",        &out_highMuon_NBPixHits,        "highMuon_NBPixHits/F");

  outTree->Branch("highMuon_NPixLayers",        &out_highMuon_NPixLayers,        "highMuon_NPixLayers/F");
  outTree->Branch("highMuon_NTraLayers",        &out_highMuon_NTraLayers,        "highMuon_NTraLayers/F");
  outTree->Branch("highMuon_NStrLayers",        &out_highMuon_NStrLayers,        "highMuon_NStrLayers/F");
  outTree->Branch("highMuon_NBPixLayers",        &out_highMuon_NBPixLayers,        "highMuon_NBPixLayers/F");

  outTree->Branch("lowMuon_type",     &out_lowMuon_type,       "lowMuon_type/F");
  outTree->Branch("highMuon_type",     &out_highMuon_type,       "highMuon_type/F");

  //Tracks Flags

  outTree->Branch("highKaon_NPixelHits",        &out_highKaon_NPixelHits,        "highKaon_NPixelHits/F");
  outTree->Branch("highKaon_NStripHits",        &out_highKaon_NStripHits,        "highKaon_NStripHits/F");
  outTree->Branch("highKaon_NTrackhits",        &out_highKaon_NTrackhits,        "highKaon_NTrackhits/F");
  outTree->Branch("highKaon_NBPixHits",        &out_highKaon_NBPixHits,        "highKaon_NBPixHits/F");

  outTree->Branch("highKaon_NPixLayers",        &out_highKaon_NPixLayers,        "highKaon_NPixLayers/F");
  outTree->Branch("highKaon_NTraLayers",        &out_highKaon_NTraLayers,        "highKaon_NTraLayers/F");
  outTree->Branch("highKaon_NStrLayers",        &out_highKaon_NStrLayers,        "highKaon_NStrLayers/F");
  outTree->Branch("highKaon_NBPixLayers",        &out_highKaon_NBPixLayers,        "highKaon_NBPixLayers/F");

  outTree->Branch("lowKaon_NPixelHits",        &out_lowKaon_NPixelHits,        "lowKaon_NPixelHits/F");
  outTree->Branch("lowKaon_NStripHits",        &out_lowKaon_NStripHits,        "lowKaon_NStripHits/F");
  outTree->Branch("lowKaon_NTrackhits",        &out_lowKaon_NTrackhits,        "lowKaon_NTrackhits/F");
  outTree->Branch("lowKaon_NBPixHits",        &out_lowKaon_NBPixHits,        "lowKaon_NBPixHits/F");

  outTree->Branch("lowKaon_NPixLayers",        &out_lowKaon_NPixLayers,        "lowKaon_NPixLayers/F");
  outTree->Branch("lowKaon_NTraLayers",        &out_lowKaon_NTraLayers,        "lowKaon_NTraLayers/F");
  outTree->Branch("lowKaon_NStrLayers",        &out_lowKaon_NStrLayers,        "lowKaon_NStrLayers/F");
  outTree->Branch("lowKaon_NBPixLayers",        &out_lowKaon_NBPixLayers,        "lowKaon_NBPixLayers/F");




}

Bool_t TwoMuTwoK::Process(Long64_t entry)
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

  bool test = (*out_run) < 320000;
  //int a = (int) (*trigger);
  //std::cout << (*trigger);
  //

  if(test)
  {
    out_run  = (*run);
    out_event        = (*event);
    out_lumi         = (*lumi);
    out_numPrimaryVertices   = (*numPrimaryVertices);
    out_trigger      = (*trigger);
    out_noXCandidates        = (*noXCandidates);

    out_dimuonditrk_m        = (*dimuonditrk_m);
    out_dimuonditrk_m_rf     = (*dimuonditrk_m_rf);
    out_dimuonditrk_m_rf_d_c         = (*dimuonditrk_m_rf_d_c);
    out_dimuonditrk_m_rf_c   = (*dimuonditrk_m_rf_c);
    out_dimuonditrk_pt       = (*dimuonditrk_pt);
    out_dimuonditrk_eta      = (*dimuonditrk_eta);
    out_dimuonditrk_phi      = (*dimuonditrk_phi);
    out_dimuonditrk_y        = (*dimuonditrk_y);
    out_dimuonditrk_vx       = (*dimuonditrk_vx);
    out_dimuonditrk_vy       = (*dimuonditrk_vy);
    out_dimuonditrk_vz       = (*dimuonditrk_vz);
    out_pv_x         = (*pv_x);
    out_pv_y         = (*pv_y);
    out_pv_z         = (*pv_z);
    out_dimuon_m     = (*dimuon_m);
    out_dimuon_pt    = (*dimuon_pt);
    out_ditrak_m     = (*ditrak_m);
    out_ditrak_pt    = (*ditrak_pt);
    out_highKaon_pt  = (*highKaon_pt);
    out_lowKaon_pt   = (*lowKaon_pt);
    out_highMuon_pt  = (*highMuon_pt);
    out_lowMuon_pt   = (*lowMuon_pt);
    out_highKaon_eta         = (*highKaon_eta);
    out_lowKaon_eta  = (*lowKaon_eta);
    out_highMuon_eta         = (*highMuon_eta);
    out_lowMuon_eta  = (*lowMuon_eta);
    out_highKaon_phi         = (*highKaon_phi);
    out_lowKaon_phi  = (*lowKaon_phi);
    out_highMuon_phi         = (*highMuon_phi);
    out_lowMuon_phi  = (*lowMuon_phi);
    out_highKaon_dxy         = (*highKaon_dxy);
    out_lowKaon_dxy  = (*lowKaon_dxy);
    out_highMuon_dxy         = (*highMuon_dxy);
    out_lowMuon_dxy  = (*lowMuon_dxy);
    out_highKaon_dz  = (*highKaon_dz);
    out_lowKaon_dz   = (*lowKaon_dz);
    out_highMuon_dz  = (*highMuon_dz);
    out_lowMuon_dz   = (*lowMuon_dz);
    out_dimuonditrk_refPK_mass       = (*dimuonditrk_refPK_mass);
    out_dimuonditrk_refKP_mass       = (*dimuonditrk_refKP_mass);
    out_dimuonditrk_refPP_mass       = (*dimuonditrk_refPP_mass);
    out_dimuonditrk_refPK_vChi2      = (*dimuonditrk_refPK_vChi2);
    out_dimuonditrk_refKP_vChi2      = (*dimuonditrk_refKP_vChi2);
    out_dimuonditrk_refPP_vChi2      = (*dimuonditrk_refPP_vChi2);
    out_dimuonditrk_refPK_nDof       = (*dimuonditrk_refPK_nDof);
    out_dimuonditrk_refKP_nDof       = (*dimuonditrk_refKP_nDof);
    out_dimuonditrk_refPP_nDof       = (*dimuonditrk_refPP_nDof);
    out_dimuonditrk_refPK_vProb      = (*dimuonditrk_refPK_vProb);
    out_dimuonditrk_refKP_vProb      = (*dimuonditrk_refKP_vProb);
    out_dimuonditrk_refPP_vProb      = (*dimuonditrk_refPP_vProb);
    out_dimuon_vProb         = (*dimuon_vProb);
    out_dimuon_vChi2         = (*dimuon_vChi2);
    out_dimuon_DCA   = (*dimuon_DCA);
    out_dimuon_ctauPV        = (*dimuon_ctauPV);
    out_dimuon_ctauErrPV     = (*dimuon_ctauErrPV);
    out_dimuon_cosAlpha      = (*dimuon_cosAlpha);
    out_dimuon_triggerMatch  = (*dimuon_triggerMatch);
    out_dimuonditrk_vProb    = (*dimuonditrk_vProb);
    out_dimuonditrk_vChi2    = (*dimuonditrk_vChi2);
    out_dimuonditrk_nDof     = (*dimuonditrk_nDof);
    out_dimuonditrk_charge   = (*dimuonditrk_charge);
    out_dimuonditrk_cosAlpha         = (*dimuonditrk_cosAlpha);
    out_dimuonditrk_ctauPV   = (*dimuonditrk_ctauPV);
    out_dimuonditrk_ctauErrPV        = (*dimuonditrk_ctauErrPV);
    out_dimuonditrk_tPFromPV         = (*dimuonditrk_tPFromPV);
    out_dimuonditrk_tMFromPV         = (*dimuonditrk_tMFromPV);
    out_dimuonditrk_cosAlphaDZ       = (*dimuonditrk_cosAlphaDZ);
    out_dimuonditrk_ctauPVDZ         = (*dimuonditrk_ctauPVDZ);
    out_dimuonditrk_ctauErrPVDZ      = (*dimuonditrk_ctauErrPVDZ);
    out_dimuonditrk_tPFromPVDZ       = (*dimuonditrk_tPFromPVDZ);
    out_dimuonditrk_tMFromPVDZ       = (*dimuonditrk_tMFromPVDZ);
    out_dimuonditrk_cosAlphaBS       = (*dimuonditrk_cosAlphaBS);
    out_dimuonditrk_ctauPVBS         = (*dimuonditrk_ctauPVBS);
    out_dimuonditrk_ctauErrPVBS      = (*dimuonditrk_ctauErrPVBS);
    out_dimuonditrk_tPFromPVBS       = (*dimuonditrk_tPFromPVBS);
    out_dimuonditrk_tMFromPVBS       = (*dimuonditrk_tMFromPVBS);
    out_dimuonditrk_dca_m1m2         = (*dimuonditrk_dca_m1m2);
    out_dimuonditrk_dca_m1t1         = (*dimuonditrk_dca_m1t1);
    out_dimuonditrk_dca_m1t2         = (*dimuonditrk_dca_m1t2);
    out_dimuonditrk_dca_m2t1         = (*dimuonditrk_dca_m2t1);
    out_dimuonditrk_dca_m2t2         = (*dimuonditrk_dca_m2t2);
    out_dimuonditrk_dca_t1t2         = (*dimuonditrk_dca_t1t2);
    out_dimuonditrk_rf_vProb         = (*dimuonditrk_rf_vProb);
    out_dimuonditrk_rf_vChi2         = (*dimuonditrk_rf_vChi2);
    out_dimuonditrk_rf_nDof  = (*dimuonditrk_rf_nDof);
    out_dimuonditrk_rf_cosAlpha      = (*dimuonditrk_rf_cosAlpha);
    out_dimuonditrk_rf_ctauPV        = (*dimuonditrk_rf_ctauPV);
    out_dimuonditrk_rf_ctauErrPV     = (*dimuonditrk_rf_ctauErrPV);
    out_dimuonditrk_rf_c_vProb       = (*dimuonditrk_rf_c_vProb);
    out_dimuonditrk_rf_c_vChi2       = (*dimuonditrk_rf_c_vChi2);
    out_dimuonditrk_rf_c_nDof        = (*dimuonditrk_rf_c_nDof);
    out_dimuonditrk_rf_c_cosAlpha    = (*dimuonditrk_rf_c_cosAlpha);
    out_dimuonditrk_rf_c_ctauPV      = (*dimuonditrk_rf_c_ctauPV);
    out_dimuonditrk_rf_c_ctauErrPV   = (*dimuonditrk_rf_c_ctauErrPV);
    out_highKaonMatch        = (*highKaonMatch);
    out_lowKaonMatch         = (*lowKaonMatch);
    out_lowMuonMatch         = (*lowMuonMatch);
    out_highMuonMatch        = (*highMuonMatch);
    out_lowMuon_isTight      = (*lowMuon_isTight);
    out_lowMuon_isLoose      = (*lowMuon_isLoose);
    out_lowMuon_isSoft       = (*lowMuon_isSoft);
    out_lowMuon_isMedium     = (*lowMuon_isMedium);
    out_lowMuon_isHighPt     = (*lowMuon_isHighPt);
    out_lowMuon_isTracker    = (*lowMuon_isTracker);
    out_lowMuon_isGlobal     = (*lowMuon_isGlobal);
    out_lowMuon_NPixelHits   = (*lowMuon_NPixelHits);
    out_lowMuon_NStripHits   = (*lowMuon_NStripHits);
    out_lowMuon_NTrackhits   = (*lowMuon_NTrackhits);
    out_lowMuon_NBPixHits    = (*lowMuon_NBPixHits);
    out_lowMuon_NPixLayers   = (*lowMuon_NPixLayers);
    out_lowMuon_NTraLayers   = (*lowMuon_NTraLayers);
    out_lowMuon_NStrLayers   = (*lowMuon_NStrLayers);
    out_lowMuon_NBPixLayers  = (*lowMuon_NBPixLayers);
    out_highMuon_isTight     = (*highMuon_isTight);
    out_highMuon_isLoose     = (*highMuon_isLoose);
    out_highMuon_isSoft      = (*highMuon_isSoft);
    out_highMuon_isMedium    = (*highMuon_isMedium);
    out_highMuon_isHighPt    = (*highMuon_isHighPt);
    out_highMuon_isTracker   = (*highMuon_isTracker);
    out_highMuon_isGlobal    = (*highMuon_isGlobal);
    out_highMuon_NPixelHits  = (*highMuon_NPixelHits);
    out_highMuon_NStripHits  = (*highMuon_NStripHits);
    out_highMuon_NTrackhits  = (*highMuon_NTrackhits);
    out_highMuon_NBPixHits   = (*highMuon_NBPixHits);
    out_highMuon_NPixLayers  = (*highMuon_NPixLayers);
    out_highMuon_NTraLayers  = (*highMuon_NTraLayers);
    out_highMuon_NStrLayers  = (*highMuon_NStrLayers);
    out_highMuon_NBPixLayers         = (*highMuon_NBPixLayers);
    out_lowMuon_type         = (*lowMuon_type);
    out_highMuon_type        = (*highMuon_type);
    out_highKaon_NPixelHits  = (*highKaon_NPixelHits);
    out_highKaon_NStripHits  = (*highKaon_NStripHits);
    out_highKaon_NTrackhits  = (*highKaon_NTrackhits);
    out_highKaon_NBPixHits   = (*highKaon_NBPixHits);
    out_highKaon_NPixLayers  = (*highKaon_NPixLayers);
    out_highKaon_NTraLayers  = (*highKaon_NTraLayers);
    out_highKaon_NStrLayers  = (*highKaon_NStrLayers);
    out_highKaon_NBPixLayers         = (*highKaon_NBPixLayers);
    out_lowKaon_NPixelHits   = (*lowKaon_NPixelHits);
    out_lowKaon_NStripHits   = (*lowKaon_NStripHits);
    out_lowKaon_NTrackhits   = (*lowKaon_NTrackhits);
    out_lowKaon_NBPixHits    = (*lowKaon_NBPixHits);
    out_lowKaon_NPixLayers   = (*lowKaon_NPixLayers);
    out_lowKaon_NTraLayers   = (*lowKaon_NTraLayers);
    out_lowKaon_NStrLayers   = (*lowKaon_NStrLayers);
    out_lowKaon_NBPixLayers  = (*lowKaon_NBPixLayers);

    out_isBestCandidate      = (*isBestCandidate);
    

  }

  return kTRUE;
}

void TwoMuTwoK::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

}

void TwoMuTwoK::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

}
