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

  std::string outputString = "2mu2k_sPlot_tree.root";
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
  outTree->Branch("dimuonditrk_m_rf_c",       &out_dimuonditrk_m_rf_c,        "dimuonditrk_m_rf_c/F");

  outTree->Branch("dimuonditrk_pt",          &out_dimuonditrk_pt,          "dimuonditrk_pt/F");
  outTree->Branch("dimuonditrk_eta",          &out_dimuonditrk_eta,          "dimuonditrk_eta/F");
  outTree->Branch("dimuonditrk_phi",          &out_dimuonditrk_phi,          "dimuonditrk_phi/F");
  outTree->Branch("dimuonditrk_y",          &out_dimuonditrk_y,          "dimuonditrk_y/F");

  outTree->Branch("dimuon_m",       &out_dimuon_m,       "dimuon_m/F");
  outTree->Branch("dimuon_pt",    &out_dimuon_pt,    "dimuon_pt/F");
  outTree->Branch("ditrak_m",     &out_ditrak_m,     "ditrak_m/F");
  outTree->Branch("ditrak_pt",       &out_ditrak_pt,        "ditrak_pt/F");

  outTree->Branch("highKaon_pt",          &out_highKaon_pt,          "highKaon_pt/F");
  outTree->Branch("lowKaon_pt",       &out_lowKaon_pt,       "lowKaon_pt/F");
  outTree->Branch("highMuon_pt",    &out_highMuon_pt,    "highMuon_pt/F");
  outTree->Branch("lowMuon_pt",     &out_lowMuon_pt,     "lowMuon_pt/F");


  //2mu vertexing
  outTree->Branch("dimuon_vProb",        &out_dimuon_vProb,        "dimuon_vProb/F");
  outTree->Branch("dimuon_vChi2",       &out_dimuon_vChi2,        "dimuon_vChi2/F");
  outTree->Branch("dimuon_triggerMatch", &out_dimuon_triggerMatch, "dimuon_triggerMatch/F");

  //2mu+2Trk vertexing
  outTree->Branch("dimuonditrk_vProb",      &out_dimuonditrk_vProb,        "dimuonditrk_vProb/F");
  outTree->Branch("dimuonditrk_vChi2",      &out_dimuonditrk_vChi2,        "dimuonditrk_vChi2/F");
  outTree->Branch("dimuonditrk_nDof",       &out_dimuonditrk_nDof,         "dimuonditrk_nDof/F");
  outTree->Branch("dimuonditrk_charge",     &out_dimuonditrk_charge,       "dimuonditrk_charge/F");

  outTree->Branch("dimuonditrk_cosAlpha",      &out_dimuonditrk_cosAlpha,        "dimuonditrk_cosAlpha/F");
  outTree->Branch("dimuonditrk_ctauPV",      &out_dimuonditrk_ctauPV,        "dimuonditrk_ctauPV/F");
  outTree->Branch("dimuonditrk_ctauErrPV",      &out_dimuonditrk_ctauErrPV,        "dimuonditrk_ctauErrPV/F");

  outTree->Branch("highKaonMatch",     &out_highKaonMatch,       "highKaonMatch/F");
  outTree->Branch("lowKaonMatch",     &out_lowKaonMatch,       "lowKaonMatch/F");
  outTree->Branch("lowMuonMatch",     &out_lowMuonMatch,       "lowMuonMatch/F");
  outTree->Branch("highMuonMatch",     &out_highMuonMatch,       "highMuonMatch/F");

  outTree->Branch("isBestCandidate",        &out_isBestCandidate,        "isBestCandidate/F");




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

  bool test = true;

  test = test && ((*highMuon_pt) >= 1.0) && ((*highMuon_pt) >= 1.0);

  //test = test && (*lowMuonMatch>0.0) && (*highMuonMatch>0.0);

  test = test && (*dimuonditrk_vProb> 0.01);

  //int a = (int) (*trigger);
  //std::cout << (*trigger);
  //

  if(test)
  {
    out_run  = (Float_t)(*run);
    out_event        = (Float_t)(*event);
    out_lumi         = (Float_t)(*lumi);
    out_numPrimaryVertices   = (Float_t)(*numPrimaryVertices);
    out_trigger      = (Float_t)(*trigger);
    out_noXCandidates        = (Float_t)(*noXCandidates);

    out_dimuonditrk_m        = (Float_t)(*dimuonditrk_m);
    out_dimuonditrk_m_rf     = (Float_t)(*dimuonditrk_m_rf);
    out_dimuonditrk_m_rf_c   = (Float_t)(*dimuonditrk_m_rf_c);
    out_dimuonditrk_pt       = (Float_t)(*dimuonditrk_pt);


    out_dimuon_m     = (Float_t)(*dimuon_m);
    out_dimuon_pt    = (Float_t)(*dimuon_pt);
    out_ditrak_m     = (Float_t)(*ditrak_m);
    out_ditrak_pt    = (Float_t)(*ditrak_pt);
    out_highKaon_pt  = (Float_t)(*highKaon_pt);
    out_lowKaon_pt   = (Float_t)(*lowKaon_pt);
    out_highMuon_pt  = (Float_t)(*highMuon_pt);
    out_lowMuon_pt   = (Float_t)(*lowMuon_pt);

    out_dimuon_vProb         = (Float_t)(*dimuon_vProb);
    out_dimuon_vChi2         = (Float_t)(*dimuon_vChi2);
    out_dimuon_triggerMatch  = (Float_t)(*dimuon_triggerMatch);

    out_dimuonditrk_vProb    = (Float_t)(*dimuonditrk_vProb);
    out_dimuonditrk_vChi2    = (Float_t)(*dimuonditrk_vChi2);
    out_dimuonditrk_nDof     = (Float_t)(*dimuonditrk_nDof);
    out_dimuonditrk_charge   = (Float_t)(*dimuonditrk_charge);

    out_dimuonditrk_cosAlpha         = (Float_t)(*dimuonditrk_cosAlpha);
    out_dimuonditrk_ctauPV   = (Float_t)(*dimuonditrk_ctauPV);
    out_dimuonditrk_ctauErrPV        = (Float_t)(*dimuonditrk_ctauErrPV);

    out_highKaonMatch        = (Float_t)(*highKaonMatch);
    out_lowKaonMatch         = (Float_t)(*lowKaonMatch);
    out_lowMuonMatch         = (Float_t)(*lowMuonMatch);
    out_highMuonMatch        = (Float_t)(*highMuonMatch);

    out_isBestCandidate      = (Float_t)(*isBestCandidate);

    outTree->Fill();
  }

  return kTRUE;
}

void TwoMuTwoK::SlaveTerminate()
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

void TwoMuTwoK::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

}
