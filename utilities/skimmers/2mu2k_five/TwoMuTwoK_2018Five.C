#define TwoMuTwoK_2018Five_cxx
// The class definition in TwoMuTwoK_2018Five.h has been generated automatically
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
// root> T->Process("TwoMuTwoK_2018Five.C")
// root> T->Process("TwoMuTwoK_2018Five.C","some options")
// root> T->Process("TwoMuTwoK_2018Five.C+")
//


#include "TwoMuTwoK_2018Five.h"
#include <TH2.h>
#include <bitset>

void TwoMuTwoK_2018Five::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
}

void TwoMuTwoK_2018Five::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).
  TString option = GetOption();

  std::string outputString = "2mu2k_tree.root";
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

  outTuple = fChain->CloneTree(0);//new TNtuple("outuple","outuple","run:evt:xM:ttM:mmM:xM_ref:ttM_ref:mmM_ref:xL:xPt:xEta:xVtx:xCos:xHlt:muonp_pT:muonn_pT:kaonn_pT:kaonp_pT:mmPt:ttPt");

}

Bool_t TwoMuTwoK_2018Five::Process(Long64_t entry)
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

  Float_t run_out, evt_out, xM_out, ttM_out, mmM_out, xM_ref_out, ttM_ref_out, mmM_ref_out;
  Float_t xL_out, xPt_out, xEta_out, xVtx_out, xCos_out, xHlt_out,muonp_pT_out, muonn_pT_out, kaonn_pT_out, kaonp_pT_out;
  Float_t mmPt, ttPt;

  fReader.SetEntry(entry);

  //std::bitset<16> tt(*trigger);
  //std::bitset<16> dimuonTmatch(*dimuon_triggerMatch);
  bool phiM = (*ditrak_p4).M() > 1.01 && (*ditrak_p4).M() < 1.03;
  bool jpsiM = (*dimuon_p4).M() > 3.00 && (*dimuon_p4).M() < 3.20;
  bool cosAlpha = (*dimuonditrk_cosAlpha) > 0.9;
  bool vertexP = (*dimuonditrk_vProb) > 0.01;
  bool jPT = (*dimuon_p4).Pt() > 7.0;
  bool pPT = (*ditrak_p4).Pt() > 1.0;
  bool theTrigger = (*trigger) > 0;
  //bool tMatchDimuon = dimuonTmatch.test(0);
  bool isMatched = (*dimuon_triggerMatch)>0 && (*highKaonMatch>0) && (*lowKaonMatch>0);
  bool triggerBit = (*trigger) > 0; //tt.test(0);
  bool isBest = (*isBestCandidate);

  if(theTrigger && phiM && jpsiM && cosAlpha && vertexP && isMatched && jPT && pPT)
  {
    // run_out =  *run;
    // evt_out =  *event;
    //
    // ttM_out = (*ditrak_p4).M();
    // mmM_out= (*dimuon_p4).M();
    // xM_out   = (*dimuonditrk_p4).M();
    //
    // ttM_ref_out = (*ditrak_rf_p4).M();
    // mmM_ref_out= (*dimuon_rf_p4).M();
    // xM_ref_out   = (*dimuonditrk_rf_p4).M();
    //
    // xL_out = (*dimuonditrk_ctauPV)/(*dimuonditrk_ctauErrPV);
    // xPt_out = (*dimuonditrk_rf_p4).Pt();
    // xEta_out = (*dimuonditrk_rf_p4).Eta();
    // xVtx_out = *dimuonditrk_vProb;
    // xCos_out = *dimuonditrk_cosAlpha;
    // xHlt_out = *trigger;
    //
    // muonp_pT_out = (*muonp_rf_p4).Pt();
    // muonn_pT_out = (*muonn_rf_p4).Pt();
    // kaonn_pT_out = (*kaonp_rf_p4).Pt();
    // kaonp_pT_out = (*kaonn_rf_p4).Pt();
    //
    // mmPt = (*dimuon_p4).Pt();
    // ttPt = (*ditrak_p4).Pt();
    //
    // Float_t params[] = {run_out,evt_out,xM_out,ttM_out,mmM_out,xM_ref_out,ttM_ref_out,mmM_ref_out,
    // xL_out,xPt_out,xEta_out,xVtx_out,xCos_out,xHlt_out,muonp_pT_out,muonn_pT_out,kaonn_pT_out,kaonp_pT_out,mmPt,ttPt};

    outTuple->Fill();
  }


  return kTRUE;
}

void TwoMuTwoK_2018Five::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  TDirectory *savedir = gDirectory;
  if (fOut)
  {
    fOut->cd();
    //gStyle->SetOptStat(111111) ;


    outTuple->Write();
    OutFile->Print();
    fOutput->Add(OutFile);
    gDirectory = savedir;
    fOut->Close();

  }

}

void TwoMuTwoK_2018Five::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

}
