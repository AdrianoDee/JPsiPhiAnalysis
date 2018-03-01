#define TwoMuTwoTrakTrigSkim_cxx
// The class definition in TwoMuTwoTrakTrigSkim.h has been generated automatically
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
// root> T->Process("TwoMuTwoTrakTrigSkim.C")
// root> T->Process("TwoMuTwoTrakTrigSkim.C","some options")
// root> T->Process("TwoMuTwoTrakTrigSkim.C+")
//


#include "TwoMuTwoTrakTrigSkim.h"
#include <TH2.h>
#include <TStyle.h>

void TwoMuTwoTrakTrigSkim::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void TwoMuTwoTrakTrigSkim::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   std::string outputString = "2Trak2MuTrig_tree.root";
   OutFile = new TProofOutputFile( outputString.data() );
   fOut = OutFile->OpenFile("RECREATE");
   if (!(fOut=OutFile->OpenFile("RECREATE")))
   {
     Warning("SlaveBegin","Problems opening file: %s%s", OutFile->GetDir(), OutFile->GetFileName() );
   }

   outTuple = new TNtuple("outuple","outuple","run:xM:ttM:mmM:xTrigM:ttTrigM:mmTrigM:muonp_pT:muonn_pT:kaonp_pT:kaonn_pT:matchMN:matchMP:matchKN:matchKP:vProb");


}

Bool_t TwoMuTwoTrakTrigSkim::Process(Long64_t entry)
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

   Float_t run_out,ttM,mmM,xM;
   Float_t xTrigM,ttTrigM,mmTrigM;
   Float_t matchKP,matchKN,matchMN,matchMP;
   Float_t muonp_pT, muonn_pT, kaonn_pT, kaonp_pT, vProb_out;

   fReader.SetEntry(entry);

   bool phiMass = (*ditrakTrigger_p4).M() > 0.9 && (*ditrakTrigger_p4).M() < 1.3;
   bool jpsiMass = (*dimuonTrigger_p4).M() > 2.88 && (*dimuonTrigger_p4).M() < 3.32;
   bool xMass = (*dimuonditrkTrigger_p4).M() > 4.0 && (*dimuonditrkTrigger_p4).M() < 6.0;

   std::bitset<16> tOne(*muonN_tMatch);
   std::bitset<16> tTwo(*muonP_tMatch);
   std::bitset<16> tThree(*trakN_tMatch);
   std::bitset<16> tFour(*trakP_tMatch);

   std::bitset<16> theTrig(*trigger);

   if(phiMass && xMass && jpsiMass)// && tOne.test(0) && tTwo.test(0)  && tThree.test(0) && tFour.test(0))
   {
     run_out = (*run);

     xM = (*dimuonditrk_p4).M();
     ttM = (*ditrak_p4).M();
     mmM = (*dimuon_p4).M();

     xTrigM = (*dimuonditrkTrigger_p4).M();
     ttTrigM = (*ditrakTrigger_p4).M();
     mmTrigM = (*dimuonTrigger_p4).M();

     muonp_pT = (*muonp_p4).Pt();
     muonn_pT = (*muonn_p4).Pt();
     kaonp_pT = (*kaonp_p4).Pt();
     kaonn_pT = (*kaonn_p4).Pt();

     matchMN = (*muonN_tMatch);
     matchMP = (*muonP_tMatch);
     matchKN = (*trakN_tMatch);
     matchKP = (*trakP_tMatch);

     vProb_out = (*dimuonditrk_vProb);

     Float_t params[16] = {run_out,xM,ttM,mmM,xTrigM,ttTrigM,mmTrigM,muonp_pT,muonn_pT,kaonp_pT,kaonn_pT,matchMN,matchMP,matchKN,matchKP,vProb_out};
     outTuple->Fill(params);
   }

   return kTRUE;
}

void TwoMuTwoTrakTrigSkim::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

   TDirectory *savedir = gDirectory;
   if (fOut)
   {
     fOut->cd();
     gStyle->SetOptStat(111111) ;


     outTuple->Write();
     OutFile->Print();
     fOutput->Add(OutFile);
     gDirectory = savedir;
     fOut->Close();

   }

}

void TwoMuTwoTrakTrigSkim::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}
