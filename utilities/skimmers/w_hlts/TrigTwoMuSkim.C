#define TrigTwoMuSkim_cxx
// The class definition in TrigTwoMuSkim.h has been generated automatically
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
// root> T->Process("TrigTwoMuSkim.C")
// root> T->Process("TrigTwoMuSkim.C","some options")
// root> T->Process("TrigTwoMuSkim.C+")
//


#include "TrigTwoMuSkim.h"
#include <TH2.h>
#include <TStyle.h>

void TrigTwoMuSkim::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void TrigTwoMuSkim::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   std::string outputString = "2Mu2Trig_tree.root";
   OutFile = new TProofOutputFile( outputString.data() );
   fOut = OutFile->OpenFile("RECREATE");
   if (!(fOut=OutFile->OpenFile("RECREATE")))
   {
     Warning("SlaveBegin","Problems opening file: %s%s", OutFile->GetDir(), OutFile->GetFileName() );
   }

   outTuple = new TNtuple("outuple","outuple","run:vProb:mmM:trigtrigM:trigp_pT:trign_pT:matchOne:matchTwo");



}

Bool_t TrigTwoMuSkim::Process(Long64_t entry)
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

   Float_t run_out, mmM, trigtrigM;
   Float_t trigp_pT, trign_pT,vProb_out;
   UInt_t matchOne, matchTwo;
   fReader.SetEntry(entry);

   bool trigMass = (*dimuonTrigger_p4).M() > 2.88 && (*dimuonTrigger_p4).M() < 3.32;

   std::bitset<16> tOne(*tMatchN);
   std::bitset<16> tTwo(*tMatchP);
   std::bitset<16> theTrig(*trigger);
   
   int triggerToTest = 0;
 
   if(trigMass && tOne.test(triggerToTest) && tTwo.test(triggerToTest) && theTrig.test(triggerToTest))
   {
     run_out = (*run);
     mmM = (*dimuon_p4).M();
     trigtrigM = (*dimuonTrigger_p4).M();
     trigp_pT = (*muonP_p4).Pt();
     trign_pT = (*muonN_p4).Pt();
     matchOne = (*tMatchN);
     matchTwo = (*tMatchP);
     vProb_out = (*vProb);

     outTuple->Fill(run_out,vProb_out,mmM,trigtrigM,trigp_pT,trign_pT,matchOne,matchTwo);
   }

   return kTRUE;
}

void TrigTwoMuSkim::SlaveTerminate()
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

void TrigTwoMuSkim::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}
