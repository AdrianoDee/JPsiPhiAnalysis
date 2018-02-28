#define TrakTrigSkim_cxx
// The class definition in TrakTrigSkim.h has been generated automatically
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
// root> T->Process("TrakTrigSkim.C")
// root> T->Process("TrakTrigSkim.C","some options")
// root> T->Process("TrakTrigSkim.C+")
//


#include "TrakTrigSkim.h"
#include <TH2.h>
#include <TStyle.h>

void TrakTrigSkim::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void TrakTrigSkim::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   std::string outputString = "2Trak2Trig_tree.root";
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

   outTuple = new TNtuple("outuple","outuple","run:ttM:trigtrigM:trigp_pT:trign_pT");


}

Bool_t TrakTrigSkim::Process(Long64_t entry)
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
   Float_t run_out, ttM, trigtrigM;
   Float_t trigp_pT, trign_pT;

   bool trigMass = (*ditrig_p4).M() < 1.31 && (*ditrig_p4).M() > 0.94;

   fReader.SetEntry(entry);

   if(trigMass)
   {
     run_out = (*run);
     ttM = (*ditrak_p4).M();
     trigtrigM = (*ditrig_p4).M();
     trigp_pT = (*trigP_p4).Pt();
     trign_pT = (*trigN_p4).Pt();

     outTuple->Fill(run_out,ttM,trigtrigM,trigp_pT,trign_pT);
   }

   return kTRUE;
}

void TrakTrigSkim::SlaveTerminate()
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

void TrakTrigSkim::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}
