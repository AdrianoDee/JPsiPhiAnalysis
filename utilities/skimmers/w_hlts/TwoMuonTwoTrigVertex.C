#define TwoMuonTwoTrigVertex_cxx
// The class definition in TwoMuonTwoTrigVertex.h has been generated automatically
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
// root> T->Process("TwoMuonTwoTrigVertex.C")
// root> T->Process("TwoMuonTwoTrigVertex.C","some options")
// root> T->Process("TwoMuonTwoTrigVertex.C+")
//


#include "TwoMuonTwoTrigVertex.h"
#include <TH2.h>
#include <TStyle.h>

void TwoMuonTwoTrigVertex::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void TwoMuonTwoTrigVertex::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   std::string outputString = "2Trak2MuonVertexTrig_tree.root";
   OutFile = new TProofOutputFile( outputString.data() );
   fOut = OutFile->OpenFile("RECREATE");
   if (!(fOut=OutFile->OpenFile("RECREATE")))
   {
     Warning("SlaveBegin","Problems opening file: %s%s", OutFile->GetDir(), OutFile->GetFileName() );
   }

   outTuple = new TNtuple("outuple","outuple","run:xM:ttM:mmM:xTrigM:ttTrigM:mmTrigM:traktrakHLT:trakHLT:phiHLT:matchMN:matchMP:matchKN:matchKP:vProb:lxysig");


}

Bool_t TwoMuonTwoTrigVertex::Process(Long64_t entry)
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

   if(run.Get())
   {
   std::bitset<16> tMuonOne(*muonN_tMatch);
   std::bitset<16> tMuonTwo(*muonP_tMatch);
   std::bitset<16> tKaonOne(*trakN_tMatch);
   std::bitset<16> tKaonTwo(*trakP_tMatch);

   bool phiHLT      = tMuonOne.test(0) && tMuonTwo.test(0) && tKaonOne.test(0) && tKaonTwo.test(0);
   bool traktrakHLT = tMuonOne.test(1) && tMuonTwo.test(1) && tKaonOne.test(1) && tKaonTwo.test(1);
   bool trakHLT     = tMuonOne.test(2) && tMuonTwo.test(2) && (tKaonOne.test(2) || tKaonTwo.test(2));

   int triggerToTest = 0;

   // // float run_out = float(*run);
   // float xM      = float((*dimuonditrak_p4).M());
   // float ttM     = float((*ditrak_p4).M());
   // float mmM     = float((*dimuon_p4).M());
   float params[16] = {float(*run),float((*dimuonditrak_p4).M()),float((*ditrak_p4).M()),float((*dimuon_p4).M()),float((*dimuonditrkTrigger_p4).M()),float((*dimuonTrigger_p4).M()),
     float((*ditrakTrigger_p4).M()),float(traktrakHLT),float(trakHLT),float(phiHLT),float(*muonN_tMatch),float(*muonP_tMatch),float(*trakN_tMatch),float(*trakP_tMatch),float(*dimuonditrk_vProb),float(*dimuonditrk_lxy)/float(*dimuonditrk_lxyErr)};

   if(phiHLT || traktrakHLT || trakHLT)
    // triggerToTest++;
     outTuple->Fill(params);
   }

   return kTRUE;
}

void TwoMuonTwoTrigVertex::SlaveTerminate()
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

void TwoMuonTwoTrigVertex::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}
