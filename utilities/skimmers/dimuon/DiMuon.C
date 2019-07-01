#define Dimuon_cxx
// The class definition in Dimuon.h has been generated automatically
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
// root> T->Process("Dimuon.C")
// root> T->Process("Dimuon.C","some options")
// root> T->Process("Dimuon.C+")
//


#include "Dimuon.h"
#include <TH2.h>
#include <TStyle.h>

void Dimuon::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void Dimuon::SlaveBegin(TTree * /*tree*/)
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

}

Bool_t Dimuon::Process(Long64_t entry)
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


   out_run =       (Float_t)(*run);
  out_event =     (Float_t)(*event);
  out_lumiblock =         (Float_t)(*lumiblock);
  out_ndimuon =   (Float_t)(*ndimuon);
  out_nmuons =    (Float_t)(*nmuons);
  out_trigger =   (Float_t)(*trigger);
  out_charge =    (Float_t)(*charge);
  out_tMatch =    (Float_t)(*tMatch);
  out_dimuon_p4 =         (Float_t)(*dimuon_p4);
  out_highMuon_p4 =       (Float_t)(*highMuon_p4);
  out_lowMuon_p4 =        (Float_t)(*lowMuon_p4);
  out_dimuon_m =  (Float_t)(*dimuon_m);
  out_dimuon_pt =         (Float_t)(*dimuon_pt);
  out_dimuon_eta =        (Float_t)(*dimuon_eta);
  out_dimuon_phi =        (Float_t)(*dimuon_phi);
  out_dimuon_p =  (Float_t)(*dimuon_p);
  out_MassErr =   (Float_t)(*MassErr);
  out_vProb =     (Float_t)(*vProb);
  out_DCA =       (Float_t)(*DCA);
  out_ppdlPV =    (Float_t)(*ppdlPV);
  out_ppdlErrPV =         (Float_t)(*ppdlErrPV);
  out_ppdlBS =    (Float_t)(*ppdlBS);
  out_ppdlErrBS =         (Float_t)(*ppdlErrBS);
  out_cosAlpha =  (Float_t)(*cosAlpha);
  out_lxyPV =     (Float_t)(*lxyPV);
  out_lxyBS =     (Float_t)(*lxyBS);
  out_numPrimaryVertices =        (Float_t)(*numPrimaryVertices);

   return kTRUE;
}

void Dimuon::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void Dimuon::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}
