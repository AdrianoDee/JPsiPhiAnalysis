#define DiTrakSkim_cxx
// The class definition in DiTrakSkim.h has been generated automatically
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
// root> T->Process("DiTrakSkim.C")
// root> T->Process("DiTrakSkim.C","some oPtions")
// root> T->Process("DiTrakSkim.C+")
//


#include "DiTrakSkim.h"
#include <TH2.h>
#include <TStyle.h>

void DiTrakSkim::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void DiTrakSkim::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   std::string outputString = "2Trak2TrigMatch_tree.root";
   OutFile = new TProofOutputFile( outputString.data() );
   fOut = OutFile->OpenFile("RECREATE");
   if (!(fOut=OutFile->OpenFile("RECREATE")))
   {
     Warning("SlaveBegin","Problems opening file: %s%s", OutFile->GetDir(), OutFile->GetFileName() );
   }

   outTuple = new TNtuple("outuple","outuple","run:ttM:trigtrigM:trakP_Pt:trakN_Pt:trigP_Pt:trigN_Pt:matchOne:matchTwo");

}

Bool_t DiTrakSkim::Process(Long64_t entry)
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

   Float_t run_out, ttM, trigtrigM;
   Float_t trigN_Pt, trigP_Pt;
   Float_t trakN_Pt, trakP_Pt;
   UInt_t matchP = 0, matchN = 0;
   fReader.SetEntry(entry);

   std::bitset<16> theTrigger(*trigger);

   int triggerToTest = 0; //trigger-filter one to one

   TLorentzVector tP = (*trakP_p4);
   TLorentzVector tN = (*trakN_p4);
   std::vector<TLorentzVector> filteredTrigs;
   std::vector<int> filteredIndex;

   if(theTrigger.test(triggerToTest))
   {
     matchP = 0, matchN = 0;
     //Filtering Trigger objects
     for (size_t i = 0; i < trigs_filters.GetSize(); i++) {
       // std::bitset<16> theFilt(trigs_filters[i]);

       // if(theFilt.test(triggerToTest))
       // {
         TLorentzVector trig;

         trig.SetPtEtaPhiM(trigs_pt[i], trigs_eta[i], trigs_phi[i], trigs_m[i]);
         filteredTrigs.push_back(trig);
         filteredIndex.push_back(trigs_filters[i]);
       // }
     }

     //Matching the two tracks
     bool matchPos = false, matchNeg = false;
     //tPos
     TLorentzVector trigPos, trigNeg;

     for (size_t i = 0; i < filteredTrigs.size(); i++)
     {

       if(MatchByDRDPt(*trakP_p4,filteredTrigs[i]))
       {
         if(matchPos)
         {
           if(DeltaR(*trakP_p4,trigPos) > DeltaR(*trakP_p4,filteredTrigs[i]))
           {
             trigPos = filteredTrigs[i];
             matchP  = filteredIndex[i];
           }
         }
         else
         {
           trigPos = filteredTrigs[i];
           matchP  = filteredIndex[i];
         }

         matchPos = true;
       }
     }
     //tNeg

     for (size_t i = 0; i < filteredTrigs.size(); i++)
     {

       if(MatchByDRDPt(*trakN_p4,filteredTrigs[i]))
       {
         if(matchNeg)
         {
           if(DeltaR(*trakN_p4,trigNeg) > DeltaR(*trakN_p4,filteredTrigs[i]))
           {
             trigNeg = filteredTrigs[i];
             matchP  = filteredIndex[i];
           }
         }
         else
         {
           trigNeg = filteredTrigs[i];
           matchP  = filteredIndex[i];
         }
         matchNeg = true;
       }
     }

     TLorentzVector ditrig_p4 = (trigNeg) + (trigPos);

     bool matched = matchPos || matchNeg;
     std::bitset<16> filtP(matchP);
     std::bitset<16> filtN(matchN);
     bool testFilter = filtP.test(triggerToTest) || filtP.test(triggerToTest);
     bool trigMass = (ditrig_p4).M() < 1.31 && (ditrig_p4).M() > 0.94;

     if(matched && trigMass && testFilter)
     {
       run_out = (*run);
       ttM = (*ditrak_p4).M();
       trigtrigM = (ditrig_p4).M();
       trakP_Pt = (*trakP_p4).Pt();
       trakN_Pt = (*trakN_p4).Pt();
       trigN_Pt = (trigNeg).Pt();
       trigP_Pt = (trigPos).Pt();

       outTuple->Fill(run_out,ttM,trigtrigM,trakP_Pt,trakN_Pt,trigP_Pt,trigN_Pt,matchP,matchN);

     }

 }

   return kTRUE;
}

void DiTrakSkim::SlaveTerminate()
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

float DiTrakSkim::DeltaR(const TLorentzVector t1, const TLorentzVector t2)
{
   float p1 = t1.Phi();
   float p2 = t2.Phi();
   float e1 = t1.Eta();
   float e2 = t2.Eta();

   auto dp=std::abs(p1-p2); if (dp>float(M_PI)) dp-=float(2*M_PI);

   return sqrt((e1-e2)*(e1-e2) + dp*dp);
}

bool DiTrakSkim::MatchByDRDPt(const TLorentzVector t1, const TLorentzVector t2)
{
  return (fabs(t1.Pt()-t2.Pt())/t2.Pt()<maxDPtRel &&
        DeltaR(t1,t2) < maxDeltaR);
}

void DiTrakSkim::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}
