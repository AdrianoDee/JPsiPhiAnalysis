#define 2mu2k_2012_cxx
// The class definition in 2mu2k_2012.h has been generated automatically
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
// root> T->Process("2mu2k_2012.C")
// root> T->Process("2mu2k_2012.C","some options")
// root> T->Process("2mu2k_2012.C+")
//


#include "2mu2k_2012.h"
#include <TH2.h>
#include <TStyle.h>

void 2mu2k_2012::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void 2mu2k_2012::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();



   std::string outputString = "2mu2k_2012_tree.root";
   OutFile = new TProofOutputFile( outputString.data() );
   fOut = OutFile->OpenFile("RECREATE");
   if (!(fOut=OutFile->OpenFile("RECREATE")))
   {
     Warning("SlaveBegin","Problems opening file: %s%s", OutFile->GetDir(), OutFile->GetFileName() );
   }


  outTree->Branch("TrigRes",      &out_TrigRes,   "TrigRes/F");
  outTree->Branch("TrigNames",    &out_TrigNames,         "TrigNames/F");
  outTree->Branch("MatchTriggerNames",    &out_MatchTriggerNames,         "MatchTriggerNames/F");
  outTree->Branch("L1TrigRes",    &out_L1TrigRes,         "L1TrigRes/F");
  outTree->Branch("evtNum",       &out_evtNum,    "evtNum/F");
  outTree->Branch("runNum",       &out_runNum,    "runNum/F");
  outTree->Branch("lumiNum",      &out_lumiNum,   "lumiNum/F");
  outTree->Branch("priVtx_n",     &out_priVtx_n,  "priVtx_n/F");
  outTree->Branch("priVtx_X",     &out_priVtx_X,  "priVtx_X/F");
  outTree->Branch("priVtx_Y",     &out_priVtx_Y,  "priVtx_Y/F");
  outTree->Branch("priVtx_Z",     &out_priVtx_Z,  "priVtx_Z/F");
  outTree->Branch("priVtx_XE",    &out_priVtx_XE,         "priVtx_XE/F");
  outTree->Branch("priVtx_YE",    &out_priVtx_YE,         "priVtx_YE/F");
  outTree->Branch("priVtx_ZE",    &out_priVtx_ZE,         "priVtx_ZE/F");
  outTree->Branch("priVtx_NormChi2",      &out_priVtx_NormChi2,   "priVtx_NormChi2/F");
  outTree->Branch("priVtx_Chi2",  &out_priVtx_Chi2,       "priVtx_Chi2/F");
  outTree->Branch("priVtx_CL",    &out_priVtx_CL,         "priVtx_CL/F");
  outTree->Branch("priVtx_tracks",        &out_priVtx_tracks,     "priVtx_tracks/F");
  outTree->Branch("priVtx_tracksPtSq",    &out_priVtx_tracksPtSq,         "priVtx_tracksPtSq/F");

  outTree->Branch("XMass",        &out_XMass,     "XMass/F");
  outTree->Branch("XPx",  &out_XPx,       "XPx/F");
  outTree->Branch("XPy",  &out_XPy,       "XPy/F");
  outTree->Branch("XPz",  &out_XPz,       "XPz/F");
  outTree->Branch("XPxE",         &out_XPxE,      "XPxE/F");
  outTree->Branch("XPyE",         &out_XPyE,      "XPyE/F");
  outTree->Branch("XPzE",         &out_XPzE,      "XPzE/F");
  outTree->Branch("XVtx_CL",      &out_XVtx_CL,   "XVtx_CL/F");
  outTree->Branch("XVtx_Chi2",    &out_XVtx_Chi2,         "XVtx_Chi2/F");
  outTree->Branch("XDecayVtx_X",  &out_XDecayVtx_X,       "XDecayVtx_X/F");
  outTree->Branch("XDecayVtx_Y",  &out_XDecayVtx_Y,       "XDecayVtx_Y/F");
  outTree->Branch("XDecayVtx_Z",  &out_XDecayVtx_Z,       "XDecayVtx_Z/F");
  outTree->Branch("XDecayVtx_XE",         &out_XDecayVtx_XE,      "XDecayVtx_XE/F");
  outTree->Branch("XDecayVtx_YE",         &out_XDecayVtx_YE,      "XDecayVtx_YE/F");
  outTree->Branch("XDecayVtx_ZE",         &out_XDecayVtx_ZE,      "XDecayVtx_ZE/F");
  outTree->Branch("XCosAlphaBS",  &out_XCosAlphaBS,       "XCosAlphaBS/F");
  outTree->Branch("XCosAlpha3DBS",        &out_XCosAlpha3DBS,     "XCosAlpha3DBS/F");
  outTree->Branch("XCTauBS",      &out_XCTauBS,   "XCTauBS/F");
  outTree->Branch("XCTauBSE",     &out_XCTauBSE,  "XCTauBSE/F");
  outTree->Branch("XLxyBS",       &out_XLxyBS,    "XLxyBS/F");
  outTree->Branch("XLxyBSE",      &out_XLxyBSE,   "XLxyBSE/F");
  outTree->Branch("XLxyzBS",      &out_XLxyzBS,   "XLxyzBS/F");
  outTree->Branch("XLxyzBSE",     &out_XLxyzBSE,  "XLxyzBSE/F");
  outTree->Branch("XCosAlphaPV",  &out_XCosAlphaPV,       "XCosAlphaPV/F");
  outTree->Branch("XCosAlpha3DPV",        &out_XCosAlpha3DPV,     "XCosAlpha3DPV/F");
  outTree->Branch("XCTauPV",      &out_XCTauPV,   "XCTauPV/F");
  outTree->Branch("XCTauPVE",     &out_XCTauPVE,  "XCTauPVE/F");
  outTree->Branch("XLxyPV",       &out_XLxyPV,    "XLxyPV/F");
  outTree->Branch("XLxyPVE",      &out_XLxyPVE,   "XLxyPVE/F");
  outTree->Branch("XLxyzPV",      &out_XLxyzPV,   "XLxyzPV/F");
  outTree->Branch("XLxyzPVE",     &out_XLxyzPVE,  "XLxyzPVE/F");
  outTree->Branch("PriVtx_XCosAlpha_n",   &out_PriVtx_XCosAlpha_n,        "PriVtx_XCosAlpha_n/F");
  outTree->Branch("PriVtx_XCosAlpha_X",   &out_PriVtx_XCosAlpha_X,        "PriVtx_XCosAlpha_X/F");
  outTree->Branch("PriVtx_XCosAlpha_Y",   &out_PriVtx_XCosAlpha_Y,        "PriVtx_XCosAlpha_Y/F");
  outTree->Branch("PriVtx_XCosAlpha_Z",   &out_PriVtx_XCosAlpha_Z,        "PriVtx_XCosAlpha_Z/F");
  outTree->Branch("PriVtx_XCosAlpha_EX",  &out_PriVtx_XCosAlpha_EX,       "PriVtx_XCosAlpha_EX/F");
  outTree->Branch("PriVtx_XCosAlpha_EY",  &out_PriVtx_XCosAlpha_EY,       "PriVtx_XCosAlpha_EY/F");
  outTree->Branch("PriVtx_XCosAlpha_EZ",  &out_PriVtx_XCosAlpha_EZ,       "PriVtx_XCosAlpha_EZ/F");
  outTree->Branch("PriVtx_XCosAlpha_Chi2",        &out_PriVtx_XCosAlpha_Chi2,     "PriVtx_XCosAlpha_Chi2/F");
  outTree->Branch("PriVtx_XCosAlpha_CL",  &out_PriVtx_XCosAlpha_CL,       "PriVtx_XCosAlpha_CL/F");
  outTree->Branch("PriVtx_XCosAlpha_tracks",      &out_PriVtx_XCosAlpha_tracks,   "PriVtx_XCosAlpha_tracks/F");
  outTree->Branch("XCosAlphaPVCosAlpha",  &out_XCosAlphaPVCosAlpha,       "XCosAlphaPVCosAlpha/F");
  outTree->Branch("XCosAlpha3DPVCosAlpha",        &out_XCosAlpha3DPVCosAlpha,     "XCosAlpha3DPVCosAlpha/F");
  outTree->Branch("XCTauPVCosAlpha",      &out_XCTauPVCosAlpha,   "XCTauPVCosAlpha/F");
  outTree->Branch("XCTauPVCosAlphaE",     &out_XCTauPVCosAlphaE,  "XCTauPVCosAlphaE/F");
  outTree->Branch("XLxyPVCosAlpha",       &out_XLxyPVCosAlpha,    "XLxyPVCosAlpha/F");
  outTree->Branch("XLxyPVCosAlphaE",      &out_XLxyPVCosAlphaE,   "XLxyPVCosAlphaE/F");
  outTree->Branch("XLxyzPVCosAlpha",      &out_XLxyzPVCosAlpha,   "XLxyzPVCosAlpha/F");
  outTree->Branch("XLxyzPVCosAlphaE",     &out_XLxyzPVCosAlphaE,  "XLxyzPVCosAlphaE/F");
  outTree->Branch("PriVtx_XCosAlpha3D_n",         &out_PriVtx_XCosAlpha3D_n,      "PriVtx_XCosAlpha3D_n/F");
  outTree->Branch("PriVtx_XCosAlpha3D_X",         &out_PriVtx_XCosAlpha3D_X,      "PriVtx_XCosAlpha3D_X/F");
  outTree->Branch("PriVtx_XCosAlpha3D_Y",         &out_PriVtx_XCosAlpha3D_Y,      "PriVtx_XCosAlpha3D_Y/F");
  outTree->Branch("PriVtx_XCosAlpha3D_Z",         &out_PriVtx_XCosAlpha3D_Z,      "PriVtx_XCosAlpha3D_Z/F");
  outTree->Branch("PriVtx_XCosAlpha3D_EX",        &out_PriVtx_XCosAlpha3D_EX,     "PriVtx_XCosAlpha3D_EX/F");
  outTree->Branch("PriVtx_XCosAlpha3D_EY",        &out_PriVtx_XCosAlpha3D_EY,     "PriVtx_XCosAlpha3D_EY/F");
  outTree->Branch("PriVtx_XCosAlpha3D_EZ",        &out_PriVtx_XCosAlpha3D_EZ,     "PriVtx_XCosAlpha3D_EZ/F");
  outTree->Branch("PriVtx_XCosAlpha3D_Chi2",      &out_PriVtx_XCosAlpha3D_Chi2,   "PriVtx_XCosAlpha3D_Chi2/F");
  outTree->Branch("PriVtx_XCosAlpha3D_CL",        &out_PriVtx_XCosAlpha3D_CL,     "PriVtx_XCosAlpha3D_CL/F");
  outTree->Branch("PriVtx_XCosAlpha3D_tracks",    &out_PriVtx_XCosAlpha3D_tracks,         "PriVtx_XCosAlpha3D_tracks/F");
  outTree->Branch("XCosAlphaPVCosAlpha3D",        &out_XCosAlphaPVCosAlpha3D,     "XCosAlphaPVCosAlpha3D/F");
  outTree->Branch("XCosAlpha3DPVCosAlpha3D",      &out_XCosAlpha3DPVCosAlpha3D,   "XCosAlpha3DPVCosAlpha3D/F");
  outTree->Branch("XCTauPVCosAlpha3D",    &out_XCTauPVCosAlpha3D,         "XCTauPVCosAlpha3D/F");
  outTree->Branch("XCTauPVCosAlpha3DE",   &out_XCTauPVCosAlpha3DE,        "XCTauPVCosAlpha3DE/F");
  outTree->Branch("XLxyPVCosAlpha3D",     &out_XLxyPVCosAlpha3D,  "XLxyPVCosAlpha3D/F");
  outTree->Branch("XLxyPVCosAlpha3DE",    &out_XLxyPVCosAlpha3DE,         "XLxyPVCosAlpha3DE/F");
  outTree->Branch("XLxyzPVCosAlpha3D",    &out_XLxyzPVCosAlpha3D,         "XLxyzPVCosAlpha3D/F");
  outTree->Branch("XLxyzPVCosAlpha3DE",   &out_XLxyzPVCosAlpha3DE,        "XLxyzPVCosAlpha3DE/F");
  outTree->Branch("XLessPV_tracksPtSq",   &out_XLessPV_tracksPtSq,        "XLessPV_tracksPtSq/F");
  outTree->Branch("XLessPV_4tracksPtSq",  &out_XLessPV_4tracksPtSq,       "XLessPV_4tracksPtSq/F");
  outTree->Branch("PriVtxXLess_n",        &out_PriVtxXLess_n,     "PriVtxXLess_n/F");
  outTree->Branch("PriVtxXLess_X",        &out_PriVtxXLess_X,     "PriVtxXLess_X/F");
  outTree->Branch("PriVtxXLess_Y",        &out_PriVtxXLess_Y,     "PriVtxXLess_Y/F");
  outTree->Branch("PriVtxXLess_Z",        &out_PriVtxXLess_Z,     "PriVtxXLess_Z/F");
  outTree->Branch("PriVtxXLess_EX",       &out_PriVtxXLess_EX,    "PriVtxXLess_EX/F");
  outTree->Branch("PriVtxXLess_EY",       &out_PriVtxXLess_EY,    "PriVtxXLess_EY/F");
  outTree->Branch("PriVtxXLess_EZ",       &out_PriVtxXLess_EZ,    "PriVtxXLess_EZ/F");
  outTree->Branch("PriVtxXLess_Chi2",     &out_PriVtxXLess_Chi2,  "PriVtxXLess_Chi2/F");
  outTree->Branch("PriVtxXLess_CL",       &out_PriVtxXLess_CL,    "PriVtxXLess_CL/F");
  outTree->Branch("PriVtxXLess_tracks",   &out_PriVtxXLess_tracks,        "PriVtxXLess_tracks/F");
  outTree->Branch("XCosAlphaXLessPV",     &out_XCosAlphaXLessPV,  "XCosAlphaXLessPV/F");
  outTree->Branch("XCosAlpha3DXLessPV",   &out_XCosAlpha3DXLessPV,        "XCosAlpha3DXLessPV/F");
  outTree->Branch("XCTauXLessPV",         &out_XCTauXLessPV,      "XCTauXLessPV/F");
  outTree->Branch("XCTauXLessPVE",        &out_XCTauXLessPVE,     "XCTauXLessPVE/F");
  outTree->Branch("XLxyXLessPV",  &out_XLxyXLessPV,       "XLxyXLessPV/F");
  outTree->Branch("XLxyXLessPVE",         &out_XLxyXLessPVE,      "XLxyXLessPVE/F");
  outTree->Branch("XLxyzXLessPV",         &out_XLxyzXLessPV,      "XLxyzXLessPV/F");
  outTree->Branch("XLxyzXLessPVE",        &out_XLxyzXLessPVE,     "XLxyzXLessPVE/F");
  outTree->Branch("PriVtxXLess_XCosAlpha_n",      &out_PriVtxXLess_XCosAlpha_n,   "PriVtxXLess_XCosAlpha_n/F");
  outTree->Branch("PriVtxXLess_XCosAlpha_X",      &out_PriVtxXLess_XCosAlpha_X,   "PriVtxXLess_XCosAlpha_X/F");
  outTree->Branch("PriVtxXLess_XCosAlpha_Y",      &out_PriVtxXLess_XCosAlpha_Y,   "PriVtxXLess_XCosAlpha_Y/F");
  outTree->Branch("PriVtxXLess_XCosAlpha_Z",      &out_PriVtxXLess_XCosAlpha_Z,   "PriVtxXLess_XCosAlpha_Z/F");
  outTree->Branch("PriVtxXLess_XCosAlpha_EX",     &out_PriVtxXLess_XCosAlpha_EX,  "PriVtxXLess_XCosAlpha_EX/F");
  outTree->Branch("PriVtxXLess_XCosAlpha_EY",     &out_PriVtxXLess_XCosAlpha_EY,  "PriVtxXLess_XCosAlpha_EY/F");
  outTree->Branch("PriVtxXLess_XCosAlpha_EZ",     &out_PriVtxXLess_XCosAlpha_EZ,  "PriVtxXLess_XCosAlpha_EZ/F");
  outTree->Branch("PriVtxXLess_XCosAlpha_Chi2",   &out_PriVtxXLess_XCosAlpha_Chi2,        "PriVtxXLess_XCosAlpha_Chi2/F");
  outTree->Branch("PriVtxXLess_XCosAlpha_CL",     &out_PriVtxXLess_XCosAlpha_CL,  "PriVtxXLess_XCosAlpha_CL/F");
  outTree->Branch("PriVtxXLess_XCosAlpha_tracks",         &out_PriVtxXLess_XCosAlpha_tracks,      "PriVtxXLess_XCosAlpha_tracks/F");
  outTree->Branch("XCosAlphaXLessPVCosAlpha",     &out_XCosAlphaXLessPVCosAlpha,  "XCosAlphaXLessPVCosAlpha/F");
  outTree->Branch("XCosAlpha3DXLessPVCosAlpha",   &out_XCosAlpha3DXLessPVCosAlpha,        "XCosAlpha3DXLessPVCosAlpha/F");
  outTree->Branch("XCTauXLessPVCosAlpha",         &out_XCTauXLessPVCosAlpha,      "XCTauXLessPVCosAlpha/F");
  outTree->Branch("XCTauXLessPVCosAlphaE",        &out_XCTauXLessPVCosAlphaE,     "XCTauXLessPVCosAlphaE/F");
  outTree->Branch("XLxyXLessPVCosAlpha",  &out_XLxyXLessPVCosAlpha,       "XLxyXLessPVCosAlpha/F");
  outTree->Branch("XLxyXLessPVCosAlphaE",         &out_XLxyXLessPVCosAlphaE,      "XLxyXLessPVCosAlphaE/F");
  outTree->Branch("XLxyzXLessPVCosAlpha",         &out_XLxyzXLessPVCosAlpha,      "XLxyzXLessPVCosAlpha/F");
  outTree->Branch("XLxyzXLessPVCosAlphaE",        &out_XLxyzXLessPVCosAlphaE,     "XLxyzXLessPVCosAlphaE/F");
  outTree->Branch("PriVtxXLess_XCosAlpha3D_n",    &out_PriVtxXLess_XCosAlpha3D_n,         "PriVtxXLess_XCosAlpha3D_n/F");
  outTree->Branch("PriVtxXLess_XCosAlpha3D_X",    &out_PriVtxXLess_XCosAlpha3D_X,         "PriVtxXLess_XCosAlpha3D_X/F");
  outTree->Branch("PriVtxXLess_XCosAlpha3D_Y",    &out_PriVtxXLess_XCosAlpha3D_Y,         "PriVtxXLess_XCosAlpha3D_Y/F");
  outTree->Branch("PriVtxXLess_XCosAlpha3D_Z",    &out_PriVtxXLess_XCosAlpha3D_Z,         "PriVtxXLess_XCosAlpha3D_Z/F");
  outTree->Branch("PriVtxXLess_XCosAlpha3D_EX",   &out_PriVtxXLess_XCosAlpha3D_EX,        "PriVtxXLess_XCosAlpha3D_EX/F");
  outTree->Branch("PriVtxXLess_XCosAlpha3D_EY",   &out_PriVtxXLess_XCosAlpha3D_EY,        "PriVtxXLess_XCosAlpha3D_EY/F");
  outTree->Branch("PriVtxXLess_XCosAlpha3D_EZ",   &out_PriVtxXLess_XCosAlpha3D_EZ,        "PriVtxXLess_XCosAlpha3D_EZ/F");
  outTree->Branch("PriVtxXLess_XCosAlpha3D_Chi2",         &out_PriVtxXLess_XCosAlpha3D_Chi2,      "PriVtxXLess_XCosAlpha3D_Chi2/F");
  outTree->Branch("PriVtxXLess_XCosAlpha3D_CL",   &out_PriVtxXLess_XCosAlpha3D_CL,        "PriVtxXLess_XCosAlpha3D_CL/F");
  outTree->Branch("PriVtxXLess_XCosAlpha3D_tracks",       &out_PriVtxXLess_XCosAlpha3D_tracks,    "PriVtxXLess_XCosAlpha3D_tracks/F");
  outTree->Branch("XCosAlphaXLessPVCosAlpha3D",   &out_XCosAlphaXLessPVCosAlpha3D,        "XCosAlphaXLessPVCosAlpha3D/F");
  outTree->Branch("XCosAlpha3DXLessPVCosAlpha3D",         &out_XCosAlpha3DXLessPVCosAlpha3D,      "XCosAlpha3DXLessPVCosAlpha3D/F");
  outTree->Branch("XCTauXLessPVCosAlpha3D",       &out_XCTauXLessPVCosAlpha3D,    "XCTauXLessPVCosAlpha3D/F");
  outTree->Branch("XCTauXLessPVCosAlpha3DE",      &out_XCTauXLessPVCosAlpha3DE,   "XCTauXLessPVCosAlpha3DE/F");
  outTree->Branch("XLxyXLessPVCosAlpha3D",        &out_XLxyXLessPVCosAlpha3D,     "XLxyXLessPVCosAlpha3D/F");
  outTree->Branch("XLxyXLessPVCosAlpha3DE",       &out_XLxyXLessPVCosAlpha3DE,    "XLxyXLessPVCosAlpha3DE/F");
  outTree->Branch("XLxyzXLessPVCosAlpha3D",       &out_XLxyzXLessPVCosAlpha3D,    "XLxyzXLessPVCosAlpha3D/F");
  outTree->Branch("XLxyzXLessPVCosAlpha3DE",      &out_XLxyzXLessPVCosAlpha3DE,   "XLxyzXLessPVCosAlpha3DE/F");
  outTree->Branch("PriVtxXCorr_n",        &out_PriVtxXCorr_n,     "PriVtxXCorr_n/F");
  outTree->Branch("PriVtxXCorr_X",        &out_PriVtxXCorr_X,     "PriVtxXCorr_X/F");
  outTree->Branch("PriVtxXCorr_Y",        &out_PriVtxXCorr_Y,     "PriVtxXCorr_Y/F");
  outTree->Branch("PriVtxXCorr_Z",        &out_PriVtxXCorr_Z,     "PriVtxXCorr_Z/F");
  outTree->Branch("PriVtxXCorr_EX",       &out_PriVtxXCorr_EX,    "PriVtxXCorr_EX/F");
  outTree->Branch("PriVtxXCorr_EY",       &out_PriVtxXCorr_EY,    "PriVtxXCorr_EY/F");
  outTree->Branch("PriVtxXCorr_EZ",       &out_PriVtxXCorr_EZ,    "PriVtxXCorr_EZ/F");
  outTree->Branch("PriVtxXCorr_Chi2",     &out_PriVtxXCorr_Chi2,  "PriVtxXCorr_Chi2/F");
  outTree->Branch("PriVtxXCorr_CL",       &out_PriVtxXCorr_CL,    "PriVtxXCorr_CL/F");
  outTree->Branch("PriVtxXCorr_tracks",   &out_PriVtxXCorr_tracks,        "PriVtxXCorr_tracks/F");
  outTree->Branch("XCosAlphaPVX",         &out_XCosAlphaPVX,      "XCosAlphaPVX/F");
  outTree->Branch("XCTauPVX",     &out_XCTauPVX,  "XCTauPVX/F");
  outTree->Branch("XCTauPVXE",    &out_XCTauPVXE,         "XCTauPVXE/F");
  outTree->Branch("XLxyPVX",      &out_XLxyPVX,   "XLxyPVX/F");
  outTree->Branch("XLxyzPVX",     &out_XLxyzPVX,  "XLxyzPVX/F");
  outTree->Branch("XCTauPVX_3D",  &out_XCTauPVX_3D,       "XCTauPVX_3D/F");
  outTree->Branch("XCTauPVX_3D_err",      &out_XCTauPVX_3D_err,   "XCTauPVX_3D_err/F");
  outTree->Branch("kaon1_dxy_PV",         &out_kaon1_dxy_PV,      "kaon1_dxy_PV/F");
  outTree->Branch("kaon1_dz_PV",  &out_kaon1_dz_PV,       "kaon1_dz_PV/F");
  outTree->Branch("kaon2_dxy_PV",         &out_kaon2_dxy_PV,      "kaon2_dxy_PV/F");
  outTree->Branch("kaon2_dz_PV",  &out_kaon2_dz_PV,       "kaon2_dz_PV/F");
  outTree->Branch("kaon1_dxy_BS",         &out_kaon1_dxy_BS,      "kaon1_dxy_BS/F");
  outTree->Branch("kaon1_dz_BS",  &out_kaon1_dz_BS,       "kaon1_dz_BS/F");
  outTree->Branch("kaon2_dxy_BS",         &out_kaon2_dxy_BS,      "kaon2_dxy_BS/F");
  outTree->Branch("kaon2_dz_BS",  &out_kaon2_dz_BS,       "kaon2_dz_BS/F");
  outTree->Branch("kaon1_dxy_XLessPV",    &out_kaon1_dxy_XLessPV,         "kaon1_dxy_XLessPV/F");
  outTree->Branch("kaon1_dz_XLessPV",     &out_kaon1_dz_XLessPV,  "kaon1_dz_XLessPV/F");
  outTree->Branch("kaon2_dxy_XLessPV",    &out_kaon2_dxy_XLessPV,         "kaon2_dxy_XLessPV/F");
  outTree->Branch("kaon2_dz_XLessPV",     &out_kaon2_dz_XLessPV,  "kaon2_dz_XLessPV/F");
  outTree->Branch("kaon1_dxyE",   &out_kaon1_dxyE,        "kaon1_dxyE/F");
  outTree->Branch("kaon1_dzE",    &out_kaon1_dzE,         "kaon1_dzE/F");
  outTree->Branch("kaon2_dxyE",   &out_kaon2_dxyE,        "kaon2_dxyE/F");
  outTree->Branch("kaon2_dzE",    &out_kaon2_dzE,         "kaon2_dzE/F");
  outTree->Branch("XMuMuIdx",     &out_XMuMuIdx,  "XMuMuIdx/F");
  outTree->Branch("XKaon1Idx",    &out_XKaon1Idx,         "XKaon1Idx/F");
  outTree->Branch("XKaon2Idx",    &out_XKaon2Idx,         "XKaon2Idx/F");
  outTree->Branch("Kaon1FromPV",  &out_Kaon1FromPV,       "Kaon1FromPV/F");
  outTree->Branch("Kaon2FromPV",  &out_Kaon2FromPV,       "Kaon2FromPV/F");
  outTree->Branch("Muon1Px_MuMuKK",       &out_Muon1Px_MuMuKK,    "Muon1Px_MuMuKK/F");
  outTree->Branch("Muon1Py_MuMuKK",       &out_Muon1Py_MuMuKK,    "Muon1Py_MuMuKK/F");
  outTree->Branch("Muon1Pz_MuMuKK",       &out_Muon1Pz_MuMuKK,    "Muon1Pz_MuMuKK/F");
  outTree->Branch("Muon1E_MuMuKK",        &out_Muon1E_MuMuKK,     "Muon1E_MuMuKK/F");
  outTree->Branch("Muon2Px_MuMuKK",       &out_Muon2Px_MuMuKK,    "Muon2Px_MuMuKK/F");
  outTree->Branch("Muon2Py_MuMuKK",       &out_Muon2Py_MuMuKK,    "Muon2Py_MuMuKK/F");
  outTree->Branch("Muon2Pz_MuMuKK",       &out_Muon2Pz_MuMuKK,    "Muon2Pz_MuMuKK/F");
  outTree->Branch("Muon2E_MuMuKK",        &out_Muon2E_MuMuKK,     "Muon2E_MuMuKK/F");
  outTree->Branch("Kaon1Px_MuMuKK",       &out_Kaon1Px_MuMuKK,    "Kaon1Px_MuMuKK/F");
  outTree->Branch("Kaon1Py_MuMuKK",       &out_Kaon1Py_MuMuKK,    "Kaon1Py_MuMuKK/F");
  outTree->Branch("Kaon1Pz_MuMuKK",       &out_Kaon1Pz_MuMuKK,    "Kaon1Pz_MuMuKK/F");
  outTree->Branch("Kaon1E_MuMuKK",        &out_Kaon1E_MuMuKK,     "Kaon1E_MuMuKK/F");
  outTree->Branch("kaon1_nsigdedx",       &out_kaon1_nsigdedx,    "kaon1_nsigdedx/F");
  outTree->Branch("kaon1_dedx",   &out_kaon1_dedx,        "kaon1_dedx/F");
  outTree->Branch("kaon1_dedxMass",       &out_kaon1_dedxMass,    "kaon1_dedxMass/F");
  outTree->Branch("kaon1_theo",   &out_kaon1_theo,        "kaon1_theo/F");
  outTree->Branch("kaon1_sigma",  &out_kaon1_sigma,       "kaon1_sigma/F");
  outTree->Branch("kaon1_dedx_byHits",    &out_kaon1_dedx_byHits,         "kaon1_dedx_byHits/F");
  outTree->Branch("kaon1_dedxErr_byHits",         &out_kaon1_dedxErr_byHits,      "kaon1_dedxErr_byHits/F");
  outTree->Branch("kaon1_saturMeas_byHits",       &out_kaon1_saturMeas_byHits,    "kaon1_saturMeas_byHits/F");
  outTree->Branch("kaon1_Meas_byHits",    &out_kaon1_Meas_byHits,         "kaon1_Meas_byHits/F");
  outTree->Branch("Kaon2Px_MuMuKK",       &out_Kaon2Px_MuMuKK,    "Kaon2Px_MuMuKK/F");
  outTree->Branch("Kaon2Py_MuMuKK",       &out_Kaon2Py_MuMuKK,    "Kaon2Py_MuMuKK/F");
  outTree->Branch("Kaon2Pz_MuMuKK",       &out_Kaon2Pz_MuMuKK,    "Kaon2Pz_MuMuKK/F");
  outTree->Branch("Kaon2E_MuMuKK",        &out_Kaon2E_MuMuKK,     "Kaon2E_MuMuKK/F");
  outTree->Branch("kaon2_nsigdedx",       &out_kaon2_nsigdedx,    "kaon2_nsigdedx/F");
  outTree->Branch("kaon2_dedx",   &out_kaon2_dedx,        "kaon2_dedx/F");
  outTree->Branch("kaon2_dedxMass",       &out_kaon2_dedxMass,    "kaon2_dedxMass/F");
  outTree->Branch("kaon2_theo",   &out_kaon2_theo,        "kaon2_theo/F");
  outTree->Branch("kaon2_sigma",  &out_kaon2_sigma,       "kaon2_sigma/F");
  outTree->Branch("kaon2_dedx_byHits",    &out_kaon2_dedx_byHits,         "kaon2_dedx_byHits/F");
  outTree->Branch("kaon2_dedxErr_byHits",         &out_kaon2_dedxErr_byHits,      "kaon2_dedxErr_byHits/F");
  outTree->Branch("kaon2_saturMeas_byHits",       &out_kaon2_saturMeas_byHits,    "kaon2_saturMeas_byHits/F");
  outTree->Branch("kaon2_Meas_byHits",    &out_kaon2_Meas_byHits,         "kaon2_Meas_byHits/F");

}

Bool_t 2mu2k_2012::Process(Long64_t entry)
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

   ////////////////// Histograms //////////////////
   JPsi_mass = 3.096916; /// pdg mass
   Phi_mass = 1.019455;

   bool HLT_4_v9 = false, HLT_4_v10 = false, HLT_4_v11 = false, HLT_4_v12 = false ;
   bool HLT_8_v3 = false, HLT_8_v4 = false, HLT_8_v5 = false, HLT_8_v6 = false, HLT_8_v7 = false ;
   bool HLT_4_vAny = false;
   bool HLT_8_vAny = false;
   bool HLT_Any = false;

   std::vector<bool> hltsFlags;

  for (Int_t i = 0; i != abs((int)(TrigRes->size())); ++i)
  {
    if ( TrigNames->at(i).find("HLT_DoubleMu4_Jpsi_Displaced_v9") != string::npos  &&  TrigRes->at(i) == 1) HLT_4_v9 = true;
    if ( TrigNames->at(i).find("HLT_DoubleMu4_Jpsi_Displaced_v10") != string::npos  &&  TrigRes->at(i) == 1) HLT_4_v10 = true;
    if ( TrigNames->at(i).find("HLT_DoubleMu4_Jpsi_Displaced_v11") != string::npos  &&  TrigRes->at(i) == 1) HLT_4_v11 = true;
    if ( TrigNames->at(i).find("HLT_DoubleMu4_Jpsi_Displaced_v12") != string::npos  &&  TrigRes->at(i) == 1) HLT_4_v12 = true;

    if ( TrigNames->at(i).find("HLT_Dimuon8_Jpsi_v3") != string::npos  &&  TrigRes->at(i) == 1) HLT_8_v3 = true;
    if ( TrigNames->at(i).find("HLT_Dimuon8_Jpsi_v4") != string::npos  &&  TrigRes->at(i) == 1) HLT_8_v4 = true;
    if ( TrigNames->at(i).find("HLT_Dimuon8_Jpsi_v5") != string::npos  &&  TrigRes->at(i) == 1) HLT_8_v5 = true;
    if ( TrigNames->at(i).find("HLT_Dimuon8_Jpsi_v6") != string::npos  &&  TrigRes->at(i) == 1) HLT_8_v6 = true;
    if ( TrigNames->at(i).find("HLT_Dimuon8_Jpsi_v7") != string::npos  &&  TrigRes->at(i) == 1) HLT_8_v7 = true;
  }

  if (HLT_4_v9 || HLT_4_v10 || HLT_4_v11 || HLT_4_v12) HLT_4_vAny = true;
  if (HLT_8_v3 || HLT_8_v4 || HLT_8_v5 || HLT_8_v6 || HLT_8_v7) HLT_8_vAny = true;
  if (HLT_8_vAny || HLT_4_vAny) HLT_Any = true;

  if(HLT_Any)
  for(int iX=0; iX<(*nX); ++iX)
  {
    int iJPsi = (*XMuMuIdx)[iX];
    int iMu1 = (*mu1Idx)[iJPsi] ; // define for original muon1
    int iMu2 = (*mu2Idx)[iJPsi] ; // define for original muon2
    int iK1 = (*ka1Idx)[iX] ; // define for original kaon1
    int iK2 = (*ka2Idx)[iX] ;

    double mu1_E = 0., mu2_E = 0., K1_E = 0., K2_E = 0.;

    TLorentzVector mu1, mu2, oMu1, oMu2;

    mu1.SetPxPyPzE((*Muon1Px_MuMuKK)[iX],(*Muon1Py_MuMuKK)[iX],(*Muon1Pz_MuMuKK)[iX],(*Muon1E_MuMuKK)[iX]);
    mu2.SetPxPyPzE((*Muon2Px_MuMuKK)[iX],(*Muon2Py_MuMuKK)[iX],(*Muon2Pz_MuMuKK)[iX],(*Muon2E_MuMuKK)[iX]);

    mu1_E = sqrt( pow((*muPx)[iMu1], 2) + pow((*muPy)[iMu1], 2) + pow((*muPz)[iMu1], 2) + pow(muon_mass, 2) ) ;
    mu2_E = sqrt( pow((*muPx)[iMu2], 2) + pow((*muPy)[iMu2], 2) + pow((*muPz)[iMu2], 2) + pow(muon_mass, 2) ) ;
    oMu1.SetPxPyPzE( (*muPx)[iMu1], (*muPy)[iMu1], (*muPz)[iMu1], mu1_E) ;
    oMu2.SetPxPyPzE( (*muPx)[iMu2], (*muPy)[iMu2], (*muPz)[iMu2], mu2_E) ;

    TLorentzVector JPsi;
    JPsi = mu1 + mu2;

    TLorentzVector JPsiOriginal;
    JPsiOriginal = oMu1 + oMu2;

    TLorentzVector kaon1,kaon2;

    K1_E=sqrt(pow((*Kaon1Px_MuMuKK)[iX],2)+pow((*Kaon1Py_MuMuKK)[iX],2)+pow((*Kaon1Pz_MuMuKK)[iX],2)+pow(kaonCh_mass,2));
    kaon1.SetPxPyPzE((*Kaon1Px_MuMuKK)[iX],(*Kaon1Py_MuMuKK)[iX],(*Kaon1Pz_MuMuKK)[iX],K1_E);
    K2_E=sqrt(pow((*Kaon2Px_MuMuKK)[iX],2)+pow((*Kaon2Py_MuMuKK)[iX],2)+pow((*Kaon2Pz_MuMuKK)[iX],2)+pow(kaonCh_mass,2));
    kaon2.SetPxPyPzE((*Kaon2Px_MuMuKK)[iX],(*Kaon2Py_MuMuKK)[iX],(*Kaon2Pz_MuMuKK)[iX],K2_E);

    TLorentzVector Phi;
    Phi = kaon1 + kaon2;

    out_TrigRes =   (Float_t)(*TrigRes);
    out_TrigNames =         (Float_t)(*TrigNames);
    out_MatchTriggerNames =         (Float_t)(*MatchTriggerNames);
    out_L1TrigRes =         (Float_t)(*L1TrigRes);
    out_evtNum =    (Float_t)(*evtNum);
    out_runNum =    (Float_t)(*runNum);
    out_lumiNum =   (Float_t)(*lumiNum);
    out_priVtx_n =  (Float_t)(*priVtx_n);
    out_priVtx_X =  (Float_t)(*priVtx_X);
    out_priVtx_Y =  (Float_t)(*priVtx_Y);
    out_priVtx_Z =  (Float_t)(*priVtx_Z);
    out_priVtx_XE =         (Float_t)(*priVtx_XE);
    out_priVtx_YE =         (Float_t)(*priVtx_YE);
    out_priVtx_ZE =         (Float_t)(*priVtx_ZE);
    out_priVtx_NormChi2 =   (Float_t)(*priVtx_NormChi2);
    out_priVtx_Chi2 =       (Float_t)(*priVtx_Chi2);
    out_priVtx_CL =         (Float_t)(*priVtx_CL);
    out_priVtx_tracks =     (Float_t)(*priVtx_tracks);
    out_priVtx_tracksPtSq =         (Float_t)(*priVtx_tracksPtSq);

    out_XMass =     (Float_t)(*XMass);
    out_JpsiMass = (Float_t)(JPsiOriginal.M());
    out_JpsiMass_ref = (Float_t)(JPsi.M());
    out_PhiMass = (Float_t)(Phi.M());
    out_XPx =       (Float_t)(*XPx);
    out_XPy =       (Float_t)(*XPy);
    out_XPz =       (Float_t)(*XPz);
    out_XPxE =      (Float_t)(*XPxE);
    out_XPyE =      (Float_t)(*XPyE);
    out_XPzE =      (Float_t)(*XPzE);
    out_XVtx_CL =   (Float_t)(*XVtx_CL);
    out_XVtx_Chi2 =         (Float_t)(*XVtx_Chi2);
    out_XDecayVtx_X =       (Float_t)(*XDecayVtx_X);
    out_XDecayVtx_Y =       (Float_t)(*XDecayVtx_Y);
    out_XDecayVtx_Z =       (Float_t)(*XDecayVtx_Z);
    out_XDecayVtx_XE =      (Float_t)(*XDecayVtx_XE);
    out_XDecayVtx_YE =      (Float_t)(*XDecayVtx_YE);
    out_XDecayVtx_ZE =      (Float_t)(*XDecayVtx_ZE);
    out_XCosAlphaBS =       (Float_t)(*XCosAlphaBS);
    out_XCosAlpha3DBS =     (Float_t)(*XCosAlpha3DBS);
    out_XCTauBS =   (Float_t)(*XCTauBS);
    out_XCTauBSE =  (Float_t)(*XCTauBSE);
    out_XLxyBS =    (Float_t)(*XLxyBS);
    out_XLxyBSE =   (Float_t)(*XLxyBSE);
    out_XLxyzBS =   (Float_t)(*XLxyzBS);
    out_XLxyzBSE =  (Float_t)(*XLxyzBSE);
    out_XCosAlphaPV =       (Float_t)(*XCosAlphaPV);
    out_XCosAlpha3DPV =     (Float_t)(*XCosAlpha3DPV);
    out_XCTauPV =   (Float_t)(*XCTauPV);
    out_XCTauPVE =  (Float_t)(*XCTauPVE);
    out_XLxyPV =    (Float_t)(*XLxyPV);
    out_XLxyPVE =   (Float_t)(*XLxyPVE);
    out_XLxyzPV =   (Float_t)(*XLxyzPV);
    out_XLxyzPVE =  (Float_t)(*XLxyzPVE);
    out_PriVtx_XCosAlpha_n =        (Float_t)(*PriVtx_XCosAlpha_n);
    out_PriVtx_XCosAlpha_X =        (Float_t)(*PriVtx_XCosAlpha_X);
    out_PriVtx_XCosAlpha_Y =        (Float_t)(*PriVtx_XCosAlpha_Y);
    out_PriVtx_XCosAlpha_Z =        (Float_t)(*PriVtx_XCosAlpha_Z);
    out_PriVtx_XCosAlpha_EX =       (Float_t)(*PriVtx_XCosAlpha_EX);
    out_PriVtx_XCosAlpha_EY =       (Float_t)(*PriVtx_XCosAlpha_EY);
    out_PriVtx_XCosAlpha_EZ =       (Float_t)(*PriVtx_XCosAlpha_EZ);
    out_PriVtx_XCosAlpha_Chi2 =     (Float_t)(*PriVtx_XCosAlpha_Chi2);
    out_PriVtx_XCosAlpha_CL =       (Float_t)(*PriVtx_XCosAlpha_CL);
    out_PriVtx_XCosAlpha_tracks =   (Float_t)(*PriVtx_XCosAlpha_tracks);
    out_XCosAlphaPVCosAlpha =       (Float_t)(*XCosAlphaPVCosAlpha);
    out_XCosAlpha3DPVCosAlpha =     (Float_t)(*XCosAlpha3DPVCosAlpha);
    out_XCTauPVCosAlpha =   (Float_t)(*XCTauPVCosAlpha);
    out_XCTauPVCosAlphaE =  (Float_t)(*XCTauPVCosAlphaE);
    out_XLxyPVCosAlpha =    (Float_t)(*XLxyPVCosAlpha);
    out_XLxyPVCosAlphaE =   (Float_t)(*XLxyPVCosAlphaE);
    out_XLxyzPVCosAlpha =   (Float_t)(*XLxyzPVCosAlpha);
    out_XLxyzPVCosAlphaE =  (Float_t)(*XLxyzPVCosAlphaE);
    out_PriVtx_XCosAlpha3D_n =      (Float_t)(*PriVtx_XCosAlpha3D_n);
    out_PriVtx_XCosAlpha3D_X =      (Float_t)(*PriVtx_XCosAlpha3D_X);
    out_PriVtx_XCosAlpha3D_Y =      (Float_t)(*PriVtx_XCosAlpha3D_Y);
    out_PriVtx_XCosAlpha3D_Z =      (Float_t)(*PriVtx_XCosAlpha3D_Z);
    out_PriVtx_XCosAlpha3D_EX =     (Float_t)(*PriVtx_XCosAlpha3D_EX);
    out_PriVtx_XCosAlpha3D_EY =     (Float_t)(*PriVtx_XCosAlpha3D_EY);
    out_PriVtx_XCosAlpha3D_EZ =     (Float_t)(*PriVtx_XCosAlpha3D_EZ);
    out_PriVtx_XCosAlpha3D_Chi2 =   (Float_t)(*PriVtx_XCosAlpha3D_Chi2);
    out_PriVtx_XCosAlpha3D_CL =     (Float_t)(*PriVtx_XCosAlpha3D_CL);
    out_PriVtx_XCosAlpha3D_tracks =         (Float_t)(*PriVtx_XCosAlpha3D_tracks);
    out_XCosAlphaPVCosAlpha3D =     (Float_t)(*XCosAlphaPVCosAlpha3D);
    out_XCosAlpha3DPVCosAlpha3D =   (Float_t)(*XCosAlpha3DPVCosAlpha3D);
    out_XCTauPVCosAlpha3D =         (Float_t)(*XCTauPVCosAlpha3D);
    out_XCTauPVCosAlpha3DE =        (Float_t)(*XCTauPVCosAlpha3DE);
    out_XLxyPVCosAlpha3D =  (Float_t)(*XLxyPVCosAlpha3D);
    out_XLxyPVCosAlpha3DE =         (Float_t)(*XLxyPVCosAlpha3DE);
    out_XLxyzPVCosAlpha3D =         (Float_t)(*XLxyzPVCosAlpha3D);
    out_XLxyzPVCosAlpha3DE =        (Float_t)(*XLxyzPVCosAlpha3DE);
    out_XLessPV_tracksPtSq =        (Float_t)(*XLessPV_tracksPtSq);
    out_XLessPV_4tracksPtSq =       (Float_t)(*XLessPV_4tracksPtSq);
    out_PriVtxXLess_n =     (Float_t)(*PriVtxXLess_n);
    out_PriVtxXLess_X =     (Float_t)(*PriVtxXLess_X);
    out_PriVtxXLess_Y =     (Float_t)(*PriVtxXLess_Y);
    out_PriVtxXLess_Z =     (Float_t)(*PriVtxXLess_Z);
    out_PriVtxXLess_EX =    (Float_t)(*PriVtxXLess_EX);
    out_PriVtxXLess_EY =    (Float_t)(*PriVtxXLess_EY);
    out_PriVtxXLess_EZ =    (Float_t)(*PriVtxXLess_EZ);
    out_PriVtxXLess_Chi2 =  (Float_t)(*PriVtxXLess_Chi2);
    out_PriVtxXLess_CL =    (Float_t)(*PriVtxXLess_CL);
    out_PriVtxXLess_tracks =        (Float_t)(*PriVtxXLess_tracks);
    out_XCosAlphaXLessPV =  (Float_t)(*XCosAlphaXLessPV);
    out_XCosAlpha3DXLessPV =        (Float_t)(*XCosAlpha3DXLessPV);
    out_XCTauXLessPV =      (Float_t)(*XCTauXLessPV);
    out_XCTauXLessPVE =     (Float_t)(*XCTauXLessPVE);
    out_XLxyXLessPV =       (Float_t)(*XLxyXLessPV);
    out_XLxyXLessPVE =      (Float_t)(*XLxyXLessPVE);
    out_XLxyzXLessPV =      (Float_t)(*XLxyzXLessPV);
    out_XLxyzXLessPVE =     (Float_t)(*XLxyzXLessPVE);
    out_PriVtxXLess_XCosAlpha_n =   (Float_t)(*PriVtxXLess_XCosAlpha_n);
    out_PriVtxXLess_XCosAlpha_X =   (Float_t)(*PriVtxXLess_XCosAlpha_X);
    out_PriVtxXLess_XCosAlpha_Y =   (Float_t)(*PriVtxXLess_XCosAlpha_Y);
    out_PriVtxXLess_XCosAlpha_Z =   (Float_t)(*PriVtxXLess_XCosAlpha_Z);
    out_PriVtxXLess_XCosAlpha_EX =  (Float_t)(*PriVtxXLess_XCosAlpha_EX);
    out_PriVtxXLess_XCosAlpha_EY =  (Float_t)(*PriVtxXLess_XCosAlpha_EY);
    out_PriVtxXLess_XCosAlpha_EZ =  (Float_t)(*PriVtxXLess_XCosAlpha_EZ);
    out_PriVtxXLess_XCosAlpha_Chi2 =        (Float_t)(*PriVtxXLess_XCosAlpha_Chi2);
    out_PriVtxXLess_XCosAlpha_CL =  (Float_t)(*PriVtxXLess_XCosAlpha_CL);
    out_PriVtxXLess_XCosAlpha_tracks =      (Float_t)(*PriVtxXLess_XCosAlpha_tracks);
    out_XCosAlphaXLessPVCosAlpha =  (Float_t)(*XCosAlphaXLessPVCosAlpha);
    out_XCosAlpha3DXLessPVCosAlpha =        (Float_t)(*XCosAlpha3DXLessPVCosAlpha);
    out_XCTauXLessPVCosAlpha =      (Float_t)(*XCTauXLessPVCosAlpha);
    out_XCTauXLessPVCosAlphaE =     (Float_t)(*XCTauXLessPVCosAlphaE);
    out_XLxyXLessPVCosAlpha =       (Float_t)(*XLxyXLessPVCosAlpha);
    out_XLxyXLessPVCosAlphaE =      (Float_t)(*XLxyXLessPVCosAlphaE);
    out_XLxyzXLessPVCosAlpha =      (Float_t)(*XLxyzXLessPVCosAlpha);
    out_XLxyzXLessPVCosAlphaE =     (Float_t)(*XLxyzXLessPVCosAlphaE);
    out_PriVtxXLess_XCosAlpha3D_n =         (Float_t)(*PriVtxXLess_XCosAlpha3D_n);
    out_PriVtxXLess_XCosAlpha3D_X =         (Float_t)(*PriVtxXLess_XCosAlpha3D_X);
    out_PriVtxXLess_XCosAlpha3D_Y =         (Float_t)(*PriVtxXLess_XCosAlpha3D_Y);
    out_PriVtxXLess_XCosAlpha3D_Z =         (Float_t)(*PriVtxXLess_XCosAlpha3D_Z);
    out_PriVtxXLess_XCosAlpha3D_EX =        (Float_t)(*PriVtxXLess_XCosAlpha3D_EX);
    out_PriVtxXLess_XCosAlpha3D_EY =        (Float_t)(*PriVtxXLess_XCosAlpha3D_EY);
    out_PriVtxXLess_XCosAlpha3D_EZ =        (Float_t)(*PriVtxXLess_XCosAlpha3D_EZ);
    out_PriVtxXLess_XCosAlpha3D_Chi2 =      (Float_t)(*PriVtxXLess_XCosAlpha3D_Chi2);
    out_PriVtxXLess_XCosAlpha3D_CL =        (Float_t)(*PriVtxXLess_XCosAlpha3D_CL);
    out_PriVtxXLess_XCosAlpha3D_tracks =    (Float_t)(*PriVtxXLess_XCosAlpha3D_tracks);
    out_XCosAlphaXLessPVCosAlpha3D =        (Float_t)(*XCosAlphaXLessPVCosAlpha3D);
    out_XCosAlpha3DXLessPVCosAlpha3D =      (Float_t)(*XCosAlpha3DXLessPVCosAlpha3D);
    out_XCTauXLessPVCosAlpha3D =    (Float_t)(*XCTauXLessPVCosAlpha3D);
    out_XCTauXLessPVCosAlpha3DE =   (Float_t)(*XCTauXLessPVCosAlpha3DE);
    out_XLxyXLessPVCosAlpha3D =     (Float_t)(*XLxyXLessPVCosAlpha3D);
    out_XLxyXLessPVCosAlpha3DE =    (Float_t)(*XLxyXLessPVCosAlpha3DE);
    out_XLxyzXLessPVCosAlpha3D =    (Float_t)(*XLxyzXLessPVCosAlpha3D);
    out_XLxyzXLessPVCosAlpha3DE =   (Float_t)(*XLxyzXLessPVCosAlpha3DE);
    out_PriVtxXCorr_n =     (Float_t)(*PriVtxXCorr_n);
    out_PriVtxXCorr_X =     (Float_t)(*PriVtxXCorr_X);
    out_PriVtxXCorr_Y =     (Float_t)(*PriVtxXCorr_Y);
    out_PriVtxXCorr_Z =     (Float_t)(*PriVtxXCorr_Z);
    out_PriVtxXCorr_EX =    (Float_t)(*PriVtxXCorr_EX);
    out_PriVtxXCorr_EY =    (Float_t)(*PriVtxXCorr_EY);
    out_PriVtxXCorr_EZ =    (Float_t)(*PriVtxXCorr_EZ);
    out_PriVtxXCorr_Chi2 =  (Float_t)(*PriVtxXCorr_Chi2);
    out_PriVtxXCorr_CL =    (Float_t)(*PriVtxXCorr_CL);
    out_PriVtxXCorr_tracks =        (Float_t)(*PriVtxXCorr_tracks);
    out_XCosAlphaPVX =      (Float_t)(*XCosAlphaPVX);
    out_XCTauPVX =  (Float_t)(*XCTauPVX);
    out_XCTauPVXE =         (Float_t)(*XCTauPVXE);
    out_XLxyPVX =   (Float_t)(*XLxyPVX);
    out_XLxyzPVX =  (Float_t)(*XLxyzPVX);
    out_XCTauPVX_3D =       (Float_t)(*XCTauPVX_3D);
    out_XCTauPVX_3D_err =   (Float_t)(*XCTauPVX_3D_err);
    out_kaon1_dxy_PV =      (Float_t)(*kaon1_dxy_PV);
    out_kaon1_dz_PV =       (Float_t)(*kaon1_dz_PV);
    out_kaon2_dxy_PV =      (Float_t)(*kaon2_dxy_PV);
    out_kaon2_dz_PV =       (Float_t)(*kaon2_dz_PV);
    out_kaon1_dxy_BS =      (Float_t)(*kaon1_dxy_BS);
    out_kaon1_dz_BS =       (Float_t)(*kaon1_dz_BS);
    out_kaon2_dxy_BS =      (Float_t)(*kaon2_dxy_BS);
    out_kaon2_dz_BS =       (Float_t)(*kaon2_dz_BS);
    out_kaon1_dxy_XLessPV =         (Float_t)(*kaon1_dxy_XLessPV);
    out_kaon1_dz_XLessPV =  (Float_t)(*kaon1_dz_XLessPV);
    out_kaon2_dxy_XLessPV =         (Float_t)(*kaon2_dxy_XLessPV);
    out_kaon2_dz_XLessPV =  (Float_t)(*kaon2_dz_XLessPV);
    out_kaon1_dxyE =        (Float_t)(*kaon1_dxyE);
    out_kaon1_dzE =         (Float_t)(*kaon1_dzE);
    out_kaon2_dxyE =        (Float_t)(*kaon2_dxyE);
    out_kaon2_dzE =         (Float_t)(*kaon2_dzE);
    out_XMuMuIdx =  (Float_t)(*XMuMuIdx);
    out_XKaon1Idx =         (Float_t)(*XKaon1Idx);
    out_XKaon2Idx =         (Float_t)(*XKaon2Idx);
    out_Kaon1FromPV =       (Float_t)(*Kaon1FromPV);
    out_Kaon2FromPV =       (Float_t)(*Kaon2FromPV);
    out_Muon1Px_MuMuKK =    (Float_t)(*Muon1Px_MuMuKK);
    out_Muon1Py_MuMuKK =    (Float_t)(*Muon1Py_MuMuKK);
    out_Muon1Pz_MuMuKK =    (Float_t)(*Muon1Pz_MuMuKK);
    out_Muon1E_MuMuKK =     (Float_t)(*Muon1E_MuMuKK);
    out_Muon2Px_MuMuKK =    (Float_t)(*Muon2Px_MuMuKK);
    out_Muon2Py_MuMuKK =    (Float_t)(*Muon2Py_MuMuKK);
    out_Muon2Pz_MuMuKK =    (Float_t)(*Muon2Pz_MuMuKK);
    out_Muon2E_MuMuKK =     (Float_t)(*Muon2E_MuMuKK);
    out_Kaon1Px_MuMuKK =    (Float_t)(*Kaon1Px_MuMuKK);
    out_Kaon1Py_MuMuKK =    (Float_t)(*Kaon1Py_MuMuKK);
    out_Kaon1Pz_MuMuKK =    (Float_t)(*Kaon1Pz_MuMuKK);
    out_Kaon1E_MuMuKK =     (Float_t)(*Kaon1E_MuMuKK);
    out_kaon1_nsigdedx =    (Float_t)(*kaon1_nsigdedx);
    out_kaon1_dedx =        (Float_t)(*kaon1_dedx);
    out_kaon1_dedxMass =    (Float_t)(*kaon1_dedxMass);
    out_kaon1_theo =        (Float_t)(*kaon1_theo);
    out_kaon1_sigma =       (Float_t)(*kaon1_sigma);
    out_kaon1_dedx_byHits =         (Float_t)(*kaon1_dedx_byHits);
    out_kaon1_dedxErr_byHits =      (Float_t)(*kaon1_dedxErr_byHits);
    out_kaon1_saturMeas_byHits =    (Float_t)(*kaon1_saturMeas_byHits);
    out_kaon1_Meas_byHits =         (Float_t)(*kaon1_Meas_byHits);
    out_Kaon2Px_MuMuKK =    (Float_t)(*Kaon2Px_MuMuKK);
    out_Kaon2Py_MuMuKK =    (Float_t)(*Kaon2Py_MuMuKK);
    out_Kaon2Pz_MuMuKK =    (Float_t)(*Kaon2Pz_MuMuKK);
    out_Kaon2E_MuMuKK =     (Float_t)(*Kaon2E_MuMuKK);
    out_kaon2_nsigdedx =    (Float_t)(*kaon2_nsigdedx);
    out_kaon2_dedx =        (Float_t)(*kaon2_dedx);
    out_kaon2_dedxMass =    (Float_t)(*kaon2_dedxMass);
    out_kaon2_theo =        (Float_t)(*kaon2_theo);
    out_kaon2_sigma =       (Float_t)(*kaon2_sigma);
    out_kaon2_dedx_byHits =         (Float_t)(*kaon2_dedx_byHits);
    out_kaon2_dedxErr_byHits =      (Float_t)(*kaon2_dedxErr_byHits);
    out_kaon2_saturMeas_byHits =    (Float_t)(*kaon2_saturMeas_byHits);
    out_kaon2_Meas_byHits =         (Float_t)(*kaon2_Meas_byHits);


  }

   return kTRUE;
}

void 2mu2k_2012::SlaveTerminate()
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

void 2mu2k_2012::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}
