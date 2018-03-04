{

#include "TTree.h"
#include "TDSet.h"
#include "TProof.h"
#include "TString.h"


  // INPUT DATA SAMPLE ON LOCAL DISK

  TDSet* dataset = new TDSet("TTree", "ditrakTree", "DiTrakRootupler");
  //
  //dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2k2Trig_Charmonium_Run2017F-17Nov2017-v1_MINIAOD_305388-309000__20180227_203341/merge.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2m2k2Trig_Charmonium_Run2017F-17Nov2017-v1_MINIAOD_305388-309000__20180303_205340/merge.root");
  // dataset->Add("/Users/adrianodiflorio/Documents/Git/X4140/ProofLite/Y4140_testrootuple.root");
  TString selector = "DiTrakSkim";
  TProof *p = TProof::Open("workers=40"); // 12 workers for qsub

  // Processing
  cout << ">> Processing " << selector << " ... " << endl;

  TString selectorplus = selector;
  selectorplus += ".C+";
  p->Process(dataset, selectorplus);

}
