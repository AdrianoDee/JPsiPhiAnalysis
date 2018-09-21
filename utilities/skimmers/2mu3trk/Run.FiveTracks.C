{

#include "TTree.h"
#include "TDSet.h"
#include "TProof.h"
#include "TString.h"


  // INPUT DATA SAMPLE ON LOCAL DISK

  TDSet* dataset = new TDSet("TTree", "FiveTracksTree", "rootupleFive");
  //
  //dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/2mu2k_miniaod_17Nov2017_BCDEF_2017.root");
  //dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/2mu2k_miniaod_07Aug17_BCDEFGH_2016.root");
  //phi mass window [0.97-1.06]
  //dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/2018/2mu2k-Run2018A-PromptReco-v3.root");
  //dataset->Add("/lustre/cms/store/user/adiflori/Charmonium//2mu2k_miniaod_17Nov2017_BCDEF_2017_phi_097_106.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018C-PromptReco-v3_MINIAOD___20180913_014624/180912_234632/0000/rootuple-2018-dimuonditrak_fivedataC2018_3_0_207.root");

  TString selector = "/lustre/home/adrianodif/jpsiphi/2018/data_2018/analysis/utilities/skimmers/2mu2k/FiveTracks";
  TProof *p = TProof::Open("workers=5"); // 12 workers for qsub
  //gProofDebugMask = TProofDebug::kAll;
  //gProofDebugLevel = 5;

  // Processing
  cout << ">> Processing " << selector << " ... " << endl;

  TString selectorplus = selector;
  selectorplus += ".C+";
  p->Process(dataset, selectorplus);

}
