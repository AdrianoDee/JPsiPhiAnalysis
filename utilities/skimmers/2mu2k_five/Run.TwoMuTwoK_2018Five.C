{

#include "TTree.h"
#include "TDSet.h"
#include "TProof.h"
#include "TString.h"


  // INPUT DATA SAMPLE ON LOCAL DISK

  TDSet* dataset = new TDSet("TTree", "JPsiPhiTree", "rootuple");
  //
  //dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/2mu2k_miniaod_17Nov2017_BCDEF_2017.root");
  //dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/2mu2k_miniaod_07Aug17_BCDEFGH_2016.root");
  //phi mass window [0.97-1.06]
  //dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/2018/2mu2k-Run2018A-PromptReco-v3.root");
  //dataset->Add("/lustre/cms/store/user/adiflori/Charmonium//2mu2k_miniaod_17Nov2017_BCDEF_2017_phi_097_106.root");
  dataset->Add("/lustre/home/adrianodif/jpsiphi/2018/CMSSW_10_2_1/src/jpsiphi/jpsiphi/test/rootuple-2018-dimuonditrak_bbbar_hard_0.root");

  TString selector = "/lustre/home/adrianodif/jpsiphi/analysis/utilities/skimmers/2mu2k_five/TwoMuTwoK_Five";
  TProof *p = TProof::Open("workers=1"); // 12 workers for qsub

  // Processing
  cout << ">> Processing " << selector << " ... " << endl;

  TString selectorplus = selector;
  selectorplus += ".C+";
  p->Process(dataset, selectorplus);

}
