{

#include "TTree.h"
#include "TDSet.h"
#include "TProof.h"
#include "TString.h"


  // INPUT DATA SAMPLE ON LOCAL DISK

  TDSet* dataset = new TDSet("TTree", "dimuonTree", "rootupleMuMu");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018A-17Sep2018-v1_MINIAOD___20190627_124818_six_five/190627_104825/0000.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018A-17Sep2018-v1_MINIAOD___20190627_124818_six_five/190627_104825/0001.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018A-17Sep2018-v1_MINIAOD___20190627_124818_six_five/190627_104825/0002.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018A-17Sep2018-v1_MINIAOD___20190627_124818_six_five/190627_104825/0003.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018A-17Sep2018-v1_MINIAOD___20190627_124818_six_five/190627_104825/0004.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018B-17Sep2018-v1_MINIAOD___20190627_112839_six_five/190627_092845/0000.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018B-17Sep2018-v1_MINIAOD___20190627_112839_six_five/190627_092845/0001.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018B-17Sep2018-v1_MINIAOD___20190627_112839_six_five/190627_092845/0002.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018B-17Sep2018-v1_MINIAOD___20190627_112839_six_five/190627_092845/0003.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018B-17Sep2018-v1_MINIAOD___20190627_112839_six_five/190627_092845/0004.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018B-17Sep2018-v1_MINIAOD___20190627_112839_six_five/190627_092845/0005.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018C-17Sep2018-v1_MINIAOD___20190627_112903_six_five/190627_092908/0000.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018C-17Sep2018-v1_MINIAOD___20190627_112903_six_five/190627_092908/0001.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018C-17Sep2018-v1_MINIAOD___20190627_112903_six_five/190627_092908/0002.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018C-17Sep2018-v1_MINIAOD___20190627_112903_six_five/190627_092908/0003.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018C-17Sep2018-v1_MINIAOD___20190627_112903_six_five/190627_092908/0004.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018C-17Sep2018-v1_MINIAOD___20190627_112903_six_five/190627_092908/0005.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20190627_125057_six_five/190627_105102/0000.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20190627_125057_six_five/190627_105102/0001.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20190627_125057_six_five/190627_105102/0002.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20190627_125057_six_five/190627_105102/0003.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20190627_125057_six_five/190627_105102/0004.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20190627_125057_six_five/190627_105102/0005.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20190627_125057_six_five/190627_105102/0006.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20190627_125057_six_five/190627_105102/0007.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20190627_125057_six_five/190627_105102/0008.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20190627_125057_six_five/190627_105102/0009.root");


  TString selector = "/lustre/home/adrianodif/jpsiphi/2018/data_2018/analysis/utilities/skimmers/dimuon/DiMuon";
  TProof *p = TProof::Open("workers=5"); // 12 workers for qsub
  //gProofDebugMask = TProofDebug::kAll;
  //gProofDebugLevel = 5;

  // Processing
  cout << ">> Processing " << selector << " ... " << endl;

  TString selectorplus = selector;
  selectorplus += ".C+";
  p->Process(dataset, selectorplus);

}
