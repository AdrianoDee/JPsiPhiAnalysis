{

#include "TTree.h"
#include "TDSet.h"
#include "TProof.h"
#include "TString.h"


  // INPUT DATA SAMPLE ON LOCAL DISK

  TDSet* dataset = new TDSet("TTree", "FourMuonTree", "rootuple");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018A-PromptReco-v1_MINIAOD___20180919_154605/180919_134618/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018A-PromptReco-v1_MINIAOD___20180919_154605/180919_134618/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018A-PromptReco-v1_MINIAOD___20180919_154605/180919_134618/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018A-PromptReco-v1_MINIAOD___20180919_154605/180919_134618/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018A-PromptReco-v1_MINIAOD___20180919_154605/180919_134618/0004/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018A-PromptReco-v1_MINIAOD___20180919_154605/180919_134618/0005/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018A-PromptReco-v1_MINIAOD___20180919_154605/180919_134618/0006/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018A-PromptReco-v2_MINIAOD___20180919_153708/180919_133721/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018A-PromptReco-v2_MINIAOD___20180919_153708/180919_133721/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018A-PromptReco-v3_MINIAOD___20180919_154229/180919_134242/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018A-PromptReco-v3_MINIAOD___20180919_154229/180919_134242/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018A-PromptReco-v3_MINIAOD___20180919_154229/180919_134242/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018A-PromptReco-v3_MINIAOD___20180919_154229/180919_134242/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018B-PromptReco-v1_MINIAOD___20180919_154322/180919_134334/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018B-PromptReco-v1_MINIAOD___20180919_154322/180919_134334/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018B-PromptReco-v1_MINIAOD___20180919_154322/180919_134334/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018B-PromptReco-v1_MINIAOD___20180919_154322/180919_134334/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018B-PromptReco-v1_MINIAOD___20180919_154322/180919_134334/0004/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018C-PromptReco-v1_MINIAOD___20180919_154021/180919_134032/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018C-PromptReco-v1_MINIAOD___20180919_154138/180919_134150/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018C-PromptReco-v2_MINIAOD___20180919_154359/180919_134409/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018C-PromptReco-v2_MINIAOD___20180919_154359/180919_134409/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018C-PromptReco-v2_MINIAOD___20180919_154359/180919_134409/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018C-PromptReco-v3_MINIAOD___20180919_154445/180919_134500/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018C-PromptReco-v3_MINIAOD___20180919_154445/180919_134500/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018C-PromptReco-v3_MINIAOD___20180919_154445/180919_134500/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20180919_163853/180919_143902/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20180919_163853/180919_143902/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20180919_163853/180919_143902/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20180919_163853/180919_143902/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20180919_163853/180919_143902/0004/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20180919_163853/180919_143902/0005/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20180919_163853/180919_143902/0006/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20180919_163853/180919_143902/0007/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20180919_163853/180919_143902/0008/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2mu_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20180919_163853/180919_143902/0009/sum.root");


  TString selector = "/lustre/home/adrianodif/jpsiphi/2018/data_2018/analysis/utilities/skimmers/2mu2mu/FourMuons";
  TProof *p = TProof::Open("workers=5"); // 12 workers for qsub
  //gProofDebugMask = TProofDebug::kAll;
  //gProofDebugLevel = 5;

  // Processing
  cout << ">> Processing " << selector << " ... " << endl;

  TString selectorplus = selector;
  selectorplus += ".C+";
  p->Process(dataset, selectorplus);

}
