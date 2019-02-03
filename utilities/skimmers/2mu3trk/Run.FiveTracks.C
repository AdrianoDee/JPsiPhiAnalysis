{

#include "TTree.h"
#include "TDSet.h"
#include "TProof.h"
#include "TString.h"


  // INPUT DATA SAMPLE ON LOCAL DISK

  TDSet* dataset = new TDSet("TTree", "FiveTracksTree", "rootupleFive");
  /// 2018

  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v1_MINIAOD___20181124_150539/181124_140547/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v1_MINIAOD___20181124_150539/181124_140547/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v1_MINIAOD___20181124_150539/181124_140547/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v1_MINIAOD___20181124_150539/181124_140547/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v2_MINIAOD___20181124_150752/181124_140759/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v3_MINIAOD___20181124_151626/181124_141633/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v3_MINIAOD___20181124_151626/181124_141633/0001/sum.root");
  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018B-PromptReco-v1_MINIAOD___20181126_111547/181126_101555/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018B-PromptReco-v1_MINIAOD___20181126_111547/181126_101555/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018B-PromptReco-v1_MINIAOD___20181126_111547/181126_101555/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018B-PromptReco-v1_MINIAOD___20181126_111547/181126_101555/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018B-PromptReco-v1_MINIAOD___20181126_111547/181126_101555/0004/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018B-PromptReco-v2_MINIAOD___20181124_151949/181124_142002/0000/sum.root");
  */
  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018C-PromptReco-v1_MINIAOD___20181124_152043/181124_142053/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018C-PromptReco-v2_MINIAOD___20181124_152128/181124_142152/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018C-PromptReco-v2_MINIAOD___20181124_152128/181124_142152/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018C-PromptReco-v3_MINIAOD___20181124_152643/181124_142709/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018C-PromptReco-v3_MINIAOD___20181124_152643/181124_142709/0001/sum.root");
  */
  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20181124_165527/181124_155537/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20181124_165527/181124_155537/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20181124_165527/181124_155537/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20181124_165527/181124_155537/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20181124_165527/181124_155537/0004/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20181124_165527/181124_155537/0005/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20181124_165527/181124_155537/0006/sum.root");
  */

  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018C-PromptReco-v3_MINIAOD___20180913_014624/180912_234632/0000/rootuple-2018-dimuonditrak_fivedataC2018_3_0_207.root");

  TString selector = "/lustre/home/adrianodif/jpsiphi/2018/data_2018/analysis/utilities/skimmers/2mu3trk/FiveTracks";
  TProof *p = TProof::Open("workers=5"); // 12 workers for qsub
  //gProofDebugMask = TProofDebug::kAll;
  //gProofDebugLevel = 5;

  // Processing
  cout << ">> Processing " << selector << " ... " << endl;

  TString selectorplus = selector;
  selectorplus += ".C+";
  p->Process(dataset, selectorplus);

}
