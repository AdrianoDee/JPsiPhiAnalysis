{

#include "TTree.h"
#include "TDSet.h"
#include "TProof.h"
#include "TString.h"


  // INPUT DATA SAMPLE ON LOCAL DISK

  TDSet* dataset = new TDSet("TTree", "JPsiPhiTree", "rootuple");

  /// 2018
  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v1_MINIAOD___20180913_015222/180912_235229/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v1_MINIAOD___20180913_015222/180912_235229/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v1_MINIAOD___20180913_015222/180912_235229/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v1_MINIAOD___20180913_015222/180912_235229/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v1_MINIAOD___20180913_015222/180912_235229/0004/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v1_MINIAOD___20180913_015222/180912_235229/0005/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v1_MINIAOD___20180913_015222/180912_235229/0006/sum.root");


  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v2_MINIAOD___20180913_015153/180912_235201/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v2_MINIAOD___20180913_015153/180912_235201/0001/sum.root");


  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v3_MINIAOD___20180913_015127/180912_235134/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v3_MINIAOD___20180913_015127/180912_235134/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v3_MINIAOD___20180913_015127/180912_235134/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v3_MINIAOD___20180913_015127/180912_235134/0003/sum.root");


  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018B-PromptReco-v1_MINIAOD___20180913_015054/180912_235101/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018B-PromptReco-v1_MINIAOD___20180913_015054/180912_235101/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018B-PromptReco-v1_MINIAOD___20180913_015054/180912_235101/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018B-PromptReco-v1_MINIAOD___20180913_015054/180912_235101/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018B-PromptReco-v1_MINIAOD___20180913_015054/180912_235101/0004/sum.root");


  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018B-PromptReco-v2_MINIAOD___20180913_014941/180912_234949/0000/sum.root");

  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018C-PromptReco-v1_MINIAOD___20180913_014813/180912_234820/0000/sum.root");

  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018C-PromptReco-v2_MINIAOD___20180913_014734/180912_234742/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018C-PromptReco-v2_MINIAOD___20180913_014734/180912_234742/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018C-PromptReco-v2_MINIAOD___20180913_014734/180912_234742/0002/sum.root");


  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018C-PromptReco-v3_MINIAOD___20180913_014624/180912_234632/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018C-PromptReco-v3_MINIAOD___20180913_014624/180912_234632/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018C-PromptReco-v3_MINIAOD___20180913_014624/180912_234632/0002/sum.root");


  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20180913_015627/180912_235635/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20180913_015627/180912_235635/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20180913_015627/180912_235635/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20180913_015627/180912_235635/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20180913_015627/180912_235635/0004/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20180913_015627/180912_235635/0005/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20180913_015627/180912_235635/0006/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20180913_015627/180912_235635/0007/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20180913_015627/180912_235635/0008/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20180913_015627/180912_235635/0009/sum.root");
  */

  /// 2017
  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017B-31Mar2018-v1_MINIAOD___20181022_235142/181022_215230/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017B-31Mar2018-v1_MINIAOD___20181022_235142/181022_215230/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017B-31Mar2018-v1_MINIAOD___20181022_235142/181022_215230/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017B-31Mar2018-v1_MINIAOD___20181022_235142/181022_215230/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017B-31Mar2018-v1_MINIAOD___20181022_235142/181022_215230/0004/sum.root");
  */

  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017C-31Mar2018-v1_MINIAOD___20181023_001009/181022_221020/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017C-31Mar2018-v1_MINIAOD___20181023_001009/181022_221020/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017C-31Mar2018-v1_MINIAOD___20181023_001009/181022_221020/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017C-31Mar2018-v1_MINIAOD___20181023_001009/181022_221020/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017C-31Mar2018-v1_MINIAOD___20181023_001009/181022_221020/0004/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017C-31Mar2018-v1_MINIAOD___20181023_001009/181022_221020/0005/sum.root");
  */
  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017D-31Mar2018-v1_MINIAOD___20181022_235327/181022_215340/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017D-31Mar2018-v1_MINIAOD___20181022_235327/181022_215340/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017D-31Mar2018-v1_MINIAOD___20181022_235327/181022_215340/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017D-31Mar2018-v1_MINIAOD___20181022_235327/181022_215340/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017D-31Mar2018-v1_MINIAOD___20181022_235327/181022_215340/0004/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017D-31Mar2018-v1_MINIAOD___20181022_235327/181022_215340/0005/sum.root");
  */
  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017E-31Mar2018-v1_MINIAOD___20181022_235401/181022_215411/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017E-31Mar2018-v1_MINIAOD___20181022_235401/181022_215411/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017E-31Mar2018-v1_MINIAOD___20181022_235401/181022_215411/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017E-31Mar2018-v1_MINIAOD___20181022_235401/181022_215411/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017E-31Mar2018-v1_MINIAOD___20181022_235401/181022_215411/0004/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017E-31Mar2018-v1_MINIAOD___20181022_235401/181022_215411/0005/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017E-31Mar2018-v1_MINIAOD___20181022_235401/181022_215411/0006/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017E-31Mar2018-v1_MINIAOD___20181022_235401/181022_215411/0007/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017E-31Mar2018-v1_MINIAOD___20181022_235401/181022_215411/0008/sum.root");
  */
  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017F-31Mar2018-v1_MINIAOD___20181023_000853/181022_220909/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017F-31Mar2018-v1_MINIAOD___20181023_000853/181022_220909/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017F-31Mar2018-v1_MINIAOD___20181023_000853/181022_220909/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017F-31Mar2018-v1_MINIAOD___20181023_000853/181022_220909/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017F-31Mar2018-v1_MINIAOD___20181023_000853/181022_220909/0004/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017F-31Mar2018-v1_MINIAOD___20181023_000853/181022_220909/0005/sum.root");
  */

  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018C-PromptReco-v3_MINIAOD___20180913_014624/180912_234632/0000/rootuple-2018-dimuonditrak_fivedataC2018_3_0_207.root");

  TString selector = "/lustre/home/adrianodif/jpsiphi/2018/data_2018/analysis/utilities/skimmers/2mu2k/TwoMuTwoK";
  TProof *p = TProof::Open("workers=5"); // 12 workers for qsub
  //gProofDebugMask = TProofDebug::kAll;
  //gProofDebugLevel = 5;

  // Processing
  cout << ">> Processing " << selector << " ... " << endl;

  TString selectorplus = selector;
  selectorplus += ".C+";
  p->Process(dataset, selectorplus);

}
