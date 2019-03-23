{

#include "TTree.h"
#include "TDSet.h"
#include "TProof.h"
#include "TString.h"


  // INPUT DATA SAMPLE ON LOCAL DISK

  TDSet* dataset = new TDSet("TTree", "JPsiPhiTree", "rootuple");
  //
  //Y MCs
  //dataset->Add("/lustre/home/adrianodif/jpsiphi/2018/CMSSW_10_2_1/src/jpsiphi/jpsiphi/test/Y4700/y4700_official_mc_2018.root");
  //
  //B MCs
  dataset->Add("/lustre/cms/store/user/adiflori/OfficialEfficiency/Bs/BsToJpsiPhi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/bs_official_bmm.root");
  /// 2018

  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v1_MINIAOD___20181124_150539/181124_140547/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v1_MINIAOD___20181124_150539/181124_140547/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v1_MINIAOD___20181124_150539/181124_140547/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v1_MINIAOD___20181124_150539/181124_140547/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v2_MINIAOD___20181124_150752/181124_140759/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v3_MINIAOD___20181124_151626/181124_141633/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018A-PromptReco-v3_MINIAOD___20181124_151626/181124_141633/0001/sum.root");
*/
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
  /// 2017
  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017B-31Mar2018-v1_MINIAOD___20181122_114838/181122_105020/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017B-31Mar2018-v1_MINIAOD___20181122_114838/181122_105020/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017B-31Mar2018-v1_MINIAOD___20181122_114838/181122_105020/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017B-31Mar2018-v1_MINIAOD___20181122_114838/181122_105020/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017B-31Mar2018-v1_MINIAOD___20181122_114838/181122_105020/0004/sum.root");
  */
  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017C-31Mar2018-v1_MINIAOD___20181122_141444/181122_131456/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017C-31Mar2018-v1_MINIAOD___20181122_141444/181122_131456/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017C-31Mar2018-v1_MINIAOD___20181122_141444/181122_131456/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017C-31Mar2018-v1_MINIAOD___20181122_141444/181122_131456/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017C-31Mar2018-v1_MINIAOD___20181122_141444/181122_131456/0004/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017C-31Mar2018-v1_MINIAOD___20181122_141444/181122_131456/0005/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017C-31Mar2018-v1_MINIAOD___20181122_141444/181122_131456/0006/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017C-31Mar2018-v1_MINIAOD___20181122_141444/181122_131456/0007/sum.root");
  */
  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017D-31Mar2018-v1_MINIAOD___20181122_115226/181122_105239/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017D-31Mar2018-v1_MINIAOD___20181122_115226/181122_105239/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017D-31Mar2018-v1_MINIAOD___20181122_115226/181122_105239/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017D-31Mar2018-v1_MINIAOD___20181122_115226/181122_105239/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017D-31Mar2018-v1_MINIAOD___20181122_115226/181122_105239/0004/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017D-31Mar2018-v1_MINIAOD___20181122_115226/181122_105239/0005/sum.root");
  */
  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017E-31Mar2018-v1_MINIAOD___20181122_115317/181122_105331/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017E-31Mar2018-v1_MINIAOD___20181122_115317/181122_105331/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017E-31Mar2018-v1_MINIAOD___20181122_115317/181122_105331/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017E-31Mar2018-v1_MINIAOD___20181122_115317/181122_105331/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017E-31Mar2018-v1_MINIAOD___20181122_115317/181122_105331/0004/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017E-31Mar2018-v1_MINIAOD___20181122_115317/181122_105331/0005/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017E-31Mar2018-v1_MINIAOD___20181122_115317/181122_105331/0006/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017E-31Mar2018-v1_MINIAOD___20181122_115317/181122_105331/0007/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017E-31Mar2018-v1_MINIAOD___20181122_115317/181122_105331/0008/sum.root");
  */
  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017F-31Mar2018-v1_MINIAOD___20181122_141620/181122_131633/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017F-31Mar2018-v1_MINIAOD___20181122_141620/181122_131633/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017F-31Mar2018-v1_MINIAOD___20181122_141620/181122_131633/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017F-31Mar2018-v1_MINIAOD___20181122_141620/181122_131633/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017F-31Mar2018-v1_MINIAOD___20181122_141620/181122_131633/0004/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017F-31Mar2018-v1_MINIAOD___20181122_141620/181122_131633/0005/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017F-31Mar2018-v1_MINIAOD___20181122_141620/181122_131633/0006/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017F-31Mar2018-v1_MINIAOD___20181122_141620/181122_131633/0007/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2017F-31Mar2018-v1_MINIAOD___20181122_141620/181122_131633/0008/sum.root");
  */

  /// 2016

  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016B-17Jul2018_ver2-v1_MINIAOD___20181023_104525/181023_084536/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016B-17Jul2018_ver2-v1_MINIAOD___20181023_104525/181023_084536/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016B-17Jul2018_ver2-v1_MINIAOD___20181023_104525/181023_084536/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016B-17Jul2018_ver2-v1_MINIAOD___20181023_104525/181023_084536/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016B-17Jul2018_ver2-v1_MINIAOD___20181023_104525/181023_084536/0004/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016B-17Jul2018_ver2-v1_MINIAOD___20181023_104525/181023_084536/0005/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016B-17Jul2018_ver2-v1_MINIAOD___20181023_104525/181023_084536/0006/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016B-17Jul2018_ver2-v1_MINIAOD___20181023_104525/181023_084536/0007/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016B-17Jul2018_ver2-v1_MINIAOD___20181023_104525/181023_084536/0008/sum.root");
  */

  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016C-17Jul2018-v1_MINIAOD___20181023_100656/181023_080711/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016C-17Jul2018-v1_MINIAOD___20181023_100656/181023_080711/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016C-17Jul2018-v1_MINIAOD___20181023_100656/181023_080711/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016C-17Jul2018-v1_MINIAOD___20181023_100656/181023_080711/0003/sum.root");
  */
  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016D-17Jul2018-v1_MINIAOD___20181023_103007/181023_083021/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016D-17Jul2018-v1_MINIAOD___20181023_103007/181023_083021/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016D-17Jul2018-v1_MINIAOD___20181023_103007/181023_083021/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016D-17Jul2018-v1_MINIAOD___20181023_103007/181023_083021/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016D-17Jul2018-v1_MINIAOD___20181023_103007/181023_083021/0004/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016D-17Jul2018-v1_MINIAOD___20181023_103007/181023_083021/0005/sum.root");
  */
  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016E-17Jul2018-v1_MINIAOD___20181023_103035/181023_083048/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016E-17Jul2018-v1_MINIAOD___20181023_103035/181023_083048/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016E-17Jul2018-v1_MINIAOD___20181023_103035/181023_083048/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016E-17Jul2018-v1_MINIAOD___20181023_103035/181023_083048/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016E-17Jul2018-v1_MINIAOD___20181023_103035/181023_083048/0004/sum.root");
  */
  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016F-17Jul2018-v1_MINIAOD___20181023_104610/181023_084620/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016F-17Jul2018-v1_MINIAOD___20181023_104610/181023_084620/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016F-17Jul2018-v1_MINIAOD___20181023_104610/181023_084620/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016F-17Jul2018-v1_MINIAOD___20181023_104610/181023_084620/0003/sum.root");
  */
  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016G-17Jul2018-v1_MINIAOD___20181023_104628/181023_084638/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016G-17Jul2018-v1_MINIAOD___20181023_104628/181023_084638/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016G-17Jul2018-v1_MINIAOD___20181023_104628/181023_084638/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016G-17Jul2018-v1_MINIAOD___20181023_104628/181023_084638/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016G-17Jul2018-v1_MINIAOD___20181023_104628/181023_084638/0004/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016G-17Jul2018-v1_MINIAOD___20181023_104628/181023_084638/0005/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016G-17Jul2018-v1_MINIAOD___20181023_104628/181023_084638/0006/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016G-17Jul2018-v1_MINIAOD___20181023_104628/181023_084638/0007/sum.root");
  */
  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016H-17Jul2018-v1_MINIAOD___20181023_104708/181023_084719/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016H-17Jul2018-v1_MINIAOD___20181023_104708/181023_084719/0001/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016H-17Jul2018-v1_MINIAOD___20181023_104708/181023_084719/0002/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016H-17Jul2018-v1_MINIAOD___20181023_104708/181023_084719/0003/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016H-17Jul2018-v1_MINIAOD___20181023_104708/181023_084719/0004/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016H-17Jul2018-v1_MINIAOD___20181023_104708/181023_084719/0005/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016H-17Jul2018-v1_MINIAOD___20181023_104708/181023_084719/0006/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016H-17Jul2018-v1_MINIAOD___20181023_104708/181023_084719/0007/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2016H-17Jul2018-v1_MINIAOD___20181023_104708/181023_084719/0008/sum.root");
  */
  /*
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20181106_172137/181106_162145/0000/sum.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_five_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20181106_172137/181106_162145/0001/sum.root");
  */
  //dataset->Add("/lustre/home/adrianodif/jpsiphi/2018/CMSSW_10_2_1/src/jpsiphi/jpsiphi/test/bujpsiphik/bu_jpsiphik_4.root");
  //dataset->Add("/lustre/home/adrianodif/jpsiphi/2018/CMSSW_10_2_1/src/jpsiphi/jpsiphi/test/bujpsiphik/bu_jpsiphik_2.root");
  //dataset->Add("/lustre/home/adrianodif/jpsiphi/2018/CMSSW_10_2_1/src/jpsiphi/jpsiphi/test/bujpsiphik/bu_jpsiphik.root");
  dataset->Add("/lustre/home/adrianodif/jpsiphi/2018/CMSSW_10_2_1/src/jpsiphi/jpsiphi/test/bspsiphi/bs_psiphi.root");

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
