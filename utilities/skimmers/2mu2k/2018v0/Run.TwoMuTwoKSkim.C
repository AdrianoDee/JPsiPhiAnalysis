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
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/2mu2k_miniaod_07Aug17_BCDEFGH_2016_phi_097_106.root");
  //dataset->Add("/lustre/cms/store/user/adiflori/Charmonium//2mu2k_miniaod_17Nov2017_BCDEF_2017_phi_097_106.root");

  TString selector = "/lustre/home/adrianodif/jpsiphi/analysis/utilities/skimmers/2mu2k/TwoMuTwoKSkim";
  TProof *p = TProof::Open("workers=40"); // 12 workers for qsub

  // Processing
  cout << ">> Processing " << selector << " ... " << endl;

  TString selectorplus = selector;
  selectorplus += ".C+";
  p->Process(dataset, selectorplus);

}
