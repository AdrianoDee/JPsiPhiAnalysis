{

#include "TTree.h"
#include "TDSet.h"
#include "TProof.h"
#include "TString.h"


  // INPUT DATA SAMPLE ON LOCAL DISK

  TDSet* dataset = new TDSet("TTree", "DiMuon PiK Tree", "rootuplePiK");
  //
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/2mukpi_miniaod_17Nov2017_BCDEF_2017.root");

  // dataset->Add("/Users/adrianodiflorio/Documents/Git/X4140/ProofLite/Y4140_testrootuple.root");
  TString selector = "/lustre/home/adrianodif/jpsiphi/analysis/utilities/skimmers/2mupik/TwoMuPiKSkim";
  TProof *p = TProof::Open("workers=40"); // 12 workers for qsub

  // Processing
  cout << ">> Processing " << selector << " ... " << endl;

  TString selectorplus = selector;
  selectorplus += ".C+";
  p->Process(dataset, selectorplus);

}
