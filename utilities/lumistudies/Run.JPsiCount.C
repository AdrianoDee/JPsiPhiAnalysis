{

#include "TTree.h"
#include "TDSet.h"
#include "TProof.h"
#include "TString.h"


  // INPUT DATA SAMPLE ON LOCAL DISK

  TDSet* dataset = new TDSet("TTree", "dimuonTree", "rootupleMuMu");
  //
  //dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/2mu2k_miniaod_17Nov2017_BCDEF_2017.root");
  dataset->Add("/lustre/cms/store/user/adiflori/Charmonium/2mu2k_miniaod_07Aug17_BCDEFGH_2016.root");

  TString selector = "/lustre/home/adrianodif/jpsiphi/analysis/utilities/lumistudies/JPsiCount";
  TProof *p = TProof::Open("workers=40"); // 12 workers for qsub

  // Processing
  cout << ">> Processing " << selector << " ... " << endl;

  TString selectorplus = selector;
  selectorplus += ".C+";
  p->Process(dataset, selectorplus);

}
