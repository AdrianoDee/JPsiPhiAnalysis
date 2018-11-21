{

  #include "TTree.h"
  #include "TDSet.h"
  #include "TProof.h"
  #include "TString.h"


  // INPUT DATA SAMPLE ON LOCAL DISK

  TDSet* dataset = new TDSet("TTree", "X_data", "mkcands");


  //2012 B
  /*
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runB_resplit_Oct17/runB_add_01.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runB_resplit_Oct17/runB_split_00.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runB_resplit_Oct17/runB_split_01.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runB_resplit_Oct17/runB_split_02.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runB_resplit_Oct17/runB_split_03.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runB_resplit_Oct17/runB_split_04.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runB_resplit_Oct17/runB_split_05.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runB_resplit_Oct17/runB_split_06.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runB_resplit_Oct17/runB_split_07.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runB_resplit_Oct17/runB_split_08.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runB_resplit_Oct17/runB_split_09.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runB_resplit_Oct17/runB_split_10.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runB_resplit_Oct17/runB_split_11.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runB_resplit_Oct17/runB_split_12.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runB_resplit_Oct17/runB_split_13.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runB_resplit_Oct17/runB_split_14.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runB_resplit_Oct17/runB_split_15.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runB_resplit_Oct17/runB_split_16.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runB_resplit_Oct17/runB_split_17.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runB_resplit_Oct17/runB_split_18.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runB_resplit_Oct17/runB_split_19.root");
  */

  //2012 C
  /*
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runC_resplit_Oct17/runC_add_01.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runC_resplit_Oct17/runC_split_00.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runC_resplit_Oct17/runC_split_01.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runC_resplit_Oct17/runC_split_02.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runC_resplit_Oct17/runC_split_03.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runC_resplit_Oct17/runC_split_04.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runC_resplit_Oct17/runC_split_05.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runC_resplit_Oct17/runC_split_06.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runC_resplit_Oct17/runC_split_07.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runC_resplit_Oct17/runC_split_08.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runC_resplit_Oct17/runC_split_09.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runC_resplit_Oct17/runC_split_10.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runC_resplit_Oct17/runC_split_11.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runC_resplit_Oct17/runC_split_12.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runC_resplit_Oct17/runC_split_13.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runC_resplit_Oct17/runC_split_14.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runC_resplit_Oct17/runC_split_15.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runC_resplit_Oct17/runC_split_16.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runC_resplit_Oct17/runC_split_17.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runC_resplit_Oct17/runC_split_18.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runC_resplit_Oct17/runC_split_19.root");
  */
  //2012 D
  /*
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runD_resplit_Oct17/runD_add_01.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runD_resplit_Oct17/runD_split_00.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runD_resplit_Oct17/runD_split_01.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runD_resplit_Oct17/runD_split_02.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runD_resplit_Oct17/runD_split_03.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runD_resplit_Oct17/runD_split_04.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runD_resplit_Oct17/runD_split_05.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runD_resplit_Oct17/runD_split_06.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runD_resplit_Oct17/runD_split_07.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runD_resplit_Oct17/runD_split_08.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runD_resplit_Oct17/runD_split_09.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runD_resplit_Oct17/runD_split_10.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runD_resplit_Oct17/runD_split_11.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runD_resplit_Oct17/runD_split_12.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runD_resplit_Oct17/runD_split_13.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runD_resplit_Oct17/runD_split_14.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runD_resplit_Oct17/runD_split_15.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runD_resplit_Oct17/runD_split_16.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runD_resplit_Oct17/runD_split_17.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runD_resplit_Oct17/runD_split_18.root");
  dataset->Add("/lustre/home/adrianodif/MuOniaParked/2012_data/runD_resplit_Oct17/runD_split_19.root");
  */
  TString selector = "/lustre/home/adrianodif/jpsiphi/2018/data_2018/analysis/utilities/skimmers/2mu2k/2012/TwoMuTwoK_2012";

  TProof *p = TProof::Open("workers=5"); // 12 workers for qsub
  //gProofDebugMask = TProofDebug::kAll;
  //gProofDebugLevel = 5;

  // Processing
  cout << ">> Processing " << selector << " ... " << endl;

  TString selectorplus = selector;
  selectorplus += ".C+";
  p->Process(dataset, selectorplus);

}
