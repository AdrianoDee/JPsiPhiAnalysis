#define SixTracks_cxx
// The class definition in SixTracks.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("SixTracks.C")
// root> T->Process("SixTracks.C","some options")
// root> T->Process("SixTracks.C+")
//


#include "SixTracks.h"
#include <TH2.h>
#include <TStyle.h>

void SixTracks::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
}

void SixTracks::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  std::string outputString = "2mu4k_six_tree.root";
  OutFile = new TProofOutputFile( outputString.data() );
  fOut = OutFile->OpenFile("RECREATE");
  if (!(fOut=OutFile->OpenFile("RECREATE")))
  {
    Warning("SlaveBegin","Problems opening file: %s%s", OutFile->GetDir(), OutFile->GetFileName() );
  }

  ////////////////// Histograms //////////////////
  JPsi_mass = 3.096916; /// pdg mass
  Phi_mass = 1.019455; /// pdg mass
  Phi_mean = 1.019723;
  Phi_sigma = 2.35607e-03;//2.28400e-03;

  outTree = new TTree("SixTrackSkimmedTree","SixTrackSkimmedTree");

  Last login: Mon Jul  1 10:00:17 on ttys006
  (base) visitor-50230508:CNNFiltering adrianodif$ ssh -X ui-centos7.recas.ba.infn.it
  Last login: Mon Jul  1 11:33:16 2019 from 2001:1458:204:1::102:3fcc
  screen -r C
  Identity added: /lustre/home/adrianodif/.ssh/id_tesla (/lustre/home/adrianodif/.ssh/id_tesla)
  [adrianodif@ui-centos7 ~]$ screen -r C
  There is a screen on:
  	55914.COH	(Attached)
  There is no screen to be resumed matching C.
  [adrianodif@ui-centos7 ~]$ screen -D
  [55914.COH power detached.]

  [adrianodif@ui-centos7 ~]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 ~]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 ~]$ cd /lustre/cms/store/user/adiflori/GEN-MINIAODSIMQCD_PhiFilter_Dau_SoftQCD_MiniTracks_/crab_GEN-MINIAODSIM_QCD_PhiFilter_Dau_SoftQCD_MiniTracks__20190701_061305/190701_041401/0001/qcd_ml
  [adrianodif@ui-centos7 qcd_ml]$ ls
  ^C
  [adrianodif@ui-centos7 qcd_ml]$ vi to^C
  [adrianodif@ui-centos7 qcd_ml]$ vi toEvent.py
  [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
  [adrianodif@ui-centos7 qcd_ml]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
  [adrianodif@ui-centos7 qcd_ml]$ screen -r C^C
  [adrianodif@ui-centos7 qcd_ml]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 qcd_ml]$ vu ^C
  [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
  [adrianodif@ui-centos7 qcd_ml]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
  [adrianodif@ui-centos7 qcd_ml]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
  [adrianodif@ui-centos7 qcd_ml]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
  [adrianodif@ui-centos7 qcd_ml]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
  [adrianodif@ui-centos7 qcd_ml]$ rm *pion*
  ^C
  [adrianodif@ui-centos7 qcd_ml]$ rm *electron* &
  [1] 32315
  [adrianodif@ui-centos7 qcd_ml]$ rm *muon*
  ^C
  [adrianodif@ui-centos7 qcd_ml]$ rm *muon* &
  [2] 34149
  [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
  [1]-  Done                    rm *electron*
  [adrianodif@ui-centos7 qcd_ml]$
  [adrianodif@ui-centos7 qcd_ml]$ screen -r C
  [detached from 55914.COH]
  [2]+  Done                    rm *muon*
  [adrianodif@ui-centos7 qcd_ml]$ vi to^C
  [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.p^C
  [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
  [adrianodif@ui-centos7 qcd_ml]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
  [adrianodif@ui-centos7 qcd_ml]$ ps
     PID TTY          TIME CMD
    8843 pts/48   00:00:01 bash
   64411 pts/48   00:00:00 ps
  [adrianodif@ui-centos7 qcd_ml]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
  [adrianodif@ui-centos7 qcd_ml]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
  ^[OA[adrianodif@ui-centos7 qcd_ml]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 qcd_ml]$ pwd
  /lustre/cms/store/user/adiflori/GEN-MINIAODSIMQCD_PhiFilter_Dau_SoftQCD_MiniTracks_/crab_GEN-MINIAODSIM_QCD_PhiFilter_Dau_SoftQCD_MiniTracks__20190701_061305/190701_041401/0001/qcd_ml
  [adrianodif@ui-centos7 qcd_ml]$ cd -
  /lustre/home/adrianodif
  [adrianodif@ui-centos7 ~]$ cd -
  /lustre/cms/store/user/adiflori/GEN-MINIAODSIMQCD_PhiFilter_Dau_SoftQCD_MiniTracks_/crab_GEN-MINIAODSIM_QCD_PhiFilter_Dau_SoftQCD_MiniTracks__20190701_061305/190701_041401/0001/qcd_ml
  [adrianodif@ui-centos7 qcd_ml]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
  [adrianodif@ui-centos7 qcd_ml]$ ^C
  [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
  c[adrianodif@ui-centos7 qcd_ml]$ screen -^C
  [adrianodif@ui-centos7 qcd_ml]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 qcd_ml]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 qcd_ml]$ screen -r ^C
  [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
  [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
  [adrianodif@ui-centos7 qcd_ml]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
  [adrianodif@ui-centos7 qcd_ml]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
  [adrianodif@ui-centos7 qcd_ml]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
  [adrianodif@ui-centos7 qcd_ml]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
  [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
  [adrianodif@ui-centos7 qcd_ml]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 qcd_ml]$ cd /lustre/cms/store/user/adiflori/GEN-MINIAODSIMQCD_PhiFilter_Dau_SoftQCD_MiniTracks_/crab_GEN-MINIAODSIM_QCD_PhiFilter_Dau_SoftQCD_MiniTracks__20190701_061305/190701_041401/0000/qcd_ml
  [adrianodif@ui-centos7 qcd_ml]$ vi pdg_0
  ^C^C^C^C^[zz^Z
  [adrianodif@ui-centos7 qcd_ml]$ cd data
  vi pdg_[adrianodif@ui-centos7 data]$ vi pdg_3
  [adrianodif@ui-centos7 data]$ vi pdg_4
  [adrianodif@ui-centos7 data]$ vi toPdgs.py^C
  [adrianodif@ui-centos7 data]$ screen -r C
  [detached from 55914.COH]
  (reverse-i-search)`v': ^C pdg_4
  [adrianodif@ui-centos7 data]$ vi toP^C
  [adrianodif@ui-centos7 data]$ vi toPdg^C
  [adrianodif@ui-centos7 data]$ vi toPdgs.py
  [adrianodif@ui-centos7 data]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 data]$ vi toP^C
  [adrianodif@ui-centos7 data]$ vi toPdgs.py
  [adrianodif@ui-centos7 data]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 data]$ vi to^C
  [adrianodif@ui-centos7 data]$ vi toPdgs.py
  [adrianodif@ui-centos7 data]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 data]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 data]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 data]$ cd /lustre/cms/store/user/adiflori/GEN-MINIAODSIMQCD_PhiFilter_Dau_SoftQCD_MiniTracks_/crab_GEN-MINIAODSIM_QCD_PhiFilter_Dau_SoftQCD_MiniTracks__20190701_061305/190701_041401/0001/qcd_ml
  [adrianodif@ui-centos7 qcd_ml]$ cp toPdgs.py toElectron.py
  [adrianodif@ui-centos7 qcd_ml]$ vi toElectron
  [adrianodif@ui-centos7 qcd_ml]$ vi toElectron.py^C
  [adrianodif@ui-centos7 qcd_ml]$ vi toElectron.py
  (reverse-i-search)`for': ^Cr f in *; do tar xvf $f -P; done
  [adrianodif@ui-centos7 qcd_ml]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 qcd_ml]$ ls -Slh *electron* | tail -10
  ^C
  [adrianodif@ui-centos7 qcd_ml]$ screen -r C
  [detached from 55914.COH]
  [adrianodif@ui-centos7 qcd_ml]$ vi toElectron.py
  [adrianodif@ui-centos7 qcd_ml]$ char^C
  [adrianodif@ui-centos7 qcd_ml]$ mystore
  [adrianodif@ui-centos7 adiflori]$ cd Charmonium/
  [adrianodif@ui-centos7 Charmonium]$ ls -tlrh
  totale 416K
  -rw-rw-r--+ 1 adrianodif cms 6,1K 14 set  2018 sumlist.txt
  -rw-rw-r--+ 1 adrianodif cms 6,0K 24 set  2018 2mu2mulist.txt
  -rw-rw-r--+ 1 adrianodif cms  24K 26 set  2018 FourMuons.h
  -rw-rw-r--+ 1 adrianodif cms 3,1K 26 set  2018 FourMuons.C
  drwxrwxr-x+ 2 adrianodif cms 4,0K  9 ott  2018 2mu2mu
  -rw-rw-r--+ 1 adrianodif cms 3,5K 27 ott  2018 2017data2mu2k
  -rw-rw-r--+ 1 adrianodif cms 177K 27 ott  2018 log
  -rw-rw-r--+ 1 adrianodif cms 7,6K 28 ott  2018 2016data
  -rw-rw-r--+ 1 adrianodif cms  847 19 nov  2018 mclist.sh
  -rw-rw-r--+ 1 adrianodif cms 6,2K 23 nov  2018 2017data
  -rw-rw-r--+ 1 adrianodif cms 4,2K 28 nov  2018 2018data
  drwxrwxr-x+ 3        497 497 4,0K 23 mar 11.27 crab_miniaod_2mu2k_five_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20190323_105234
  -rw-rw-r--+ 1 adrianodif cms  29K  9 apr 18.12 dummy.h
  -rw-rw-r--+ 1 adrianodif cms 3,1K  9 apr 18.12 dummy.C
  drwxrwxr-x+ 3        497 497 4,0K 26 giu 18.50 crab_miniaod_2mu2k_Charmonium_Run2018B-17Sep2018-v1_MINIAOD___20190626_182314
  drwxrwxr-x+ 3        497 497 4,0K 27 giu 09.02 crab_miniaod_2mu2k_Charmonium_Run2018C-17Sep2018-v1_MINIAOD___20190627_084018
  drwxrwxr-x+ 3        497 497 4,0K 27 giu 10.25 crab_miniaod_2mu2k_Charmonium_Run2018A-17Sep2018-v1_MINIAOD___20190627_095241
  drwxrwxr-x+ 3        497 497 4,0K 27 giu 11.51 crab_miniaod_2mu2k_Charmonium_Run2018C-17Sep2018-v1_MINIAOD___20190627_112903_six_five
  drwxrwxr-x+ 3        497 497 4,0K 27 giu 12.37 crab_miniaod_2mu2k_Charmonium_Run2018B-17Sep2018-v1_MINIAOD___20190627_112839_six_five
  drwxrwxr-x+ 3        497 497 4,0K 27 giu 13.30 crab_miniaod_2mu2k_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20190627_125057_six_five
  drwxrwxr-x+ 3        497 497 4,0K 27 giu 13.32 crab_miniaod_2mu2k_Charmonium_Run2018A-17Sep2018-v1_MINIAOD___20190627_124818_six_five
  [adrianodif@ui-centos7 Charmonium]$ cd crab_miniaod_2mu2k_Charmonium_Run2018A-17Sep2018-v1_MINIAOD___20190627_124818_six_five/
  [adrianodif@ui-centos7 crab_miniaod_2mu2k_Charmonium_Run2018A-17Sep2018-v1_MINIAOD___20190627_124818_six_five]$ ls -tlrh
  totale 0
  drwxrwxr-x+ 7 497 497 4,0K 30 giu 14.05 190627_104825
  [adrianodif@ui-centos7 crab_miniaod_2mu2k_Charmonium_Run2018A-17Sep2018-v1_MINIAOD___20190627_124818_six_five]$ cd 190627_104825/
  [adrianodif@ui-centos7 190627_104825]$ ls -tlrh
  totale 110G
  drwxrwxr-x+ 2        497 497 128K 28 giu 16.04 0002
  drwxrwxr-x+ 2        497 497 128K 28 giu 16.04 0001
  drwxrwxr-x+ 2        497 497 128K 28 giu 16.04 0003
  drwxrwxr-x+ 2        497 497  96K 28 giu 16.04 0004
  drwxrwxr-x+ 2        497 497 128K 28 giu 16.04 0000
  -rw-r--r--+ 1 adrianodif cms  17G 28 giu 16.34 0000.root
  -rw-r--r--+ 1 adrianodif cms  22G 28 giu 16.38 0004.root
  -rw-r--r--+ 1 adrianodif cms  23G 28 giu 16.39 0001.root
  -rw-r--r--+ 1 adrianodif cms  23G 28 giu 16.41 0003.root
  -rw-r--r--+ 1 adrianodif cms  26G 28 giu 16.44 0002.root
  -rw-r-----+ 1 adrianodif cms 1,6K 28 giu 17.00 toHdF.py
  -rw-rw-r--+ 1 adrianodif cms  137 28 giu 19.20 DU
  -rw-rw-r--+ 1 adrianodif cms  34K 29 giu 09.38 SixTracks.h
  -rw-rw-r--+ 1 adrianodif cms 3,1K 29 giu 09.38 SixTracks.C
  [adrianodif@ui-centos7 190627_104825]$ pwd
  /lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018A-17Sep2018-v1_MINIAOD___20190627_124818_six_five/190627_104825
  [adrianodif@ui-centos7 190627_104825]$ cd --
  [adrianodif@ui-centos7 ~]$ cd jpsiphi/2018/data_2018/analysis/utilities/skimmers/
  [adrianodif@ui-centos7 skimmers]$ ls -tlrh
  totale 128K
  drwxr-x---+  2 adrianodif cms 4,0K 13 set  2018 2muBkg
  drwxr-x---+  2 adrianodif cms 4,0K 13 set  2018 2mukpi
  drwxr-x---+  2 adrianodif cms 4,0K 13 set  2018 2mupik
  drwxr-x---+  2 adrianodif cms 4,0K 13 set  2018 4mu
  -rw-r-----+  1 adrianodif cms 3,1K 13 set  2018 FourMuSkim.C
  drwxr-x---+  2 adrianodif cms 4,0K 13 set  2018 legacy
  -rw-r-----+  1 adrianodif cms  54K 13 set  2018 skimRunII_xmass.C
  -rw-r-----+  1 adrianodif cms 5,1K 13 set  2018 skimRunII_xmass.ipynb
  -rw-r-----+  1 adrianodif cms 1,9K 13 set  2018 skimRunII_xmass.py
  drwxr-x---+  2 adrianodif cms 4,0K 13 set  2018 w_hlts
  drwxr-x---+  2 adrianodif cms 4,0K 13 set  2018 wo_hlts
  -rw-r-----+  1 adrianodif cms 1,4K 25 set  2018 skimvars.sh
  drwxr-x---+  2 adrianodif cms 4,0K 27 ott  2018 2mu2mu
  drwxr-x---+  3 adrianodif cms 4,0K 31 ott  2018 splotskimmers
  drwxr-x---+  2 adrianodif cms 4,0K 21 nov  2018 2mu2k_five
  drwxr-x---+ 12 adrianodif cms 4,0K 23 mar 16.34 2mu2k
  drwxr-x---+  2 adrianodif cms 4,0K 23 mar 16.43 2mu4trk
  -rw-r-----+  1 adrianodif cms  336  6 apr 23.54 DU
  drwxr-x---+  8 adrianodif cms 4,0K  6 apr 23.55 2mu3trk
  -rw-r-----+  1 adrianodif cms  11K 10 apr 15.45 merging.py
  drwxr-x---+  3 adrianodif cms 4,0K 11 apr 04.03 fourtracks
  drwxr-x---+  3 adrianodif cms 4,0K 11 apr 04.05 fivetracks
  drwxr-x---+  5 adrianodif cms 4,0K 11 apr 10.11 sixtracks
  -rw-r-----+  1 adrianodif cms 1,6K 22 mag 00.38 toHdF.py
  drwxr-x---+  2 adrianodif cms 4,0K 29 mag 07.28 2012
  drwxr-x---+  2 adrianodif cms 4,0K 29 giu 09.40 sixtracks_new
  [adrianodif@ui-centos7 skimmers]$ cd sixtracks_new/
  [adrianodif@ui-centos7 sixtracks_new]$ ls -tlrh
  totale 0
  [adrianodif@ui-centos7 sixtracks_new]$ ls -tlr^C
  [adrianodif@ui-centos7 sixtracks_new]$ cp /lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018A-17Sep2018-v1_MINIAOD___20190627_124818_six_five/190627_104825/SixTracks.* .
  [adrianodif@ui-centos7 sixtracks_new]$ cp ../skimvars.sh .
  [adrianodif@ui-centos7 sixtracks_new]$ . skimvars.sh SixTracks.h
  Use this script with the header generated by MakeSelector
  sed: espressione -e #1, carattere 2: ci sono altri caratteri dopo il comando
  [adrianodif@ui-centos7 sixtracks_new]$ vi skimvars.sh
  [adrianodif@ui-centos7 sixtracks_new]$ vi lis^C
  [adrianodif@ui-centos7 sixtracks_new]$ ls
  assign.txt  branches.txt  SixTracks.C  SixTracks.h  skimvars.sh  vars.txt
  [adrianodif@ui-centos7 sixtracks_new]$ vi assign.txt
  [adrianodif@ui-centos7 sixtracks_new]$ vi SixTracks.h









































































































































































































































































  [adrianodif@ui-centos7 sixtracks_new]$ vi SixTracks.C
  [adrianodif@ui-centos7 sixtracks_new]$ vi assign.txt
  [adrianodif@ui-centos7 sixtracks_new]$ ls
  assign.txt  branches.txt  SixTracks.C  SixTracks.h  skimvars.sh  vars.txt
  [adrianodif@ui-centos7 sixtracks_new]$ vi branches.txt

  outTree->Branch("run",  &out_run,       "run/F");
  outTree->Branch("event",        &out_event,     "event/F");
  outTree->Branch("lumi",         &out_lumi,      "lumi/F");
  outTree->Branch("numPrimaryVertices",   &out_numPrimaryVertices,        "numPrimaryVertices/F");
  outTree->Branch("trigger",      &out_trigger,   "trigger/F");
  outTree->Branch("noSixCandidates",      &out_noSixCandidates,   "noSixCandidates/F");
  outTree->Branch("five_id",      &out_five_id,   "five_id/F");
  outTree->Branch("dimuon_id",    &out_dimuon_id,         "dimuon_id/F");
  outTree->Branch("p1_id",        &out_p1_id,     "p1_id/F");
  outTree->Branch("m1_id",        &out_m1_id,     "m1_id/F");
  outTree->Branch("p2_id",        &out_p2_id,     "p2_id/F");
  outTree->Branch("m2_id",        &out_m2_id,     "m2_id/F");
  outTree->Branch("swapped",      &out_swapped,   "swapped/F");
  outTree->Branch("sameSign",     &out_sameSign,  "sameSign/F");
  outTree->Branch("sameSign_mmtt",        &out_sameSign_mmtt,     "sameSign_mmtt/F");
  outTree->Branch("six_p4",       &out_six_p4,    "six_p4/F");
  outTree->Branch("five_p4",      &out_five_p4,   "five_p4/F");
  outTree->Branch("mumukk_p4",    &out_mumukk_p4,         "mumukk_p4/F");
  outTree->Branch("ditrack_p4",   &out_ditrack_p4,        "ditrack_p4/F");
  outTree->Branch("dimuon_p4",    &out_dimuon_p4,         "dimuon_p4/F");
  outTree->Branch("lowMuon_p4",   &out_lowMuon_p4,        "lowMuon_p4/F");
  outTree->Branch("highMuon_p4",  &out_highMuon_p4,       "highMuon_p4/F");
  outTree->Branch("highKaon_p4",  &out_highKaon_p4,       "highKaon_p4/F");
  outTree->Branch("lowKaon_p4",   &out_lowKaon_p4,        "lowKaon_p4/F");
  outTree->Branch("thirdKaon_p4",         &out_thirdKaon_p4,      "thirdKaon_p4/F");
  outTree->Branch("fourthKaon_p4",        &out_fourthKaon_p4,     "fourthKaon_p4/F");
  outTree->Branch("highPion_p4",  &out_highPion_p4,       "highPion_p4/F");
  outTree->Branch("lowPion_p4",   &out_lowPion_p4,        "lowPion_p4/F");
  outTree->Branch("thirdPion_p4",         &out_thirdPion_p4,      "thirdPion_p4/F");
  outTree->Branch("fourthPion_p4",        &out_fourthPion_p4,     "fourthPion_p4/F");
  outTree->Branch("highProton_p4",        &out_highProton_p4,     "highProton_p4/F");
  outTree->Branch("lowProton_p4",         &out_lowProton_p4,      "lowProton_p4/F");
  outTree->Branch("thirdProton_p4",       &out_thirdProton_p4,    "thirdProton_p4/F");
  outTree->Branch("fourthProton_p4",      &out_fourthProton_p4,   "fourthProton_p4/F");
  outTree->Branch("mumukk_m",     &out_mumukk_m,  "mumukk_m/F");
  outTree->Branch("mumukk_pt",    &out_mumukk_pt,         "mumukk_pt/F");
  outTree->Branch("mumukk_eta",   &out_mumukk_eta,        "mumukk_eta/F");
  outTree->Branch("mumukk_phi",   &out_mumukk_phi,        "mumukk_phi/F");
  outTree->Branch("mumukk_p",     &out_mumukk_p,  "mumukk_p/F");
  outTree->Branch("dimuon_m",     &out_dimuon_m,  "dimuon_m/F");
  outTree->Branch("dimuon_pt",    &out_dimuon_pt,         "dimuon_pt/F");
  outTree->Branch("dimuon_eta",   &out_dimuon_eta,        "dimuon_eta/F");
  outTree->Branch("dimuon_phi",   &out_dimuon_phi,        "dimuon_phi/F");
  outTree->Branch("dimuon_p",     &out_dimuon_p,  "dimuon_p/F");
  outTree->Branch("highTrackMatch",       &out_highTrackMatch,    "highTrackMatch/F");
  outTree->Branch("lowTrackMatch",        &out_lowTrackMatch,     "lowTrackMatch/F");
  outTree->Branch("lowMuonMatch",         &out_lowMuonMatch,      "lowMuonMatch/F");
  outTree->Branch("highMuonMatch",        &out_highMuonMatch,     "highMuonMatch/F");
  outTree->Branch("thirdTrackMatch",      &out_thirdTrackMatch,   "thirdTrackMatch/F");
  outTree->Branch("fourthTrackMatch",     &out_fourthTrackMatch,  "fourthTrackMatch/F");
  outTree->Branch("ditrack_m",    &out_ditrack_m,         "ditrack_m/F");
  outTree->Branch("diTrackOne_pt",        &out_diTrackOne_pt,     "diTrackOne_pt/F");
  outTree->Branch("diTrackOne_eta",       &out_diTrackOne_eta,    "diTrackOne_eta/F");
  outTree->Branch("diTrackOne_phi",       &out_diTrackOne_phi,    "diTrackOne_phi/F");
  outTree->Branch("diTrackOne_p",         &out_diTrackOne_p,      "diTrackOne_p/F");
  outTree->Branch("diTrackTwo_pt",        &out_diTrackTwo_pt,     "diTrackTwo_pt/F");
  outTree->Branch("diTrackTwo_eta",       &out_diTrackTwo_eta,    "diTrackTwo_eta/F");
  outTree->Branch("diTrackTwo_phi",       &out_diTrackTwo_phi,    "diTrackTwo_phi/F");
  outTree->Branch("diTrackTwo_p",         &out_diTrackTwo_p,      "diTrackTwo_p/F");
  outTree->Branch("diTrackThree_pt",      &out_diTrackThree_pt,   "diTrackThree_pt/F");
  outTree->Branch("diTrackThree_eta",     &out_diTrackThree_eta,  "diTrackThree_eta/F");
  outTree->Branch("diTrackThree_phi",     &out_diTrackThree_phi,  "diTrackThree_phi/F");
  outTree->Branch("diTrackThree_p",       &out_diTrackThree_p,    "diTrackThree_p/F");
  outTree->Branch("diTrackFour_pt",       &out_diTrackFour_pt,    "diTrackFour_pt/F");
  outTree->Branch("diTrackFour_eta",      &out_diTrackFour_eta,   "diTrackFour_eta/F");
  outTree->Branch("diTrackFour_phi",      &out_diTrackFour_phi,   "diTrackFour_phi/F");
  outTree->Branch("diTrackFour_p",        &out_diTrackFour_p,     "diTrackFour_p/F");
  outTree->Branch("diTrackFive_pt",       &out_diTrackFive_pt,    "diTrackFive_pt/F");
  outTree->Branch("diTrackFive_eta",      &out_diTrackFive_eta,   "diTrackFive_eta/F");
  outTree->Branch("diTrackFive_phi",      &out_diTrackFive_phi,   "diTrackFive_phi/F");
  outTree->Branch("diTrackFive_p",        &out_diTrackFive_p,     "diTrackFive_p/F");
  outTree->Branch("diTrackSix_pt",        &out_diTrackSix_pt,     "diTrackSix_pt/F");
  outTree->Branch("diTrackSix_eta",       &out_diTrackSix_eta,    "diTrackSix_eta/F");
  outTree->Branch("diTrackSix_phi",       &out_diTrackSix_phi,    "diTrackSix_phi/F");
  outTree->Branch("diTrackSix_p",         &out_diTrackSix_p,      "diTrackSix_p/F");
  outTree->Branch("dimuonDiTrkOne_mmpp",  &out_dimuonDiTrkOne_mmpp,       "dimuonDiTrkOne_mmpp/F");
  outTree->Branch("dimuonDiTrkTwo_mmpp",  &out_dimuonDiTrkTwo_mmpp,       "dimuonDiTrkTwo_mmpp/F");
  outTree->Branch("dimuonDiTrkThree_mmpp",        &out_dimuonDiTrkThree_mmpp,     "dimuonDiTrkThree_mmpp/F");
  outTree->Branch("dimuonDiTrkFour_mmpp",         &out_dimuonDiTrkFour_mmpp,      "dimuonDiTrkFour_mmpp/F");
  outTree->Branch("dimuonDiTrkOne_mmkk",  &out_dimuonDiTrkOne_mmkk,       "dimuonDiTrkOne_mmkk/F");
  outTree->Branch("dimuonDiTrkTwo_mmkk",  &out_dimuonDiTrkTwo_mmkk,       "dimuonDiTrkTwo_mmkk/F");
  outTree->Branch("dimuonDiTrkThree_mmkk",        &out_dimuonDiTrkThree_mmkk,     "dimuonDiTrkThree_mmkk/F");
  outTree->Branch("dimuonDiTrkFour_mmkk",         &out_dimuonDiTrkFour_mmkk,      "dimuonDiTrkFour_mmkk/F");
  outTree->Branch("dimuonDiTrkOne_mmpk",  &out_dimuonDiTrkOne_mmpk,       "dimuonDiTrkOne_mmpk/F");
  outTree->Branch("dimuonDiTrkTwo_mmpk",  &out_dimuonDiTrkTwo_mmpk,       "dimuonDiTrkTwo_mmpk/F");
  outTree->Branch("dimuonDiTrkThree_mmpk",        &out_dimuonDiTrkThree_mmpk,     "dimuonDiTrkThree_mmpk/F");
  outTree->Branch("dimuonDiTrkFour_mmpk",         &out_dimuonDiTrkFour_mmpk,      "dimuonDiTrkFour_mmpk/F");
  outTree->Branch("dimuonDiTrkOne_mmkp",  &out_dimuonDiTrkOne_mmkp,       "dimuonDiTrkOne_mmkp/F");
  outTree->Branch("dimuonDiTrkTwo_mmkp",  &out_dimuonDiTrkTwo_mmkp,       "dimuonDiTrkTwo_mmkp/F");
  outTree->Branch("dimuonDiTrkThree_mmkp",        &out_dimuonDiTrkThree_mmkp,     "dimuonDiTrkThree_mmkp/F");
  outTree->Branch("dimuonDiTrkFour_mmkp",         &out_dimuonDiTrkFour_mmkp,      "dimuonDiTrkFour_mmkp/F");
  outTree->Branch("diTrackOne_kk",        &out_diTrackOne_kk,     "diTrackOne_kk/F");
  outTree->Branch("diTrackTwo_kk",        &out_diTrackTwo_kk,     "diTrackTwo_kk/F");
  outTree->Branch("diTrackThree_kk",      &out_diTrackThree_kk,   "diTrackThree_kk/F");
  outTree->Branch("diTrackFour_kk",       &out_diTrackFour_kk,    "diTrackFour_kk/F");
  outTree->Branch("diTrackFive_kk",       &out_diTrackFive_kk,    "diTrackFive_kk/F");
  outTree->Branch("diTrackSix_kk",        &out_diTrackSix_kk,     "diTrackSix_kk/F");
  outTree->Branch("diTrackOne_pp",        &out_diTrackOne_pp,     "diTrackOne_pp/F");
  outTree->Branch("diTrackTwo_pp",        &out_diTrackTwo_pp,     "diTrackTwo_pp/F");
  outTree->Branch("diTrackThree_pp",      &out_diTrackThree_pp,   "diTrackThree_pp/F");
  outTree->Branch("diTrackFour_pp",       &out_diTrackFour_pp,    "diTrackFour_pp/F");
  outTree->Branch("diTrackFive_pp",       &out_diTrackFive_pp,    "diTrackFive_pp/F");
  outTree->Branch("diTrackSix_pp",        &out_diTrackSix_pp,     "diTrackSix_pp/F");
  outTree->Branch("diTrackOne_pk",        &out_diTrackOne_pk,     "diTrackOne_pk/F");
  outTree->Branch("diTrackTwo_pk",        &out_diTrackTwo_pk,     "diTrackTwo_pk/F");
  outTree->Branch("diTrackThree_pk",      &out_diTrackThree_pk,   "diTrackThree_pk/F");
  outTree->Branch("diTrackFour_pk",       &out_diTrackFour_pk,    "diTrackFour_pk/F");
  outTree->Branch("diTrackFive_pk",       &out_diTrackFive_pk,    "diTrackFive_pk/F");
  outTree->Branch("diTrackSix_pk",        &out_diTrackSix_pk,     "diTrackSix_pk/F");
  outTree->Branch("diTrackOne_kp",        &out_diTrackOne_kp,     "diTrackOne_kp/F");
  outTree->Branch("diTrackTwo_kp",        &out_diTrackTwo_kp,     "diTrackTwo_kp/F");
  outTree->Branch("diTrackThree_kp",      &out_diTrackThree_kp,   "diTrackThree_kp/F");
  outTree->Branch("diTrackFour_kp",       &out_diTrackFour_kp,    "diTrackFour_kp/F");
  outTree->Branch("diTrackFive_kp",       &out_diTrackFive_kp,    "diTrackFive_kp/F");
  outTree->Branch("diTrackSix_kp",        &out_diTrackSix_kp,     "diTrackSix_kp/F");
  outTree->Branch("highMuon_pt",  &out_highMuon_pt,       "highMuon_pt/F");
  outTree->Branch("highMuon_eta",         &out_highMuon_eta,      "highMuon_eta/F");
  outTree->Branch("highMuon_phi",         &out_highMuon_phi,      "highMuon_phi/F");
  outTree->Branch("highMuon_charge",      &out_highMuon_charge,   "highMuon_charge/F");
  outTree->Branch("highMuon_dz",  &out_highMuon_dz,       "highMuon_dz/F");
  outTree->Branch("highMuon_dxy",         &out_highMuon_dxy,      "highMuon_dxy/F");
  outTree->Branch("lowMuon_pt",   &out_lowMuon_pt,        "lowMuon_pt/F");
  outTree->Branch("lowMuon_eta",  &out_lowMuon_eta,       "lowMuon_eta/F");
  outTree->Branch("lowMuon_phi",  &out_lowMuon_phi,       "lowMuon_phi/F");
  outTree->Branch("lowMuon_charge",       &out_lowMuon_charge,    "lowMuon_charge/F");
  outTree->Branch("lowMuon_dz",   &out_lowMuon_dz,        "lowMuon_dz/F");
  outTree->Branch("lowMuon_dxy",  &out_lowMuon_dxy,       "lowMuon_dxy/F");
  outTree->Branch("highTrack_pt",         &out_highTrack_pt,      "highTrack_pt/F");
  outTree->Branch("highTrack_eta",        &out_highTrack_eta,     "highTrack_eta/F");
  outTree->Branch("highTrack_phi",        &out_highTrack_phi,     "highTrack_phi/F");
  outTree->Branch("highTrack_charge",     &out_highTrack_charge,  "highTrack_charge/F");
  outTree->Branch("highTrack_dz",         &out_highTrack_dz,      "highTrack_dz/F");
  outTree->Branch("highTrack_dxy",        &out_highTrack_dxy,     "highTrack_dxy/F");
  outTree->Branch("lowTrack_pt",  &out_lowTrack_pt,       "lowTrack_pt/F");
  outTree->Branch("lowTrack_eta",         &out_lowTrack_eta,      "lowTrack_eta/F");
  outTree->Branch("lowTrack_phi",         &out_lowTrack_phi,      "lowTrack_phi/F");
  outTree->Branch("lowTrack_charge",      &out_lowTrack_charge,   "lowTrack_charge/F");
  outTree->Branch("lowTrack_dz",  &out_lowTrack_dz,       "lowTrack_dz/F");
  outTree->Branch("lowTrack_dxy",         &out_lowTrack_dxy,      "lowTrack_dxy/F");
  outTree->Branch("thirdTrack_pt",        &out_thirdTrack_pt,     "thirdTrack_pt/F");
  outTree->Branch("thirdTrack_eta",       &out_thirdTrack_eta,    "thirdTrack_eta/F");
  outTree->Branch("thirdTrack_phi",       &out_thirdTrack_phi,    "thirdTrack_phi/F");
  outTree->Branch("thirdTrack_charge",    &out_thirdTrack_charge,         "thirdTrack_charge/F");
  outTree->Branch("thirdTrack_dz",        &out_thirdTrack_dz,     "thirdTrack_dz/F");
  outTree->Branch("thirdTrack_dxy",       &out_thirdTrack_dxy,    "thirdTrack_dxy/F");
  outTree->Branch("dimuonDiTrkOne_pt",    &out_dimuonDiTrkOne_pt,         "dimuonDiTrkOne_pt/F");
  outTree->Branch("dimuonDiTrkOne_eta",   &out_dimuonDiTrkOne_eta,        "dimuonDiTrkOne_eta/F");
  outTree->Branch("dimuonDiTrkOne_phi",   &out_dimuonDiTrkOne_phi,        "dimuonDiTrkOne_phi/F");
  outTree->Branch("dimuonDiTrkOne_charge",        &out_dimuonDiTrkOne_charge,     "dimuonDiTrkOne_charge/F");
  outTree->Branch("dimuonDiTrkOne_p",     &out_dimuonDiTrkOne_p,  "dimuonDiTrkOne_p/F");
  outTree->Branch("dimuonDiTrkTwo_pt",    &out_dimuonDiTrkTwo_pt,         "dimuonDiTrkTwo_pt/F");
  outTree->Branch("dimuonDiTrkTwo_eta",   &out_dimuonDiTrkTwo_eta,        "dimuonDiTrkTwo_eta/F");
  outTree->Branch("dimuonDiTrkTwo_phi",   &out_dimuonDiTrkTwo_phi,        "dimuonDiTrkTwo_phi/F");
  outTree->Branch("dimuonDiTrkTwo_charge",        &out_dimuonDiTrkTwo_charge,     "dimuonDiTrkTwo_charge/F");
  outTree->Branch("dimuonDiTrkTwo_p",     &out_dimuonDiTrkTwo_p,  "dimuonDiTrkTwo_p/F");
  outTree->Branch("dimuonDiTrkThree_pt",  &out_dimuonDiTrkThree_pt,       "dimuonDiTrkThree_pt/F");
  outTree->Branch("dimuonDiTrkThree_eta",         &out_dimuonDiTrkThree_eta,      "dimuonDiTrkThree_eta/F");
  outTree->Branch("dimuonDiTrkThree_phi",         &out_dimuonDiTrkThree_phi,      "dimuonDiTrkThree_phi/F");
  outTree->Branch("dimuonDiTrkThree_charge",      &out_dimuonDiTrkThree_charge,   "dimuonDiTrkThree_charge/F");
  outTree->Branch("dimuonDiTrkThree_p",   &out_dimuonDiTrkThree_p,        "dimuonDiTrkThree_p/F");
  outTree->Branch("dimuonDiTrkFour_pt",   &out_dimuonDiTrkFour_pt,        "dimuonDiTrkFour_pt/F");
  outTree->Branch("dimuonDiTrkFour_eta",  &out_dimuonDiTrkFour_eta,       "dimuonDiTrkFour_eta/F");
  outTree->Branch("dimuonDiTrkFour_phi",  &out_dimuonDiTrkFour_phi,       "dimuonDiTrkFour_phi/F");
  outTree->Branch("dimuonDiTrkFour_charge",       &out_dimuonDiTrkFour_charge,    "dimuonDiTrkFour_charge/F");
  outTree->Branch("dimuonDiTrkFour_p",    &out_dimuonDiTrkFour_p,         "dimuonDiTrkFour_p/F");
  outTree->Branch("dimuonDiTrkFive_pt",   &out_dimuonDiTrkFive_pt,        "dimuonDiTrkFive_pt/F");
  outTree->Branch("dimuonDiTrkFive_eta",  &out_dimuonDiTrkFive_eta,       "dimuonDiTrkFive_eta/F");
  outTree->Branch("dimuonDiTrkFive_phi",  &out_dimuonDiTrkFive_phi,       "dimuonDiTrkFive_phi/F");
  outTree->Branch("dimuonDiTrkFive_charge",       &out_dimuonDiTrkFive_charge,    "dimuonDiTrkFive_charge/F");
  outTree->Branch("dimuonDiTrkFive_p",    &out_dimuonDiTrkFive_p,         "dimuonDiTrkFive_p/F");
  outTree->Branch("dimuonDiTrkSix_pt",    &out_dimuonDiTrkSix_pt,         "dimuonDiTrkSix_pt/F");
  outTree->Branch("dimuonDiTrkSix_eta",   &out_dimuonDiTrkSix_eta,        "dimuonDiTrkSix_eta/F");
  outTree->Branch("dimuonDiTrkSix_phi",   &out_dimuonDiTrkSix_phi,        "dimuonDiTrkSix_phi/F");
  outTree->Branch("dimuonDiTrkSix_charge",        &out_dimuonDiTrkSix_charge,     "dimuonDiTrkSix_charge/F");
  outTree->Branch("dimuonDiTrkSix_p",     &out_dimuonDiTrkSix_p,  "dimuonDiTrkSix_p/F");
  outTree->Branch("dimuon_vProb",         &out_dimuon_vProb,      "dimuon_vProb/F");
  outTree->Branch("dimuon_vChi2",         &out_dimuon_vChi2,      "dimuon_vChi2/F");
  outTree->Branch("dimuon_DCA",   &out_dimuon_DCA,        "dimuon_DCA/F");
  outTree->Branch("dimuon_ctauPV",        &out_dimuon_ctauPV,     "dimuon_ctauPV/F");
  outTree->Branch("dimuon_ctauErrPV",     &out_dimuon_ctauErrPV,  "dimuon_ctauErrPV/F");
  outTree->Branch("dimuon_cosAlpha",      &out_dimuon_cosAlpha,   "dimuon_cosAlpha/F");
  outTree->Branch("triTrackOne_kkk",      &out_triTrackOne_kkk,   "triTrackOne_kkk/F");
  outTree->Branch("triTrackTwo_kkk",      &out_triTrackTwo_kkk,   "triTrackTwo_kkk/F");
  outTree->Branch("triTrackThree_kkk",    &out_triTrackThree_kkk,         "triTrackThree_kkk/F");
  outTree->Branch("triTrackFour_kkk",     &out_triTrackFour_kkk,  "triTrackFour_kkk/F");
  outTree->Branch("triTrackOne_kkp",      &out_triTrackOne_kkp,   "triTrackOne_kkp/F");
  outTree->Branch("triTrackTwo_kkp",      &out_triTrackTwo_kkp,   "triTrackTwo_kkp/F");
  outTree->Branch("triTrackThree_kkp",    &out_triTrackThree_kkp,         "triTrackThree_kkp/F");
  outTree->Branch("triTrackOne_kpp",      &out_triTrackOne_kpp,   "triTrackOne_kpp/F");
  outTree->Branch("triTrackTwo_kpp",      &out_triTrackTwo_kpp,   "triTrackTwo_kpp/F");
  outTree->Branch("triTrackThree_kpp",    &out_triTrackThree_kpp,         "triTrackThree_kpp/F");
  outTree->Branch("triTrackOne_ppp",      &out_triTrackOne_ppp,   "triTrackOne_ppp/F");
  outTree->Branch("triTrackTwo_ppp",      &out_triTrackTwo_ppp,   "triTrackTwo_ppp/F");
  outTree->Branch("triTrackThree_ppp",    &out_triTrackThree_ppp,         "triTrackThree_ppp/F");
  outTree->Branch("triTrackOne_pt",       &out_triTrackOne_pt,    "triTrackOne_pt/F");
  outTree->Branch("triTrackOne_eta",      &out_triTrackOne_eta,   "triTrackOne_eta/F");
  outTree->Branch("triTrackOne_phi",      &out_triTrackOne_phi,   "triTrackOne_phi/F");
  outTree->Branch("triTrackOne_charge",   &out_triTrackOne_charge,        "triTrackOne_charge/F");
  outTree->Branch("triTrackTwo_pt",       &out_triTrackTwo_pt,    "triTrackTwo_pt/F");
  outTree->Branch("triTrackTwo_eta",      &out_triTrackTwo_eta,   "triTrackTwo_eta/F");
  outTree->Branch("triTrackTwo_phi",      &out_triTrackTwo_phi,   "triTrackTwo_phi/F");
  outTree->Branch("triTrackTwo_charge",   &out_triTrackTwo_charge,        "triTrackTwo_charge/F");
  outTree->Branch("triTrackThree_pt",     &out_triTrackThree_pt,  "triTrackThree_pt/F");
  outTree->Branch("triTrackThree_eta",    &out_triTrackThree_eta,         "triTrackThree_eta/F");
  outTree->Branch("triTrackThree_phi",    &out_triTrackThree_phi,         "triTrackThree_phi/F");
  outTree->Branch("triTrackThree_charge",         &out_triTrackThree_charge,      "triTrackThree_charge/F");
  outTree->Branch("triTrackFour_pt",      &out_triTrackFour_pt,   "triTrackFour_pt/F");
  outTree->Branch("triTrackFour_eta",     &out_triTrackFour_eta,  "triTrackFour_eta/F");
  outTree->Branch("triTrackFour_phi",     &out_triTrackFour_phi,  "triTrackFour_phi/F");
  outTree->Branch("triTrackFour_charge",  &out_triTrackFour_charge,       "triTrackFour_charge/F");
  outTree->Branch("mumukk_vProb",         &out_mumukk_vProb,      "mumukk_vProb/F");
  outTree->Branch("mumukk_vChi2",         &out_mumukk_vChi2,      "mumukk_vChi2/F");
  outTree->Branch("mumukk_nDof",  &out_mumukk_nDof,       "mumukk_nDof/F");
  outTree->Branch("mumukk_charge",        &out_mumukk_charge,     "mumukk_charge/F");
  outTree->Branch("mumukk_cosAlpha",      &out_mumukk_cosAlpha,   "mumukk_cosAlpha/F");
  outTree->Branch("mumukk_ctauPV",        &out_mumukk_ctauPV,     "mumukk_ctauPV/F");
  outTree->Branch("mumukk_ctauErrPV",     &out_mumukk_ctauErrPV,  "mumukk_ctauErrPV/F");
  outTree->Branch("mumukk_cosAlphaCA",    &out_mumukk_cosAlphaCA,         "mumukk_cosAlphaCA/F");
  outTree->Branch("mumukk_ctauPVCA",      &out_mumukk_ctauPVCA,   "mumukk_ctauPVCA/F");
  outTree->Branch("mumukk_ctauErrPVCA",   &out_mumukk_ctauErrPVCA,        "mumukk_ctauErrPVCA/F");
  outTree->Branch("mumukk_cosAlphaDZ",    &out_mumukk_cosAlphaDZ,         "mumukk_cosAlphaDZ/F");
  outTree->Branch("mumukk_ctauPVDZ",      &out_mumukk_ctauPVDZ,   "mumukk_ctauPVDZ/F");
  outTree->Branch("mumukk_ctauErrPVDZ",   &out_mumukk_ctauErrPVDZ,        "mumukk_ctauErrPVDZ/F");
  outTree->Branch("mumukk_cosAlphaBS",    &out_mumukk_cosAlphaBS,         "mumukk_cosAlphaBS/F");
  outTree->Branch("mumukk_ctauPVBS",      &out_mumukk_ctauPVBS,   "mumukk_ctauPVBS/F");
  outTree->Branch("mumukk_ctauErrPVBS",   &out_mumukk_ctauErrPVBS,        "mumukk_ctauErrPVBS/F");
  outTree->Branch("mumukk_vx",    &out_mumukk_vx,         "mumukk_vx/F");
  outTree->Branch("mumukk_vy",    &out_mumukk_vy,         "mumukk_vy/F");
  outTree->Branch("mumukk_vz",    &out_mumukk_vz,         "mumukk_vz/F");
  outTree->Branch("dca_m1m2",     &out_dca_m1m2,  "dca_m1m2/F");
  outTree->Branch("dca_m1t1",     &out_dca_m1t1,  "dca_m1t1/F");
  outTree->Branch("dca_m1t2",     &out_dca_m1t2,  "dca_m1t2/F");
  outTree->Branch("dca_m2t1",     &out_dca_m2t1,  "dca_m2t1/F");
  outTree->Branch("dca_m2t2",     &out_dca_m2t2,  "dca_m2t2/F");
  outTree->Branch("dca_t1t2",     &out_dca_t1t2,  "dca_t1t2/F");
  outTree->Branch("dca_m1t3",     &out_dca_m1t3,  "dca_m1t3/F");
  outTree->Branch("dca_m2t3",     &out_dca_m2t3,  "dca_m2t3/F");
  outTree->Branch("dca_t1t3",     &out_dca_t1t3,  "dca_t1t3/F");
  outTree->Branch("dca_t2t3",     &out_dca_t2t3,  "dca_t2t3/F");
  outTree->Branch("dca_m1t4",     &out_dca_m1t4,  "dca_m1t4/F");
  outTree->Branch("dca_m2t4",     &out_dca_m2t4,  "dca_m2t4/F");
  outTree->Branch("dca_t1t4",     &out_dca_t1t4,  "dca_t1t4/F");
  outTree->Branch("dca_t2t4",     &out_dca_t2t4,  "dca_t2t4/F");
  outTree->Branch("dca_t3t4",     &out_dca_t3t4,  "dca_t3t4/F");
  outTree->Branch("highTrackMuonDR",      &out_highTrackMuonDR,   "highTrackMuonDR/F");
  outTree->Branch("highTrackMuonDP",      &out_highTrackMuonDP,   "highTrackMuonDP/F");
  outTree->Branch("highTrackMuonDPt",     &out_highTrackMuonDPt,  "highTrackMuonDPt/F");
  outTree->Branch("lowTrackMuonDR",       &out_lowTrackMuonDR,    "lowTrackMuonDR/F");
  outTree->Branch("lowTrackMuonDP",       &out_lowTrackMuonDP,    "lowTrackMuonDP/F");
  outTree->Branch("lowTrackMuonDPt",      &out_lowTrackMuonDPt,   "lowTrackMuonDPt/F");
  outTree->Branch("thirdTrackMuonDR",     &out_thirdTrackMuonDR,  "thirdTrackMuonDR/F");
  outTree->Branch("thirdTrackMuonDP",     &out_thirdTrackMuonDP,  "thirdTrackMuonDP/F");
  outTree->Branch("thirdTrackMuonDPt",    &out_thirdTrackMuonDPt,         "thirdTrackMuonDPt/F");
  outTree->Branch("fourthTrackMuonDR",    &out_fourthTrackMuonDR,         "fourthTrackMuonDR/F");
  outTree->Branch("fourthTrackMuonDP",    &out_fourthTrackMuonDP,         "fourthTrackMuonDP/F");
  outTree->Branch("fourthTrackMuonDPt",   &out_fourthTrackMuonDPt,        "fourthTrackMuonDPt/F");
  outTree->Branch("tPFromPV",     &out_tPFromPV,  "tPFromPV/F");
  outTree->Branch("tMFromPV",     &out_tMFromPV,  "tMFromPV/F");
  outTree->Branch("tTFromPV",     &out_tTFromPV,  "tTFromPV/F");
  outTree->Branch("tFFromPV",     &out_tFFromPV,  "tFFromPV/F");
  outTree->Branch("tPFromPVCA",   &out_tPFromPVCA,        "tPFromPVCA/F");
  outTree->Branch("tMFromPVCA",   &out_tMFromPVCA,        "tMFromPVCA/F");
  outTree->Branch("tTFromPVCA",   &out_tTFromPVCA,        "tTFromPVCA/F");
  outTree->Branch("tFFromPVCA",   &out_tFFromPVCA,        "tFFromPVCA/F");
  outTree->Branch("tPFromPVDZ",   &out_tPFromPVDZ,        "tPFromPVDZ/F");
  outTree->Branch("tMFromPVDZ",   &out_tMFromPVDZ,        "tMFromPVDZ/F");
  outTree->Branch("tTFromPVDZ",   &out_tTFromPVDZ,        "tTFromPVDZ/F");
  outTree->Branch("tFFromPVDZ",   &out_tFFromPVDZ,        "tFFromPVDZ/F");
  outTree->Branch("five_m",       &out_five_m,    "five_m/F");
  outTree->Branch("five_m_ref",   &out_five_m_ref,        "five_m_ref/F");
  outTree->Branch("five_mass_ppk",        &out_five_mass_ppk,     "five_mass_ppk/F");
  outTree->Branch("five_mass_kpp",        &out_five_mass_kpp,     "five_mass_kpp/F");
  outTree->Branch("five_mass_pkp",        &out_five_mass_pkp,     "five_mass_pkp/F");
  outTree->Branch("five_mass_ppp",        &out_five_mass_ppp,     "five_mass_ppp/F");
  outTree->Branch("fiveOne_pt",   &out_fiveOne_pt,        "fiveOne_pt/F");
  outTree->Branch("fiveOne_eta",  &out_fiveOne_eta,       "fiveOne_eta/F");
  outTree->Branch("fiveOne_phi",  &out_fiveOne_phi,       "fiveOne_phi/F");
  outTree->Branch("fiveOne_p",    &out_fiveOne_p,         "fiveOne_p/F");
  outTree->Branch("fiveTwo_pt",   &out_fiveTwo_pt,        "fiveTwo_pt/F");
  outTree->Branch("fiveTwo_eta",  &out_fiveTwo_eta,       "fiveTwo_eta/F");
  outTree->Branch("fiveTwo_phi",  &out_fiveTwo_phi,       "fiveTwo_phi/F");
  outTree->Branch("fiveTwo_p",    &out_fiveTwo_p,         "fiveTwo_p/F");
  outTree->Branch("fiveThree_pt",         &out_fiveThree_pt,      "fiveThree_pt/F");
  outTree->Branch("fiveThree_eta",        &out_fiveThree_eta,     "fiveThree_eta/F");
  outTree->Branch("fiveThree_phi",        &out_fiveThree_phi,     "fiveThree_phi/F");
  outTree->Branch("fiveThree_p",  &out_fiveThree_p,       "fiveThree_p/F");
  outTree->Branch("fiveFour_pt",  &out_fiveFour_pt,       "fiveFour_pt/F");
  outTree->Branch("fiveFour_eta",         &out_fiveFour_eta,      "fiveFour_eta/F");
  outTree->Branch("fiveFour_phi",         &out_fiveFour_phi,      "fiveFour_phi/F");
  outTree->Branch("fiveFour_p",   &out_fiveFour_p,        "fiveFour_p/F");
  outTree->Branch("fiveFive_pt",  &out_fiveFive_pt,       "fiveFive_pt/F");
  outTree->Branch("fiveFive_eta",         &out_fiveFive_eta,      "fiveFive_eta/F");
  outTree->Branch("fiveFive_phi",         &out_fiveFive_phi,      "fiveFive_phi/F");
  outTree->Branch("fiveFive_p",   &out_fiveFive_p,        "fiveFive_p/F");
  outTree->Branch("five_cosAlpha",        &out_five_cosAlpha,     "five_cosAlpha/F");
  outTree->Branch("five_ctauPV",  &out_five_ctauPV,       "five_ctauPV/F");
  outTree->Branch("five_ctauErrPV",       &out_five_ctauErrPV,    "five_ctauErrPV/F");
  outTree->Branch("five_cosAlphaCA",      &out_five_cosAlphaCA,   "five_cosAlphaCA/F");
  outTree->Branch("five_ctauPVCA",        &out_five_ctauPVCA,     "five_ctauPVCA/F");
  outTree->Branch("five_ctauErrPVCA",     &out_five_ctauErrPVCA,  "five_ctauErrPVCA/F");
  outTree->Branch("five_cosAlphaDZ",      &out_five_cosAlphaDZ,   "five_cosAlphaDZ/F");
  outTree->Branch("five_ctauPVDZ",        &out_five_ctauPVDZ,     "five_ctauPVDZ/F");
  outTree->Branch("five_ctauErrPVDZ",     &out_five_ctauErrPVDZ,  "five_ctauErrPVDZ/F");
  outTree->Branch("five_cosAlphaBS",      &out_five_cosAlphaBS,   "five_cosAlphaBS/F");
  outTree->Branch("five_ctauPVBS",        &out_five_ctauPVBS,     "five_ctauPVBS/F");
  outTree->Branch("five_ctauErrPVBS",     &out_five_ctauErrPVBS,  "five_ctauErrPVBS/F");
  outTree->Branch("five_vProb",   &out_five_vProb,        "five_vProb/F");
  outTree->Branch("five_nDof",    &out_five_nDof,         "five_nDof/F");
  outTree->Branch("five_vChi2",   &out_five_vChi2,        "five_vChi2/F");
  outTree->Branch("five_vx",      &out_five_vx,   "five_vx/F");
  outTree->Branch("five_vy",      &out_five_vy,   "five_vy/F");
  outTree->Branch("five_vz",      &out_five_vz,   "five_vz/F");
  outTree->Branch("five_charge",  &out_five_charge,       "five_charge/F");
  outTree->Branch("bestPV_X",     &out_bestPV_X,  "bestPV_X/F");
  outTree->Branch("bestPV_Y",     &out_bestPV_Y,  "bestPV_Y/F");
  outTree->Branch("bestPV_Z",     &out_bestPV_Z,  "bestPV_Z/F");
  outTree->Branch("cosAlphaPV_X",         &out_cosAlphaPV_X,      "cosAlphaPV_X/F");
  outTree->Branch("cosAlphaPV_Y",         &out_cosAlphaPV_Y,      "cosAlphaPV_Y/F");
  outTree->Branch("cosAlphaPV_Z",         &out_cosAlphaPV_Z,      "cosAlphaPV_Z/F");
  outTree->Branch("bS_X",         &out_bS_X,      "bS_X/F");
  outTree->Branch("bS_Y",         &out_bS_Y,      "bS_Y/F");
  outTree->Branch("bS_Z",         &out_bS_Z,      "bS_Z/F");
  outTree->Branch("zPV_X",        &out_zPV_X,     "zPV_X/F");
  outTree->Branch("zPV_Y",        &out_zPV_Y,     "zPV_Y/F");
  outTree->Branch("zPV_Z",        &out_zPV_Z,     "zPV_Z/F");
  outTree->Branch("lowMuon_isTight",      &out_lowMuon_isTight,   "lowMuon_isTight/F");
  outTree->Branch("lowMuon_isLoose",      &out_lowMuon_isLoose,   "lowMuon_isLoose/F");
  outTree->Branch("lowMuon_isSoft",       &out_lowMuon_isSoft,    "lowMuon_isSoft/F");
  outTree->Branch("lowMuon_isMedium",     &out_lowMuon_isMedium,  "lowMuon_isMedium/F");
  outTree->Branch("lowMuon_isHighPt",     &out_lowMuon_isHighPt,  "lowMuon_isHighPt/F");
  outTree->Branch("lowMuon_isTracker",    &out_lowMuon_isTracker,         "lowMuon_isTracker/F");
  outTree->Branch("lowMuon_isGlobal",     &out_lowMuon_isGlobal,  "lowMuon_isGlobal/F");
  outTree->Branch("lowMuon_NPixelHits",   &out_lowMuon_NPixelHits,        "lowMuon_NPixelHits/F");
  outTree->Branch("lowMuon_NStripHits",   &out_lowMuon_NStripHits,        "lowMuon_NStripHits/F");
  outTree->Branch("lowMuon_NTrackhits",   &out_lowMuon_NTrackhits,        "lowMuon_NTrackhits/F");
  outTree->Branch("lowMuon_NBPixHits",    &out_lowMuon_NBPixHits,         "lowMuon_NBPixHits/F");
  outTree->Branch("lowMuon_NPixLayers",   &out_lowMuon_NPixLayers,        "lowMuon_NPixLayers/F");
  outTree->Branch("lowMuon_NTraLayers",   &out_lowMuon_NTraLayers,        "lowMuon_NTraLayers/F");
  outTree->Branch("lowMuon_NStrLayers",   &out_lowMuon_NStrLayers,        "lowMuon_NStrLayers/F");
  outTree->Branch("lowMuon_NBPixLayers",  &out_lowMuon_NBPixLayers,       "lowMuon_NBPixLayers/F");
  outTree->Branch("highMuon_isTight",     &out_highMuon_isTight,  "highMuon_isTight/F");
  outTree->Branch("highMuon_isLoose",     &out_highMuon_isLoose,  "highMuon_isLoose/F");
  outTree->Branch("highMuon_isSoft",      &out_highMuon_isSoft,   "highMuon_isSoft/F");
  outTree->Branch("highMuon_isMedium",    &out_highMuon_isMedium,         "highMuon_isMedium/F");
  outTree->Branch("highMuon_isHighPt",    &out_highMuon_isHighPt,         "highMuon_isHighPt/F");
  outTree->Branch("highMuon_isTracker",   &out_highMuon_isTracker,        "highMuon_isTracker/F");
  outTree->Branch("highMuon_isGlobal",    &out_highMuon_isGlobal,         "highMuon_isGlobal/F");
  outTree->Branch("highMuon_NPixelHits",  &out_highMuon_NPixelHits,       "highMuon_NPixelHits/F");
  outTree->Branch("highMuon_NStripHits",  &out_highMuon_NStripHits,       "highMuon_NStripHits/F");
  outTree->Branch("highMuon_NTrackhits",  &out_highMuon_NTrackhits,       "highMuon_NTrackhits/F");
  outTree->Branch("highMuon_NBPixHits",   &out_highMuon_NBPixHits,        "highMuon_NBPixHits/F");
  outTree->Branch("highMuon_NPixLayers",  &out_highMuon_NPixLayers,       "highMuon_NPixLayers/F");
  outTree->Branch("highMuon_NTraLayers",  &out_highMuon_NTraLayers,       "highMuon_NTraLayers/F");
  outTree->Branch("highMuon_NStrLayers",  &out_highMuon_NStrLayers,       "highMuon_NStrLayers/F");
  outTree->Branch("highMuon_NBPixLayers",         &out_highMuon_NBPixLayers,      "highMuon_NBPixLayers/F");
  outTree->Branch("lowMuon_type",         &out_lowMuon_type,      "lowMuon_type/F");
  outTree->Branch("highMuon_type",        &out_highMuon_type,     "highMuon_type/F");
  outTree->Branch("highTrack_NPixelHits",         &out_highTrack_NPixelHits,      "highTrack_NPixelHits/F");
  outTree->Branch("highTrack_NStripHits",         &out_highTrack_NStripHits,      "highTrack_NStripHits/F");
  outTree->Branch("highTrack_NTrackhits",         &out_highTrack_NTrackhits,      "highTrack_NTrackhits/F");
  outTree->Branch("highTrack_NBPixHits",  &out_highTrack_NBPixHits,       "highTrack_NBPixHits/F");
  outTree->Branch("highTrack_NPixLayers",         &out_highTrack_NPixLayers,      "highTrack_NPixLayers/F");
  outTree->Branch("highTrack_NTraLayers",         &out_highTrack_NTraLayers,      "highTrack_NTraLayers/F");
  outTree->Branch("highTrack_NStrLayers",         &out_highTrack_NStrLayers,      "highTrack_NStrLayers/F");
  outTree->Branch("highTrack_NBPixLayers",        &out_highTrack_NBPixLayers,     "highTrack_NBPixLayers/F");
  outTree->Branch("lowTrack_NPixelHits",  &out_lowTrack_NPixelHits,       "lowTrack_NPixelHits/F");
  outTree->Branch("lowTrack_NStripHits",  &out_lowTrack_NStripHits,       "lowTrack_NStripHits/F");
  outTree->Branch("lowTrack_NTrackhits",  &out_lowTrack_NTrackhits,       "lowTrack_NTrackhits/F");
  outTree->Branch("lowTrack_NBPixHits",   &out_lowTrack_NBPixHits,        "lowTrack_NBPixHits/F");
  outTree->Branch("lowTrack_NPixLayers",  &out_lowTrack_NPixLayers,       "lowTrack_NPixLayers/F");
  outTree->Branch("lowTrack_NTraLayers",  &out_lowTrack_NTraLayers,       "lowTrack_NTraLayers/F");
  outTree->Branch("lowTrack_NStrLayers",  &out_lowTrack_NStrLayers,       "lowTrack_NStrLayers/F");
  outTree->Branch("lowTrack_NBPixLayers",         &out_lowTrack_NBPixLayers,      "lowTrack_NBPixLayers/F");
  outTree->Branch("thirdTrack_NPixelHits",        &out_thirdTrack_NPixelHits,     "thirdTrack_NPixelHits/F");
  outTree->Branch("thirdTrack_NStripHits",        &out_thirdTrack_NStripHits,     "thirdTrack_NStripHits/F");
  outTree->Branch("thirdTrack_NTrackhits",        &out_thirdTrack_NTrackhits,     "thirdTrack_NTrackhits/F");
  outTree->Branch("thirdTrack_NBPixHits",         &out_thirdTrack_NBPixHits,      "thirdTrack_NBPixHits/F");
  outTree->Branch("thirdTrack_NPixLayers",        &out_thirdTrack_NPixLayers,     "thirdTrack_NPixLayers/F");
  outTree->Branch("thirdTrack_NTraLayers",        &out_thirdTrack_NTraLayers,     "thirdTrack_NTraLayers/F");
  outTree->Branch("thirdTrack_NStrLayers",        &out_thirdTrack_NStrLayers,     "thirdTrack_NStrLayers/F");
  outTree->Branch("thirdTrack_NBPixLayers",       &out_thirdTrack_NBPixLayers,    "thirdTrack_NBPixLayers/F");
  outTree->Branch("fourthTrack_NPixelHits",       &out_fourthTrack_NPixelHits,    "fourthTrack_NPixelHits/F");
  outTree->Branch("fourthTrack_NStripHits",       &out_fourthTrack_NStripHits,    "fourthTrack_NStripHits/F");
  outTree->Branch("fourthTrack_NTrackhits",       &out_fourthTrack_NTrackhits,    "fourthTrack_NTrackhits/F");
  outTree->Branch("fourthTrack_NBPixHits",        &out_fourthTrack_NBPixHits,     "fourthTrack_NBPixHits/F");
  outTree->Branch("fourthTrack_NPixLayers",       &out_fourthTrack_NPixLayers,    "fourthTrack_NPixLayers/F");
  outTree->Branch("fourthTrack_NTraLayers",       &out_fourthTrack_NTraLayers,    "fourthTrack_NTraLayers/F");
  outTree->Branch("fourthTrack_NStrLayers",       &out_fourthTrack_NStrLayers,    "fourthTrack_NStrLayers/F");
  outTree->Branch("fourthTrack_NBPixLayers",      &out_fourthTrack_NBPixLayers,   "fourthTrack_NBPixLayers/F");
  outTree->Branch("six_m_kkpp",   &out_six_m_kkpp,        "six_m_kkpp/F");
  outTree->Branch("six_m_ref_kkpp",       &out_six_m_ref_kkpp,    "six_m_ref_kkpp/F");
  outTree->Branch("six_mass_ppkk",        &out_six_mass_ppkk,     "six_mass_ppkk/F");
  outTree->Branch("six_mass_pkpk",        &out_six_mass_pkpk,     "six_mass_pkpk/F");
  outTree->Branch("six_mass_pppp",        &out_six_mass_pppp,     "six_mass_pppp/F");
  outTree->Branch("six_mass_kpkp",        &out_six_mass_kpkp,     "six_mass_kpkp/F");
  outTree->Branch("six_mass_kppk",        &out_six_mass_kppk,     "six_mass_kppk/F");
  outTree->Branch("six_mass_kkkk",        &out_six_mass_kkkk,     "six_mass_kkkk/F");
  outTree->Branch("six_pt",       &out_six_pt,    "six_pt/F");
  outTree->Branch("six_eta",      &out_six_eta,   "six_eta/F");
  outTree->Branch("six_phi",      &out_six_phi,   "six_phi/F");
  outTree->Branch("six_p",        &out_six_p,     "six_p/F");
  outTree->Branch("six_cosAlpha",         &out_six_cosAlpha,      "six_cosAlpha/F");
  outTree->Branch("six_ctauPV",   &out_six_ctauPV,        "six_ctauPV/F");
  outTree->Branch("six_ctauErrPV",        &out_six_ctauErrPV,     "six_ctauErrPV/F");
  outTree->Branch("six_cosAlphaCA",       &out_six_cosAlphaCA,    "six_cosAlphaCA/F");
  outTree->Branch("six_ctauPVCA",         &out_six_ctauPVCA,      "six_ctauPVCA/F");
  outTree->Branch("six_ctauErrPVCA",      &out_six_ctauErrPVCA,   "six_ctauErrPVCA/F");
  outTree->Branch("six_cosAlphaDZ",       &out_six_cosAlphaDZ,    "six_cosAlphaDZ/F");
  outTree->Branch("six_ctauPVDZ",         &out_six_ctauPVDZ,      "six_ctauPVDZ/F");
  outTree->Branch("six_ctauErrPVDZ",      &out_six_ctauErrPVDZ,   "six_ctauErrPVDZ/F");
  outTree->Branch("six_cosAlphaBS",       &out_six_cosAlphaBS,    "six_cosAlphaBS/F");
  outTree->Branch("six_ctauPVBS",         &out_six_ctauPVBS,      "six_ctauPVBS/F");
  outTree->Branch("six_ctauErrPVBS",      &out_six_ctauErrPVBS,   "six_ctauErrPVBS/F");
  outTree->Branch("six_vProb",    &out_six_vProb,         "six_vProb/F");
  outTree->Branch("six_nDof",     &out_six_nDof,  "six_nDof/F");
  outTree->Branch("six_vChi2",    &out_six_vChi2,         "six_vChi2/F");
  outTree->Branch("six_vx",       &out_six_vx,    "six_vx/F");
  outTree->Branch("six_vy",       &out_six_vy,    "six_vy/F");
  outTree->Branch("six_vz",       &out_six_vz,    "six_vz/F");
  outTree->Branch("six_charge",   &out_six_charge,        "six_charge/F");


}

Bool_t SixTracks::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  fReader.SetEntry(entry);

  bool test = true;

  test = test && ((*lowMuon_pt) >= 2.0);

  test = test && (*lowTrack_pt >= 0.7);

  //
  test = test && (*lowMuonMatch>0.0) && (*highMuonMatch>0.0);
  //
  test = test && (*dimuonditrk_vProb> 0.015);

  test = test && (fabs(*dimuonditrk_cosAlpha)> 0.90);

  test = test && ((*dimuon_pt)> 3.5);

  test = test && (*highTrack_NBPixHits>1) && (*lowTrack_NPixelHits>1);

  // test = test && ((*dimuonditrk_pt)> 6.0);

  //int a = (int) (*trigger);
  //std::cout << (*trigger);

  if(test)
  {

    Last login: Mon Jul  1 10:00:17 on ttys006
    (base) visitor-50230508:CNNFiltering adrianodif$ ssh -X ui-centos7.recas.ba.infn.it
    Last login: Mon Jul  1 11:33:16 2019 from 2001:1458:204:1::102:3fcc
    screen -r C
    Identity added: /lustre/home/adrianodif/.ssh/id_tesla (/lustre/home/adrianodif/.ssh/id_tesla)
    [adrianodif@ui-centos7 ~]$ screen -r C
    There is a screen on:
    	55914.COH	(Attached)
    There is no screen to be resumed matching C.
    [adrianodif@ui-centos7 ~]$ screen -D
    [55914.COH power detached.]

    [adrianodif@ui-centos7 ~]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 ~]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 ~]$ cd /lustre/cms/store/user/adiflori/GEN-MINIAODSIMQCD_PhiFilter_Dau_SoftQCD_MiniTracks_/crab_GEN-MINIAODSIM_QCD_PhiFilter_Dau_SoftQCD_MiniTracks__20190701_061305/190701_041401/0001/qcd_ml
    [adrianodif@ui-centos7 qcd_ml]$ ls
    ^C
    [adrianodif@ui-centos7 qcd_ml]$ vi to^C
    [adrianodif@ui-centos7 qcd_ml]$ vi toEvent.py
    [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
    [adrianodif@ui-centos7 qcd_ml]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
    [adrianodif@ui-centos7 qcd_ml]$ screen -r C^C
    [adrianodif@ui-centos7 qcd_ml]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 qcd_ml]$ vu ^C
    [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
    [adrianodif@ui-centos7 qcd_ml]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
    [adrianodif@ui-centos7 qcd_ml]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
    [adrianodif@ui-centos7 qcd_ml]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
    [adrianodif@ui-centos7 qcd_ml]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
    [adrianodif@ui-centos7 qcd_ml]$ rm *pion*
    ^C
    [adrianodif@ui-centos7 qcd_ml]$ rm *electron* &
    [1] 32315
    [adrianodif@ui-centos7 qcd_ml]$ rm *muon*
    ^C
    [adrianodif@ui-centos7 qcd_ml]$ rm *muon* &
    [2] 34149
    [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
    [1]-  Done                    rm *electron*
    [adrianodif@ui-centos7 qcd_ml]$
    [adrianodif@ui-centos7 qcd_ml]$ screen -r C
    [detached from 55914.COH]
    [2]+  Done                    rm *muon*
    [adrianodif@ui-centos7 qcd_ml]$ vi to^C
    [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.p^C
    [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
    [adrianodif@ui-centos7 qcd_ml]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
    [adrianodif@ui-centos7 qcd_ml]$ ps
       PID TTY          TIME CMD
      8843 pts/48   00:00:01 bash
     64411 pts/48   00:00:00 ps
    [adrianodif@ui-centos7 qcd_ml]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
    [adrianodif@ui-centos7 qcd_ml]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
    ^[OA[adrianodif@ui-centos7 qcd_ml]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 qcd_ml]$ pwd
    /lustre/cms/store/user/adiflori/GEN-MINIAODSIMQCD_PhiFilter_Dau_SoftQCD_MiniTracks_/crab_GEN-MINIAODSIM_QCD_PhiFilter_Dau_SoftQCD_MiniTracks__20190701_061305/190701_041401/0001/qcd_ml
    [adrianodif@ui-centos7 qcd_ml]$ cd -
    /lustre/home/adrianodif
    [adrianodif@ui-centos7 ~]$ cd -
    /lustre/cms/store/user/adiflori/GEN-MINIAODSIMQCD_PhiFilter_Dau_SoftQCD_MiniTracks_/crab_GEN-MINIAODSIM_QCD_PhiFilter_Dau_SoftQCD_MiniTracks__20190701_061305/190701_041401/0001/qcd_ml
    [adrianodif@ui-centos7 qcd_ml]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
    [adrianodif@ui-centos7 qcd_ml]$ ^C
    [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
    c[adrianodif@ui-centos7 qcd_ml]$ screen -^C
    [adrianodif@ui-centos7 qcd_ml]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 qcd_ml]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 qcd_ml]$ screen -r ^C
    [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
    [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
    [adrianodif@ui-centos7 qcd_ml]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
    [adrianodif@ui-centos7 qcd_ml]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
    [adrianodif@ui-centos7 qcd_ml]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
    [adrianodif@ui-centos7 qcd_ml]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
    [adrianodif@ui-centos7 qcd_ml]$ vi toPdgs.py
    [adrianodif@ui-centos7 qcd_ml]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 qcd_ml]$ cd /lustre/cms/store/user/adiflori/GEN-MINIAODSIMQCD_PhiFilter_Dau_SoftQCD_MiniTracks_/crab_GEN-MINIAODSIM_QCD_PhiFilter_Dau_SoftQCD_MiniTracks__20190701_061305/190701_041401/0000/qcd_ml
    [adrianodif@ui-centos7 qcd_ml]$ vi pdg_0
    ^C^C^C^C^[zz^Z
    [adrianodif@ui-centos7 qcd_ml]$ cd data
    vi pdg_[adrianodif@ui-centos7 data]$ vi pdg_3
    [adrianodif@ui-centos7 data]$ vi pdg_4
    [adrianodif@ui-centos7 data]$ vi toPdgs.py^C
    [adrianodif@ui-centos7 data]$ screen -r C
    [detached from 55914.COH]
    (reverse-i-search)`v': ^C pdg_4
    [adrianodif@ui-centos7 data]$ vi toP^C
    [adrianodif@ui-centos7 data]$ vi toPdg^C
    [adrianodif@ui-centos7 data]$ vi toPdgs.py
    [adrianodif@ui-centos7 data]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 data]$ vi toP^C
    [adrianodif@ui-centos7 data]$ vi toPdgs.py
    [adrianodif@ui-centos7 data]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 data]$ vi to^C
    [adrianodif@ui-centos7 data]$ vi toPdgs.py
    [adrianodif@ui-centos7 data]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 data]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 data]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 data]$ cd /lustre/cms/store/user/adiflori/GEN-MINIAODSIMQCD_PhiFilter_Dau_SoftQCD_MiniTracks_/crab_GEN-MINIAODSIM_QCD_PhiFilter_Dau_SoftQCD_MiniTracks__20190701_061305/190701_041401/0001/qcd_ml
    [adrianodif@ui-centos7 qcd_ml]$ cp toPdgs.py toElectron.py
    [adrianodif@ui-centos7 qcd_ml]$ vi toElectron
    [adrianodif@ui-centos7 qcd_ml]$ vi toElectron.py^C
    [adrianodif@ui-centos7 qcd_ml]$ vi toElectron.py
    (reverse-i-search)`for': ^Cr f in *; do tar xvf $f -P; done
    [adrianodif@ui-centos7 qcd_ml]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 qcd_ml]$ ls -Slh *electron* | tail -10
    ^C
    [adrianodif@ui-centos7 qcd_ml]$ screen -r C
    [detached from 55914.COH]
    [adrianodif@ui-centos7 qcd_ml]$ vi toElectron.py
    [adrianodif@ui-centos7 qcd_ml]$ char^C
    [adrianodif@ui-centos7 qcd_ml]$ mystore
    [adrianodif@ui-centos7 adiflori]$ cd Charmonium/
    [adrianodif@ui-centos7 Charmonium]$ ls -tlrh
    totale 416K
    -rw-rw-r--+ 1 adrianodif cms 6,1K 14 set  2018 sumlist.txt
    -rw-rw-r--+ 1 adrianodif cms 6,0K 24 set  2018 2mu2mulist.txt
    -rw-rw-r--+ 1 adrianodif cms  24K 26 set  2018 FourMuons.h
    -rw-rw-r--+ 1 adrianodif cms 3,1K 26 set  2018 FourMuons.C
    drwxrwxr-x+ 2 adrianodif cms 4,0K  9 ott  2018 2mu2mu
    -rw-rw-r--+ 1 adrianodif cms 3,5K 27 ott  2018 2017data2mu2k
    -rw-rw-r--+ 1 adrianodif cms 177K 27 ott  2018 log
    -rw-rw-r--+ 1 adrianodif cms 7,6K 28 ott  2018 2016data
    -rw-rw-r--+ 1 adrianodif cms  847 19 nov  2018 mclist.sh
    -rw-rw-r--+ 1 adrianodif cms 6,2K 23 nov  2018 2017data
    -rw-rw-r--+ 1 adrianodif cms 4,2K 28 nov  2018 2018data
    drwxrwxr-x+ 3        497 497 4,0K 23 mar 11.27 crab_miniaod_2mu2k_five_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20190323_105234
    -rw-rw-r--+ 1 adrianodif cms  29K  9 apr 18.12 dummy.h
    -rw-rw-r--+ 1 adrianodif cms 3,1K  9 apr 18.12 dummy.C
    drwxrwxr-x+ 3        497 497 4,0K 26 giu 18.50 crab_miniaod_2mu2k_Charmonium_Run2018B-17Sep2018-v1_MINIAOD___20190626_182314
    drwxrwxr-x+ 3        497 497 4,0K 27 giu 09.02 crab_miniaod_2mu2k_Charmonium_Run2018C-17Sep2018-v1_MINIAOD___20190627_084018
    drwxrwxr-x+ 3        497 497 4,0K 27 giu 10.25 crab_miniaod_2mu2k_Charmonium_Run2018A-17Sep2018-v1_MINIAOD___20190627_095241
    drwxrwxr-x+ 3        497 497 4,0K 27 giu 11.51 crab_miniaod_2mu2k_Charmonium_Run2018C-17Sep2018-v1_MINIAOD___20190627_112903_six_five
    drwxrwxr-x+ 3        497 497 4,0K 27 giu 12.37 crab_miniaod_2mu2k_Charmonium_Run2018B-17Sep2018-v1_MINIAOD___20190627_112839_six_five
    drwxrwxr-x+ 3        497 497 4,0K 27 giu 13.30 crab_miniaod_2mu2k_Charmonium_Run2018D-PromptReco-v2_MINIAOD___20190627_125057_six_five
    drwxrwxr-x+ 3        497 497 4,0K 27 giu 13.32 crab_miniaod_2mu2k_Charmonium_Run2018A-17Sep2018-v1_MINIAOD___20190627_124818_six_five
    [adrianodif@ui-centos7 Charmonium]$ cd crab_miniaod_2mu2k_Charmonium_Run2018A-17Sep2018-v1_MINIAOD___20190627_124818_six_five/
    [adrianodif@ui-centos7 crab_miniaod_2mu2k_Charmonium_Run2018A-17Sep2018-v1_MINIAOD___20190627_124818_six_five]$ ls -tlrh
    totale 0
    drwxrwxr-x+ 7 497 497 4,0K 30 giu 14.05 190627_104825
    [adrianodif@ui-centos7 crab_miniaod_2mu2k_Charmonium_Run2018A-17Sep2018-v1_MINIAOD___20190627_124818_six_five]$ cd 190627_104825/
    [adrianodif@ui-centos7 190627_104825]$ ls -tlrh
    totale 110G
    drwxrwxr-x+ 2        497 497 128K 28 giu 16.04 0002
    drwxrwxr-x+ 2        497 497 128K 28 giu 16.04 0001
    drwxrwxr-x+ 2        497 497 128K 28 giu 16.04 0003
    drwxrwxr-x+ 2        497 497  96K 28 giu 16.04 0004
    drwxrwxr-x+ 2        497 497 128K 28 giu 16.04 0000
    -rw-r--r--+ 1 adrianodif cms  17G 28 giu 16.34 0000.root
    -rw-r--r--+ 1 adrianodif cms  22G 28 giu 16.38 0004.root
    -rw-r--r--+ 1 adrianodif cms  23G 28 giu 16.39 0001.root
    -rw-r--r--+ 1 adrianodif cms  23G 28 giu 16.41 0003.root
    -rw-r--r--+ 1 adrianodif cms  26G 28 giu 16.44 0002.root
    -rw-r-----+ 1 adrianodif cms 1,6K 28 giu 17.00 toHdF.py
    -rw-rw-r--+ 1 adrianodif cms  137 28 giu 19.20 DU
    -rw-rw-r--+ 1 adrianodif cms  34K 29 giu 09.38 SixTracks.h
    -rw-rw-r--+ 1 adrianodif cms 3,1K 29 giu 09.38 SixTracks.C
    [adrianodif@ui-centos7 190627_104825]$ pwd
    /lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018A-17Sep2018-v1_MINIAOD___20190627_124818_six_five/190627_104825
    [adrianodif@ui-centos7 190627_104825]$ cd --
    [adrianodif@ui-centos7 ~]$ cd jpsiphi/2018/data_2018/analysis/utilities/skimmers/
    [adrianodif@ui-centos7 skimmers]$ ls -tlrh
    totale 128K
    drwxr-x---+  2 adrianodif cms 4,0K 13 set  2018 2muBkg
    drwxr-x---+  2 adrianodif cms 4,0K 13 set  2018 2mukpi
    drwxr-x---+  2 adrianodif cms 4,0K 13 set  2018 2mupik
    drwxr-x---+  2 adrianodif cms 4,0K 13 set  2018 4mu
    -rw-r-----+  1 adrianodif cms 3,1K 13 set  2018 FourMuSkim.C
    drwxr-x---+  2 adrianodif cms 4,0K 13 set  2018 legacy
    -rw-r-----+  1 adrianodif cms  54K 13 set  2018 skimRunII_xmass.C
    -rw-r-----+  1 adrianodif cms 5,1K 13 set  2018 skimRunII_xmass.ipynb
    -rw-r-----+  1 adrianodif cms 1,9K 13 set  2018 skimRunII_xmass.py
    drwxr-x---+  2 adrianodif cms 4,0K 13 set  2018 w_hlts
    drwxr-x---+  2 adrianodif cms 4,0K 13 set  2018 wo_hlts
    -rw-r-----+  1 adrianodif cms 1,4K 25 set  2018 skimvars.sh
    drwxr-x---+  2 adrianodif cms 4,0K 27 ott  2018 2mu2mu
    drwxr-x---+  3 adrianodif cms 4,0K 31 ott  2018 splotskimmers
    drwxr-x---+  2 adrianodif cms 4,0K 21 nov  2018 2mu2k_five
    drwxr-x---+ 12 adrianodif cms 4,0K 23 mar 16.34 2mu2k
    drwxr-x---+  2 adrianodif cms 4,0K 23 mar 16.43 2mu4trk
    -rw-r-----+  1 adrianodif cms  336  6 apr 23.54 DU
    drwxr-x---+  8 adrianodif cms 4,0K  6 apr 23.55 2mu3trk
    -rw-r-----+  1 adrianodif cms  11K 10 apr 15.45 merging.py
    drwxr-x---+  3 adrianodif cms 4,0K 11 apr 04.03 fourtracks
    drwxr-x---+  3 adrianodif cms 4,0K 11 apr 04.05 fivetracks
    drwxr-x---+  5 adrianodif cms 4,0K 11 apr 10.11 sixtracks
    -rw-r-----+  1 adrianodif cms 1,6K 22 mag 00.38 toHdF.py
    drwxr-x---+  2 adrianodif cms 4,0K 29 mag 07.28 2012
    drwxr-x---+  2 adrianodif cms 4,0K 29 giu 09.40 sixtracks_new
    [adrianodif@ui-centos7 skimmers]$ cd sixtracks_new/
    [adrianodif@ui-centos7 sixtracks_new]$ ls -tlrh
    totale 0
    [adrianodif@ui-centos7 sixtracks_new]$ ls -tlr^C
    [adrianodif@ui-centos7 sixtracks_new]$ cp /lustre/cms/store/user/adiflori/Charmonium/crab_miniaod_2mu2k_Charmonium_Run2018A-17Sep2018-v1_MINIAOD___20190627_124818_six_five/190627_104825/SixTracks.* .
    [adrianodif@ui-centos7 sixtracks_new]$ cp ../skimvars.sh .
    [adrianodif@ui-centos7 sixtracks_new]$ . skimvars.sh SixTracks.h
    Use this script with the header generated by MakeSelector
    sed: espressione -e #1, carattere 2: ci sono altri caratteri dopo il comando
    [adrianodif@ui-centos7 sixtracks_new]$ vi skimvars.sh
    [adrianodif@ui-centos7 sixtracks_new]$ vi lis^C
    [adrianodif@ui-centos7 sixtracks_new]$ ls
    assign.txt  branches.txt  SixTracks.C  SixTracks.h  skimvars.sh  vars.txt
    [adrianodif@ui-centos7 sixtracks_new]$ vi assign.txt
    [adrianodif@ui-centos7 sixtracks_new]$ vi SixTracks.h









































































































































































































































































    [adrianodif@ui-centos7 sixtracks_new]$ vi SixTracks.C
    [adrianodif@ui-centos7 sixtracks_new]$ vi assign.txt
    [adrianodif@ui-centos7 sixtracks_new]$ ls
    assign.txt  branches.txt  SixTracks.C  SixTracks.h  skimvars.sh  vars.txt
    [adrianodif@ui-centos7 sixtracks_new]$ vi branches.txt
    [adrianodif@ui-centos7 sixtracks_new]$ vi assign.txt

    out_run =       (Float_t)(*run);
    out_event =     (Float_t)(*event);
    out_lumi =      (Float_t)(*lumi);
    out_numPrimaryVertices =        (Float_t)(*numPrimaryVertices);
    out_trigger =   (Float_t)(*trigger);
    out_noSixCandidates =   (Float_t)(*noSixCandidates);
    out_five_id =   (Float_t)(*five_id);
    out_dimuon_id =         (Float_t)(*dimuon_id);
    out_p1_id =     (Float_t)(*p1_id);
    out_m1_id =     (Float_t)(*m1_id);
    out_p2_id =     (Float_t)(*p2_id);
    out_m2_id =     (Float_t)(*m2_id);
    out_swapped =   (Float_t)(*swapped);
    out_sameSign =  (Float_t)(*sameSign);
    out_sameSign_mmtt =     (Float_t)(*sameSign_mmtt);
    out_six_p4 =    (Float_t)(*six_p4);
    out_five_p4 =   (Float_t)(*five_p4);
    out_mumukk_p4 =         (Float_t)(*mumukk_p4);
    out_ditrack_p4 =        (Float_t)(*ditrack_p4);
    out_dimuon_p4 =         (Float_t)(*dimuon_p4);
    out_lowMuon_p4 =        (Float_t)(*lowMuon_p4);
    out_highMuon_p4 =       (Float_t)(*highMuon_p4);
    out_highKaon_p4 =       (Float_t)(*highKaon_p4);
    out_lowKaon_p4 =        (Float_t)(*lowKaon_p4);
    out_thirdKaon_p4 =      (Float_t)(*thirdKaon_p4);
    out_fourthKaon_p4 =     (Float_t)(*fourthKaon_p4);
    out_highPion_p4 =       (Float_t)(*highPion_p4);
    out_lowPion_p4 =        (Float_t)(*lowPion_p4);
    out_thirdPion_p4 =      (Float_t)(*thirdPion_p4);
    out_fourthPion_p4 =     (Float_t)(*fourthPion_p4);
    out_highProton_p4 =     (Float_t)(*highProton_p4);
    out_lowProton_p4 =      (Float_t)(*lowProton_p4);
    out_thirdProton_p4 =    (Float_t)(*thirdProton_p4);
    out_fourthProton_p4 =   (Float_t)(*fourthProton_p4);
    out_mumukk_m =  (Float_t)(*mumukk_m);
    out_mumukk_pt =         (Float_t)(*mumukk_pt);
    out_mumukk_eta =        (Float_t)(*mumukk_eta);
    out_mumukk_phi =        (Float_t)(*mumukk_phi);
    out_mumukk_p =  (Float_t)(*mumukk_p);
    out_dimuon_m =  (Float_t)(*dimuon_m);
    out_dimuon_pt =         (Float_t)(*dimuon_pt);
    out_dimuon_eta =        (Float_t)(*dimuon_eta);
    out_dimuon_phi =        (Float_t)(*dimuon_phi);
    out_dimuon_p =  (Float_t)(*dimuon_p);
    out_highTrackMatch =    (Float_t)(*highTrackMatch);
    out_lowTrackMatch =     (Float_t)(*lowTrackMatch);
    out_lowMuonMatch =      (Float_t)(*lowMuonMatch);
    out_highMuonMatch =     (Float_t)(*highMuonMatch);
    out_thirdTrackMatch =   (Float_t)(*thirdTrackMatch);
    out_fourthTrackMatch =  (Float_t)(*fourthTrackMatch);
    out_ditrack_m =         (Float_t)(*ditrack_m);
    out_diTrackOne_pt =     (Float_t)(*diTrackOne_pt);
    out_diTrackOne_eta =    (Float_t)(*diTrackOne_eta);
    out_diTrackOne_phi =    (Float_t)(*diTrackOne_phi);
    out_diTrackOne_p =      (Float_t)(*diTrackOne_p);
    out_diTrackTwo_pt =     (Float_t)(*diTrackTwo_pt);
    out_diTrackTwo_eta =    (Float_t)(*diTrackTwo_eta);
    out_diTrackTwo_phi =    (Float_t)(*diTrackTwo_phi);
    out_diTrackTwo_p =      (Float_t)(*diTrackTwo_p);
    out_diTrackThree_pt =   (Float_t)(*diTrackThree_pt);
    out_diTrackThree_eta =  (Float_t)(*diTrackThree_eta);
    out_diTrackThree_phi =  (Float_t)(*diTrackThree_phi);
    out_diTrackThree_p =    (Float_t)(*diTrackThree_p);
    out_diTrackFour_pt =    (Float_t)(*diTrackFour_pt);
    out_diTrackFour_eta =   (Float_t)(*diTrackFour_eta);
    out_diTrackFour_phi =   (Float_t)(*diTrackFour_phi);
    out_diTrackFour_p =     (Float_t)(*diTrackFour_p);
    out_diTrackFive_pt =    (Float_t)(*diTrackFive_pt);
    out_diTrackFive_eta =   (Float_t)(*diTrackFive_eta);
    out_diTrackFive_phi =   (Float_t)(*diTrackFive_phi);
    out_diTrackFive_p =     (Float_t)(*diTrackFive_p);
    out_diTrackSix_pt =     (Float_t)(*diTrackSix_pt);
    out_diTrackSix_eta =    (Float_t)(*diTrackSix_eta);
    out_diTrackSix_phi =    (Float_t)(*diTrackSix_phi);
    out_diTrackSix_p =      (Float_t)(*diTrackSix_p);
    out_dimuonDiTrkOne_mmpp =       (Float_t)(*dimuonDiTrkOne_mmpp);
    out_dimuonDiTrkTwo_mmpp =       (Float_t)(*dimuonDiTrkTwo_mmpp);
    out_dimuonDiTrkThree_mmpp =     (Float_t)(*dimuonDiTrkThree_mmpp);
    out_dimuonDiTrkFour_mmpp =      (Float_t)(*dimuonDiTrkFour_mmpp);
    out_dimuonDiTrkOne_mmkk =       (Float_t)(*dimuonDiTrkOne_mmkk);
    out_dimuonDiTrkTwo_mmkk =       (Float_t)(*dimuonDiTrkTwo_mmkk);
    out_dimuonDiTrkThree_mmkk =     (Float_t)(*dimuonDiTrkThree_mmkk);
    out_dimuonDiTrkFour_mmkk =      (Float_t)(*dimuonDiTrkFour_mmkk);
    out_dimuonDiTrkOne_mmpk =       (Float_t)(*dimuonDiTrkOne_mmpk);
    out_dimuonDiTrkTwo_mmpk =       (Float_t)(*dimuonDiTrkTwo_mmpk);
    out_dimuonDiTrkThree_mmpk =     (Float_t)(*dimuonDiTrkThree_mmpk);
    out_dimuonDiTrkFour_mmpk =      (Float_t)(*dimuonDiTrkFour_mmpk);
    out_dimuonDiTrkOne_mmkp =       (Float_t)(*dimuonDiTrkOne_mmkp);
    out_dimuonDiTrkTwo_mmkp =       (Float_t)(*dimuonDiTrkTwo_mmkp);
    out_dimuonDiTrkThree_mmkp =     (Float_t)(*dimuonDiTrkThree_mmkp);
    out_dimuonDiTrkFour_mmkp =      (Float_t)(*dimuonDiTrkFour_mmkp);
    out_diTrackOne_kk =     (Float_t)(*diTrackOne_kk);
    out_diTrackTwo_kk =     (Float_t)(*diTrackTwo_kk);
    out_diTrackThree_kk =   (Float_t)(*diTrackThree_kk);
    out_diTrackFour_kk =    (Float_t)(*diTrackFour_kk);
    out_diTrackFive_kk =    (Float_t)(*diTrackFive_kk);
    out_diTrackSix_kk =     (Float_t)(*diTrackSix_kk);
    out_diTrackOne_pp =     (Float_t)(*diTrackOne_pp);
    out_diTrackTwo_pp =     (Float_t)(*diTrackTwo_pp);
    out_diTrackThree_pp =   (Float_t)(*diTrackThree_pp);
    out_diTrackFour_pp =    (Float_t)(*diTrackFour_pp);
    out_diTrackFive_pp =    (Float_t)(*diTrackFive_pp);
    out_diTrackSix_pp =     (Float_t)(*diTrackSix_pp);
    out_diTrackOne_pk =     (Float_t)(*diTrackOne_pk);
    out_diTrackTwo_pk =     (Float_t)(*diTrackTwo_pk);
    out_diTrackThree_pk =   (Float_t)(*diTrackThree_pk);
    out_diTrackFour_pk =    (Float_t)(*diTrackFour_pk);
    out_diTrackFive_pk =    (Float_t)(*diTrackFive_pk);
    out_diTrackSix_pk =     (Float_t)(*diTrackSix_pk);
    out_diTrackOne_kp =     (Float_t)(*diTrackOne_kp);
    out_diTrackTwo_kp =     (Float_t)(*diTrackTwo_kp);
    out_diTrackThree_kp =   (Float_t)(*diTrackThree_kp);
    out_diTrackFour_kp =    (Float_t)(*diTrackFour_kp);
    out_diTrackFive_kp =    (Float_t)(*diTrackFive_kp);
    out_diTrackSix_kp =     (Float_t)(*diTrackSix_kp);
    out_highMuon_pt =       (Float_t)(*highMuon_pt);
    out_highMuon_eta =      (Float_t)(*highMuon_eta);
    out_highMuon_phi =      (Float_t)(*highMuon_phi);
    out_highMuon_charge =   (Float_t)(*highMuon_charge);
    out_highMuon_dz =       (Float_t)(*highMuon_dz);
    out_highMuon_dxy =      (Float_t)(*highMuon_dxy);
    out_lowMuon_pt =        (Float_t)(*lowMuon_pt);
    out_lowMuon_eta =       (Float_t)(*lowMuon_eta);
    out_lowMuon_phi =       (Float_t)(*lowMuon_phi);
    out_lowMuon_charge =    (Float_t)(*lowMuon_charge);
    out_lowMuon_dz =        (Float_t)(*lowMuon_dz);
    out_lowMuon_dxy =       (Float_t)(*lowMuon_dxy);
    out_highTrack_pt =      (Float_t)(*highTrack_pt);
    out_highTrack_eta =     (Float_t)(*highTrack_eta);
    out_highTrack_phi =     (Float_t)(*highTrack_phi);
    out_highTrack_charge =  (Float_t)(*highTrack_charge);
    out_highTrack_dz =      (Float_t)(*highTrack_dz);
    out_highTrack_dxy =     (Float_t)(*highTrack_dxy);
    out_lowTrack_pt =       (Float_t)(*lowTrack_pt);
    out_lowTrack_eta =      (Float_t)(*lowTrack_eta);
    out_lowTrack_phi =      (Float_t)(*lowTrack_phi);
    out_lowTrack_charge =   (Float_t)(*lowTrack_charge);
    out_lowTrack_dz =       (Float_t)(*lowTrack_dz);
    out_lowTrack_dxy =      (Float_t)(*lowTrack_dxy);
    out_thirdTrack_pt =     (Float_t)(*thirdTrack_pt);
    out_thirdTrack_eta =    (Float_t)(*thirdTrack_eta);
    out_thirdTrack_phi =    (Float_t)(*thirdTrack_phi);
    out_thirdTrack_charge =         (Float_t)(*thirdTrack_charge);
    out_thirdTrack_dz =     (Float_t)(*thirdTrack_dz);
    out_thirdTrack_dxy =    (Float_t)(*thirdTrack_dxy);
    out_dimuonDiTrkOne_pt =         (Float_t)(*dimuonDiTrkOne_pt);
    out_dimuonDiTrkOne_eta =        (Float_t)(*dimuonDiTrkOne_eta);
    out_dimuonDiTrkOne_phi =        (Float_t)(*dimuonDiTrkOne_phi);
    out_dimuonDiTrkOne_charge =     (Float_t)(*dimuonDiTrkOne_charge);
    out_dimuonDiTrkOne_p =  (Float_t)(*dimuonDiTrkOne_p);
    out_dimuonDiTrkTwo_pt =         (Float_t)(*dimuonDiTrkTwo_pt);
    out_dimuonDiTrkTwo_eta =        (Float_t)(*dimuonDiTrkTwo_eta);
    out_dimuonDiTrkTwo_phi =        (Float_t)(*dimuonDiTrkTwo_phi);
    out_dimuonDiTrkTwo_charge =     (Float_t)(*dimuonDiTrkTwo_charge);
    out_dimuonDiTrkTwo_p =  (Float_t)(*dimuonDiTrkTwo_p);
    out_dimuonDiTrkThree_pt =       (Float_t)(*dimuonDiTrkThree_pt);
    out_dimuonDiTrkThree_eta =      (Float_t)(*dimuonDiTrkThree_eta);
    out_dimuonDiTrkThree_phi =      (Float_t)(*dimuonDiTrkThree_phi);
    out_dimuonDiTrkThree_charge =   (Float_t)(*dimuonDiTrkThree_charge);
    out_dimuonDiTrkThree_p =        (Float_t)(*dimuonDiTrkThree_p);
    out_dimuonDiTrkFour_pt =        (Float_t)(*dimuonDiTrkFour_pt);
    out_dimuonDiTrkFour_eta =       (Float_t)(*dimuonDiTrkFour_eta);
    out_dimuonDiTrkFour_phi =       (Float_t)(*dimuonDiTrkFour_phi);
    out_dimuonDiTrkFour_charge =    (Float_t)(*dimuonDiTrkFour_charge);
    out_dimuonDiTrkFour_p =         (Float_t)(*dimuonDiTrkFour_p);
    out_dimuonDiTrkFive_pt =        (Float_t)(*dimuonDiTrkFive_pt);
    out_dimuonDiTrkFive_eta =       (Float_t)(*dimuonDiTrkFive_eta);
    out_dimuonDiTrkFive_phi =       (Float_t)(*dimuonDiTrkFive_phi);
    out_dimuonDiTrkFive_charge =    (Float_t)(*dimuonDiTrkFive_charge);
    out_dimuonDiTrkFive_p =         (Float_t)(*dimuonDiTrkFive_p);
    out_dimuonDiTrkSix_pt =         (Float_t)(*dimuonDiTrkSix_pt);
    out_dimuonDiTrkSix_eta =        (Float_t)(*dimuonDiTrkSix_eta);
    out_dimuonDiTrkSix_phi =        (Float_t)(*dimuonDiTrkSix_phi);
    out_dimuonDiTrkSix_charge =     (Float_t)(*dimuonDiTrkSix_charge);
    out_dimuonDiTrkSix_p =  (Float_t)(*dimuonDiTrkSix_p);
    out_dimuon_vProb =      (Float_t)(*dimuon_vProb);
    out_dimuon_vChi2 =      (Float_t)(*dimuon_vChi2);
    out_dimuon_DCA =        (Float_t)(*dimuon_DCA);
    out_dimuon_ctauPV =     (Float_t)(*dimuon_ctauPV);
    out_dimuon_ctauErrPV =  (Float_t)(*dimuon_ctauErrPV);
    out_dimuon_cosAlpha =   (Float_t)(*dimuon_cosAlpha);
    out_triTrackOne_kkk =   (Float_t)(*triTrackOne_kkk);
    out_triTrackTwo_kkk =   (Float_t)(*triTrackTwo_kkk);
    out_triTrackThree_kkk =         (Float_t)(*triTrackThree_kkk);
    out_triTrackFour_kkk =  (Float_t)(*triTrackFour_kkk);
    out_triTrackOne_kkp =   (Float_t)(*triTrackOne_kkp);
    out_triTrackTwo_kkp =   (Float_t)(*triTrackTwo_kkp);
    out_triTrackThree_kkp =         (Float_t)(*triTrackThree_kkp);
    out_triTrackOne_kpp =   (Float_t)(*triTrackOne_kpp);
    out_triTrackTwo_kpp =   (Float_t)(*triTrackTwo_kpp);
    out_triTrackThree_kpp =         (Float_t)(*triTrackThree_kpp);
    out_triTrackOne_ppp =   (Float_t)(*triTrackOne_ppp);
    out_triTrackTwo_ppp =   (Float_t)(*triTrackTwo_ppp);
    out_triTrackThree_ppp =         (Float_t)(*triTrackThree_ppp);
    out_triTrackOne_pt =    (Float_t)(*triTrackOne_pt);
    out_triTrackOne_eta =   (Float_t)(*triTrackOne_eta);
    out_triTrackOne_phi =   (Float_t)(*triTrackOne_phi);
    out_triTrackOne_charge =        (Float_t)(*triTrackOne_charge);
    out_triTrackTwo_pt =    (Float_t)(*triTrackTwo_pt);
    out_triTrackTwo_eta =   (Float_t)(*triTrackTwo_eta);
    out_triTrackTwo_phi =   (Float_t)(*triTrackTwo_phi);
    out_triTrackTwo_charge =        (Float_t)(*triTrackTwo_charge);
    out_triTrackThree_pt =  (Float_t)(*triTrackThree_pt);
    out_triTrackThree_eta =         (Float_t)(*triTrackThree_eta);
    out_triTrackThree_phi =         (Float_t)(*triTrackThree_phi);
    out_triTrackThree_charge =      (Float_t)(*triTrackThree_charge);
    out_triTrackFour_pt =   (Float_t)(*triTrackFour_pt);
    out_triTrackFour_eta =  (Float_t)(*triTrackFour_eta);
    out_triTrackFour_phi =  (Float_t)(*triTrackFour_phi);
    out_triTrackFour_charge =       (Float_t)(*triTrackFour_charge);
    out_mumukk_vProb =      (Float_t)(*mumukk_vProb);
    out_mumukk_vChi2 =      (Float_t)(*mumukk_vChi2);
    out_mumukk_nDof =       (Float_t)(*mumukk_nDof);
    out_mumukk_charge =     (Float_t)(*mumukk_charge);
    out_mumukk_cosAlpha =   (Float_t)(*mumukk_cosAlpha);
    out_mumukk_ctauPV =     (Float_t)(*mumukk_ctauPV);
    out_mumukk_ctauErrPV =  (Float_t)(*mumukk_ctauErrPV);
    out_mumukk_cosAlphaCA =         (Float_t)(*mumukk_cosAlphaCA);
    out_mumukk_ctauPVCA =   (Float_t)(*mumukk_ctauPVCA);
    out_mumukk_ctauErrPVCA =        (Float_t)(*mumukk_ctauErrPVCA);
    out_mumukk_cosAlphaDZ =         (Float_t)(*mumukk_cosAlphaDZ);
    out_mumukk_ctauPVDZ =   (Float_t)(*mumukk_ctauPVDZ);
    out_mumukk_ctauErrPVDZ =        (Float_t)(*mumukk_ctauErrPVDZ);
    out_mumukk_cosAlphaBS =         (Float_t)(*mumukk_cosAlphaBS);
    out_mumukk_ctauPVBS =   (Float_t)(*mumukk_ctauPVBS);
    out_mumukk_ctauErrPVBS =        (Float_t)(*mumukk_ctauErrPVBS);
    out_mumukk_vx =         (Float_t)(*mumukk_vx);
    out_mumukk_vy =         (Float_t)(*mumukk_vy);
    out_mumukk_vz =         (Float_t)(*mumukk_vz);
    out_dca_m1m2 =  (Float_t)(*dca_m1m2);
    out_dca_m1t1 =  (Float_t)(*dca_m1t1);
    out_dca_m1t2 =  (Float_t)(*dca_m1t2);
    out_dca_m2t1 =  (Float_t)(*dca_m2t1);
    out_dca_m2t2 =  (Float_t)(*dca_m2t2);
    out_dca_t1t2 =  (Float_t)(*dca_t1t2);
    out_dca_m1t3 =  (Float_t)(*dca_m1t3);
    out_dca_m2t3 =  (Float_t)(*dca_m2t3);
    out_dca_t1t3 =  (Float_t)(*dca_t1t3);
    out_dca_t2t3 =  (Float_t)(*dca_t2t3);
    out_dca_m1t4 =  (Float_t)(*dca_m1t4);
    out_dca_m2t4 =  (Float_t)(*dca_m2t4);
    out_dca_t1t4 =  (Float_t)(*dca_t1t4);
    out_dca_t2t4 =  (Float_t)(*dca_t2t4);
    out_dca_t3t4 =  (Float_t)(*dca_t3t4);
    out_highTrackMuonDR =   (Float_t)(*highTrackMuonDR);
    out_highTrackMuonDP =   (Float_t)(*highTrackMuonDP);
    out_highTrackMuonDPt =  (Float_t)(*highTrackMuonDPt);
    out_lowTrackMuonDR =    (Float_t)(*lowTrackMuonDR);
    out_lowTrackMuonDP =    (Float_t)(*lowTrackMuonDP);
    out_lowTrackMuonDPt =   (Float_t)(*lowTrackMuonDPt);
    out_thirdTrackMuonDR =  (Float_t)(*thirdTrackMuonDR);
    out_thirdTrackMuonDP =  (Float_t)(*thirdTrackMuonDP);
    out_thirdTrackMuonDPt =         (Float_t)(*thirdTrackMuonDPt);
    out_fourthTrackMuonDR =         (Float_t)(*fourthTrackMuonDR);
    out_fourthTrackMuonDP =         (Float_t)(*fourthTrackMuonDP);
    out_fourthTrackMuonDPt =        (Float_t)(*fourthTrackMuonDPt);
    out_tPFromPV =  (Float_t)(*tPFromPV);
    out_tMFromPV =  (Float_t)(*tMFromPV);
    out_tTFromPV =  (Float_t)(*tTFromPV);
    out_tFFromPV =  (Float_t)(*tFFromPV);
    out_tPFromPVCA =        (Float_t)(*tPFromPVCA);
    out_tMFromPVCA =        (Float_t)(*tMFromPVCA);
    out_tTFromPVCA =        (Float_t)(*tTFromPVCA);
    out_tFFromPVCA =        (Float_t)(*tFFromPVCA);
    out_tPFromPVDZ =        (Float_t)(*tPFromPVDZ);
    out_tMFromPVDZ =        (Float_t)(*tMFromPVDZ);
    out_tTFromPVDZ =        (Float_t)(*tTFromPVDZ);
    out_tFFromPVDZ =        (Float_t)(*tFFromPVDZ);
    out_five_m =    (Float_t)(*five_m);
    out_five_m_ref =        (Float_t)(*five_m_ref);
    out_five_mass_ppk =     (Float_t)(*five_mass_ppk);
    out_five_mass_kpp =     (Float_t)(*five_mass_kpp);
    out_five_mass_pkp =     (Float_t)(*five_mass_pkp);
    out_five_mass_ppp =     (Float_t)(*five_mass_ppp);
    out_fiveOne_pt =        (Float_t)(*fiveOne_pt);
    out_fiveOne_eta =       (Float_t)(*fiveOne_eta);
    out_fiveOne_phi =       (Float_t)(*fiveOne_phi);
    out_fiveOne_p =         (Float_t)(*fiveOne_p);
    out_fiveTwo_pt =        (Float_t)(*fiveTwo_pt);
    out_fiveTwo_eta =       (Float_t)(*fiveTwo_eta);
    out_fiveTwo_phi =       (Float_t)(*fiveTwo_phi);
    out_fiveTwo_p =         (Float_t)(*fiveTwo_p);
    out_fiveThree_pt =      (Float_t)(*fiveThree_pt);
    out_fiveThree_eta =     (Float_t)(*fiveThree_eta);
    out_fiveThree_phi =     (Float_t)(*fiveThree_phi);
    out_fiveThree_p =       (Float_t)(*fiveThree_p);
    out_fiveFour_pt =       (Float_t)(*fiveFour_pt);
    out_fiveFour_eta =      (Float_t)(*fiveFour_eta);
    out_fiveFour_phi =      (Float_t)(*fiveFour_phi);
    out_fiveFour_p =        (Float_t)(*fiveFour_p);
    out_fiveFive_pt =       (Float_t)(*fiveFive_pt);
    out_fiveFive_eta =      (Float_t)(*fiveFive_eta);
    out_fiveFive_phi =      (Float_t)(*fiveFive_phi);
    out_fiveFive_p =        (Float_t)(*fiveFive_p);
    out_five_cosAlpha =     (Float_t)(*five_cosAlpha);
    out_five_ctauPV =       (Float_t)(*five_ctauPV);
    out_five_ctauErrPV =    (Float_t)(*five_ctauErrPV);
    out_five_cosAlphaCA =   (Float_t)(*five_cosAlphaCA);
    out_five_ctauPVCA =     (Float_t)(*five_ctauPVCA);
    out_five_ctauErrPVCA =  (Float_t)(*five_ctauErrPVCA);
    out_five_cosAlphaDZ =   (Float_t)(*five_cosAlphaDZ);
    out_five_ctauPVDZ =     (Float_t)(*five_ctauPVDZ);
    out_five_ctauErrPVDZ =  (Float_t)(*five_ctauErrPVDZ);
    out_five_cosAlphaBS =   (Float_t)(*five_cosAlphaBS);
    out_five_ctauPVBS =     (Float_t)(*five_ctauPVBS);
    out_five_ctauErrPVBS =  (Float_t)(*five_ctauErrPVBS);
    out_five_vProb =        (Float_t)(*five_vProb);
    out_five_nDof =         (Float_t)(*five_nDof);
    out_five_vChi2 =        (Float_t)(*five_vChi2);
    out_five_vx =   (Float_t)(*five_vx);
    out_five_vy =   (Float_t)(*five_vy);
    out_five_vz =   (Float_t)(*five_vz);
    out_five_charge =       (Float_t)(*five_charge);
    out_bestPV_X =  (Float_t)(*bestPV_X);
    out_bestPV_Y =  (Float_t)(*bestPV_Y);
    out_bestPV_Z =  (Float_t)(*bestPV_Z);
    out_cosAlphaPV_X =      (Float_t)(*cosAlphaPV_X);
    out_cosAlphaPV_Y =      (Float_t)(*cosAlphaPV_Y);
    out_cosAlphaPV_Z =      (Float_t)(*cosAlphaPV_Z);
    out_bS_X =      (Float_t)(*bS_X);
    out_bS_Y =      (Float_t)(*bS_Y);
    out_bS_Z =      (Float_t)(*bS_Z);
    out_zPV_X =     (Float_t)(*zPV_X);
    out_zPV_Y =     (Float_t)(*zPV_Y);
    out_zPV_Z =     (Float_t)(*zPV_Z);
    out_lowMuon_isTight =   (Float_t)(*lowMuon_isTight);
    out_lowMuon_isLoose =   (Float_t)(*lowMuon_isLoose);
    out_lowMuon_isSoft =    (Float_t)(*lowMuon_isSoft);
    out_lowMuon_isMedium =  (Float_t)(*lowMuon_isMedium);
    out_lowMuon_isHighPt =  (Float_t)(*lowMuon_isHighPt);
    out_lowMuon_isTracker =         (Float_t)(*lowMuon_isTracker);
    out_lowMuon_isGlobal =  (Float_t)(*lowMuon_isGlobal);
    out_lowMuon_NPixelHits =        (Float_t)(*lowMuon_NPixelHits);
    out_lowMuon_NStripHits =        (Float_t)(*lowMuon_NStripHits);
    out_lowMuon_NTrackhits =        (Float_t)(*lowMuon_NTrackhits);
    out_lowMuon_NBPixHits =         (Float_t)(*lowMuon_NBPixHits);
    out_lowMuon_NPixLayers =        (Float_t)(*lowMuon_NPixLayers);
    out_lowMuon_NTraLayers =        (Float_t)(*lowMuon_NTraLayers);
    out_lowMuon_NStrLayers =        (Float_t)(*lowMuon_NStrLayers);
    out_lowMuon_NBPixLayers =       (Float_t)(*lowMuon_NBPixLayers);
    out_highMuon_isTight =  (Float_t)(*highMuon_isTight);
    out_highMuon_isLoose =  (Float_t)(*highMuon_isLoose);
    out_highMuon_isSoft =   (Float_t)(*highMuon_isSoft);
    out_highMuon_isMedium =         (Float_t)(*highMuon_isMedium);
    out_highMuon_isHighPt =         (Float_t)(*highMuon_isHighPt);
    out_highMuon_isTracker =        (Float_t)(*highMuon_isTracker);
    out_highMuon_isGlobal =         (Float_t)(*highMuon_isGlobal);
    out_highMuon_NPixelHits =       (Float_t)(*highMuon_NPixelHits);
    out_highMuon_NStripHits =       (Float_t)(*highMuon_NStripHits);
    out_highMuon_NTrackhits =       (Float_t)(*highMuon_NTrackhits);
    out_highMuon_NBPixHits =        (Float_t)(*highMuon_NBPixHits);
    out_highMuon_NPixLayers =       (Float_t)(*highMuon_NPixLayers);
    out_highMuon_NTraLayers =       (Float_t)(*highMuon_NTraLayers);
    out_highMuon_NStrLayers =       (Float_t)(*highMuon_NStrLayers);
    out_highMuon_NBPixLayers =      (Float_t)(*highMuon_NBPixLayers);
    out_lowMuon_type =      (Float_t)(*lowMuon_type);
    out_highMuon_type =     (Float_t)(*highMuon_type);
    out_highTrack_NPixelHits =      (Float_t)(*highTrack_NPixelHits);
    out_highTrack_NStripHits =      (Float_t)(*highTrack_NStripHits);
    out_highTrack_NTrackhits =      (Float_t)(*highTrack_NTrackhits);
    out_highTrack_NBPixHits =       (Float_t)(*highTrack_NBPixHits);
    out_highTrack_NPixLayers =      (Float_t)(*highTrack_NPixLayers);
    out_highTrack_NTraLayers =      (Float_t)(*highTrack_NTraLayers);
    out_highTrack_NStrLayers =      (Float_t)(*highTrack_NStrLayers);
    out_highTrack_NBPixLayers =     (Float_t)(*highTrack_NBPixLayers);
    out_lowTrack_NPixelHits =       (Float_t)(*lowTrack_NPixelHits);
    out_lowTrack_NStripHits =       (Float_t)(*lowTrack_NStripHits);
    out_lowTrack_NTrackhits =       (Float_t)(*lowTrack_NTrackhits);
    out_lowTrack_NBPixHits =        (Float_t)(*lowTrack_NBPixHits);
    out_lowTrack_NPixLayers =       (Float_t)(*lowTrack_NPixLayers);
    out_lowTrack_NTraLayers =       (Float_t)(*lowTrack_NTraLayers);
    out_lowTrack_NStrLayers =       (Float_t)(*lowTrack_NStrLayers);
    out_lowTrack_NBPixLayers =      (Float_t)(*lowTrack_NBPixLayers);
    out_thirdTrack_NPixelHits =     (Float_t)(*thirdTrack_NPixelHits);
    out_thirdTrack_NStripHits =     (Float_t)(*thirdTrack_NStripHits);
    out_thirdTrack_NTrackhits =     (Float_t)(*thirdTrack_NTrackhits);
    out_thirdTrack_NBPixHits =      (Float_t)(*thirdTrack_NBPixHits);
    out_thirdTrack_NPixLayers =     (Float_t)(*thirdTrack_NPixLayers);
    out_thirdTrack_NTraLayers =     (Float_t)(*thirdTrack_NTraLayers);
    out_thirdTrack_NStrLayers =     (Float_t)(*thirdTrack_NStrLayers);
    out_thirdTrack_NBPixLayers =    (Float_t)(*thirdTrack_NBPixLayers);
    out_fourthTrack_NPixelHits =    (Float_t)(*fourthTrack_NPixelHits);
    out_fourthTrack_NStripHits =    (Float_t)(*fourthTrack_NStripHits);
    out_fourthTrack_NTrackhits =    (Float_t)(*fourthTrack_NTrackhits);
    out_fourthTrack_NBPixHits =     (Float_t)(*fourthTrack_NBPixHits);
    out_fourthTrack_NPixLayers =    (Float_t)(*fourthTrack_NPixLayers);
    out_fourthTrack_NTraLayers =    (Float_t)(*fourthTrack_NTraLayers);
    out_fourthTrack_NStrLayers =    (Float_t)(*fourthTrack_NStrLayers);
    out_fourthTrack_NBPixLayers =   (Float_t)(*fourthTrack_NBPixLayers);
    out_six_m_kkpp =        (Float_t)(*six_m_kkpp);
    out_six_m_ref_kkpp =    (Float_t)(*six_m_ref_kkpp);
    out_six_mass_ppkk =     (Float_t)(*six_mass_ppkk);
    out_six_mass_pkpk =     (Float_t)(*six_mass_pkpk);
    out_six_mass_pppp =     (Float_t)(*six_mass_pppp);
    out_six_mass_kpkp =     (Float_t)(*six_mass_kpkp);
    out_six_mass_kppk =     (Float_t)(*six_mass_kppk);
    out_six_mass_kkkk =     (Float_t)(*six_mass_kkkk);
    out_six_pt =    (Float_t)(*six_pt);
    out_six_eta =   (Float_t)(*six_eta);
    out_six_phi =   (Float_t)(*six_phi);
    out_six_p =     (Float_t)(*six_p);
    out_six_cosAlpha =      (Float_t)(*six_cosAlpha);
    out_six_ctauPV =        (Float_t)(*six_ctauPV);
    out_six_ctauErrPV =     (Float_t)(*six_ctauErrPV);
    out_six_cosAlphaCA =    (Float_t)(*six_cosAlphaCA);
    out_six_ctauPVCA =      (Float_t)(*six_ctauPVCA);
    out_six_ctauErrPVCA =   (Float_t)(*six_ctauErrPVCA);
    out_six_cosAlphaDZ =    (Float_t)(*six_cosAlphaDZ);
    out_six_ctauPVDZ =      (Float_t)(*six_ctauPVDZ);
    out_six_ctauErrPVDZ =   (Float_t)(*six_ctauErrPVDZ);
    out_six_cosAlphaBS =    (Float_t)(*six_cosAlphaBS);
    out_six_ctauPVBS =      (Float_t)(*six_ctauPVBS);
    out_six_ctauErrPVBS =   (Float_t)(*six_ctauErrPVBS);
    out_six_vProb =         (Float_t)(*six_vProb);
    out_six_nDof =  (Float_t)(*six_nDof);
    out_six_vChi2 =         (Float_t)(*six_vChi2);
    out_six_vx =    (Float_t)(*six_vx);
    out_six_vy =    (Float_t)(*six_vy);
    out_six_vz =    (Float_t)(*six_vz);
    out_six_charge =        (Float_t)(*six_charge);


    outTree->Fill();

  }

  return kTRUE;
}

void SixTracks::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  TDirectory *savedir = gDirectory;
  if (fOut)
  {
    fOut->cd();
    gStyle->SetOptStat(111111) ;


    outTree->Write();
    OutFile->Print();
    fOutput->Add(OutFile);
    gDirectory = savedir;
    fOut->Close();

  }

}

void SixTracks::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

}
