#define mumukk_cxx
// The class definition in Y4140.h has been generated automatically
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
// Root > T->Process("Y4140.C")
// Root > T->Process("Y4140.C","some options")
// Root > T->Process("Y4140.C+")
//

#include "mumukk.h"
#include <TH2.h>
#include <TStyle.h>

/// SEMRA added
#include <TCanvas.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TSystem.h>
#include <TTree.h>
#include <TBranch.h>
//#include <TCint.h>
#include <TRandom.h>
#include <TMath.h>
#include <TDirectory.h>
#include "TEnv.h"
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TString.h>
#include <TProof.h>
#include <TProofOutputFile.h>
#include "TLorentzVector.h"
#include "TPoint.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <map>
#include "TSelectorCint.h"

#include <iostream>
#include <fstream>

void mumukk::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

}


void mumukk::SlaveBegin(TTree * /*tree*/)
{
  //std::cout << "SLAVE BEGIN" << std::endl;

  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  std::string outputString = "mumukk_tree.root";
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

  outTuple = new TNtuple("outuple","outuple","run:evt:lum:xM:kkM:mumuM:xL:xPt:xEta:xVtx:xCos:xHlt");


}

bool mumukk::Process(Long64_t entry)
{
  //std::cout << "PROCESS" << std::endl;

  GetEntry(entry); /// ROOT event loop

  // The entry argument specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either Y4140::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing

  double muon_mass = 0.1056583715;
  double kaonCh_mass = 0.493677;
  double Bs0_Low_Mass = 5.1;
  double Bs0_High_Mass = 5.6;
  double Y_High_Mass = 4.35;


  ////////////////// Bs0 & X(4140) Loop //////////////////
  bool HLT_4_v9 = false, HLT_4_v10 = false, HLT_4_v11 = false, HLT_4_v12 = false ;
  bool HLT_8_v3 = false, HLT_8_v4 = false, HLT_8_v5 = false, HLT_8_v6 = false, HLT_8_v7 = false ;
  bool HLT_4_vAny = false;
  bool HLT_8_vAny = false;
  bool HLT_Any = false;

  std::vector<bool> hltsFlags;

  for (Int_t i = 0; i != abs((int)(TrigRes->size())); ++i)
  {
    if ( TrigNames->at(i).find("HLT_DoubleMu4_Jpsi_Displaced_v9") != string::npos  &&  TrigRes->at(i) == 1) HLT_4_v9 = true;
    if ( TrigNames->at(i).find("HLT_DoubleMu4_Jpsi_Displaced_v10") != string::npos  &&  TrigRes->at(i) == 1) HLT_4_v10 = true;
    if ( TrigNames->at(i).find("HLT_DoubleMu4_Jpsi_Displaced_v11") != string::npos  &&  TrigRes->at(i) == 1) HLT_4_v11 = true;
    if ( TrigNames->at(i).find("HLT_DoubleMu4_Jpsi_Displaced_v12") != string::npos  &&  TrigRes->at(i) == 1) HLT_4_v12 = true;

    if ( TrigNames->at(i).find("HLT_Dimuon8_Jpsi_v3") != string::npos  &&  TrigRes->at(i) == 1) HLT_8_v3 = true;
    if ( TrigNames->at(i).find("HLT_Dimuon8_Jpsi_v4") != string::npos  &&  TrigRes->at(i) == 1) HLT_8_v4 = true;
    if ( TrigNames->at(i).find("HLT_Dimuon8_Jpsi_v5") != string::npos  &&  TrigRes->at(i) == 1) HLT_8_v5 = true;
    if ( TrigNames->at(i).find("HLT_Dimuon8_Jpsi_v6") != string::npos  &&  TrigRes->at(i) == 1) HLT_8_v6 = true;
    if ( TrigNames->at(i).find("HLT_Dimuon8_Jpsi_v7") != string::npos  &&  TrigRes->at(i) == 1) HLT_8_v7 = true;
  }

  if (HLT_4_v9 || HLT_4_v10 || HLT_4_v11 || HLT_4_v12) HLT_4_vAny = true;
  if (HLT_8_v3 || HLT_8_v4 || HLT_8_v5 || HLT_8_v6 || HLT_8_v7) HLT_8_vAny = true;
  if (HLT_8_vAny || HLT_4_vAny) HLT_Any = true;

  if(!HLT_Any) return kTRUE;

  int muonQual[] = {1,3,4,12};

  // for (std::map < std::pair <int,int>, int>::iterator lumiMapIt=lumiMap.begin(); lumiMapIt!=lumiMap.end(); ++lumiMapIt)
  //  runlumi << lumiMapIt->first.first << " : " << lumiMapIt->first.second << std::endl;


  // runLumi->Fill(runNum,lumiNum);
  std::vector<double> xCandsB0sMasses,B0_Cand_Masses_L_Sides,B0_Cand_Masses_R_Sides,B0_Cand_Masses,B0_phi_masses;
  B0_phi_masses.clear();
  xCandsB0sMasses.clear();
  B0_Cand_Masses_R_Sides.clear();
  B0_Cand_Masses_L_Sides.clear();
  B0_Cand_Masses.clear();

  Float_t          run_out;
  Float_t          evt_out;
  Float_t          lum_out;
  Float_t         X_mass;
  Float_t         kk_mass;
  Float_t         mumu_mass;
  Float_t         X_LFly;
  Float_t         X_pt;
  Float_t         X_eta;
  Float_t         X_vtx;
  Float_t         X_cosAlpha;
  Float_t         X_hlt;

  for(unsigned int iX=0; iX<nX; ++iX)
  {

    bool muonQualityCut = false, muonChiCut = false, muonPhitsCut = false, muonShitsCut = false;
    bool muonDZPVCut= false, muonDXYPVCut = false, muonSoftCuts = false, muonsCuts = false;

    bool jPsiPtCut = false,jPsiMassCut = false, jPsiVtxCut = false, jPsiMuEtaPtCut = false, jPsiMusPtCut = false, jPsiCuts = false;

    bool kaonOneChiCut = false, kaonOnePhitsCut = false, kaonOneShitsCut = false, kaonTwoChiCut = false;
    bool kaonTwoPhitsCut = false, kaonTwoShitsCut = false, kaonsPt = false;

    bool kaonOneCuts = false, kaonTwoCuts = false, kaonsCuts = false;
    bool cosAlphaCut = false, vtxCLCut = false;

    bool CWMass = false, SWMass = false;
    bool promptRegion = false, mixedRegion = false, nonPromptRegion = false;

    bool extraCuts = false;

    int iJPsi = (*XMuMuIdx)[iX];

    //doneJPsiIt = doneJPsi.find(iJPsi);

    //if(doneJPsiIt!=doneJPsi.end())
    //continue;
    //else
    //doneJPsi[iJPsi] = 1.0;

    int iMu1 = (*mu1Idx)[iJPsi] ; // define for original muon1
    int iMu2 = (*mu2Idx)[iJPsi] ; // define for original muon2
    int iK1 = (*ka1Idx)[iX] ; // define for original kaon1
    int iK2 = (*ka2Idx)[iX] ;

    double mu1_E = 0., mu2_E = 0., K1_E = 0., K2_E = 0.;

    TLorentzVector mu1, mu2, refMu1, refMu2;

    mu1_E = sqrt( pow((*muPx)[iMu1], 2) + pow((*muPy)[iMu1], 2) + pow((*muPz)[iMu1], 2) + pow(muon_mass, 2) ) ;
    mu2_E = sqrt( pow((*muPx)[iMu2], 2) + pow((*muPy)[iMu2], 2) + pow((*muPz)[iMu2], 2) + pow(muon_mass, 2) ) ;
    refMu1.SetPxPyPzE((*Muon1Px_MuMuKK)[iX],(*Muon1Py_MuMuKK)[iX],(*Muon1Pz_MuMuKK)[iX],(*Muon1E_MuMuKK)[iX]);
    refMu2.SetPxPyPzE((*Muon2Px_MuMuKK)[iX],(*Muon2Py_MuMuKK)[iX],(*Muon2Pz_MuMuKK)[iX],(*Muon2E_MuMuKK)[iX]);
    mu1.SetPxPyPzE( (*muPx)[iMu1], (*muPy)[iMu1], (*muPz)[iMu1], mu1_E) ;
    mu2.SetPxPyPzE( (*muPx)[iMu2], (*muPy)[iMu2], (*muPz)[iMu2], mu2_E) ;

    TLorentzVector JPsi;
    JPsi = refMu1 + refMu2;

    TLorentzVector JPsiOriginal;
    JPsiOriginal = mu1 + mu2;

    TLorentzVector kaon1,kaon2;

    K1_E=sqrt(pow((*Kaon1Px_MuMuKK)[iX],2)+pow((*Kaon1Py_MuMuKK)[iX],2)+pow((*Kaon1Pz_MuMuKK)[iX],2)+pow(kaonCh_mass,2));
    kaon1.SetPxPyPzE((*Kaon1Px_MuMuKK)[iX],(*Kaon1Py_MuMuKK)[iX],(*Kaon1Pz_MuMuKK)[iX],K1_E);
    K2_E=sqrt(pow((*Kaon2Px_MuMuKK)[iX],2)+pow((*Kaon2Py_MuMuKK)[iX],2)+pow((*Kaon2Pz_MuMuKK)[iX],2)+pow(kaonCh_mass,2));
    kaon2.SetPxPyPzE((*Kaon2Px_MuMuKK)[iX],(*Kaon2Py_MuMuKK)[iX],(*Kaon2Pz_MuMuKK)[iX],K2_E);

    TLorentzVector Phi;
    Phi = kaon1 + kaon2;

    // Muon1_Mass->Fill(mu1.M());
    // Muon2_Mass->Fill(mu2.M());

    TLorentzVector XCand;
    XCand = JPsi + Phi;

    SWMass = (((XCand.M() > 4.05) && (XCand.M() < 4.8)));
    CWMass = ((XCand.M() > 5.2) && (XCand.M() < 5.55));

    mixedRegion     = (((*XLxyPV)[iX] / (*XLxyPVE)[iX]) >= 2. && ((*XLxyPV)[iX] / (*XLxyPVE)[iX])  <= 3.5);
    nonPromptRegion = (((*XLxyPV)[iX] / (*XLxyPVE)[iX]) > 3.5);
    promptRegion    = (((*XLxyPV)[iX] / (*XLxyPVE)[iX]) < 2.);

    // muonQualityCut = ( ((*muQual)[iMu1]) & (1 << muonQual[3]) ) && ( ((*muQual)[iMu2]) & (1 << muonQual[3]) );
    // muonChiCut     = (( ( (*muChi2)[iMu1] / (*muNDF)[iMu1] ) < 3 ) && ( ( (*muChi2)[iMu2] / (*muNDF)[iMu2] ) < 3 ));
    // muonPhitsCut   = ((*muPhits)[iMu1] > 0 && (*muPhits)[iMu2] > 0);
    // muonShitsCut   = ((*muShits)[iMu1] > 5 && (*muShits)[iMu2] > 5);
    muonDZPVCut    = (fabs((*muDzVtx)[iMu1]) < 20.0 && fabs((*muDzVtx)[iMu2]) < 20.0);
    muonDXYPVCut   = (fabs((*muDxyVtx)[iMu1]) < 0.3 && fabs((*muDxyVtx)[iMu2]) < 0.3);

    //muonsCuts = muonQualityCut && muonChiCut && muonShitsCut && muonPhitsCut && muonDZPVCut && muonDXYPVCut;
    muonsCuts = muonDZPVCut && muonDXYPVCut;

    if(!muonsCuts) break;

    jPsiPtCut      = (JPsi.Pt() > 7.0);
    jPsiMassCut    =  (JPsi.M()<3.15 && JPsi.M()>3.04);
    jPsiVtxCut     = ((*MuMuVtx_CL)[iJPsi]) > 0.1;
    jPsiMuEtaPtCut = (fabs(mu1.Eta()) < 2.2) && (fabs(mu2.Eta()) < 2.2);
    jPsiMusPtCut   = ((mu1.Pt() > 4.0) && (mu2.Pt() > 4.0));

    jPsiCuts = jPsiPtCut && jPsiMassCut && jPsiVtxCut  && jPsiMuEtaPtCut && jPsiMuEtaPtCut && jPsiMusPtCut;

    if(!jPsiCuts) break;

    // kaonOneChiCut    = (((*trackChi2)[iK1] / (*trackNDF)[iK1]) < 5.0);
    // kaonOnePhitsCut  = ((*trackPhits)[iK1] > 0);
    // kaonOneShitsCut  = ((*trackShits)[iK1] >= 7);
    // kaonTwoChiCut    = (((*trackChi2)[iK2] / (*trackNDF)[iK2]) < 5.0);
    // kaonTwoPhitsCut  = ((*trackPhits)[iK2] > 0);
    // kaonTwoShitsCut  = ((*trackShits)[iK2] >= 7);
    kaonsPt = ((kaon1.Pt()>0.7) && (kaon2.Pt()>0.7));

    kaonOneCuts = kaonOneChiCut && kaonOnePhitsCut && kaonOneShitsCut;
    kaonTwoCuts = kaonTwoChiCut && kaonTwoPhitsCut && kaonTwoShitsCut;
    kaonsCuts = kaonOneCuts && kaonTwoCuts && kaonsPt;
    kaonsCuts = kaonsPt;

    if(!kaonsCuts) break;

    cosAlphaCut = (fabs((*XCosAlphaPV)[iX]) > 0.99);
    vtxCLCut =  (((*XVtx_CL)[iX]) > 0.01);

    extraCuts = vtxCLCut && cosAlphaCut;

    mixedRegion     = (((*XLxyPV)[iX] / (*XLxyPVE)[iX]) >= 2. && ((*XLxyPV)[iX] / (*XLxyPVE)[iX])  <= 3.0);
    nonPromptRegion = (((*XLxyPV)[iX] / (*XLxyPVE)[iX]) > 3.0);
    promptRegion    = (((*XLxyPV)[iX] / (*XLxyPVE)[iX]) < 2.);

    if(!extraCuts) break;

    if(muonsCuts && kaonsCuts && jPsiCuts && extraCuts && HLT_Any)
    //if(HLT_Any)
    {
      X_hlt = 0.0;
      if (HLT_4_vAny) X_hlt += 4.0;
      if (HLT_8_vAny) X_hlt += 8.0;

      X_mass  = XCand.M();
      kk_mass = Phi.M();
      mumu_mass = JPsi.M();

      X_LFly = ((*XLxyPV)[iX] / (*XLxyPVE)[iX]);
      X_pt = XCand.Pt();
      X_eta = XCand.Eta();
      X_vtx = ((*XVtx_CL)[iX]);
      X_cosAlpha = fabs((*XCosAlphaPV)[iX]);

      run_out = runNum;
      evt_out = evtNum;
      lum_out = lumiNum;

      outTuple->Fill(run_out,evt_out,lum_out,X_mass,kk_mass,mumu_mass,X_LFly,X_pt,X_eta,X_vtx,X_cosAlpha,X_hlt);
    }


  }

  return kTRUE;
}


void mumukk::SlaveTerminate()
{

  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  TDirectory *savedir = gDirectory;
  if (fOut)
  {
    fOut->cd();
    gStyle->SetOptStat(111111) ;


    outTuple->Write();
    OutFile->Print();
    fOutput->Add(OutFile);
    gDirectory = savedir;
    fOut->Close();

  }


}


void mumukk::Terminate()
{


}
