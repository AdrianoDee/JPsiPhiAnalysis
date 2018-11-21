#define TwoMuTwoK_2012_cxx
// The class definition in TwoMuTwoK_2012.h has been generated automatically
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
// root> T->Process("TwoMuTwoK_2012.C")
// root> T->Process("TwoMuTwoK_2012.C","some options")
// root> T->Process("TwoMuTwoK_2012.C+")
//


#include "TwoMuTwoK_2012.h"
#include <TH2.h>
#include <TStyle.h>

void TwoMuTwoK_2012::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
}

void TwoMuTwoK_2012::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).
  TString option = GetOption();

  std::string outputString = "TwoMuTwoK_five_tree.root";
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

  outTree = new TTree("2mu2kSkimmedTree","2mu2kSkimmedTree");

  outTree->Branch("event",      &event,   "event/F");
  outTree->Branch("run",      &run,   "run/F");
  outTree->Branch("lumi",      &lumi,   "lumi/F");

  outTree->Branch("HLT_Dimuon8_Jpsi", &hlt8, "HLT_Dimuon8_Jpsi/F");
  outTree->Branch("HLT_DoubleMu4_Jpsi_Displaced", &hlt4, "HLT_DoubleMu4_Jpsi_Displaced/F");

  outTree->Branch("priVtx_n",     &out_priVtx_n,  "priVtx_n/F");
  outTree->Branch("priVtx_X",     &out_priVtx_X,  "priVtx_X/F");
  outTree->Branch("priVtx_Y",     &out_priVtx_Y,  "priVtx_Y/F");
  outTree->Branch("priVtx_Z",     &out_priVtx_Z,  "priVtx_Z/F");

  outTree->Branch("muOnePx",         &out_muOnePx,      "muOnePx/F");
  outTree->Branch("muOnePy",         &out_muOnePy,      "muOnePy/F");
  outTree->Branch("muOnePz",         &out_muOnePz,      "muOnePz/F");
  outTree->Branch("muOneCharge",     &out_muOneCharge,  "muOneCharge/F");
  outTree->Branch("muOneD0",         &out_muOneD0,      "muOneD0/F");
  outTree->Branch("muOneDz",         &out_muOneDz,      "muOneDz/F");
  outTree->Branch("muOneChi2",       &out_muOneChi2,    "muOneChi2/F");
  outTree->Branch("muOneNDF",        &out_muOneNDF,     "muOneNDF/F");
  outTree->Branch("muOnePhits",      &out_muOnePhits,   "muOnePhits/F");
  outTree->Branch("muOneShits",      &out_muOneShits,   "muOneShits/F");
  outTree->Branch("muOneLayersTr",   &out_muOneLayersTr,        "muOneLayersTr/F");
  outTree->Branch("muOneLayersPix",  &out_muOneLayersPix,       "muOneLayersPix/F");
  outTree->Branch("muOneD0E",        &out_muOneD0E,     "muOneD0E/F");
  outTree->Branch("muOneDzVtxErr",   &out_muOneDzVtxErr,        "muOneDzVtxErr/F");
  outTree->Branch("muOneIsGlobal",   &out_muOneIsGlobal,        "muOneIsGlobal/F");
  outTree->Branch("muOneIsPF",       &out_muOneIsPF,    "muOneIsPF/F");
  outTree->Branch("muOneIsSoft",     &out_muOneIsSoft,    "muOneIsSoft/F");

  outTree->Branch("muTwoPx",         &out_muTwoPx,      "muTwoPx/F");
  outTree->Branch("muTwoPy",         &out_muTwoPy,      "muTwoPy/F");
  outTree->Branch("muTwoPz",         &out_muTwoPz,      "muTwoPz/F");
  outTree->Branch("muTwoCharge",     &out_muTwoCharge,  "muTwoCharge/F");
  outTree->Branch("muTwoD0",         &out_muTwoD0,      "muTwoD0/F");
  outTree->Branch("muTwoDz",         &out_muTwoDz,      "muTwoDz/F");
  outTree->Branch("muTwoChi2",       &out_muTwoChi2,    "muTwoChi2/F");
  outTree->Branch("muTwoNDF",        &out_muTwoNDF,     "muTwoNDF/F");
  outTree->Branch("muTwoPhits",      &out_muTwoPhits,   "muTwoPhits/F");
  outTree->Branch("muTwoShits",      &out_muTwoShits,   "muTwoShits/F");
  outTree->Branch("muTwoLayersTr",   &out_muTwoLayersTr,        "muTwoLayersTr/F");
  outTree->Branch("muTwoLayersPix",  &out_muTwoLayersPix,       "muTwoLayersPix/F");
  outTree->Branch("muTwoD0E",        &out_muTwoD0E,     "muTwoD0E/F");
  outTree->Branch("muTwoDzVtxErr",   &out_muTwoDzVtxErr,        "muTwoDzVtxErr/F");
  outTree->Branch("muTwoIsGlobal",   &out_muTwoIsGlobal,        "muTwoIsGlobal/F");
  outTree->Branch("muTwoIsPF",       &out_muTwoIsPF,    "muTwoIsPF/F");
  outTree->Branch("muTwoIsSoft",     &out_muTwoIsSoft,    "muTwoIsSoft/F");

  outTree->Branch("kaonOnePx",          &out_kaonOnePx,           "kaonOnePx/F");
  outTree->Branch("kaonOnePy",          &out_kaonOnePy,           "kaonOnePy/F");
  outTree->Branch("kaonOnePz",          &out_kaonOnePz,           "kaonOnePz/F");
  outTree->Branch("kaonOneE",           &out_kaonOneE,            "kaonOneE/F");
  outTree->Branch("kaonOneChi2",      &out_kaonOneChi2,           "kaonOneChi2/F");
  outTree->Branch("kaonOneNDF",         &out_kaonOneNDF,          "kaonOneNDF/F");
  outTree->Branch("kaonOneD0",          &out_kaonOneD0,           "kaonOneChiOne/F");
  outTree->Branch("kaonOneDzVtx",       &out_kaonOneDzVtx,        "kaonOneDzVtx/F");
  outTree->Branch("kaonOneDxyVtx",      &out_kaonOneDxyVtx,       "kaonOneDxyVtx/F");
  outTree->Branch("kaonOnetrackPhits",  &out_kaonOnetrackPhits,   "kaonOnetrackPhits/F");
  outTree->Branch("kaonOnetrackShits",  &out_kaonOnetrackShits,   "kaonOnetrackShits/F");
  outTree->Branch("kaonOnetrackCharge", &out_kaonOnetrackCharge,  "kaonOnetrackCharge/F");
  outTree->Branch("kaonOneHighPurity",  &out_kaonOneHighPurity,   "kaonOneHighPurity/F");


  outTree->Branch("kaonTwoPx",          &out_kaonTwoPx,          "kaonTwoPx/F");
  outTree->Branch("kaonTwoPy",          &out_kaonTwoPy,          "kaonTwoPy/F");
  outTree->Branch("kaonTwoPz",          &out_kaonTwoPz,          "kaonTwoPz/F");
  outTree->Branch("kaonTwoE",           &out_kaonTwoE,           "kaonTwoE/F");
  outTree->Branch("kaonTwoChi2",      &out_kaonTwoChi2,          "kaonTwoChi2/F");
  outTree->Branch("kaonTwoNDF",         &out_kaonTwoNDF,         "kaonTwoNDF/F");
  outTree->Branch("kaonTwoD0",          &out_kaonTwoD0,          "kaonTwoChiTwo/F");
  outTree->Branch("kaonTwoDzVtx",       &out_kaonTwoDzVtx,       "kaonTwoDzVtx/F");
  outTree->Branch("kaonTwoDxyVtx",      &out_kaonTwoDxyVtx,      "kaonTwoDxyVtx/F");
  outTree->Branch("kaonTwotrackPhits",  &out_kaonTwotrackPhits,  "kaonTwotrackPhits/F");
  outTree->Branch("kaonTwotrackShits",  &out_kaonTwotrackShits,  "kaonTwotrackShits/F");
  outTree->Branch("kaonTwotrackCharge", &out_kaonTwotrackCharge,        "kaonTwotrackCharge/F");
  outTree->Branch("kaonTwoHighPurity",  &out_kaonTwoHighPurity,  "kaonTwoHighPurity/F");

  outTree->Branch("MuMuMass",     &out_MuMuMass,  "MuMuMass/F");
  outTree->Branch("MuMuMass_original",     &out_MuMuMas_original,  "MuMuMass/F");
  outTree->Branch("MuMuPx",       &out_MuMuPx,    "MuMuPx/F");
  outTree->Branch("MuMuPy",       &out_MuMuPy,    "MuMuPy/F");
  outTree->Branch("MuMuPz",       &out_MuMuPz,    "MuMuPz/F");
  outTree->Branch("MuMuVtx_CL",   &out_MuMuVtx_CL,        "MuMuVtx_CL/F");
  outTree->Branch("MuMuVtx_Chi2",         &out_MuMuVtx_Chi2,      "MuMuVtx_Chi2/F");
  outTree->Branch("MuMuDecayVtx_X",       &out_MuMuDecayVtx_X,    "MuMuDecayVtx_X/F");
  outTree->Branch("MuMuDecayVtx_Y",       &out_MuMuDecayVtx_Y,    "MuMuDecayVtx_Y/F");
  outTree->Branch("MuMuDecayVtx_Z",       &out_MuMuDecayVtx_Z,    "MuMuDecayVtx_Z/F");

  outTree->Branch("KKMass",       &out_KKMass,    "KKMass/F");
  outTree->Branch("KKPx",         &out_KKPx,      "KKPx/F");
  outTree->Branch("KKPy",         &out_KKPy,      "KKPy/F");
  outTree->Branch("KKPz",         &out_KKPz,      "KKPz/F");
  outTree->Branch("KKPt",         &out_KKPt,      "KKPz/F");

  outTree->Branch("nX",   &out_nX,        "nX/F");

  outTree->Branch("XMass_original",        &out_XMass_original,     "XMass/F");

  outTree->Branch("XMass",        &out_XMass,     "XMass/F");
  outTree->Branch("XPx",  &out_XPx,       "XPx/F");
  outTree->Branch("XPy",  &out_XPy,       "XPy/F");
  outTree->Branch("XPz",  &out_XPz,       "XPz/F");

  outTree->Branch("XVtx_CL",      &out_XVtx_CL,   "XVtx_CL/F");
  outTree->Branch("XVtx_Chi2",    &out_XVtx_Chi2,         "XVtx_Chi2/F");
  outTree->Branch("XDecayVtx_X",  &out_XDecayVtx_X,       "XDecayVtx_X/F");
  outTree->Branch("XDecayVtx_Y",  &out_XDecayVtx_Y,       "XDecayVtx_Y/F");
  outTree->Branch("XDecayVtx_Z",  &out_XDecayVtx_Z,       "XDecayVtx_Z/F");
  outTree->Branch("XDecayVtx_XE",         &out_XDecayVtx_XE,      "XDecayVtx_XE/F");
  outTree->Branch("XDecayVtx_YE",         &out_XDecayVtx_YE,      "XDecayVtx_YE/F");
  outTree->Branch("XDecayVtx_ZE",         &out_XDecayVtx_ZE,      "XDecayVtx_ZE/F");
  outTree->Branch("XCosAlphaBS",  &out_XCosAlphaBS,       "XCosAlphaBS/F");
  outTree->Branch("XCosAlpha3DBS",        &out_XCosAlpha3DBS,     "XCosAlpha3DBS/F");
  outTree->Branch("XCTauBS",      &out_XCTauBS,   "XCTauBS/F");
  outTree->Branch("XCTauBSE",     &out_XCTauBSE,  "XCTauBSE/F");
  outTree->Branch("XLxyBS",       &out_XLxyBS,    "XLxyBS/F");
  outTree->Branch("XLxyBSE",      &out_XLxyBSE,   "XLxyBSE/F");
  outTree->Branch("XLxyzBS",      &out_XLxyzBS,   "XLxyzBS/F");
  outTree->Branch("XLxyzBSE",     &out_XLxyzBSE,  "XLxyzBSE/F");
  outTree->Branch("XCosAlphaPV",  &out_XCosAlphaPV,       "XCosAlphaPV/F");
  outTree->Branch("XCosAlpha3DPV",        &out_XCosAlpha3DPV,     "XCosAlpha3DPV/F");
  outTree->Branch("XCTauPV",      &out_XCTauPV,   "XCTauPV/F");
  outTree->Branch("XCTauPVE",     &out_XCTauPVE,  "XCTauPVE/F");
  outTree->Branch("XLxyPV",       &out_XLxyPV,    "XLxyPV/F");
  outTree->Branch("XLxyPVE",      &out_XLxyPVE,   "XLxyPVE/F");
  outTree->Branch("XLxyzPV",      &out_XLxyzPV,   "XLxyzPV/F");
  outTree->Branch("XLxyzPVE",     &out_XLxyzPVE,  "XLxyzPVE/F");


}

Bool_t TwoMuTwoK_2012::Process(Long64_t entry)
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

  float muon_mass = 0.1056583745;
  float kaon_mass = 0.493677;
  ////////////////// Bs0 & X(4140) Loop //////////////////
  bool HLT_4_v9 = false, HLT_4_v10 = false, HLT_4_v11 = false, HLT_4_v12 = false ;
  bool HLT_8_v3 = false, HLT_8_v4 = false, HLT_8_v5 = false, HLT_8_v6 = false, HLT_8_v7 = false ;
  bool HLT_4_vAny = false;
  bool HLT_8_vAny = false;
  bool HLT_Any = false;

  std::vector<bool> hltsFlags;

  for (Int_t i = 0; i != abs((int)(TrigRes.GetSize())); ++i)
  {
    if ( TrigNames[i].find("HLT_DoubleMu4_Jpsi_Displaced_v9") != string::npos  &&  TrigRes[i] == 1) HLT_4_v9 = true;
    if ( TrigNames[i].find("HLT_DoubleMu4_Jpsi_Displaced_v10") != string::npos  &&  TrigRes[i] == 1) HLT_4_v10 = true;
    if ( TrigNames[i].find("HLT_DoubleMu4_Jpsi_Displaced_v11") != string::npos  &&  TrigRes[i] == 1) HLT_4_v11 = true;
    if ( TrigNames[i].find("HLT_DoubleMu4_Jpsi_Displaced_v12") != string::npos  &&  TrigRes[i] == 1) HLT_4_v12 = true;

    if ( TrigNames[i].find("HLT_Dimuon8_Jpsi_v3") != string::npos  &&  TrigRes[i] == 1) HLT_8_v3 = true;
    if ( TrigNames[i].find("HLT_Dimuon8_Jpsi_v4") != string::npos  &&  TrigRes[i] == 1) HLT_8_v4 = true;
    if ( TrigNames[i].find("HLT_Dimuon8_Jpsi_v5") != string::npos  &&  TrigRes[i] == 1) HLT_8_v5 = true;
    if ( TrigNames[i].find("HLT_Dimuon8_Jpsi_v6") != string::npos  &&  TrigRes[i] == 1) HLT_8_v6 = true;
    if ( TrigNames[i].find("HLT_Dimuon8_Jpsi_v7") != string::npos  &&  TrigRes[i] == 1) HLT_8_v7 = true;
  }

  if (HLT_4_v9 || HLT_4_v10 || HLT_4_v11 || HLT_4_v12) HLT_4_vAny = true;
  if (HLT_8_v3 || HLT_8_v4 || HLT_8_v5 || HLT_8_v6 || HLT_8_v7) HLT_8_vAny = true;
  if (HLT_8_vAny || HLT_4_vAny) HLT_Any = true;

  int muonQual[4] = {1,3,4,12};


  if(HLT_Any)
  for(int iX=0; iX<*nX; ++iX)
  {

    out_nX = *nX;

    hlt8 = (Float_t) HLT_8_vAny;
    hlt4 = (Float_t) HLT_4_vAny;


    int iJPsi = XMuMuIdx[iX];

    event = (Float_t) *evtNum;
    run = (Float_t) *runNum;
    lumi = (Float_t) *lumiNum;

    priVtx_n = (Float_t)(*out_priVtx_n);
    priVtx_X (Float_t)(*out_priVtx_X);
    priVtx_Y = (Float_t)(*out_priVtx_Y);
    priVtx_Z = (Float_t)(*out_priVtx_Z);

    int iMu1 = (mu1Idx)[iJPsi] ; // define for original muon1
    int iMu2 = (mu2Idx)[iJPsi] ; // define for original muon2
    int iK1 = (ka1Idx)[iX] ; // define for original kaon1
    int iK2 = (ka2Idx)[iX] ;
    //
    double mu1_E = 0., mu2_E = 0., K1_E = 0., K2_E = 0.;
    //
    TLorentzVector mu1, mu2, oMu1, oMu2;
    //
    mu1.SetPxPyPzE((Muon1Px_MuMuKK)[iX],(Muon1Py_MuMuKK)[iX],(Muon1Pz_MuMuKK)[iX],(Muon1E_MuMuKK)[iX]);
    mu2.SetPxPyPzE((Muon2Px_MuMuKK)[iX],(Muon2Py_MuMuKK)[iX],(Muon2Pz_MuMuKK)[iX],(Muon2E_MuMuKK)[iX]);
    //
    mu1_E = sqrt( pow((muPx)[iMu1], 2) + pow((muPy)[iMu1], 2) + pow((muPz)[iMu1], 2) + pow(muon_mass, 2) ) ;
    mu2_E = sqrt( pow((muPx)[iMu2], 2) + pow((muPy)[iMu2], 2) + pow((muPz)[iMu2], 2) + pow(muon_mass, 2) ) ;
    oMu1.SetPxPyPzE( (muPx)[iMu1], (muPy)[iMu1], (muPz)[iMu1], mu1_E) ;
    oMu2.SetPxPyPzE( (muPx)[iMu2], (muPy)[iMu2], (muPz)[iMu2], mu2_E) ;
    //
    TLorentzVector JPsi;
    JPsi = mu1 + mu2;
    //
    TLorentzVector JPsiOriginal;
    JPsiOriginal = oMu1 + oMu2;
    //
    TLorentzVector kaon1,kaon2;
    //
    K1_E=sqrt(pow((Kaon1Px_MuMuKK)[iX],2)+pow((Kaon1Py_MuMuKK)[iX],2)+pow((Kaon1Pz_MuMuKK)[iX],2)+pow(kaon_mass,2));
    kaon1.SetPxPyPzE((Kaon1Px_MuMuKK)[iX],(Kaon1Py_MuMuKK)[iX],(Kaon1Pz_MuMuKK)[iX],K1_E);
    K2_E=sqrt(pow((Kaon2Px_MuMuKK)[iX],2)+pow((Kaon2Py_MuMuKK)[iX],2)+pow((Kaon2Pz_MuMuKK)[iX],2)+pow(kaon_mass,2));
    kaon2.SetPxPyPzE((Kaon2Px_MuMuKK)[iX],(Kaon2Py_MuMuKK)[iX],(Kaon2Pz_MuMuKK)[iX],K2_E);
    //
    TLorentzVector Phi;
    Phi = kaon1 + kaon2;

    //
    // // Muon1_Mass->Fill(mu1.M());
    // // Muon2_Mass->Fill(mu2.M());
    //i

    TLorentzVector XCand;
    XCand = JPsi + Phi;

    out_XMass_original = XCand.M();
    // std::cout << out_XMass_original << std::endl;
    //
    // SWMass = (((XCand.M() > 4.0) && (XCand.M() < 5.0)));
    // CWMass = ((XCand.M() > 5.15) && (XCand.M() < 5.55));
    //
    //
    // mixedRegion     = (((*XLxyPV)[iX] / (*XLxyPVE)[iX]) >= 2. && ((*XLxyPV)[iX] / (*XLxyPVE)[iX])  <= 3.0);
    // nonPromptRegion = (((*XLxyPV)[iX] / (*XLxyPVE)[iX]) > 3.0);
    // promptRegion    = (((*XLxyPV)[iX] / (*XLxyPVE)[iX]) < 2.);
    //
    bool muonOneIsSoft = true;
    muonOneIsSoft = muonOneIsSoft && ( ((muQual)[iMu1]) & (1 << muonQual[3]) );
    muonOneIsSoft = muonOneIsSoft && ( ( (muChi2)[iMu1] / (muNDF)[iMu1] ) < 3 );
    muonOneIsSoft = muonOneIsSoft && ((muPhits)[iMu1] > 0);
    muonOneIsSoft = muonOneIsSoft && ((muShits)[iMu1] > 5);
    muonOneIsSoft = muonOneIsSoft && (fabs((muDzVtx)[iMu1]) < 20.0);
    muonOneIsSoft = muonOneIsSoft && (fabs((muDxyVtx)[iMu1]) < 0.3);

    bool muonTwoIsSoft = true;
    muonTwoIsSoft = muonTwoIsSoft && ( ((muQual)[iMu2]) & (1 << muonQual[3]) );
    muonTwoIsSoft = muonTwoIsSoft && ( ( (muChi2)[iMu2] / (muNDF)[iMu2] ) < 3 );
    muonTwoIsSoft = muonTwoIsSoft && ((muPhits)[iMu2] > 0);
    muonTwoIsSoft = muonTwoIsSoft && ((muShits)[iMu2] > 5);
    muonTwoIsSoft = muonTwoIsSoft && (fabs((muDzVtx)[iMu2]) < 20.0);
    muonTwoIsSoft = muonTwoIsSoft && (fabs((muDxyVtx)[iMu2]) < 0.3);

    out_MuMuMas_original = (Float_t)  JPsi.M();

    out_muOnePx =      (Float_t)(muPx[iMu1]);
    out_muOnePy =      (Float_t)(muPy[iMu1]);
    out_muOnePz =      (Float_t)(muPz[iMu1]);
    out_muOneCharge =  (Float_t)(muCharge[iMu1]);
    out_muOneD0 =      (Float_t)(muD0[iMu1]);
    out_muOneDz =      (Float_t)(muDz[iMu1]);
    out_muOneChi2 =    (Float_t)(muChi2[iMu1]);
    out_muOneNDF =     (Float_t)(muNDF[iMu1]);
    out_muOnePhits =   (Float_t)(muPhits[iMu1]);
    out_muOneShits =   (Float_t)(muShits[iMu1]);
    out_muOneLayersTr =        (Float_t)(muLayersTr[iMu1]);
    out_muOneLayersPix =       (Float_t)(muLayersPix[iMu1]);
    out_muOneD0E =     (Float_t)(muD0E[iMu1]);
    out_muOneDzVtxErr =        (Float_t)(muDzVtxErr[iMu1]);
    out_muOneIsGlobal =        (Float_t)((*muIsGlobal)[iMu1]);
    out_muOneIsPF =    (Float_t)((*muIsPF)[iMu1]);
    out_muOneIsSoft  = muonOneIsSoft;

    out_kaonOnePx  = ( trackPx[iK1]);
    out_kaonOnePy  = ( trackPy[iK1]);
    out_kaonOnePz  = ( trackPz[iK1]);
    out_kaonOneE  = ( trackEnergy[iK1]);
    out_kaonOneChi2  = ( trackChi2[iK1]);
    out_kaonOneNDF  = ( trackNDF[iK1]);
    out_kaonOneD0  = ( trackD0[iK1]);
    out_kaonOneDzVtx  = ( trackDzVtx[iK1]);
    out_kaonOneDxyVtx  = ( trackDxyVtx[iK1]);
    out_kaonOnetrackPhits  = ( trackPhits[iK1]);
    out_kaonOnetrackShits  = ( trackShits[iK1]);
    out_kaonOneHighPurity = (TrackHighPurity[iK1]);
    out_kaonOnetrackCharge  = ( trackCharge[iK1]);

    out_kaonTwoPx  = ( trackPx[iK2]);
    out_kaonTwoPy  = ( trackPy[iK2]);
    out_kaonTwoPz  = ( trackPz[iK2]);
    out_kaonTwoE  = ( trackEnergy[iK2]);
    out_kaonTwoChi2  = ( trackChi2[iK2]);
    out_kaonTwoNDF  = ( trackNDF[iK2]);
    out_kaonTwoD0  = ( trackD0[iK2]);
    out_kaonTwoDzVtx  = ( trackDzVtx[iK2]);
    out_kaonTwoDxyVtx  = ( trackDxyVtx[iK2]);
    out_kaonTwotrackPhits  = ( trackPhits[iK2]);
    out_kaonTwotrackShits  = ( trackShits[iK2]);
    out_kaonTwoHighPurity = (TrackHighPurity[iK2]);
    out_kaonTwotrackCharge  = (Float_t)  ( trackCharge[iK2]);

    out_muTwoPx =      (Float_t)(muPx[iMu2]);
    out_muTwoPy =      (Float_t)(muPy[iMu2]);
    out_muTwoPz =      (Float_t)(muPz[iMu2]);
    out_muTwoCharge =  (Float_t)(muCharge[iMu2]);
    out_muTwoD0 =      (Float_t)(muD0[iMu2]);
    out_muTwoDz =      (Float_t)(muDz[iMu2]);
    out_muTwoChi2 =    (Float_t)(muChi2[iMu2]);
    out_muTwoNDF =     (Float_t)(muNDF[iMu2]);
    out_muTwoPhits =   (Float_t)(muPhits[iMu2]);
    out_muTwoShits =   (Float_t)(muShits[iMu2]);
    out_muTwoLayersTr =        (Float_t)(muLayersTr[iMu2]);
    out_muTwoLayersPix =       (Float_t)(muLayersPix[iMu2]);
    out_muTwoD0E =     (Float_t)(muD0E[iMu2]);
    out_muTwoDzVtxErr =        (Float_t)(muDzVtxErr[iMu2]);
    out_muTwoIsGlobal =        (Float_t)((*muIsGlobal)[iMu2]);
    out_muTwoIsPF =    (Float_t)((*muIsPF)[iMu2]);
    out_muTwoIsSoft  = (Float_t) muonTwoIsSoft;

    out_MuMuMass =  (Float_t)(MuMuMass[iJPsi]);
    out_MuMuPx =    (Float_t)(MuMuPx[iJPsi]);
    out_MuMuPy =    (Float_t)(MuMuPy[iJPsi]);
    out_MuMuPz =    (Float_t)(MuMuPz[iJPsi]);
    out_MuMuVtx_CL =        (Float_t)(MuMuVtx_CL[iJPsi]);
    out_MuMuVtx_Chi2 =      (Float_t)(MuMuVtx_Chi2[iJPsi]);
    out_MuMuDecayVtx_X =    (Float_t)(MuMuDecayVtx_X[iJPsi]);
    out_MuMuDecayVtx_Y =    (Float_t)(MuMuDecayVtx_Y[iJPsi]);
    out_MuMuDecayVtx_Z =    (Float_t)(MuMuDecayVtx_Z[iJPsi]);
    out_MuMuDecayVtx_XE =   (Float_t)(MuMuDecayVtx_XE[iJPsi]);
    out_MuMuDecayVtx_YE =   (Float_t)(MuMuDecayVtx_YE[iJPsi]);
    out_MuMuDecayVtx_ZE =   (Float_t)(MuMuDecayVtx_ZE[iJPsi]);

    out_KKPx = Phi.Px();
    out_KKPy = Phi.Py();
    out_KKPz = Phi.Pz();
    out_KKPt = Phi.Pt();
    out_KKMass = Phi.M();

    out_XMass =     (Float_t)(XMass[iX]);
    out_XPx =       (Float_t)(XPx[iX]);
    out_XPy =       (Float_t)(XPy[iX]);
    out_XPz =       (Float_t)(XPz[iX]);

    out_XPxE =      (Float_t)(XPxE[iX]);
    out_XPyE =      (Float_t)(XPyE[iX]);
    out_XPzE =      (Float_t)(XPzE[iX]);


    out_XVtx_CL =   (Float_t)(XVtx_CL[iX]);
    out_XVtx_Chi2 =         (Float_t)(XVtx_Chi2[iX]);
    out_XDecayVtx_X =       (Float_t)(XDecayVtx_X[iX]);
    out_XDecayVtx_Y =       (Float_t)(XDecayVtx_Y[iX]);
    out_XDecayVtx_Z =       (Float_t)(XDecayVtx_Z[iX]);
    out_XDecayVtx_XE =      (Float_t)(XDecayVtx_XE[iX]);
    out_XDecayVtx_YE =      (Float_t)(XDecayVtx_YE[iX]);
    out_XDecayVtx_ZE =      (Float_t)(XDecayVtx_ZE[iX]);
    out_XCosAlphaBS =       (Float_t)(XCosAlphaBS[iX]);
    out_XCosAlpha3DBS =     (Float_t)(XCosAlpha3DBS[iX]);
    out_XCTauBS =   (Float_t)(XCTauBS[iX]);
    out_XCTauBSE =  (Float_t)(XCTauBSE[iX]);
    out_XLxyBS =    (Float_t)(XLxyBS[iX]);
    out_XLxyBSE =   (Float_t)(XLxyBSE[iX]);
    out_XLxyzBS =   (Float_t)(XLxyzBS[iX]);
    out_XLxyzBSE =  (Float_t)(XLxyzBSE[iX]);
    out_XCosAlphaPV =       (Float_t)(XCosAlphaPV[iX]);
    out_XCosAlpha3DPV =     (Float_t)(XCosAlpha3DPV[iX]);
    out_XCTauPV =   (Float_t)(XCTauPV[iX]);
    out_XCTauPVE =  (Float_t)(XCTauPVE[iX]);
    out_XLxyPV =    (Float_t)(XLxyPV[iX]);
    out_XLxyPVE =   (Float_t)(XLxyPVE[iX]);
    out_XLxyzPV =   (Float_t)(XLxyzPV[iX]);
    out_XLxyzPVE =  (Float_t)(XLxyzPVE[iX]);


    out_kaon1_dxy_PV =      (Float_t)(kaon1_dxy_PV[iX]);
    out_kaon1_dz_PV =       (Float_t)(kaon1_dz_PV[iX]);
    out_kaon2_dxy_PV =      (Float_t)(kaon2_dxy_PV[iX]);
    out_kaon2_dz_PV =       (Float_t)(kaon2_dz_PV[iX]);
    out_kaon1_dxy_BS =      (Float_t)(kaon1_dxy_BS[iX]);
    out_kaon1_dz_BS =       (Float_t)(kaon1_dz_BS[iX]);
    out_kaon2_dxy_BS =      (Float_t)(kaon2_dxy_BS[iX]);
    out_kaon2_dz_BS =       (Float_t)(kaon2_dz_BS[iX]);


    out_Kaon1FromPV =       (Float_t)((*Kaon1FromPV)[iX]);
    out_Kaon2FromPV =       (Float_t)((*Kaon2FromPV)[iX]);

    out_Muon1Px_MuMuKK =    (Float_t)(Muon1Px_MuMuKK[iX]);
    out_Muon1Py_MuMuKK =    (Float_t)(Muon1Py_MuMuKK[iX]);
    out_Muon1Pz_MuMuKK =    (Float_t)(Muon1Pz_MuMuKK[iX]);
    out_Muon1E_MuMuKK =     (Float_t)(Muon1E_MuMuKK[iX]);
    out_Muon2Px_MuMuKK =    (Float_t)(Muon2Px_MuMuKK[iX]);
    out_Muon2Py_MuMuKK =    (Float_t)(Muon2Py_MuMuKK[iX]);
    out_Muon2Pz_MuMuKK =    (Float_t)(Muon2Pz_MuMuKK[iX]);
    out_Muon2E_MuMuKK =     (Float_t)(Muon2E_MuMuKK[iX]);
    out_Kaon1Px_MuMuKK =    (Float_t)(Kaon1Px_MuMuKK[iX]);
    out_Kaon1Py_MuMuKK =    (Float_t)(Kaon1Py_MuMuKK[iX]);
    out_Kaon1Pz_MuMuKK =    (Float_t)(Kaon1Pz_MuMuKK[iX]);
    out_Kaon1E_MuMuKK =     (Float_t)(Kaon1E_MuMuKK[iX]);

    outTree->Fill();

    //
    // muonsCuts = muonQualityCut && muonChiCut && muonShitsCut && muonPhitsCut && muonDZPVCut && muonDXYPVCut;
    // muonsCuts = muonDZPVCut && muonDXYPVCut;
    //
    // jPsiPtCut      = (JPsi.Pt() > 7.0);
    // jPsiMassCut    =  (JPsiOriginal.M()<3.4 && JPsiOriginal.M()>2.8);
    // jPsiVtxCut     = ((*MuMuVtx_CL)[iJPsi]) > 0.1;
    // jPsiMuEtaPtCut = (fabs(mu1.Eta()) < 2.2) && (fabs(mu2.Eta()) < 2.2);
    // jPsiMusPtCut   = ((mu1.Pt() > 4.0) && (mu2.Pt() > 4.0));
    //
    // jPsiCuts = jPsiPtCut && jPsiMassCut && jPsiVtxCut  && jPsiMuEtaPtCut && jPsiMuEtaPtCut && jPsiMusPtCut;
    //
    // // kaonOneChiCut    = (((*trackChi2)[iK1] / (*trackNDF)[iK1]) < 5.0);
    // // kaonOnePhitsCut  = ((*trackPhits)[iK1] > 0);
    // // kaonOneShitsCut  = ((*trackShits)[iK1] >= 7);
    // // kaonTwoChiCut    = (((*trackChi2)[iK2] / (*trackNDF)[iK2]) < 5.0);
    // // kaonTwoPhitsCut  = ((*trackPhits)[iK2] > 0);
    // // kaonTwoShitsCut  = ((*trackShits)[iK2] >= 7);
    // kaonsPt = ((kaon1.Pt()>0.7) && (kaon2.Pt()>0.7));
    //
    // kaonOneCuts = kaonOneChiCut && kaonOnePhitsCut && kaonOneShitsCut;
    // kaonTwoCuts = kaonTwoChiCut && kaonTwoPhitsCut && kaonTwoShitsCut;
    // kaonsCuts = kaonOneCuts && kaonTwoCuts && kaonsPt;
    // kaonsCuts = kaonsPt;
    //
    // cosAlphaCut = (fabs((*XCosAlphaPV)[iX]) > 0.99);
    // vtxCLCut =  (((*XVtx_CL)[iX]) > 0.01);
    //
    // extraCuts = vtxCLCut && cosAlphaCut;

  }

  return kTRUE;
}

void TwoMuTwoK_2012::SlaveTerminate()
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

void TwoMuTwoK_2012::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

}
