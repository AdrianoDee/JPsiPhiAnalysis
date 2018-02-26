#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TDirectory.h>
#include <iostream>
#include <bitset>
#include <TCanvas.h>
#include <algorithm>
#include <TLegend.h>
#include <TStyle.h>
#include <string>
#include <TColor.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <vector>

int noHlts = 13;

double pi = 3.14159265358979323846;
double pdg_Phi_mass = 1.019455;

std::string hltsName[13] = {
  "HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi", //0
"HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi", //1
"HLT_Mu20_TkMu0_Phi", //2
"HLT_Dimuon14_Phi_Barrel_Seagulls", //3
"HLT_Mu25_TkMu0_Phi", //4
"HLT_Dimuon24_Phi_noCorrL1", //5
"HLT_DoubleMu4_JpsiTrkTrk_Displaced", //6
"HLT_DoubleMu4_JpsiTrk_Displaced", //7
"HLT_DoubleMu4_Jpsi_Displaced", //8
"HLT_DoubleMu4_3_Jpsi_Displaced", //9
"HLT_Dimuon20_Jpsi_Barrel_Seagulls", //10
"HLT_Dimuon25_Jpsi", //11
"HLT_Dimuon0_Jpsi" //12
};


int skimXTree(std::string path, std::string filename, std::string treename = "xTree", std::string dirname = "rootuple")
{

  TFile *oldfile = TFile::Open((path+filename).data());
  TDirectory *directory = (TDirectory*)oldfile->Get(dirname.data());
  TTree *oldtree = (TTree*)directory->Get(treename.data());
  Long64_t nentries = oldtree->GetEntries();
  ULong64_t event   = 0;
  oldtree->SetBranchAddress("event",&event);
  //Create a new file + a clone of old tree in new file
  TFile *newfile = new TFile((treename + "_skim_" + filename).data(),"RECREATE");
  TTree *newtree = oldtree->CloneTree();

  newtree->Print();
  newtree->Write();

  return 0;

}



int selectXTree()
{

  TFile *oldfile = TFile::Open("/Users/adrianodiflorio/Documents/Git/X4140/iPythons/skimmedNP.root");

  TTree *oldtree = (TTree*)oldfile->Get("jTree");
  Long64_t nentries = oldtree->GetEntries();
  Double_t xyl   = 0.0;
  Double_t xylErr   = 0.0;
  Double_t cosA  = 0.0;

  Double_t phiM  = 0.0;
  Double_t jPsiM  = 0.0;

  Double_t ctau  = 0.0;
  Double_t ctauErr  = 0.0;

  Double_t vProb  = 0.0;

  Int_t phiMType = 0, phiPType = 0;
  UInt_t phi_trigger = 0, jpsi_trigger = 0;

  oldtree->SetBranchAddress("vProb",&vProb);
  oldtree->SetBranchAddress("l_xy",&xyl);
  oldtree->SetBranchAddress("lErr_xy",&xylErr);
  oldtree->SetBranchAddress("cosAlpha",&cosA);
  oldtree->SetBranchAddress("phi_M",&phiM);
  oldtree->SetBranchAddress("jpsi_M",&jPsiM);

  oldtree->SetBranchAddress("phi_muonM_type",&phiMType);
  oldtree->SetBranchAddress("phi_muonP_type",&phiPType);

  oldtree->SetBranchAddress("ctauPV",&ctau);
  oldtree->SetBranchAddress("ctauErrPV",&ctauErr);

  oldtree->SetBranchAddress("phi_trigger",&phi_trigger);
  oldtree->SetBranchAddress("jpsi_trigger",&jpsi_trigger);
  //Create a new file + a clone of old tree in new file
  TFile *newfile = new TFile("skimmedNPCos.root","RECREATE");
  TTree *newtree = oldtree->CloneTree(0);

  for (Long64_t i=0;i<nentries; i++) {
    oldtree->GetEntry(i);
    std::bitset<16> pM(phiMType);
    std::bitset<16> pP(phiPType);
    // std::cout << phiMType << "-" << binary << std::endl;
    // std::cout<<vProb <<std::endl;
    //if (vProb > 0.0) newtree->Fill();
    if (jPsiM > 2.8 && phiM < 1.1 && cosA > 0.995 && phiM > 0.95 && pP.test(1) && pM.test(1) && vProb > 0.1  ) newtree->Fill();
  }

  newtree->Print();
  newtree->Write();

  return 0;

}


int drawPTree(std::string path = "/Users/adrianodiflorio/Documents/Git/X4140/iPythons/xTree.root",std::string treename = "xTree")
{

  UInt_t colors[13] = {1,2,3,6,7,8,30,40,46,38,29,34,9};

  TFile *oldfile = TFile::Open(path.data());
  TTree *oldtree = (TTree*)oldfile->Get(treename.data());

  Long64_t nentries = oldtree->GetEntries();

  Double_t xyl   = 0.0;
  Double_t xylErr   = 0.0;
  Double_t cosA  = 0.0;

  Double_t phiM  = 0.0;
  Double_t jPsiM  = 0.0;

  Double_t ctau  = 0.0;
  Double_t ctauErr  = 0.0;

  Double_t vProb  = 0.0;

  Int_t phiMType = 0, phiPType = 0;
  UInt_t jpsi_trigger = 0, trigger = 0;
  Int_t phi_trigger = 0;

  Int_t isTrackerM, isTrackerP;

  oldtree->SetBranchAddress("pM",&phiM);
  oldtree->SetBranchAddress("p_vProb",&vProb);
  oldtree->SetBranchAddress("trigger",&trigger);
  oldtree->SetBranchAddress("p_triggerMatch",&phi_trigger);

  oldtree->SetBranchAddress("p_muonP_type",&phiPType);
  oldtree->SetBranchAddress("p_muonM_type",&phiMType);

  oldtree->SetBranchAddress("p_muonP_isTracker",&isTrackerP);
  oldtree->SetBranchAddress("p_muonM_isTracker",&isTrackerM);

  //Create a new file + a clone of old tree in new file
  TCanvas c("c","c",1200,1600);

  TFile *newfile = new TFile("drawSkim.root","RECREATE");
  // for(int j = 0; j < 1; j++)
  // {


  Double_t xmin = 5.0, xmax = 6.0;
  Int_t xBin = ((xmax - xmin)/0.01);

  TTree *newtree = oldtree->CloneTree(0);
  // TH1F* phi_triggrHist = new TH1F("phi_triggrHist","phi_triggrHist",600,0.6,1.2);
  TH1F* phiHist = new TH1F("phiHist","phiHist",250,0.0,1.25);

  std::vector<TH1F*> phiHists;


  for (int i = 0; i < 13; i++)
  phiHists.push_back(new TH1F((hltsName[i] + "_phi").data(),(hltsName[i] + "_phi").data(),200,0.25,1.25));


  for (Long64_t i=0;i<nentries; i++) {
    oldtree->GetEntry(i);
    std::bitset<16> tB(trigger);
    std::bitset<16> pM(phiMType);
    std::bitset<16> pP(phiPType);
    bool test = false;
    // bool jpsimass = jPsiM < 3.15 && jPsiM > 3.0;
    // bool phimass = phiM < 1.06 && phiM > 0.98;
    for (int j = 0; j < 13; j++){
      // if (tB.test(j) && cosA > 0.995 && vProb > 0.01 && xyl/xylErr > 2.0 && trigger > 0)
      if (tB.test(j) && vProb > 0.05 && phi_trigger > 0 )
      {
        test = true;
        phiHists[j]->Fill(phiM);

      }
    }
    // if(cosA > 0.995 && phimass && jpsimass && vProb>0.02 && xyl/xylErr > 2.0)
    if (test)
    phiHist->Fill(phiM);

  }


  //newtree->Draw("phi_M","","same");
  // }
  phiHist->SetMinimum(1.0);
  phiHist->SetMaximum(phiHist->GetMaximum()*5.0);
  //oldtree->Draw("phi_M");
  // TH1F* phi_triggrHist = (TH1F*)gDirectory->Get("phi_triggrHist");

  phiHist->SetLineColor(kBlue);
  phiHist->Write();
  phiHist->Draw();

  TLegend* leg = new TLegend(0.1,0.5,0.45,0.9);
  leg->AddEntry(phiHist,(phiHist->GetName()),"l");
  for (int i = 0; i < 13; i++)
  {
    phiHists[i]->SetLineColor(colors[i]);
    phiHists[i]->SetLineWidth(2);
    if(i>5) phiHists[i]->SetLineStyle(kDashed);
    phiHists[i]->Draw("same");
    leg->AddEntry(phiHists[i],(phiHists[i]->GetName()),"l");
    phiHists[i]->Write();
  }

  leg->Draw();
  c.SetLogy(1);
  c.SaveAs("phitriggerCheck.png");
  c.SaveAs("phitriggerCheck.eps");
  c.SaveAs("phitriggerCheck.root");

  return 0;

}


int drawXXTree(std::string path = "/Users/adrianodiflorio/Documents/Git/X4140/iPythons/xTree.root",std::string treename = "xTree")
{

  UInt_t colors[13] = {1,2,3,6,7,8,30,40,46,38,29,34,9};

  TFile *oldfile = TFile::Open(path.data());
  TTree *oldtree = (TTree*)oldfile->Get(treename.data());

  Long64_t nentries = oldtree->GetEntries();

  Double_t xyl   = 0.0;
  Double_t xylErr   = 0.0;
  Double_t cosA  = 0.0;

  Double_t xM  = 0.0;
  Double_t jPsiM  = 0.0, phiM = 0.0;

  Double_t ctau  = 0.0;
  Double_t ctauErr  = 0.0;

  Double_t vProb  = 0.0;

  Int_t xMType = 0, xPType = 0;
  UInt_t trigger = 0;
  UInt_t phi_trigger = 0, jpsi_trigger = 0;

  Int_t isTrackerM, isTrackerP;

  oldtree->SetBranchAddress("xM",&xM);
  oldtree->SetBranchAddress("vProb",&vProb);
  oldtree->SetBranchAddress("trigger",&trigger);
  oldtree->SetBranchAddress("phi_trigger",&phi_trigger);
  oldtree->SetBranchAddress("jpsi_trigger",&jpsi_trigger);

  oldtree->SetBranchAddress("phi_muonP_type",&xPType);
  oldtree->SetBranchAddress("phi_muonM_type",&xMType);

  // oldtree->SetBranchAddress("p_muonP_isTracker",&isTrackerP);
  // oldtree->SetBranchAddress("p_muonM_isTracker",&isTrackerM);

  oldtree->SetBranchAddress("cosAlpha",&isTrackerP);
  oldtree->SetBranchAddress("p_muonP_isTracker",&isTrackerP);
  oldtree->SetBranchAddress("p_muonP_isTracker",&isTrackerP);

  oldtree->SetBranchAddress("l_xy",&xyl);
  oldtree->SetBranchAddress("lErr_xy",&xylErr);
  oldtree->SetBranchAddress("cosAlpha",&cosA);
  oldtree->SetBranchAddress("phi_M",&phiM);
  oldtree->SetBranchAddress("jpsi_M",&jPsiM);

  //Create a new file + a clone of old tree in new file
  TCanvas c("c","c",1200,1600);

  TFile *newfile = new TFile("drawSkim.root","RECREATE");
  // for(int j = 0; j < 1; j++)
  // {


  Double_t xmin = 5.0, xmax = 6.0;
  Int_t xBin = ((xmax - xmin)/0.01);

  TTree *newtree = oldtree->CloneTree(0);
  // TH1F* x_triggrHist = new TH1F("x_triggrHist","x_triggrHist",600,0.6,1.2);
  TH1F* xHist = new TH1F("xHist","xHist",600,3.9,6.1);

  std::vector<TH1F*> xHists;

  for (int i = 0; i < 13; i++)
  xHists.push_back(new TH1F((hltsName[i] + "_x").data(),(hltsName[i] + "_x").data(),600,3.9,6.1));



  for (Long64_t i=0;i<nentries; i++) {
    oldtree->GetEntry(i);
    std::bitset<16> tB(trigger);
    // std::bitset<16> pM(xMType);
    // std::bitset<16> pP(xPType);
    bool test = false;
    // bool jpsimass = jPsiM < 3.15 && jPsiM > 3.0;
    // bool xmass = xM < 1.06 && xM > 0.98;
    for (int j = 0; j < 13; j++){
      if (tB.test(j) && cosA > 0.995 && vProb > 0.01 && xyl/xylErr > 2.0 && trigger > 0)
      //if (xM < 5.4 && xM > 5.3 && tB.test(j) && vProb > 0.1 )
      {
        test = true;
        xHists[j]->Fill(xM);
      }
    }
    // if(cosA > 0.995 && xmass && jpsimass && vProb>0.02 && xyl/xylErr > 2.0)
    if (test)
    xHist->Fill(xM);

  }


  //newtree->Draw("x_M","","same");
  // }
  xHist->SetMinimum(1.0);
  xHist->SetMaximum(xHist->GetMaximum()*5.0);
  //oldtree->Draw("x_M");
  // TH1F* x_triggrHist = (TH1F*)gDirectory->Get("x_triggrHist");

  xHist->SetLineColor(kBlue);
  xHist->Write();
  xHist->Draw();

  TLegend* leg = new TLegend(0.1,0.5,0.45,0.9);
  leg->AddEntry(xHist,(xHist->GetName()),"l");
  for (int i = 0; i < 13; i++)
  {
    xHists[i]->SetLineColor(colors[i]);
    xHists[i]->SetLineWidth(2);
    if(i>5) xHists[i]->SetLineStyle(kDashed);
    xHists[i]->Draw("same");
    leg->AddEntry(xHists[i],(xHists[i]->GetName()),"l");
    xHists[i]->Write();
  }

  leg->Draw();
  c.SetLogy(1);
  c.SaveAs("xtriggerCheck.png");
  c.SaveAs("xtriggerCheck.eps");
  c.SaveAs("xtriggerCheck.root");

  return 0;

}



int drawXTree(std::string path = "/Users/adrianodiflorio/Documents/Git/X4140/iPythons/xTree.root",std::string treename = "xTree")
{

  UInt_t colors[13] = {1,2,3,6,7,8,30,40,46,38,29,34,9};

  TFile *oldfile = TFile::Open(path.data());
  TTree *oldtree = (TTree*)oldfile->Get(treename.data());

  Long64_t nentries = oldtree->GetEntries();
  Double_t xM = 0.0;
  Double_t xyl   = 0.0;
  Double_t xylErr   = 0.0;
  Double_t cosA  = 0.0;

  Double_t phiM  = 0.0;
  Double_t jPsiM  = 0.0;

  Double_t ctau  = 0.0;
  Double_t ctauErr  = 0.0;

  Double_t vProb  = 0.0;

  Int_t phiMType = 0, phiPType = 0;
  UInt_t phi_trigger = 0, jpsi_trigger = 0, trigger = 0;

  TLorentzVector *xP4 = 0, *jP4 = 0, *pP4 = 0;
  TLorentzVector *mM_jpsi_P4 = 0, *mP_jpsi_P4 = 0, *mM_phi_P4 = 0, *mP_phi_P4 = 0;

  oldtree->SetBranchAddress("xM",&xM);
  oldtree->SetBranchAddress("vProb",&vProb);
  oldtree->SetBranchAddress("trigger",&trigger);
  oldtree->SetBranchAddress("l_xy",&xyl);
  oldtree->SetBranchAddress("lErr_xy",&xylErr);
  oldtree->SetBranchAddress("cosAlpha",&cosA);
  oldtree->SetBranchAddress("phi_M",&phiM);
  oldtree->SetBranchAddress("jpsi_M",&jPsiM);

  oldtree->SetBranchAddress("phi_muonM_type",&phiMType);
  oldtree->SetBranchAddress("phi_muonP_type",&phiPType);

  oldtree->SetBranchAddress("ctauPV",&ctau);
  oldtree->SetBranchAddress("ctauErrPV",&ctauErr);

  oldtree->SetBranchAddress("phi_trigger",&phi_trigger);
  oldtree->SetBranchAddress("jpsi_trigger",&jpsi_trigger);

  oldtree->SetBranchAddress("x_p4",&xP4);

  oldtree->SetBranchAddress("phi_p4",&pP4);
  oldtree->SetBranchAddress("muonM_phi_p4",&mM_jpsi_P4);
  oldtree->SetBranchAddress("muonP_phi_p4",&mP_jpsi_P4);

  oldtree->SetBranchAddress("jpsi_p4",&jP4);
  oldtree->SetBranchAddress("muonM_jpsi_p4",&mM_phi_P4);
  oldtree->SetBranchAddress("muonP_jpsi_p4",&mP_phi_P4);
  //Create a new file + a clone of old tree in new file
  TCanvas c("c","c",1200,1600);

  TFile *newfile = new TFile("drawSkim.root","RECREATE");
  // for(int j = 0; j < 1; j++)
  // {


  Double_t xmin = 4.0, xmax = 6.0;
  Int_t xBin = ((xmax - xmin)/0.01);

  TTree *newtree = oldtree->CloneTree(0);
  // TH1F* phi_triggrHist = new TH1F("phi_triggrHist","phi_triggrHist",600,0.6,1.2);
  TH1F* phiHist = new TH1F("phiHist","phiHist",250,0.0,1.25);
  TH1F* jpsiHist = new TH1F("jpsiHist","jpsiHist",140,2.6,3.3);
  TH1F* xHist = new TH1F("xHist","xHist",xBin,xmin,xmax);

  TH1F* x_ptHist = new TH1F("x_ptHist","x_ptHist",1000,0.0,100.0);
  TH1F* jpsi_ptHist = new TH1F("jpsi_ptHist","jpsi_ptHist",1000,0.0,100.0);
  TH1F* phi_ptHist = new TH1F("phi_ptHist","phi_ptHist",1000,0.0,100.0);

  TH1F* jpsiMP_ptHist = new TH1F("jpsiMP_ptHist","jpsiMP_ptHist",1000,0.0,100.0);
  TH1F* jpsiMM_ptHist = new TH1F("jpsiMM_ptHist","jpsiMM_ptHist",1000,0.0,100.0);
  TH1F* jpsiMHig_ptHist = new TH1F("jpsiMHig_ptHist","jpsiMHig_ptHist",1000,0.0,100.0);
  TH1F* jpsiMLow_ptHist = new TH1F("jpsiMLow_ptHist","jpsiMLow_ptHist",1000,0.0,100.0);

  TH1F* phiMP_ptHist = new TH1F("phiMP_ptHist","phiMP_ptHist",1000,0.0,100.0);
  TH1F* phiMM_ptHist = new TH1F("phiMM_ptHist","phiMM_ptHist",1000,0.0,100.0);
  TH1F* phiMHig_ptHist = new TH1F("phiMHig_ptHist","phiMHig_ptHist",1000,0.0,100.0);
  TH1F* phiMLow_ptHist = new TH1F("phiMLow_ptHist","phiMLow_ptHist",1000,0.0,100.0);

  TH2F* phiPts = new TH2F("phiPts","phiPts",1000,0.0,100.0,1000,0.0,100.0);
  TH2F* jpsPts = new TH2F("jpsPts","jpsPts",1000,0.0,100.0,1000,0.0,100.0);

  TH1F* dRJpsiPhi = new TH1F("dRJpsiPhi","dRJpsiPhi",1000,-10.0,10.0);

  std::vector<TH1F*> phiHists;
  std::vector<TH1F*> jpsiHists;
  std::vector<TH1F*> xHists;


  for (int i = 0; i < noHlts; i++)
  {
    phiHists.push_back(new TH1F((hltsName[i] + "_phi").data(),(hltsName[i] + "_phi").data(),500,0.25,1.25));
    jpsiHists.push_back(new TH1F((hltsName[i] + "_jpsi").data(),(hltsName[i] + "_jpsi").data(),140,2.6,3.3));
    xHists.push_back(new TH1F((hltsName[i] + "_x").data(),(hltsName[i] + "_x").data(),xBin,xmin,xmax));
  }


  for (Long64_t i=0;i<nentries; i++) {
    oldtree->GetEntry(i);
    std::bitset<16> tB(trigger);
    std::bitset<16> pM(phiMType);
    std::bitset<16> pP(phiPType);

    float deltaEta = jP4->Eta() - pP4->Eta();
    float deltaPhi = std::fabs(jP4->Phi() - pP4->Phi());
    if (deltaPhi>pi) deltaPhi-=pi;
    float deltaR = sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);

    dRJpsiPhi->Fill(deltaR); //cut < 1

    x_ptHist->Fill(xP4->Pt());
    jpsi_ptHist->Fill(jP4->Pt());

    jpsiMP_ptHist->Fill(mP_jpsi_P4->Pt());
    jpsiMM_ptHist->Fill(mM_jpsi_P4->Pt());

    jpsiMHig_ptHist->Fill(std::max(mP_jpsi_P4->Pt(),mM_jpsi_P4->Pt()));
    jpsiMLow_ptHist->Fill(-std::max(-mP_jpsi_P4->Pt(),-mM_jpsi_P4->Pt()));

    phiMHig_ptHist->Fill(std::max(mP_phi_P4->Pt(),mM_phi_P4->Pt()));
    phiMLow_ptHist->Fill(-std::max(-mP_phi_P4->Pt(),-mM_phi_P4->Pt()));

    phiMP_ptHist->Fill(mP_phi_P4->Pt());
    phiMM_ptHist->Fill(mM_phi_P4->Pt());

    phiPts->Fill(-std::max(-mP_phi_P4->Pt(),-mM_phi_P4->Pt()),std::max(mP_phi_P4->Pt(),mM_phi_P4->Pt()));
    jpsPts->Fill(-std::max(-mP_jpsi_P4->Pt(),-mM_jpsi_P4->Pt()), std::max(mP_jpsi_P4->Pt(),mM_jpsi_P4->Pt()));

    bool test = false;
    bool jpsimass = jPsiM < 3.2 && jPsiM > 3.0;
    bool phimass = phiM > 1.005 && phiM < 1.03;
    for (int j = 0; j < 13; j++){
      // if (tB.test(j) && cosA > 0.995 && vProb > 0.01 && xyl/xylErr > 2.0 && trigger > 0)
      // if (tB.test(j))
      if (xyl/xylErr > 0.0 && xP4->Pt() > 6.0 && cosA > 0.997 && pP4->Pt() > 8.0 && jP4->Pt() > 5.0 && mP_phi_P4->Pt() > 2.0 && mM_phi_P4->Pt() > 2.0 && tB.test(j) && vProb > 0.05 && deltaR < 1.0 && deltaR > 0.0 && jpsimass && phimass)
      {
        phi_ptHist->Fill(pP4->Pt());
        test = true;
        phiHists[j]->Fill(phiM);
        jpsiHists[j]->Fill(jPsiM);
        //  if(phimass && jpsimass && cosA > 0.995 && vProb>0.02 && xyl/xylErr > 2.0)

        xHists[j]->Fill(xM);
      }
    }
    // if(cosA > 0.995 && phimass && jpsimass && vProb>0.02 && xyl/xylErr > 2.0)
    if (test){
      xHist->Fill(xM);
      phiHist->Fill(phiM);
      jpsiHist->Fill(jPsiM);
    }
  }

  dRJpsiPhi->Write();
  x_ptHist->Write();
  jpsi_ptHist->Write();
  phi_ptHist->Write();
  jpsiMP_ptHist->Write();
  jpsiMM_ptHist->Write();
  jpsiMHig_ptHist->Write();
  jpsiMLow_ptHist->Write();
  phiMHig_ptHist->Write();
  phiMLow_ptHist->Write();
  phiMP_ptHist->Write();
  phiMM_ptHist->Write();

  phiPts->Write();
  jpsPts->Write();

  //newtree->Draw("phi_M","","same");
  // }
  phiHist->SetMinimum(1.0);
  phiHist->SetMaximum(phiHist->GetMaximum()*5.0);
  //oldtree->Draw("phi_M");
  // TH1F* phi_triggrHist = (TH1F*)gDirectory->Get("phi_triggrHist");

  phiHist->SetLineColor(kBlue);
  phiHist->Write();
  phiHist->Draw();

  TLegend* leg = new TLegend(0.1,0.5,0.45,0.9);
  leg->AddEntry(phiHist,(phiHist->GetName()),"l");
  for (int i = 0; i < 13; i++)
  {
    phiHists[i]->SetLineColor(colors[i]);
    phiHists[i]->SetLineWidth(2);
    if(i>5) phiHists[i]->SetLineStyle(kDashed);
    phiHists[i]->Draw("same");
    leg->AddEntry(phiHists[i],(phiHists[i]->GetName()),"l");
    phiHists[i]->Write();
  }

  leg->Draw();
  c.SetLogy(1);
  c.SaveAs("phitriggerCheck.png");
  c.SaveAs("phitriggerCheck.eps");
  c.SaveAs("phitriggerCheck.root");

  jpsiHist->SetMinimum(1.0);
  jpsiHist->SetMaximum(jpsiHist->GetMaximum()*5.0);
  //oldtree->Draw("phi_M");
  // TH1F* phi_triggrHist = (TH1F*)gDirectory->Get("phi_triggrHist");

  jpsiHist->SetLineColor(kBlue);
  jpsiHist->Write();
  jpsiHist->Draw();

  leg = new TLegend(0.1,0.5,0.4,0.9);
  leg->AddEntry(jpsiHist,(jpsiHist->GetName()),"l");
  for (int i = 0; i < 13; i++)
  {
    jpsiHists[i]->SetLineColor(colors[i]);
    jpsiHists[i]->SetLineWidth(2);
    if(i>5) jpsiHists[i]->SetLineStyle(kDashed);
    jpsiHists[i]->Draw("same");
    leg->AddEntry(jpsiHists[i],(jpsiHists[i]->GetName()),"l");
    jpsiHists[i]->Write();
  }


  // phi_triggrHist->Draw("same");
  leg->Draw();
  c.SetLogy(1);
  c.SaveAs("jpsitriggerCheck.png");
  c.SaveAs("jpsitriggerCheck.eps");
  c.SaveAs("jpsitriggerCheck.root");

  xHist->SetMinimum(0.1);
  xHist->SetMaximum(xHist->GetMaximum());//*5.0);

  TLine line (5.35,0.1,5.35,xHist->GetMaximum());
  line.SetLineColor(kRed);
  line.SetLineWidth(2);

  xHist->SetLineColor(kBlue);
  xHist->Write();
  xHist->Draw();

  leg = new TLegend(0.1,0.55,0.35,0.9);
  leg->AddEntry(xHist,(xHist->GetName()),"l");
  for (int i = 0; i < 13; i++)
  {
    xHists[i]->SetLineColor(colors[i]);
    xHists[i]->SetLineWidth(2);
    if(i>5) xHists[i]->SetLineStyle(kDashed);
    xHists[i]->Draw("same");
    leg->AddEntry(xHists[i],(xHists[i]->GetName()),"l");
    xHists[i]->Write();
  }

  // line.Draw();
  // phi_triggrHist->Draw("same");
  leg->Draw();
  c.SetLogy(0);
  c.SaveAs("xtriggerCheck.png");
  c.SaveAs("xtriggerCheck.eps");
  c.SaveAs("xtriggerCheck.root");

  return 0;

}

int jpsiRuns(std::string path, std::string filename, std::string treename, std::string dirname)
{
  UInt_t colors[13] = {1,2,3,6,7,8,30,40,46,38,29,34,9};

  TFile *oldfile = TFile::Open((path+"/"+filename).data());
  TDirectory *directory = (TDirectory*)oldfile->Get(dirname.data());
  TTree *oldtree = (TTree*)directory->Get((treename).data());


  Long64_t nentries = oldtree->GetEntries();
  Double_t xM = 0.0, xDeltaM = 0.0;

  Float_t cosA  = 0.0;

  Double_t jPsiM  = 0.0;

  Float_t ctau  = 0.0;
  Float_t ctauErr  = 0.0;

  Float_t vProb  = 0.0;
  Double_t deltaR = 0.0;

  UInt_t n_jpsi = 0;
  UInt_t trigger = 0;
  UInt_t run = 0;

  TLorentzVector *xP4 = 0, *jP4 = 0, *pP4 = 0;
  TLorentzVector *muonp_p4 = 0, *muonn_p4 = 0, *kaonp_p4 = 0, *kaonn_p4 = 0;

  //Output Variables
  Int_t t = 0, evt = 0, r = 0;
  Double_t jPt, xPt, mNPt, mPPt, kNPt, kPPt, cosAlpha, vP, pPt;
  Double_t xyl   = 0.0;
  Double_t xylErr   = 0.0;

  oldtree->SetBranchAddress("vProb",&vProb);
  oldtree->SetBranchAddress("trigger",&trigger);
  oldtree->SetBranchAddress("run",&run);
  oldtree->SetBranchAddress("nonia",&n_jpsi);
  oldtree->SetBranchAddress("cosAlpha",&cosA);
  oldtree->SetBranchAddress("ppdlPV",&ctau);
  oldtree->SetBranchAddress("ppdlErrPV",&ctauErr);
  oldtree->SetBranchAddress("dimuon_p4",&jP4);

  // std::string hltsName[13] = {
  //   "HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi", //0
  // "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi", //1
  // "HLT_Mu20_TkMu0_Phi", //2
  // "HLT_Dimuon14_Phi_Barrel_Seagulls", //3
  // "HLT_Mu25_TkMu0_Phi", //4
  // "HLT_Dimuon24_Phi_noCorrL1", //5
  // "HLT_DoubleMu4_JpsiTrkTrk_Displaced", //6
  // "HLT_DoubleMu4_JpsiTrk_Displaced", //7
  // "HLT_DoubleMu4_Jpsi_Displaced", //8
  // "HLT_DoubleMu4_3_Jpsi_Displaced", //9
  // "HLT_Dimuon20_Jpsi_Barrel_Seagulls", //10
  // "HLT_Dimuon25_Jpsi", //11
  // "HLT_Dimuon0_Jpsi" //12
  // };

  std::vector <int> triggersToTest;

  //OUR TRIGGERs
  triggersToTest.push_back(0);
  triggersToTest.push_back(1);

  //PHI INCLUSIVE TRIGGERs
  triggersToTest.push_back(2);
  triggersToTest.push_back(3);
  triggersToTest.push_back(4);
  triggersToTest.push_back(5);

  //JPSI INCLUSIVE TRIGGERs
  triggersToTest.push_back(10);
  triggersToTest.push_back(11);
  triggersToTest.push_back(12);

  //DISPLACED JPSI TRIGGERs
  triggersToTest.push_back(6); //dis
  triggersToTest.push_back(7); //dis
  triggersToTest.push_back(8); //dis
  triggersToTest.push_back(9); //dis



  //oldtree->SetBranchAddress("oniat_rf_p4",&xP4);

  TFile *newfile = new TFile((treename + "_jpsisRun_" + filename).data(),"RECREATE");

  std::vector <TH1F*> JPsi_vs_run_hists;

  for (size_t i = 0; i < triggersToTest.size(); i++)
    JPsi_vs_run_hists.push_back(new TH1F (("JPsi_vs_run_" + std::to_string(triggersToTest[i])).data(), "JPsi_vs_run; Run[#];J/Psi[#]",40000, 280000, 320000));

  TH1F* JPsi_vs_run = new TH1F ("JPsi_vs_run", "JPsi_vs_run; Run[#];J/Psi[#]",20000, 190000, 210000);



  for (Long64_t i=0;i<nentries; i++) {

    int barWidth = 70;

    float progress = float(i)/float(nentries);

    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();

    oldtree->GetEntry(i);

    std::bitset<16> tB(trigger);

    bool tested = false;

    for (size_t j = 0; j < triggersToTest.size(); j++){

      int testingTrigger = triggersToTest[j];

      if (tB.test(testingTrigger) && vProb > 0.0)
      {
        tested = true;
        JPsi_vs_run_hists[j]->Fill(run);
      }
    }

      if(tested)
        JPsi_vs_run->Fill(run);

  }

  for (size_t j = 0; j < triggersToTest.size(); j++)
    JPsi_vs_run_hists[j]->Write();

  JPsi_vs_run->Write();

  return 0;

}


int jpsiRunsMMKK(std::string path, std::string filename, std::string treename, std::string dirname)
{
  UInt_t colors[13] = {1,2,3,6,7,8,30,40,46,38,29,34,9};

  TFile *oldfile = TFile::Open((path+"/"+filename).data());
  // TDirectory *directory = (TDirectory*)oldfile->Get(dirname.data());
  TTree *oldtree = (TTree*)oldfile->Get((treename).data());


  Long64_t nentries = oldtree->GetEntries();
  Double_t xM = 0.0, xDeltaM = 0.0;

  Double_t cosA  = 0.0;

  Double_t jPsiM  = 0.0, phiM = 0.0;

  Float_t ctau  = 0.0;
  Float_t ctauErr  = 0.0;

  Double_t vProb  = 0.0;
  Double_t deltaR = 0.0;

  Int_t trigger = 0;
  Int_t run = 0;
  Int_t ev = 0;

  TLorentzVector *xP4 = 0, *jP4 = 0, *pP4 = 0;
  TLorentzVector *muonp_p4 = 0, *muonn_p4 = 0, *kaonp_p4 = 0, *kaonn_p4 = 0;

  //Output Variables
  Double_t jPt, xPt, mNPt, mPPt, kNPt, kPPt, cosAlpha, vP, pPt;
  Double_t xyl   = 0.0;
  Double_t xylErr   = 0.0;

  oldtree->SetBranchAddress("trigger",&trigger);
  oldtree->SetBranchAddress("run",&run);
  oldtree->SetBranchAddress("event",&ev);
  // oldtree->SetBranchAddress("nonia",&n_jpsi);
  // oldtree->SetBranchAddress("cosAlpha",&cosA);
  // oldtree->SetBranchAddress("ppdlPV",&ctau);
  // oldtree->SetBranchAddress("ppdlErrPV",&ctauErr);
  oldtree->SetBranchAddress("psi_p4",&jP4);

  oldtree->SetBranchAddress("oniat_ctauPV",&xyl);
  oldtree->SetBranchAddress("oniat_ctauErrPV",&xylErr);
  oldtree->SetBranchAddress("oniat_cosAlpha",&cosA);
  oldtree->SetBranchAddress("oniat_vProb",&vProb);

  oldtree->SetBranchAddress("phi_p4",&pP4);
  oldtree->SetBranchAddress("oniat_rf_p4",&xP4);
  // std::string hltsName[13] = {
  //   "HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi", //0
  // "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi", //1
  // "HLT_Mu20_TkMu0_Phi", //2
  // "HLT_Dimuon14_Phi_Barrel_Seagulls", //3
  // "HLT_Mu25_TkMu0_Phi", //4
  // "HLT_Dimuon24_Phi_noCorrL1", //5
  // "HLT_DoubleMu4_JpsiTrkTrk_Displaced", //6
  // "HLT_DoubleMu4_JpsiTrk_Displaced", //7
  // "HLT_DoubleMu4_Jpsi_Displaced", //8
  // "HLT_DoubleMu4_3_Jpsi_Displaced", //9
  // "HLT_Dimuon20_Jpsi_Barrel_Seagulls", //10
  // "HLT_Dimuon25_Jpsi", //11
  // "HLT_Dimuon0_Jpsi" //12
  // };

  std::vector <int> triggersToTest;

  //OUR TRIGGERs
  // triggersToTest.push_back(0);
  triggersToTest.push_back(1);

  // //PHI INCLUSIVE TRIGGERs
  // triggersToTest.push_back(2);
  // triggersToTest.push_back(3);
  // triggersToTest.push_back(4);
  // triggersToTest.push_back(5);

  //JPSI INCLUSIVE TRIGGERs
  triggersToTest.push_back(10);
  triggersToTest.push_back(11);
  triggersToTest.push_back(12);

  //DISPLACED JPSI TRIGGERs
  triggersToTest.push_back(6); //dis
  triggersToTest.push_back(7); //dis
  triggersToTest.push_back(8); //dis
  triggersToTest.push_back(9); //dis

  //oldtree->SetBranchAddress("oniat_rf_p4",&xP4);

  TFile *newfile = new TFile((treename + "_jpsisRun_" + filename).data(),"RECREATE");

  std::vector <TH1F*> JPsi_vs_run_hists;

  for (size_t i = 0; i < triggersToTest.size(); i++)
    JPsi_vs_run_hists.push_back(new TH1F (("JPsi_vs_run_" + std::to_string(triggersToTest[i])).data(), "JPsi_vs_run; Run[#];J/Psi[#]",40000, 280000, 320000));

  TH1F* JPsi_vs_run = new TH1F ("JPsi_vs_run", "JPsi_vs_run; Run[#];J/Psi[#]",40000, 280000, 320000);

  std::map<Int_t,int> eventMap;
  std::map<Int_t,int> runMap;

  for (Long64_t i=0;i<nentries; i++) {

    int barWidth = 70;

    float progress = float(i)/float(nentries);

    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();

    oldtree->GetEntry(i);

    if(eventMap.find(ev) != eventMap.end())
      if(runMap.find(run) != runMap.end())
        continue;

    std::bitset<16> tB(trigger);

    bool tested = false;
    bool phimass = pP4->M() > 1.015 && pP4->M() < 1.025;
    bool xmass = xP4->M() > 5.15 && xP4->M() < 5.55;

    for (size_t j = 0; j < triggersToTest.size(); j++){

      int testingTrigger = triggersToTest[j];

      //if (tB.test(testingTrigger) && vProb > 0.0)

      if (tB.test(testingTrigger) && xyl/xylErr > 3.0 && cosA > 0.997 && jP4->Pt() > 7.0 && vProb > 0.0 && phimass && xmass)
      {
        tested = true;
        JPsi_vs_run_hists[j]->Fill(run);
        eventMap[ev] = 1;
        runMap[run] = 1;
        std::cout<<run<<std::endl;
      }
    }

      if(tested)
        JPsi_vs_run->Fill(run);

  }

  for (size_t j = 0; j < triggersToTest.size(); j++)
    JPsi_vs_run_hists[j]->Write();

  JPsi_vs_run->Write();

  return 0;

}

int skimMMKKTree(std::string path, std::string filename, std::string treename)
{
  UInt_t colors[13] = {1,2,3,6,7,8,30,40,46,38,29,34,9};

  TFile *oldfile = TFile::Open((path+"/"+filename).data());
  TTree *oldtree = (TTree*)oldfile->Get((treename).data());


  Long64_t nentries = oldtree->GetEntries();
  Double_t xM = 0.0, xDeltaM = 0.0;

  Double_t cosA  = 0.0;

  Double_t phiM  = 0.0;
  Double_t jPsiM  = 0.0;

  Double_t ctau  = 0.0;
  Double_t ctauErr  = 0.0;

  Double_t vProb  = 0.0;
  Double_t deltaR = 0.0;

  Int_t phiMType = 0, phiPType = 0;
  Int_t run = 0;
  Int_t phi_trigger = 0, jpsi_trigger = 0, trigger = 0;

  TLorentzVector *xP4 = 0, *jP4 = 0, *pP4 = 0;
  TLorentzVector *muonp_p4 = 0, *muonn_p4 = 0, *kaonp_p4 = 0, *kaonn_p4 = 0;

  //Output Variables
  Int_t t = 0, evt = 0, r = 0;
  Double_t jPt, xPt, mNPt, mPPt, kNPt, kPPt, cosAlpha, vP, pPt;
  Double_t xyl   = 0.0;
  Double_t xylErr   = 0.0;

  oldtree->SetBranchAddress("oniat_vProb",&vProb);
  oldtree->SetBranchAddress("trigger",&trigger);
  oldtree->SetBranchAddress("run",&run);
  // oldtree->SetBranchAddress("l_xy",&xyl);
  // oldtree->SetBranchAddress("lErr_xy",&xylErr);
  oldtree->SetBranchAddress("oniat_cosAlpha",&cosA);
  // oldtree->SetBranchAddress("phi_M",&phiM);
  // oldtree->SetBranchAddress("jpsi_M",&jPsiM);

  // oldtree->SetBranchAddress("phi_muonM_type",&phiMType);
  // oldtree->SetBranchAddress("phi_muonP_type",&phiPType);
  //
  oldtree->SetBranchAddress("oniat_ctauPV",&ctau);
  oldtree->SetBranchAddress("oniat_ctauErrPV",&ctauErr);
  //
  // oldtree->SetBranchAddress("phi_trigger",&phi_trigger);
  // oldtree->SetBranchAddress("jpsi_trigger",&jpsi_trigger);

  oldtree->SetBranchAddress("oniat_p4",&xP4);

  //oldtree->SetBranchAddress("oniat_rf_p4",&xP4);


  oldtree->SetBranchAddress("phi_p4",&pP4);
  oldtree->SetBranchAddress("kaonn_p4",&kaonn_p4);
  oldtree->SetBranchAddress("kaonp_p4",&kaonp_p4);

  oldtree->SetBranchAddress("psi_p4",&jP4);
  oldtree->SetBranchAddress("muonp_p4",&muonp_p4);
  oldtree->SetBranchAddress("muonn_p4",&muonn_p4);

  oldtree->SetBranchAddress("oniat_p4",&xP4);

  //oldtree->SetBranchAddress("oniat_rf_p4",&xP4);

  TFile *newfile = new TFile((treename + "_skimMMKK_" + filename).data(),"RECREATE");

  TTree *newtree = new TTree("skimmed_x_tree","skimmed_x_tree");

  newtree->Branch("oniat_vProb",&vProb);
  newtree->Branch("trigger",&trigger);
  newtree->Branch("run",&run);

  newtree->Branch("oniat_p4",&xP4);

  newtree->Branch("phi_p4",&pP4);
  newtree->Branch("kaonn_p4",&kaonn_p4);
  newtree->Branch("kaonp_p4",&kaonp_p4);

  newtree->Branch("psi_p4",&jP4);
  newtree->Branch("muonp_p4",&muonp_p4);
  newtree->Branch("muonn_p4",&muonn_p4);

  newtree->Branch("evt",&evt,"evt/I");
  newtree->Branch("run",&r,"run/I");
  newtree->Branch("trigger",&t,"trigger/I");
  newtree->Branch("xM",&xM,"xM/D");
  newtree->Branch("xDeltaM",&xDeltaM,"xDeltaM/D");
  newtree->Branch("phiM",&phiM,"phiM/D");
  newtree->Branch("jPsiM",&jPsiM,"jPsiM/D");
  newtree->Branch("deltaR",&deltaR,"deltaR/D");

  newtree->Branch("jPt",&jPt,"jPt/D");
  newtree->Branch("xPt",&xPt,"xPt/D");
  newtree->Branch("pPt",&pPt,"pPt/D");
  newtree->Branch("mNPt",&mNPt,"mNPt/D");
  newtree->Branch("mPPt",&mPPt,"mPPt/D");
  newtree->Branch("kNPt",&kNPt,"kNPt/D");
  newtree->Branch("kPPt",&kPPt,"kPPt/D");
  newtree->Branch("vProb",&vP,"vProb/D");
  newtree->Branch("cosAlpha",&cosAlpha,"cosAlpha/D");
  newtree->Branch("xyl",&xyl,"xyl/D");
  newtree->Branch("xylErr",&xylErr,"xylErr/D");
  newtree->Branch("xylErr",&xylErr,"xylErr/D");


  for (Long64_t i=0;i<nentries; i++) {
    oldtree->GetEntry(i);
    std::bitset<16> tB(trigger);
    // std::bitset<16> pM(phiKType);
    // std::bitset<16> pP(phiPType);

    float deltaEta = jP4->Eta() - pP4->Eta();
    float deltaPhi = std::fabs(jP4->Phi() - pP4->Phi());
    if (deltaPhi>pi) deltaPhi-=pi;

    deltaR = sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);

    phiM = pP4->M();
    jPsiM = jP4->M();
    xM = xP4->M();
    xDeltaM = xP4->M() - pP4->M() + pdg_Phi_mass;

    // evt = event;
    r = run;
    t = trigger;

    jPt  = jP4->Pt();
    xPt  = xP4->Pt();
    xPt  = pP4->Pt();
    mNPt = muonn_p4->Pt();
    mPPt = muonp_p4->Pt();
    kNPt = kaonn_p4->Pt();
    kPPt = kaonp_p4->Pt();
    vP       = vProb;
    cosAlpha = cosA;
    xyl      = ctau;
    xylErr   = ctauErr;

    if(ctau/ctauErr > 3.0 && vProb > 0.1)
      newtree->Fill();


  }


  return 0;

}

int drawMMKKTree(std::string path, std::string filename, std::string treename)
{

  UInt_t colors[13] = {1,2,3,6,7,8,30,40,46,38,29,34,9};

  TFile *oldfile = TFile::Open((path+"/"+filename).data());
  TTree *oldtree = (TTree*)oldfile->Get((treename).data());

  Long64_t nentries = oldtree->GetEntries();
  Double_t xM = 0.0, xDeltaM = 0.0;
  Double_t xyl   = 0.0;
  Double_t xylErr   = 0.0;
  Double_t cosA  = 0.0;

  Double_t phiM  = 0.0;
  Double_t jPsiM  = 0.0;

  Double_t ctau  = 0.0;
  Double_t ctauErr  = 0.0;

  Double_t vProb  = 0.0;

  Int_t phiMType = 0, phiPType = 0;
  Int_t run = 0;
  Int_t phi_trigger = 0, jpsi_trigger = 0, trigger = 0;

  TLorentzVector *xP4 = 0, *jP4 = 0, *pP4 = 0;
  TLorentzVector *muonp_p4 = 0, *muonn_p4 = 0, *kaonp_p4 = 0, *kaonn_p4 = 0;


  oldtree->SetBranchAddress("oniat_vProb",&vProb);
  oldtree->SetBranchAddress("trigger",&trigger);
  oldtree->SetBranchAddress("run",&run);
  // oldtree->SetBranchAddress("l_xy",&xyl);
  // oldtree->SetBranchAddress("lErr_xy",&xylErr);
  oldtree->SetBranchAddress("oniat_cosAlpha",&cosA);
  // oldtree->SetBranchAddress("phi_M",&phiM);
  // oldtree->SetBranchAddress("jpsi_M",&jPsiM);

  // oldtree->SetBranchAddress("phi_muonM_type",&phiMType);
  // oldtree->SetBranchAddress("phi_muonP_type",&phiPType);
  //
  oldtree->SetBranchAddress("oniat_ctauPV",&ctau);
  oldtree->SetBranchAddress("oniat_ctauErrPV",&ctauErr);
  //
  // oldtree->SetBranchAddress("phi_trigger",&phi_trigger);
  // oldtree->SetBranchAddress("jpsi_trigger",&jpsi_trigger);

  //oldtree->SetBranchAddress("oniat_p4",&xP4);

  oldtree->SetBranchAddress("oniat_rf_p4",&xP4);


  oldtree->SetBranchAddress("phi_p4",&pP4);
  oldtree->SetBranchAddress("kaonn_p4",&kaonn_p4);
  oldtree->SetBranchAddress("kaonp_p4",&kaonp_p4);

  oldtree->SetBranchAddress("psi_p4",&jP4);
  oldtree->SetBranchAddress("muonp_p4",&muonp_p4);
  oldtree->SetBranchAddress("muonn_p4",&muonn_p4);

  // oldtree->SetBranchAddress("muonM_jpsi_p4",&mM_phi_P4);
  // oldtree->SetBranchAddress("muonP_jpsi_p4",&mP_phi_P4);
  //Create a new file + a clone of old tree in new file
  TCanvas c("c","c",1200,1600);

  TFile *newfile = new TFile(("DrawSkim_"+ filename).data(),"RECREATE");
  // for(int j = 0; j < 1; j++)
  // {


  Double_t xmin = 4.0, xmax = 6.0;
  Int_t xBin = ((xmax - xmin)/0.001);

  TTree *newtree = oldtree->CloneTree(0);
  // TH1F* phi_triggrHist = new TH1F("phi_triggrHist","phi_triggrHist",600,0.6,1.2);
  TH1F* phiHist = new TH1F("phiHist","phiHist",1000,0.9,1.15);
  TH1F* jpsiHist = new TH1F("jpsiHist","jpsiHist",1000,2.9,3.3);
  TH1F* xHist = new TH1F("xHist","xHist",xBin,xmin,xmax);
  TH1F* xHistDeltaM = new TH1F("xHistDeltaM","xHistDeltaM",xBin,xmin,xmax);

  TH1F* x_ptHist = new TH1F("x_ptHist","x_ptHist",1000,0.0,100.0);
  TH1F* jpsi_ptHist = new TH1F("jpsi_ptHist","jpsi_ptHist",1000,0.0,100.0);
  TH1F* phi_ptHist = new TH1F("phi_ptHist","phi_ptHist",1000,0.0,100.0);

  TH1F* jpsiMP_ptHist = new TH1F("jpsiMP_ptHist","jpsiMP_ptHist",1000,0.0,100.0);
  TH1F* jpsiMM_ptHist = new TH1F("jpsiMM_ptHist","jpsiMM_ptHist",1000,0.0,100.0);
  TH1F* jpsiMHig_ptHist = new TH1F("jpsiMHig_ptHist","jpsiMHig_ptHist",1000,0.0,100.0);
  TH1F* jpsiMLow_ptHist = new TH1F("jpsiMLow_ptHist","jpsiMLow_ptHist",1000,0.0,100.0);

  TH1F* phiKP_ptHist = new TH1F("phiKP_ptHist","phiKP_ptHist",1000,0.0,100.0);
  TH1F* phiKM_ptHist = new TH1F("phiKM_ptHist","phiKM_ptHist",1000,0.0,100.0);
  TH1F* phiKHig_ptHist = new TH1F("phiKHig_ptHist","phiKHig_ptHist",1000,0.0,100.0);
  TH1F* phiKLow_ptHist = new TH1F("phiKLow_ptHist","phiKLow_ptHist",1000,0.0,100.0);

  TH1F* dRJpsiPhi = new TH1F("dRJpsiPhi","dRJpsiPhi",1000,-10.0,10.0);

  TH2F* mpts = new TH2F("mpts","mpts",1000,0.0,20.0,1000,0.0,100.0);
  TH2F* kpts = new TH2F("kpts","kpts",1000,0.0,20.0,1000,0.0,100.0);

  std::vector<TH1F*> phiHists;
  std::vector<TH1F*> jpsiHists;
  std::vector<TH1F*> xHists;
  std::vector<TH1F*> ptJHists;
  std::vector<TH1F*> xHistsDeltaM;

  std::vector<TH2F*> psiMuonsPts;
  std::vector<TH2F*> phiKaonsPts;

  for (int i = 0; i < noHlts; i++)
  {
    phiHists.push_back(new TH1F((hltsName[i] + "_phi").data(),(hltsName[i] + "_phi").data(),1000,0.9,1.15));
    jpsiHists.push_back(new TH1F((hltsName[i] + "_jpsi").data(),(hltsName[i] + "_jpsi").data(),1000,2.9,3.3));
    xHists.push_back(new TH1F((hltsName[i] + "_x").data(),(hltsName[i] + "_x").data(),xBin,xmin,xmax));
    ptJHists.push_back(new TH1F((hltsName[i] + "_jpsi_pt").data(),(hltsName[i] + "_jpsi_pt").data(),1000,0.0,100.0));
    xHistsDeltaM.push_back(new TH1F((hltsName[i] + "_x_deltam").data(),(hltsName[i] + "_x_deltam").data(),xBin,xmin,xmax));
    psiMuonsPts.push_back(new TH2F((hltsName[i] + "_mpts").data(),(hltsName[i] + "_mpts").data(),1000,0.0,100.0,1000,0.0,100.0));
    phiKaonsPts.push_back(new TH2F((hltsName[i] + "_kpts").data(),(hltsName[i] + "_kpts").data(),1000,0.0,100.0,1000,0.0,100.0));
  }

  // std::string hltsName[13] = {
  //   "HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi", //0
  // "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi", //1
  // "HLT_Mu20_TkMu0_Phi", //2
  // "HLT_Dimuon14_Phi_Barrel_Seagulls", //3
  // "HLT_Mu25_TkMu0_Phi", //4
  // "HLT_Dimuon24_Phi_noCorrL1", //5
  // "HLT_DoubleMu4_JpsiTrkTrk_Displaced", //6
  // "HLT_DoubleMu4_JpsiTrk_Displaced", //7
  // "HLT_DoubleMu4_Jpsi_Displaced", //8
  // "HLT_DoubleMu4_3_Jpsi_Displaced", //9
  // "HLT_Dimuon20_Jpsi_Barrel_Seagulls", //10
  // "HLT_Dimuon25_Jpsi", //11
  // "HLT_Dimuon0_Jpsi" //12
  // };

  std::vector <int> triggersToTest;

  //OUR TRIGGERs
  // triggersToTest.push_back(0);
  triggersToTest.push_back(1);

  //PHI INCLUSIVE TRIGGERs
  // triggersToTest.push_back(2);
  // triggersToTest.push_back(3);
  // triggersToTest.push_back(4);
  // triggersToTest.push_back(5);

  //JPSI INCLUSIVE TRIGGERs
  // triggersToTest.push_back(10);
  // triggersToTest.push_back(11);
  // triggersToTest.push_back(12);

  //DISPLACED JPSI TRIGGERs
  triggersToTest.push_back(6); //dis
  triggersToTest.push_back(7); //dis
  // triggersToTest.push_back(8); //dis
  // triggersToTest.push_back(9); //dis

  for (Long64_t i=0;i<nentries; i++) {
    oldtree->GetEntry(i);
    std::bitset<16> tB(trigger);
    // std::bitset<16> pM(phiKType);
    // std::bitset<16> pP(phiPType);

    phiM = pP4->M();
    jPsiM = jP4->M();
    xM = xP4->M();
    xDeltaM = xP4->M() - pP4->M() + pdg_Phi_mass;

    float deltaEta = jP4->Eta() - pP4->Eta();
    float deltaPhi = std::fabs(jP4->Phi() - pP4->Phi());
    if (deltaPhi>pi) deltaPhi-=pi;
    float deltaR = sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);

    // jpsiMP_ptHist->Fill(mP_jpsi_P4->Pt());
    // jpsiMM_ptHist->Fill(mM_jpsi_P4->Pt());
    //
    // jpsiMHig_ptHist->Fill(std::max(mP_jpsi_P4->Pt(),mM_jpsi_P4->Pt()));
    // jpsiMLow_ptHist->Fill(-std::max(-mP_jpsi_P4->Pt(),-mM_jpsi_P4->Pt()));
    //
    // phiKHig_ptHist->Fill(std::max(mP_phi_P4->Pt(),mM_phi_P4->Pt()));
    // phiKLow_ptHist->Fill(-std::max(-mP_phi_P4->Pt(),-mM_phi_P4->Pt()));
    //
    // phiKP_ptHist->Fill(mP_phi_P4->Pt());
    // phiKM_ptHist->Fill(mM_phi_P4->Pt());

    bool test = false;
    bool jpsimass = jPsiM < 3.19 && jPsiM > 3.0;
    bool phimass = phiM > 0.95 && phiM < 1.3;
    for (size_t j = 0; j < triggersToTest.size(); j++){

      int testingTrigger = triggersToTest[j];
      // if (tB.test(j) && cosA > 0.995 && vProb > 0.01 && xyl/xylErr > 2.0 && trigger > 0)
      // if (tB.test(j))
      //if (xyl/xylErr > 0.0 && xP4->Pt() > 6.0 && cosA > 0.997 && pP4->Pt() > 8.0 && jP4->Pt() > 5.0 && mP_phi_P4->Pt() > 2.0 && mM_phi_P4->Pt() > 2.0 && tB.test(j) && vProb > 0.05 && deltaR < 1.0 && deltaR > 0.0 && jpsimass && phiKass)
// if (tB.test(testingTrigger) && run > 305388 && vProb > 0.5 && cosA > 0.997 && deltaR < 0.8  && jpsimass && phimass && ctau/ctauErr > 3.0)
      // if (tB.test(testingTrigger) && run > 305388 && ctau/ctauErr > 3.0 && phimass && std::max(kaonn_p4->Pt(),kaonp_p4->Pt())>1.2 && -std::max(-kaonn_p4->Pt(),-kaonp_p4->Pt())>1.0 && jP4->Pt() > 4.0)
      // if (tB.test(testingTrigger) && run > 305388 && ctau/ctauErr > 3.0 && phimass && kaonn_p4->Pt() >1.0 && kaonp_p4->Pt()>1.0)
      // if (tB.test(testingTrigger) && ctau/ctauErr < 2.0 && phimass && vProb > 0.2 && cosA > 0.997 && deltaR < 2.0  && jpsimass)
      if (tB.test(testingTrigger))
      {

        test = true;
        phiHists[testingTrigger]->Fill(phiM);
        jpsiHists[testingTrigger]->Fill(jPsiM);
        //  if(phiKass && jpsimass && cosA > 0.995 && vProb>0.02 && xyl/xylErr > 2.0)
        ptJHists[testingTrigger]->Fill(jP4->Pt());
        xHists[testingTrigger]->Fill(xM);
        xHistsDeltaM[testingTrigger]->Fill(xDeltaM);
        psiMuonsPts[testingTrigger]->Fill(-std::max(-muonp_p4->Pt(),-muonn_p4->Pt()),std::max(muonp_p4->Pt(),muonn_p4->Pt()));
        phiKaonsPts[testingTrigger]->Fill(-std::max(-kaonn_p4->Pt(),-kaonp_p4->Pt()),std::max(kaonn_p4->Pt(),kaonp_p4->Pt()));
      }
    }
    // if(cosA > 0.995 && phiKass && jpsimass && vProb>0.02 && xyl/xylErr > 2.0)
    if (test){
      xHist->Fill(xM);
      xHistDeltaM->Fill(xDeltaM);
      phiHist->Fill(phiM);
      jpsiHist->Fill(jPsiM);
      phi_ptHist->Fill(pP4->Pt());
      dRJpsiPhi->Fill(deltaR); //cut < 1
      x_ptHist->Fill(xP4->Pt());
      jpsi_ptHist->Fill(jP4->Pt());

      jpsiMP_ptHist->Fill(muonp_p4->Pt());
      jpsiMM_ptHist->Fill(muonn_p4->Pt());

      jpsiMHig_ptHist->Fill(std::max(muonp_p4->Pt(),muonn_p4->Pt()));
      jpsiMLow_ptHist->Fill(-std::max(-muonp_p4->Pt(),-muonn_p4->Pt()));

      phiKHig_ptHist->Fill(std::max(kaonn_p4->Pt(),kaonp_p4->Pt()));
      phiKLow_ptHist->Fill(-std::max(-kaonn_p4->Pt(),-kaonp_p4->Pt()));

      phiKP_ptHist->Fill(kaonp_p4->Pt());
      phiKM_ptHist->Fill(kaonn_p4->Pt());

      kpts->Fill(-std::max(-muonp_p4->Pt(),-muonn_p4->Pt()),std::max(muonp_p4->Pt(),muonn_p4->Pt()));
      mpts->Fill(-std::max(-kaonp_p4->Pt(),-kaonn_p4->Pt()), std::max(kaonp_p4->Pt(),kaonn_p4->Pt()));


    }
  }

  xHistDeltaM->Write();
  dRJpsiPhi->Write();
  x_ptHist->Write();
  jpsi_ptHist->Write();
  phi_ptHist->Write();
  jpsiMP_ptHist->Write();
  jpsiMM_ptHist->Write();
  jpsiMHig_ptHist->Write();
  jpsiMLow_ptHist->Write();
  phiKHig_ptHist->Write();
  phiKLow_ptHist->Write();
  phiKP_ptHist->Write();
  phiKM_ptHist->Write();
  kpts->Write();
  mpts->Write();
  //newtree->Draw("phi_M","","same");
  // }
  phiHist->SetMinimum(1.0);
  phiHist->SetMaximum(phiHist->GetMaximum()*5.0);
  //oldtree->Draw("phi_M");
  // TH1F* phi_triggrHist = (TH1F*)gDirectory->Get("phi_triggrHist");

  phiHist->SetLineColor(kBlue);
  phiHist->Write();
  phiHist->Draw();

  TLegend* leg = new TLegend(0.1,0.5,0.45,0.9);
  leg->AddEntry(phiHist,(phiHist->GetName()),"l");
  for (int i = 0; i < 13; i++)
  {
    phiHists[i]->SetLineColor(colors[i]);
    phiHists[i]->SetLineWidth(2);
    if(i>5) phiHists[i]->SetLineStyle(kDashed);
    phiHists[i]->Draw("same");
    leg->AddEntry(phiHists[i],(phiHists[i]->GetName()),"l");
    phiHists[i]->Write();
    ptJHists[i]->Write();
    xHistsDeltaM[i]->Write();
    psiMuonsPts[i]->Write();
    phiKaonsPts[i]->Write();
  }

  leg->Draw();
  c.SetLogy(1);
  c.SaveAs("phitriggerCheck.png");
  c.SaveAs("phitriggerCheck.eps");
  c.SaveAs("phitriggerCheck.root");

  jpsiHist->SetMinimum(1.0);
  jpsiHist->SetMaximum(jpsiHist->GetMaximum()*5.0);
  //oldtree->Draw("phi_M");
  // TH1F* phi_triggrHist = (TH1F*)gDirectory->Get("phi_triggrHist");

  jpsiHist->SetLineColor(kBlue);
  jpsiHist->Write();
  jpsiHist->Draw();

  leg = new TLegend(0.1,0.5,0.4,0.9);
  leg->AddEntry(jpsiHist,(jpsiHist->GetName()),"l");
  for (int i = 0; i < 13; i++)
  {
    jpsiHists[i]->SetLineColor(colors[i]);
    jpsiHists[i]->SetLineWidth(2);
    if(i>5) jpsiHists[i]->SetLineStyle(kDashed);
    jpsiHists[i]->Draw("same");
    leg->AddEntry(jpsiHists[i],(jpsiHists[i]->GetName()),"l");
    jpsiHists[i]->Write();
  }


  // phi_triggrHist->Draw("same");
  leg->Draw();
  c.SetLogy(1);
  c.SaveAs("jpsitriggerCheck.png");
  c.SaveAs("jpsitriggerCheck.eps");
  c.SaveAs("jpsitriggerCheck.root");

  xHist->SetMinimum(0.1);
  xHist->SetMaximum(xHist->GetMaximum());//*5.0);

  TLine line (5.35,0.1,5.35,xHist->GetMaximum());
  line.SetLineColor(kRed);
  line.SetLineWidth(2);

  xHist->SetLineColor(kBlue);
  xHist->Write();
  xHist->Draw();

  leg = new TLegend(0.1,0.55,0.35,0.9);
  leg->AddEntry(xHist,(xHist->GetName()),"l");
  for (int i = 0; i < 13; i++)
  {
    xHists[i]->SetLineColor(colors[i]);
    xHists[i]->SetLineWidth(2);
    if(i>5) xHists[i]->SetLineStyle(kDashed);
    xHists[i]->Draw("same");
    leg->AddEntry(xHists[i],(xHists[i]->GetName()),"l");
    xHists[i]->Write();
  }

  // line.Draw();
  // phi_triggrHist->Draw("same");
  leg->Draw();
  c.SetLogy(0);
  c.SaveAs("xtriggerCheck.png");
  c.SaveAs("xtriggerCheck.eps");
  c.SaveAs("xtriggerCheck.root");

  return 0;

}

int selectXTreeHLT()
{

  TFile *oldfile = TFile::Open("/Users/adrianodiflorio/Documents/Git/X4140/iPythons/skimmed_vProb.root");

  TTree *oldtree = (TTree*)oldfile->Get("xTree");
  Long64_t nentries = oldtree->GetEntries();
  Double_t vProb   = 0.0;
  UInt_t trigger = 0;
  oldtree->SetBranchAddress("vProb",&vProb);
  oldtree->SetBranchAddress("trigger",&trigger);
  //Create a new file + a clone of old tree in new file
  TFile *newfile = new TFile("skimmed.root","RECREATE");
  TTree *newtree = oldtree->CloneTree(0);

  for (Long64_t i=0;i<nentries; i++) {
    oldtree->GetEntry(i);
    // std::bitset<16> binary(trigger);
    // std::cout << trigger << "-" << binary << std::endl;
    if(trigger>0)
    newtree->Fill();
  }

  newtree->Print();
  newtree->Write();

  return 0;

}

/*
# coding: utf-8

# In[1]:


import ROOT
from ROOT import TFile,TH1,TH1F,TCanvas,TNtuple,TTreeReader,TTreeReaderValue
from ROOT import RooFit
from ROOT.RooFit import Layout
from ROOT import RooStats
from ROOT import RooAbsData
RooAbsData.setDefaultStorageType ( RooAbsData.Tree )
from array import array
import sys


# In[3]:


from ROOT import RooRealVar,RooAbsPdf,RooChebychev,RooExponential,RooGaussian,RooAbsPdf,RooPlot,RooAddPdf,RooDataHist,RooArgSet,RooArgList
from ROOT import kGreen,kRed,kBlack,kBlue,kDashed,kDotted,kMagenta,RooVoigtian
from ROOT.RooFit import Components,LineColor,LineStyle,Name,Normalization,Range,SelectVars
from ROOT import RooDataSet,RooFormulaVar,RooLinkedList


# In[7]:


no_hlts = 13


# In[8]:


#rootfile = "/Users/adrianodiflorio/Desktop/mmkk2017/09Jan2017.root"
#rootfile = "/Users/adrianodiflorio/Desktop/mmkk2017/allphi_DataF.root"
rootfile = "/Users/adrianodiflorio/Desktop/mmkk2017/phiJpsiTriggersBCDEF.root"
inputfile = TFile(rootfile,"READ")
inputfile.ls()
xTupleDir = (inputfile.Get("rootuple"))
xTupleDir.ls()
#pTuple = (xTupleDir.Get("pTree"))
xTuple = (xTupleDir.Get("xTree"))
#jTuple = (xTupleDir.Get("jTree"))

event = 2

xTuple.SetBranchAddress("event",event);
newfile = TFile("small.root","recreate");
newtree = xTuple.CloneTree();
newtree.CopyEntries(xTuple);

newtree.Write()

sys.exit()

# In[16]:

file = TFile("newFile.root","RECREATE")
canvas = TCanvas("canvas","canvas",1200,1000)
mass = RooRealVar("xM","M(#mu#mu#mu#mu)[GeV]",5.15,5.55)
trigger = RooRealVar("trigger","trigger",0.0,10000)
vProb = RooRealVar("vProb","vProb",-1.0,1.0)
alldata = RooDataSet("alldata","alldata",xTuple,RooArgSet(mass), RooFormulaVar("vProb","vProb","vProb>0.01",RooArgList(vProb)))#,cutFormula)
frame = mass.frame(Range(5.15,5.55))
alldata.plotOn(frame,RooLinkedList())
alldata.Write()
frame.Draw()


# In[ ]:


canvas.SaveAs("testCanvas.eps")
*/
