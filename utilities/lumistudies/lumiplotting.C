#include <TH2F.h>
#include <TH1F.h>
#include <TFile.h>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <TCanvas.h>
#include "TStyle.h"
#include <TLine.h>
#include <TLegend.h>

// B dataset: 193833-196531
// C dataset: 198022-203742
// D dataset: 203777-208686

void lumiplotting(TFile* f)
{

	Int_t colors[13] = {1,2,3,6,7,8,30,40,46,38,29,34,9};

	gStyle->SetOptStat(000000000);

	int bdown = 193833, bup = 196531, cdown = 198022, cup = 203742, ddown =203777, dup = 208686;

	double yMax = 0.04, yMin = 0.0;

	// TLine blinedown(bdown-2,yMin,bdown-2,yMax);
	// TLine blineup(bup-2,yMin,bup-2,yMax);
  //
	// TLine clinedown(cdown-2,yMin,cdown-2,yMax);
	// TLine clineup(cup-2,yMin,cup-2,yMax);
  //
	// TLine dlinedown(ddown-2,yMin,ddown-2,yMax);
	// TLine dlineup(dup-2,yMin,dup-2,yMax);
  //
	// blinedown.SetLineColor(kRed);
	// blinedown.SetLineWidth(1);
  //
	// blineup.SetLineColor(kRed);
	// blineup.SetLineWidth(1);
  //
	// clinedown.SetLineColor(kRed);
	// clinedown.SetLineWidth(1);
  //
	// clineup.SetLineColor(kRed);
	// clineup.SetLineWidth(1);
  //
	// dlinedown.SetLineColor(kRed);
	// dlinedown.SetLineWidth(1);
  //
	// dlineup.SetLineColor(kRed);
	// dlineup.SetLineWidth(1);

	std::vector< std::string > datasetNames {"B","C","D","E","F"};

	std::vector< TH1F* > hltHists;
	std::vector< TH1F* > hltHistsAllRange;
	for (size_t i = 0; i < 15; i++) {
		TH1F* hltHist = (TH1F*)f->Get(("lumiJPsirate" + std::to_string(i)).data());
		if(hltHist)
		{
			std::cout<< "Found hlt hist for HLT no. " << i << std::endl;
			hltHists.push_back(hltHist);
		}

	}

	hltHistsAllRange.push_back(new TH1F("allruns_hlt_all", "JPsis from B^{0}_{s} candidate; Run Number;no. of J/Psi per lumi (mb)",40000, 280000, 320000));

	for (size_t i = 0; i < hltHists.size(); i++)
	{
		// for (size_t j = 0; j < datasetNames.size(); j++) {
			hltHistsAllRange.push_back(new TH1F(("allruns_hlt_" + std::to_string(i)).data(), "JPsis from B^{0}_{s} candidate; Run Number;no. of J/Psi per lumi (mb)",40000, 280000, 320000));
		// }
	}

	TCanvas canvas("canvas","canvas",1200,800);

	TLegend legend(0.6,0.8,0.9,0.9);

	for(int i = 0; i<=hltHistsAllRange[0]->GetNbinsX();++i)
		hltHistsAllRange[0]->SetBinContent(i,0.0);

	for (size_t j = 1; j < hltHistsAllRange.size(); j++)
	{
		for(int i = 0; i<hltHistsAllRange[j]->GetNbinsX();++i)
		{
			double center = hltHistsAllRange[j]->GetBinCenter(i);

			hltHistsAllRange[j]->SetBinContent(hltHists[j-1]->FindBin(center),hltHists[j-1]->GetBinContent(i));
			hltHistsAllRange[j]->SetBinError(hltHists[j-1]->FindBin(center),hltHists[j-1]->GetBinError(i));
		}
		hltHistsAllRange[0]->Add(hltHistsAllRange[j]);
	}

	for (size_t i = 0; i < hltHistsAllRange.size(); i++)
	{
		// for (size_t j = 0; j < datasetNames.size(); j++) {
			hltHistsAllRange[i]->SetLineColor(1);
			hltHistsAllRange[i]->SetLineStyle(0);
			hltHistsAllRange[i]->SetLineWidth(1);
			hltHistsAllRange[i]->SetMarkerStyle(20);
			hltHistsAllRange[i]->SetMarkerColor(colors[i]);
			hltHistsAllRange[i]->SetMaximum(yMax);
			hltHistsAllRange[i]->SetMinimum(yMin);
		// }
	}

	hltHistsAllRange[0]->Draw("PE");

	for (size_t j = 1; j < hltHistsAllRange.size(); j++)
	{
		hltHistsAllRange[j]->Draw("PESame");
	}

		// b_runs_hlt4->Draw("PE");
		// b_runs_hlt8->Draw("PEsame");

		// legend.AddEntry("b_runs_hlt4","HLT_4 + pt cut","lep");
		// legend.AddEntry("b_runs_hlt8","HLT_8 + Lxy cut","lep");
		// legend.Draw();

	canvas.SaveAs("jspiTrendRun.root");

	for (size_t j = 1; j < hltHistsAllRange.size(); j++)
	{
		hltHistsAllRange[j]->Draw();
		canvas.SaveAs(("jspiTrendRun_" + std::to_string(j) + ".root").data());
	}

		// c_runs_hlt4->Draw("PE");
		// c_runs_hlt8->Draw("PEsame");
		// legend.Draw();
		// canvas.SaveAs("cRun.eps");
    //
		// d_runs_hlt4->Draw("PE");
		// d_runs_hlt8->Draw("PEsame");
		// legend.Draw();
		// canvas.SaveAs("dRuns.eps");
    //
		// all_runs_hlt4->Draw("PE");
		// //all_runs_hlt4->SetOptStat(0000);
		// all_runs_hlt8->Draw("PEsame");
    //
		// blinedown.Draw();
		// blineup.Draw();
    //
		// dlinedown.Draw();
		// dlineup.Draw();
    //
		// clinedown.Draw();
		// legend.Draw();
		// // clineup.Draw();
    //
		// canvas.SaveAs("allRuns.eps");

}
