#include         <iostream>
#include         <fstream>
#include         <TGaxis.h>
#include         <sstream>
#include         <TGraph.h>
#include         <TGraphErrors.h>
#include         <TCanvas.h>
#include         <TApplication.h>
#include         <stdlib.h>
#include         <TF1.h>
#include         <TH1.h>
#include         <TProfile.h>
#include         <TObjArray.h>
#include         <TStyle.h>
#include         <TMarker.h>
#include         <math.h>
#include         <TStyle.h>
#include         <TPaveStats.h>
#include         <TPaveText.h>
#include         <vector>
#include         <string.h>
#include         <fstream>
#include         <TROOT.h>
#include         <TFile.h>
#include         <TLegend.h>
#include         <TLegendEntry.h>
#include         <time.h>
#include         <TH2F.h>
#include         <assert.h>
#include         <string>
#include         <TRandom.h>
#include         <TTree.h>
#include         <TChain.h>
#include         <TVector.h>
#include         <vector>
#include         <utility>
#include         <TLeaf.h>
#include         <math.h>
#include	 <TCut.h>
#include	 <TMath.h>
#include	 <TRandom3.h>

using            namespace std;

// some forward declarations of useful functions
double SetPoissonErrors(int counts);
void FillInAcceptancesFromFile(TString fileName, double totalCounts, int flag);
void FillInCountsFromFile(TString fileName, TH1D* hToFill);
TGraph* CreateAcceptancesGraph(double timeLower, double timeUpper);

// some global access vectors to read in acceptance files
vector <double> t094;
vector <double> t144;
vector <double> t244;
vector <double> t344;
vector <double> t444;
vector <double> t544;
vector <double> t644;
vector <double> n094;
vector <double> n144;
vector <double> n244;
vector <double> n344;
vector <double> n444;
vector <double> n544;
vector <double> n644;


struct DataEvent
{
  double TDCE;
  double TDCW;
  int PID;
  int Type;
  int Side;
  double Erecon;
  double Erecon_ee;
  int badTimeFlag;

};

struct SimEvent
{
  double Edep_EdepE;
  double Edep_EdepW;
  double EdepQ_EdepQE;
  double EdepQ_EdepQW;
  double MWPCEnergy_MWPCEnergyE;
  double MWPCEnergy_MWPCEnergyW;
  double time_timeE;
  double time_timeW;

};

// Used for visualization, keeps the graph on screen.
TApplication plot_program("FADC_readin",0,0,0,0);

//-------------------------------------------------//
//------------ Start of Program -------------------//
//-------------------------------------------------//

int main(int argc, char* argv[])
{
  gRandom->SetSeed(0);
  gROOT->SetStyle("Pub");

  // set time window variables.
  double timeLowerEdgeE1 = -10;
  double timeUpperEdgeE1 = 140;
  double timeLowerEdgeW1 = -10;
  double timeUpperEdgeW1 = 3;
  double timeLowerEdgeE2 = -10;
  double timeUpperEdgeE2 = 5;
  double timeLowerEdgeW2 = -10;
  double timeUpperEdgeW2 = 140;

  double windowLower = 0;
  double windowUpper = 12;

  TChain *bgchain = new TChain("pass3");
  TChain *fgchain = new TChain("pass3");

  bgchain->Add("TimeCalibrated_BGRuns_type1_fixed_v3.root");
  fgchain->Add("TimeCalibrated_FGRuns_type1_fixed_v3.root");

  // define all the cuts we will use later.
  TCut basicCut = "(PID == 1 && badTimeFlag == 0)";
  TCut energyCut = "(Erecon_ee > 0 && Erecon_ee < 800)";
  TCut fiducialCut = "(((xE.center)*(xE.center) + (yE.center)*(yE.center) < 49*49) && ((xW.center)*(xW.center) + (yW.center)*(yW.center) < 49*49))";
  TCut simulationCut = "EdepQ.EdepQE > 0 && EdepQ.EdepQW > 0 && MWPCEnergy.MWPCEnergyE > 0 && MWPCEnergy.MWPCEnergyW > 0 && time.timeE < 200 && time.timeW < 200 && timeE < timeW";
  TCut inESTPCut = Form("newTDC2TimeE > -1 && newTDC2TimeE < 5 && newTDC2TimeW > %f", windowLower);
  TCut inWSTPCut = Form("newTDC2TimeW > -4 && newTDC2TimeW < 2 && newTDC2TimeE > %f", windowLower);
  TCut inESTPAndTimeWinCut = Form("newTDC2TimeE > -1 && newTDC2TimeE < 5 && newTDC2TimeW > %f && newTDC2TimeW < %f", windowLower, windowUpper);
  TCut inWSTPAndTimeWinCut = Form("newTDC2TimeW > -4 && newTDC2TimeW < 2 && newTDC2TimeE > %f && newTDC2TimeE < %f", windowLower, windowUpper);


  TCut time1STPCut = Form("newTDC2TimeE < %f && newTDC2TimeE > %f && newTDC2TimeW > %f && newTDC2TimeW < %f", timeUpperEdgeE1, timeLowerEdgeE1, timeLowerEdgeW1, timeUpperEdgeW1);
  TCut time2STPCut = Form("(newTDC2TimeE < %f && newTDC2TimeE > %f && newTDC2TimeW > %f && newTDC2TimeW < %f) || (newTDC2TimeE < %f && newTDC2TimeE > %f && newTDC2TimeW > %f && newTDC2TimeW < %f)",
			timeUpperEdgeE1, timeLowerEdgeE1, timeLowerEdgeW1, timeUpperEdgeW1,
			timeUpperEdgeE2, timeLowerEdgeE2, timeLowerEdgeW2, timeUpperEdgeW2);

  // first canvas, timing spectra with MPM simulation.
  // first canvas fills in our bg and fg spectra. It is a temp canvas effectively
  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->Divide(2,1);

  TH1D *hTDCEfg = new TH1D("hTDCEfg", "FG TDC Spectra", 300, -10, 140);
  hTDCEfg->GetXaxis()->SetTitle("Time (ns)");
  hTDCEfg->SetLineColor(2);

  TH1D *hTDCWfg = new TH1D("hTDCWfg", "FG TDC Spectra", 300, -10, 140);
  hTDCWfg->SetLineColor(4);

  TH1D *hTDCEbg = new TH1D("hTDCEbg", "BG TDC Spectra", 300, -10, 140);
  hTDCEbg->GetXaxis()->SetTitle("Time (ns)");
  hTDCEbg->SetLineColor(2);

  TH1D *hTDCWbg = new TH1D("hTDCWbg", "BG TDC Spectra", 300, -10, 140);
  hTDCWbg->SetLineColor(4);

  c1->cd(1);
  // i.e. if you are in the West STP, look at East events outside of East STP. For hTDCEbg.
  bgchain->Draw("newTDC2TimeE >> hTDCEbg", basicCut && fiducialCut && energyCut && inWSTPCut);
  bgchain->Draw("newTDC2TimeW >> hTDCWbg", basicCut && fiducialCut && energyCut && inESTPCut, "SAME");

  TLegend *l1a = new TLegend(0.1,0.8,0.4,0.9);
  l1a->AddEntry(hTDCEbg,"Calibrated TDCE","l");
  l1a->AddEntry(hTDCWbg,"Calibrated TDCW","l");
  l1a->Draw();
  gPad->SetLogy();

  c1->cd(2);
  fgchain->Draw("newTDC2TimeE >> hTDCEfg", basicCut && fiducialCut && energyCut && inWSTPCut);
  fgchain->Draw("newTDC2TimeW >> hTDCWfg", basicCut && fiducialCut && energyCut && inESTPCut, "SAME");

  TLegend *l1b = new TLegend(0.1,0.8,0.4,0.9);
  l1b->AddEntry(hTDCEfg,"Calibrated TDCE","l");
  l1b->AddEntry(hTDCWfg,"Calibrated TDCW","l");
  l1b->Draw();
  gPad->SetLogy();


  TCanvas *c2 = new TCanvas("c2","c2");
  c2->cd();

  TChain *mpm_sim_chain = new TChain("anaTree");
  mpm_sim_chain->Add("mpm_sim_betas_500EWFiles_2012-2013.root");

  TRandom3 *engine = new TRandom3(0);

  TH1D *hSim0 = new TH1D("hSim0sigma", "hSim", 300, -10, 140);
  hSim0->SetLineColor(6);
  TH1D *hSim2 = new TH1D("hSim2sigma", "hSim", 300, -10, 140);
  hSim2->SetLineColor(3);
  TH1D *hSim3 = new TH1D("hSim3sigma", "hSim", 300, -10, 140);
  hSim3->SetLineColor(2);

  SimEvent evt;
  mpm_sim_chain->GetBranch("Edep")->GetLeaf("EdepE")->SetAddress(&evt.Edep_EdepE);
  mpm_sim_chain->GetBranch("Edep")->GetLeaf("EdepW")->SetAddress(&evt.Edep_EdepW);
  mpm_sim_chain->GetBranch("EdepQ")->GetLeaf("EdepQE")->SetAddress(&evt.EdepQ_EdepQE);
  mpm_sim_chain->GetBranch("EdepQ")->GetLeaf("EdepQW")->SetAddress(&evt.EdepQ_EdepQW);
  mpm_sim_chain->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyE")->SetAddress(&evt.MWPCEnergy_MWPCEnergyE);
  mpm_sim_chain->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyW")->SetAddress(&evt.MWPCEnergy_MWPCEnergyW);
  mpm_sim_chain->GetBranch("time")->GetLeaf("timeE")->SetAddress(&evt.time_timeE);
  mpm_sim_chain->GetBranch("time")->GetLeaf("timeW")->SetAddress(&evt.time_timeW);

  for(unsigned int i = 0; i < 500000 /*mpm_sim_chain->GetEntries()*/; i++)
  {
    mpm_sim_chain->GetEntry(i);

    if(evt.EdepQ_EdepQE > 0 && evt.EdepQ_EdepQW > 0 && evt.MWPCEnergy_MWPCEnergyE > 0 && evt.MWPCEnergy_MWPCEnergyW > 0
        && evt.time_timeE < 200 && evt.time_timeW < 200 && evt.time_timeE < evt.time_timeW)
    {
      hSim0->Fill(evt.time_timeW - evt.time_timeE);
      hSim2->Fill(engine->Gaus(evt.time_timeW - evt.time_timeE, 2));
      hSim3->Fill(engine->Gaus(evt.time_timeW - evt.time_timeE, 3));
    }

    if(i % 100000 == 0)
    {
      cout << "Filling event " << i << "/" << mpm_sim_chain->GetEntries() << endl;
    }
  }

  TH1D *hTimeE_bgSub = new TH1D("BGSubE", "BG Subtracted Time East", 300, -10, 140);
  hTimeE_bgSub->SetLineColor(1);
  hTimeE_bgSub->Add(hTDCEfg, hTDCEbg, 1, -5.07);
  hTimeE_bgSub->SetStats(kFALSE);
  hTimeE_bgSub->GetXaxis()->SetTitle("#left| T_{East} - T_{West} #right| (ns)");
  hTimeE_bgSub->GetYaxis()->SetTitle("Counts");
  hTimeE_bgSub->SetTitle("");
  hTimeE_bgSub->Draw();

  TH1D *hTimeW_bgSub = new TH1D("BGSubW", "BG Subtracted Time West", 300, -10, 140);
  hTimeW_bgSub->SetLineColor(4);
  hTimeW_bgSub->Add(hTDCWfg, hTDCWbg, 1, -5.07);
  hTimeW_bgSub->Draw("SAME");

  hSim0->Scale((double)hTimeW_bgSub->GetEntries() / (hSim0->GetEntries()));
  hSim2->Scale((double)hTimeW_bgSub->GetEntries() / (hSim2->GetEntries()));
  hSim3->Scale((double)hTimeW_bgSub->GetEntries() / (hSim3->GetEntries()));
  hSim3->Draw("SAME");

  gPad->Update();

  TLine *l0 = new TLine(0, gPad->GetUymin(), 0, 1.2*gPad->GetUymax());
  l0->SetLineStyle(9);
  l0->Draw();

  TLine *l12 = new TLine(12, gPad->GetUymin(), 12, 1.2*gPad->GetUymax());
  l12->SetLineStyle(9);
  l12->Draw();

  TLegend *l2 = new TLegend(0.6,0.7,0.85,0.85);
  l2->AddEntry(hSim3, "GEANT4", "l");
  l2->AddEntry(hTimeE_bgSub, "T_{East} > T_{West}", "l");
  l2->AddEntry(hTimeW_bgSub, "T_{West} > T_{East}", "l");
  l2->SetTextSize(0.05);
  l2->SetBorderSize(0);
  l2->Draw();
  gPad->SetLogy();

  c2->Print("Figures/2_MPMTimingCompare.pdf");
  c2->Print("Figures/2_MPMTimingCompare.eps");


  // third canvas, let's begin looking at energies
  // this will be our time window upper edge, energy spectra plot
  TCanvas *c31 = new TCanvas("c31", "c31");
  TCanvas *c32 = new TCanvas("c32", "c32");

  vector <TH1D*> hbgErecon;
  vector <TH1D*> hfgErecon;
  vector <TH1D*> hbgSubErecon;
  vector <TCut> timeWindowCuts;

  vector <double> timeUpper;
  timeUpper.push_back(140);
  timeUpper.push_back(20);
  timeUpper.push_back(12);

  TLegend *l3a = new TLegend(0.6,0.7,0.85,0.85);
  l3a->SetBorderSize(0);
  l3a->SetTextSize(0.05);
  TLegend *l3b = new TLegend(0.6,0.7,0.85,0.85);
  l3b->SetBorderSize(0);
  l3b->SetTextSize(0.05);

  for(unsigned int i = 0; i < timeUpper.size(); i++)
  {
    hbgErecon.push_back(new TH1D(Form("bgErecon_%i", i), "BG Erecon", 160, 0, 4000));
    hfgErecon.push_back(new TH1D(Form("fgErecon_%i", i), "FG Erecon", 160, 0, 4000));
    hbgSubErecon.push_back(new TH1D(Form("bgSubErecon_%i", i), "FG - 5.07BG Erecon", 160, 0, 4000));
    hbgErecon[i]->SetStats(kFALSE);
//    hbgErecon[i]->SetBinErrorOption(TH1::kPoisson);
    hbgErecon[i]->Sumw2();
    hfgErecon[i]->SetStats(kFALSE);
//    hfgErecon[i]->SetBinErrorOption(TH1::kPoisson);
    hfgErecon[i]->Sumw2();
    hbgSubErecon[i]->SetStats(kFALSE);
//    hbgSubErecon[i]->SetBinErrorOption(TH1::kPoisson);
    hbgSubErecon[i]->Sumw2();
    timeWindowCuts.push_back(Form("(newTDC2TimeE > -1 && newTDC2TimeE < 5 && newTDC2TimeW > 0 && newTDC2TimeW < %f) || (newTDC2TimeW > -4 && newTDC2TimeW < 2 && newTDC2TimeE > 0 && newTDC2TimeE < %f)", timeUpper[i], timeUpper[i]));

    if(i == 0)
    {
      hbgErecon[i]->SetLineColor(1);
      hfgErecon[i]->SetLineColor(1);
      hbgSubErecon[i]->SetLineColor(1);
    }
    if(i == 1)
    {
      hbgErecon[i]->SetLineColor(2);
      hfgErecon[i]->SetLineColor(2);
      hbgSubErecon[i]->SetLineColor(2);
    }
    if(i == 2)
    {
      hbgErecon[i]->SetLineColor(4);
      hfgErecon[i]->SetLineColor(4);
      hbgSubErecon[i]->SetLineColor(4);
    }

    c31->cd();
    if(i == 0)
    {
      hbgErecon[i]->SetTitle("");
      hbgErecon[i]->GetXaxis()->SetTitle("E_{e^{+}e^{-}} (keV)");
      hbgErecon[i]->GetXaxis()->CenterTitle();
      hbgErecon[i]->GetYaxis()->SetTitle("Counts");
      hbgErecon[i]->GetYaxis()->CenterTitle();
      hbgErecon[i]->GetYaxis()->SetRangeUser(0.1, 5000);
      bgchain->Draw(Form("Erecon_ee >> bgErecon_%i", i), basicCut && fiducialCut && "Erecon_ee > 0" && timeWindowCuts[i], "HIST");
    }
    else
    {
      bgchain->Draw(Form("Erecon_ee >> bgErecon_%i", i), basicCut && fiducialCut && "Erecon_ee > 0" && timeWindowCuts[i], "HIST SAME");
    }
    l3a->AddEntry(hbgErecon[i], Form("(0, %i) ns", (int)timeUpper[i]), "l");

    c32->cd();
    if(i == 0)
    {
      hfgErecon[i]->SetTitle("");
      hfgErecon[i]->GetXaxis()->SetTitle("E_{e^{+}e^{-}} (keV)");
      hfgErecon[i]->GetXaxis()->CenterTitle();
      hfgErecon[i]->GetYaxis()->SetTitle("");
      hfgErecon[i]->GetYaxis()->SetRangeUser(1, 50000);
      fgchain->Draw(Form("Erecon_ee >> fgErecon_%i", i), basicCut && fiducialCut && "Erecon_ee > 0" && timeWindowCuts[i], "HIST");
    }
    else
    {
      fgchain->Draw(Form("Erecon_ee >> fgErecon_%i", i), basicCut && fiducialCut && "Erecon_ee > 0" && timeWindowCuts[i], "HIST SAME");
    }
    l3b->AddEntry(hfgErecon[i], Form("(0, %i) ns", (int)timeUpper[i]), "l");

  }

  for(unsigned int i = 0; i < timeUpper.size(); i++)
  {
    for(int j = 0; j < hbgSubErecon[i]->GetNbinsX(); j++)
    {
      hbgErecon[i]->SetBinError(j, SetPoissonErrors(hbgErecon[i]->GetBinContent(j)));
      hfgErecon[i]->SetBinError(j, SetPoissonErrors(hfgErecon[i]->GetBinContent(j)));
/*
      cout << "At hist for time cut i = " << i
	   << ", we have at bin j = " << j
	   << ", BG counts = " << hbgErecon[i]->GetBinContent(j) << " +/- " << hbgErecon[i]->GetBinError(j)
	   << ", FG counts = " << hfgErecon[i]->GetBinContent(j) << " +/- " << hfgErecon[i]->GetBinError(j) << endl;
*/
    }
  }

  for(unsigned int i = 0; i < timeUpper.size(); i++)
  {
    hbgSubErecon[i]->Add(hfgErecon[i], hbgErecon[i], 1, -5.07);
  }

  c31->cd();
  gPad->Update();

  TLine *lE0a = new TLine(0, 0, 0, 5000);
  lE0a->SetLineStyle(9);

  TLine *lE800a = new TLine(800, 0, 800, 5000);
  lE800a->SetLineStyle(9);

  gPad->SetLogy();
  lE0a->Draw();
  lE800a->Draw();
  l3a->Draw();


  c32->cd();
  gPad->Update();

  TLine *lE0b = new TLine(0, gPad->GetUymin(), 0, 50000);
  lE0b->SetLineStyle(9);

  TLine *lE800b = new TLine(800, gPad->GetUymin(), 800, 50000);
  lE800b->SetLineStyle(9);

  gPad->SetLogy();
  lE0b->Draw();
  lE800b->Draw();
  l3b->Draw();



  c31->Print("Figures/31_Erecon_TimeWindowEndpoint_BG.pdf");
  c31->Print("Figures/31_Erecon_TimeWindowEndpoint_BG.eps");
  c32->Print("Figures/32_Erecon_TimeWindowEndpoint_FG.pdf");
  c32->Print("Figures/32_Erecon_TimeWindowEndpoint_FG.eps");

  // fourth canvas, energies with background subtraction and statistics propagated and acceptances from Brad F.
  TCanvas *c4 = new TCanvas("c4","c4");
  c4->Divide(2,1);

  TH1D *hbgErecon_withCuts = new TH1D("bgEreconfull", Form("BG Erecon: non-STP Events (%f, %f) ns", windowLower, windowUpper), 8, -6, 794);
  hbgErecon_withCuts->Sumw2();
  hbgErecon_withCuts->GetXaxis()->SetTitle("Erecon_ee (keV)");
  TH1D *hfgErecon_withCuts = new TH1D("fgEreconfull", Form("FG Erecon: non-STP Events (%f, %f) ns", windowLower, windowUpper), 8, -6, 794);
  hfgErecon_withCuts->Sumw2();
  hfgErecon_withCuts->GetXaxis()->SetTitle("Erecon_ee (keV)");


  c4->cd(1);
  bgchain->Draw("Erecon_ee >> bgEreconfull", basicCut && fiducialCut && (inESTPAndTimeWinCut || inWSTPAndTimeWinCut) && energyCut);
  c4->cd(2);
  fgchain->Draw("Erecon_ee >> fgEreconfull", basicCut && fiducialCut && (inESTPAndTimeWinCut || inWSTPAndTimeWinCut) && energyCut);

  cout << "Setting Poisson error bars..." << endl;
  for(int i = 0; i < hfgErecon_withCuts->GetNbinsX(); i++)
  {
    hbgErecon_withCuts->SetBinError(i, SetPoissonErrors(hbgErecon_withCuts->GetBinContent(i)));
    hfgErecon_withCuts->SetBinError(i, SetPoissonErrors(hfgErecon_withCuts->GetBinContent(i)));

//    cout << "At bin i = " << i << " BG counts = " << hbgErecon_withCuts->GetBinContent(i) << ", FG counts = " << hfgErecon_withCuts->GetBinContent(i) << endl;
  }


//  c4->Print("4_BG_FG_histograms_withErrors.pdf");
//  c4->Print("4_BG_FG_histograms_withErrors.eps");

  // intermediate 45 canvas
  TCanvas *c45 = new TCanvas ("c45","c45");
  c45->cd();

  TLegend *l45 = new TLegend(0.2,0.65,0.6,0.85);

  TH1D* hTestSignal = new TH1D("hTestSignal", "hTestSignal", 160, 0, 4000);
  FillInCountsFromFile("test_signal_644KeV.dat", hTestSignal);
  hTestSignal->Scale(0.1834);
  hTestSignal->SetLineColor(2);
  hTestSignal->SetStats(kFALSE);
  hTestSignal->SetTitle("");
  hTestSignal->GetXaxis()->SetTitle("E_{e^{+}e^{-}} (keV)");
  hTestSignal->GetXaxis()->CenterTitle();
  hTestSignal->GetXaxis()->SetRangeUser(0, 800);
  hTestSignal->GetYaxis()->SetTitle("Counts");
  hTestSignal->GetYaxis()->CenterTitle();
  hTestSignal->GetYaxis()->SetRangeUser(1, 100000);
  gPad->SetLogy();
  hTestSignal->Draw();

  TH1D* hTestSignal322 = new TH1D("hTestSignal322", "hTestSignal322", 160, 0, 4000);
  FillInCountsFromFile("anom_peak_322.dat", hTestSignal322);
  hTestSignal322->Scale(0.1124);
  hTestSignal322->SetLineColor(6);
  hTestSignal322->Draw("SAME");

  for(int i = 0; i <= hbgSubErecon[2]->GetNbinsX(); i++)
  {
    if(hbgSubErecon[2]->GetBinContent(i) < 0)
    {
      hbgSubErecon[2]->SetBinContent(i, 0);
    }
/*
    cout << "Final i = 2 BG sub histogram has at bin j = " << i
	 << ", FG - 5.07BG counts = " << hbgSubErecon[2]->GetBinContent(i)
	 << " +/- " << hbgSubErecon[2]->GetBinError(i) << endl;
*/
  }

//  hbgSubErecon[0]->Draw("HIST SAME");
//  hbgSubErecon[1]->Draw("HIST SAME");
  gStyle->SetErrorX(0);
  hbgSubErecon[2]->SetMarkerStyle(21);
  hbgSubErecon[2]->SetMarkerColor(4);
  hbgSubErecon[2]->Draw("E0 SAME");

  l45->AddEntry(hTestSignal322, "Simulated Signal at 322 keV", "l");
  l45->AddEntry(hTestSignal, "Simulated Signal at 644 keV", "l");
  l45->AddEntry(hbgSubErecon[2], "(0, 12) ns Data", "l");
  l45->SetTextSize(0.05);
  l45->SetBorderSize(0);
  l45->Draw();

  c45->Print("Figures/45_BGSub_andSignal.pdf");
  c45->Print("Figures/45_BGSub_andSignal.png");
  c45->Print("Figures/45_BGSub_andSignal.eps");


  // fifth canvas, does background subtraction and acceptances, because ROOT is tough to work with.
  TCanvas *c5 = new TCanvas("c5","c5");
  c5->cd();

  gStyle->SetOptStat(0000);
  TH1D *hErecon_bgSub = new TH1D("fullCuts", "BG subtracted Erecon_ee", 8, -6, 794);
  hErecon_bgSub->SetLineColor(2);
  hErecon_bgSub->Sumw2();
  hErecon_bgSub->Add(hfgErecon_withCuts, hbgErecon_withCuts, 1, -5.07);	// 5.07 comes from 860262/169717 live time ratio

  double totalEntries = 0;
  for(int i = 0; i < hErecon_bgSub->GetNbinsX(); i++)
  {
    totalEntries = totalEntries + hErecon_bgSub->GetBinContent(i);
  }
  hErecon_bgSub->SetEntries(totalEntries);

  // sixth canvas, for acceptances
  TCanvas *c6 = new TCanvas("c6","c6");
  FillInAcceptancesFromFile("m094002MeV_4502123TotalCounts.dat", 4502123, 94);
  FillInAcceptancesFromFile("m144002MeV_455639TotalCounts.dat", 455639, 144);
  FillInAcceptancesFromFile("m244002MeV_460822TotalCounts.dat", 460822, 244);
  FillInAcceptancesFromFile("m344002MeV_466564TotalCounts.dat", 466564, 344);
  FillInAcceptancesFromFile("m444002MeV_471023TotalCounts.dat", 471023, 444);
  FillInAcceptancesFromFile("m544002MeV_475630TotalCounts.dat", 475630, 544);
  FillInAcceptancesFromFile("m644002MeV_478192TotalCounts.dat", 478192, 644);

  TGraph *gAccept = CreateAcceptancesGraph(0.0, 12.0);
  double xAcceptKin = 0;
  double yAcceptKin = 0;
  // adjust the kinematic acceptance for G4 sim, for total acceptance
  for(int i = 0; i < gAccept->GetN(); i++)
  {
    gAccept->GetPoint(i, xAcceptKin, yAcceptKin);
    gAccept->SetPoint(i, xAcceptKin, yAcceptKin*0.845);
  }
  gAccept->SetMarkerStyle(21);
  gAccept->SetMarkerColor(1);
  gAccept->SetLineColor(1);
  gAccept->GetHistogram()->SetTitle("");
  gAccept->GetHistogram()->GetXaxis()->SetTitle("E_{e^{+}e^{-}} (keV)");
  gAccept->GetHistogram()->GetXaxis()->CenterTitle();
  gAccept->GetHistogram()->GetYaxis()->SetTitle("Fraction Accepted");
  gAccept->GetHistogram()->GetYaxis()->CenterTitle();
  gAccept->GetHistogram()->GetXaxis()->SetRangeUser(0, 800);
  gAccept->GetHistogram()->GetYaxis()->SetRangeUser(0, 0.22);
  gAccept->Draw("AC");

  c6->Print("Figures/6_Acceptance.pdf");
  c6->Print("Figures/6_Acceptance.png");
  c6->Print("Figures/6_Acceptance.eps");

  double G4x[4] = {200, 321, 480, 640};
  double G4y[4] = {82560.0/96558.0, 82212.0/96721.0, 80955.0/96271.0, 80283.0/96067.0};
  TGraph *gG4 = new TGraph(4, G4x, G4y);

  double gEx[7] = {94, 144, 244, 344, 444, 544, 644};
  // 50.0 KeV is our bin width on either side of the center
  double sigmaEKeV[7] = {50.0/15.33, 50.0/18.97, 50.0/24.70, 50.0/29.32, 50.0/33.32, 50.0/36.88, 50.0/40.12};
  // calculated by taking a Gaussian probability calculator and plugging in (denom) as SD, range from -50 to 50
  double gEy[7] = {0.9989, 0.9916, 0.9571, 0.9119, 0.8665, 0.8248, 0.7873};

  TGraph *gEResolution = new TGraph(7, gEx, gEy);

  c5->cd();
  double xAccept = 0;
  double yAccept = 0;
  double xEResolution = 0;
  double yEResolution = 0;

  double finalEntries = 0;

  TH1D* hFinalNumbers = new TH1D("hFinal", "Final: Cuts and scaled for Acceptances", 8, -6, 794);
  hFinalNumbers->SetStats(kFALSE);
  hFinalNumbers->SetTitle("");
  hFinalNumbers->GetXaxis()->SetTitle("E_{e^{+}e^{-}} (keV)");
  hFinalNumbers->GetYaxis()->SetTitle("Counts");
  hFinalNumbers->GetYaxis()->SetRangeUser(0, 1000);
  hFinalNumbers->SetLineColor(1);

  for(int i = 0; i < hFinalNumbers->GetNbinsX(); i++)
  {
    gAccept->GetPoint(i, xAccept, yAccept);
    cout << "At i = " << i << " we have xAccept = " << xAccept << " and yAccept = " << yAccept << endl;

    gEResolution->GetPoint(i, xEResolution, yEResolution);
    cout << "At i = 0" << i << " we have xEResolution = " << xEResolution << " and yEResolution = " << yEResolution << endl;

    hFinalNumbers->SetBinContent(i, (double)hErecon_bgSub->GetBinContent(i) / (yAccept*yEResolution));
    hFinalNumbers->SetBinError(i, (double)hErecon_bgSub->GetBinError(i) / (yAccept*yEResolution));

    finalEntries = finalEntries + ((double)hErecon_bgSub->GetBinContent(i) / (yAccept*yEResolution));
  }

  hFinalNumbers->SetEntries(finalEntries);

  hFinalNumbers->Draw();
  hErecon_bgSub->SetBinContent(hErecon_bgSub->GetNbinsX(), 0);
  hErecon_bgSub->SetBinError(hErecon_bgSub->GetNbinsX(), 0);
  hErecon_bgSub->Draw("SAME");


  cout << "Total counts in histogram after final acceptances propagated: " << hFinalNumbers->GetEntries() << endl;
  cout << "Total counts in histogram after background subtraction: " << hErecon_bgSub->GetEntries() << endl;

  for(int i = 0; i < hFinalNumbers->GetNbinsX(); i++)
  {
    cout << "hFinalNumbers at x = " << hFinalNumbers->GetBinCenter(i)
	 << " has counts = " << hFinalNumbers->GetBinContent(i)
	 << " and associated symmetric error = " << hFinalNumbers->GetBinError(i) << endl;
  }

  c5->Print("Figures/5_BGSub_WithAcceptances_finalNumbers.pdf");
  c5->Print("Figures/5_BGSub_WithAcceptances_finalNumbers.eps");

  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}


double SetPoissonErrors(int counts)
{
  double upperErrBar = 0;

  // all values taken from PDG paper rpp2017-rev-statistics.pdf
  // Using Table 40.3, for 95% upper one-sided limit
  if(counts == 0)
  {
    upperErrBar = 3;
  }
  else if(counts == 1)
  {
    upperErrBar = 4.74 - 1;
  }
  else if(counts == 2)
  {
    upperErrBar = 6.30 - 2;
  }
  else if(counts == 3)
  {
    upperErrBar = 7.75 - 3;
  }
  else if(counts == 4)
  {
    upperErrBar = 9.15 - 4;
  }
  else if(counts == 5)
  {
    upperErrBar = 10.51 - 5;
  }
  else if(counts == 6)
  {
    upperErrBar = 11.84 - 6;
  }
  else if(counts == 7)
  {
    upperErrBar = 13.15 - 7;
  }
  else if(counts == 8)
  {
    upperErrBar = 14.43 - 8;
  }
  else if(counts == 9)
  {
    upperErrBar = 15.71 - 9;
  }
  else if(counts == 10)
  {
    upperErrBar = 16.96 - 10;
  }
  else if(counts > 10)
  {
    upperErrBar = 2*sqrt(counts);
  }

  return (upperErrBar / 2.0);
}

void FillInCountsFromFile(TString fileName, TH1D* hToFill)
{
  double tmpEnergy;
  double tmpCounts;

  string buf;
  ifstream infile;
  cout << "The file being opened is: " << fileName << endl;
  infile.open(Form("Acceptances/%s", fileName.Data()));

  if(!infile.is_open())
    cout << "Problem opening " << fileName << endl;

  while(getline(infile, buf))
  {
    istringstream bufstream(buf);
    if(!bufstream.eof())
    {
      bufstream >> tmpEnergy >> tmpCounts;

      for(int i = 0; i < tmpCounts; i++)
      {
        hToFill->Fill(tmpEnergy*1000);
      }
    }
  }

  cout << "Done filling counts into signal histogram." << endl;
}


void FillInAcceptancesFromFile(TString fileName, double totalCounts, int flag)
{
  vector <double> time;
  vector <double> counts;

  double tmpTime;
  double tmpCounts;

  string buf;
  ifstream infile;
  cout << "The file being opened is: " << fileName << endl;
  infile.open(Form("Acceptances/%s", fileName.Data()));

  if(!infile.is_open())
    cout << "Problem opening " << fileName << endl;

  while(getline(infile, buf))
  {
    istringstream bufstream(buf);
    if(!bufstream.eof())
    {
      bufstream >> tmpTime >> tmpCounts;

      time.push_back(tmpTime);
      counts.push_back(tmpCounts / (totalCounts));

    }
  }

  // a bunch of different cases because I don't want to use addreses
  if(flag == 94)
  {
    t094 = time;
    n094 = counts;
  }
  else if(flag == 144)
  {
    t144 = time;
    n144 = counts;
  }
  else if(flag == 244)
  {
    t244 = time;
    n244 = counts;
  }
  else if(flag == 344)
  {
    t344 = time;
    n344 = counts;
  }
  else if(flag == 444)
  {
    t444 = time;
    n444 = counts;
  }
  else if(flag == 544)
  {
    t544 = time;
    n544 = counts;
  }
  else if(flag == 644)
  {
    t644 = time;
    n644 = counts;
  }
  else
  {
    cout << "Flag makes no sense, not sure which vectors to fill!" << endl;
  }
}

TGraph* CreateAcceptancesGraph(double timeLower, double timeUpper)
{
  int nPoints = 0;
  vector <double> energy;
  vector <double> acceptance;

  energy.push_back(0);
  acceptance.push_back(0);
  nPoints = nPoints + 1;

  double accept094 = 0;
  for(unsigned int i = 0; i < t094.size(); i++)
  {
    if(t094[i] < timeLower || t094[i] > timeUpper)
    {
      continue;
    }
    accept094 = accept094 + n094[i];
  }
  energy.push_back(94);
  acceptance.push_back(accept094);
  nPoints = nPoints + 1;

  double accept144 = 0;
  for(unsigned int i = 0; i < t144.size(); i++)
  {
    if(t144[i] < timeLower || t144[i] > timeUpper)
    {
      continue;
    }
    accept144 = accept144 + n144[i];
  }
  energy.push_back(144);
  acceptance.push_back(accept144);
  nPoints = nPoints + 1;

  double accept244 = 0;
  for(unsigned int i = 0; i < t244.size(); i++)
  {
    if(t244[i] < timeLower || t244[i] > timeUpper)
    {
      continue;
    }
    accept244 = accept244 + n244[i];
  }
  energy.push_back(244);
  acceptance.push_back(accept244);
  nPoints = nPoints + 1;

  double accept344 = 0;
  for(unsigned int i = 0; i < t344.size(); i++)
  {
    if(t344[i] < timeLower || t344[i] > timeUpper)
    {
      continue;
    }
    accept344 = accept344 + n344[i];
  }
  energy.push_back(344);
  acceptance.push_back(accept344);
  nPoints = nPoints + 1;

  double accept444 = 0;
  for(unsigned int i = 0; i < t444.size(); i++)
  {
    if(t444[i] < timeLower || t444[i] > timeUpper)
    {
      continue;
    }
    accept444 = accept444 + n444[i];
  }
  energy.push_back(444);
  acceptance.push_back(accept444);
  nPoints = nPoints + 1;

  double accept544 = 0;
  for(unsigned int i = 0; i < t544.size(); i++)
  {
    if(t544[i] < timeLower || t544[i] > timeUpper)
    {
      continue;
    }
    accept544 = accept544 + n544[i];
  }
  energy.push_back(544);
  acceptance.push_back(accept544);
  nPoints = nPoints + 1;

  double accept644 = 0;
  for(unsigned int i = 0; i < t644.size(); i++)
  {
    if(t644[i] < timeLower || t644[i] > timeUpper)
    {
      continue;
    }
    accept644 = accept644 + n644[i];
  }
  energy.push_back(644);
  acceptance.push_back(accept644);
  nPoints = nPoints + 1;

  energy.push_back(694);
  acceptance.push_back(0.222);
  nPoints = nPoints + 1;

  energy.push_back(744);
  acceptance.push_back(0.229);
  nPoints = nPoints + 1;

  energy.push_back(794);
  acceptance.push_back(0.233);
  nPoints = nPoints + 1;


  for(int i = 0; i < nPoints; i++)
  {
    cout << "Points being fed into the graph: (" << energy[i] << ", " << acceptance[i] << ")" << endl;

  }

  TGraph *g = new TGraph(nPoints, &(energy[0]), &(acceptance[0]));

  return g;
}
