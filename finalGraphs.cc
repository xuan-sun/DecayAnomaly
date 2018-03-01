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

double SetPoissonErrors(int counts);

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

  TH1D *hTime_sim = new TH1D("hSim", "MPM Simulated Timing, East Hits First", 320, -20, 140);
  TH1D *hTime_tdcE = new TH1D("hTDC2TimeE", "TDC converted to Time East", 320, -20, 140);
  TH1D *hTime_tdcW = new TH1D("hTDC2TimeW", "TDC converted to Time West", 320, -20, 140);


  TH1D *hTDCEfg = new TH1D("hTDCEfg", "FG TDC Spectra", 320, -20, 140);
  hTDCEfg->GetXaxis()->SetTitle("Time (ns)");
  hTDCEfg->SetLineColor(2);

  TH1D *hTDCWfg = new TH1D("hTDCWfg", "FG TDC Spectra", 320, -20, 140);
  hTDCWfg->SetLineColor(4);

  TH1D *hTDCEbg = new TH1D("hTDCEbg", "BG TDC Spectra", 320, -20, 140);
  hTDCEbg->GetXaxis()->SetTitle("Time (ns)");
  hTDCEbg->SetLineColor(2);

  TH1D *hTDCWbg = new TH1D("hTDCWbg", "BG TDC Spectra", 320, -20, 140);
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

  TH1D *hSim0 = new TH1D("hSim0sigma", "hSim", 320, -20, 140);
  hSim0->SetLineColor(6);
  TH1D *hSim2 = new TH1D("hSim2sigma", "hSim", 320, -20, 140);
  hSim2->SetLineColor(3);

  SimEvent evt;
  mpm_sim_chain->GetBranch("Edep")->GetLeaf("EdepE")->SetAddress(&evt.Edep_EdepE);
  mpm_sim_chain->GetBranch("Edep")->GetLeaf("EdepW")->SetAddress(&evt.Edep_EdepW);
  mpm_sim_chain->GetBranch("EdepQ")->GetLeaf("EdepQE")->SetAddress(&evt.EdepQ_EdepQE);
  mpm_sim_chain->GetBranch("EdepQ")->GetLeaf("EdepQW")->SetAddress(&evt.EdepQ_EdepQW);
  mpm_sim_chain->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyE")->SetAddress(&evt.MWPCEnergy_MWPCEnergyE);
  mpm_sim_chain->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyW")->SetAddress(&evt.MWPCEnergy_MWPCEnergyW);
  mpm_sim_chain->GetBranch("time")->GetLeaf("timeE")->SetAddress(&evt.time_timeE);
  mpm_sim_chain->GetBranch("time")->GetLeaf("timeW")->SetAddress(&evt.time_timeW);

  for(unsigned int i = 0; i < mpm_sim_chain->GetEntries(); i++)
  {
    mpm_sim_chain->GetEntry(i);

    if(evt.EdepQ_EdepQE > 0 && evt.EdepQ_EdepQW > 0 && evt.MWPCEnergy_MWPCEnergyE > 0 && evt.MWPCEnergy_MWPCEnergyW > 0
        && evt.time_timeE < 200 && evt.time_timeW < 200 && evt.time_timeE < evt.time_timeW)
    {
      hSim0->Fill(evt.time_timeW - evt.time_timeE);
      hSim2->Fill(engine->Gaus(evt.time_timeW - evt.time_timeE, 2));
    }

    if(i % 100000 == 0)
    {
      cout << "Filling event " << i << "/" << mpm_sim_chain->GetEntries() << endl;
    }
  }

  TH1D *hTimeE_bgSub = new TH1D("BGSubE", "BG Subtracted Time East", 320, -20, 140);
  hTimeE_bgSub->SetLineColor(2);
  hTimeE_bgSub->Add(hTDCEfg, hTDCEbg, 1, -5.07);
  hTimeE_bgSub->SetStats(kFALSE);
  hTimeE_bgSub->Draw();

  TH1D *hTimeW_bgSub = new TH1D("BGSubW", "BG Subtracted Time West", 320, -20, 140);
  hTimeW_bgSub->SetLineColor(4);
  hTimeW_bgSub->Add(hTDCWfg, hTDCWbg, 1, -5.07);
  hTimeW_bgSub->Draw("SAME");

  hSim0->Scale((double)hTimeW_bgSub->GetEntries() / (hSim0->GetEntries()));
  hSim0->Draw("SAME");
  hSim2->Scale((double)hTimeW_bgSub->GetEntries() / (hSim2->GetEntries()));
  hSim2->Draw("SAME");

  TLegend *l2 = new TLegend(0.6,0.75,0.9,0.9);
  l2->AddEntry(hSim0, "MPM G4 Sim no smearing", "l");
  l2->AddEntry(hSim2, "MPM G4 Sim 2ns sigma", "l");
  l2->AddEntry(hTimeE_bgSub, "FG - 5.07BG East", "l");
  l2->AddEntry(hTimeW_bgSub, "FG - 5.07BG West", "l");
  l2->Draw();
  gPad->SetLogy();

  c2->Print("2_MPMTimingCompare.pdf");

  // third canvas, let's begin looking at energies
  // this will be our time window upper edge, energy spectra plot
  TCanvas *c3 = new TCanvas("c3", "c3");
  c3->Divide(2,1);

  vector <TH1D*> hbgErecon;
  vector <TH1D*> hfgErecon;
  vector <TCut> timeWindowCuts;

  vector <double> timeUpper;
  timeUpper.push_back(140);
  timeUpper.push_back(20);
  timeUpper.push_back(16);
  timeUpper.push_back(14);
  timeUpper.push_back(12);

  TLegend *l3a = new TLegend(0.6,0.6,0.9,0.9);
  TLegend *l3b = new TLegend(0.6,0.6,0.9,0.9);

  for(unsigned int i = 0; i < timeUpper.size(); i++)
  {
    hbgErecon.push_back(new TH1D(Form("bgErecon_%i", i), "BG Erecon", 160, 0, 4000));
    hfgErecon.push_back(new TH1D(Form("fgErecon_%i", i), "FG Erecon", 160, 0, 4000));
    hbgErecon[i]->SetLineColor(i+1);
    hbgErecon[i]->SetStats(kFALSE);
    hfgErecon[i]->SetLineColor(i+1);
    hfgErecon[i]->SetStats(kFALSE);
    timeWindowCuts.push_back(Form("(newTDC2TimeE > -1 && newTDC2TimeE < 5 && newTDC2TimeW > 0 && newTDC2TimeW < %f) || (newTDC2TimeW > -4 && newTDC2TimeW < 2 && newTDC2TimeE > 0 && newTDC2TimeE < %f)", timeUpper[i], timeUpper[i]));

    c3->cd(1);
    if(i == 0)
    {
      bgchain->Draw(Form("Erecon_ee >> bgErecon_%i", i), basicCut && fiducialCut && "Erecon_ee > 0" && timeWindowCuts[i]);
    }
    else
    {
      bgchain->Draw(Form("Erecon_ee >> bgErecon_%i", i), basicCut && fiducialCut && "Erecon_ee > 0" && timeWindowCuts[i], "SAME");
    }
    l3a->AddEntry(hbgErecon[i], Form("(0, %f) ns acceptance", timeUpper[i]), "l");

    c3->cd(2);
    if(i == 0)
    {
      fgchain->Draw(Form("Erecon_ee >> fgErecon_%i", i), basicCut && fiducialCut && "Erecon_ee > 0" && timeWindowCuts[i]);
    }
    else
    {
      fgchain->Draw(Form("Erecon_ee >> fgErecon_%i", i), basicCut && fiducialCut && "Erecon_ee > 0" && timeWindowCuts[i], "SAME");
    }
    l3b->AddEntry(hfgErecon[i], Form("(0, %f) ns acceptance", timeUpper[i]), "l");
  }

  c3->cd(1);
  l3a->Draw();
  c3->cd(2);
  l3b->Draw();
  gPad->SetLogy();

  c3->Print("3_Erecon_TimeWindowEndpointCompare.pdf");

  // fourth canvas, energies with background subtraction and statstics.
  TCanvas *c4 = new TCanvas("c4","c4");
  c4->Divide(3,1);

  TH1D *hbgErecon_withCuts = new TH1D("bgEreconfull", Form("BG Erecon: non-STP Events (%f, %f) ns", windowLower, windowUpper), 10, 0, 1000);
  hbgErecon_withCuts->Sumw2();
  hbgErecon_withCuts->GetXaxis()->SetTitle("Erecon_ee (KeV)");
  TH1D *hfgErecon_withCuts = new TH1D("fgEreconfull", Form("FG Erecon: non-STP Events (%f, %f) ns", windowLower, windowUpper), 10, 0, 1000);
  hfgErecon_withCuts->Sumw2();
  hfgErecon_withCuts->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  c4->cd(1);
  bgchain->Draw("Erecon_ee >> bgEreconfull", basicCut && fiducialCut && (inESTPAndTimeWinCut || inWSTPAndTimeWinCut) && energyCut);
  c4->cd(2);
  fgchain->Draw("Erecon_ee >> fgEreconfull", basicCut && fiducialCut && (inESTPAndTimeWinCut || inWSTPAndTimeWinCut) && energyCut);

  cout << "Setting Poisson error bars..." << endl;
  for(int i = 0; i < hfgErecon_withCuts->GetNbinsX(); i++)
  {
    hbgErecon_withCuts->SetBinError(i, SetPoissonErrors(hbgErecon_withCuts->GetBinContent(i)));
    hfgErecon_withCuts->SetBinError(i, SetPoissonErrors(hfgErecon_withCuts->GetBinContent(i)));
  }

  gStyle->SetOptStat(11);
  TH1D *hErecon_bgSub = new TH1D("fullCuts", "BG subtracted Erecon_ee", 10, 0, 1000);
  hErecon_bgSub->Sumw2();
  hErecon_bgSub->Add(hfgErecon_withCuts, hbgErecon_withCuts, 1, -5.07);	// 5.07 comes from 860262/169717 live time ratio

  c4->cd(3);
  hErecon_bgSub->Draw();

  double totalEntries = 0;
  for(int i = 0; i < hErecon_bgSub->GetNbinsX(); i++)
  {
    totalEntries = totalEntries + hErecon_bgSub->GetBinContent(i);
  }
  hErecon_bgSub->SetEntries(totalEntries);

  c4->Print("4_Erecon_BG_FG_finalCuts.pdf");


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
