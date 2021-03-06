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
  if(argc < 9)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable) (E Low 1) (E High 1) (W Low 1) (W High 1) (E Low 2) (E High 2) (W High 2) (W High 2)" << endl;
    return 0;
  }

  gRandom->SetSeed(0);

  // read in the arguments
  double timeLowerEdgeE1 = atof(argv[1]);
  double timeUpperEdgeE1 = atof(argv[2]);
  double timeLowerEdgeW1 = atof(argv[3]);
  double timeUpperEdgeW1 = atof(argv[4]);
  double timeLowerEdgeE2 = atof(argv[5]);
  double timeUpperEdgeE2 = atof(argv[6]);
  double timeLowerEdgeW2 = atof(argv[7]);
  double timeUpperEdgeW2 = atof(argv[8]);

  TString stp = "E";
  TString roi = "W";

  TChain *bgchain = new TChain("pass3");
  TChain *fgchain = new TChain("pass3");

  bgchain->Add("TimeCalibrated_BGRuns_type1_fixed_v3.root");
  fgchain->Add("TimeCalibrated_FGRuns_type1_fixed_v3.root");

  // define all the cuts we will use later.
  TCut basicCut = "(PID == 1 && badTimeFlag == 0)";
  TCut energyCut = "(Erecon_ee > 0 && Erecon_ee < 800)";
  TCut fiducialCut = "(((xE.center)*(xE.center) + (yE.center)*(yE.center) < 49*49) && ((xW.center)*(xW.center) + (yW.center)*(yW.center) < 49*49))";
  TCut simulationCut = "EdepQ.EdepQE > 0 && EdepQ.EdepQW > 0 && MWPCEnergy.MWPCEnergyE > 0 && MWPCEnergy.MWPCEnergyW > 0 && time.timeE < 200 && time.timeW < 200 && timeE < timeW";
  TCut time1STPCut = Form("newTDC2TimeE < %f && newTDC2TimeE > %f && newTDC2TimeW > %f && newTDC2TimeW < %f", timeUpperEdgeE1, timeLowerEdgeE1, timeLowerEdgeW1, timeUpperEdgeW1);
  TCut time2STPCut = Form("(newTDC2TimeE < %f && newTDC2TimeE > %f && newTDC2TimeW > %f && newTDC2TimeW < %f) || (newTDC2TimeE < %f && newTDC2TimeE > %f && newTDC2TimeW > %f && newTDC2TimeW < %f)",
			timeUpperEdgeE1, timeLowerEdgeE1, timeLowerEdgeW1, timeUpperEdgeW1,
			timeUpperEdgeE2, timeLowerEdgeE2, timeLowerEdgeW2, timeUpperEdgeW2);


  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->Divide(2,1);
  c1->cd(1);

  TH1D *hbgErecon_timeWin = new TH1D("bgErecon_timeWin", Form("BG Erecon_ee global time: %f < E < %f ns, %f < W < %f ns",
				timeLowerEdgeE1, timeUpperEdgeE1, timeLowerEdgeW1, timeUpperEdgeW1), 160, 0, 4000);
  hbgErecon_timeWin->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  TH1D *hfgErecon_timeWin = new TH1D("fgErecon_timeWin", Form("FG Erecon_ee global time: %f < E < %f ns, %f < W < %f ns",
				timeLowerEdgeE1, timeUpperEdgeE1, timeLowerEdgeW1, timeUpperEdgeW1), 160, 0, 4000);
  hfgErecon_timeWin->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  bgchain->Draw("Erecon_ee >> bgErecon_timeWin", basicCut && fiducialCut && time2STPCut);

  c1->cd(2);
  fgchain->Draw("Erecon_ee >> fgErecon_timeWin", basicCut && fiducialCut && time2STPCut);

//  c1->Print("1_Erecon_timingWindow.pdf");

  // second canvas for plots
  TCanvas *c2 = new TCanvas("c2", "c2");
  c2->Divide(3,1);

  TH1D *hbgSpectra = new TH1D("bgfull", Form("BG Erecon_ee, all cuts: %f < E < %f, %f < W < %f",
				timeLowerEdgeE1, timeUpperEdgeE1, timeLowerEdgeW1, timeUpperEdgeW1), 10, 0, 1000);
  hbgSpectra->Sumw2();
  hbgSpectra->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  TH1D *hfgSpectra = new TH1D("fgfull", Form("FG Erecon_ee, all cuts: %f < E < %f, %f < W < %f",
				timeLowerEdgeE1, timeUpperEdgeE1, timeLowerEdgeW1, timeUpperEdgeW1), 10, 0, 1000);
  hfgSpectra->Sumw2();
  hfgSpectra->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  c2->cd(1);
  bgchain->Draw("Erecon_ee >> bgfull", basicCut && fiducialCut && time2STPCut && energyCut);

  c2->cd(2);
  fgchain->Draw("Erecon_ee >> fgfull", basicCut && fiducialCut && time2STPCut && energyCut);

  // setting Poisson error bars for the BG and FG histograms
  cout << "Setting Poisson error bars..." << endl;
  for(int i = 0; i < hfgSpectra->GetNbinsX(); i++)
  {
    hbgSpectra->SetBinError(i, SetPoissonErrors(hbgSpectra->GetBinContent(i)));
    hfgSpectra->SetBinError(i, SetPoissonErrors(hfgSpectra->GetBinContent(i)));
  }

  TH1D *h_fullCuts = new TH1D("fullCuts", "BG subtracted Erecon_ee", 10, 0, 1000);
  h_fullCuts->Sumw2();
  h_fullCuts->Add(hfgSpectra, hbgSpectra, 1, -5.07);	// 5.07 comes from 860262/169717 live time ratio

  c2->cd(3);
  h_fullCuts->Draw();

//  c2->Print("2_BGsubtracted_timeWindow.pdf");

  // fifth canvas
  TCanvas *c5 = new TCanvas("c5", "c5");
  c5->Divide(2,1);

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

  c5->cd(1);
//  bgchain->Draw("(-1)*(182.27/4096.0)*(TDCE - 38 - 2917) >> hTDCEbg", basicCut && fiducialCut && time2STPCut && energyCut);
//  bgchain->Draw("(-1)*(182.27/4096.0)*(TDCW + 38 - 3138) >> hTDCWbg", basicCut && fiducialCut && time2STPCut && energyCut, "SAME");
  bgchain->Draw("newTDC2TimeE >> hTDCEbg", basicCut && fiducialCut && time2STPCut && energyCut);
  bgchain->Draw("newTDC2TimeW >> hTDCWbg", basicCut && fiducialCut && time2STPCut && energyCut, "SAME");

  TLegend *lbg = new TLegend(0.1,0.8,0.4,0.9);
  lbg->AddEntry(hTDCEbg,"Calibrated TDCE","l");
  lbg->AddEntry(hTDCWbg,"Calibrated TDCW","l");
  lbg->Draw();
  gPad->SetLogy();

  c5->cd(2);
//  fgchain->Draw("(-1)*(182.27/4096.0)*(TDCE - 38 - 2917) >> hTDCEfg", basicCut && fiducialCut && time2STPCut && energyCut);
//  fgchain->Draw("(-1)*(182.27/4096.0)*(TDCW + 38 - 3138) >> hTDCWfg", basicCut && fiducialCut && time2STPCut && energyCut, "SAME");
  fgchain->Draw("newTDC2TimeE >> hTDCEfg", basicCut && fiducialCut && time2STPCut && energyCut);
  fgchain->Draw("newTDC2TimeW >> hTDCWfg", basicCut && fiducialCut && time2STPCut && energyCut, "SAME");
  TLegend *lfg = new TLegend(0.1,0.8,0.4,0.9);
  lfg->AddEntry(hTDCEfg,"Calibrated TDCE","l");
  lfg->AddEntry(hTDCWfg,"Calibrated TDCW","l");
  lfg->Draw();
  gPad->SetLogy();

//  c5->Print("5_TDCShiftedTimes_fullRange.pdf");

  // seventh canvas, overlay MPM simulations on these roughly calibrated data plots
  TCanvas *c7 = new TCanvas("c7", "c7");
  c7->cd();

  TChain *mpm_sim_chain = new TChain("anaTree");
  mpm_sim_chain->Add("mpm_sim_betas_100EWFiles_2012-2013.root");

  TRandom3 *engine = new TRandom3(0);
  TH1D *hSim0 = new TH1D("hSim0sigma", "hSim", 320, -20, 140);
  hSim0->SetLineColor(6);
  TH1D *hSim2 = new TH1D("hSim2sigma", "hSim", 320, -20, 140);
  hSim2->SetLineColor(3);
  TH1D *hSim4 = new TH1D("hSim4sigma", "hSim", 320, -20, 140);
  hSim4->SetLineColor(5);

  SimEvent evt;
  mpm_sim_chain->GetBranch("Edep")->GetLeaf("EdepE")->SetAddress(&evt.Edep_EdepE);
  mpm_sim_chain->GetBranch("Edep")->GetLeaf("EdepW")->SetAddress(&evt.Edep_EdepW);
  mpm_sim_chain->GetBranch("EdepQ")->GetLeaf("EdepQE")->SetAddress(&evt.EdepQ_EdepQE);
  mpm_sim_chain->GetBranch("EdepQ")->GetLeaf("EdepQW")->SetAddress(&evt.EdepQ_EdepQW);
  mpm_sim_chain->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyE")->SetAddress(&evt.MWPCEnergy_MWPCEnergyE);
  mpm_sim_chain->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyW")->SetAddress(&evt.MWPCEnergy_MWPCEnergyW);
  mpm_sim_chain->GetBranch("time")->GetLeaf("timeE")->SetAddress(&evt.time_timeE);
  mpm_sim_chain->GetBranch("time")->GetLeaf("timeW")->SetAddress(&evt.time_timeW);

  for(unsigned int i = 0; i < 1000000 /*mpm_sim_chain->GetEntries()*/; i++)
  {
    mpm_sim_chain->GetEntry(i);

    if(evt.EdepQ_EdepQE > 0 && evt.EdepQ_EdepQW > 0 && evt.MWPCEnergy_MWPCEnergyE > 0 && evt.MWPCEnergy_MWPCEnergyW > 0
	&& evt.time_timeE < 200 && evt.time_timeW < 200 && evt.time_timeE < evt.time_timeW)
    {
      hSim0->Fill(evt.time_timeW - evt.time_timeE);
      hSim2->Fill(engine->Gaus(evt.time_timeW - evt.time_timeE, 2));
      hSim4->Fill(engine->Gaus(evt.time_timeW - evt.time_timeE, 4));
    }

    if(i % 100000 == 0)
    {
      cout << "Filling event " << i << "/" << mpm_sim_chain->GetEntries() << endl;
    }
  }

//  hTDCEfg->Draw();
//  hTDCWfg->Draw("SAME");

  TH1D *hTimeE_bgSub = new TH1D("BGSubE", "BG Subtracted Time East", 320, -20, 140);
  hTimeE_bgSub->SetLineColor(4);
  hTimeE_bgSub->Add(hTDCEfg, hTDCEbg, 1, -5.07);
  hTimeE_bgSub->Draw();

//  cout << "value of get entries for hTimeE_bgSub is " << hTimeE_bgSub->GetEntries() << endl;

  TH1D *hTimeW_bgSub = new TH1D("BGSubW", "BG Subtracted Time West", 320, -20, 140);
  hTimeW_bgSub->SetLineColor(2);
  hTimeW_bgSub->Add(hTDCWfg, hTDCWbg, 1, -5.07);
  hTimeW_bgSub->Draw("SAME");

  hSim0->Scale((double)hTimeW_bgSub->GetEntries() / (hSim0->GetEntries()));
  hSim0->Draw("SAME");
  hSim2->Scale((double)hTimeW_bgSub->GetEntries() / (hSim2->GetEntries()));
  hSim2->Draw("SAME");
  hSim4->Scale((double)hTimeW_bgSub->GetEntries() / (hSim4->GetEntries()));
//  hSim4->Draw("SAME");

  TLegend *l7 = new TLegend(0.6,0.75,0.9,0.9);
//  l7->AddEntry(hTDCEbg,"Calibrated TDCE","l");
//  l7->AddEntry(hTDCWbg,"Calibrated TDCW","l");
  l7->AddEntry(hSim0, "MPM G4 Sim no smearing", "l");
  l7->AddEntry(hSim2, "MPM G4 Sim 2ns sigma", "l");
//  l7->AddEntry(hSim4, "MPM G4 Sim 4ns sigma", "l");
  l7->AddEntry(hTimeE_bgSub, "FG - 5.07BG East", "l");
  l7->AddEntry(hTimeW_bgSub, "FG - 5.07BG West", "l");
  l7->Draw();
  gPad->SetLogy();

//  c7->Print("7_MPMSimCompare_fullRange.pdf");

  // eighth canvas
  TCanvas *c8 = new TCanvas("c8","c8");
  c8->cd();

  hTimeW_bgSub->Draw();
  hTimeE_bgSub->Draw("SAME");
  hTimeW_bgSub->SetTitle(Form("FG - 5.07BG: %f < E < %f, %f < W < %f", timeLowerEdgeE1, timeUpperEdgeE1, timeLowerEdgeW1, timeUpperEdgeW1));
  hTimeW_bgSub->GetXaxis()->SetTitle("Time (ns)");
  hTimeW_bgSub->GetXaxis()->SetRangeUser(-10, 30);
  hTimeW_bgSub->GetYaxis()->SetTitle("Counts");
  gPad->SetLogy();

  TLegend *l8 = new TLegend(0.3, 0.75, 0.6, 0.9);
  l8->AddEntry(hTimeW_bgSub, "FG - 5.07BG West", "l");
  l8->AddEntry(hTimeE_bgSub, "FG - 5.07BG East", "l");
  l8->Draw();

//  c8->Print("8_TimingSpectra_fullRange.pdf");

  // print out all the stats that we'll use
  cout << "For " << timeLowerEdgeE1 << " < E < " << timeUpperEdgeE1 << ", " << timeLowerEdgeW1 << " < W < " << timeUpperEdgeW1 << endl;
  cout << "We have full background spectrum counts: " << hbgErecon_timeWin->GetEntries() << endl;
  cout << "And full foreground spectrum counts: " << hfgErecon_timeWin->GetEntries() << endl;
  cout << "Restricted energy range background spectrum counts: " << hbgSpectra->GetEntries() << endl;
  cout << "And the corresponding foreground spectrum counts: " << hfgSpectra->GetEntries() << endl;

/*
  ofstream outfile;
  outfile.open("FinalCounts_variousTimeWindows_allForegroundRuns.txt", ios::app);

  outfile << timeLowerEdgeE << "\t"
	  << timeUpperEdgeE << "\t"
	  << timeLowerEdgeW << "\t"
	  << timeUpperEdgeW << "\t"
	  << hbgErecon_timeWin->GetEntries() << "\t"
	  << hfgErecon_timeWin->GetEntries() << "\t"
	  << hbgSpectra->GetEntries() << "\t"
	  << hfgSpectra->GetEntries() << "\n";

  outfile.close();
*/

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
