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

using            namespace std;

double SetPoissonErrors(int counts);

struct Event
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

// Used for visualization, keeps the graph on screen.
TApplication plot_program("FADC_readin",0,0,0,0);

//-------------------------------------------------//
//------------ Start of Program -------------------//
//-------------------------------------------------//

int main(int argc, char* argv[])
{
  if(argc < 5)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable) (E Low) (E High) (W Low) (W High)" << endl;
    return 0;
  }

  // read in the arguments
  double timeLowerEdgeE = atof(argv[1]);
  double timeUpperEdgeE = atof(argv[2]);
  double timeLowerEdgeW = atof(argv[3]);
  double timeUpperEdgeW = atof(argv[4]);

  TString stp = "E";
  TString roi = "W";


  TChain *bgchain = new TChain("pass3");
  TChain *fgchain = new TChain("pass3");

  bgchain->Add("TimeCalibrated_BGRuns_type1_fixed_v2.root");
  fgchain->Add("TimeCalibrated_FGRuns_type1_fixed_v2.root");

  // define all the cuts we will use later.
  TCut basicCut = "(PID == 1 && badTimeFlag == 0)";
  TCut energyCut = "(Erecon_ee > 0 && Erecon_ee < 640)";
  TCut fiducialCut = "(((xE.center)*(xE.center) + (yE.center)*(yE.center) < 49*49) && ((xW.center)*(xW.center) + (yW.center)*(yW.center) < 49*49))";

  // bad cuts down here
//  TCut scaledTimeCut = Form("((newTimeScaledW < 4 && newTimeScaledE > 4 && newTimeScaledE < %f) || (newTimeScaledE < 4 && newTimeScaledW > 4 && newTimeScaledW < %f)) && newTimeScaledW > -2 && newTimeScaledE > -2", timeWindowUpperEdge, timeWindowUpperEdge);
//  TCut shiftedTimeCut = Form("((newTimeShiftedW < 4 && newTimeShiftedE > 4 && newTimeShiftedE < %f) || (newTimeShiftedE < 4 && newTimeShiftedW > 4 && newTimeShiftedW < %f)) && newTimeShiftedW > -2 && newTimeShiftedE > -2", timeWindowUpperEdge, timeWindowUpperEdge);
//  TCut tdcTimeCut = "((TDCE > 2850 && TDCW > 2850 && TDCW < 3050) || (TDCW > 3050 && TDCE > 2650 && TDCE < 2850))";

  TCut globalTimeCut = "((newTimeGlobalShiftW < 4 && newTimeGlobalShiftE > 4 && newTimeGlobalShiftE < 16) || (newTimeGlobalShiftE < 4 && newTimeGlobalShiftW > 4 && newTimeGlobalShiftW < 16)) && newTimeGlobalShiftW > -2 && newTimeGlobalShiftE > -2";

  TCut time1STPCut = Form("newTimeGlobalShiftE < %f && newTimeGlobalShiftE > %f && newTimeGlobalShiftW > %f && newTimeGlobalShiftW < %f", timeUpperEdgeE, timeLowerEdgeE, timeLowerEdgeW, timeUpperEdgeW);

//  TCut timeWCut = Form("newTimeGlobalShiftW < %f && newTimeGlobalShiftE > %f && newTimeGlobalShiftE < %f && newTimeGlobalShiftE > -10 && newTimeGlobalShiftW > -10", timeWindowLowerEdge, timeWindowLowerEdge, timeWindowUpperEdge);

  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->Divide(2,1);
  c1->cd(1);

  TH1D *hbgErecon_timeWin = new TH1D("bgErecon_timeWin", Form("BG Erecon_ee global time: %f < %s < %f ns, %f < %s < %f ns",
				timeLowerEdgeE, stp.Data(), timeUpperEdgeE, timeLowerEdgeW, roi.Data(), timeUpperEdgeW), 160, 0, 4000);
  hbgErecon_timeWin->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  TH1D *hfgErecon_timeWin = new TH1D("fgErecon_timeWin", Form("FG Erecon_ee global time: %f < %s < %f ns, %f < %s < %f ns",
				timeLowerEdgeE, stp.Data(), timeUpperEdgeE, timeLowerEdgeW, roi.Data(), timeUpperEdgeW), 160, 0, 4000);
  hfgErecon_timeWin->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  bgchain->Draw("Erecon_ee >> bgErecon_timeWin", basicCut && fiducialCut && time1STPCut);

  c1->cd(2);
  fgchain->Draw("Erecon_ee >> fgErecon_timeWin", basicCut && fiducialCut && time1STPCut);

  c1->Print("1_Erecon_timingWindow.pdf");

  // second canvas for plots
  TCanvas *c2 = new TCanvas("c2", "c2");
  c2->Divide(3,1);

  TH1D *hbgSpectra = new TH1D("bgfull", Form("BG Erecon_ee, all cuts: %f < %s < %f, %f < %s < %f",
				timeLowerEdgeE, stp.Data(), timeUpperEdgeE, timeLowerEdgeW, roi.Data(), timeUpperEdgeW), 40, 0, 1000);
  hbgSpectra->Sumw2();
  hbgSpectra->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  TH1D *hfgSpectra = new TH1D("fgfull", Form("FG Erecon_ee, all cuts: %f < %s < %f, %f < %s < %f",
				timeLowerEdgeE, stp.Data(), timeUpperEdgeE, timeLowerEdgeW, roi.Data(), timeUpperEdgeW), 40, 0, 1000);
  hfgSpectra->Sumw2();
  hfgSpectra->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  c2->cd(1);
  bgchain->Draw("Erecon_ee >> bgfull", basicCut && fiducialCut && time1STPCut && energyCut);

  c2->cd(2);
  fgchain->Draw("Erecon_ee >> fgfull", basicCut && fiducialCut && time1STPCut && energyCut);

  // setting Poisson error bars for the BG and FG histograms
  cout << "Setting Poisson error bars..." << endl;
  for(int i = 0; i < hfgSpectra->GetNbinsX(); i++)
  {
    hbgSpectra->SetBinError(i, SetPoissonErrors(hbgSpectra->GetBinContent(i)));
    hfgSpectra->SetBinError(i, SetPoissonErrors(hfgSpectra->GetBinContent(i)));
  }

  TH1D *h_fullCuts = new TH1D("fullCuts", "BG subtracted Erecon_ee", 40, 0, 1000);
  h_fullCuts->Sumw2();
  h_fullCuts->Add(hfgSpectra, hbgSpectra, 1, -5.07);	// 5.07 comes from 860262/169717 live time ratio

  c2->cd(3);
  h_fullCuts->Draw();

  c2->Print("2_BGsubtracted_timeWindow.pdf");

  // third canvas, starting time plots
/*  TCanvas *c3 = new TCanvas("c3", "c3");
  c3->cd();
  TH2D *hTime2D = new TH2D("time2D", "time2D", 160, -20, 140, 160, -20, 140);
  hTime2D->GetXaxis()->SetTitle("newTimeGlobalShiftW (ns)");
  hTime2D->GetYaxis()->SetTitle("newTimeGlobalShiftE (ns)");
  fgchain->Draw("newTimeGlobalShiftE:newTimeGlobalShiftW >> time2D", basicCut && fiducialCut, "COLZ");
  gPad->SetLogz();
  c3->Print("3_2DTimePlots_noTimingCuts.pdf");
*/

  // fourth canvas
/*  TCanvas *c4 = new TCanvas("c4", "c4");
  c4->cd();
  TH2D *hTime2D_subrange = new TH2D("time2D_subrange", "time2D subrange", 25, -5, 20, 25, -5, 20);
  hTime2D_subrange->GetXaxis()->SetTitle("newTimeGlobalShiftW (ns)");
  hTime2D_subrange->GetYaxis()->SetTitle("newTimeGlobalShiftE (ns)");
  fgchain->Draw("newTimeGlobalShiftE:newTimeGlobalShiftW >> time2D_subrange", basicCut && fiducialCut && globalTimeCut && energyCut, "COLZ");
  gPad->SetLogz();
  c4->Print("4_2DTimePlots_withTimeCuts.pdf");
*/
  // fifth canvas
/*  TCanvas *c5 = new TCanvas("c5","c5");
  c5->Divide(3,1);

  TChain *bgchain_8ns = new TChain("pass3");
  TChain *fgchain_8ns = new TChain("pass3");

  bgchain_8ns->Add("TimeCalibrated_BGRuns_type1_fixed.root");
  fgchain_8ns->Add("TimeCalibrated_FGRuns_type1_fixed.root");

  TH1D *hbg_8ns_scaled = new TH1D("bg_8ns_scaled", Form("Background %f ns window, scaled to counts of full window", modelTimeUpperEdge - timeWindowLowerEdge), 160, 0, 4000);
  hbgErecon_timeWin->GetXaxis()->SetTitle("Erecon_ee (KeV)");
  TH1D *hfg_8ns_scaled = new TH1D("fg_8ns_scaled", Form("Foreground %f ns window, scaled to counts of full window", modelTimeUpperEdge - timeWindowLowerEdge), 160, 0, 4000);
  hfgErecon_timeWin->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  TH1D *hfg_realWindow = new TH1D("fgreal", "Foreground of full window", 160, 0, 4000);
  hfgErecon_timeWin->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  c5->cd(1);
  bgchain_8ns->Draw("Erecon_ee >> bg_8ns_scaled", basicCut && fiducialCut && scaledTimeCut_8ns);

  c5->cd(2);
  hfg_realWindow->SetFillColor(38);
  hfg_realWindow->SetFillStyle(3005);
  fgchain->Draw("Erecon_ee >> fgreal", basicCut && fiducialCut && scaledTimeCut);

  c5->cd(3);
  fgchain_8ns->Draw("Erecon_ee >> fg_8ns_scaled", basicCut && fiducialCut && scaledTimeCut_8ns);

  c5->cd(2);

  double fullwindow_bgCounts = hbgErecon_timeWin->GetEntries();
  double window8ns_bgCounts = hbg_8ns_scaled->GetEntries();
  hfg_8ns_scaled->Scale(fullwindow_bgCounts / window8ns_bgCounts);
  hfg_8ns_scaled->SetFillColor(46);
  hfg_8ns_scaled->SetFillStyle(3004);
  hfg_8ns_scaled->Draw("SAME");


  c5->cd(3);
//  hbgErecon_timeWin->SetFillColor(30);
//  hbgErecon_timeWin->Draw();
  hfgErecon_timeWin->SetFillColor(38);
  hfgErecon_timeWin->Draw();
*/
  // sixth canvas
  TCanvas *c6 = new TCanvas("c6", "c6");
  c6->Divide(2,1);

  TH1D *hTDCEfg = new TH1D("hTDCEfg", "FG TDC Spectra", 4000, 0, 4000);
  hTDCEfg->GetXaxis()->SetTitle("Channels");
  hTDCEfg->SetLineColor(2);

  TH1D *hTDCWfg = new TH1D("hTDCWfg", "FG TDCW", 4000, 0, 4000);
  hTDCWfg->GetXaxis()->SetTitle("Channels");
  hTDCWfg->SetLineColor(4);

  TH1D *hTDCEbg = new TH1D("hTDCEbg", "BG TDC Spectra", 4000, 0, 4000);
  hTDCEbg->GetXaxis()->SetTitle("Channels");
  hTDCEbg->SetLineColor(2);

  TH1D *hTDCWbg = new TH1D("hTDCWbg", "BG TDCW", 4000, 0, 4000);
  hTDCWbg->GetXaxis()->SetTitle("Channels");
  hTDCWbg->SetLineColor(4);

  c6->cd(1);
  bgchain->Draw("TDCE + 150 >> hTDCEbg", basicCut && fiducialCut /*&& time1STPCut && energyCut*/);
  bgchain->Draw("TDCW >> hTDCWbg", basicCut && fiducialCut /*&& time1STPCut && energyCut*/, "SAME");
  TLegend *lbg = new TLegend(0.1,0.7,0.48,0.9);
  lbg->AddEntry(hTDCEbg,"TDCE + 150","l");
  lbg->AddEntry(hTDCWbg,"TDCW","l");
  lbg->Draw();
  gPad->SetLogy();

  c6->cd(2);
  fgchain->Draw("TDCE + 150 >> hTDCEfg", basicCut && fiducialCut /*&& time1STPCut && energyCut*/);
  fgchain->Draw("TDCW >> hTDCWfg", basicCut && fiducialCut /*&& time1STPCut && energyCut*/, "SAME");
  TLegend *lfg = new TLegend(0.1,0.7,0.48,0.9);
  lfg->AddEntry(hTDCEfg,"TDCE + 150","l");
  lfg->AddEntry(hTDCWfg,"TDCW","l");
  lfg->Draw();
  gPad->SetLogy();


  c6->Print("6_globalShiftedTimes_fullRange.pdf");


  // print out all the stats that we'll use
  cout << "For " << timeLowerEdgeE << " < E < " << timeUpperEdgeE << ", " << timeLowerEdgeW << " < W < " << timeUpperEdgeW << endl;
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
