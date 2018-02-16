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
    upperErrBar = sqrt(counts);
  }

  return upperErrBar;
}

bgfgGraphs()
{
  double timeWindowUpperEdge = 18;
  double timeWindowLowerEdge = 4;

  double modelTimeUpperEdge = 12;

  TChain *bgchain = new TChain("pass3");
  TChain *fgchain = new TChain("pass3");

  bgchain->Add("TimeCalibrated_BGRuns_type1_fixed.root");
  fgchain->Add("TimeCalibrated_FGRuns_type1_fixed.root");

  // define all the cuts we will use later.
  TCut basicCut = "(PID == 1 && badTimeFlag == 0)";
  TCut energyCut = "(Erecon_ee > 0 && Erecon_ee < 640)";
  // assuming TDCE self-timing peak at 2850 channels, TDCW at 3050
/*  TCut oneSTPeak_100channels = "((TDCE > 2850 && TDCW > 2950 && TDCW < 3050) || (TDCW > 3050 && TDCE > 2750 && TDCE < 2850))";
  TCut oneSTPeak_200channels = "((TDCE > 2850 && TDCW > 2850 && TDCW < 3050) || (TDCW > 3050 && TDCE > 2650 && TDCE < 2850))";
  TCut oneSTPeak_300channels = "((TDCE > 2850 && TDCW > 2750 && TDCW < 3050) || (TDCW > 3050 && TDCE > 2550 && TDCE < 2850))";
  TCut oneSTPeak_400channels = "((TDCE > 2850 && TDCW > 2650 && TDCW < 3050) || (TDCW > 3050 && TDCE > 2450 && TDCW < 2850))";
  TCut oneSTPeak_500channels = "((TDCE > 2850 && TDCW > 2550 && TDCW < 3050) || (TDCW > 3050 && TDCE > 2350 && TDCW < 2850))";
*/  TCut fiducialCut = "(((xE.center)*(xE.center) + (yE.center)*(yE.center) < 49*49) && ((xW.center)*(xW.center) + (yW.center)*(yW.center) < 49*49))";
  TCut scaledTimeCut = Form("((newTimeScaledW < 4 && newTimeScaledE > 4 && newTimeScaledE < %f) || (newTimeScaledE < 4 && newTimeScaledW > 4 && newTimeScaledW < %f)) && newTimeScaledW > -2 && newTimeScaledE > -2", timeWindowUpperEdge, timeWindowUpperEdge);
  TCut scaledTimeCut_8ns = Form("((newTimeScaledW < 4 && newTimeScaledE > 4 && newTimeScaledE < %f) || (newTimeScaledE < 4 && newTimeScaledW > 4 && newTimeScaledW < %f)) && newTimeScaledW > -2 && newTimeScaledE > -2", modelTimeUpperEdge, modelTimeUpperEdge);
  TCut shiftedTimeCut = Form("((newTimeShiftedW < 4 && newTimeShiftedE > 4 && newTimeShiftedE < %f) || (newTimeShiftedE < 4 && newTimeShiftedW > 4 && newTimeShiftedW < %f)) && newTimeShiftedW > -2 && newTimeShiftedE > -2", timeWindowUpperEdge, timeWindowUpperEdge);
  TCut tdcTimeCut = "((TDCE > 2850 && TDCW > 2850 && TDCW < 3050) || (TDCW > 3050 && TDCE > 2650 && TDCE < 2850))";

  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->Divide(2,1);
  c1->cd(1);

  TH1D *hbgErecon_timeWin = new TH1D("bgErecon_timeWin", Form("BG Erecon_ee scaled time %f ns, %f ns low edge", timeWindowUpperEdge - timeWindowLowerEdge, timeWindowLowerEdge), 160, 0, 4000);
  hbgErecon_timeWin->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  TH1D *hfgErecon_timeWin = new TH1D("fgErecon_timeWin", Form("FG Erecon_ee scaled time %f ns, %f ns low edge", timeWindowUpperEdge - timeWindowLowerEdge, timeWindowLowerEdge), 160, 0, 4000);
  hfgErecon_timeWin->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  bgchain->Draw("Erecon_ee >> bgErecon_timeWin", basicCut && fiducialCut && scaledTimeCut);

  c1->cd(2);
  fgchain->Draw("Erecon_ee >> fgErecon_timeWin", basicCut && fiducialCut && scaledTimeCut);

  c1->Print("1_Erecon_timingWindow.pdf");

  // fourth canvas for plots
  TCanvas *c2 = new TCanvas("c2", "c2");
  c2->Divide(3,1);

  TH1D *hbgSpectra = new TH1D("bgfull", Form("BG Erecon_ee, all cuts, scaled time %f ns, %f ns low edge", timeWindowUpperEdge - timeWindowLowerEdge, timeWindowLowerEdge), 40, 0, 1000);
  hbgSpectra->Sumw2();
  hbgSpectra->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  TH1D *hfgSpectra = new TH1D("fgfull", Form("FG Erecon_ee, all cuts, scaled time %f ns, %f ns low edge", timeWindowUpperEdge - timeWindowLowerEdge, timeWindowLowerEdge), 40, 0, 1000);
  hfgSpectra->Sumw2();
  hfgSpectra->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  c2->cd(1);
  bgchain->Draw("Erecon_ee >> bgfull", basicCut && fiducialCut && scaledTimeCut && energyCut);

  c2->cd(2);
  fgchain->Draw("Erecon_ee >> fgfull", basicCut && fiducialCut && scaledTimeCut && energyCut);

  // setting Poisson error bars for the BG and FG histograms
/*  cout << "Setting Poisson error bars..." << endl;
  for(int i = 0; i < hfgSpectra->GetNbinsX(); i++)
  {
    hbgSpectra->SetBinError(i, SetPoissonErrors(hbgSpectra->GetBinContent(i)));
    hfgSpectra->SetBinError(i, SetPoissonErrors(hfgSpectra->GetBinContent(i)));
  }
*/
  TH1D *h_fullCuts = new TH1D("fullCuts", "BG subtracted Erecon_ee", 40, 0, 1000);
  h_fullCuts->Sumw2();
  h_fullCuts->Add(hfgSpectra, hbgSpectra, 1, -5.07);	// 5.07 comes from 860262/169717 live time ratio

/*
  for(int i = 0; i < hfgSpectra->GetNbinsX(); i++)
  {
    cout << "BG Hist bin " << i << " has value " << hbgSpectra->GetBinContent(i) << " and error " << hbgSpectra->GetBinError(i) << endl;
    cout << "FG Hist bin " << i << " has value " << hfgSpectra->GetBinContent(i) << " and error " << hfgSpectra->GetBinError(i) << endl;
    cout << "Sum Hist bin " << i << " has value " << h_fullCuts->GetBinContent(i) << " and error " << h_fullCuts->GetBinError(i) << endl;
  }
*/

  c2->cd(3);
  h_fullCuts->Draw();

  // fit this with a flat line to get a measure of the beta contamination
/*  double backgroundFitMin = 200;
  double backgroundFitMax = 640;
  TF1 *fBetas = new TF1("beta contamination", "[0]", backgroundFitMin, backgroundFitMax);
  fBetas->FixParameter(0, 0);
  h_fullCuts->Fit(fBetas, "RP", "", backgroundFitMin, backgroundFitMax);
  TF1 *fBetaFit = h_fullCuts->GetFunction("beta contamination");
  double betaFit = fBetaFit->GetParameter(0);
  double betaFitErr = fBetaFit->GetParError(0);
  double chi2 = fBetaFit->GetChisquare();
  double ndf = fBetaFit->GetNDF();

  cout << "Flat-line fit gives value " << betaFit << " with error " << betaFitErr << " and chi-squared " << chi2 << " with " << ndf << " degrees of freedom." << endl;
*/
  c2->Print("2_BGsubtracted_timeWindow.pdf");

  // third canvas, starting time plots
  TCanvas *c3 = new TCanvas("c3", "c3");
  c3->cd();
  TH2D *hTime2D = new TH2D("time2D", "time2D", 160, -20, 140, 160, -20, 140);
  hTime2D->GetXaxis()->SetTitle("newTimeScaledW (ns)");
  hTime2D->GetYaxis()->SetTitle("newTimeScaledE (ns)");
  fgchain->Draw("newTimeScaledE:newTimeScaledW >> time2D", basicCut && fiducialCut, "COLZ");
  gPad->SetLogz();
  c3->Print("3_2DTimePlots_noTimingCuts.pdf");


  // fourth canvas
  TCanvas *c4 = new TCanvas("c4", "c4");
  c4->cd();
  TH2D *hTime2D_subrange = new TH2D("time2D_subrange", "time2D subrange", 25, -5, 20, 25, -5, 20);
  hTime2D_subrange->GetXaxis()->SetTitle("newTimeScaledW (ns)");
  hTime2D_subrange->GetXaxis()->SetTitle("newTimeScaledE (ns)");
  fgchain->Draw("newTimeScaledE:newTimeScaledW >> time2D_subrange", basicCut && fiducialCut && scaledTimeCut && energyCut, "COLZ");
  gPad->SetLogz();
  c4->Print("4_2DTimePlots_withTimeCuts.pdf");

  // fifth canvas
  TCanvas *c5 = new TCanvas("c5","c5");
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


/*  TCanvas *c5 = new TCanvas("c5", "c5");
  c5->Divide(2,1);

  TH1D *hnewTimeE_bg = new TH1D("htimeEbg", "BG TDCE", 4000, 0, 4000);
  hnewTimeE_bg->GetXaxis()->SetTitle("Channel");
  TH1D *hnewTimeW_bg = new TH1D("htimeWbg", "BG TDCW", 4000, 0, 4000);
  hnewTimeW_bg->GetXaxis()->SetTitle("Channel");

  c5->cd(1);
  bgchain->Draw("TDCE >> htimeEbg", basicCut && Erecon_ee_range && fiducialCut);
  gPad->SetLogy();

  c5->cd(2);
  bgchain->Draw("TDCW >> htimeWbg", basicCut && Erecon_ee_range && fiducialCut);
  gPad->SetLogy();

  c5->Print("3_ScaledTimeBG_runbyrun.pdf");

  // sixth canvas for time FG plots
  TCanvas *c6 = new TCanvas("c6", "c6");
  c6->Divide(2,1);

  TH1D *hnewTimeE_fg = new TH1D("htimeEfg", "FG TDCE", 4000, 0, 4000);
  hnewTimeE_fg->GetXaxis()->SetTitle("Channel");
  TH1D *hnewTimeW_fg = new TH1D("htimeWfg", "FG TDCW", 4000, 0, 4000);
  hnewTimeW_fg->GetXaxis()->SetTitle("Channel");

  c6->cd(1);
  fgchain->Draw("TDCE >> htimeEfg", basicCut && Erecon_ee_range && fiducialCut);
  gPad->SetLogy();

  c6->cd(2);
  fgchain->Draw("TDCW >> htimeWfg", basicCut && Erecon_ee_range && fiducialCut);
  gPad->SetLogy();

  c6->Print("4_ScaledTimeFG_runbyrun.pdf");
*/

  cout << "End of program." << endl;

}
