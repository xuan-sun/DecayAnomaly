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
  TCut scaledTimeCut = "((newTimeScaledW < 4 && newTimeScaledE > 4 && newTimeScaledE < 17) || (newTimeScaledE < 4 && newTimeScaledW > 4 && newTimeScaledW < 17)) && newTimeScaledW > -2 && newTimeScaledE > -2";
  TCut shiftedTimeCut = "((newTimeShiftedW < 4 && newTimeShiftedE > 4 && newTimeShiftedE < 17) || (newTimeShiftedE < 4 && newTimeShiftedW > 4 && newTimeShiftedW < 17)) && newTimeShiftedW > -2 && newTimeShiftedE > -2";
  TCut tdcTimeCut = "((TDCE > 2850 && TDCW > 2715 && TDCW < 3050) || (TDCW > 3050 && TDCE > 2515 && TDCE < 2850))";

  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->Divide(2,1);
  c1->cd(1);

  TH1D *hbgErecon_timeWin = new TH1D("bgErecon_timeWin", "BG Erecon_ee scaled time 13ns, 4ns low edge", 160, 0, 4000);
  hbgErecon_timeWin->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  TH1D *hfgErecon_timeWin = new TH1D("fgErecon_timeWin", "FG Erecon_ee scaled time 13ns, 4ns low edge", 160, 0, 4000);
  hfgErecon_timeWin->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  bgchain->Draw("Erecon_ee >> bgErecon_timeWin", basicCut && fiducialCut && scaledTimeCut);

  c1->cd(2);
  fgchain->Draw("Erecon_ee >> fgErecon_timeWin", basicCut && fiducialCut && scaledTimeCut);

  c1->Print("1_Erecon_timingWindow.pdf");

  // fourth canvas for plots
  TCanvas *c2 = new TCanvas("c2", "c2");
  c2->Divide(3,1);

  TH1D *hbgSpectra = new TH1D("bgfull", "BG full Erecon_ee spectra with all cuts, scaled time 13ns, 4ns low edge", 40, 0, 1000);
  hbgSpectra->Sumw2();
  hbgSpectra->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  TH1D *hfgSpectra = new TH1D("fgfull", "FG full Erecon_ee spectra with all cuts, scaled time 13ns, 4ns low edge", 40, 0, 1000);
  hfgSpectra->Sumw2();
  hfgSpectra->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  c2->cd(1);
  bgchain->Draw("Erecon_ee >> bgfull", basicCut && fiducialCut && scaledTimeCut && energyCut);

  c2->cd(2);
  fgchain->Draw("Erecon_ee >> fgfull", basicCut && fiducialCut && scaledTimeCut && energyCut);

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
  c4->Print("4_2DTimePlots_13nsTimeCuts.pdf");


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
