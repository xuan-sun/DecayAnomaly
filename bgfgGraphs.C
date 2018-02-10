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

  bgchain->Add("TimeCalibrated_BGRuns_type1.root");
  fgchain->Add("TimeCalibrated_FGRuns_type1.root");


  // define all the cuts we will use later.
  TCut basicCut = "(PID == 1 && badTimeFlag == 0)";
  TCut Erecon_ee_range = "(Erecon_ee > 0 && Erecon_ee < 640)";
  // assuming TDCE self-timing peak at 2850 channels, TDCW at 3050
  TCut oneSTPeak_100channels = "((TDCE > 2850 && TDCW > 2950 && TDCW < 3050) || (TDCW > 3050 && TDCE > 2750 && TDCE < 2850))";
  TCut oneSTPeak_200channels = "((TDCE > 2850 && TDCW > 2850 && TDCW < 3050) || (TDCW > 3050 && TDCE > 2650 && TDCE < 2850))";
  TCut oneSTPeak_300channels = "((TDCE > 2850 && TDCW > 2750 && TDCW < 3050) || (TDCW > 3050 && TDCE > 2550 && TDCE < 2850))";
  TCut oneSTPeak_400channels = "((TDCE > 2850 && TDCW > 2650 && TDCW < 3050) || (TDCW > 3050 && TDCE > 2450 && TDCW < 2850))";
  TCut oneSTPeak_500channels = "((TDCE > 2850 && TDCW > 2550 && TDCW < 3050) || (TDCW > 3050 && TDCE > 2350 && TDCW < 2850))";
  TCut radialCutE = "((xE.center)*(xE.center) + (yE.center)*(yE.center) < 49*49)";
  TCut radialCutW = "((xW.center)*(xW.center) + (yW.center)*(yW.center) < 49*49)";
  TCut fiducialCut = "(((xE.center)*(xE.center) + (yE.center)*(yE.center) < 49*49) && ((xW.center)*(xW.center) + (yW.center)*(yW.center) < 49*49))";
  TCut tennsTimeCut = "((newTimeW < 4 && newTimeE > 4 && newTimeE < 14) || (newTimeE < 4 && newTimeW > 4 && newTimeW < 14))";

  // first canvas for plots
  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->Divide(2,1);
  c1->cd(1);

  TH1D* hbgErecon = new TH1D("BGErecon", "BG Erecon_ee", 160, 0, 4000);
  hbgErecon->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  TH1D* hfgErecon = new TH1D("FGErecon", "FG Erecon_ee", 160, 0, 4000);
  hfgErecon->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  bgchain->Draw("Erecon_ee", basicCut);

  c1->cd(2);
  fgchain->Draw("Erecon_ee", basicCut);

  c1->Print("1_Erecon_basicCut.pdf");

  // second canvas for plots
  TCanvas *c2 = new TCanvas("c2", "c2");
  c2->Divide(2,1);
  c2->cd(1);

  TH1D *hbgEreconFid = new TH1D("BGEreconFid", "BG Erecon_ee + Fiducial Cut", 160, 0, 4000);
  hbgEreconFid->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  TH1D *hfgEreconFid = new TH1D("FGEreconFid", "FG Erecon_ee + Fiducial Cut", 160, 0, 4000);
  hfgEreconFid->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  bgchain->Draw("Erecon_ee", basicCut && fiducialCut);

  c2->cd(2);
  fgchain->Draw("Erecon_ee", basicCut && fiducialCut);

  c2->Print("2_Erecon_basic+radialCut.pdf");

  // third canvas for plots
  TCanvas *c3 = new TCanvas("c3", "c3");
  c3->Divide(2,1);
  c3->cd(1);

  TH1D *hbgErecon_timeWin = new TH1D("bgErecon_timeWin", "BG Erecon_ee 10ns window", 160, 0, 4000);
  hbgErecon_timeWin->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  TH1D *hfgErecon_timeWin = new TH1D("fgErecon_timeWin", "FG Erecon_ee 10ns window", 160, 0, 4000);
  hfgErecon_timeWin->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  bgchain->Draw("Erecon_ee >> bgErecon_timeWin", basicCut && fiducialCut && oneSTPeak_200channels);

  c3->cd(2);
  fgchain->Draw("Erecon_ee >> fgErecon_timeWin", basicCut && fiducialCut && oneSTPeak_200channels);

  c3->Print("3_Erecon_timeWindow.pdf");

  // fourth canvas for plots
  TCanvas *c4 = new TCanvas("c4", "c4");
  c4->Divide(3,1);

  TH1D *hbgSpectra = new TH1D("bgfull", "BG full Erecon_ee spectra with all cuts, 10ns", 40, 0, 1000);
  hbgSpectra->Sumw2();
  hbgSpectra->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  TH1D *hfgSpectra = new TH1D("fgfull", "FG full Erecon_ee spectra with all cuts, 10ns", 40, 0, 1000);
  hfgSpectra->Sumw2();
  hfgSpectra->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  c4->cd(1);
  bgchain->Draw("Erecon_ee >> bgfull", basicCut && fiducialCut && oneSTPeak_200channels && Erecon_ee_range);

  c4->cd(2);
  fgchain->Draw("Erecon_ee >> fgfull", basicCut && fiducialCut && oneSTPeak_200channels && Erecon_ee_range);

  // setting Poisson error bars for the BG and FG histograms
  cout << "Setting Poisson error bars..." << endl;
  for(int i = 0; i < hfgSpectra->GetNbinsX(); i++)
  {
    hbgSpectra->SetBinError(i, SetPoissonErrors(hbgSpectra->GetBinContent(i)));
    hfgSpectra->SetBinError(i, SetPoissonErrors(hfgSpectra->GetBinContent(i)));
  }

  TH1D *h_fullCuts = new TH1D("fullCuts", "BG subtracted Erecon_ee, 10ns window", 40, 0, 1000);
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

  c4->cd(3);
  h_fullCuts->Draw();

  c4->Print("4_BGsubtracted_timeWindow.pdf");

  // fifth canvas for time BG plots
  TCanvas *c5 = new TCanvas("c5", "c5");
  c5->Divide(3,1);

  TH1D *hnewTimeE_bg = new TH1D("htimeEbg", "BG newTimeE", 1600, -20, 140);
  hnewTimeE_bg->GetXaxis()->SetTitle("Time of Flight (ns)");
  TH1D *hnewTimeW_bg = new TH1D("htimeWbg", "BG newTimeW", 1600, -20, 140);
  hnewTimeW_bg->GetXaxis()->SetTitle("Time of Flight (ns)");

  c5->cd(1);
  bgchain->Draw("newTimeE >> htimeEbg", basicCut && Erecon_ee_range && fiducialCut);
  gPad->SetLogy();

  c5->cd(2);
  bgchain->Draw("newTimeW >> htimeWbg", basicCut && Erecon_ee_range && fiducialCut);
  gPad->SetLogy();

  c5->cd(3);


  // sixth canvas for time FG plots
  TCanvas *c6 = new TCanvas("c6", "c6");
  c6->Divide(3,1);

  TH1D *hnewTimeE_fg = new TH1D("htimeEfg", "FG newTimeE", 1600, -20, 140);
  hnewTimeE_fg->GetXaxis()->SetTitle("Time of Flight (ns)");
  TH1D *hnewTimeW_fg = new TH1D("htimeWfg", "FG newTimeW", 1600, -20, 140);
  hnewTimeW_fg->GetXaxis()->SetTitle("Time of Flight (ns)");

  c6->cd(1);
  fgchain->Draw("newTimeE >> htimeEfg", basicCut && Erecon_ee_range && fiducialCut);
  gPad->SetLogy();

  c6->cd(2);
  fgchain->Draw("newTimeW >> htimeWfg", basicCut && Erecon_ee_range && fiducialCut);
  gPad->SetLogy();



  cout << "End of program." << endl;

}
