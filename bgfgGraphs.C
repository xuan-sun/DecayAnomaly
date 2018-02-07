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

double east_channelToTime = 140.0 / 2850.0;
double west_channelToTime = 140.0 / 3050.0;

bgfgGraphs()
{
  TChain *bgchain = new TChain("pass3");
  TChain *fgchain = new TChain("pass3");

  bgchain->Add("AllBGRuns_Octets_80-120_type1.root");
  fgchain->Add("AllFGRuns_Octets_80-120_type1.root");



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

  // first canvas for plots
  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->Divide(2,1);
  c1->cd(1);

  TH1D* hbgErecon = new TH1D("BGErecon", "BG Erecon_ee", 400, 0, 4000);
  hbgErecon->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  TH1D* hfgErecon = new TH1D("FGErecon", "FG Erecon_ee", 400, 0, 4000);
  hfgErecon->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  bgchain->Draw("Erecon_ee", basicCut);

  c1->cd(2);
  fgchain->Draw("Erecon_ee", basicCut);

  c1->Print("1_Erecon_basicCut.pdf");

  // second canvas for plots
  TCanvas *c2 = new TCanvas("c2", "c2");
  c2->Divide(2,1);
  c2->cd(1);

  TH1D *hbgEreconFid = new TH1D("BGEreconFid", "BG Erecon_ee + Fiducial Cut", 400, 0, 4000);
  hbgEreconFid->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  TH1D *hfgEreconFid = new TH1D("FGEreconFid", "FG Erecon_ee + Fiducial Cut", 400, 0, 4000);
  hfgEreconFid->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  bgchain->Draw("Erecon_ee", basicCut && radialCutE && radialCutW);

  c2->cd(2);
  fgchain->Draw("Erecon_ee", basicCut && radialCutE && radialCutW);

  c2->Print("2_Erecon_basic+radialCut.pdf");

  // third canvas for plots
  TCanvas *c3 = new TCanvas("c3", "c3");
  c3->Divide(2,1);
  c3->cd(1);

  TH1D *hbgErecon_200chan = new TH1D("bgErecon_200chan", "BG Erecon_ee 200 channel window", 400, 0, 4000);
  hbgErecon_200chan->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  TH1D *hfgErecon_200chan = new TH1D("fgErecon_200chan", "FG Erecon_ee 200 channel window", 400, 0, 4000);
  hfgErecon_200chan->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  bgchain->Draw("Erecon_ee >> bgErecon_200chan", basicCut && radialCutE && radialCutW && oneSTPeak_200channels);

  c3->cd(2);
  fgchain->Draw("Erecon_ee >> fgErecon_200chan", basicCut && radialCutE && radialCutW && oneSTPeak_200channels);

  c3->Print("3_Erecon_200chanWindow.pdf");

  // fourth canvas for plots
  TCanvas *c4 = new TCanvas("c4", "c4");
  c4->Divide(3,1);

  TH1D *hbgSpectra = new TH1D("bgfull", "BG full Erecon_ee spectra with all cuts, 200 chan", 100, 0, 1000);
  hbgSpectra->Sumw2();
  hbgSpectra->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  TH1D *hfgSpectra = new TH1D("fgfull", "FG full Erecon_ee spectra with all cuts, 200 chan", 100, 0, 1000);
  hfgSpectra->Sumw2();
  hfgSpectra->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  c4->cd(1);
  bgchain->Draw("Erecon_ee >> bgfull", basicCut && radialCutE && radialCutW && oneSTPeak_200channels && Erecon_ee_range);

  c4->cd(2);
  fgchain->Draw("Erecon_ee >> fgfull", basicCut && radialCutE && radialCutW && oneSTPeak_200channels && Erecon_ee_range);


  TH1D *h200chan_full = new TH1D("200chanfull", "BG subtracted Erecon_ee, 200 channel window", 100, 0, 1000);
  h200chan_full->Sumw2();
  h200chan_full->Add(hfgSpectra, hbgSpectra, 1, -5.07);	// 5.07 comes from 860262/169717 live time ratio

  c4->cd(3);
  h200chan_full->Draw();

  c4->Print("4_BGsubtracted_200chanWindow.pdf");

}
