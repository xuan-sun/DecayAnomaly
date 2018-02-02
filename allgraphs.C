#define	path	"Type1_Octets_80-120"

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

allgraphs()
{
  cout << "Chaining together all our TFiles... ";

  // chain together all our octets of interest
  TChain *chain = new TChain("pass3");
  chain->Add(Form("%s/Octet_80_type1.root", path));
  chain->Add(Form("%s/Octet_81_type1.root", path));
  chain->Add(Form("%s/Octet_82_type1.root", path));
  chain->Add(Form("%s/Octet_83_type1.root", path));
  chain->Add(Form("%s/Octet_84_type1.root", path));
  chain->Add(Form("%s/Octet_85_type1.root", path));
  chain->Add(Form("%s/Octet_86_type1.root", path));
  chain->Add(Form("%s/Octet_87_type1.root", path));
  chain->Add(Form("%s/Octet_88_type1.root", path));
  chain->Add(Form("%s/Octet_89_type1.root", path));
  chain->Add(Form("%s/Octet_90_type1.root", path));
  chain->Add(Form("%s/Octet_92_type1.root", path));
  chain->Add(Form("%s/Octet_94_type1.root", path));
  chain->Add(Form("%s/Octet_95_type1.root", path));
  chain->Add(Form("%s/Octet_96_type1.root", path));
  chain->Add(Form("%s/Octet_97_type1.root", path));
  chain->Add(Form("%s/Octet_98_type1.root", path));
  chain->Add(Form("%s/Octet_99_type1.root", path));
  chain->Add(Form("%s/Octet_100_type1.root", path));
  chain->Add(Form("%s/Octet_102_type1.root", path));
  chain->Add(Form("%s/Octet_103_type1.root", path));
  chain->Add(Form("%s/Octet_104_type1.root", path));
  chain->Add(Form("%s/Octet_105_type1.root", path));
  chain->Add(Form("%s/Octet_106_type1.root", path));
  chain->Add(Form("%s/Octet_108_type1.root", path));
  chain->Add(Form("%s/Octet_109_type1.root", path));
  chain->Add(Form("%s/Octet_110_type1.root", path));
  chain->Add(Form("%s/Octet_111_type1.root", path));
  chain->Add(Form("%s/Octet_112_type1.root", path));
  chain->Add(Form("%s/Octet_113_type1.root", path));
  chain->Add(Form("%s/Octet_114_type1.root", path));
  chain->Add(Form("%s/Octet_115_type1.root", path));
  chain->Add(Form("%s/Octet_116_type1.root", path));
  chain->Add(Form("%s/Octet_117_type1.root", path));
  chain->Add(Form("%s/Octet_118_type1.root", path));
  chain->Add(Form("%s/Octet_119_type1.root", path));
  chain->Add(Form("%s/Octet_120_type1.root", path));

  cout << "Done." << endl;

  // define all the cuts we will use later.
  TCut basicCut = "(PID == 1 && badTimeFlag == 0)";
  TCut Erecon_ee_range = "(Erecon_ee > 0 && Erecon_ee < 640)";
  // assuming TDCE self-timing peak at 2850 channels, TDCW at 3050
  TCut oneSTPeak_100channels = "((TDCE > 2850 && TDCW > 2950 && TDCW < 3050) || (TDCW > 3050 && TDCE > 2750 && TDCE < 2850))";
  TCut oneSTPeak_200channels = "((TDCE > 2850 && TDCW > 2850 && TDCW < 3050) || (TDCW > 3050 && TDCE > 2650 && TDCE < 2850))";
  TCut oneSTPeak_300channels = "((TDCE > 2850 && TDCW > 2750 && TDCW < 3050) || (TDCW > 3050 && TDCE > 2550 && TDCE < 2850))";

  cout << "Making our first canvas' plots of TDC East and West values... ";

  // make plot of the TDCE and TDCW full range values
  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->Divide(2,1);
  c1->cd(1);

  // create histograms (our choice of binning) for TDCE, TDCW totals.
  TH1D *h1tdce = new TH1D("h1tdce", "h1 tdce", 800, 0, 4000);
  TH1D *h1tdcw = new TH1D("h1tdcw", "h1 tdcw", 800, 0, 4000);

  chain->Draw("TDCE >> h1tdce", basicCut);
  gPad->SetLogy();

  c1->cd(2);
  chain->Draw("TDCW >> h1tdcw", basicCut);
  gPad->SetLogy();

  cout << "Done." << endl;

  cout << "Making second canvas with subrange TDC plots... ";

  // make second plot of TDCE/TDCW values, sub range to look at "gap"
  TCanvas *c2 = new TCanvas("c2", "c2");
  c2->Divide(2,1);
  c2->cd(1);

  TH1D *h1tdce_subrange = new TH1D("h1tdce_subrange", "h1 tdce subrange", 160, 2400, 3200);
  TH1D *h1tdcw_subrange = new TH1D("h1tdcw_subrange", "h1 tdcw subrange", 160, 2400, 3200);

  h1tdce_subrange->GetYaxis()->SetRangeUser(1, 5000);
  chain->Draw("TDCE >> h1tdce_subrange", basicCut);

  c2->cd(2);
  h1tdcw_subrange->GetYaxis()->SetRangeUser(1, 5000);
  chain->Draw("TDCW >> h1tdce_subrange", basicCut);

  cout << "Done." << endl;

  cout << "Creating a third canvas with a true time conversion using the self-timing peak... " << endl;

  // convert both East and West TDC values into "time" using crude self-timing peak to 140ns conversion
  TCanvas *c3 = new TCanvas("c3", "c3");
  c3->Divide(3,1);
  c3->cd(1);

  TH1D *h2timeE = new TH1D("timeE", "time E", 180, 0, 180);
  h2timeE->GetXaxis()->SetTitle("Time (ns)");
  h2timeE->GetYaxis()->SetTitle("Counts");
  TH1D *h2timeW = new TH1D("timeW", "time W", 180, 0, 180);
  h2timeW->GetXaxis()->SetTitle("Time (ns)");
  h2timeW->GetYaxis()->SetTitle("Counts");

  TH1D *h2totalTime = new TH1D("total time", "total time", 180, 0, 180);
  h2totalTime->GetXaxis()->SetTitle("Time (ns)");
  h2totalTime->GetYaxis()->SetTitle("Counts");

  Event evt;
  chain->SetBranchAddress("TDCE", &evt.TDCE);
  chain->SetBranchAddress("TDCW", &evt.TDCW);
  chain->SetBranchAddress("PID", &evt.PID);
  chain->SetBranchAddress("Type", &evt.Type);
  chain->SetBranchAddress("Side", &evt.Side);
  chain->SetBranchAddress("Erecon", &evt.Erecon);
  chain->SetBranchAddress("Erecon_ee", &evt.Erecon_ee);
  chain->SetBranchAddress("badTimeFlag", &evt.badTimeFlag);


  for(unsigned int i = 0; i < chain->GetEntries(); i++)
  {
    chain->GetEntry(i);
    if(evt.PID == 1 && evt.badTimeFlag == 0)
    {
      h2timeE->Fill(evt.TDCE * 140.0 / 2850.0);	// converts to ns
      h2timeW->Fill(evt.TDCW * 140.0 / 3050.0);
    }

    if(evt.PID == 1 && evt.badTimeFlag == 0 && (evt.TDCE < 2850 && evt.TDCW > 3050))
    {
      h2totalTime->Fill(evt.TDCE * 140.0 / 2850.0);
    }
    if(evt.PID == 1 && evt.badTimeFlag == 0 && (evt.TDCE > 2850 && evt.TDCW < 3050))
    {
      h2totalTime->Fill(evt.TDCW * 140.0 / 3050.0);
    }


    if(i%100000 == 0)
    {
      cout << "Completed loading " << i << " events out of " << chain->GetEntries() << endl;
    }
  }

  h2timeE->Draw();

  c3->cd(2);
  h2timeW->Draw();

  c3->cd(3);
  h2totalTime->Draw();

  cout << "Creating a fourth canvas for Erecon_ee plots... ";

  // Now for the good stuff: Erecon_ee for different channel windows
  TCanvas *c4 = new TCanvas("c4", "c4");
  c4->Divide(3,1);

  TH1D *h3100chan = new TH1D("100channels", "100 channel window", 500, 0, 5000);
  h3100chan->GetXaxis()->SetTitle("Erecon_ee (KeV)");
  TH1D *h3200chan = new TH1D("200channels", "200 channel window", 500, 0, 5000);
  h3200chan->GetXaxis()->SetTitle("Erecon_ee (KeV)");
  TH1D *h3300chan = new TH1D("300channels", "300 channel window", 500, 0, 5000);
  h3300chan->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  c4->cd(1);
  chain->Draw("Erecon_ee >> 100channels", basicCut && oneSTPeak_100channels);

  c4->cd(2);
  chain->Draw("Erecon_ee >> 200channels", basicCut && oneSTPeak_200channels);

  c4->cd(3);
  chain->Draw("Erecon_ee >> 300channels", basicCut && oneSTPeak_300channels);

  cout << "Done." << endl;

  cout << "Creating a fifth canvas to look at sub-range of allowed KE for Erecon_ee...";

  // The final coarse results: number of events in allowed subwindow
  TCanvas *c5 = new TCanvas("c5", "c5");
  c5->Divide(3,1);

  TH1D *h3100chan_subrange = new TH1D("100channels_subrange", "100 channel window", 80, 0, 800);
  h3100chan_subrange->GetXaxis()->SetTitle("Erecon_ee (KeV)");
  TH1D *h3200chan_subrange = new TH1D("200channels_subrange", "200 channel window", 80, 0, 800);
  h3200chan_subrange->GetXaxis()->SetTitle("Erecon_ee (KeV)");
  TH1D *h3300chan_subrange = new TH1D("300channels_subrange", "300 channel window", 80, 0, 800);
  h3300chan_subrange->GetXaxis()->SetTitle("Erecon_ee (KeV)");

  c5->cd(1);
  chain->Draw("Erecon_ee >> 100channels_subrange", basicCut && Erecon_ee_range && oneSTPeak_100channels);
  c5->cd(2);
  chain->Draw("Erecon_ee >> 200channels_subrange", basicCut && Erecon_ee_range && oneSTPeak_200channels);
  c5->cd(3);
  chain->Draw("Erecon_ee >> 300channels_subrange", basicCut && Erecon_ee_range && oneSTPeak_300channels);

  cout << "Done." << endl;

  // Making 2D histogram to look at Erecon_ee vs "true time" contour plot.

}
