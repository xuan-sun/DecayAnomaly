#define	path	"Type1_Octets_80-120"

{
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

  // make plot of the TDCE and TDCW full range values
  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->Divide(2,1);
  c1->cd(1);

  // create histograms (our choice of binning) for TDCE, TDCW totals.
  TH1D *h1tdce = new TH1D("h1tdce", "h1 tdce", 800, 0, 4000);
  TH1D *h1tdcw = new TH1D("h1tdcw", "h1 tdcw", 800, 0, 4000);

  chain->Draw("TDCE >> h1tdce", "PID == 1 && badTimeFlag == 0");
  gPad->SetLogy();

  c1->cd(2);
  chain->Draw("TDCW >> h1tdcw", "PID == 1 && badTimeFlag == 0");
  gPad->SetLogy();

  // make second plot of TDCE/TDCW values, sub range to look at "gap"
  TCanvas *c2 = new TCanvas("c2", "c2");
  c2->Divide(2,1);
  c2->cd(1);

  TH1D *h1tdce_subrange = new TH1D("h1tdce_subrange", "h1 tdce subrange", 160, 2400, 3200);
  TH1D *h1tdcw_subrange = new TH1D("h1tdcw_subrange", "h1 tdcw subrange", 160, 2400, 3200);

  h1tdce_subrange->GetYaxis()->SetRangeUser(1, 5000);
  chain->Draw("TDCE >> h1tdce_subrange", "PID == 1 && badTimeFlag == 0");

  c2->cd(2);
  h1tdcw_subrange->GetYaxis()->SetRangeUser(1, 5000);
  chain->Draw("TDCW >> h1tdce_subrange", "PID == 1 && badTimeFlag == 0");


}
