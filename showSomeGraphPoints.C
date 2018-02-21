{
  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->cd();

  double x[5] = {10, 11, 12, 13, 14};
  double y[5] = {9723.0/1863.0, 15935.0/3081.0, 23853.0/4540.0, 34313.0/6371.0, 47670.0/8600.0};

  TGraph *g = new TGraph(5, x, y);
  g->SetMarkerStyle(8);
  g->Draw("AP");


  c1->Print("Ratio_FG2BG_variousTimeWindows.pdf");
}
