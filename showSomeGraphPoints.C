{
  TString fileName = "FinalCounts_variousTimeWindows_allForegroundRuns.txt";

  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->cd();
/*
  double x[5] = {10, 11, 12, 13, 14};
  double y[5] = {9723.0/1863.0, 15935.0/3081.0, 23853.0/4540.0, 34313.0/6371.0, 47670.0/8600.0};
  TGraph *g = new TGraph(5, x, y);
  g->SetMarkerStyle(8);
  g->Draw("AP");
*/

  vector <double> yValWn2;
  vector <double> yErrWn2;
  vector <double> yValW0;
  vector <double> yErrW0;
  vector <double> yValW2;
  vector <double> yErrW2;
  vector <double> xValWn2;
  vector <double> xErrWn2;
  vector <double> xValW0;
  vector <double> xErrW0;
  vector <double> xValW2;
  vector <double> xErrW2;

  double eLow, eHigh, wLow, wHigh;
  double bgAll, fgAll, bgValid, fgValid;

  string buf;
  ifstream infile;
  cout << "The file being opened is: " << fileName << endl;
  infile.open(fileName);

  if(!infile.is_open())
    cout << "Problem opening " << fileName << endl;

  while(getline(infile, buf))
  {
    istringstream bufstream(buf);
    if(!bufstream.eof())
    {
      bufstream >> eLow >> eHigh >> wLow >> wHigh
		>> bgAll >> fgAll >> bgValid >> fgValid;

      if(wLow == -2)
      {
	yValWn2.push_back(fgValid);
	yErrWn2.push_back(sqrt(fgValid));
	xValWn2.push_back(wHigh);
	xErrWn2.push_back(1);
      }
      else if(wLow == 0)
      {
	yValW0.push_back(fgValid);
	yErrW0.push_back(sqrt(fgValid));
	xValW0.push_back(wHigh);
	xErrW0.push_back(1);
      }
      else if(wLow == 2)
      {
	yValW2.push_back(fgValid);
	yErrW2.push_back(sqrt(fgValid));
	xValW2.push_back(wHigh);
	xErrW2.push_back(1);
      }

    }
  }
  infile.close();

  TGraphErrors *gn2 = new TGraphErrors(xValWn2.size(), &(xValWn2[0]), &(yValWn2[0]), &(xErrWn2[0]), &(yErrWn2[0]));
  gn2->SetMarkerStyle(22);
  gn2->SetMarkerColor(1);
  gn2->Draw("AP");

  TGraphErrors *g0 = new TGraphErrors(xValW0.size(), &(xValW0[0]), &(yValW0[0]), &(xErrW0[0]), &(yErrW0[0]));
  g0->SetMarkerStyle(22);
  g0->SetMarkerColor(2);
  g0->Draw("PSAME");

  TGraphErrors *g2 = new TGraphErrors(xValW2.size(), &(xValW2[0]), &(yValW2[0]), &(xErrW2[0]), &(yErrW2[0]));
  g2->SetMarkerStyle(22);
  g2->SetMarkerColor(3);
  g2->Draw("PSAME");


  TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry(gn2,"West Low Edge -2ns","lep");
  legend->AddEntry(g0,"West Low Edge 0ns","lep");
  legend->AddEntry(g2,"West Low Edge 2ns","lep");
  legend->Draw();

  gn2->SetTitle("Counts as a function of West High Edge of Accepted Timing Window");
  gn2->GetXaxis()->SetTitle("West Upper Edge of Window (ns)");
  gn2->GetYaxis()->SetTitle("Counts (N)");

  c1->Print("Counts_variousTimeWindows_lowAndHighEdges_allCuts.pdf");
}
