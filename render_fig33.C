{
  gROOT->SetStyle("Pub");

  TString fileName = "FinalCounts_variousTimeWindows_allCuts_bgfg.txt";

  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->cd();

  vector <double> yValW0;
  vector <double> yErrW0;
  vector <double> xValW0;
  vector <double> xErrW0;

  double windowLower, windowUpper;
  double bgValid, fgValid;
  double bgErr, fgErr;

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
      bufstream >> windowLower >> windowUpper
		>> bgValid >> fgValid;


	bgErr = sqrt(bgValid);
	fgErr = sqrt(fgValid);

	yValW0.push_back(fgValid / (bgValid*5.07));
	yErrW0.push_back((fgValid/(bgValid*5.07)) * sqrt((bgErr/bgValid)*(bgErr/bgValid) + (fgErr/fgValid)*(fgErr/fgValid)));
	xValW0.push_back(windowUpper);
	xErrW0.push_back(0.05);
    }
  }
  infile.close();

  TGraphErrors *g0 = new TGraphErrors(xValW0.size(), &(xValW0[0]), &(yValW0[0]), &(xErrW0[0]), &(yErrW0[0]));
  g0->SetMarkerStyle(21);
  g0->SetMarkerColor(1);
  g0->GetHistogram()->GetYaxis()->SetRangeUser(0, 3);
  g0->Draw("AP");


  TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry(g0,"West Low Edge 0ns","lep");
//  legend->Draw();

  g0->SetTitle("Various Time Windows");
  g0->GetXaxis()->SetTitle("Upper Edge of Window (ns)");
  g0->GetXaxis()->CenterTitle();
  g0->GetYaxis()->SetTitle("N_{FG}/5.07N_{BG}");
  g0->GetYaxis()->CenterTitle();

  c1->Print("33_Ratio_bgfg_variousTimeWindows_allCuts.pdf");
}
