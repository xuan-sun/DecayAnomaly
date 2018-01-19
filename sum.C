{
  TFile f;

  TH1D *hTotal = new TH1D("total", "Erecon Octets 80-120", 150, 0, 1500);


  for(int i = 80; i < 122; i++)
  {
    if(i == 91 || i == 93 || i == 101 || i == 107 || i == 121)
    {
      continue;
    }

    f.Open(Form("Octet_%i_TDC_selfTimingCuts_Erecon_type1.root", i));
    TH1D *hTemp = (TH1D*)f.Get("Total hist");
    hTotal->Add(hTemp, 1.0);
  }

  hTotal->Draw();



}
