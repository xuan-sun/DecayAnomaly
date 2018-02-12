#define	path	"/mnt/Data/xuansun/newReplayData_ee"
#define	fileName	"Foreground_run_numbers_octets_80-120.txt"

void Simons_runbyrun_TDC_peaks()
{

  // define all the cuts we will use later.
  TCut basicCut = "(PID == 1 && badTimeFlag == 0)";
  TCut fiducialCut = "(((xE.center)*(xE.center) + (yE.center)*(yE.center) < 49*49) && ((xW.center)*(xW.center) + (yW.center)*(yW.center) < 49*49))";
  TCut typeCut = "Type == 1";


  double EmaxBin = -1;
  double EsigmaLoBin = -1;
  double EsigmaHiBin = -1;
  double EsigmaAvg = -1;
  double WmaxBin = -1;
  double WsigmaLoBin = -1;
  double WsigmaHiBin = -1;
  double WsigmaAvg = -1;

  ofstream outfile;
  outfile.open("Foreground_TDC_STPeak_runbyrun_v2.txt");


  // first read in and save the right events to file
  string buf;
  ifstream infile;
  cout << "The file being opened is: " << fileName << endl;
  infile.open(fileName);

  if(!infile.is_open())
    cout << "Problem opening " << fileName << endl;

  int runIndex;


  while(getline(infile, buf))
  {
    istringstream bufstream(buf);
    if(!bufstream.eof())
    {
      bufstream >> runIndex;


      TChain *chain = new TChain("pass3");
      chain->Add(Form("%s/replay_pass3_%i.root", path, runIndex));

      TH1D *hTDCE = new TH1D("hTDCE", "hTDCE", 5000, 0, 5000);
      TH1D *hTDCW = new TH1D("hTDCW", "hTDCW", 5000, 0, 5000);


      chain->Draw("TDCE >> hTDCE", basicCut && fiducialCut && typeCut);
      chain->Draw("TDCW >> hTDCW", basicCut && fiducialCut && typeCut);

      EmaxBin = hTDCE->GetBinLowEdge(hTDCE->GetMaximumBin());
      EsigmaLoBin = hTDCE->GetBinLowEdge(hTDCE->FindFirstBinAbove(0.6*hTDCE->GetBinContent(hTDCE->GetMaximumBin())));
      EsigmaHiBin = hTDCE->GetBinLowEdge(hTDCE->FindLastBinAbove(0.6*hTDCE->GetBinContent(hTDCE->GetMaximumBin())));
      EsigmaAvg = ( (EmaxBin-EsigmaLoBin) + (EsigmaHiBin-EmaxBin) )/2;

      WmaxBin = hTDCW->GetBinLowEdge(hTDCW->GetMaximumBin());
      WsigmaLoBin = hTDCW->GetBinLowEdge(hTDCW->FindFirstBinAbove(0.6*hTDCW->GetBinContent(hTDCW->GetMaximumBin())));
      WsigmaHiBin = hTDCW->GetBinLowEdge(hTDCW->FindLastBinAbove(0.6*hTDCW->GetBinContent(hTDCW->GetMaximumBin())));
      WsigmaAvg = ( (WmaxBin-WsigmaLoBin) + (WsigmaHiBin-WmaxBin) )/2;

      outfile << runIndex << "\t"
	      << EmaxBin << "\t"
	      << EsigmaLoBin << "\t"
	      << EsigmaHiBin << "\t"
	      << EsigmaAvg << "\t"
	      << WmaxBin << "\t"
	      << WsigmaLoBin << "\t"
	      << WsigmaHiBin << "\t"
	      << WsigmaAvg << "\n";

      cout << "Finished saving run index " << runIndex << " to file." << endl;

      delete chain;
      delete hTDCE;
      delete hTDCW;
    }
  }

  infile.close();
  outfile.close();
}
