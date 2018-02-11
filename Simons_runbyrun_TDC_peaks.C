#define	path	"/mnt/Data/xuansun/newReplayData_ee"
#define	fileName	"Background_run_numbers_octets_80-120.txt"

void runbyrun_TDC_peaks()
{
  // define all the cuts we will use later.
  TCut basicCut = "(PID == 1 && badTimeFlag == 0 && Type == 1)";
  TCut Erecon_ee_range = "(Erecon_ee > 0 && Erecon_ee < 640)";
  // assuming TDCE self-timing peak at 2850 channels, TDCW at 3050
  TCut oneSTPeak_100channels = "((TDCE > 2850 && TDCW > 2950 && TDCW < 3050) || (TDCW > 3050 && TDCE > 2750 && TDCE < 2850))";
  TCut oneSTPeak_200channels = "((TDCE > 2850 && TDCW > 2850 && TDCW < 3050) || (TDCW > 3050 && TDCE > 2650 && TDCE < 2850))";
  TCut oneSTPeak_300channels = "((TDCE > 2850 && TDCW > 2750 && TDCW < 3050) || (TDCW > 3050 && TDCE > 2550 && TDCE < 2850))";
  TCut oneSTPeak_400channels = "((TDCE > 2850 && TDCW > 2650 && TDCW < 3050) || (TDCW > 3050 && TDCE > 2450 && TDCW < 2850))";
  TCut oneSTPeak_500channels = "((TDCE > 2850 && TDCW > 2550 && TDCW < 3050) || (TDCW > 3050 && TDCE > 2350 && TDCW < 2850))";
  TCut radialCutE = "((xE.center)*(xE.center) + (yE.center)*(yE.center) < 49*49)";
  TCut radialCutW = "((xW.center)*(xW.center) + (yW.center)*(yW.center) < 49*49)";

  double EmaxBin = -1;
  double EsigmaLoBin = -1;
  double EsigmaHiBin = -1;
  double EsigmaAvg = -1;
  double WmaxBin = -1;
  double WsigmaLoBin = -1;
  double WsigmaHiBin = -1;
  double WsigmaAvg = -1;

  ofstream outfile;
  outfile.open("TDC_STPeak_runbyrun.txt");


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


      chain->Draw("TDCE >> hTDCE", basicCut && radialCutE && radialCutW);
      chain->Draw("TDCW >> hTDCW", basicCut && radialCutE && radialCutW);

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

      cout << "Finshed saving run index " << runIndex << " to file." << endl;

      delete chain;
      delete hTDCE;
      delete hTDCW;
    }
  }

  infile.close();
  outfile.close();
}
