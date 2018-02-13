#define	path	"Type1_runs_BG_FG"

TChain* MergeFileInChain(TString fileName, int flag)
{
  TChain *mainChain = new TChain("pass3");

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

      if(flag == 0)	// background runs
      {
        if(runIndex == 23017 || runIndex == 23021)
	{
	  continue;
	}
      }
      else if(flag == 1)	// foreground runs
      {
	if(runIndex == 22199)
	{
	  continue;
	}
      }

      mainChain->Add(Form("%s/newTimeCalib_replay_pass3_%i_type1.root", path, runIndex));

    }

  }
  infile.close();

  cout << "Done chaining the files together." << endl;

  return mainChain;
}

runbyrun_chainmaker()
{
  TChain *mainBGChain = MergeFileInChain("Background_run_numbers_octets_80-120.txt", 0);
  cout << "Merging all chains into one file..." << endl;
  mainBGChain->Merge("TimeCalibrated_BGRuns_type1_fixed.root", "RECREATE");
  cout << "Done." << endl;


  TChain *mainFGChain = MergeFileInChain("Foreground_run_numbers_octets_80-120.txt", 1);
  cout << "Merging all chains into one file..." << endl;
  mainFGChain->Merge("TimeCalibrated_FGRuns_type1_fixed.root", "RECREATE");
  cout << "Done." << endl;

}


