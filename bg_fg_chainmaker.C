#define	path	"Type1_runs_BG_FG"

{
  TChain *mainFGChain = new TChain("pass3");
  for(int i = 0; i < 299; i++)
  {
    if(i == 107)	// this is for run number 22199 which is FG and peak finder is wrong
    {
      continue;
    }
    mainFGChain->Add(Form("%s/test_%i.root", path, i));
  }
  mainFGChain->Merge("TimeCalibrated_FGRuns_type1.root");

  TChain *mainBGChain = new TChain("pass3");
  for(int i = 299; i < 595; i++)
  {
    if(i == 554 || i == 556)	// run numbers 23017 and 23021
    {
      continue;
    }

    mainBGChain->Add(Form("%s/test_%i.root", path, i));
  }
  mainBGChain->Merge("TimeCalibrated_BGRuns_type1.root");

}
