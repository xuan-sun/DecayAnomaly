#define	path	"/mnt/Data/xuansun/newReplayData_ee"
#define	filebgName	"Foreground_run_numbers_octets_80-120.txt"

{
  Double_t time = 0;
  Int_t pid = -1;
  Int_t badTimeFlag = -1;
  Int_t Type = -1;


  // first read in and save the right events to file
  string buf;
  ifstream infile;
  cout << "The file being opened is: " << filebgName << endl;
  infile.open(filebgName);

  if(!infile.is_open())
    cout << "Problem opening " << filebgName << endl;

  int runIndex;


  while(getline(infile, buf))
  {
    istringstream bufstream(buf);
    if(!bufstream.eof())
    {
      bufstream >> runIndex;


      TFile f(Form("replay_pass3_%i_type1.root", runIndex), "RECREATE");

      TChain *chain = new TChain("pass3");
      chain->SetBranchAddress("Type", &Type);


      chain->Add(Form("%s/replay_pass3_%i.root", path, runIndex));


      TTree *subtree = chain->CloneTree(0);

      for(unsigned int i = 0; i < chain->GetEntries(); i++)
      {
	chain->GetEntry(i);
	if(Type == 1)
	{
	  subtree->Fill();
	}
      }

      subtree->Write();


      delete subtree;
      delete chain;

      f.Close();

      cout << "Finshed saving run index " << runIndex << " to file." << endl;
    }
  }

  infile.close();

  // now open the run index files again and merge the chains
  // first read in and save the right events to file
  TChain *mainChain = new TChain("pass3");

  string buf2;
  ifstream infile2;
  cout << "The file being opened is: " << filebgName << endl;
  infile2.open(filebgName);

  if(!infile2.is_open())
    cout << "Problem opening " << filebgName << endl;

  int runIndex2;


  while(getline(infile2, buf2))
  {
    istringstream bufstream(buf2);
    if(!bufstream.eof())
    {
      bufstream >> runIndex2;

      mainChain->Add(Form("replay_pass3_%i_type1.root", runIndex2));

    }
  }

  infile.close();

  mainChain->Merge("AllFGRuns_Octets_80-120_type1.root");
}
