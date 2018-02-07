#define	path	"/mnt/Data/xuansun/newReplayData_ee"
#define	filefgName	"Foreground_run_numbers_octets_80-120.txt"

{
  long double totalLiveTime = 0;
  TChain *mainfgChain = new TChain("pass3");


  Double_t time = 0;
  Int_t pid = -1;
  Int_t badTimeFlag = -1;

  string buf;
  ifstream infile;
  cout << "The file being opened is: " << filefgName << endl;
  infile.open(filefgName);

  if(!infile.is_open())
    cout << "Problem opening " << filefgName << endl;

  int fgRunIndex;

  int lastGoodEventCounter = 0;

  while(getline(infile, buf))
  {
    istringstream bufstream(buf);
    if(!bufstream.eof())
    {
      bufstream >> fgRunIndex;

/*      TChain *liveTimeChain = new TChain("pass3");
      liveTimeChain->SetBranchAddress("Time", &time);
      liveTimeChain->SetBranchAddress("PID", &pid);
      liveTimeChain->SetBranchAddress("badTimeFlag", &badTimeFlag);


      liveTimeChain->Add(Form("%s/replay_pass3_%i.root", path, bgRunIndex));


      for(unsigned int i = 0; i < liveTimeChain->GetEntries(); i++)
      {

        liveTimeChain->GetEntry(i);

        if(pid == 1 && badTimeFlag == 0)
        {
          lastGoodEventCounter = i;
        }
      }

      liveTimeChain->GetEntry(lastGoodEventCounter);

      totalLiveTime = totalLiveTime + time;

      delete liveTimeChain;
*/
      cout << "Adding run index " << fgRunIndex << " to our main chain." << endl;
      mainfgChain->Add(Form("%s/replay_pass3_%i.root", path, fgRunIndex));

      lastGoodEventCounter = 0;
    }

  }

//  cout << "THE TOTAL LIVE TIME IS: " << totalLiveTime << endl;

  infile.close();


}
