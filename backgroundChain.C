#define	path	"/mnt/Data/xuansun/newReplayData_ee"
#define	filebgName	"Background_run_numbers_octets_80-120.txt"

{
  long double totalLiveTime = 0;
  TChain *mainbgChain = new TChain("pass3");


  Double_t time = 0;
  Int_t pid = -1;
  Int_t badTimeFlag = -1;

  string buf;
  ifstream infile;
  cout << "The file being opened is: " << filebgName << endl;
  infile.open(filebgName);

  if(!infile.is_open())
    cout << "Problem opening " << filebgName << endl;

  int bgRunIndex;

  int lastGoodEventCounter = 0;

  while(getline(infile, buf))
  {
    istringstream bufstream(buf);
    if(!bufstream.eof())
    {
      bufstream >> bgRunIndex;

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
      cout << "Adding run index " << bgRunIndex << " to our main chain." << endl;
      mainbgChain->Add(Form("%s/replay_pass3_%i.root", path, bgRunIndex));

      lastGoodEventCounter = 0;
    }

  }

  infile.close();


}
