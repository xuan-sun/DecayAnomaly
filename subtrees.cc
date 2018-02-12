#include	 <iostream>
#include	 <fstream>
#include	 <TGaxis.h>
#include	 <sstream>
#include	 <TGraph.h>
#include	 <TGraphErrors.h>
#include	 <TCanvas.h>
#include	 <TApplication.h>
#include	 <stdlib.h>
#include	 <TF1.h>
#include	 <TH1.h>
#include	 <TProfile.h>
#include	 <TObjArray.h>
#include	 <TStyle.h>
#include	 <TMarker.h>
#include	 <math.h>
#include	 <TStyle.h>
#include	 <TPaveStats.h>
#include	 <TPaveText.h>
#include	 <vector>
#include	 <string.h>
#include	 <fstream>
#include	 <TROOT.h>
#include	 <TFile.h>
#include	 <TLegend.h>
#include         <TLegendEntry.h>
#include	 <time.h>
#include	 <TH2F.h>
#include         <assert.h>
#include	 <string>
#include	 <TRandom.h>
#include 	 <TTree.h>
#include	 <TChain.h>
#include	 <TVector.h>
#include	 <vector>
#include	 <utility>
#include	 <TLeaf.h>
#include	 <math.h>

using		 namespace std;

struct Event
{
  int TriggerNum;
  int EvtN;
  double DeltaT;
  double Tof;
  double TimeE;
  double TimeW;
  double Time;
  double TDCE;
  double TDCW;
  double TDCE1;
  double TDCE2;
  double TDCE3;
  double TDCE4;
  double TDCW1;
  double TDCW2;
  double TDCW3;
  double TDCW4;
  double EvisE;
  double EvisW;
  double CathSumE;
  double CathSumW;
  double CathMaxE;
  double CathMaxW;
  double EMWPC_E;
  double EMWPC_W;
  double AnodeE;
  double AnodeW;
  int PID;
  int Side;
  int Type;
  double Erecon;
  double Erecon_ee;
  int badTimeFlag;
  double xE_center;
  double yE_center;
  double xW_center;
  double yW_center;
};

// forward declarations for useful functions
vector < pair <string,int> >  LoadOctetList(TString fileName);
vector < TChain* > GetChainsOfRuns(vector < pair <string,int> > octetList, TString dataPath);
vector < TTree* > CreateOctetTrees(vector <TChain*> runsChains);

// these are actual beta run indices
const int index_A2 = 0;
const int index_A5 = 1;
const int index_A7 = 2;
const int index_A10 = 3;
const int index_B2 = 4;
const int index_B5 = 5;
const int index_B7 = 6;
const int index_B10 = 7;

// these are the background runs
// they correspond to beta run index + 8
const int index_A1 = 8;
const int index_A4 = 9;
const int index_A9 = 10;
const int index_A12 = 11;
const int index_B1 = 12;
const int index_B4 = 13;
const int index_B9 = 14;
const int index_B12 = 15;

// Used for visualization, keeps the graph on screen.
//TApplication plot_program("FADC_readin",0,0,0,0);

//-------------------------------------------------//
//------------ Start of Program -------------------//
//-------------------------------------------------//

int main(int argc, char* argv[])
{
  if(argc < 2)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable) (octet #)" << endl;
    return 0;
  }

  // creating canvas for plotting
  TCanvas *C = new TCanvas("canvas", "canvas", 800, 400);

  // read in arguments.
  Int_t octNb = atoi(argv[1]);

  // Reads in the octet list and saves the run files indices corresponding to an octet number
  vector < pair <string,int> > octetIndices = LoadOctetList(TString::Format("%s/octet_list_%i.dat", "OctetLists", octNb));
  // Points TChains at the run files idenified in the octet lists above
  vector < TChain* > runFiles = GetChainsOfRuns(octetIndices, "/mnt/Data/xuansun/newReplayData_ee/");

  cout << "Chained all the runs together... " << endl;

  // save subtrees with the data listed in Event
  vector < TTree* > contents = CreateOctetTrees(runFiles);

  cout << "vector of TTree* called contents made..." << endl;

  cout << "Now writing to file... " << endl;

  TList *allTreesList = new TList();
  for(unsigned int i = 0; i < contents.size(); i++)
  {
    if(contents[i] == NULL)
    {
      cout << "contents[" << i << "] is NULL" << endl;
    }

    TFile f(TString::Format("Octet_%i_type1_runIndex_%i.root", octNb, i), "RECREATE");
    contents[i]->Write();
    f.Close();



  }

  cout << "Finished saving TFiles for every run... " << endl;

  cout << "Loading and merging all the TFiles into one TChain... " << endl;

  TChain* allRunsChainsInOctet = new TChain("pass3");
  for(unsigned int i = 0; i < contents.size(); i++)
  {
    allRunsChainsInOctet->Add(TString::Format("Octet_%i_type1_runIndex_%i.root", octNb, i));
  }
  allRunsChainsInOctet->Merge(TString::Format("Octet_%i_type1.root", octNb));




//  f.Close();
  cout << "-------------- End of Program ---------------" << endl;
//  plot_program.Run();

  return 0;
}

vector < pair <string,int> >  LoadOctetList(TString fileName)
{
  string buf;
  ifstream infile;
  cout << "The file being opened is: " << fileName << endl;
  infile.open(fileName);

  if(!infile.is_open())
    cout << "Problem opening " << fileName << endl;

  string runType;
  int runIndex;
  vector <pair <string, int> > pairs;

  while(getline(infile, buf))
  {
    istringstream bufstream(buf);
    if(!bufstream.eof())
    {
      bufstream >> runType >> runIndex;
      if(runType == "A2")
      {
        pairs.push_back(make_pair("A2", runIndex));
      }
      else if(runType == "A5")
      {
        pairs.push_back(make_pair("A5", runIndex));
      }
      else if(runType == "A7")
      {
        pairs.push_back(make_pair("A7", runIndex));
      }
      else if(runType == "A10")
      {
        pairs.push_back(make_pair("A10", runIndex));
      }
      else if(runType == "B2")
      {
        pairs.push_back(make_pair("B2", runIndex));
      }
      else if(runType == "B5")
      {
        pairs.push_back(make_pair("B5", runIndex));
      }
      else if(runType == "B7")
      {
        pairs.push_back(make_pair("B7", runIndex));
      }
      else if(runType == "B10")
      {
        pairs.push_back(make_pair("B10", runIndex));
      }
      // Save the background runs as well
      else if(runType == "A1")
      {
	pairs.push_back(make_pair("A1", runIndex));
      }
      else if(runType == "A4")
      {
        pairs.push_back(make_pair("A4", runIndex));
      }
      else if(runType == "A9")
      {
        pairs.push_back(make_pair("A9", runIndex));
      }
      else if(runType == "A12")
      {
        pairs.push_back(make_pair("A12", runIndex));
      }
      else if(runType == "B1")
      {
        pairs.push_back(make_pair("B1", runIndex));
      }
      else if(runType == "B4")
      {
        pairs.push_back(make_pair("B4", runIndex));
      }
      else if(runType == "B9")
      {
        pairs.push_back(make_pair("B9", runIndex));
      }
      else if(runType == "B12")
      {
        pairs.push_back(make_pair("B12", runIndex));
      }
    }
  }

  infile.close();       // ensure you're closing the read-in file.

  return pairs;
}

vector < TChain* > GetChainsOfRuns(vector < pair <string,int> > octetList, TString dataPath)
{
  vector < TChain* > runs;

  for(unsigned int i = 0; i < 16; i++)
  {
    runs.push_back(new TChain("pass3"));
  }

  for(unsigned int l = 0; l < octetList.size(); l++)
  {
    if(octetList[l].first == "A2")
    {
      runs[index_A2] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "A5")
    {
      runs[index_A5] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "A7")
    {
      runs[index_A7] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "A10")
    {
      runs[index_A10] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B2")
    {
      runs[index_B2] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B5")
    {
      runs[index_B5] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B7")
    {
      runs[index_B7] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B10")
    {
      runs[index_B10] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "A1")
    {
      runs[index_A1] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "A4")
    {
      runs[index_A4] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "A9")
    {
      runs[index_A9] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "A12")
    {
      runs[index_A12] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B1")
    {
      runs[index_B1] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B4")
    {
      runs[index_B4] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B9")
    {
      runs[index_B9] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B12")
    {
      runs[index_B12] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
  }

  return runs;
}

vector <TTree*> CreateOctetTrees(vector <TChain*> runsChains)
{
  // reminder: each evt[i] index corresponds to the indices noted at global scope.
  vector <Event*> evt;
  vector <TTree*> subtrees;

  for(unsigned int i = 0; i < runsChains.size(); i++)
  {
    evt.push_back(new Event);
    subtrees.push_back(runsChains[i]->CloneTree(0));

    runsChains[i]->SetBranchAddress("TriggerNum", &evt[i]->TriggerNum);
    runsChains[i]->SetBranchAddress("EvtN", &evt[i]->EvtN);
    runsChains[i]->SetBranchAddress("DeltaT", &evt[i]->DeltaT);
    runsChains[i]->SetBranchAddress("Tof", &evt[i]->Tof);
    runsChains[i]->SetBranchAddress("Time", &evt[i]->Time);
    runsChains[i]->SetBranchAddress("TimeE", &evt[i]->TimeE);
    runsChains[i]->SetBranchAddress("TimeW", &evt[i]->TimeW);
    runsChains[i]->SetBranchAddress("TDCE", &evt[i]->TDCE);
    runsChains[i]->SetBranchAddress("TDCW", &evt[i]->TDCW);
    runsChains[i]->SetBranchAddress("TDCE1", &evt[i]->TDCE1);
    runsChains[i]->SetBranchAddress("TDCE2", &evt[i]->TDCE2);
    runsChains[i]->SetBranchAddress("TDCE3", &evt[i]->TDCE3);
    runsChains[i]->SetBranchAddress("TDCE4", &evt[i]->TDCE4);
    runsChains[i]->SetBranchAddress("TDCW1", &evt[i]->TDCW1);
    runsChains[i]->SetBranchAddress("TDCW2", &evt[i]->TDCW2);
    runsChains[i]->SetBranchAddress("TDCW3", &evt[i]->TDCW3);
    runsChains[i]->SetBranchAddress("TDCW4", &evt[i]->TDCW4);
    runsChains[i]->SetBranchAddress("EvisE", &evt[i]->EvisE);
    runsChains[i]->SetBranchAddress("EvisW", &evt[i]->EvisW);
    runsChains[i]->SetBranchAddress("CathSumE", &evt[i]->CathSumE);
    runsChains[i]->SetBranchAddress("CathSumW", &evt[i]->CathSumW);
    runsChains[i]->SetBranchAddress("EMWPC_E", &evt[i]->EMWPC_E);
    runsChains[i]->SetBranchAddress("EMWPC_W", &evt[i]->EMWPC_W);
    runsChains[i]->SetBranchAddress("AnodeE", &evt[i]->AnodeE);
    runsChains[i]->SetBranchAddress("AnodeW", &evt[i]->AnodeW);
    runsChains[i]->SetBranchAddress("PID", &evt[i]->PID);
    runsChains[i]->SetBranchAddress("Side", &evt[i]->Side);
    runsChains[i]->SetBranchAddress("Type", &evt[i]->Type);
    runsChains[i]->SetBranchAddress("Erecon", &evt[i]->Erecon);
    runsChains[i]->SetBranchAddress("Erecon_ee", &evt[i]->Erecon_ee);
    runsChains[i]->SetBranchAddress("badTimeFlag", &evt[i]->badTimeFlag);
    // this additional syntax is needed to get the right leaf inside branch inside tree named "pass3"
    runsChains[i]->GetBranch("xE")->GetLeaf("center")->SetAddress(&evt[i]->xE_center);
    runsChains[i]->GetBranch("yE")->GetLeaf("center")->SetAddress(&evt[i]->yE_center);
    runsChains[i]->GetBranch("xW")->GetLeaf("center")->SetAddress(&evt[i]->xW_center);
    runsChains[i]->GetBranch("yW")->GetLeaf("center")->SetAddress(&evt[i]->yW_center);
  }


  for(unsigned int j = 0; j < runsChains.size(); j++)
  {
    for(unsigned int i = 0; i < runsChains[j]->GetEntries(); i++)
    {
      runsChains[j]->GetEntry(i);	/* this is where cuts happen */
      if(evt[j]->Type == 1)
      {
        subtrees[j]->Fill();
      }
    }
  }

  return subtrees;

}

