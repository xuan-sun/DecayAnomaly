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

#define		 PATH		"/mnt/Data/xuansun/newReplayData_ee"

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
//  double xE_center;
//  double yE_center;
//  double xW_center;
//  double yW_center;
};

vector <int> runIndices;
vector <double> eSTP;	// east self-timing peak channel
vector <double> wSTP;	// west self-timing peak channel


// forward declarations for useful functions
void FillGlobalVectorsFromFile(TString fileName);
TTree* AddBranchToClonedTree(int runNumber, int lineIndex);

// Used for visualization, keeps the graph on screen.
//TApplication plot_program("FADC_readin",0,0,0,0);

//-------------------------------------------------//
//------------ Start of Program -------------------//
//-------------------------------------------------//

int main(int argc, char* argv[])
{
  // creating canvas for plotting
  TCanvas *C = new TCanvas("canvas", "canvas", 800, 400);

  FillGlobalVectorsFromFile("Foreground_TDC_STPeak_runbyrun.txt");
  FillGlobalVectorsFromFile("Background_TDC_STPeak_runbyrun.txt");

  for(unsigned int i = 0; i < runIndices.size(); i++)
  {
    TFile f(Form("newTimeCalib_replay_pass3_%i_type1.root", runIndices[i]), "RECREATE");
    cout << "runIndices[" << i << "] = " << runIndices[i] << endl;
    TTree *t = AddBranchToClonedTree(runIndices[i], i);
    cout << "Writing to file now.." << endl;

    t->Write();
    cout << "Write completed" << endl;

    delete t;

    f.Close();
    cout << "File closed" << endl;


  }


  cout << "-------------- End of Program ---------------" << endl;
//  plot_program.Run();

  return 0;
}

void FillGlobalVectorsFromFile(TString fileName)
{
  // note we call it Ebin and Wbin but the bins are 1 channel width so it's also units of channel
  int runIndex = -1;
  double EbinMax, ElowerEdge, EUpperEdge, Eavg;
  double WbinMax, WlowerEdge, WUpperEdge, Wavg;


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
      bufstream >> runIndex
	        >> EbinMax >> ElowerEdge >> EUpperEdge >> Eavg
	        >> WbinMax >> WlowerEdge >> WUpperEdge >> Wavg;

      runIndices.push_back(runIndex);
      eSTP.push_back(EbinMax);
      wSTP.push_back(WbinMax);

    }
  }

  infile.close();       // ensure you're closing the read-in file.

  cout << "Done filling in data from " << fileName <<endl;
}


//vector <TTree*> CreateOctetTrees(vector <TChain*> runsChains)
TTree* AddBranchToClonedTree(int runNumber, int lineIndex)
{
  TChain *chain = new TChain("pass3");
  chain->Add(Form("%s/replay_pass3_%i.root", PATH, runNumber));

  double newTimeE = -1;
  double newTimeW = -1;

  Event* evt = new Event;
  TTree* calibratedSubTree = chain->CloneTree(0);

  chain->SetBranchAddress("TriggerNum", &evt->TriggerNum);
  chain->SetBranchAddress("EvtN", &evt->EvtN);
  chain->SetBranchAddress("DeltaT", &evt->DeltaT);
  chain->SetBranchAddress("Tof", &evt->Tof);
  chain->SetBranchAddress("Time", &evt->Time);
  chain->SetBranchAddress("TimeE", &evt->TimeE);
  chain->SetBranchAddress("TimeW", &evt->TimeW);
  chain->SetBranchAddress("TDCE", &evt->TDCE);
  chain->SetBranchAddress("TDCW", &evt->TDCW);
  chain->SetBranchAddress("TDCE1", &evt->TDCE1);
  chain->SetBranchAddress("TDCE2", &evt->TDCE2);
  chain->SetBranchAddress("TDCE3", &evt->TDCE3);
  chain->SetBranchAddress("TDCE4", &evt->TDCE4);
  chain->SetBranchAddress("TDCW1", &evt->TDCW1);
  chain->SetBranchAddress("TDCW2", &evt->TDCW2);
  chain->SetBranchAddress("TDCW3", &evt->TDCW3);
  chain->SetBranchAddress("TDCW4", &evt->TDCW4);
  chain->SetBranchAddress("EvisE", &evt->EvisE);
  chain->SetBranchAddress("EvisW", &evt->EvisW);
  chain->SetBranchAddress("CathSumE", &evt->CathSumE);
  chain->SetBranchAddress("CathSumW", &evt->CathSumW);
  chain->SetBranchAddress("EMWPC_E", &evt->EMWPC_E);
  chain->SetBranchAddress("EMWPC_W", &evt->EMWPC_W);
  chain->SetBranchAddress("AnodeE", &evt->AnodeE);
  chain->SetBranchAddress("AnodeW", &evt->AnodeW);
  chain->SetBranchAddress("PID", &evt->PID);
  chain->SetBranchAddress("Side", &evt->Side);
  chain->SetBranchAddress("Type", &evt->Type);
  chain->SetBranchAddress("Erecon", &evt->Erecon);
  chain->SetBranchAddress("Erecon_ee", &evt->Erecon_ee);
  chain->SetBranchAddress("badTimeFlag", &evt->badTimeFlag);
/*  chain->GetBranch("xE")->GetLeaf("center")->SetAddress(&evt->xE_center);
  chain->GetBranch("yE")->GetLeaf("center")->SetAddress(&evt->yE_center);
  chain->GetBranch("xW")->GetLeaf("center")->SetAddress(&evt->xW_center);
  chain->GetBranch("yW")->GetLeaf("center")->SetAddress(&evt->yW_center);
*/
  TBranch *bTimeE = calibratedSubTree->Branch("newTimeE", &newTimeE, "newTimeE/D");
  TBranch *bTimeW = calibratedSubTree->Branch("newTimeW", &newTimeW, "newTimeW/D");


  for(unsigned int i = 0; i < chain->GetEntries(); i++)
  {
    chain->GetEntry(i);
    if(evt->Type == 1)
    {
      newTimeE = 140 - (evt->TDCE * 140.0 / eSTP[lineIndex]);
      newTimeW = 140 - (evt->TDCW * 140.0 / wSTP[lineIndex]);

      bTimeE->Fill();
      bTimeW->Fill();
      calibratedSubTree->Fill();
    }

  }

  delete evt;
  delete chain;

  return calibratedSubTree;
}

