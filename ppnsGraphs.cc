#include         <iostream>
#include         <fstream>
#include         <TGaxis.h>
#include         <sstream>
#include         <TGraph.h>
#include         <TGraphErrors.h>
#include         <TCanvas.h>
#include         <TApplication.h>
#include         <stdlib.h>
#include         <TF1.h>
#include         <TH1.h>
#include         <TProfile.h>
#include         <TObjArray.h>
#include         <TStyle.h>
#include         <TMarker.h>
#include         <math.h>
#include         <TStyle.h>
#include         <TPaveStats.h>
#include         <TPaveText.h>
#include         <vector>
#include         <string.h>
#include         <fstream>
#include         <TROOT.h>
#include         <TFile.h>
#include         <TLegend.h>
#include         <TLegendEntry.h>
#include         <time.h>
#include         <TH2F.h>
#include         <assert.h>
#include         <string>
#include         <TRandom.h>
#include         <TTree.h>
#include         <TChain.h>
#include         <TVector.h>
#include         <vector>
#include         <utility>
#include         <TLeaf.h>
#include         <math.h>
#include	 <TCut.h>
#include	 <TMath.h>
#include	 <TRandom3.h>

using            namespace std;

// some forward declarations of useful functions
double SetPoissonErrors(int counts);
TH1D* FillInTimingCounts(TString fileName);

struct SimEvent
{
  double Edep_EdepE;
  double Edep_EdepW;
  double EdepQ_EdepQE;
  double EdepQ_EdepQW;
  double MWPCEnergy_MWPCEnergyE;
  double MWPCEnergy_MWPCEnergyW;
  double time_timeE;
  double time_timeW;

};

// Used for visualization, keeps the graph on screen.
TApplication plot_program("FADC_readin",0,0,0,0);

//-------------------------------------------------//
//------------ Start of Program -------------------//
//-------------------------------------------------//

int main(int argc, char* argv[])
{
  gRandom->SetSeed(0);
  gROOT->SetStyle("Pub");

  // set time window variables.
  double timeLowerEdgeE1 = -10;
  double timeUpperEdgeE1 = 140;
  double timeLowerEdgeW1 = -10;
  double timeUpperEdgeW1 = 3;
  double timeLowerEdgeE2 = -10;
  double timeUpperEdgeE2 = 5;
  double timeLowerEdgeW2 = -10;
  double timeUpperEdgeW2 = 140;

  double windowLower = 0;
  double windowUpper = 12;

  // define all the cuts we will use later.
  TCut basicCut = "(PID == 1 && badTimeFlag == 0)";
  TCut energyCut = "(Erecon_ee > 0 && Erecon_ee < 800)";
  TCut fiducialCut = "(((xE.center)*(xE.center) + (yE.center)*(yE.center) < 49*49) && ((xW.center)*(xW.center) + (yW.center)*(yW.center) < 49*49))";
  TCut simulationCut = "EdepQ.EdepQE > 0 && EdepQ.EdepQW > 0 && MWPCEnergy.MWPCEnergyE > 0 && MWPCEnergy.MWPCEnergyW > 0 && time.timeE < 200 && time.timeW < 200 && timeE < timeW";
  TCut inESTPCut = Form("newTDC2TimeE > -1 && newTDC2TimeE < 5 && newTDC2TimeW > %f", windowLower);
  TCut inWSTPCut = Form("newTDC2TimeW > -4 && newTDC2TimeW < 2 && newTDC2TimeE > %f", windowLower);
  TCut inESTPAndTimeWinCut = Form("newTDC2TimeE > -1 && newTDC2TimeE < 5 && newTDC2TimeW > %f && newTDC2TimeW < %f", windowLower, windowUpper);
  TCut inWSTPAndTimeWinCut = Form("newTDC2TimeW > -4 && newTDC2TimeW < 2 && newTDC2TimeE > %f && newTDC2TimeE < %f", windowLower, windowUpper);


  TCut time1STPCut = Form("newTDC2TimeE < %f && newTDC2TimeE > %f && newTDC2TimeW > %f && newTDC2TimeW < %f", timeUpperEdgeE1, timeLowerEdgeE1, timeLowerEdgeW1, timeUpperEdgeW1);
  TCut time2STPCut = Form("(newTDC2TimeE < %f && newTDC2TimeE > %f && newTDC2TimeW > %f && newTDC2TimeW < %f) || (newTDC2TimeE < %f && newTDC2TimeE > %f && newTDC2TimeW > %f && newTDC2TimeW < %f)",
			timeUpperEdgeE1, timeLowerEdgeE1, timeLowerEdgeW1, timeUpperEdgeW1,
			timeUpperEdgeE2, timeLowerEdgeE2, timeLowerEdgeW2, timeUpperEdgeW2);

  TCanvas *c2 = new TCanvas("c2","c2");
  c2->cd();

  TChain *mpm_sim_chain = new TChain("anaTree");
  mpm_sim_chain->Add("mpm_sim_betas_500EWFiles_2012-2013.root");

  TRandom3 *engine = new TRandom3(0);

  TH1D *hSim0 = new TH1D("hSim0sigma", "hSim", 280, 0, 140);
  hSim0->SetLineColor(6);
  TH1D *hSim2 = new TH1D("hSim2sigma", "hSim", 280, 0, 140);
  hSim2->SetLineColor(3);
  TH1D *hSim3 = new TH1D("hSim3sigma", "hSim", 280, 0, 140);
  hSim3->SetLineColor(2);

  SimEvent evt;
  mpm_sim_chain->GetBranch("Edep")->GetLeaf("EdepE")->SetAddress(&evt.Edep_EdepE);
  mpm_sim_chain->GetBranch("Edep")->GetLeaf("EdepW")->SetAddress(&evt.Edep_EdepW);
  mpm_sim_chain->GetBranch("EdepQ")->GetLeaf("EdepQE")->SetAddress(&evt.EdepQ_EdepQE);
  mpm_sim_chain->GetBranch("EdepQ")->GetLeaf("EdepQW")->SetAddress(&evt.EdepQ_EdepQW);
  mpm_sim_chain->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyE")->SetAddress(&evt.MWPCEnergy_MWPCEnergyE);
  mpm_sim_chain->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyW")->SetAddress(&evt.MWPCEnergy_MWPCEnergyW);
  mpm_sim_chain->GetBranch("time")->GetLeaf("timeE")->SetAddress(&evt.time_timeE);
  mpm_sim_chain->GetBranch("time")->GetLeaf("timeW")->SetAddress(&evt.time_timeW);

  for(unsigned int i = 0; i < 500000 /*mpm_sim_chain->GetEntries()*/; i++)
  {
    mpm_sim_chain->GetEntry(i);

    if(evt.EdepQ_EdepQE > 0 && evt.EdepQ_EdepQW > 0 && evt.MWPCEnergy_MWPCEnergyE > 0 && evt.MWPCEnergy_MWPCEnergyW > 0
        && evt.time_timeE < 200 && evt.time_timeW < 200 && evt.time_timeE < evt.time_timeW)
    {
      hSim0->Fill(evt.time_timeW - evt.time_timeE);
      hSim2->Fill(engine->Gaus(evt.time_timeW - evt.time_timeE, 2));
      hSim3->Fill(engine->Gaus(evt.time_timeW - evt.time_timeE, 3));
    }

    if(i % 100000 == 0)
    {
      cout << "Filling event " << i << "/" << mpm_sim_chain->GetEntries() << endl;
    }
  }

  hSim3->GetXaxis()->SetTitle("| T_{East} - T_{West} | (ns)");
  hSim3->GetXaxis()->CenterTitle();
  hSim3->GetYaxis()->SetTitle("Counts");
  hSim3->GetYaxis()->CenterTitle();
  hSim3->Draw();

  gPad->Update();

  TLine *l0 = new TLine(0, gPad->GetUymin(), 0, 1.5*gPad->GetUymax());
  l0->SetLineStyle(9);
  l0->Draw();

  TLine *l12 = new TLine(12, gPad->GetUymin(), 12, 1.5*gPad->GetUymax());
  l12->SetLineStyle(9);
  l12->Draw();

  double totalSimCounts = 0;
  for(int i = 0; i <= hSim3->GetNbinsX(); i++)
  {
    totalSimCounts = totalSimCounts + hSim3->GetBinContent(i);
  }

  // read in Brad's timing counts and plot them too
  TH1D* timingCounts = FillInTimingCounts("Acceptances/timing_644.txt");

  double totalBradSimCounts = 0;
  for(int i = 0; i <= timingCounts->GetNbinsX(); i++)
  {
    totalBradSimCounts = totalBradSimCounts + timingCounts->GetBinContent(i);
  }

  // this 0.30 is the 1% branching ratio to the dark decay channel to resolve the lifetime anomaly
  // relative to the population of Type 1 events that we are simulating above
  timingCounts->Scale(0.30 * totalSimCounts / totalBradSimCounts);

  timingCounts->SetLineColor(1);
  timingCounts->Draw("SAME");

  TLegend *l2 = new TLegend(0.55,0.7,0.85,0.85);
  l2->AddEntry(hSim3, "GEANT4 Type 1", "l");
  l2->AddEntry(timingCounts, "Monte Carlo e^{+}e^{-}", "l");
  l2->SetTextSize(0.05);
  l2->SetBorderSize(0);
  l2->Draw();
  gPad->SetLogy();

//  c2->Print("PresentationFigures/timingSpectraSimmed.pdf");
//  c2->Print("PresentationFigures/timingSpectraSimmed.eps");
//  c2->Print("PresentationFigures/timingSpectraSimmed.png");


  //------------------------------------------------------//
  // Creating another publication graph of TDC channels
  // spectrum, for PPNS conference proceedings paper.
  //------------------------------------------------------//




  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}


double SetPoissonErrors(int counts)
{
  double upperErrBar = 0;

  // all values taken from PDG paper rpp2017-rev-statistics.pdf
  // Using Table 40.3, for 95% upper one-sided limit
  if(counts == 0)
  {
    upperErrBar = 3;
  }
  else if(counts == 1)
  {
    upperErrBar = 4.74 - 1;
  }
  else if(counts == 2)
  {
    upperErrBar = 6.30 - 2;
  }
  else if(counts == 3)
  {
    upperErrBar = 7.75 - 3;
  }
  else if(counts == 4)
  {
    upperErrBar = 9.15 - 4;
  }
  else if(counts == 5)
  {
    upperErrBar = 10.51 - 5;
  }
  else if(counts == 6)
  {
    upperErrBar = 11.84 - 6;
  }
  else if(counts == 7)
  {
    upperErrBar = 13.15 - 7;
  }
  else if(counts == 8)
  {
    upperErrBar = 14.43 - 8;
  }
  else if(counts == 9)
  {
    upperErrBar = 15.71 - 9;
  }
  else if(counts == 10)
  {
    upperErrBar = 16.96 - 10;
  }
  else if(counts > 10)
  {
    upperErrBar = 2*sqrt(counts);
  }

  return (upperErrBar / 2.0);
}

TH1D* FillInTimingCounts(TString fileName)
{
  TH1D *h = new TH1D("timing", "timing", 280, 0, 140);

  double tmpTime = 0;
  double tmpCounts = 0;

  string buf;
  ifstream infile;
  cout << "The file being opened is: " << fileName << endl;
  infile.open(fileName.Data());

  if(!infile.is_open())
    cout << "Problem opening " << fileName << endl;

  while(getline(infile, buf))
  {
    istringstream bufstream(buf);
    if(!bufstream.eof())
    {
      bufstream >> tmpTime >> tmpCounts;

      if(tmpTime == 139.75)
      {
        continue;
      }
      h->SetBinContent(h->FindBin(tmpTime) , tmpCounts);
    }
  }

  return h;
}
