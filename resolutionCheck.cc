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
void FillInAcceptancesFromFile(TString fileName, double totalCounts, int flag);
void FillInCountsFromFile(TString fileName, TH1D* hToFill);
TGraph* CreateAcceptancesGraph(double timeLower, double timeUpper);

// some global access vectors to read in acceptance files
vector <double> t094;
vector <double> t144;
vector <double> t244;
vector <double> t344;
vector <double> t444;
vector <double> t544;
vector <double> t644;
vector <double> n094;
vector <double> n144;
vector <double> n244;
vector <double> n344;
vector <double> n444;
vector <double> n544;
vector <double> n644;


struct DataEvent
{
  double TDCE;
  double TDCW;
  int PID;
  int Type;
  int Side;
  double Erecon;
  double Erecon_ee;
  int badTimeFlag;

};

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

  TChain *bgchain = new TChain("pass3");
  TChain *fgchain = new TChain("pass3");

  bgchain->Add("TimeCalibrated_BGRuns_type1_fixed_v3.root");
  fgchain->Add("TimeCalibrated_FGRuns_type1_fixed_v3.root");

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

  // fourth canvas, energies with background subtraction and statistics propagated and acceptances from Brad F.
  TCanvas *c4 = new TCanvas("c4","c4");
  c4->Divide(2,1);

  TH1D *hbgErecon_withCuts = new TH1D("bgEreconfull", Form("BG Erecon: non-STP Events (%f, %f) ns", windowLower, windowUpper), 8, 44, 844);
  hbgErecon_withCuts->Sumw2();
  hbgErecon_withCuts->GetXaxis()->SetTitle("Erecon_ee (keV)");
  TH1D *hfgErecon_withCuts = new TH1D("fgEreconfull", Form("FG Erecon: non-STP Events (%f, %f) ns", windowLower, windowUpper), 8, 44, 844);
  hfgErecon_withCuts->Sumw2();
  hfgErecon_withCuts->GetXaxis()->SetTitle("Erecon_ee (keV)");


  c4->cd(1);
  bgchain->Draw("Erecon_ee >> bgEreconfull", basicCut && fiducialCut && (inESTPAndTimeWinCut || inWSTPAndTimeWinCut) && energyCut);
  c4->cd(2);
  fgchain->Draw("Erecon_ee >> fgEreconfull", basicCut && fiducialCut && (inESTPAndTimeWinCut || inWSTPAndTimeWinCut) && energyCut);

  cout << "Setting Poisson error bars..." << endl;
  for(int i = 0; i < hfgErecon_withCuts->GetNbinsX(); i++)
  {
    hbgErecon_withCuts->SetBinError(i, SetPoissonErrors(hbgErecon_withCuts->GetBinContent(i)));
    hfgErecon_withCuts->SetBinError(i, SetPoissonErrors(hfgErecon_withCuts->GetBinContent(i)));

//    cout << "At bin i = " << i << " BG counts = " << hbgErecon_withCuts->GetBinContent(i) << ", FG counts = " << hfgErecon_withCuts->GetBinContent(i) << endl;
  }


//  c4->Print("4_BG_FG_histograms_withErrors.pdf");
//  c4->Print("4_BG_FG_histograms_withErrors.eps");

  // fifth canvas, does background subtraction and acceptances, because ROOT is tough to work with.
  TCanvas *c5 = new TCanvas("c5","c5");
  c5->cd();

  gStyle->SetOptStat(0000);
  TH1D *hErecon_bgSub = new TH1D("fullCuts", "BG subtracted Erecon_ee", 8, 44, 844);
  hErecon_bgSub->SetLineColor(2);
  hErecon_bgSub->Sumw2();
  hErecon_bgSub->Add(hfgErecon_withCuts, hbgErecon_withCuts, 1, -5.07);	// 5.07 comes from 860262/169717 live time ratio

  double totalEntries = 0;
  for(int i = 0; i < hErecon_bgSub->GetNbinsX(); i++)
  {
    totalEntries = totalEntries + hErecon_bgSub->GetBinContent(i);
  }
  hErecon_bgSub->SetEntries(totalEntries);

  // sixth canvas, for acceptances
  // This is the old stuff which does the proper bin centers
  TCanvas *c6 = new TCanvas("c6","c6");
  FillInAcceptancesFromFile("m094002MeV_4502123TotalCounts.dat", 4502123, 94);
  FillInAcceptancesFromFile("m144002MeV_455639TotalCounts.dat", 455639, 144);
  FillInAcceptancesFromFile("m244002MeV_460822TotalCounts.dat", 460822, 244);
  FillInAcceptancesFromFile("m344002MeV_466564TotalCounts.dat", 466564, 344);
  FillInAcceptancesFromFile("m444002MeV_471023TotalCounts.dat", 471023, 444);
  FillInAcceptancesFromFile("m544002MeV_475630TotalCounts.dat", 475630, 544);
  FillInAcceptancesFromFile("m644002MeV_478192TotalCounts.dat", 478192, 644);

  // now we are doing acceptances for a shifted bin center
  vector <double> acceptanceShift;
  vector <double> energyShift;
  energyShift.push_back(94);
  acceptanceShift.push_back(0.000592);
  energyShift.push_back(194);
  acceptanceShift.push_back(0.0501791);
  energyShift.push_back(294);
  acceptanceShift.push_back(0.1132812);
  energyShift.push_back(394);
  acceptanceShift.push_back(0.1592465);
  energyShift.push_back(494);
  acceptanceShift.push_back(0.1891385);
  energyShift.push_back(594);
  acceptanceShift.push_back(0.208866);
  energyShift.push_back(694);
  acceptanceShift.push_back(0.222);

//  TGraph *gAccept = CreateAcceptancesGraph(0.0, 12.0);

  TGraph *gAccept = new TGraph(energyShift.size(), &(energyShift[0]), &(acceptanceShift[0]));

  double xAcceptKin = 0;
  double yAcceptKin = 0;
  // adjust the kinematic acceptance for G4 sim, for total acceptance
  for(int i = 0; i < gAccept->GetN(); i++)
  {
    gAccept->GetPoint(i, xAcceptKin, yAcceptKin);
    gAccept->SetPoint(i, xAcceptKin, yAcceptKin*0.845);
  }
  gAccept->SetMarkerStyle(21);
  gAccept->SetMarkerColor(1);
  gAccept->SetLineColor(1);
  gAccept->GetHistogram()->SetTitle("");
  gAccept->GetHistogram()->GetXaxis()->SetTitle("E_{e^{+}e^{-}} (keV)");
  gAccept->GetHistogram()->GetXaxis()->CenterTitle();
  gAccept->GetHistogram()->GetYaxis()->SetTitle("Fraction Accepted");
  gAccept->GetHistogram()->GetYaxis()->CenterTitle();
  gAccept->GetHistogram()->GetXaxis()->SetRangeUser(0, 800);
  gAccept->GetHistogram()->GetYaxis()->SetRangeUser(0, 0.22);
  gAccept->Draw("AC");

  c6->Print("Figures/6_Acceptance.pdf");
  c6->Print("Figures/6_Acceptance.eps");

  double G4x[4] = {200, 321, 480, 640};
  double G4y[4] = {82560.0/96558.0, 82212.0/96721.0, 80955.0/96271.0, 80283.0/96067.0};
  TGraph *gG4 = new TGraph(4, G4x, G4y);

  // this is the previous, centered at the right bins
  double gEx_old[7] = {94, 144, 244, 344, 444, 544, 644};
  // 50.0 KeV is our bin width on either side of the center
  double sigmaEKeV_old[7] = {50.0/15.33, 50.0/18.97, 50.0/24.70, 50.0/29.32, 50.0/33.32, 50.0/36.88, 50.0/40.12};
  // calculated by taking a Gaussian probability calculator and plugging in (denom) as SD, range from -50 to 50
  double gEy_old[7] = {0.9989, 0.9916, 0.9571, 0.9119, 0.8665, 0.8248, 0.7873};



  // this resolution is for a shifted bin center
  double gEx[7] = {94, 194, 294, 394, 494, 594, 694};
  // 50.0 KeV is our bin width on either side of the center
  double sigmaEKeV[7] = {50.0/15.33, 50.0/22.02, 50.0/27.11, 50.0/31.38, 50.0/35.14, 50.0/38.53, 50.0/41.65};
  // calculated by taking a Gaussian probability calculator and plugging in (denom) as SD, range from -50 to 50
  double gEy[7] = {0.9989, 0.9768, 0.9349, 0.8889, 0.8452, 0.8056, 0.77};


  TGraph *gEResolution = new TGraph(7, gEx, gEy);

  c5->cd();
  double xAccept = 0;
  double yAccept = 0;
  double xEResolution = 0;
  double yEResolution = 0;

  double finalEntries = 0;

  TH1D* hFinalNumbers = new TH1D("hFinal", "Final: Cuts and scaled for Acceptances", 8, 44, 844);
  hFinalNumbers->SetStats(kFALSE);
  hFinalNumbers->SetTitle("");
  hFinalNumbers->GetXaxis()->SetTitle("E_{e^{+}e^{-}} (keV)");
  hFinalNumbers->GetYaxis()->SetTitle("Counts");
  hFinalNumbers->GetYaxis()->SetRangeUser(0, 1000);
  hFinalNumbers->SetLineColor(1);

  for(int i = 0; i < hFinalNumbers->GetNbinsX(); i++)
  {
    gAccept->GetPoint(i, xAccept, yAccept);
    cout << "At i = " << i << " we have xAccept = " << xAccept << " and yAccept = " << yAccept << endl;

    gEResolution->GetPoint(i, xEResolution, yEResolution);
    cout << "At i = 0" << i << " we have xEResolution = " << xEResolution << " and yEResolution = " << yEResolution << endl;

    hFinalNumbers->SetBinContent(i, (double)hErecon_bgSub->GetBinContent(i) / (yAccept*yEResolution));
    hFinalNumbers->SetBinError(i, (double)hErecon_bgSub->GetBinError(i) / (yAccept*yEResolution));

    finalEntries = finalEntries + ((double)hErecon_bgSub->GetBinContent(i) / (yAccept*yEResolution));
  }

  hFinalNumbers->SetEntries(finalEntries);

  hFinalNumbers->Draw();
  hErecon_bgSub->SetBinContent(hErecon_bgSub->GetNbinsX(), 0);
  hErecon_bgSub->SetBinError(hErecon_bgSub->GetNbinsX(), 0);
  hErecon_bgSub->Draw("SAME");


  cout << "Total counts in histogram after final acceptances propagated: " << hFinalNumbers->GetEntries() << endl;
  cout << "Total counts in histogram after background subtraction: " << hErecon_bgSub->GetEntries() << endl;

  for(int i = 0; i < hFinalNumbers->GetNbinsX(); i++)
  {
    cout << "hFinalNumbers at x = " << hFinalNumbers->GetBinCenter(i)
	 << " has counts = " << hFinalNumbers->GetBinContent(i)
	 << " and associated symmetric error = " << hFinalNumbers->GetBinError(i) << endl;
  }

//  c5->Print("Figures/5_BGSub_WithAcceptances_finalNumbers.pdf");
//  c5->Print("Figures/5_BGSub_WithAcceptances_finalNumbers.eps");

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

void FillInCountsFromFile(TString fileName, TH1D* hToFill)
{
  double tmpEnergy;
  double tmpCounts;

  string buf;
  ifstream infile;
  cout << "The file being opened is: " << fileName << endl;
  infile.open(Form("Acceptances/%s", fileName.Data()));

  if(!infile.is_open())
    cout << "Problem opening " << fileName << endl;

  while(getline(infile, buf))
  {
    istringstream bufstream(buf);
    if(!bufstream.eof())
    {
      bufstream >> tmpEnergy >> tmpCounts;

      for(int i = 0; i < tmpCounts; i++)
      {
        hToFill->Fill(tmpEnergy*1000);
      }
    }
  }

  cout << "Done filling counts into signal histogram." << endl;
}


void FillInAcceptancesFromFile(TString fileName, double totalCounts, int flag)
{
  vector <double> time;
  vector <double> counts;

  double tmpTime;
  double tmpCounts;

  string buf;
  ifstream infile;
  cout << "The file being opened is: " << fileName << endl;
  infile.open(Form("Acceptances/%s", fileName.Data()));

  if(!infile.is_open())
    cout << "Problem opening " << fileName << endl;

  while(getline(infile, buf))
  {
    istringstream bufstream(buf);
    if(!bufstream.eof())
    {
      bufstream >> tmpTime >> tmpCounts;

      time.push_back(tmpTime);
      counts.push_back(tmpCounts / (totalCounts));

    }
  }

  // a bunch of different cases because I don't want to use addreses
  if(flag == 94)
  {
    t094 = time;
    n094 = counts;
  }
  else if(flag == 144)
  {
    t144 = time;
    n144 = counts;
  }
  else if(flag == 244)
  {
    t244 = time;
    n244 = counts;
  }
  else if(flag == 344)
  {
    t344 = time;
    n344 = counts;
  }
  else if(flag == 444)
  {
    t444 = time;
    n444 = counts;
  }
  else if(flag == 544)
  {
    t544 = time;
    n544 = counts;
  }
  else if(flag == 644)
  {
    t644 = time;
    n644 = counts;
  }
  else
  {
    cout << "Flag makes no sense, not sure which vectors to fill!" << endl;
  }
}

TGraph* CreateAcceptancesGraph(double timeLower, double timeUpper)
{
  int nPoints = 0;
  vector <double> energy;
  vector <double> acceptance;

  energy.push_back(0);
  acceptance.push_back(0);
  nPoints = nPoints + 1;

  double accept094 = 0;
  for(unsigned int i = 0; i < t094.size(); i++)
  {
    if(t094[i] < timeLower || t094[i] > timeUpper)
    {
      continue;
    }
    accept094 = accept094 + n094[i];
  }
  energy.push_back(94);
  acceptance.push_back(accept094);
  nPoints = nPoints + 1;

  double accept144 = 0;
  for(unsigned int i = 0; i < t144.size(); i++)
  {
    if(t144[i] < timeLower || t144[i] > timeUpper)
    {
      continue;
    }
    accept144 = accept144 + n144[i];
  }
  energy.push_back(144);
  acceptance.push_back(accept144);
  nPoints = nPoints + 1;

  double accept244 = 0;
  for(unsigned int i = 0; i < t244.size(); i++)
  {
    if(t244[i] < timeLower || t244[i] > timeUpper)
    {
      continue;
    }
    accept244 = accept244 + n244[i];
  }
  energy.push_back(244);
  acceptance.push_back(accept244);
  nPoints = nPoints + 1;

  double accept344 = 0;
  for(unsigned int i = 0; i < t344.size(); i++)
  {
    if(t344[i] < timeLower || t344[i] > timeUpper)
    {
      continue;
    }
    accept344 = accept344 + n344[i];
  }
  energy.push_back(344);
  acceptance.push_back(accept344);
  nPoints = nPoints + 1;

  double accept444 = 0;
  for(unsigned int i = 0; i < t444.size(); i++)
  {
    if(t444[i] < timeLower || t444[i] > timeUpper)
    {
      continue;
    }
    accept444 = accept444 + n444[i];
  }
  energy.push_back(444);
  acceptance.push_back(accept444);
  nPoints = nPoints + 1;

  double accept544 = 0;
  for(unsigned int i = 0; i < t544.size(); i++)
  {
    if(t544[i] < timeLower || t544[i] > timeUpper)
    {
      continue;
    }
    accept544 = accept544 + n544[i];
  }
  energy.push_back(544);
  acceptance.push_back(accept544);
  nPoints = nPoints + 1;

  double accept644 = 0;
  for(unsigned int i = 0; i < t644.size(); i++)
  {
    if(t644[i] < timeLower || t644[i] > timeUpper)
    {
      continue;
    }
    accept644 = accept644 + n644[i];
  }
  energy.push_back(644);
  acceptance.push_back(accept644);
  nPoints = nPoints + 1;

  energy.push_back(694);
  acceptance.push_back(0.222);
  nPoints = nPoints + 1;

  energy.push_back(744);
  acceptance.push_back(0.229);
  nPoints = nPoints + 1;

  energy.push_back(794);
  acceptance.push_back(0.233);
  nPoints = nPoints + 1;


  for(int i = 0; i < nPoints; i++)
  {
    cout << "Points being fed into the graph: (" << energy[i] << ", " << acceptance[i] << ")" << endl;

  }

  TGraph *g = new TGraph(nPoints, &(energy[0]), &(acceptance[0]));

  return g;
}
