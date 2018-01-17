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

#define		TYPE	"REPLACEWITHTYPE"

using		 namespace std;

// plotting functions
void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xtitle, TString ytitle, TString command);

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
  C->Divide(2,1);

  // read in arguments.
  Int_t octNb = atoi(argv[1]);

  TFile *fE = TFile::Open(TString::Format("Octet_%i_TDCE_%s.root", octNb, TYPE));
  TFile *fW = TFile::Open(TString::Format("Octet_%i_TDCW_%s.root", octNb, TYPE));
  TH1D* hE = (TH1D*)fE->Get("Total hist");
  TH1D* hW = (TH1D*)fW->Get("Total hist");

  PlotHist(C, 1, 1, hE, Form("TDCE Octet %i, %s", octNb, TYPE), "Channels", "Counts", "");
  PlotHist(C, 1, 2, hW, Form("TDCW Octet %i, %s", octNb, TYPE), "Channels", "Counts", "");

  // Save our plot and print it out as a pdf.
  C -> Print(Form("outputHist_octet_%03i_%s.pdf", octNb, TYPE));
  cout << "-------------- End of Program ---------------" << endl;
//  plot_program.Run();

  return 0;
}

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xtitle, TString ytitle, TString command)
{
  C -> cd(canvasIndex);
  gPad->SetLogy();
  hPlot -> SetTitle(title);
  hPlot -> GetXaxis() -> SetTitle(xtitle);
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle(ytitle);
  hPlot -> GetYaxis() -> CenterTitle();
//  hPlot -> GetYaxis() -> SetRangeUser(0, 0.000004);

  if(styleIndex == 1)
  {
    hPlot -> SetFillColor(46);
    hPlot -> SetFillStyle(3004);
//    hPlot -> SetFillStyle(3001);
  }
  if(styleIndex == 2)
  {
    hPlot -> SetFillColor(38);
    hPlot -> SetFillStyle(3005);
//    hPlot -> SetFillStyle(3001);
  }
  if(styleIndex == 3)
  {
    hPlot -> SetFillColor(29);
//    hPlot -> SetFillStyle(3005);
    hPlot -> SetFillStyle(3001);
  }

  hPlot -> Draw(command);
  C -> Update();
}
