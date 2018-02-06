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

#define		path	"Type1_Octets_80-120"

using namespace std;

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
  int octNb = atoi(argv[1]);


  TChain *chain = new TChain("pass3");
  chain->Add(Form("%s/Octet_%i_type1.root", path, octNb));

  TH2D *h = new TH2D("TDC", Form("TDCE vs TDCW, Octet %i", octNb), 200, 0, 4000, 200, 0, 4000);
  h->GetXaxis()->SetTitle("TDCW");
  h->GetYaxis()->SetTitle("TDCE");

  chain->Draw("TDCE:TDCW >> TDC", "PID == 1 && badTimeFlag == 0", "COLZ");

  C->Print(Form("2Devents_TDCE_TDCW_Octet_%03i.pdf", octNb));

}
