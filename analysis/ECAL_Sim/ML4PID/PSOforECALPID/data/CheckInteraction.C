#include <TPaveStats.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <TFitResult.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TH1F.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystemFile.h"
#include "TGraph.h"
#include "TAxis.h"

void CheckInteraction(TString filename) {
  std::cout<<"Loading file : "<<filename<<std::endl;
  TFile *f = new TFile(filename+".root","update");
  TTree *T = (TTree*)f->Get("ntp"); 
  
  Float_t ecal_interaction = 0;
  T->SetBranchAddress("ecal_interaction",&ecal_interaction);
  Float_t hcal_interaction = 0;
  T->SetBranchAddress("hcal_interaction",&hcal_interaction);
  Float_t total_interaction = 0;
  T->SetBranchAddress("total_interaction",&total_interaction);

  Long64_t nentries = T->GetEntries();

  //Tables for check interactions
  int count_e = 0;
  int count_h = 0;
  int count_t = 0;
  for (int i = 0; i<nentries; i++) {
    T->GetEntry(i);
    if(ecal_interaction==1) {
      count_e += 1;
    }
    if(hcal_interaction==1) {
      count_h += 1;
    }
    if(total_interaction==1) {
      count_t += 1;
    }

  }
  cout<<"Total events: "<<nentries<<endl;
  cout<<"(ecal, hcal, total) = ("<<count_e<<", "<<count_h<<", "<<count_t<<")"<<endl;
}
