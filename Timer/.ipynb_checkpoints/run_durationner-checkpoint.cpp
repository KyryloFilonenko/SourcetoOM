#include <string>
#include "TCanvas.h"
#include "TPaveText.h"
#include <iostream>
#include <memory>
using namespace RooFit;
#define _USE_MATH_DEFINES
#include <fstream>
#include <sstream>
#include <vector>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TStyle.h>
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TROOT.h"
#include <TText.h>
#include <TLatex.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TMinuit.h>



long long get_run_unixtime(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;

    while (std::getline(file, line)) {

        if (line.find("run.run_unixtime_ms=") == 0) {
            std::string value = line.substr(line.find("=") + 1);
            return std::stoll(value); 
        }
    }

    return -1;
}


void durationner(int run_number){

  std::vector<long> *timestamp = new std::vector<long>;
  long t0 = 0;
  long t1 = 0;
  int run_duration = 0;
  TFile *file = new TFile(Form("/sps/nemo/snemo/snemo_data/reco_data/UDD/delta-tdc-10us-v3/snemo_run-%d_udd.root", run_number), "READ");
  if (!file || file->IsZombie()) return;  
  TTree *tree = (TTree *)file->Get("SimData");

  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("digicalo.timestamp", 1);
  tree->SetBranchAddress("digicalo.timestamp", &timestamp);

  tree->GetEntry(0);
  t0 = timestamp->at(0);

  tree->GetEntry(tree->GetEntries()-1);
  t1 = timestamp->at(timestamp->size()-1);

  run_duration = (t1 - t0)*6.25e-9;

  long long runstart = get_run_unixtime(Form("/sps/nemo/snemo/snemo_data/raw_data/CBD/run-%d/snemo_crate-0_run-%d.log", run_number, run_number));

  cout << run_number << " " << runstart << " " << run_duration << std::endl;
  
  std::ofstream out("/sps/nemo/scratch/kfilonen/Calibration/SourcetoOM/Timer/run_durations_betabeta.txt", std::ios::app);
  out << run_number << " " << runstart << " " << run_duration << std::endl;
  out.close();


  file->Close();

  
  }


void sort_and_replace_input() { // to sort runs and supress double entries

    std::ifstream in("/sps/nemo/scratch/kfilonen/Calibration/SourcetoOM/Timer/run_durations.txt");
    std::ofstream tmp("/sps/nemo/scratch/kfilonen/Calibration/SourcetoOM/Timer/run_durations_tmp.txt");

    std::string header;
    std::getline(in, header);
    tmp << header << std::endl;

    struct RunInfo {
      long long unixtime;
      double duration;
    };
    
    std::map<int, RunInfo> runs;

    int run_number;
    double duration;
    long long unixtime;
      
    while (in >> run_number >> unixtime >> duration) {
      runs[run_number] = {unixtime, duration};
    }
    
    for (const auto& [rn, dur] : runs) {
      tmp << rn << " " << dur.unixtime << " " << dur.duration <<  std::endl;
    }

    in.close();
    tmp.close();
    
    std::remove("/sps/nemo/scratch/kfilonen/Calibration/SourcetoOM/Timer/run_durations.txt");
    std::rename("/sps/nemo/scratch/kfilonen/Calibration/SourcetoOM/Timer/run_durations_tmp.txt", "/sps/nemo/scratch/kfilonen/Calibration/SourcetoOM/Timer/run_durations.txt");

}

void run_durationner(int run_number){
  durationner(run_number);
  sort_and_replace_input();
}
