//
//Written by: B. Olaizola and A. Illana 2022-05-29
//Experiment: IS622
//Root version: 6.14/04
//This script has a TChain, so you can give a list of files and it will loop through all of them
//
/////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// To compile:
//
// g++ ids_histo.C -Wall `root-config --cflags` `root-config --glibs` -o ids_histo
//
// ./ids_histo Run362.root
//
/////////////////////////////////////////////////////////////////////////////////////////////////

// ToDo:
//  - Adding a quick addback function
//  - Improving the input and output files


//C,C++ Libraries
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>

// ROOT libraries
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TColor.h>
#include <TROOT.h>
#include "TChain.h"
#include <iomanip>

using namespace std;

// NUMBER OF HISTOGRAMS FOR ANGULAR CORRELATIONS, GIVEN BY ANGULAR RESOLUTION

  double range = 180.;
  double angular_resolution = 12.9; // Follow IS709 Crystal pairing spreadsheet
  const int angular_binning = round(range/angular_resolution);

// FUNCTION TO GET THE CENTER OF THE ANGULAR BINNING
double get_bin_center(double angle) {
    if (angle < 0) angle = 0;
    if (angle > 180) angle = 180;

    int bin_index = static_cast<int>(angle / angular_resolution);
    if (bin_index >= angular_binning) bin_index = angular_binning - 1;

    double bin_center = (bin_index + 0.5) * angular_resolution;
    return bin_center;
}

void PrintBar(float _progress){

  int barWidth = 70;
  cout << "[";
  int pos = barWidth * (_progress/100.);
  for (int i = 0; i < barWidth; ++i) {
      if (i < pos) cout << "=";
      else if (i == pos) cout << ">";
      else cout << " ";
  }
  if (_progress > 100.) _progress = 100.;
  cout << "] " << int(_progress) << " %\r";
  cout.flush();

  return;
}


int main(int argc,char* argv[]){
  if(argc == 1) {
    std::cout<<" Usage: "<<argv[0]<<" <analysis tree file(s)>"<<std::endl;
    return 1;
  }

  TChain* tree = new TChain("ids");

  for(int i = 1; i < argc; ++i) {
    std::cout<<"'"<<argv[i]<<"': added "<<tree->Add(argv[i])<<" files"<<std::endl;
  }
  string inputFileName = argv[1];

  string outputFileName = inputFileName.substr(0, inputFileName.size()-5);
  outputFileName = outputFileName + "_histograms.root" ;
  cout << "" <<  endl;
  cout << " Output file: " << outputFileName <<  endl;
  cout << "" <<  endl;

///////////////////Number of each type of detectors////////////////

  //const int num_clov_crys = 24; //This is number of crystals, not clovers, so usually crystals=4*clovers
  const int num_labr      = 2;  //This is number of LaBr crystals
  const int num_tac       = 4;  //This is number of tacs
  const int num_beta      = 3;  //This is number of beta plastics
  const int num_veto      = 4;  //This is number of VETO plastics

// ============================== 
// DEFINITION OF THE CRYSTALS
// ==============================

    const int n_gantries = 5;
    const int n_clovers_per_gantry = 2;
    const int n_colors = 4;

    const std::string gantries[n_gantries] = { "g1", "g2", "g3", "g4", "g5" };
    const char clover_labels[n_clovers_per_gantry] = { 'A', 'B' };
    const std::string colors[n_colors] = { "red", "green", "blue", "white" };

    const int num_clov_crys = n_gantries * n_clovers_per_gantry * n_colors;

    std::vector<std::string> clover(num_clov_crys);
    std::vector<std::string> color(num_clov_crys);

    int idx = 0;
    for(int i=0; i<n_gantries; i++) {
        for(int j=0; j<n_clovers_per_gantry; j++) {
            std::string name = gantries[i] + "_" + clover_labels[j];
            for(int c=0; c<n_colors; c++) {
                clover[idx] = name;
                color[idx] = colors[c];
                idx++;
            }
        }
    }


// ==============================
// CREATION OF THE MAP
// ==============================

// Reads the .txt angle pairing file and creates the map
    std::map<std::string, double> angulos;
    std::ifstream archivo("angles_and_combinations_both_10_clovers.txt");

    std::string linea;
    std::getline(archivo, linea); // skip header

    while (std::getline(archivo, linea)) {
        size_t coma_pos = linea.rfind(',');
        if (coma_pos == std::string::npos) continue;
        std::string clave = linea.substr(0, coma_pos);
        std::string angulo_str = linea.substr(coma_pos + 1);

        try {
            double angulo = std::stod(angulo_str);
            angulos[clave] = angulo;
        } catch (...) {
            continue;
        }
    }
    archivo.close();

  //Config Index
  string TacName[4]  = {"SiPM1", "SiPM2", "LaBrLaBr","PMTLaBr"};
  string BetaName[3] = {"PMT", "SiPM1", "SiPM2"};
  string VetoName[6] = {"Clov0", "Clov1", "Clov2", "Clov3", "Clov4", "Clov5"};

  int tac_labr      = 2;
  int tac_betapmt   = 3;
  int tac_betasipm1 = 0;
  int tac_betasipm2 = 1;
  int beta_pmt   = 0;
  int beta_sipm1 = 1;
  int beta_sipm2 = 2;

  ///////////////////////Time differences, you probably want to edit/////////////////////////
  //Clover-clover time windows
  static double clov_t_coinc            = 60.;//True time coincidence
  static double clov_rand_coinc[2]      = {100., 500.};//Time-random coincidences
  //LaBr-LaBr time windows
  static double labr_t_coinc            = 60.;
  //Clov-LaBr time windows
  static double clov_labr_t_coinc[2]    = {1040., 600.};
  static double clov_beta_t_coinc[2]    = {1040., 600.};
  //Beta Energy windows
  static double betapmt_energy_gate[2]  = {20., 4000.};
  static double betasipm_energy_gate[2] = {20., 40000.};

  //Declaration of leaves types
  //Int_t           MULT;
  ULong64_t       TIME_REF;
  ULong64_t       TIMESTAMP;
  Double_t        E_Clov[num_clov_crys];
  Int_t           T_Clov[num_clov_crys];
  Int_t           M_Clov;
  Double_t        E_Beta[num_beta+num_veto];
  Int_t           T_Beta[num_beta+num_veto];
  Int_t           M_Beta;
  Double_t        E_LaBr[num_labr];
  Int_t           T_LaBr[num_labr];
  Int_t           M_LaBr;
  Double_t        E_TAC[num_tac];
  Int_t           T_TAC[num_tac];
  Int_t           M_TAC;

  // Set branch addresses.
  tree->SetBranchAddress("Time_vs_ref",&TIME_REF);
  tree->SetBranchAddress("Timestamp",&TIMESTAMP);

  bool gotClov = false;
  if(tree->FindBranch("Energy_Clov") == nullptr) { // We check to see if we have a clover branch in the analysis tree
    gotClov = false;
  } else {
    tree->SetBranchAddress("Energy_Clov",&E_Clov);
    tree->SetBranchAddress("Time_Clov",&T_Clov);
    tree->SetBranchAddress("Mult_Clov",&M_Clov);
    gotClov = true;
    cout<<" Setting the Clover branch "<<endl;
  }
/*
  bool gotLaBr = false;
  if(tree->FindBranch("Energy_LaBr") == nullptr) { // We check to see if we hbiuave a LaBr branch in the analysis tree
    gotLaBr = false;
  } else {
    tree->SetBranchAddress("Energy_LaBr",&E_LaBr);
    tree->SetBranchAddress("Time_LaBr",&T_LaBr);
    tree->SetBranchAddress("Mult_LaBr",&M_LaBr);
    gotLaBr = true;
    cout<<" Setting the LaBr branch "<<endl;
  }

  bool gotTAC = false;
  if(tree->FindBranch("Energy_TAC") == nullptr) { // We check to see if we have a TAC branch in the analysis tree
    gotTAC = false;
  } else {
    tree->SetBranchAddress("Energy_TAC",&E_TAC);
    tree->SetBranchAddress("Time_TAC",&T_TAC);
    tree->SetBranchAddress("Mult_TAC",&M_TAC);
    gotTAC = true;
    cout<<" Setting the TAC branch "<<endl;
  }
*/
  bool gotBeta = false;
  if(tree->FindBranch("Energy_Beta") == nullptr) { // We check to see if we have a beta branch in the analysis tree
      gotBeta = false;
  } else {
    tree->SetBranchAddress("Energy_Beta",&E_Beta);
    tree->SetBranchAddress("Time_Beta",&T_Beta);
    tree->SetBranchAddress("Mult_Beta",&M_Beta);
    gotBeta = true;
    cout<<" Setting the beta branch "<<endl;
  }

  //General variables
  double time_diff;
  double time_ref_c = 1200.;     // 101Sr 200 ms, 101Y 600 ms, 101Zr 1200 ms
  int event_separation = 10;


  //Here we create the histograms. We will add them to a TList, so they are all written at the same time.
  //Clover histograms
  auto* list_clov = new TList; //list of the clov histograms
  TH1F* clov_time                 = new TH1F("clov_time","Clover-clover time difference; Time 4 ns/Chan",1000,-500.,500.); list_clov->Add(clov_time); /////
  TH1F* clov_singles              = new TH1F("clov_singles","Clover singles; E_{gamma} [keV]; Counts/0.25 keV",16000,0.,4000.); list_clov->Add(clov_singles); /////
  TH1F* clov_singles_betagated    = new TH1F("clov_singles_betagated","Clover singles beta gated; E_{gamma} [keV]; Counts/0.25 keV",16000,0.,4000.); list_clov->Add(clov_singles_betagated);
  TH2F* clov_gg                   = new TH2F("clov_gg","Clover-Clover with time random subtracted; E_{gamma} [1.0 keV/Chan]; E_{gamma} [1.0 keV/Chan]",4000,0.,4000.,4000,0.,4000.); list_clov->Add(clov_gg);/////////
  TH2F* clov_gg_Notrand           = new TH2F("clov_gg_Notrand","Clover-Clover without time-random coincidences subtraction; E_{gamma} [1.0 keV/Chan]; E_{gamma} [1.0 keV/Chan]",4000,0.,4000.,4000,0.,4000.); list_clov->Add(clov_gg_Notrand);//////
  TH2F* clov_gg_trand             = new TH2F("clov_gg_trand","Clover-Clover time-random coincidences; E_{gamma} [1.0 keV/Chan]; E_{gamma} [1.0 keV/Chan]",4000,0.,4000.,4000,0.,4000.); list_clov->Add(clov_gg_trand);///////
  TH2F* clov_gg_betagated         = new TH2F("clov_gg_betagated","Clover-Clover with time random subtracted beta gated; E_{gamma} [1.0 keV/Chan]; E_{gamma} [1.0 keV/Chan]",4000,0.,4000.,4000,0.,4000.); list_clov->Add(clov_gg_betagated);//////
  TH2F* clov_gg_Notrand_betagated = new TH2F("clov_gg_Notrand_betagated","Clover-Clover without time-random coincidences subtraction beta gated;  E_{gamma} [1.0 keV/Chan]; E_{gamma} [1.0 keV/Chan]",4000,0.,4000.,4000,0.,4000.); list_clov->Add(clov_gg_Notrand_betagated);//////
  TH2F* clov_gg_trand_betagated   = new TH2F("clov_gg_trand_betagated","Clover-Clover time-random coincidences beta gated;  E_{gamma} [1.0 keV/Chan]; E_{gamma} [1.0 keV/Chan]",4000,0.,4000.,4000,0.,4000.); list_clov->Add(clov_gg_trand_betagated); ///////
  TH2F* clov_vs_TRef              = new TH2F("clov_vs_TRef","TimeRef vs Clovers; Time [10 ms]; E_{gamma} [keV]",4000,0,4000, 4000,0.,4000.); list_clov->Add(clov_vs_TRef);
  TH2F* clov_vs_TRef_betagated    = new TH2F("clov_vs_TRef_betagated","TimeRef vs Clovers beta gated; Time [10 ms]; E_{gamma} [keV]",4000,0,4000, 4000,0.,4000.); list_clov->Add(clov_vs_TRef_betagated); ////// binning
  TH2F* clov_Id_vs_Ener           = new TH2F("clov_Id_vs_Ener", "Clover ID vs Energy", 24, 0, 24, 8000, 0, 2000); list_clov->Add(clov_Id_vs_Ener);
  TH2F* histo                     = new TH2F("histo", "histo", 4000,0,4000, 4000,0.,4000.); list_clov->Add(histo);   //clov_vs_TRef_betagated (gg filter)
  TH2F* histo2                    = new TH2F("histo2", "histo2", 4000,0,4000, 4000,0.,4000.); list_clov->Add(histo2);   //clov_vs_TRef_betagated (gg filter)
  TH2F* histo_bkg                 = new TH2F("histo_bkg", "histo_bkg", 4000,0,4000, 4000,0.,4000.); list_clov->Add(histo_bkg); // background (gg filter)


    //LaBr histograms
  auto* list_labr = new TList; //list of the LaBr histograms
  TH1F* labr_time              = new TH1F("labr_time","LaBr-LaBr time difference; Time 4 ns/Chan",600, -300., 300.); list_labr->Add(labr_time);
  TH1F* labr_singles           = new TH1F("labr_singles","LaBr singles; E_{gamma} [keV]; Counts/2.0 keV",2000,0.,4000.); list_labr->Add(labr_singles);
  TH1F* labr_singles_betagated = new TH1F("labr_singles_betagated","LaBr singles beta gated; E_{gamma} [keV]; Counts/2.0 keV",2000,0.,4000.); list_labr->Add(labr_singles_betagated);  
  TH2F* labr_gg       = new TH2F("labr_gg","LaBr-LaBr withOUT time random subtraction; E_{gamma} [keV/Chan]; E_{gamma} [keV/Chan]",1000,0.,2000.,1000,0.,2000.); list_labr->Add(labr_gg);
  TH2F* labr_gg_betagated      = new TH2F("labr_gg_betagated","LaBr-LaBr withOUT time random subtraction and beta gated; E_{gamma} [keV/Chan]; E_{gamma} [keV/Chan]",1000,0.,2000.,1000,0.,2000.); list_labr->Add(labr_gg_betagated);

  //Clov vs LaBr histograms
  auto* list_clovlabr = new TList; //list of the LaBr histograms
  TH1F* clov_labr_time         = new TH1F("Clov_labr_time","Clover-LaBr time difference; Time 4 ns/Chan",1000,-500.,500.); list_labr->Add(clov_labr_time);
  TH2F* gg_clov_labr           = new TH2F("gg_clov_labr","Clover-LaBr; E_{gamma} Clov [1.0 keV/Chan]; E_{gamma} LaBr [2.0 keV/Chan]",2000,0.,2000.,1000,0.,2000.); list_clovlabr->Add(gg_clov_labr);
  TH2F* gg_clov_labr_betagated = new TH2F("gg_clov_labr_betagated","Clover-LaBr beta gated; E_{gamma} Clov [1.0 keV/Chan]; E_{gamma} LaBr [2.0 keV/Chan]",2000,0.,2000.,1000,0.,2000.); list_clovlabr->Add(gg_clov_labr);

  TH2F* beta_LaBr_clov = new TH2F("beta_LaBr_clov","Beta_LaBr_Clover; Time [ns]; E_{gamma} Clov [1.0 keV/Chan]; Counts/10 ps",500,0.,500.,2000,0.,2000.); list_clovlabr->Add(beta_LaBr_clov);
  TH2F* beta_LaBr_clov_bkg = new TH2F("beta_LaBr_clov_bkg","Beta_LaBr_Clover_bkg; Time [ns]; E_{gamma} Clov [1.0 keV/Chan]; Counts/10 ps",500,0.,500.,2000,0.,2000.); list_clovlabr->Add(beta_LaBr_clov_bkg);

  //Clov vs Beta histograms

  TH1F* clov_beta_time         = new TH1F("Clov_beta_time","Clover-Beta time difference; Time 4 ns/Chan",1000,-500.,500.); list_clov->Add(clov_beta_time);
  TH2F* clov_timebeta = new TH2F("clov_timebeta","START: beta - STOP: Clov; Time 4 ns/Chan; E_{gamma} [keV/Chan]; Counts/10 ps",500,0.,500.,2000,0.,2000.); list_clov->Add(clov_timebeta);
  TH2F* beta_timeclov = new TH2F("beta_timeclov","START: beta - STOP: Clov; Time 4 ns/Chan; E_{gamma} [keV/Chan]; Counts/10 ps",500,0.,500.,2000,0.,4000.); list_clov->Add(beta_timeclov);

  TH2F* beta_timeclov_pmt = new TH2F("beta_timeclov_pmt","START: beta_pmt - STOP: Clov; Time 4 ns/Chan; E_{gamma} [keV/Chan]; Counts/10 ps",1000,0.,1000.,4000,0.,8000.); list_clov->Add(beta_timeclov_pmt);
  TH2F* beta_timeclov_sipm1 = new TH2F("beta_timeclov_sipm1","START: beta_sipm1 - STOP: Clov; Time 4 ns/Chan; E_{gamma} [keV/Chan]; Counts/10 ps",1000,0.,1000.,20000,0.,80000.); list_clov->Add(beta_timeclov_sipm1);
  TH2F* beta_timeclov_sipm2 = new TH2F("beta_timeclov_sipm2","START: beta_sipm2 - STOP: Clov; Time 4 ns/Chan; E_{gamma} [keV/Chan]; Counts/10 ps",1000,0.,1000.,20000,0.,80000.); list_clov->Add(beta_timeclov_sipm2);


  TH1F* Ebeta_pmt              = new TH1F("Ebeta_pmt","Ebetah_pmt; E_{gamma} [keV]; Counts keV",1000,0.,20000.); list_clov->Add(Ebeta_pmt);
  TH1F* Ebeta_sipm1              = new TH1F("Ebeta_sipm1","Ebeta_sipm1; E_{gamma} [keV]; Counts keV",4000,0.,80000.); list_clov->Add(Ebeta_sipm1);
  TH1F* Ebeta_sipm2              = new TH1F("Ebeta_sipm2","Ebeta_sipm2; E_{gamma} [keV]; Counts keV",4000,0.,80000.); list_clov->Add(Ebeta_sipm2);


  TH2F* gg_clov_beta = new TH2F("gg_clov_beta","Clover-Beta; E_{gamma} Clov [1.0 keV/Chan]; E_{gamma} Beta [1.0 keV/Chan]",2000,0.,2000.,8000,0.,8000.); list_clov->Add(gg_clov_beta); 

  TH2F* gg_clov_beta_pmt = new TH2F("gg_clov_beta_pmt","Clover-Beta_pmt; E_{gamma} Clov [1.0 keV/Chan]; E_{gamma} Beta [1.0 keV/Chan]",2000,0.,2000.,10000,0.,10000.); list_clov->Add(gg_clov_beta_pmt);
  TH2F* gg_clov_beta_sipm1 = new TH2F("gg_clov_beta_sipm1","Clover-Beta_sipm1; E_{gamma} Clov [1.0 keV/Chan]; E_{gamma} Beta [1.0 keV/Chan]",2000,0.,2000.,32000,0.,32000.); list_clov->Add(gg_clov_beta_sipm1);
  TH2F* gg_clov_beta_sipm2 = new TH2F("gg_clov_beta_sipm2","Clover-Beta_sipm2; E_{gamma} Clov [1.0 keV/Chan]; E_{gamma} Beta [1.0 keV/Chan]",2000,0.,2000.,32000,0.,32000.); list_clov->Add(gg_clov_beta_sipm2);


  //TAC histograms
  auto* list_tac = new TList; //list of the TAC histograms
  TH2F* tac_labr1_labr2     = new TH2F("tac_labr1_labr2","TAC START: LaBr1 - STOP: LaBr2; Time [ps]; E_{gamma} [keV/Chan]; Counts/ps",10000,0,100000.,1000,0.,4000.); list_tac->Add(tac_labr1_labr2);
  TH2F* tac_labr1_labr2_gate     = new TH2F("tac_labr1_labr2_gate","TAC START: LaBr1 - STOP: LaBr2; Time [ps]; E_{gamma} [keV/Chan]; Counts/ps",10000,0,100000.,1000,0.,4000.); list_tac->Add(tac_labr1_labr2_gate);
  TH2F* tac_betaPMT_labr1   = new TH2F("tac_betaPMT_labr1","TAC START: betaPMT - STOP: LaBr1; Time [ps]; E_{gamma} [keV/Chan]; Counts/10 ps",10000,0.,100000.,2000,0.,2000.); list_tac->Add(tac_betaPMT_labr1);
  TH2F* tac_betaPMT_labr2   = new TH2F("tac_betaPMT_labr2","TAC START: betaPMT - STOP: LaBr2; Time [ps]; E_{gamma} [keV/Chan]; Counts/10 ps",10000,0.,100000.,2000,0.,2000.); list_tac->Add(tac_betaPMT_labr2);
  TH2F* tac_betaSiPM1_labr1 = new TH2F("tac_betaSiPM1_labr1","TAC START: betaSiPM1 - STOP: LaBr1; Time [ps]; E_{gamma} [keV/Chan]; Counts/10 ps",10000,0.,100000.,2000,0.,2000.); list_tac->Add(tac_betaSiPM1_labr1);
  TH2F* tac_betaSiPM1_labr2 = new TH2F("tac_betaSiPM1_labr2","TAC START: betaSiPM1 - STOP: LaBr2; Time [ps]; E_{gamma} [keV/Chan]; Counts/10 ps",10000,0.,100000.,2000,0.,2000.); list_tac->Add(tac_betaSiPM1_labr2);
  TH2F* tac_betaSiPM2_labr1 = new TH2F("tac_betaSiPM2_labr1","TAC START: betaSiPM2 - STOP: LaBr1; Time [ps]; E_{gamma} [keV/Chan]; Counts/10 ps",10000,0.,100000.,2000,0.,2000.); list_tac->Add(tac_betaSiPM2_labr1);
  TH2F* tac_betaSiPM2_labr2 = new TH2F("tac_betaSiPM2_labr2","TAC START: betaSiPM2 - STOP: LaBr2; Time [ps]; E_{gamma} [keV/Chan]; Counts/10 ps",10000,0.,100000.,2000,0.,2000.); list_tac->Add(tac_betaSiPM2_labr2);

  //Beta histograms
  auto* list_beta = new TList; //list of the Beta histograms
  TH1F* betaDet[4];
  for(Int_t i=0; i<4; i++){
    betaDet[i]  = new TH1F(Form("betaDet_%d",i),Form("Beta %s; E_{beta} [Chan]; Counts/10 Chan",BetaName[i].c_str()),100000,0.,100000.); list_beta->Add(betaDet[i]);
  }

  //Others
  auto* list_others = new TList; //list of the extra histograms
  TH1F* TRef = new TH1F("TRef","Time difference; Time 10 ms/Chan",400000,0,400000); list_others->Add(TRef);

  // ==============================
  // CRYSTAL PAIRING HISTOGRAMS 
  // ==============================

  // Function to create the crystal correlations

  auto createHistograms = [&](TList* list, const std::string& prefix, const std::string& title_prefix) {
      for (int i = 0; i < angular_binning; i++) {
        double bin_center = get_bin_center(i * angular_resolution);
        std::stringstream ss;
        ss << std::fixed << std::setprecision(2) << bin_center;
        std::string name = prefix + std::to_string(i);
        std::string title = title_prefix + ss.str() + "; Energy [keV]; Energy [keV]";
        TH2F* crystal_correlations = new TH2F(name.c_str(), title.c_str(), 4000, 0, 4000, 4000, 0, 4000);
        list->Add(crystal_correlations);
      }
  };

  // LISTS FOR ANGULAR CORRELATIONS FOR DIFFERENT CONDITIONS
  auto* list_crys = new TList;
  auto* list_crys_bckg = new TList;
  auto* list_crys_true = new TList;
  auto* list_crys_mixed = new TList;
  auto* list_crys_beta_gated = new TList;
  auto* list_crys_bckg_beta_gated = new TList;
  auto* list_crys_true_beta_gated = new TList;
  auto* list_crys_mixed_beta_gated = new TList;

  createHistograms(list_crys, "Angular_correlation_hist_", "Histogram for angle ");
  createHistograms(list_crys_bckg, "Angular_correlation_hist_bckg_", "Histogram background for ");
  createHistograms(list_crys_mixed, "Angular_correlation_hist_mixed_", "Histogram mixed for angle ");
  createHistograms(list_crys_beta_gated, "Angular_correlation_hist_beta_gated_", "Histogram beta gated for angle ");
  createHistograms(list_crys_bckg_beta_gated, "Angular_correlation_hist_bckg_beta_gated_", "Histogram background beta gated for ");
  createHistograms(list_crys_mixed_beta_gated, "Angular_correlation_hist_mixed_beta_gated_", "Histogram mixed beta gated for angle ");

  //Here is where we start looping over the events
  Long64_t nentries = tree->GetEntries();
  cout<< " Filling Clovers, Betas, TACs and Vetoes" << endl;
  cout<< " Total events: "<< nentries/1e6 << " millions"<<endl;

  //Long64_t totalEvents=0;
  float progress = 0.0, progress_buffer = 0.0;

  int beta_count=0;
  int beta_indv_count[num_beta];
  //ULong64_t TimeRef_value = 0;

  // ==================
  // FILLING HISTOGRAMS 
  // ==================

  for (Long64_t i=0; i<nentries; i++) {//Loops over all the events (entries) of the tree
    tree->GetEntry(i);

    beta_count = 0;
    for(Int_t j = 0; j < num_beta; j++)
      beta_indv_count[j] = 0;

    if(gotBeta){
      for(Int_t j = 0; j < num_beta; j++){//Loop over beta event
        if(E_Beta[j] > 0.) {
          betaDet[j]->Fill(E_Beta[j]);
          beta_count+=1;
        }
        if(j == beta_pmt && E_Beta[j] > betapmt_energy_gate[0] && E_Beta[j] < betapmt_energy_gate[1]){
          beta_indv_count[beta_pmt]+=1;
        }
        else if((j == beta_sipm1 || j == beta_sipm2) && E_Beta[j] > betasipm_energy_gate[0] && E_Beta[j] < betasipm_energy_gate[1]){
          beta_indv_count[j]+=1;
        }
      }
    }//Beta boolean check

    if(gotClov){//Check for the very unlikely event we run without clovers
      for(Int_t j=0; j<num_clov_crys; j++){//Loop over first clover event
        if(E_Clov[j]>1. ){   // && TIME_REF < time_ref_c
          clov_singles->Fill(E_Clov[j]);
          //if (TIME_REF < time_ref_c){                ///////////////////////////////////////////////////////////////////////////////////
            clov_vs_TRef->Fill(TIME_REF,E_Clov[j]);
            TRef->Fill(TIME_REF);
          //}
          if(beta_count>0 ){ /////////////////////////////////////////////////////////////////////////////////// && TIME_REF < time_ref_c
            clov_singles_betagated->Fill(E_Clov[j]);//If some detector found a beta, we write the clover energy
            clov_Id_vs_Ener->Fill(j,E_Clov[j]);
            //if (TIME_REF < time_ref_c){
              clov_vs_TRef_betagated->Fill(TIME_REF,E_Clov[j]);

             for(Int_t k=0; k<num_beta; k++){//Loop over first labr event
              if(E_Beta[k]>0.){    ////////////////////////////////////////////////////////////////////////// k=0/                  1
                time_diff = T_Clov[j]-T_Beta[k];
                clov_beta_time->Fill(time_diff);
                if(clov_beta_t_coinc[0] > time_diff && time_diff < clov_beta_t_coinc[1]){  //////// gg clov condition
                  gg_clov_beta->Fill(E_Clov[j], E_Beta[k]);

                  if(E_Beta[0]>0){
                    gg_clov_beta_pmt->Fill(E_Clov[j], E_Beta[0]);
                  }
                  if(E_Beta[1]>1){
                    gg_clov_beta_sipm1->Fill(E_Clov[j], E_Beta[1]);
                  }
                  if(E_Beta[2]>2){
                    gg_clov_beta_sipm2->Fill(E_Clov[j], E_Beta[2]);
                  }

                  clov_timebeta->Fill(time_diff, E_Clov[j]);
                  if(E_Beta[k]>0){
                    beta_timeclov->Fill(time_diff, E_Beta[k]);
                  }

                  if(E_Beta[0]>0){
                    Ebeta_pmt->Fill(E_Beta[0]);
                    if(E_Clov[j]>127. && E_Clov[j]<131.){
                      beta_timeclov_pmt->Fill(time_diff, E_Beta[0]);
                    }
                  }

                  if(E_Beta[1]>0){
                    if(E_Clov[j]>127. && E_Clov[j]<131.){
                      beta_timeclov_sipm1->Fill(time_diff, E_Beta[1]);                                   //// if(E_Clov[j]>1200. && E_Clov[j]<1204.){   //// 100Rb (1201.7 keV) --> ~85 ns
                    }
                    Ebeta_sipm1->Fill(E_Beta[1]);
                  }

                  if(E_Beta[2]>0){
                    if(E_Clov[j]>127. && E_Clov[j]<131.){  /// 100Rb (1201.7 keV)                    //// if(E_Clov[j]>894. && E_Clov[j]<901.){   //// 88Rb (898.03 keV)
                      beta_timeclov_sipm2->Fill(time_diff, E_Beta[2]);
                    }
                    Ebeta_sipm2->Fill(E_Beta[2]);
                  }
                }
              }
            }
          }

          /*
          if(gotLaBr){//Check to see if we run with LaBr
            for(Int_t k=0; k<num_labr; k++){//Loop over first labr event
              if(E_LaBr[k]>1.){
                time_diff = T_Clov[j]-T_LaBr[k];
                clov_labr_time->Fill(time_diff);
                if(clov_labr_t_coinc[0] > time_diff && time_diff < clov_labr_t_coinc[1]){
                  gg_clov_labr->Fill(E_Clov[j], E_LaBr[k]);
                  if(beta_indv_count[beta_pmt]>0) {   //   if(beta_count>0)    Only beta_pmt contribution
                    gg_clov_labr_betagated->Fill(E_LaBr[k], E_Clov[j]);
                    //if(E_Clov[j]>1922. && E_Clov[j]<1928.){                                           //// 100Sr (129.2 keV), LT ~ 4 ns   Bgg (1926.8 keV-129.2 keV)
                    //if(E_Clov[j]>1090. && E_Clov[j]<1094.){                                           //// 101Sr (271.2 keV), LT ~ 0.1 ns   Bgg (1091.8 keV-271.2 keV)
                    //if(E_Clov[j]>269. && E_Clov[j]<273.){                                           //// 101Sr (271.2 keV), LT ~ 0.1 ns   Bgg (271.2 keV-1091.8 keV)
                    //if(E_Clov[j]>250. && E_Clov[j]<254.){                                           //// 101Sr (111.6 keV), LT ~ 0.1 ns   Bgg (251.6 keV-111.6 keV)
                    //if(E_Clov[j]>230. && E_Clov[j]<234.){                                           //// 101Sr (232.7 keV), LT ~ 0.4 ns   Bgg (232.7 keV-251.6 keV)
                    //if(E_Clov[j]>230. && E_Clov[j]<235.){                                           //// 101Sr (363.2 keV), Bgg (232.7 keV-251.6/363.2 keV)

                    //if(E_Clov[j]>995. && E_Clov[j]<999.){                                           //// 101Y (996.53 keV) LT ~ 0.1 ns   Bgg (996.53 keV-128.34 keV)
                    //if(E_Clov[j]>127. && E_Clov[j]<131.){                                           //// 101Y (128.34), Bgg (128.34 keV-996.53 keV)
                    //if(E_Clov[j]>162. && E_Clov[j]<165.){                                           //// 101Y (163.35), Bgg (163.35 keV-128.34 keV)
                    //if(E_Clov[j]>421. && E_Clov[j]<425.){                                           //// 101Y (422.84), Bgg (422.84 keV-163.35 keV)

                    //if(E_Clov[j]>132. && E_Clov[j]<136.){                                           //// 101Zr (133.67), Bgg (133.67 keV-98.21 keV)
                    //if(E_Clov[j]>103. && E_Clov[j]<106.){                                           //// 101Zr (104.43), Bgg (104.43 keV-216.68 keV)
                    //if(E_Clov[j]>174. && E_Clov[j]<178.){                                           //// 101Zr (176.44), Bgg (176.44 keV-133.65/231.91 keV)
                    //if(E_Clov[j]>230. && E_Clov[j]<233.){                                           //// 101Zr (231.91), Bgg (231.91 keV-176.44 keV)
                    if(E_Clov[j]>144. && E_Clov[j]<148.){                                           //// 101Zr (146.64), Bgg (146.64 keV-104.43 keV)

                    beta_LaBr_clov->Fill(time_diff, E_LaBr[k]);                                       //// Original: beta_LaBr_clov->Fill(time_diff, E_Clov[j]);   /// E_Clov[j] -->E_LaBr[k]
                    }
                    //else if(E_Clov[j]>276. && E_Clov[j]<280.){          // same width  !!!!
                    //else if(E_Clov[j]>258. && E_Clov[j]<262.){          // same width  !!!!
                    //else if(E_Clov[j]>237. && E_Clov[j]<241.){          // same width  !!!!    101Sr (232.7 keV)
                    //else if(E_Clov[j]>1101. && E_Clov[j]<1105.){        // same width  !!!!

                    //else if(E_Clov[j]>1020. && E_Clov[j]<1024.){        // same width  !!!!   101Y (996.53 keV) LT ~ 0.1 ns   Bgg (996.53 keV-128.34 keV)
                    //else if(E_Clov[j]>139. && E_Clov[j]<143.){        // same width  !!!!     101Y (128.34) LT ~ 0.1 ns   Bgg (128.34 keV-996.53 keV)
                    //else if(E_Clov[j]>169. && E_Clov[j]<172.){        // same width  !!!!       101Y (163.35), Bgg (163.35 keV-128.34 keV)
                    //else if(E_Clov[j]>426. && E_Clov[j]<430.){        // same width  !!!!       101Y (422.84), Bgg (422.84 keV-163.35 keV)

                    //else if(E_Clov[j]>138. && E_Clov[j]<142.){        // same width  !!!!       101Zr (133.67), Bgg (133.67 keV-98.21 keV)
                    //else if(E_Clov[j]>106. && E_Clov[j]<109.){        // same width  !!!!       101Zr (104.43), Bgg (104.43 keV-216.68 keV)
                    //else if(E_Clov[j]>190. && E_Clov[j]<194.){        // same width  !!!!       101Zr (176.44), Bgg (176.44 keV-133.65/231.91 keV)
                    //else if(E_Clov[j]>235. && E_Clov[j]<238.){          // same width  !!!!         101Zr (231.91), Bgg (231.91 keV-176.44 keV)
                    else if(E_Clov[j]>149. && E_Clov[j]<153.){          // same width  !!!!         101Zr (146.64), Bgg (146.64 keV-104.43 keV)

                    beta_LaBr_clov_bkg->Fill(time_diff, E_LaBr[k]);       // background histogram
                    beta_LaBr_clov->Fill(time_diff, E_LaBr[k], -1.);      // subtract the background contribution
                    }
                    }
                }
              }
            }
          }
          */

          for(Int_t k=0; k<num_clov_crys; k++){//Loop over second and beyond clover envents
            if (TIME_REF < time_ref_c){
              if(E_Clov[k] > 1. && j!=k ){
                time_diff = T_Clov[j]-T_Clov[k];
                clov_time->Fill(time_diff);
                if(abs(time_diff) < clov_t_coinc){//True time coincidences
                  clov_gg->Fill(E_Clov[j], E_Clov[k]);
                  std::string key = "('" + clover[j] + "', '" + color[j] + "') - ('" + clover[k] + "', '" + color[k] + "')";
                  try {
                    double angulo = angulos.at(key);
                    double angulo_binned = get_bin_center(angulo);
                    int bin_index = static_cast<int>((angulo_binned / angular_resolution)-0.5);
                    TH2D* hist_event_correlation = (TH2D*)list_crys->At(bin_index);
                    hist_event_correlation->Fill(E_Clov[j], E_Clov[k]);
                    } catch (const std::out_of_range& e) {}

                  if(beta_indv_count[0]>0){ // For gamma-gamma events only pmt beta detections
                    clov_gg_betagated->Fill(E_Clov[j], E_Clov[k]);
                  
                    // Angular correlation Beta gated
                    std::string key_beta_gated = "('" + clover[j] + "', '" + color[j] + "') - ('" + clover[k] + "', '" + color[k] + "')";
                    try {
                      double angulo = angulos.at(key_beta_gated);
                      double angulo_binned = get_bin_center(angulo);
                      int bin_index = static_cast<int>((angulo_binned / angular_resolution)-0.5);
                      TH2D* hist_event_correlation_beta_gated = (TH2D*)list_crys_beta_gated->At(bin_index);
                      hist_event_correlation_beta_gated->Fill(E_Clov[j], E_Clov[k]);
                    } catch (const std::out_of_range& e) {}
                  }
                    if(E_Clov[j]>285. && E_Clov[j]<291.){ // gg filter
                      histo->Fill(E_Clov[k], TIME_REF); // similar to clov_vs_TRef_betagated
                      histo2->Fill(E_Clov[k], TIME_REF);
                    }
                    else if(E_Clov[j]>302. && E_Clov[j]<308.){ // same width  !!!!
                      histo_bkg->Fill(E_Clov[k], TIME_REF); // background histogram
                      histo->Fill(E_Clov[k], TIME_REF, -1.); // subtract the background contribution
                    }
                }
              else if(clov_rand_coinc[0] < abs(time_diff) && abs(time_diff) < clov_rand_coinc[1]){//Time-random coincidences
                clov_gg_trand->Fill(E_Clov[j], E_Clov[k]);

                // Fill the background angular correlation histograms between a defined window
                std::string key_bckg = "('" + clover[j] + "', '" + color[j] + "') - ('" + clover[k] + "', '" + color[k] + "')";
                try {
                  double angulo = angulos.at(key_bckg);
                  double angulo_binned = get_bin_center(angulo);
                  int bin_index = static_cast<int>((angulo_binned / angular_resolution)-0.5);
                  TH2D* hist_event_correlation_bckg = (TH2D*)list_crys_bckg->At(bin_index);
                  hist_event_correlation_bckg->Fill(E_Clov[j], E_Clov[k]);
                  } catch (const std::out_of_range& e) {}
                if(beta_count>0){
                  clov_gg_trand_betagated->Fill(E_Clov[j], E_Clov[k]);
                   // Angular correlation background beta gated
                  std::string key_bckg_beta_gated = "('" + clover[j] + "', '" + color[j] + "') - ('" + clover[k] + "', '" + color[k] + "')";
                  try {
                    double angulo = angulos.at(key_bckg_beta_gated);
                    double angulo_binned = get_bin_center(angulo);
                    int bin_index = static_cast<int>((angulo_binned / angular_resolution)-0.5);
                    TH2D* hist_event_correlation_bckg_beta_gated = (TH2D*)list_crys_bckg_beta_gated->At(bin_index);
                    hist_event_correlation_bckg_beta_gated->Fill(E_Clov[j], E_Clov[k]);
                  } catch (const std::out_of_range& e) {}
                }            
              }
              } 
            }//TIME_REF
            if (abs(k - j) == event_separation) {
              // ================================
              // Event mixing (no beta gating)
              // ================================
              std::string key_mixed = "('" + clover[j] + "', '" + color[j] + "') - ('" + clover[k] + "', '" + color[k] + "')";
              try {
                double angulo = angulos.at(key_mixed);
                double angulo_binned = get_bin_center(angulo);
                int bin_index = static_cast<int>((angulo_binned / angular_resolution) - 0.5);
                TH2D* hist_event_correlation_mixed = (TH2D*)list_crys_mixed->At(bin_index);
                hist_event_correlation_mixed->Fill(E_Clov[j], E_Clov[k]);
                } catch (const std::out_of_range& e) {}

              // ===============================
              // Event mixing with beta gating PMT
              // ===============================
              if (gotBeta) {
                // event j: first beta gated 
                bool beta_j = (beta_indv_count[beta_pmt] > 0);

                // Evento k: second beta gated boolean logic
                bool beta_k = false;
                if (E_Beta[beta_pmt] > betapmt_energy_gate[0] &&
                  E_Beta[beta_pmt] < betapmt_energy_gate[1]) {
                  beta_k = true;
                }

                // If both beta gated, fill mixed beta gated histograms
                if (beta_j && beta_k) {
                  std::string key_mixed_beta = "('" + clover[j] + "', '" + color[j] + "') - ('" + clover[k] + "', '" + color[k] + "')";
                  try {
                    double angulo = angulos.at(key_mixed_beta);
                    double angulo_binned = get_bin_center(angulo);
                    int bin_index = static_cast<int>((angulo_binned / angular_resolution) - 0.5);
                    TH2D* hist_event_correlation_mixed_beta_gated =
                    (TH2D*)list_crys_mixed_beta_gated->At(bin_index);
                    hist_event_correlation_mixed_beta_gated->Fill(E_Clov[j], E_Clov[k]);
                  } catch (const std::out_of_range& e) {}
                }       
              }
            }
          }
        }
      }
    }//Clover boolean check
/*
    if(gotLaBr){//Check to see if we run with LaBr
    if (TIME_REF < time_ref_c){
      for(Int_t j=0; j<num_labr; j++){//Loop over first labr event
        if(E_LaBr[j]>1.){
          labr_singles->Fill(E_LaBr[j]);
          if(beta_count>0)
            labr_singles_betagated->Fill(E_LaBr[j]);
          for(Int_t k=0; k<num_labr; k++){//Loop over second and beyond clover envents
            if(E_LaBr[k]>1. && j!=k ){
              time_diff=T_LaBr[j]-T_LaBr[k];
              labr_time->Fill(time_diff);
              if(abs(time_diff)<labr_t_coinc){//True time coincidences
                labr_gg->Fill(E_LaBr[j], E_LaBr[k]);
                if(beta_count>0)
                  labr_gg_betagated->Fill(E_LaBr[j], E_LaBr[k]);
              }
            }
          }
        }
      }
      }//TIME_REF
    }//LaBr boolean check
*/
/*
    if(gotTAC){//Check to see if we run with TAC
      if (TIME_REF < time_ref_c){
      for(Int_t j=0; j<num_tac; j++){//Loop over tac events
        if(E_TAC[j]>1.){
          if(j==tac_labr){
          tac_labr1_labr2->Fill(E_TAC[j], E_LaBr[1]);
          }

          if(j==tac_labr){//You need to edit the TAC indexes at the begining of the script
            //100Rb_run409_run410
            //if(E_LaBr[0]>255. && E_LaBr[0]<285. && E_LaBr[1]>0){        // 101Sr 271.0->134.8... keV
            //if(E_LaBr[0]>90. && E_LaBr[0]<120. && E_LaBr[1]>0){         // 101Sr 111.6->251.6... keV
            //if(E_LaBr[0]>240. && E_LaBr[0]<270. && E_LaBr[1]>0){        // 101Sr 251.6->111.6... keV
            //if(E_LaBr[0]>125. && E_LaBr[0]<145. && E_LaBr[1]>0){        // 101Sr 134.8->271.0... keV
            //if(E_LaBr[0]>1190. && E_LaBr[0]<1225. && E_LaBr[1]>0){      // 100Sr 1201.7->287... keV

            //101Rb_Results

            //if(E_LaBr[0]>1050. && E_LaBr[0]<1200. && E_LaBr[1]>0){        // 101Sr 1091.8->271.2... keV
            //if(E_LaBr[0]>255. && E_LaBr[0]<285. && E_LaBr[1]>0){        // 101Sr 271.2->1091.8... keV
            //if(E_LaBr[0]>240. && E_LaBr[0]<270. && E_LaBr[1]>0){        // 101Sr 251.6->111.6... keV
            //if(E_LaBr[0]>100. && E_LaBr[0]<120. && E_LaBr[1]>0){        // 101Sr 111.6->251.6... keV

            //if(E_LaBr[0]>155. && E_LaBr[0]<180. && E_LaBr[1]>0){        // 101Y 163.35->128.34... keV
            //if(E_LaBr[0]>440. && E_LaBr[0]<475. && E_LaBr[1]>0){        // 101Y 462.14->128.34... keV
            //if(E_LaBr[0]>980. && E_LaBr[0]<1030. && E_LaBr[1]>0){        // 101Y 996.53->128.34... keV
            //if(E_LaBr[0]>120. && E_LaBr[0]<150. && E_LaBr[1]>0){        // 101Y 128.34->163.35/462.14/996.53... keV (Centroid shift)

            //if(E_LaBr[0]>120. && E_LaBr[0]<145. && E_LaBr[1]>0){        // 101Zr 133.67->98.21/176.44... keV
            //if(E_LaBr[0]>90. && E_LaBr[0]<115. && E_LaBr[1]>0){        // 101Zr 98.21->133.67... keV
            //if(E_LaBr[0]>200. && E_LaBr[0]<225. && E_LaBr[1]>0){        // 101Zr 216.68->104.43... keV
            if(E_LaBr[0]>160. && E_LaBr[0]<180. && E_LaBr[1]>0){        // 101Zr 176.44->133.67... keV

            //152Eu_run444   Time-ref
            //if(E_LaBr[0]>105. && E_LaBr[0]<135. && E_LaBr[1]>0){   // 152Sm 121.7817->244, 964, 1112, 1408... keV
            //if(E_LaBr[0]>225. && E_LaBr[0]<265. && E_LaBr[1]>0){   // 152Sm 244.6974->121... keV
            //if(E_LaBr[0]>835. && E_LaBr[0]<900. && E_LaBr[1]>0){   // 152Sm 867.38->244... keV
            //if(E_LaBr[0]>325. && E_LaBr[0]<365. && E_LaBr[1]>0){   // 152Gd 344.2789->779... keV
            //if(E_LaBr[0]>740. && E_LaBr[0]<810. && E_LaBr[1]>0){   // 152Gd 778.9045->344... keV
            //if(E_LaBr[0]>940. && E_LaBr[0]<980. && E_LaBr[1]>0){   // 152Sm 964.057->121... keV
            //if(E_LaBr[0]>1090. && E_LaBr[0]<1130. && E_LaBr[1]>0){   // 152Sm 1112.076->121... keV
            //if(E_LaBr[0]>1385. && E_LaBr[0]<1425. && E_LaBr[1]>0){   // 152Sm 1408.013->121... keV
            //if(E_LaBr[0]>540. && E_LaBr[0]<585. && E_LaBr[1]>0){   // 152Sm 563.986->1085... keV
            //if(E_LaBr[0]>1065. && E_LaBr[0]<1105. && E_LaBr[1]>0){   ////////////////////////////////////////////////////////////////////////152Sm 1085.837->563... keV
            //if(E_LaBr[0]>420 && E_LaBr[0]<465. && E_LaBr[1]>0){   // 152Sm 443.9606->1086... keV
            //if(E_LaBr[0]>1325. && E_LaBr[0]<1435. && E_LaBr[1]>0){   // 152Gd 1299.142->344... keV
            //if(E_LaBr[0]>950. && E_LaBr[0]<1050. && E_LaBr[1]>0){   // 152Gd 964.057->444... keV
            //if(E_LaBr[0]>1090. && E_LaBr[0]<1180. && E_LaBr[1]>0){   // 152Gd 1085.84->444... keV
            //if(E_LaBr[0]>1450. && E_LaBr[0]<1550. && E_LaBr[1]>0){   // 152Gd 563.986->1085... keV




            // Time-walk 2.0
            //if(E_LaBr[0]>425. && E_LaBr[0]<450. && E_LaBr[1]>0){   // 152Sm 443.9606->964.057, 1085.84... keV
            //if(E_LaBr[0]>955. && E_LaBr[0]<1015. && E_LaBr[1]>0){   // 152Sm 964.057->443.9606... keV
            //if(E_LaBr[0]>1110. && E_LaBr[0]<1180. && E_LaBr[1]>0){   // 152Sm 1085.84->443.9606... keV
            //if(E_LaBr[0]>230. && E_LaBr[0]<260. && E_LaBr[1]>0){   // 152Sm 244.6974->867.38... keV
            //if(E_LaBr[0]>840. && E_LaBr[0]<920. && E_LaBr[1]>0){   // 152Sm 867.38->244.6974... keV
            //if(E_LaBr[0]>315. && E_LaBr[0]<370. && E_LaBr[1]>0){   // 152Sm 344.2789->778.9045, 1299.142... keV
            //if(E_LaBr[0]>740. && E_LaBr[0]<810. && E_LaBr[1]>0){   // 152Sm 778.9045->344.2789... keV
            //if(E_LaBr[0]>1340. && E_LaBr[0]<1400. && E_LaBr[1]>0){   // 152Sm 1299.344->344.2789... keV

            tac_labr1_labr2_gate->Fill(E_TAC[j], E_LaBr[1]); /////////////////////////////////////////////////////////////////////////////////////
            }
          }
          else if(j==tac_betapmt){//You need to edit the TAC indexes at the begining of the script
            if(E_LaBr[0]>1 && E_LaBr[1] ==0){
              tac_betaPMT_labr1->Fill(E_TAC[j],E_LaBr[0]);
            }else if(E_LaBr[1]>1 && E_LaBr[0] ==0){
              tac_betaPMT_labr2->Fill(E_TAC[j],E_LaBr[1]);
            }
          }
          else if(j==tac_betasipm1){//You need to edit the TAC indexes at the begining of the script
            if(E_LaBr[0]>1 && E_LaBr[1] ==0){
              tac_betaSiPM1_labr1->Fill(E_TAC[j],E_LaBr[0]);
            }else if(E_LaBr[1]>1 && E_LaBr[0] ==0){
              tac_betaSiPM1_labr2->Fill(E_TAC[j],E_LaBr[1]);
            }
          }
          else if(j==tac_betasipm2){//You need to edit the TAC indexes at the begining of the script
            if(E_LaBr[0]>1 && E_LaBr[1] ==0){
              tac_betaSiPM2_labr1->Fill(E_TAC[j],E_LaBr[0]);
            }else if(E_LaBr[1]>1 && E_LaBr[0] ==0){
              tac_betaSiPM2_labr2->Fill(E_TAC[j],E_LaBr[1]);
            }
          }
        }
      }
      }//TIME_REF
    }//TAC boolean check
*/
    progress = ((float)i/(float)nentries)*100.0;
    if (progress >= progress_buffer || i == nentries){
      PrintBar(progress);
      cout.flush();
      progress_buffer = progress + 1.0;
    }

  }// End LOOP OVER TTREE EVENTS

  progress = progress + 1.0;
  PrintBar(progress);
  cout.flush();
  cout << endl << endl << " Sorting finished!" << endl << " Saving the histograms...";

  //We subtract the time random coincidences from the true ones:
  clov_gg_Notrand->Add(clov_gg,1);
  clov_gg_Notrand_betagated->Add(clov_gg_betagated,1);
  clov_gg->Add(clov_gg_trand,-1.*clov_t_coinc/(clov_rand_coinc[1] - clov_rand_coinc[0]));
  clov_gg_betagated->Add(clov_gg_trand_betagated,-1.*clov_t_coinc/(clov_rand_coinc[1] - clov_rand_coinc[0]));

// Angular Correlation True Coincidences
  
for (int i = 0; i < angular_binning; i++) {
  TH2D* true_corr = (TH2D*)(((TH2D*)list_crys->At(i))->Clone());
  true_corr->SetName(Form("Angular_correlation_hist_true_%d", i));
  true_corr->Add((TH2F*)list_crys_bckg->At(i), -1.*clov_t_coinc/(clov_rand_coinc[1] - clov_rand_coinc[0]));
  list_crys_true->Add(true_corr);
}

// Angular Correlation True Coincidences Beta Gated
  
for (int i = 0; i < angular_binning; i++) {
    TH2D* true_beta_gated = (TH2D*)(((TH2D*)list_crys_beta_gated->At(i))->Clone());
    true_beta_gated->SetName(Form("Angular_correlation_hist_true_beta_gated_%d", i));
    true_beta_gated->Add((TH2D*)list_crys_bckg_beta_gated->At(i), -1.*clov_t_coinc/(clov_rand_coinc[1] - clov_rand_coinc[0]));
    list_crys_true_beta_gated->Add(true_beta_gated);
}

// ==================
// WRITING HISTOGRAMS
// ==================

  //We create the output file and write the histograms
  TFile *fout = new TFile(outputFileName.c_str(),"RECREATE");

  //list_others->Write();
  if(gotClov) list_clov->Write(); //We only write the histos if we had clovers
  if(gotClov) list_crys->Write(); // Angular correlations
  if(gotClov) list_crys_bckg->Write(); // Angular correlations background
  if(gotClov) list_crys_true->Write(); // Angular correlations true = list_crys - list_crys_bckg
  if(gotClov) list_crys_mixed->Write(); // Angular Correlations mixed 
  if(gotBeta) list_beta->Write(); //We only write the histos if we had Betas
  if(gotBeta && gotClov) list_crys_beta_gated->Write(); // Angular correlations beta gated
  if(gotBeta && gotClov) list_crys_bckg_beta_gated->Write(); // Angular correlations background beta gated
  if(gotBeta && gotClov) list_crys_true_beta_gated->Write(); // Angular correlations true beta gated
  if(gotBeta && gotClov) list_crys_mixed_beta_gated->Write(); //Angular correlations mixed beta gated
  //if(gotLaBr) list_labr->Write(); //We only write the histos if we had LaBr
  //if(gotClov && gotLaBr) list_clovlabr->Write();
  //if(gotTAC)  list_tac->Write();  //We only write the histos if we had TAC

  fout->Close();

  cout << " done!" << endl << endl;

  return 0;
}