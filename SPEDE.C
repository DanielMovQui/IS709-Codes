//****************************************************************************************************************//
//                                                                                                                //
//                          Online/Offline Code for IDS/Fast-timing experiments                                   //
//                                                                                                                //
// Written by:    A. Illana (UCM)                                                                                 //
// Date:          2025-09-09                                                                                      //
// Experiment:    IS709                                                                                           //
// Root version:  6.30/06                                                                                         //
// Info:          This script has a TChain, so you can give a list of files and it will loop through all of them  //
//                                                                                                                //
// Latest Version: A. Illana (andres.illana@cern.ch) -> IS705 in May 2025                                         //
//                                                                                                                //
// To compile: g++ ids_2025SPEDEConf.C `root-config --cflags` `root-config --glibs` -o ids_2025SPEDEConf          //
//                                                                                                                //
//****************************************************************************************************************//
//
// ToDo:
//  - Parallelization...?


//C,C++ Libraries
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <time.h>
#include <vector>

// ROOT libraries
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TColor.h>
#include <TROOT.h>
#include <TChain.h>

// Global variables
TChain* tree;
std::string inputFileName, outputFileName;
static int TimeCoincAddback = 120;

// Running flags
const bool DoAddback  = true;
const bool Do60CoCase = false;

////// Starting the configuration for the experiment //////
// Number of each type of detectors
const int num_Clov_crys  = 52;  //This is number of crystals, not clovers, so usually crystals = 4*clovers
const int num_spede      = 24;   //This is number of SPEDE segments
const int num_beta       = 5;   //This is number of beta plastics

//Config Index to Names
//SiPM(Number)(Front/Back)(Top/Bottom)
//const std::string BetaName[num_beta] = {"SiPM1FT", "SiPM1FB", "SiPM1BT", "SiPM1BB", "SiPM2FT", "SiPM2FB", "SiPM2BT", "SiPM2BB", "SiPM3FT", "SiPM3FB", "SiPM3BT", "SiPM3BB"};

//For a quick analysis I will only consider one layer and one cable.
//const bool Beta_Good[num_beta] = {true, false, false, false, true, false, false, false, true, false, false, false };

// Reference ID crystal for timing evaluation
const int RefIDCrys = 0;

// Time differences, you probably want to edit
// Clover-clover time windows
const double Clov_t_coinc          = 100.;
const double Clov_rand_coinc[2]    = {  400., 1000.};
// Clover-Beta Time window
const double Clov_Beta_t_coinc[2]  = { -200.,   275.};
// Clover-SPEDE Time window
const double Clov_SPEDE_t_coinc[2] = { -200.,   60.};
const double SPEDE_SPEDE_t_coinc[2] = { -200.,   60.};

// Clover Energy min
const double EClov_Min = 20.;

// Beta Energy windows
const double BetaPMT_energy_gate[2]  = {10., 60000.};

// Beta Energy windows
//const double SPEDE_energy_gate      = 200.; //In channels
//const double SPEDE_energy_range[2]  = { 1000., 8000.};
const double SPEDE_energy_gate       = 30.;         //In keV
const double SPEDE_energy_range[2]   = {  4000., 4000.};


// Info time range used for the Xia4IDS
const double      RunUnit_val = 1.; 
const std::string RunUnit_str = "s"; 
const double      RefUnit_val = 1.; 
const std::string RefUnit_str = "ms";
const int         MaxRefUnit  = 4000; // MaxRefUnit*RefUnit_val in RefUnit_str

////// Here finish the configuration for the experiment //////


void PrintBar(float _progress){

  int barWidth = 70;
  std::cout << "[";
  int pos = barWidth * (_progress/100.);
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) std::cout << "=";
    else if (i == pos) std::cout << ">";
    else std::cout << " ";
  }
  if (_progress > 100.) _progress = 100.;
  std::cout << "] " << int(_progress) << " %\r";
  std::cout.flush();

  return;
}

// Addback algorithm
void addback(double Energy[num_Clov_crys], double Time[num_Clov_crys], int mult_clover){
  
  int mult_clover_new = 0;
  std::vector<int>    tmp_ClovSg;
  std::vector<double> tmp_ClovEn;
  std::vector<double> tmp_ClovEn_AB;
  std::vector<double> tmp_ClovTi;
  
  //Clovers (6clovers * 4crystals)
  for (int i = 0; i < num_Clov_crys/4; i++){   // for each clover
    if (tmp_ClovSg.size() > 0){
      tmp_ClovSg.clear();
      tmp_ClovEn.clear();
      tmp_ClovEn_AB.clear();
      tmp_ClovTi.clear();
    }
    for(int j = 0; j < 4; j++){ // for each crystal
      if (Energy[(4*i)+j] > 0){
        tmp_ClovSg.push_back(j);
        tmp_ClovEn.push_back(Energy[(4*i)+j]);
        tmp_ClovEn_AB.push_back(Energy[(4*i)+j]);
        tmp_ClovTi.push_back(Time[(4*i)+j]);
      }
    }
    
    //New Adddback algorithm
    for(int j = 0; j < tmp_ClovSg.size(); j++){
      for (int k = j+1; k < tmp_ClovSg.size(); k++){
        if ((tmp_ClovEn.at(j) >= tmp_ClovEn.at(k)) && (TMath::Abs(tmp_ClovTi.at(j) - tmp_ClovTi.at(k)) <= TimeCoincAddback) && tmp_ClovEn.at(j) > 0 && tmp_ClovEn.at(k) > 0 ){
          tmp_ClovEn_AB.at(j) = tmp_ClovEn_AB.at(j) + tmp_ClovEn_AB.at(k);
          tmp_ClovEn_AB.at(k) = - 1.;
          tmp_ClovTi.at(k) = 0;
        }
        else if ((tmp_ClovEn.at(j) < tmp_ClovEn.at(k)) && (TMath::Abs(tmp_ClovTi.at(j) - tmp_ClovTi.at(k)) <= TimeCoincAddback) && tmp_ClovEn.at(j) > 0 && tmp_ClovEn.at(k) > 0 ){
          tmp_ClovEn_AB.at(k) = tmp_ClovEn_AB.at(j) + tmp_ClovEn_AB.at(k);
          tmp_ClovEn_AB.at(j) = - 1.;
          tmp_ClovTi.at(j) = 0;
        }
      }
    }
    // End of addback
    for(int j = 0; j < tmp_ClovSg.size(); j++){
      int tempID = (4*i)+tmp_ClovSg.at(j);
      if (tmp_ClovEn_AB.at(j) > 0){
        Energy[tempID] = tmp_ClovEn_AB.at(j);
        Time  [tempID] = tmp_ClovTi.at(j);
        mult_clover_new++;
      }
      else{
        Energy[tempID] = 0.;
        Time  [tempID] = 0;
      }
    }
  }
  mult_clover = mult_clover_new;
  
  return;
}

inline bool CheckIfFileExist(std::string FileName, bool CheckInputs){
  
  std::ifstream InputFile;
  InputFile.open(FileName.c_str());
  if (!InputFile.good() && CheckInputs) {
    std::cerr << " The file: " << FileName.c_str() << " doesn't exist!. Please check the input file(s)!" << std::endl;
    return true;
  }
  else if (InputFile.good() && !CheckInputs) {
    std::cerr << " The output FileName, " << outputFileName << ", already exists! Please use another name!" << std::endl;
    return true;//ToBeReviewed!
  }
  InputFile.close();
  return false;
}

inline bool CheckIfTheExtIsRoot(std::string FileName){
  
  // If the file extension is not "root" the function will return false
  if(FileName.substr(FileName.find_last_of(".") + 1) == "root")
    return true;
  
  return false;
}

inline bool CheckingInputs(int argc,char* argv[]){
  
  tree = new TChain("ids");
  
  if (argc == 2){
    inputFileName  = argv[1];
    if ( CheckIfFileExist(inputFileName, true) ) return false;
    outputFileName = inputFileName.substr(0, inputFileName.size()-5);
    outputFileName = outputFileName + "_histograms.root" ;
    tree->Add(argv[1]);
    std::cout << std::endl << " Input file:  " << argv[1] << std::endl;
  }
  else if (argc > 2){
    outputFileName = argv[1];
    if ( !CheckIfTheExtIsRoot(outputFileName) )
      outputFileName = outputFileName + ".root";
    if ( CheckIfFileExist(outputFileName, false) ){
      std::cout << std::endl << " Wrong output! -> It exists an output with the same name." << std::endl;
      return false;
    }
    for(int i = 2; i < argc; ++i) {
      inputFileName = argv[i];
      if ( i == 2){
        std::cout << std::endl << " Input files: " << inputFileName << std::endl;
      }
      else{
        std::cout << "            : " << inputFileName << std::endl;
      }
      if ( CheckIfFileExist(inputFileName, true) ){
        std::cout << " Wrong imput! -> Last input does not exist." << std::endl;
        return false;
      }
      tree->Add(inputFileName.c_str());
    }
  }
  else{
    return false;
  }
  
  std::cout << " Output file: " << outputFileName << std::endl;
  std::cout << "" << std::endl;
  
  return true;
}


// Main code
int main(int argc,char* argv[]){
  
  clock_t tStart = clock();
  
  if(!CheckingInputs(argc, argv)){
    std::cout << std::endl << " Examples:" << std::endl;
    std::cout << "   " << argv[0] << " <analysis tree file>" << std::endl;
    std::cout << "   " << argv[0] << " <OutputFileName.root> <analysis tree file(s)>" << std::endl << std::endl;
    return 0;
  }
  
  std::vector<std::string> ListOfElements;
  
  // ids tree info
  // Declaration of leaves types
  Int_t           MULT;
  ULong64_t       TIME_REF;
  ULong64_t       TIMESTAMP;
  Double_t        E_Clover[num_Clov_crys];
  Double_t        T_Clover[num_Clov_crys];
  Int_t           M_Clover;
  Double_t        E_Beta[num_beta];
  Double_t        T_Beta[num_beta];
  Int_t           M_Beta;
  Double_t        E_SPEDE[num_spede];
  Double_t        T_SPEDE[num_spede];
  Int_t           M_SPEDE;

  // Set branch addresses.
  tree->SetBranchAddress("Time_vs_ref",&TIME_REF);
  tree->SetBranchAddress("Timestamp",  &TIMESTAMP);

  bool gotClov = false;
  if(tree->FindBranch("E_Clov") == nullptr) { // We check to see if we have a clover branch in the analysis tree
    gotClov = false;
  } else {
    tree->SetBranchAddress("E_Clov", &E_Clover);
    tree->SetBranchAddress("T_Clov", &T_Clover);
    tree->SetBranchAddress("M_Clov", &M_Clover);
    gotClov = true;
    std::cout << " Setting the Clover branch " << std::endl;
    ListOfElements.push_back("Clovers");
  }
  
  bool gotSPEDE = false;
  if(tree->FindBranch("E_SPEDE") == nullptr) { // We check to see if we have a SPEDE branch in the analysis tree
    gotSPEDE = false;
  } else {
    tree->SetBranchAddress("E_SPEDE", &E_SPEDE);
    tree->SetBranchAddress("T_SPEDE", &T_SPEDE);
    tree->SetBranchAddress("M_SPEDE", &M_SPEDE);
    gotSPEDE = true;
    std::cout << " Setting the SPEDE branch " << std::endl;
    ListOfElements.push_back("SPEDE");
  }
  
  bool gotBeta = false;
  if(tree->FindBranch("E_Beta") == nullptr) { // We check to see if we have a beta branch in the analysis tree
      gotBeta = false;
  } else {
    tree->SetBranchAddress("E_Beta",&E_Beta);
    tree->SetBranchAddress("T_Beta",  &T_Beta);
    tree->SetBranchAddress("M_Beta",  &M_Beta);
    gotBeta = true;
    std::cout << " Setting the beta branch " << std::endl;
    ListOfElements.push_back("Betas");
  }

  
  //Global variables
  int time_diff, time_diff2;

  //Here we create the histograms. We will add them to a TList, so they are all written at the same time.
  // Clover histograms
  auto* list_Clov = new TList; //list of the clov histograms
  TH1F* Clov_time                    = new TH1F("Clov_time","Clover-clover time difference; Time 4 ns/Chan; Counts/4 ns",800,-400.,400.); list_Clov->Add(Clov_time);
  TH1F* Clov_singles                 = new TH1F("Clov_singles","Clover singles; E_{#gamma} [keV]; Counts/0.25 keV",32000,0.,8000.); list_Clov->Add(Clov_singles);
  TH1F* Clov_singles_betagated       = new TH1F("Clov_singles_betagated","Clover singles beta gated; E_{#gamma} [keV]; Counts/0.25 keV",32000,0.,8000.); list_Clov->Add(Clov_singles_betagated);
  TH2F* Clov_Id_vs_Ener              = new TH2F("Clov_Id_vs_Ener","Clover ID Vs Clover energy; ID number; E_{#gamma} [1.0 keV/Chan]",30,0,30,8000,0.,8000.); list_Clov->Add(Clov_Id_vs_Ener);
  TH2F* Clov_gg                      = new TH2F("Clov_gg","Clover-Clover with time random subtracted; E_{#gamma} [1.0 keV/Chan]; E_{#gamma} [1.0 keV/Chan]",8000,0.,8000.,8000,0.,8000.); list_Clov->Add(Clov_gg);
  TH2F* Clov_gg_Notrand              = new TH2F("Clov_gg_Notrand","Clover-Clover without time-random coincidences; E_{#gamma} [1.0 keV/Chan]; E_{#gamma} [1.0 keV/Chan]",8000,0.,8000.,8000,0.,8000.); list_Clov->Add(Clov_gg_Notrand);
  TH2F* Clov_gg_trand                = new TH2F("Clov_gg_trand","Clover-Clover time-random coincidences; E_{#gamma} [1.0 keV/Chan]; E_{#gamma} [1.0 keV/Chan]",8000,0.,8000.,8000,0.,8000.); list_Clov->Add(Clov_gg_trand);
  TH2F* Clov_gg_betagated            = new TH2F("Clov_gg_betagated","Clover-Clover with time random subtracted beta gated; E_{#gamma} [1.0 keV/Chan]; E_{#gamma} [1.0 keV/Chan]",8000,0.,8000.,8000,0.,8000.); list_Clov->Add(Clov_gg_betagated);
  TH2F* Clov_gg_Notrand_betagated    = new TH2F("Clov_gg_Notrand_betagated","Clover-Clover without time-random coincidences beta gated;  E_{#gamma} [1.0 keV/Chan]; E_{#gamma} [1.0 keV/Chan]",8000,0.,8000.,8000,0.,8000.); list_Clov->Add(Clov_gg_Notrand_betagated);
  TH2F* Clov_gg_trand_betagated      = new TH2F("Clov_gg_trand_betagated","Clover-Clover time-random coincidences beta gated;  E_{#gamma} [1.0 keV/Chan]; E_{#gamma} [1.0 keV/Chan]",8000,0.,8000.,8000,0.,8000.); list_Clov->Add(Clov_gg_trand_betagated);
  TH2F* Clov_vs_TRef                 = new TH2F("Clov_vs_TRef",Form("TimeRef vs Clovers; Time [%d %s/Chan]; E_{#gamma} [keV]",(int)RefUnit_val,RefUnit_str.c_str()),MaxRefUnit,0,MaxRefUnit, 8000,0.,8000.); list_Clov->Add(Clov_vs_TRef);
  TH2F* Clov_vs_TRef_betagated       = new TH2F("Clov_vs_TRef_betagated",Form("TimeRef vs Clovers beta gated; Time [%d %s/Chan]; E_{#gamma} [keV]",(int)RefUnit_val,RefUnit_str.c_str()),MaxRefUnit,0,MaxRefUnit, 8000,0.,8000.); list_Clov->Add(Clov_vs_TRef_betagated);
  
  TH2F* Clov_Beta_tdiff_vs_Ener;
  if(gotClov && gotBeta){
    Clov_Beta_tdiff_vs_Ener = new TH2F("Clov_Beta_tdiff_vs_Ener","Time difference Beta-Clover Vs Clover energy; Time 4 ns/Chan; E_{#gamma} [1.0 keV/Chan]",2000,-1000.,1000.,3000,0.,3000.);
    list_Clov->Add(Clov_Beta_tdiff_vs_Ener);
  }
  
  TH2F *Clov_SPEDE_Energy_ge, *Clov_SPEDE_tdiff_vs_Ener;
  if(gotClov && gotSPEDE){
    Clov_SPEDE_tdiff_vs_Ener = new TH2F("Clov_SPEDE_tdiff_vs_Ener","Time difference SPEDE-Clover Vs Clover energy; Time 4 ns/Chan; E_{#gamma} [1.0 keV/Chan]",800,-400.,400.,4000,0.,4000.);
    list_Clov->Add(Clov_SPEDE_tdiff_vs_Ener);
    Clov_SPEDE_Energy_ge     = new TH2F("Clov_SPEDE_Energy_ge","Clover Vs SPEDE energy {#gamma}e; E_{SPEDE} [1.0 keV/Chan]; E_{#gamma} [1.0 keV/Chan]",1000,0.,4000.,4000,0.,4000.);
    list_Clov->Add(Clov_SPEDE_Energy_ge);
  }
  
  // SPEDE histograms
  auto* list_SPEDE = new TList; //list of the Beta histograms
  auto* list_SPEDE2 = new TList; //list of the Beta histograms
  TH1F* SPEDEDet[num_spede];
  TH1F *SPEDE_Energy;
  TH2F *SPEDE_EnergyVsEnergy;
  SPEDE_Energy         = new TH1F("SPEDE_Energy","SPEDE Energy; E_{SPEDE} [keV]; Counts/4.0 keV",SPEDE_energy_range[0],0.,(double)SPEDE_energy_range[1]);
  list_SPEDE->Add(SPEDE_Energy);
  SPEDE_EnergyVsEnergy = new TH2F("SPEDE_EnergyVsEnergy","SPEDE Vs SPEDE Energy; E_{SPEDE} [4.0 keV/bin]; E_{SPEDE} [4.0 keV/bin]",(double)SPEDE_energy_range[0],0.,(double)SPEDE_energy_range[1],SPEDE_energy_range[0],0.,(double)SPEDE_energy_range[1]);
  list_SPEDE->Add(SPEDE_EnergyVsEnergy);
  for(int i = 0; i < num_spede; i++){
    SPEDEDet[i] = new TH1F(Form("SPEDEDet_%d",i),Form("SPEDE %d; E_{SPEDE} [keV]; Counts/4.0 keV",i),SPEDE_energy_range[0],0.,(double)SPEDE_energy_range[1]);
    list_SPEDE2->Add(SPEDEDet[i]);
  }
  
  // Beta histograms
  auto* list_Beta = new TList; //list of the Beta histograms
  TH1F* betaDet[num_beta];
  for(int i = 0; i < num_beta+1; i++){
    betaDet[i] = new TH1F(Form("BetaDet_%d",i),Form("Beta channel %d; E_{beta} [Chan]; Counts/10 Chan",i),100000,0.,100000.);
    list_Beta->Add(betaDet[i]);
  }
  
  // TimeAlignment histograms
  auto* list_Time = new TList; //list of the TimeAlignment histograms
  TH1F *Clov_Beta_Time, *Clov_SPEDE_Time, *SPEDE_SPEDE_Time;
  TH2F *Clov_Id_vs_TimeAlignment, *Beta_Id_vs_TimeAlignment, *SPEDE_Id_vs_TimeAlignment;
  
  if(gotClov){
    Clov_Id_vs_TimeAlignment = new TH2F("Clov_Id_vs_TimeAlignment","Clover ID Vs Time Alignment; ID number; Time diference [4 ns/Chan]",30,0,30,2000,-1000,1000); 
    list_Time->Add(Clov_Id_vs_TimeAlignment);
  }
  if(gotClov && gotSPEDE){
    Clov_SPEDE_Time = new TH1F("Clov_SPEDE_time","Clover-SPEDE time difference; Time 1 ns/Chan; Counts/1 ns",20000,-10000.,10000.); list_Time->Add(Clov_SPEDE_Time);
    SPEDE_SPEDE_Time = new TH1F("SPEDE_SPEDE_Time","SPEDE-SPEDE time difference; Time 1 ns/Chan; Counts/1 ns",20000,-10000.,10000.); list_Time->Add(SPEDE_SPEDE_Time);
    SPEDE_Id_vs_TimeAlignment = new TH2F("SPEDE_Id_vs_TimeAlignment","SPEDE ID Vs Time Alignment; ID number; Time diference [1 ns/Chan]",num_spede,0,num_spede,2000,-10000,10000); list_Time->Add(SPEDE_Id_vs_TimeAlignment);
  }
  if(gotClov && gotBeta){
    Clov_Beta_Time = new TH1F("Clov_Beta_time","Clover-Beta time difference; Time 1 ns/Chan; Counts/1 ns",2000,-1000.,1000.);
    list_Time->Add(Clov_Beta_Time);
    Beta_Id_vs_TimeAlignment = new TH2F("Beta_Id_vs_TimeAlignment","Beta ID Vs Time Alignment; ID number; Time diference [4 ns/Chan]",num_beta,0,num_beta,200,-1000,1000); list_Time->Add(Beta_Id_vs_TimeAlignment);
  }
  
  // Others
  auto* list_others = new TList; //list of the extra histograms
  TH1F* TRef = new TH1F("TRef",Form("Time difference; Time [%d %s/Chan]; Counts/Chan",(int)RefUnit_val,RefUnit_str.c_str()),MaxRefUnit,0,MaxRefUnit); list_others->Add(TRef);
  
  //60Co analysis
  auto* list_60CoAna = new TList; //list of the extra histograms
  TH2F *Clov_Id_vs_TimeAlignment60Co, *Clov_TimeAlignment_Vs_Energy;
  if (Do60CoCase){
    Clov_Id_vs_TimeAlignment60Co = new TH2F("Clov_Id_vs_TimeAlignment60Co","Clover ID Vs Time Alignment - Gate on 1332 keV gamma line; ID number; Time diference [4 ns/Chan]",30,0,30,2000,-1000,1000); list_60CoAna->Add(Clov_Id_vs_TimeAlignment60Co);
    Clov_TimeAlignment_Vs_Energy = new TH2F("Clov_TimeAlignment_Vs_Energy","Clover Time Alignment Vs Energy - Gate on 1332 keV gamma line; ID number; Time diference [4 ns/Chan]",2000,-1000,1000,4000,0,4000); list_60CoAna->Add(Clov_TimeAlignment_Vs_Energy);
  }
  
  //Here is where we start looping over the events
  Long64_t nentries = tree->GetEntries();
  std::cout << std::endl;
  std::cout << " Filling ";
  for (int i = 0; i < ListOfElements.size(); i++){
    if (ListOfElements.size() == 1)
      std::cout << ListOfElements[i];
    else if (i == ListOfElements.size() - 1)
      std::cout << "and " << ListOfElements[i];
    else
      std::cout << ListOfElements[i] << ", ";
  }
  std::cout << std::endl;
  std::cout << " Total events: " << nentries/1e6 << " millions" << std::endl;
  
  Long64_t totalEvents=0;
  float progress = 0.0, progress_buffer = 0.0;
  
  int Beta_count=0;
  ULong64_t TimeRef_value = 0;
  
  std::vector<double> tmp_BetaEn;
  std::vector<double> tmp_BetaTi;
  
  tree->SetImplicitMT(true);
  for (Long64_t i = 0; i < nentries; i++) {//Loops over all the events (entries) of the tree
    tree->GetEntry(i);
    
    Beta_count = 0;
    tmp_BetaEn.clear();
    tmp_BetaTi.clear();
    
    if(gotBeta && M_Beta > 0){
      for(Int_t j = 0; j < num_beta; j++){
        if (E_Beta[j] > BetaPMT_energy_gate[0] && E_Beta[j] < BetaPMT_energy_gate[1]){
          tmp_BetaEn.push_back(E_Beta[j]);
          tmp_BetaTi.push_back(T_Beta[j]);
          betaDet[j]->Fill(E_Beta[j]);
          Beta_count += 1;
        }
      }
    }//Beta boolean check
    
    
    if(gotSPEDE && M_SPEDE > 0){
      for(Int_t j = 0; j < num_spede; j++){
        if (E_SPEDE[j] > SPEDE_energy_gate){
          SPEDE_Energy->Fill(E_SPEDE[j]);
          SPEDEDet[j]->Fill(E_SPEDE[j]);
        }
      }
    }//SPEDE boolean check
    
    
    if(gotClov && M_Clover > 0){//Check for the very unlikely event we run without clovers
      
      if(DoAddback) addback(E_Clover, T_Clover, M_Clover);
      
      for(Int_t j = 0; j < num_Clov_crys; j++){//Loop over first clover event
        if(E_Clover[j] > EClov_Min){
          Clov_singles->Fill(E_Clover[j]);
          Clov_Id_vs_Ener->Fill(j,E_Clover[j]);
          if (TIME_REF < MaxRefUnit){
            Clov_vs_TRef->Fill(TIME_REF,E_Clover[j]);
            TRef->Fill(TIME_REF);
          }
          if(Beta_count > 0)
            for(Int_t k = 0; k < tmp_BetaEn.size(); k++){
              time_diff = tmp_BetaTi.at(k) - T_Clover[j];
              Clov_Beta_Time->Fill(time_diff);
              if( tmp_BetaEn.at(k) > BetaPMT_energy_gate[0] && tmp_BetaEn.at(k) < BetaPMT_energy_gate[1]) Clov_Beta_tdiff_vs_Ener->Fill(-1.*time_diff, E_Clover[j]);
              if( tmp_BetaEn.at(k) > BetaPMT_energy_gate[0] && tmp_BetaEn.at(k) < BetaPMT_energy_gate[1] && (time_diff > Clov_Beta_t_coinc[0] && time_diff < Clov_Beta_t_coinc[1])){
                Clov_singles_betagated->Fill(E_Clover[j]);
                if (TIME_REF < MaxRefUnit)
                  Clov_vs_TRef_betagated->Fill(TIME_REF,E_Clover[j]);
              }
            }
          
          if (gotBeta && M_Beta > 0)
            for (Int_t k = 0; k < num_beta; k++)
              if (j == RefIDCrys && E_Beta[k] > 0) Beta_Id_vs_TimeAlignment->Fill(k, T_Beta[k] - T_Clover[j]);
          
          if(M_Clover > 1){
            for(Int_t k = 0; k < num_Clov_crys; k++){//Loop over second and beyond clover envents
              if(E_Clover[k] > EClov_Min && j != k ){
                time_diff = T_Clover[k] - T_Clover[j];
                Clov_time->Fill(time_diff);
                if (j == RefIDCrys) Clov_Id_vs_TimeAlignment->Fill(k, time_diff);
                if (Do60CoCase && j == RefIDCrys && E_Clover[j] > 1328 && 1336 > E_Clover[j]){
                  Clov_Id_vs_TimeAlignment60Co->Fill(k, time_diff);
                  Clov_TimeAlignment_Vs_Energy->Fill(time_diff,  E_Clover[k]);
                }
                if(abs(time_diff) < Clov_t_coinc){//True time coincidences
                  Clov_gg->Fill(E_Clover[j], E_Clover[k]);
                  if(Beta_count > 0)
                    for(Int_t l = 0; l < tmp_BetaEn.size(); l++){
                      time_diff2 = tmp_BetaTi.at(l) - T_Clover[j];
                      if( tmp_BetaEn.at(l) > BetaPMT_energy_gate[0] && tmp_BetaEn.at(l) < BetaPMT_energy_gate[1] && (time_diff2 >= Clov_Beta_t_coinc[0] && time_diff2 <= Clov_Beta_t_coinc[1]))
                        Clov_gg_betagated->Fill(E_Clover[j], E_Clover[k]);
                    }
                }
                else if(Clov_rand_coinc[0] < abs(time_diff) && abs(time_diff) < Clov_rand_coinc[1]){//Time-random coincidences
                  Clov_gg_trand->Fill(E_Clover[j], E_Clover[k]);
                  if(Beta_count > 0)
                    for(Int_t l = 0; l < tmp_BetaEn.size(); l++){
                      time_diff2 = tmp_BetaTi.at(l) - T_Clover[j];
                      if( tmp_BetaEn.at(l) > BetaPMT_energy_gate[0] && tmp_BetaEn.at(l) < BetaPMT_energy_gate[1] && (time_diff2 >= Clov_Beta_t_coinc[0] && time_diff2 <= Clov_Beta_t_coinc[1]))
                        Clov_gg_trand_betagated->Fill(E_Clover[j], E_Clover[k]);
                    }
                }
              }
            }
            if (M_SPEDE > 0){
              for(Int_t k = 0; k < num_spede; k++){//Loop over spede
                if(E_SPEDE[k] > SPEDE_energy_gate){
                  time_diff = T_SPEDE[k] - T_Clover[j];
                  Clov_SPEDE_Time->Fill(time_diff);
                  if (j == RefIDCrys) SPEDE_Id_vs_TimeAlignment->Fill(k, time_diff);
                  Clov_SPEDE_tdiff_vs_Ener->Fill(time_diff, E_SPEDE[k]);
                  if(time_diff >= Clov_SPEDE_t_coinc[0] && time_diff <= Clov_SPEDE_t_coinc[1]){//True time coincidences
                    Clov_SPEDE_Energy_ge->Fill(E_Clover[j], E_SPEDE[k]);
                  }
                }
              }
            }
          }
        }
      }
    }//Clover boolean check
    
    for(Int_t j = 0; j < num_spede; j++){//Loop over spede
      if(E_SPEDE[j] > SPEDE_energy_gate){
        for(Int_t k = 0; k < num_spede; k++){//Loop over spede
          if(j != k && E_SPEDE[k] > SPEDE_energy_gate){
            time_diff = T_SPEDE[k] - T_SPEDE[j];
            SPEDE_SPEDE_Time->Fill(time_diff);
            if(time_diff >= SPEDE_SPEDE_t_coinc[0] && time_diff <= SPEDE_SPEDE_t_coinc[1]){//True time coincidences
              SPEDE_EnergyVsEnergy->Fill(E_SPEDE[j], E_SPEDE[k]);
	    }
          }
	}
      }
    }
              

    
    progress = ((float)i/(float)nentries)*100.0;
    if (progress >= progress_buffer || i == nentries){
      PrintBar(progress);
      std::cout.flush();
      progress_buffer = progress + 1.0;
    }
    
  }// End LOOP OVER TTREE EVENTS
  
  progress = progress + 1.0;
  PrintBar(progress);
  std::cout.flush();
  std::cout << std::endl << std::endl << " Sorting finished!" << std::endl;
  std::cout << " Saving the histograms...";
  
  //We subtract the time random coincidences from the true ones:
  if(gotClov){
    Clov_gg_Notrand->Add(Clov_gg,1);
    if(gotBeta) Clov_gg_Notrand_betagated->Add(Clov_gg_betagated,1);
    Clov_gg->Add(Clov_gg_trand,-1.*Clov_t_coinc/(Clov_rand_coinc[1] - Clov_rand_coinc[0]));
    if(gotBeta) Clov_gg_betagated->Add(Clov_gg_trand_betagated,-1.*Clov_t_coinc/(Clov_rand_coinc[1] - Clov_rand_coinc[0]));
  }
  
  //We create the output file and write the histograms.
  TFile *fout = new TFile(outputFileName.c_str(),"RECREATE");
  
  list_others->Write();
  if(gotClov) list_Clov->Write();                 // We only write the histos if we had clovers
  if(gotBeta) list_Beta->Write();                 // We only write the histos if we had Betas
  if(gotSPEDE) {
    list_SPEDE->Write();
    TDirectory *SPEDEDir_dir = fout->mkdir("Idiv_SPEDE");
    SPEDEDir_dir->cd();
    list_SPEDE2->Write();                           // We only write the histos if we had SPEDE
    gDirectory->cd("/");
  }
  if (list_Time->GetSize() > 0){
    TDirectory *TimeDir_dir = fout->mkdir("TimeAlignment");
    TimeDir_dir->cd();
    list_Time->Write();                           // We only write the histos if we had All
    gDirectory->cd("/");
  }
  if(gotClov && Do60CoCase) list_60CoAna->Write();  

  fout->Close();
	
  std::cout << " done!" << std::endl << std::endl;
  
  printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  
  return 0;
}
