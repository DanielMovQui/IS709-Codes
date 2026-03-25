//
//Written by: D.Movilla 2025-10-01
//Experiment: IS622
//Root version: 6.14/04
//This script has a TChain, so you can give a list of files and it will loop through all of them
//
/////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// To compile:
//
// g++ histo_angular_corr_def10.C -Wall `root-config --cflags` `root-config --glibs` -o histo_angular_corr_def10
//
// ./histo_angular_corr_def2 Run362.root
//
/////////////////////////////////////////////////////////////////////////////////////////////////


//C,C++ Libraries
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <set>
#include <map>
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
#include <TChain.h>
#include <iomanip>
#include "TVector3.h"
#include "TMath.h"



using namespace std;

// Global variables
TChain* tree;
std::string inputFileName, outputFileName;
//static int TimeCoincAddback = 120;

// Running flags
const bool DoAddback  = false;
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
const double Clov_t_coinc          = 240.;
const double Clov_rand_coinc[2]    = { 260., 1500.};
// Clover-Beta Time window
const double Clov_Beta_t_coinc[2]  = { -200.,   275.};
// Clover-SPEDE Time window
const double Clov_SPEDE_t_coinc[2] = { -200.,   60.};
const double SPEDE_SPEDE_t_coinc[2] = { -200.,   60.};

// Clover Energy min
const double EClov_Min = 1.;

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
const int         MaxRefUnit  = 20000; // MaxRefUnit*RefUnit_val in RefUnit_str

////// Here finish the configuration for the experiment //////

// NUMBER OF HISTOGRAMS FOR ANGULAR CORRELATIONS, GIVEN BY ANGULAR RESOLUTION

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
    //outputFileName = outputFileName + "_angular_correlations_list_crys10_287_129_0.33.root" ;
    outputFileName = outputFileName + "_angular_correlations_list_crys10_0.5_t_ref_corr_new_positions.root" ;
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
  //Int_t           MULT;
  ULong64_t       TIME_REF;
  ULong64_t       TIMESTAMP;
  Double_t        E_Clover[num_Clov_crys];
  //Double_t        E_Clover_Mixed[num_Clov_crys];
  Double_t        T_Clover[num_Clov_crys];
  Int_t           M_Clover;
  //Double_t        E_Beta[num_beta];
  //Double_t        T_Beta[num_beta];
  //Int_t           M_Beta;
  //Double_t        E_SPEDE[num_spede];
  //Double_t        T_SPEDE[num_spede];
 // Int_t           M_SPEDE;

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

  //Global variables

// ============================== 
// DEFINITION OF THE CRYSTALS
// ==============================

const int n_gantries = 5;
const int n_colors = 4;

const std::string gantries[n_gantries] = { "g1", "g2", "g3", "g4", "g5" };
const std::string colors[n_colors] = { "yellow", "green", "blue", "red" };

std::vector<std::string> clover;
std::vector<std::string> color;

for(int i = 0; i < n_gantries; ++i){
    std::vector<char> labels;
    if(i == 0 || i == 4){
        labels = {'A', 'B'};
    }
    else {
        labels = {'A', 'B', 'C'};
    }

    for(char clover_label : labels) {
        char label_to_use = clover_label;

        if (i == 1) {
            if (clover_label == 'A') label_to_use = 'C';
            else if (clover_label == 'C') label_to_use = 'A';
        }

        std::string name = gantries[i] + "_" + label_to_use;
        for(int c = 0; c < n_colors; ++c){
            clover.push_back(name);
            color.push_back(colors[c]);
        }
    }
}

// =========================================
// LECTURA DE POSICIONES 3D DE LOS CRISTALES
// =========================================
std::map<int, TVector3> crystal_pos;
{
    std::ifstream file_pos("crystal_positions_def_final.txt");
  //  std::ifstream file_pos("crystals.txt");

    if (!file_pos.is_open()) {
        std::cerr << "Error: no se pudo abrir crystal_positions_def_final.txt\n";
        exit(1);
    }

    std::string line;
    int idx = 0;
    while (std::getline(file_pos, line)) {
        size_t colon = line.find(':');
        if (colon == std::string::npos) continue;

        std::string coords = line.substr(colon + 1);
        double x, y, z;
        std::stringstream ss(coords);
        ss >> x; ss.ignore(1, ','); ss >> y; ss.ignore(1, ','); ss >> z;
        crystal_pos[idx] = TVector3(x, y, z);
        idx++;
    }
    std::cout << "Loaded " << crystal_pos.size() << " crystal positions.\n";
}

double Emax = 2000.;
double Emin = 0.; 
double binning = 4000;

auto* list_Clov = new TList; //list of the clov histograms
  TH1F* Clov_time                    = new TH1F("Clov_time","Clover-clover time difference; Time 4 ns/Chan; Counts/4 ns",4000,-2000.,2000.); list_Clov->Add(Clov_time);
  TH1F* Clov_singles                 = new TH1F("Clov_singles","Clover singles; E_{#gamma} [keV]; Counts/0.5 keV",binning,Emin,Emax); list_Clov->Add(Clov_singles);
  TH2F* Clov_Id_vs_Ener              = new TH2F("Clov_Id_vs_Ener","Clover ID Vs Clover energy; ID number; E_{#gamma} [0.5 keV/Chan]",52,0,52,binning,Emin,Emax); list_Clov->Add(Clov_Id_vs_Ener);
  TH2F* Clov_gg                      = new TH2F("Clov_gg","Clover-Clover with time random subtracted; E_{#gamma} [0.5 keV/Chan]; E_{#gamma} [0.5 keV/Chan]",binning,Emin,Emax,binning,Emin,Emax); list_Clov->Add(Clov_gg);
  TH2F* Clov_gg_Notrand              = new TH2F("Clov_gg_Notrand","Clover-Clover without time-random coincidences; E_{#gamma} [0.5 keV/Chan]; E_{#gamma} [0.5 keV/Chan]",binning,Emin,Emax,binning,Emin,Emax); list_Clov->Add(Clov_gg_Notrand);
  TH2F* Clov_gg_trand                = new TH2F("Clov_gg_trand","Clover-Clover time-random coincidences; E_{#gamma} [0.5 keV/Chan]; E_{#gamma} [0.5 keV/Chan]",binning,Emin,Emax,binning,Emin,Emax); list_Clov->Add(Clov_gg_trand);
  TH2F* Clov_vs_TRef                 = new TH2F("Clov_vs_TRef",Form("TimeRef vs Clovers; Time [%d %s/Chan]; E_{#gamma} [keV]",(int)RefUnit_val,RefUnit_str.c_str()),MaxRefUnit,0,MaxRefUnit, binning,Emin,Emax); list_Clov->Add(Clov_vs_TRef);



// =========================================
// CONFIGURACIÓN DE BINS EN COS(THETA)
// =========================================
const int n_cos_bins = 20;         // número de bins
const double cos_min = -1.0;
const double cos_max = 1.0;
const double cos_bin_width = (cos_max - cos_min) / n_cos_bins;

auto* list_crys = new TList();
auto* list_crys_bckg = new TList();
auto* list_crys_true = new TList();
auto* list_crys_mixed = new TList();

// Crear vectores de histogramas: True, Background y Mixed
std::vector<TH2F*> h_cosbins(n_cos_bins);
std::vector<TH2F*> h_cosbins_bckg(n_cos_bins);
std::vector<TH2F*> h_cosbins_mixed(n_cos_bins);

for (int b = 0; b < n_cos_bins; ++b) {
    double low = cos_min + b * cos_bin_width;
    double high = low + cos_bin_width;

    // ---------------- TRUE ----------------
    TString name_true = Form("Angular_correlation_hist_True_binned_%01d", b);
    TString title_true = Form("True coincidences: E_Clover[j] vs E_Clover[k] for cosθ ∈ [%.2f, %.2f)", low, high);
    h_cosbins[b] = new TH2F(name_true, title_true, 4000, 0, 2000, 4000, 0, 2000);
    list_crys->Add(h_cosbins[b]);

    // ---------------- BACKGROUND ----------------
    TString name_bckg = Form("Angular_correlation_hist_%01d", b);
    TString title_bckg = Form("Random coincidences: E_Clover[j] vs E_Clover[k] for cosθ ∈ [%.2f, %.2f)", low, high);
    h_cosbins_bckg[b] = new TH2F(name_bckg, title_bckg, 4000, 0, 2000, 4000, 0, 2000);
    list_crys_bckg->Add(h_cosbins_bckg[b]);

    // ---------------- MIXED ----------------
    TString name_mixed = Form("Angular_correlation_hist_Mixed_binned_%01d", b);
    TString title_mixed = Form("Mixed events: E_Clover[j] vs E_Clover[k] for cosθ ∈ [%.2f, %.2f)", low, high);
    h_cosbins_mixed[b] = new TH2F(name_mixed, title_mixed, 4000, 0, 2000, 4000, 0, 2000);
    list_crys_mixed->Add(h_cosbins_mixed[b]);
}




  //Here is where we start looping over the events
  Long64_t nentries = tree->GetEntries();
  std::cout << std::endl;
  std::cout << " Filling ";
  
  std::cout << std::endl;
  std::cout << " Total events: " << nentries/1e6 << " millions" << std::endl;
  
  float progress = 0.0, progress_buffer = 0.0;
  
  std::vector<double> tmp_BetaEn;
  std::vector<double> tmp_BetaTi;
  
  tree->SetImplicitMT(true);
  int event_separation = 10;

//Co60 Gates
double E_low1 = 1155., E_high1 = 1195.; // rango 1
double E_low2 = 1310., E_high2 = 1350.; // rango 2



//Rb100 Gates
/*
double E_low1 = 100., E_high1 = 150.; // rango 1
double E_low2 = 250., E_high2 = 325.; // rango 2
*/
/*
//double E_low2 = 750., E_high2 = 850.; // rango 2
double E_low1 = 0., E_high1 = 2000.; // rango 1
double E_low2 = 0., E_high2 = 2000.; // rango 2
*/
  // ==================
  // FILLING HISTOGRAMS 
  // ==================

for (Long64_t i = 0; i < nentries - event_separation; ++i) {

  // --- Read entry i and save copies ---
  tree->GetEntry(i);

  ULong64_t TIME_REF_original = TIME_REF;  // <<<<< guardar antes de modificar entrada

for (int j = 36; j <= 43; ++j) {
    T_Clover[j] -= (-226.0);  // restamos el offset (equivalente a sumar 226 ns)
}
  const std::vector<double> E_Clover_original(E_Clover, E_Clover + num_Clov_crys);
  const std::vector<double> T_Clover_original(T_Clover, T_Clover + num_Clov_crys);
  int M_Clover_original = M_Clover;   // <-- guardar multiplicidad original
  bool gotClov_original = gotClov;

  // --- Read entry i + event_separation and save copies ---
  tree->GetEntry(i + event_separation);
  const std::vector<double> E_Clover_next(E_Clover, E_Clover + num_Clov_crys);
  // No necesitamos T_Clover_next para el mixed simple, pero puedes guardarlo si quieres:
  // const std::vector<double> T_Clover_next(T_Clover, T_Clover + num_Clov_crys);
  int M_Clover_next = M_Clover;
  bool gotClov_next = gotClov;
  ULong64_t TIME_REF_mixed = TIME_REF;

  tmp_BetaEn.clear();
  tmp_BetaTi.clear();

  // --- Procesado del evento 'original' (singles + coincidencias dentro del mismo evento) ---
  if (gotClov_original && M_Clover_original > 0) {
bool tref_ok_original = true;
bool tref_ok_mixed   = true;
/*
// Primer salto: 0 → 2400
if (TIME_REF_original > 15 && TIME_REF_original < 300)
    tref_ok_original = true;
if (TIME_REF_mixed > 15 && TIME_REF_mixed < 300)
    tref_ok_mixed = true;

// Saltos restantes: 2400, 3600, ..., 8400
for (ULong64_t base = 2400; base <= 8400; base += 1200) {
    if (TIME_REF_original > base + 15 && TIME_REF_original < base + 300)
        tref_ok_original = true;
    if (TIME_REF_mixed > base + 15 && TIME_REF_mixed < base + 300)
        tref_ok_mixed = true;
}
*/


    if (tref_ok_original){

    for (Int_t j = 0; j < num_Clov_crys; ++j) {
      if (E_Clover_original[j] <= EClov_Min) continue;
        Clov_singles->Fill(E_Clover_original[j]);
	      Clov_Id_vs_Ener->Fill(j,E_Clover_original[j]);
        if (TIME_REF < MaxRefUnit){                
            Clov_vs_TRef->Fill(TIME_REF,E_Clover_original[j]);
          }

      // TRUE / RANDOM coincidences: sólo si hay al menos 2 cristales en el evento original
      if (M_Clover_original > 1) {
        for (Int_t k = 0; k < num_Clov_crys; ++k) {
          if (j == k) continue;
          if (E_Clover_original[k] <= EClov_Min) continue;
          double time_diff = T_Clover_original[j] - T_Clover_original[k];
          Clov_time->Fill(time_diff);
          if (std::abs(time_diff) < Clov_t_coinc) {
            // true coincidence
            bool inRange1_j = (E_Clover_original[j] >= E_low1 && E_Clover_original[j] <= E_high1);
            bool inRange2_j = (E_Clover_original[j] >= E_low2 && E_Clover_original[j] <= E_high2);
            bool inRange1_k = (E_Clover_original[k] >= E_low1 && E_Clover_original[k] <= E_high1);
            bool inRange2_k = (E_Clover_original[k] >= E_low2 && E_Clover_original[k] <= E_high2);
            bool cross_event = (inRange1_j && inRange2_k) || (inRange2_j && inRange1_k);

            if (cross_event) {
            Clov_gg->Fill(E_Clover_original[j], E_Clover_original[k]);

        // ======================
        // Calcular cos(θ)
        // ======================
        if (crystal_pos.count(j) && crystal_pos.count(k)) {
            TVector3 vj = crystal_pos[j];
            TVector3 vk = crystal_pos[k];

            double cos_theta = vj.Dot(vk) / (vj.Mag() * vk.Mag());
            cos_theta = std::max(-1.0, std::min(1.0, cos_theta)); // evitar errores numéricos

            // ======================
            // Determinar bin
            // ======================
            int bin_index = int((cos_theta - cos_min) / cos_bin_width);
            if (bin_index >= 0 && bin_index < n_cos_bins) {
                h_cosbins[bin_index]->Fill(E_Clover_original[j], E_Clover_original[k]);
            }
        }
            }
          }
          else if (Clov_rand_coinc[0] < std::abs(time_diff) && std::abs(time_diff) < Clov_rand_coinc[1]) {
            Clov_gg_trand->Fill(E_Clover_original[j], E_Clover_original[k]);
            bool inRange1_j = (E_Clover_original[j] >= E_low1 && E_Clover_original[j] <= E_high1);
            bool inRange2_j = (E_Clover_original[j] >= E_low2 && E_Clover_original[j] <= E_high2);
            bool inRange1_k = (E_Clover_original[k] >= E_low1 && E_Clover_original[k] <= E_high1);
            bool inRange2_k = (E_Clover_original[k] >= E_low2 && E_Clover_original[k] <= E_high2);
            bool cross_event = (inRange1_j && inRange2_k) || (inRange2_j && inRange1_k);

            if (cross_event) {
              // ======================
        // Calcular cos(θ)
        // ======================
        if (crystal_pos.count(j) && crystal_pos.count(k)) {
            TVector3 vj = crystal_pos[j];
            TVector3 vk = crystal_pos[k];

            double cos_theta = vj.Dot(vk) / (vj.Mag() * vk.Mag());
            cos_theta = std::max(-1.0, std::min(1.0, cos_theta)); // evitar errores numéricos

            // ======================
            // Determinar bin
            // ======================
            int bin_index = int((cos_theta - cos_min) / cos_bin_width);
            if (bin_index >= 0 && bin_index < n_cos_bins) {
                h_cosbins_bckg[bin_index]->Fill(E_Clover_original[j], E_Clover_original[k]);
            }
        }
            }
          }
        } // end loop k (true/random)
      } // end if M_Clover_original > 1

      // MIXED EVENTS: emparejar gamma j del evento original con gamma k del evento i+Δ
      // Para mixed no requerimos M_Clover_original>1 ni M_Clover_next>1; basta que exista energía en cada lado.
      if (gotClov_next && M_Clover_next > 0 && tref_ok_mixed) {
        for (Int_t k = 0; k < num_Clov_crys; ++k) {
          if (E_Clover_next[k] <= EClov_Min) continue;
          if (j == k) continue; // si quieres permitir j==k entre eventos distintos, quita esta línea

          bool inRange1_j = (E_Clover_original[j] >= E_low1 && E_Clover_original[j] <= E_high1);
          bool inRange2_j = (E_Clover_original[j] >= E_low2 && E_Clover_original[j] <= E_high2);
          bool inRange1_k = (E_Clover_next[k] >= E_low1 && E_Clover_next[k] <= E_high1);
          bool inRange2_k = (E_Clover_next[k] >= E_low2 && E_Clover_next[k] <= E_high2);
          bool cross_event = (inRange1_j && inRange2_k) || (inRange2_j && inRange1_k);

          if (cross_event) {

            // ======================
        // Calcular cos(θ)
        // ======================
        if (crystal_pos.count(j) && crystal_pos.count(k)) {
            TVector3 vj = crystal_pos[j];
            TVector3 vk = crystal_pos[k];

            double cos_theta = vj.Dot(vk) / (vj.Mag() * vk.Mag());
            cos_theta = std::max(-1.0, std::min(1.0, cos_theta)); // evitar errores numéricos

            // ======================
            // Determinar bin
            // ======================
            int bin_index = int((cos_theta - cos_min) / cos_bin_width);
            if (bin_index >= 0 && bin_index < n_cos_bins) {
                h_cosbins_mixed[bin_index]->Fill(E_Clover_original[j], E_Clover_next[k]);
            }
        }
          }
        } // end loop k (mixed)
      } // end if gotClov_next && M_Clover_next > 0

    } // end loop j
  } // TIME REF original
  } // end if gotClov_original && M_Clover_original > 0

  // Progress bar (ajustado para no pasarse)
  progress = ((float)i / (float)(nentries - event_separation)) * 100.0f;
  if (progress >= progress_buffer || i == (nentries - event_separation - 1)) {
    PrintBar(progress);
    std::cout.flush();
    progress_buffer = progress + 1.0;
  }

} // end main loop

if(gotClov){
    Clov_gg_Notrand->Add(Clov_gg,1);
    //if(gotBeta) Clov_gg_Notrand_betagated->Add(Clov_gg_betagated,1);
    Clov_gg->Add(Clov_gg_trand,-1.*Clov_t_coinc/(Clov_rand_coinc[1] - Clov_rand_coinc[0]));
    //if(gotBeta) Clov_gg_betagated->Add(Clov_gg_trand_betagated,-1.*Clov_t_coinc/(Clov_rand_coinc[1] - Clov_rand_coinc[0]));
  }

int nHists = list_crys->GetSize();

for (int i = 0; i < nHists; ++i) {
    // Obtenemos los histogramas correspondientes de las dos listas
    TH2F* h_true = (TH2F*)list_crys->At(i);
    TH2F* h_bckg = (TH2F*)list_crys_bckg->At(i);

    if (!h_true || !h_bckg) continue;

    // Creamos una copia del histograma de coincidencias verdaderas
    TH2F* h_corr = (TH2F*)h_true->Clone(Form("%s", h_true->GetName()));
    h_corr->SetTitle(Form("%s (true - bckg scaled)", h_true->GetTitle()));

    // Restamos el fondo escalado por la ventana temporal
    h_corr->Add(h_bckg, -Clov_t_coinc / (Clov_rand_coinc[1] - Clov_rand_coinc[0]));

    // Lo añadimos a la lista final
    list_crys_true->Add(h_corr);
}




  progress = progress + 1.0;
  PrintBar(progress);
  std::cout.flush();
  std::cout << std::endl << std::endl << " Sorting finished!" << std::endl;
  std::cout << " Saving the histograms...";

// ==================
// WRITING HISTOGRAMS
// ==================

  //We create the output file and write the histograms
  TFile *fouty = new TFile(outputFileName.c_str(),"RECREATE");
    //list_others->Write();
  if(gotClov) list_Clov->Write(); //We only write the histos if we had clovers
  if(gotClov) list_crys->Write(); // Angular correlations
  if(gotClov) list_crys_bckg->Write(); // Angular correlations background
  if(gotClov) list_crys_true->Write(); // Angular correlations true = list_crys - list_crys_bckg
  if(gotClov) list_crys_mixed->Write(); // Angular Correlations mixed 
  //if(gotBeta && gotClov) list_crys_beta_gated->Write(); // Angular correlations beta gated
  //if(gotBeta && gotClov) list_crys_bckg_beta_gated->Write(); // Angular correlations background beta gated
  //if(gotBeta && gotClov) list_crys_true_beta_gated->Write(); // Angular correlations true beta gated
  //if(gotBeta && gotClov) list_crys_mixed_beta_gated->Write(); //Angular correlations mixed beta gated
/*
  if (list_Time->GetSize() > 0){
    TDirectory *TimeDir_dir = fouty->mkdir("TimeAlignment");
    TimeDir_dir->cd();
    list_Time->Write();                           // We only write the histos if we had All
    gDirectory->cd("/");
  }
*/
  fouty->Close();
	
  std::cout << " done!" << std::endl << std::endl;
  
  printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

  return 0;
}
