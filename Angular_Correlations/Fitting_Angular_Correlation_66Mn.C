
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
#include <TLine.h>
#include <TApplication.h> 


// ================
// FITTING FUNCTION
// ================

// Taken from GRSISORT TRWFITPEAK
/*
Double_t PeakFunction(Double_t* x, Double_t* par)
{
   Double_t xx = x[0];
   Double_t height = par[0];   
   Double_t centroid = par[1];  
   Double_t sigma = par[2];     
   Double_t beta = par[3];      
   Double_t relative = par[4];  

   Double_t gauss = height * (1.0 - relative / 100.0) * TMath::Gaus(xx, centroid, sigma);

   return gauss + relative * height / 100.0 * TMath::Exp((xx - centroid) / beta) *
          TMath::Erfc( ((xx - centroid) / (TMath::Sqrt(2.0) * sigma)) + 
                       sigma / (TMath::Sqrt(2.0) * beta) );
}*/

Double_t PeakFunction(Double_t* x, Double_t* par)
{
   Double_t xx = x[0];

   Double_t height     = par[0];
   Double_t centroid   = par[1];
   Double_t sigma      = par[2];

   Double_t beta_l     = par[3];
   Double_t rel_l      = par[4];

   Double_t beta_r     = par[5];
   Double_t rel_r      = par[6];

   Double_t gauss = height *
       (1.0 - rel_l/100.0 - rel_r/100.0) *
       TMath::Gaus(xx, centroid, sigma, true);

   Double_t left_tail =
       rel_l * height / 100.0 *
       TMath::Exp((xx - centroid)/beta_l) *
       TMath::Erfc((xx - centroid)/(TMath::Sqrt(2)*sigma)
                    + sigma/(TMath::Sqrt(2)*beta_l));

   Double_t right_tail =
       rel_r * height / 100.0 *
       TMath::Exp(-(xx - centroid)/beta_r) *
       TMath::Erfc((centroid - xx)/(TMath::Sqrt(2)*sigma)
                    + sigma/(TMath::Sqrt(2)*beta_r));

   return gauss + left_tail + right_tail;
}



/*
Double_t PeakFunction(Double_t* x, Double_t* par)
{
   Double_t xx = x[0];
   Double_t height = par[0];   
   Double_t centroid = par[1];  
   Double_t sigma = par[2];     
   Double_t beta  = par[3];      
   Double_t relative = par[4];  
   Double_t beta2 = par[5];       // << nuevo
   Double_t relative2 = par[6];   // << nuevo

   Double_t gauss = height * (1.0 - (relative + relative2) / 100.0) * 
                    TMath::Gaus(xx, centroid, sigma);

   Double_t tail1 = relative * height / 100.0 * 
                    TMath::Exp((xx - centroid) / beta) *
                    TMath::Erfc(((xx - centroid) / (TMath::Sqrt(2.0) * sigma)) +
                                sigma / (TMath::Sqrt(2.0) * beta));

   Double_t tail2 = relative2 * height / 100.0 * 
                    TMath::Exp((xx - centroid) / beta2) *
                    TMath::Erfc(((xx - centroid) / (TMath::Sqrt(2.0) * sigma)) +
                                sigma / (TMath::Sqrt(2.0) * beta2));

   return gauss + tail1 + tail2;
}
*/
/*
Double_t BackgroundFunction(Double_t* dim, Double_t* par)
{
   Double_t x        = dim[0];   // channel number used for fitting
   Double_t height   = par[0];   // height of photopeak
   Double_t centroid = par[1];   // Peak Centroid of non skew gaus
   Double_t sigma    = par[2];   // standard deviation of gaussian
   Double_t step     = par[5];   // Size of the step function;
   Double_t a = par[6];  // slope 
   Double_t b = par[7];  // offset
   Double_t step_func = TMath::Abs(step) * height / 100.0 * TMath::Erfc((x - centroid) / (TMath::Sqrt(2.0) * sigma));
   Double_t linear_bg = a*x + b;

   return step_func + linear_bg;
}*/

Double_t BackgroundFunction(Double_t* x, Double_t* par)
{
   Double_t xx       = x[0];

   Double_t step_amp = par[0];
   Double_t centroid = par[1];
   Double_t sigma    = par[2];

   Double_t slope    = par[3];
   Double_t offset   = par[4];

   Double_t step =
       step_amp *
       TMath::Erfc((xx - centroid)/(TMath::Sqrt(2)*sigma));

   Double_t linear = slope*xx + offset;

   return step + linear;
}

/*
Double_t TotalFitFunction(Double_t* x, Double_t* par) {

    Double_t peakPar[5] = {par[0], par[1], par[2], par[3], par[4]};
    Double_t bgPar[8] = {par[0], par[1], par[2], 0, 0, par[5], par[6], par[7]}; // reusar p0,p1,p2 y step para background
    Double_t peak = PeakFunction(x, peakPar);
    Double_t background = BackgroundFunction(x, bgPar);

    return peak + background;
}*/
/*
Double_t TotalFitFunction(Double_t* x, Double_t* par) {

    Double_t peakPar[10] = {par[0], par[1], par[2], par[3], par[4], 0, 0, 0, par[8], par[9]};
    Double_t bgPar[8] = {par[0], par[1], par[2], 0, 0, par[5], par[6], par[7]};
    Double_t peak = PeakFunction(x, peakPar);
    Double_t background = BackgroundFunction(x, bgPar);

    return peak + background;
}*/

Double_t TotalFitFunction(Double_t* x, Double_t* par)
{
   Double_t peakPar[7] = {
      par[0], par[1], par[2],
      par[3], par[4],
      par[5], par[6]
   };

   Double_t bgPar[5] = {
      par[7], par[1], par[2],
      par[8], par[9]
   };

   return PeakFunction(x, peakPar) +
          BackgroundFunction(x, bgPar);
}


void Fitting_Angular_Correlation_66Mn() {

    bool useTrue = true;

   // const char* filename = "/data1/daniel/TRIUMF/S1723_2022/S1723-Aug2022/66Mn/Data/Angular_Correlation_Histograms/AngularCorrelation21441-21533.root";
   // const char* filename = "/data1/daniel/TRIUMF/S1723_2022/S1723-Aug2022/66Mn/Data/Angular_Correlation_Histograms/0.25_keV_binning_50_depth/ROOT/AngularCorrelation_0.25_full.root";
   // const char* filename = "/data1/daniel/TRIUMF/S1723_2022/S1723-Aug2022/66Mn/Data/Angular_Correlation_Histograms/0.25_keV_binning_150_depth/ROOT/AngularCorrelation_0.25_150_full.root"; 
   //const char* filename = "/data1/daniel/TRIUMF/S1723_2022/S1723-Aug2022/66Mn/Data/Angular_Correlation_Histograms/0.25_keV_binning_150_depth/110cm/ROOT/AngularCorrelation_66Mn_full.root";
     // const char* filename = "/data1/daniel/TRIUMF/S1723_2022/S1723-Aug2022/66Mn/Data/Angular_Correlation_Histograms/0.5_keV_binning_150_depth_full/ROOT/AngularCorrelation_66Mn_full.root";
     const char* filename = "/data1/daniel/TRIUMF/S1723_2022/S1723-Aug2022/66Mn/Data/Angular_Correlation_Histograms/1.0_keV_binning_3kabove_150_depth_full/ROOT/AngularCorrelation_66Mn_3kabove_full.root";
     //  const char* filename = "/data1/daniel/TRIUMF/S1723_2022/S1723-Aug2022/66Mn/Data/Angular_Correlation_Histograms/0.5_keV_binning_150_depth_full_time_corrected/ROOT/AngularCorrelation_66Mn_def_corrected.root";
    // const double binning_E = 0.25;
   const double binning_E = 1.0;

    TFile* file = TFile::Open(filename, "READ");

    // Definir los rangos de gate y fondos
    double gate_low = 572, gate_high = 576;
    double bg1_low = 545, bg1_high = 565;
    double bg2_low = 578, bg2_high = 595;


/*
    // Definir los rangos de gate y fondos
    double gate_low = 837, gate_high = 842.5;
    double bg1_low = 775, bg1_high = 790;
    double bg2_low = 925, bg2_high = 955;
*/
/*
    double gate_low = 838, gate_high = 842.5;
    double bg1_low = 720, bg1_high = 730;
    double bg2_low = 850, bg2_high = 870;
*/
/*
    double gate_low = 831.5, gate_high = 836.5;
    double bg1_low = 789, bg1_high = 797;
    double bg2_low = 845, bg2_high = 855;
*/
/*
    double gate_low = 1305, gate_high = 1310.5;
    double bg1_low = 1282, bg1_high = 1289;
    double bg2_low = 1360, bg2_high = 1367;
*/
/*
    double gate_low = 1545.5, gate_high = 1550;
    double bg1_low = 1531, bg1_high = 1539;
    double bg2_low = 1576, bg2_high = 1584;
*/
/*
    double gate_low = 2677.5, gate_high = 2684.5;
    double bg1_low = 2642, bg1_high = 2656;
    double bg2_low = 2782, bg2_high = 2798;
*/
/*
    double gate_low = 1775.5, gate_high = 1780.5;
    double bg1_low = 1755, bg1_high = 1766;
    double bg2_low = 1812, bg2_high = 1826;
*/
/*
    double gate_low = 2087, gate_high = 2093.5;
    double bg1_low = 2027.5, bg1_high = 2042;
    double bg2_low = 2155, bg2_high = 2164;
*/
/*
    double gate_low = 2708, gate_high = 2715;
    double bg1_low = 2689, bg1_high = 2703;
    double bg2_low = 2730, bg2_high = 2745;
*/
/*
    double gate_low = 2127.5, gate_high = 2136.5;
    double bg1_low = 2094, bg1_high = 2107;
    double bg2_low = 2207, bg2_high = 2226;
*/
/*
    double gate_low = 1061, gate_high = 1066.5;
    double bg1_low = 1052.5, bg1_high = 1059.5;
    double bg2_low = 1067.5, bg2_high = 1074;
*/
/*
    double gate_low = 2462.5, gate_high = 2467.5;
    double bg1_low = 2449, bg1_high = 2460;
    double bg2_low = 2470, bg2_high = 2479;
    */
   /*
    double gate_low = 1276.5, gate_high = 1281.5;
    double bg1_low = 1256, bg1_high = 1262;
    double bg2_low = 1283, bg2_high = 1289;
    */
/*
    double gate_low = 3277, gate_high = 3294;
    double bg1_low = 3300, bg1_high = 3328;
    double bg2_low = 3240, bg2_high = 3265;
*/
/*
    double gate_low = 3070, gate_high = 3080;
    double bg1_low = 3020, bg1_high = 3035;
    double bg2_low = 3085, bg2_high = 3100;
*/
/*
    double gate_low = 2297, gate_high = 2303;
    double bg1_low = 2271, bg1_high = 2291;
    double bg2_low = 2321, bg2_high = 2351;
*/
/*
    double gate_low = 2894, gate_high = 2902;
    double bg1_low = 2845, bg1_high = 2860;
    double bg2_low = 2910, bg2_high = 2925;
*/
/*
    double gate_low = 2894, gate_high = 2902;
    double bg1_low = 2845, bg1_high = 2860;
    double bg2_low = 2910, bg2_high = 2925;    
*/
    // Listas para guardar áreas y errores
    std::vector<double> areas;
    std::vector<double> area_errors;

    for (int i = 0; i <51; ++i) {

        std::string histname;
        if (useTrue) {
            histname = "AngularCorrelationTrueCoincidences" + std::to_string(i);
           // histname = "Angular_correlation_hist_" + std::to_string(i) + "True_binned_";

        } else {
            histname = "AngularCorrelationMixed" + std::to_string(i);
        }
        TH2* h = (TH2*)file->Get(histname.c_str());

        int binYmin = 1;
        int binYmax = h->GetYaxis()->GetNbins(); 
        TH1* h_projX = h->ProjectionX(("h_projX_" + std::to_string(i)).c_str(), binYmin, binYmax);
          
        // Plot projection
        TCanvas* c0 = new TCanvas(("c0_" + std::to_string(i)).c_str(), ("Projection X with gates and backgrounds " + std::to_string(i)).c_str(), 900, 600);
        h_projX->SetLineColor(kBlack);
        h_projX->SetTitle(("Projection X with gate and backgrounds [" + std::to_string(i) + "];Energy (keV);Counts").c_str());
        h_projX->Draw("hist");

        // Plot of the gates and background with lines
        TLine *l_gate_low = new TLine(gate_low, 0, gate_low, h_projX->GetMaximum());
        TLine *l_gate_high = new TLine(gate_high, 0, gate_high, h_projX->GetMaximum());
        TLine *l_bg1_low = new TLine(bg1_low, 0, bg1_low, h_projX->GetMaximum());
        TLine *l_bg1_high = new TLine(bg1_high, 0, bg1_high, h_projX->GetMaximum());
        TLine *l_bg2_low = new TLine(bg2_low, 0, bg2_low, h_projX->GetMaximum());
        TLine *l_bg2_high = new TLine(bg2_high, 0, bg2_high, h_projX->GetMaximum());

        l_gate_low->SetLineColor(kBlue); l_gate_low->SetLineStyle(3); l_gate_low->Draw();
        l_gate_high->SetLineColor(kBlue); l_gate_high->SetLineStyle(3); l_gate_high->Draw();
        l_bg1_low->SetLineColor(kGreen+2); l_bg1_low->SetLineStyle(3); l_bg1_low->Draw();
        l_bg1_high->SetLineColor(kGreen+2); l_bg1_high->SetLineStyle(3); l_bg1_high->Draw();
        l_bg2_low->SetLineColor(kRed); l_bg2_low->SetLineStyle(3); l_bg2_low->Draw();
        l_bg2_high->SetLineColor(kRed); l_bg2_high->SetLineStyle(3); l_bg2_high->Draw();

        c0->Update();

        int bin_gate_low = h->GetYaxis()->FindBin(gate_low);
        int bin_gate_high = h->GetYaxis()->FindBin(gate_high);
        int bin_bg1_low = h->GetYaxis()->FindBin(bg1_low);
        int bin_bg1_high = h->GetYaxis()->FindBin(bg1_high);
        int bin_bg2_low = h->GetYaxis()->FindBin(bg2_low);
        int bin_bg2_high = h->GetYaxis()->FindBin(bg2_high);
        
        // Gates background and signal
        TH1* h_gate = h->ProjectionX(("h_gate_" + std::to_string(i)).c_str(), bin_gate_low, bin_gate_high);
        TH1* h_bg1 = h->ProjectionX(("h_bg1_" + std::to_string(i)).c_str(), bin_bg1_low, bin_bg1_high);
        TH1* h_bg2 = h->ProjectionX(("h_bg2_" + std::to_string(i)).c_str(), bin_bg2_low, bin_bg2_high);

        double width_gate = gate_high - gate_low;
        double width_bg1 = bg1_high - bg1_low;
        double width_bg2 = bg2_high - bg2_low;
        double scale_bg1 = width_gate / width_bg1;
        double scale_bg2 = width_gate / width_bg2;

        h_bg1->Scale(scale_bg1);
        h_bg2->Scale(scale_bg2);
        
        // Average background calculation
        TH1* h_bg = (TH1*)h_bg1->Clone(("h_bg_" + std::to_string(i)).c_str());
        h_bg->Add(h_bg2);
        h_bg->Scale(0.5);

        TH1* h_signal = (TH1*)h_gate->Clone(("h_signal_" + std::to_string(i)).c_str());
        h_signal->Add(h_bg, -1);

        // ===============
        // FIT AND PLOT
        // ===============        

        // Plot final coincidences from signal gated before
        TCanvas* c2 = new TCanvas(("c2_" + std::to_string(i)).c_str(), ("Coincidences " + std::to_string(i)).c_str(), 900, 600);
        h_signal->SetLineColor(kBlack); h_signal->SetLineWidth(2);
        h_signal->SetTitle(("Coincidences [" + std::to_string(i) + "];Energy (keV);Counts").c_str());
        h_signal->Draw("HIST");
       // h_signal->GetXaxis()->SetRangeUser(520, 620);
/*
        double Emin = 815;
        double Emax = 870;
        double height_init = h_signal->GetMaximum();
        double centroid_init = 840.5;
        double sigma_init = 1;       // standard deviation
        double beta_init = 1;     // Inicial asymmetry
        double relative_init = 50;    // Inicial relative height 
        double step_init = 1;
        double a_init = 1;
        double b_init = 1;
*/
/*
        double Emin = 550;
        double Emax = 590;
        double height_init = h_signal->GetMaximum();
        double centroid_init = 573;
        double sigma_init = 0.9;       // standard deviation
        double beta_init = 1;     // Inicial asymmetry
        double relative_init = 50;    // Inicial relative height 
        double step_init = 0.1;
        double a_init = 1;
        double b_init = 1;
        double beta_r = 1;      // cola derecha
        double relative_r = 50;*/

        double Emin = 3420;
        double Emax = 3480;
        double height_init = h_signal->GetMaximum();
        double centroid_init = 3449.5;
        double sigma_init = 0.9;       // standard deviation
        double beta_init = 1;     // Inicial asymmetry
        double relative_init = 50;    // Inicial relative height 
        double step_init = 0.1;
        double a_init = 1;
        double b_init = 1;
        double beta_r = 1;      // cola derecha
        double relative_r = 50;

/*
        TF1* totalFit = new TF1("totalFit", TotalFitFunction, Emin, Emax, 10);
        totalFit->SetParameters(height_init, centroid_init, sigma_init, beta_init, relative_init, step_init, a_init, b_init, beta_r, relative_r);
    */
        /* TF1* totalFit = new TF1("totalFit", TotalFitFunction, Emin, Emax, 10);
totalFit->SetParameters(height_init, centroid_init, sigma_init,
                        beta_init, relative_init, step_init,
                        a_init, b_init,
                        beta_init_l, relative_init_l);  // beta2 y relative2 INICIALES*/
          TF1* totalFit = new TF1("totalFit", TotalFitFunction, Emin, Emax, 10);  
totalFit->SetParameters(
   height_init, centroid_init, sigma_init,
   1, 5,      // beta_l, rel_l
   1.0, 10,       // beta_r, rel_r
   1,           // step_amp
   0.0, 10        // slope, offset
);

totalFit->SetParLimits(2, 0.3, 5.0);    // sigma
totalFit->SetParLimits(4, 0.0, 30.0);   // rel_l
totalFit->SetParLimits(6, 0.0, 30.0);   // rel_r
totalFit->SetParLimits(3, 0.1, 10.0);   // beta_l
totalFit->SetParLimits(5, 0.1, 10.0);   // beta_r


        TFitResultPtr fitResult = h_signal->Fit(totalFit, "RQS");
        totalFit->SetLineColor(kRed);
        totalFit->SetNpx(10000);
        totalFit->SetLineWidth(2);
        totalFit->Draw("same");
        c2->Update();


        // ========
        // RESULTS 
        // ========       

        std::cout << "===== Results of fit TRWPeak for " << histname << " =====" << std::endl;
        std::cout << "Amplitude (p0): " << totalFit->GetParameter(0) << " ± " << totalFit->GetParError(0) << std::endl;
        std::cout << "Centroid   (p1): " << totalFit->GetParameter(1) << " ± " << totalFit->GetParError(1) << std::endl;
        std::cout << "Sigma      (p2): " << totalFit->GetParameter(2) << " ± " << totalFit->GetParError(2) << std::endl;
        std::cout << "beta_l     (p3): " << totalFit->GetParameter(3) << " ± " << totalFit->GetParError(3) << std::endl;
        std::cout << "rel_l      (p4): " << totalFit->GetParameter(4) << " ± " << totalFit->GetParError(4) << std::endl;
        std::cout << "beta_r     (p5): " << totalFit->GetParameter(5) << " ± " << totalFit->GetParError(5) << std::endl;
        std::cout << "rel_r      (p6): " << totalFit->GetParameter(6) << " ± " << totalFit->GetParError(6) << std::endl;
        std::cout << "step_amp   (p7): " << totalFit->GetParameter(7) << " ± " << totalFit->GetParError(7) << std::endl;
        std::cout << "slope      (p8): " << totalFit->GetParameter(8) << " ± " << totalFit->GetParError(8) << std::endl;
        std::cout << "offset     (p9): " << totalFit->GetParameter(9) << " ± " << totalFit->GetParError(9) << std::endl;
        std::cout << "Fit status = " << fitResult->Status() << std::endl;
        std::cout << "Cov status = " << fitResult->CovMatrixStatus() << std::endl;

        std::cout << "Chi2/NDF: " << totalFit->GetChisquare() << " / " << totalFit->GetNDF() << std::endl;
/*
        double area = totalFit->Integral(Emin, Emax)/binning_E;
        double area_error = totalFit->IntegralError(Emin, Emax)/binning_E;
        std::cout << "Area: " << area << " ± " << area_error << std::endl;

// Cálculo del área solo del fondo (lineal + step)
TF1* backgroundFit = new TF1("backgroundFit", BackgroundFunction, Emin, Emax, 8);
backgroundFit->SetParameters(
    totalFit->GetParameter(0),  // height
    totalFit->GetParameter(1),  // centroid
    totalFit->GetParameter(2),  // sigma
    0, 0,                       // par[3], par[4] (no se usan en fondo)
    totalFit->GetParameter(5),  // step
    totalFit->GetParameter(6),  // a
    totalFit->GetParameter(7)   // b
);

double area_bg = backgroundFit->Integral(Emin, Emax)/binning_E;
double area_bg_error = backgroundFit->IntegralError(Emin, Emax)/binning_E;

// Área neta del pico (restando el fondo)
double area_peak = area - area_bg;
double area_peak_error = std::sqrt(area_error * area_error + area_bg_error * area_bg_error);
*/

TF1* peakOnly = new TF1("peakOnly", PeakFunction, Emin, Emax, 7);
peakOnly->SetParameters(
   totalFit->GetParameter(0),
   totalFit->GetParameter(1),
   totalFit->GetParameter(2),
   totalFit->GetParameter(3),
   totalFit->GetParameter(4),
   totalFit->GetParameter(5),
   totalFit->GetParameter(6)
);
TMatrixDSym cov = fitResult->GetCovarianceMatrix();

TMatrixDSym cov7(7);
for(int i=0;i<7;i++)
  for(int j=0;j<7;j++)
    cov7(i,j) = cov(i,j);

double area_peak =
   peakOnly->Integral(Emin, Emax) / binning_E;

double area_peak_error =
   peakOnly->IntegralError(
      Emin, Emax,
      peakOnly->GetParameters(),
      cov7.GetMatrixArray()
   ) / binning_E;


std::cout << "\n=== ÁREA CORREGIDA (solo pico) ===" << std::endl;
//std::cout << "Área total (pico + fondo): " << area << std::endl;
//std::cout << "Área fondo: " << area_bg << std::endl;
std::cout << "Área neta (solo pico): " << area_peak << " ± " << area_peak_error << std::endl;

// Guardar el área corregida
areas.push_back(area_peak);
area_errors.push_back(area_peak_error);
        
    }

// ================
// SAVING AREA DATA 
// ================


// Saving the data in a .C filled with areas an uncertainties, if wanted, modify to save true data if the code ruuns with true data. 
std::string output_name;
    std::string prefix;
    if (useTrue) {
       // output_name = "True_Areas_Rb100_Sum_binned10_0.5_808.C";
       output_name = "/home/daniel/Desktop/Analisis_Fe/66Mn/Definitive_Areas/True_Areas_66Mn_110mm_1.0_150_3450_573.C";
        //output_name = "True_Areas_Co60_prueba_305_binned10_0.5.C";
       // output_name = "True_Areas_Rb98_Sum_binned10_0.5_1079.C";

        prefix = "True";
    } else {
       // output_name = "Mixed_Areas_Rb100_Sum_binned10_0.5_808.C";
       output_name = "/home/daniel/Desktop/Analisis_Fe/66Mn/Definitive_Areas/Mixed_Areas_66Mn_110mm_1.0_150_3450_573.C";
       // output_name = "Mixed_Areas_Co60_prueba_305_binned10_0.5.C";
      // output_name = "Mixed_Areas_Rb98_Sum_binned10_0.5_1079.C";
        prefix = "Mixed";
    }

    std::ofstream fout(output_name);
    fout << "std::vector<double> areas_" << prefix << " = {";
    for (size_t j = 0; j < areas.size(); ++j) {
        fout << areas[j];
        if (j < areas.size() - 1) fout << ", ";
    }
    fout << "};\n";

    fout << "std::vector<double> area_errors_" << prefix << " = {";
    for (size_t j = 0; j < area_errors.size(); ++j) {
        fout << area_errors[j];
        if (j < area_errors.size() - 1) fout << ", ";
    }
    fout << "};\n";
    fout.close();

}
