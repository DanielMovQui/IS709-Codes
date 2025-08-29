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

Double_t PeakFunction(Double_t* x, Double_t* par)
{
   Double_t xx = x[0];
   Double_t height = par[0];   
   Double_t centroid = par[1];  
   Double_t sigma = par[2];     
   Double_t beta = par[3];      
   Double_t relative = par[4];  

   Double_t gauss = height * (1.0 - relative / 100.0) * TMath::Gaus(xx, centroid, sigma);

   if(beta == 0.0) {
       return gauss;
   }

   return gauss + relative * height / 100.0 * TMath::Exp((xx - centroid) / beta) *
          TMath::Erfc( ((xx - centroid) / (TMath::Sqrt(2.0) * sigma)) + 
                       sigma / (TMath::Sqrt(2.0) * beta) );
}



void Fitting_Angular_Correlations() {
    // Path to the file made by histo_angular_corr.C 
   // const char* filename = "/home/daniel/Desktop/IS622/98Rb/Run421_histograms.root";
    const char* filename = "/data1/daniel/TRIUMF/S1723_2022/S1723-Aug2022/Calibrations/Co60/AngularCorrelation21643_000.root";

    TFile* file = TFile::Open(filename, "READ");

    // Gates and backgrounds for the proyection
    double gate_low = 1170, gate_high = 1176;
    double bg1_low = 1140, bg1_high = 1165;
    double bg2_low = 1180, bg2_high = 1200;

    // Areas and errors vectors
    std::vector<double> areas;
    std::vector<double> area_errors;

    for (int i = 0; i <1; ++i) { // if wanted, change de i variable of the loop to make the analysis for multiple histograms
        // Data selection, only one type, or true or mixed
        // True coincidences name
        // std::string histname = "Angular_correlation_hist_true_" + std::to_string(i);
        // Mixing data name
        std::string histname = "AngularCorrelationMixed" + std::to_string(i);
        TH2* h = (TH2*)file->Get(histname.c_str());

        // ===============
        // PROJECTION IN X 
        // ===============
        int binYmin = 1;
        int binYmax = h->GetYaxis()->GetNbins(); 
        TH1* h_projX = h->ProjectionX(("h_projX_" + std::to_string(i)).c_str(), binYmin, binYmax);
        h_projX->GetXaxis()->SetRangeUser(0, 3000);

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

        // Parameters of TRWPeak Function

        double Emin = 1330;
        double Emax = 1336;
        double height_init = h_signal->GetMaximum();
        double centroid_init = 1334;
        double sigma_init = 1;       // standard deviation
        double beta_init = 0.01;     // Inicial asymmetry
        double relative_init = 1;    // Inicial relative height 

        // Fitting and plot
        TF1* trwPeakFit = new TF1(("trwPeakFit_" + std::to_string(i)).c_str(), PeakFunction, Emin, Emax, 5);
        trwPeakFit->SetParameters(height_init, centroid_init, sigma_init, beta_init, relative_init);
        h_signal->Fit(trwPeakFit, "RQ");
        trwPeakFit->SetLineColor(kRed);
        trwPeakFit->SetLineWidth(2);
        trwPeakFit->Draw("same");
        c2->Update();

        // ========
        // RESULTS 
        // ========       

        std::cout << "===== Results of fit TRWPeak for " << histname << " =====" << std::endl;
        std::cout << "Amplitude (p0): " << trwPeakFit->GetParameter(0) << " ± " << trwPeakFit->GetParError(0) << std::endl;
        std::cout << "Centroid   (p1): " << trwPeakFit->GetParameter(1) << " ± " << trwPeakFit->GetParError(1) << std::endl;
        std::cout << "Sigma      (p2): " << trwPeakFit->GetParameter(2) << " ± " << trwPeakFit->GetParError(2) << std::endl;
        std::cout << "Lambda     (p3): " << trwPeakFit->GetParameter(3) << " ± " << trwPeakFit->GetParError(3) << std::endl;
        std::cout << "Relative   (p4): " << trwPeakFit->GetParameter(4) << " ± " << trwPeakFit->GetParError(4) << std::endl;
        std::cout << "Chi2/NDF: " << trwPeakFit->GetChisquare() << " / " << trwPeakFit->GetNDF() << std::endl;

        double area = trwPeakFit->Integral(Emin, Emax);
        double area_error = trwPeakFit->IntegralError(Emin, Emax);
        std::cout << "Area: " << area << " ± " << area_error << std::endl;

        areas.push_back(area);
        area_errors.push_back(area_error);
        
    }


// ================
// SAVING AREA DATA 
// ================


// Saving the data in a .C filled with areas an uncertainties, if wanted, modify to save true data if the code ruuns with true data. 
    std::ofstream fout("Mixed_Areas_Co60_prueba.C");
    fout << "std::vector<double> areas_Mixed = {";
    for (size_t j = 0; j < areas.size(); ++j) {
        fout << areas[j];
        if (j < areas.size() - 1) fout << ", ";
    }
    fout << "};\n";
    fout << "std::vector<double> area_errors_Mixed = {";
    for (size_t j = 0; j < area_errors.size(); ++j) {
        fout << area_errors[j];
        if (j < area_errors.size() - 1) fout << ", ";
    }
    fout << "};\n";
    fout.close();

    std::cout << "Results saved in Mixed_Areas_Co60_prueba.C" << std::endl;

}

