//#include "Mixed_Areas_Co60_prueba_305.C"
//#include "True_Areas_Co60_prueba_305_binned.C"
//#include "Mixed_Areas_Co60_prueba_302_binned7.C"
//#include "True_Areas_Co60_prueba_302_binned7.C"
//#include "True_Areas_Co60_prueba_305_binned_15.C"
//#include "Mixed_Areas_Co60_prueba_305_binned_15.C"
//#include "True_Areas_Co60_prueba_305_binned7.1_10.C"
//#include "Mixed_Areas_Co60_prueba_305_binned7.1_10.C"

//#include "True_Areas_Co60_prueba_305_binned_0.1.C"
//#include "Mixed_Areas_Co60_prueba_305_binned_0.1.C"


//Co60 good data
#include "True_Areas_Co60_prueba_305_binned10.C"
#include "Mixed_Areas_Co60_prueba_305_binned10.C"

/*
#include "True_Areas_Co60_prueba_202_206_binned10_0.5.C"
#include "Mixed_Areas_Co60_prueba_202_206_binned10_0.5.C"
*/
/*
#include "True_Areas_Co60_prueba_305_new_positions.C"
#include "Mixed_Areas_Co60_prueba_305_new_positions.C"
*/
/*
#include "True_Areas_Rb100_prueba_202_206_binned10.C"
#include "Mixed_Areas_Rb100_prueba_202_206_binned10.C"*/
/*
// Rb100 808_129 good data
#include "True_Areas_Rb100_prueba_202_206_binned10_808.C"
#include "Mixed_Areas_Rb100_prueba_202_206_binned10_808.C"
*/
/*
#include "True_Areas_Rb100_prueba_202_206_binned10_287_129.C"
#include "Mixed_Areas_Rb100_prueba_202_206_binned10_287_129.C"
*/
//#include "True_Areas_Rb100_prueba_202.C"
//#include "Mixed_Areas_Rb100_prueba_202.C"
//#include "Mixed_Areas_Rb100_prueba_202_206_1201.C"
//#include "True_Areas_Rb100_prueba_202_206_1201.C"
//#include "True_Areas_Co60_prueba_305_binned9_0.25.C"
//#include "Mixed_Areas_Co60_prueba_305_binned9_0.25.C"
void Test_Ang_Corr_60Co()
{
//std::vector<double> angles_degree ={7.5, 22.5, 37.5, 52.5, 67.5, 82.5, 97.5, 112.5, 127.5, 142.5, 157.5, 172.5};
//std::vector<double> angles_degree ={17, 83, 98, 170};

//std::vector<double> angles_degree = {13.65,19.53,      29.38	    ,      35.05	   ,       44.39	  ,      52.34	   ,      60.20	  ,      68.47	  ,       75.16	 ,      84.12	  ,      91.06	   ,      99.84	   ,     107.84	   ,     116.83	  ,     124.22	   ,     133.00	  ,     139.23	   ,      148.72	  ,     156.67	   ,     164.31	   ,    170.25	   ,      179.69	  };
//std::vector<double> angles_degree = {13.65,17.79,20.40,28.87,32.53,38.68,43.09,47.16, 52.26,57.48,61.55,67.97,73.08,77.03,82.35,88.14,91.73,97.91,102.63,107.35,112.25,118.04,122.7,127.6,132.77,137.43,143.02,147.71,151.58,158.54,162.22,165.88,170.33,179.69};
//td::vector<double> angles_degree ={13.65, 22.78, 39.26, 52.15, 67.97, 83.67, 96.91, 112.29, 127.52, 141.01, 156.00, 170.35};
//std::vector<double> angles_degree = {14.3,24.16,33.90,44.77,54.85,64.48,74.89,85.79,94.85,105.33,115.63,124.91,135.56,146.03,154.18,164.44,177.81};
//std::vector<double> pairs = {104, 60, 50, 60, 170, 94, 150, 92, 206, 210, 270, 194, 194, 152, 130, 146, 84, 114, 60, 80, 8, 24};
//std::vector<double> pairs = {36,32,32,32};
//std::vector<double> angles_degree ={10, 30, 50, 70, 90, 110, 130, 150, 170};
std::vector<double> ratios, ratio_errors;
std::vector<double> cos_angles = {-0.96,-0.86,-0.74,-0.65,-0.55,-0.46,-0.34,-0.26,-0.15,-0.04,0.04,0.15,0.26,0.35,0.47,0.54,0.65,0.73,0.85,0.96};
//std::vector<double> cos_angles = {-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9};

   // std::vector<double> cos_angles = {-0.93, -0.64, -0.36, -0.1, 0.1,0.37,0.64,0.92 };

    for (size_t i = 0; i < areas_True.size(); ++i) {
        double T = areas_True[i];
        double M = areas_Mixed[i];
        //double sigma_T = area_errors_True[i];
        //double sigma_M = area_errors_Mixed[i];
        double sigma_T = sqrt(T);
        double sigma_M = sqrt(M);

        double ratio =  /*1.58**/T / M;
        double error = /*1.58**/std::sqrt(std::pow(sigma_T / M, 2) + std::pow(sigma_M * T / (M*M), 2));

        ratios.push_back(ratio);
        ratio_errors.push_back(error);
    }

    TCanvas* c = new TCanvas("c", "Normalized Counts vs cos(theta)", 800, 600);
    TGraphErrors* gr = new TGraphErrors(cos_angles.size(), &cos_angles[0], &ratios[0], 0, &ratio_errors[0]);
    gr->SetTitle("Normalized Counts vs cos(#theta) ;cos(#theta);Normalized Counts");
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->Draw("AP");

    TF1* f_legendre = new TF1("f_legendre",
        "[0]*(1 + [1]*(1.5*x*x - 0.5) + [2]*(4.375*x*x*x*x - 3.75*x*x + 0.375))",
        -1.0, 1.0);
/*
    TF1* theo_legendre = new TF1("theo_legendre",
        "[0]*(1 + 0.102041*[1]*(1.5*x*x - 0.5) + 0.0090703*[2]*(4.375*x*x*x*x - 3.75*x*x + 0.375))",
        -1.0, 1.0);
*/

// 4->2->0 pura
TF1* theo_legendre = new TF1("theo_legendre",
        "[0]*(1 + 0.102041*0.94*(1.5*x*x - 0.5) + 0.0090703*0.83*(4.375*x*x*x*x - 3.75*x*x + 0.375))",
        -1.0, 1.0);
/*
// 0->2->0 pura
TF1* theo_legendre = new TF1("theo_legendre",
        "[0]*(1 + 0.35714*0.94*(1.5*x*x - 0.5) + 1.142857*0.84*(4.375*x*x*x*x - 3.75*x*x + 0.375))",
        -1.0, 1.0);
*/
    f_legendre->SetParameter(0, *std::max_element(ratios.begin(), ratios.end()));
    f_legendre->SetParameter(1, 0.0);
    f_legendre->SetParameter(2, 0.0);

    theo_legendre->SetParameter(0, *std::max_element(ratios.begin(), ratios.end()));
    theo_legendre->SetParameter(1, 0.0);
    theo_legendre->SetParameter(2, 0.0);

    gr->Fit(f_legendre, "RQ"); 
    gr->Fit(theo_legendre, "RQ");

    f_legendre->SetLineColor(kRed);
    f_legendre->SetLineWidth(2);
    f_legendre->Draw("same");
    theo_legendre->SetLineColor(kGreen+2);
    theo_legendre->SetLineWidth(2);
    theo_legendre->Draw("same");

    TLegend *legend = new TLegend(0.5,0.5,0.75,0.75);
    legend->AddEntry(gr, "Experimental Data");
    legend->AddEntry(f_legendre, "Experimental Fit");
    legend->AddEntry(theo_legendre, "Theoretical curve 0^{+}->2^{+}->0{+}");
    legend->SetBorderSize(0); 
    legend->Draw(); 
    

    std::cout << "\n==== Resultado del ajuste con Legendre ====\n";
    std::cout << "A00: " << f_legendre->GetParameter(0) << " ± " << f_legendre->GetParError(0) << std::endl;
    std::cout << "a2:  " << f_legendre->GetParameter(1) << " ± " << f_legendre->GetParError(1) << std::endl;
    std::cout << "a4:  " << f_legendre->GetParameter(2) << " ± " << f_legendre->GetParError(2) << std::endl;

    std::cout << "\n==== Resultado para obtener los parámetros de absorción ====\n";
    std::cout << "A00: " << theo_legendre->GetParameter(0) << " ± " << theo_legendre->GetParError(0) << std::endl;
    std::cout << "gamma:  " << theo_legendre->GetParameter(1) << " ± " << theo_legendre->GetParError(1) << std::endl;
    std::cout << "beta:  " << theo_legendre->GetParameter(2) << " ± " << theo_legendre->GetParError(2) << std::endl;

    c->Update();

    // Chi squared experimental vs fit
    double chi2_exp_fit = 0.0;
    int n_points = cos_angles.size();
    int n_params = 3;

    for (size_t i = 0; i < n_points; ++i) {
        double xi = cos_angles[i];
        double yi_obs = ratios[i];
        double sigma_i = ratio_errors[i];
        double yi_fit = f_legendre->Eval(xi);
        if (sigma_i > 0) chi2_exp_fit += std::pow((yi_obs - yi_fit) / sigma_i, 2);
    }
    double chi2ndf_exp_fit = chi2_exp_fit / (n_points - n_params);

    std::cout << "\nχ² (exp vs f_legendre): " << chi2_exp_fit
              << "   χ²/NDF = " << chi2ndf_exp_fit << std::endl;

TLegend* legend_chi = new TLegend(0.25, 0.65, 0.45, 0.85);
legend_chi->SetBorderSize(0);
legend_chi->SetTextSize(0.035);
legend_chi->AddEntry((TObject*)0, Form("#chi^{2}_{exp/fit} / NDF = %.3f", chi2ndf_exp_fit), "");


legend_chi->Draw();

    // Chi squared experimental vs theo
    int n_params_theo = 1;
    double chi2_exp_theo = 0.0;
    for (size_t i = 0; i < n_points; ++i) {
        double xi = cos_angles[i];
        double yi_obs = ratios[i];
        double sigma_i = ratio_errors[i];
        double yi_theo = theo_legendre->Eval(xi);
        if (sigma_i > 0) chi2_exp_theo += std::pow((yi_obs - yi_theo) / sigma_i, 2);
    }
    double chi2ndf_exp_theo = chi2_exp_theo / (n_points - n_params_theo);

    std::cout << "χ² (exp vs theo_legendre): " << chi2_exp_theo
              << "   χ²/NDF = " << chi2ndf_exp_theo << std::endl;

    // Chi squared fit vs theo
    double chi2_fit_theo = 0.0;
    int n_params_compare = 0;
    for (size_t i = 0; i < n_points; ++i) {
        double xi = cos_angles[i];
        double yi_fit = f_legendre->Eval(xi);
        double yi_theo = theo_legendre->Eval(xi);
        double sigma_i = ratio_errors[i];
        if (sigma_i > 0) chi2_fit_theo += std::pow((yi_fit - yi_theo) / sigma_i, 2);
    }
    double chi2ndf_fit_theo = chi2_fit_theo / (n_points - n_params_compare);

    std::cout << "χ² (f_legendre vs theo_legendre): " << chi2_fit_theo
              << "   χ²/NDF = " << chi2ndf_fit_theo << std::endl;

    std::cout << "\n==========================================" << std::endl;
}
