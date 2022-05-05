#include <memory>
#include <vector>
#include <iostream>

#include "TAxis.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"

inline std::shared_ptr<TGraph> GetTGraphFromFile(const char* path, const char* name, const bool verbose) {
    TFile* infile {TFile::Open(path, "READ")};
    if (infile->IsZombie()) {
        std::cerr << "\n\nError opening input file [" << path << "]" << std::endl;
        exit(100);
    }

    std::shared_ptr<TGraph> gr = std::shared_ptr<TGraph>(static_cast<TGraph*>(infile->Get(name)));

    if (verbose)
        std::cout << "\nFound TGraph in input file [" << gr->GetName() << "] --> [" << path << "]";

    return gr;
}

inline std::shared_ptr<TGraphErrors> GetTGraphErrorsFromFile(const char* path, const char* name, const bool verbose) {
    TFile* infile {TFile::Open(path, "READ")};
    if (infile->IsZombie()) {
        std::cerr << "\n\nError opening input file [" << path << "]" << std::endl;
        exit(100);
    }

    std::shared_ptr<TGraphErrors> gr = std::shared_ptr<TGraphErrors>(static_cast<TGraphErrors*>(infile->Get(name)));

    if (verbose)
        std::cout << "\nFound TGraphErrors in input file [" << gr->GetName() << "] --> [" << path << "]";

    return gr;
}

void fluxRatio(
    const char* reference, 
    const char* sample,
    const char* output_file, 
    const bool verbose) {

    auto reference_gr = GetTGraphFromFile(reference, "convertedCSV", verbose);
    auto sample_gr    = GetTGraphFromFile(sample, "gr_flux_E3", verbose);

    std::vector<double> energy                      (sample_gr->GetN(), 0);
    std::vector<double> simple_ratio                (sample_gr->GetN(), 0);
    std::vector<double> simple_ratio_percentage     (sample_gr->GetN(), 0);
    std::vector<double> ratio                       (sample_gr->GetN(), 0);
    std::vector<double> ratio_percentage            (sample_gr->GetN(), 0);

    for (int idx=0; idx<sample_gr->GetN(); ++idx) {
        energy[idx]                     = sample_gr->GetPointX(idx);
        simple_ratio[idx]               = sample_gr->GetPointY(idx) / reference_gr->Eval(energy[idx]);
        simple_ratio_percentage[idx]    = simple_ratio[idx] * 100;
        ratio[idx]                      = fabs(sample_gr->GetPointY(idx) - reference_gr->Eval(energy[idx]))/reference_gr->Eval(energy[idx]);
        ratio_percentage[idx]           = ratio[idx] * 100;
    }

    TGraph gr_ratio(sample_gr->GetN(), &energy[0], &ratio[0]);
    TGraph gr_simple_ratio(sample_gr->GetN(), &energy[0], &simple_ratio[0]);
    TGraph gr_ratio_percentage(sample_gr->GetN(), &energy[0], &ratio_percentage[0]);
    TGraph gr_simple_ratio_percentage(sample_gr->GetN(), &energy[0], &simple_ratio_percentage[0]);

    gr_simple_ratio.SetName("gr_simple_ratio");  
    gr_simple_ratio.SetTitle("Ratio");
    gr_simple_ratio.GetXaxis()->SetTitle("Energy [GeV]");
    gr_simple_ratio.GetYaxis()->SetTitle("sample / rerefernce");

    gr_ratio.SetName("gr_ratio");  
    gr_ratio.SetTitle("Ratio");
    gr_ratio.GetXaxis()->SetTitle("Energy [GeV]");
    gr_ratio.GetYaxis()->SetTitle("(sample-reference) / rerefernce");

    gr_ratio_percentage.SetName("gr_ratio_percentage");  
    gr_ratio_percentage.SetTitle("Ratio");
    gr_ratio_percentage.GetXaxis()->SetTitle("Energy [GeV]");
    gr_ratio_percentage.GetYaxis()->SetTitle("(sample-reference) / rerefernce [%]");

    gr_simple_ratio_percentage.SetName("gr_simple_ratio_percentage");  
    gr_simple_ratio_percentage.SetTitle("Ratio");
    gr_simple_ratio_percentage.GetXaxis()->SetTitle("Energy [GeV]");
    gr_simple_ratio_percentage.GetYaxis()->SetTitle("sample / rerefernce [%]");  

    TFile* outfile {TFile::Open(output_file, "RECREATE")};
    if (outfile->IsZombie()) {
        std::cerr << "\n\nError writing output ROOT file [" << output_file << "]" << std::endl;
        exit(100);
    }

    gr_ratio.Write();
    gr_simple_ratio.Write();
    gr_ratio_percentage.Write();
    gr_simple_ratio_percentage.Write();

    outfile->Close();

}