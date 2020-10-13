#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

void text_converter(const char* in_path, const char* out_path)
{
    std::ifstream input_file(in_path);
    if (!input_file.is_open())
    {
        std::cerr << "\nERROR 100! File not open " << in_path << "\n\n";
        exit(100);
    }
    std::string input_string((std::istreambuf_iterator<char>(input_file)), (std::istreambuf_iterator<char>()));
    input_file.close();
    std::string tmp_str;
    std::istringstream input_stream(input_string);
    std::string::size_type sz;
    std::vector<double> data_x;
    std::vector<double> data_y;
    while (input_stream >> tmp_str)
    {
        data_x.push_back(stod(tmp_str, &sz));
        input_stream >> tmp_str;
        data_y.push_back(stod(tmp_str, &sz));
    }
    
    TFile outFile(out_path, "RECREATE");
    if (outFile.IsZombie())
    {
        std::cerr << "\nERROR 100! Output file not created " << "\n\n";
        exit(100); 
    }

    TGraph gr(data_x.size(), &(data_x[0]), &(data_y[0]));
    gr.SetName("convertedCSV");
    gr.Write();
    outFile.Close();
}

void text_converter_full(const char* in_path, const char* out_path)
{
    std::ifstream input_file(in_path);
    if (!input_file.is_open())
    {
        std::cerr << "\nERROR 100! File not open " << in_path << "\n\n";
        exit(100);
    }
    std::string input_string((std::istreambuf_iterator<char>(input_file)), (std::istreambuf_iterator<char>()));
    input_file.close();
    std::string tmp_str;
    std::istringstream input_stream(input_string);
    std::string::size_type sz;
    std::vector<double> data_x;
    std::vector<double> data_y;
    std::vector<double> y_x_3;
    std::vector<double> x_min;
    std::vector<double> x_max;
    std::vector<double> x_range;
    std::vector<double> stat_err;
    std::vector<double> sys_err;
    std::vector<double> tot_err;
    std::vector<double> tot_err_3;
    int counter=0;
    while (input_stream >> tmp_str)
    {
        // Get X minimum
        x_min.push_back(stod(tmp_str, &sz));
        // Get X maximum
        input_stream >> tmp_str;
        x_max.push_back(stod(tmp_str, &sz));
        // Get X range
        x_range.push_back(x_max.at(counter)-x_min.at(counter));
        // Get X
        input_stream >> tmp_str;
        data_x.push_back(stod(tmp_str, &sz));
        // Get Y
        input_stream >> tmp_str;
        data_y.push_back(stod(tmp_str, &sz));
        // Get statistical error
        input_stream >> tmp_str;
        stat_err.push_back(stod(tmp_str, &sz));
        // Get systematic error
        input_stream >> tmp_str;
        sys_err.push_back(stod(tmp_str, &sz));
        // Get total error
        tot_err.push_back(stat_err.at(counter)+sys_err.at(counter));
        // Build y_x_3
        y_x_3.push_back(data_y.at(counter)*pow(data_x.at(counter), 3));
        tot_err_3.push_back(pow(data_x.at(counter), 3) * tot_err.at(counter));
        // Update counter
        ++counter;
    }

    TFile outFile(out_path, "RECREATE");
    if (outFile.IsZombie())
    {
        std::cerr << "\nERROR 100! Output file not created " << "\n\n";
        exit(100); 
    }

    TGraphErrors gr(data_x.size(), &(data_x[0]), &(data_y[0]), &(x_range[0]), &(tot_err[0]));
    TGraphErrors gr3(data_x.size(), &(data_x[0]), &(y_x_3[0]), &(x_range[0]), &(tot_err_3[0]));
    gr.SetName("gr_all_electron_flux");
    gr3.SetName("gr_all_electron_flux_E3");
    gr.Write();
    gr3.Write();
    outFile.Close();
}

void flux_ratio(const char* pg_analysis, const char* geneva_analysis, const char* result)
{
    TFile pg_file(pg_analysis, "READ");
    if (!pg_file.IsOpen())
    {
        std::cerr << "\n\nError reading input file: [" << pg_analysis << "]\n\n";
        exit(100);
    }
    auto pg_flux = static_cast<TGraphErrors*>(pg_file.Get("gr_all_electron_flux_E3"));
    pg_file.Close();
    TFile geneva_file(geneva_analysis, "READ");
    if (!geneva_file.IsOpen())
    {
        std::cerr << "\n\nError reading input file: [" << geneva_analysis << "]\n\n";
        exit(100);
    }
    auto geneva_flux = static_cast<TGraphErrors*>(geneva_file.Get("gr_all_electron_flux_E3"));
    geneva_file.Close();
    
    auto bins = std::min(pg_flux->GetN(), geneva_flux->GetN());
    std::vector<double> diff_err (bins, 0);
    std::vector<double> diff_mean (bins, 0);
    std::vector<double> ratio (bins, 0);
    std::vector<double> energy (bins, 0);

    for (int idx=0; idx<bins; ++idx)
    {
        diff_err[idx] = abs(geneva_flux->GetPointY(idx)-pg_flux->GetPointY(idx))/geneva_flux->GetErrorY(idx);
        diff_mean[idx] = abs(geneva_flux->GetPointY(idx)-pg_flux->GetPointY(idx))/(0.5*(geneva_flux->GetPointY(idx)+pg_flux->GetPointY(idx)));
        ratio[idx] = geneva_flux->GetPointY(idx)/pg_flux->GetPointY(idx);
        energy[idx] = geneva_flux->GetPointX(idx);
    }

    TGraph gr_diff(bins, &(energy[0]), &(diff_err[0]));
    TGraph gr_diff_mean(bins, &(energy[0]), &(diff_mean[0]));
    TGraph gr_ratio(bins, &(energy[0]), &(ratio[0]));

    gr_diff.SetName("gr_diff_err");
    gr_diff_mean.SetName("gr_diff_mean");
    gr_ratio.SetName("gr_ratio");

    TFile outfile(result, "RECREATE");
    if (!outfile.IsOpen())
    {
        std::cerr << "\n\nError cerating output file: [" << result << "]\n\n";
        exit(100);
    }
    gr_diff.Write();
    gr_diff_mean.Write();
    gr_ratio.Write();
    outfile.Close();

}