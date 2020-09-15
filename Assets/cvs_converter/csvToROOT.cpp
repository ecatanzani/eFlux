#include "TFile.h"
#include "TGraph.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

void fromCSVtoROOT(const char* inCSVpath, const char* outPath)
{
    std::ifstream input_file(inCSVpath);
    if (!input_file.is_open())
    {
        std::cerr << "\nERROR 100! File not open " << inCSVpath << "\n\n";
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
    
    TFile outFile(outPath, "RECREATE");
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