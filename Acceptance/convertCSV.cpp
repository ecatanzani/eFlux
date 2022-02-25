#include <tuple>
#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>

#include "TFile.h"
#include "TGraph.h"

std::tuple<std::vector<double>, std::vector<double>> extract_data(const char* input_file) {
    
    std::ifstream file(input_file);
	if (!file.is_open())
	{
		std::cerr << "\nInput config file not found [" << input_file << "]\n\n";
		exit(100);
	}

    std::string line, item;
    std::vector<double> x_values, y_values;
    unsigned int line_counter {0};
    bool x {true};

    while (getline(file, line)) {
        istringstream lineReader(line);
        while (getline(lineReader, item, ';')) {
            std::replace(std::begin(item), std::end(item), ',', '.');
            if (x) {
                x_values.push_back(stod(item));
                x = false;
            }
            else {
                y_values.push_back(stod(item));
                x = true;
            }
        }
        ++line_counter;
    }

    std::cout << "\nNumber of lines extracted from file: " << line_counter << std::endl;
    return std::make_tuple(x_values, y_values);
}   

void convertCSV(const char* input_file, const char* output_file) {
    auto extracted_data {extract_data(input_file)};

    TFile *outfile = TFile::Open(output_file, "RECREATE");
    if (outfile->IsZombie()) {
        std::cerr << "\n\nError opening ROOT output file\n\n";
        exit(100);
    }

    TGraph gr(std::get<0>(extracted_data).size(), &(std::get<0>(extracted_data))[0], &(std::get<1>(extracted_data))[0]);
    gr.SetName("convertedCSV");
    gr.Write();
    outfile->Close();
}