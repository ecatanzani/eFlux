#include <iostream>
#include <vector>

#include "TFile.h"
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/RInterface.hxx>

void editOnTheFlight(const char* input_file, const char* output_file, const char* tree_name)
{
    TFile *inputfile = TFile::Open(input_file, "READ");
    if (!inputfile->IsOpen())
    {
        std::cerr << "\n\nError opening input file [" << input_file << "]";
        exit(100);
    }
    auto tree = static_cast<TTree*>(inputfile->Get(tree_name));

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame _fr(*tree);

    std::cout << "\n\nRDataFrame comumn names...\n\n";
    auto colNames = _fr.GetColumnNames();
    for (auto &&colName : colNames) std::cout << colName << std::endl;
    std::cout << "\n*********\n\n";

    // Create missing NUD ADC leafs
    auto build_adc_vector = [] (const int adc_1, const int adc_2, const int adc_3, const int adc_4) -> std::vector<int>
    {
        std::vector<int> NUD_ADC(4);
        NUD_ADC[0] = adc_1;
        NUD_ADC[1] = adc_2;
        NUD_ADC[2] = adc_3;
        NUD_ADC[3] = adc_4;
        return NUD_ADC;
    };

    auto _fr_preselected = _fr.Define("NUD_ADC", build_adc_vector, {"NUD_ADC_1", "NUD_ADC_2", "NUD_ADC_3", "NUD_ADC_4"});
    _fr_preselected = _fr_preselected.Define("NUD_ADC_min", [](std::vector<int> adc) -> int { return *std::min_element(adc.begin(), adc.end()); }, {"NUD_ADC"});
    _fr_preselected = _fr_preselected.Define("NUD_ADC_max", [](std::vector<int> adc) -> int { return *std::max_element(adc.begin(), adc.end()); }, {"NUD_ADC"});
    _fr_preselected = _fr_preselected.Define("NUD_ADC_rms",
                                             [](std::vector<int> adc) -> double {
                                                 double sumq = 0;
                                                 for (const auto &elm : adc)
                                                     sumq += pow(elm, 2);
                                                 return sqrt(sumq / adc.size());
                                             }, {"NUD_ADC"});

    _fr_preselected.Snapshot(
        tree_name,
        output_file,
        {"STK_bestTrack_npoints", "energy", "energy_corr", "sumRms", "sumRms_reg", "fracLast", "fracLast_reg",
         "rmsLayer_1", "rmsLayer_2", "rmsLayer_3", "rmsLayer_4", "rmsLayer_5", "rmsLayer_6", "rmsLayer_7", "rmsLayer_8",
         "rmsLayer_9", "rmsLayer_10", "rmsLayer_11", "rmsLayer_12", "rmsLayer_13", "rmsLayer_14", "fracLayer_1", "fracLayer_2",
         "fracLayer_3", "fracLayer_4", "fracLayer_5", "fracLayer_6", "fracLayer_7", "fracLayer_8", "fracLayer_9", "fracLayer_10",
         "fracLayer_11", "fracLayer_12", "fracLayer_13", "fracLayer_14", "lastBGOLayer", "nBGOentries",
         "energy_1R_radius_1", "energy_1R_radius_2", "energy_1R_radius_3", "energy_1R_radius_4", "energy_1R_radius_5", "energy_1R_radius_6",
         "energy_1R_radius_7", "energy_1R_radius_8", "energy_1R_radius_9", "energy_1R_radius_10", "energy_1R_radius_11", "energy_1R_radius_12",
         "energy_1R_radius_13", "energy_1R_radius_14",
         "energy_2R_radius_1", "energy_2R_radius_2", "energy_2R_radius_3", "energy_2R_radius_4", "energy_2R_radius_5", "energy_2R_radius_6",
         "energy_2R_radius_7", "energy_2R_radius_8", "energy_2R_radius_9", "energy_2R_radius_10", "energy_2R_radius_11", "energy_2R_radius_12",
         "energy_2R_radius_13", "energy_2R_radius_14",
         "energy_3R_radius_1", "energy_3R_radius_2", "energy_3R_radius_3", "energy_3R_radius_4", "energy_3R_radius_5", "energy_3R_radius_6",
         "energy_3R_radius_7", "energy_3R_radius_8", "energy_3R_radius_9", "energy_3R_radius_10", "energy_3R_radius_11", "energy_3R_radius_12",
         "energy_3R_radius_13", "energy_3R_radius_14",
         "energy_5R_radius_1", "energy_5R_radius_2", "energy_5R_radius_3", "energy_5R_radius_4", "energy_5R_radius_5", "energy_5R_radius_6",
         "energy_5R_radius_7", "energy_5R_radius_8", "energy_5R_radius_9", "energy_5R_radius_10", "energy_5R_radius_11", "energy_5R_radius_12",
         "energy_5R_radius_13", "energy_5R_radius_14",
         "xtrl", "NUD_ADC_1", "NUD_ADC_2", "NUD_ADC_3", "NUD_ADC_4", "NUD_total_ADC_nud_total_adc", "NUD_max_ADC_nud_max_adc", 
         "NUD_ADC_min", "NUD_ADC_max", "NUD_ADC_rms"});
}