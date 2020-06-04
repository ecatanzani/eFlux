#include "energy_cuts.h"

void load_energy_struct(energy_cuts &e_cuts, const std::string wd)
{
    //std::string cwd = GetCurrentWorkingDir();
    std::string cwd = wd;
    std::size_t index = cwd.find("eFlux");
    std::string configPath = cwd.substr(0, index + 5);
    configPath += "/config/energyConfig.txt";
    std::ifstream input_file(configPath.c_str());
    if (!input_file.is_open())
    {
        std::cerr << "\nERROR 100! File not open " << configPath << "\n\n";
        exit(100);
    }
    std::string input_string((std::istreambuf_iterator<char>(input_file)), (std::istreambuf_iterator<char>()));
    input_file.close();
    std::string tmp_str;
    std::istringstream input_stream(input_string);
    std::string::size_type sz;
    while (input_stream >> tmp_str)
    {
        if (!strcmp(tmp_str.c_str(), "cut_min_energy"))
        {
            input_stream >> tmp_str;
            e_cuts.cut_min_energy = stod(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "cut_max_energy"))
        {
            input_stream >> tmp_str;
            e_cuts.cut_max_energy = stod(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "e_bins"))
        {
            input_stream >> tmp_str;
            e_cuts.e_bins = stoi(tmp_str, &sz);
        }
    }
}
