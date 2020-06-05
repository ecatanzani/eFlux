#include "read_sets_config_file.h"

void load_input_dsets_config(data_set_conf &input_sets, const std::string wd)
{
    std::string cwd = wd;
    std::size_t index = cwd.find("eFlux");
    std::string configPath = cwd.substr(0, index + 5);
    configPath += "/config/powerlowMCsets.txt";
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
        if (!strcmp(tmp_str.c_str(), "Set:"))
        {
            // Read dataset name
            input_stream >> tmp_str;
            input_sets.sets_name.push_back(tmp_str);
            // Read dataset min energy
            input_stream >> tmp_str;
            input_sets.sets_eMin.push_back(stod(tmp_str, &sz));
            // Read dataset max energy
            input_stream >> tmp_str;
            input_sets.sets_eMax.push_back(stod(tmp_str, &sz));
            // Read dataset powerlaw input energy spectrum
            input_stream >> tmp_str;
            input_sets.sets_powerlaw_idx.push_back(stoi(tmp_str, &sz));
        }
    }
}