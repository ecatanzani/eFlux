#include "flux.h"
#include "myHeader.h"

void load_flux_struct(flux_conf &flux_params)
{
    std::string cwd = GetCurrentWorkingDir();
    std::size_t index = cwd.find("eFlux");
    std::string configPath = cwd.substr(0, index + 5);
    configPath += "/config/fluxConfig.txt";
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
        if (!strcmp(tmp_str.c_str(),"min_energy"))
        {
            input_stream >> tmp_str;
            flux_params.min_energy = stod(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(),"max_energy"))
        {
            input_stream >> tmp_str;
            flux_params.max_energy = stod(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(),"number_energy_bins"))
        {
            input_stream >> tmp_str;
            flux_params.n_energy_bins = stoi(tmp_str, &sz);
        }
    }
}