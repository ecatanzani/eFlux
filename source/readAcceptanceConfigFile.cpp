#include "acceptance.h"
#include "myHeader.h"

std::string GetCurrentWorkingDir( void ) 
{
    char buff[FILENAME_MAX];
    GetCurrentDir( buff, FILENAME_MAX );
    std::string current_working_dir(buff);
    return current_working_dir;
}

void load_acceptance_struct(acceptance_conf &acceptance_cuts)
{
    std::string cwd = GetCurrentWorkingDir();
    std::size_t index = cwd.find("eFlux");
    std::string configPath = cwd.substr(0,index+5);
    configPath += "/config/acceptanceConfig.txt";
    std::ifstream input_file(configPath.c_str());
    if(!input_file.is_open()) {
        std::cerr << "\nERROR 100! File not open " << configPath << "\n\n";
        exit(100);
    }
    std::string input_string((std::istreambuf_iterator< char >(input_file)), (std::istreambuf_iterator< char >()));
    input_file.close();
    std::string tmp_str;
    std::istringstream input_stream(input_string);
    std::string::size_type sz;
    while(input_stream>>tmp_str)
    {
        if(!strcmp(tmp_str.c_str(),"event_energy"))
        {
            input_stream>>tmp_str;
            acceptance_cuts.event_energy = stod(tmp_str,&sz);
        }
        if(!strcmp(tmp_str.c_str(),"energy_lRatio"))
        {
            input_stream>>tmp_str;
            acceptance_cuts.energy_lRatio = stod(tmp_str,&sz);
        }
        if(!strcmp(tmp_str.c_str(),"shower_axis_delta"))
        {
            input_stream>>tmp_str;
            acceptance_cuts.shower_axis_delta = stoi(tmp_str,&sz);
        }
    }
}

