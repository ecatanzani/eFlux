#include "acceptance.h"

#include <fstream>
#include <sstream>

void load_acceptance_struct(
    cuts_conf &acceptance_cuts,
    data_active_cuts &active_cuts,
    const std::string wd)
{
    //std::string cwd = GetCurrentWorkingDir();
    std::string cwd = wd;
    std::size_t index = cwd.find("eFlux");
    std::string configPath = cwd.substr(0, index + 5);
    configPath += "/config/acceptanceConfig.txt";
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
        // Load cuts variables
        if (!strcmp(tmp_str.c_str(), "min_event_energy"))
        {
            input_stream >> tmp_str;
            acceptance_cuts.min_event_energy = stod(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "max_event_energy"))
        {
            input_stream >> tmp_str;
            acceptance_cuts.max_event_energy = stod(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "energy_lRatio"))
        {
            input_stream >> tmp_str;
            acceptance_cuts.energy_lRatio = stod(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "shower_axis_delta"))
        {
            input_stream >> tmp_str;
            acceptance_cuts.shower_axis_delta = stoi(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "max_rms_shower_width"))
        {
            input_stream >> tmp_str;
            acceptance_cuts.max_rms_shower_width = stoi(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "generation_vertex_radius"))
        {
            input_stream >> tmp_str;
            acceptance_cuts.vertex_radius = stod(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "track_X_clusters"))
        {
            input_stream >> tmp_str;
            acceptance_cuts.track_X_clusters = stoi(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "track_Y_clusters"))
        {
            input_stream >> tmp_str;
            acceptance_cuts.track_Y_clusters = stoi(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "track_missingHit_X"))
        {
            input_stream >> tmp_str;
            acceptance_cuts.track_missingHit_X = stoi(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "track_missingHit_Y"))
        {
            input_stream >> tmp_str;
            acceptance_cuts.track_missingHit_Y = stoi(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "STK_BGO_delta_track"))
        {
            input_stream >> tmp_str;
            acceptance_cuts.STK_BGO_delta_track = stoi(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "STK_BGO_delta_position"))
        {
            input_stream >> tmp_str;
            acceptance_cuts.STK_BGO_delta_position = stoi(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "xtrl"))
        {
            input_stream >> tmp_str;
            acceptance_cuts.xtrl = stod(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "STK_PSD_delta_position"))
        {
            input_stream >> tmp_str;
            acceptance_cuts.STK_PSD_delta_position = stoi(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "PSD_bar_min_energy_release"))
        {
            input_stream >> tmp_str;
            acceptance_cuts.PSD_bar_min_energy_release = stod(tmp_str, &sz);
        }

        // Load cuts
        if (!strcmp(tmp_str.c_str(), "geometry"))
        {
            input_stream >> tmp_str;
            if (!strcmp(tmp_str.c_str(), "YES") || !strcmp(tmp_str.c_str(), "yes"))
                active_cuts.geometry = true;
        }
        if (!strcmp(tmp_str.c_str(), "BGO_fiducial"))
        {
            input_stream >> tmp_str;
            if (!strcmp(tmp_str.c_str(), "YES") || !strcmp(tmp_str.c_str(), "yes"))
                active_cuts.BGO_fiducial = true;
        }
        if (!strcmp(tmp_str.c_str(), "nBarLayer13"))
        {
            input_stream >> tmp_str;
            if (!strcmp(tmp_str.c_str(), "YES") || !strcmp(tmp_str.c_str(), "yes"))
            {
                active_cuts.nBarLayer13 = true;
                ++active_cuts.nActiveCuts;
            }
        }
        if (!strcmp(tmp_str.c_str(), "maxRms"))
        {
            input_stream >> tmp_str;
            if (!strcmp(tmp_str.c_str(), "YES") || !strcmp(tmp_str.c_str(), "yes"))
            {
                active_cuts.maxRms = true;
                ++active_cuts.nActiveCuts;
            }
        }
        if (!strcmp(tmp_str.c_str(), "track_selection"))
        {
            input_stream >> tmp_str;
            if (!strcmp(tmp_str.c_str(), "YES") || !strcmp(tmp_str.c_str(), "yes"))
            {
                active_cuts.track_selection = true;
                ++active_cuts.nActiveCuts;
            }
        }
        if (!strcmp(tmp_str.c_str(), "xtrl_selection"))
        {
            input_stream >> tmp_str;
            if (!strcmp(tmp_str.c_str(), "YES") || !strcmp(tmp_str.c_str(), "yes"))
            {
                active_cuts.xtrl = true;
                ++active_cuts.nActiveCuts;
            }
        }
        if (!strcmp(tmp_str.c_str(), "psd_charge"))
        {
            input_stream >> tmp_str;
            if (!strcmp(tmp_str.c_str(), "YES") || !strcmp(tmp_str.c_str(), "yes"))
            {
                active_cuts.psd_charge = true;
                ++active_cuts.nActiveCuts;
            }
        }
    }
}
