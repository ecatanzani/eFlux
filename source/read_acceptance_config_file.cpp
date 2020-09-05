#include "acceptance.h"
#include "mc_ancillary.h"
#include "binning.h"

#include <fstream>
#include <sstream>
 
std::vector<float> load_acceptance_struct(
    cuts_conf &acceptance_cuts,
    data_active_cuts &active_cuts,
    mc_ancillary_cuts &ancillary_cuts,
    const std::string wd)
{
    std::string configPath = wd;
    configPath += "/acceptanceConfig.txt";
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
    std::size_t n_bins;

    while (input_stream >> tmp_str)
    {
        // Load acceptance params
        if (!strcmp(tmp_str.c_str(), "n_energy_bins"))
            input_stream >> n_bins;

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
        if (!strcmp(tmp_str.c_str(), "layer_min_energy"))
        {
            input_stream >> tmp_str;
            acceptance_cuts.layer_min_energy = stod(tmp_str, &sz);
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
        if (!strcmp(tmp_str.c_str(), "PSD_charge"))
        {
            input_stream >> tmp_str;
            acceptance_cuts.PSD_charge = stod(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "PSD_charge_sum"))
        {
            input_stream >> tmp_str;
            acceptance_cuts.PSD_charge_sum = stod(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "STK_charge"))
        {
            input_stream >> tmp_str;
            acceptance_cuts.STK_charge = stod(tmp_str, &sz);
        }

        // Load cuts
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
        if (!strcmp(tmp_str.c_str(), "psd_stk_match"))
        {
            input_stream >> tmp_str;
            if (!strcmp(tmp_str.c_str(), "YES") || !strcmp(tmp_str.c_str(), "yes"))
            {
                active_cuts.psd_stk_match = true;
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
        if (!strcmp(tmp_str.c_str(), "stk_charge"))
        {
            input_stream >> tmp_str;
            if (!strcmp(tmp_str.c_str(), "YES") || !strcmp(tmp_str.c_str(), "yes"))
            {
                active_cuts.stk_charge = true;
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

        // Load ancillary cuts
        if (!strcmp(tmp_str.c_str(), "compute_proton_background"))
        {
            input_stream >> tmp_str;
            if (!strcmp(tmp_str.c_str(), "YES") || !strcmp(tmp_str.c_str(), "yes"))
                ancillary_cuts.compute_proton_background = true;
        }
    }

    // Build energy binning vector
    return createLogBinning(
        acceptance_cuts.min_event_energy, 
        acceptance_cuts.max_event_energy, 
        n_bins);
}
