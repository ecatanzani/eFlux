#ifndef CONFIG_H
#define CONFIG_H

#include <string>

struct cuts_config {
    double bgo_layer_min_energy     {0};    // Minimum energy per BGO layer
    double bgo_max_energy_ratio     {0};    // Maximum energy ratio per layer
    double bgo_shower_axis_delta    {0};    // BGO maximum shower axis delta (mm)
    double bgo_shower_width         {0};    // BGO maximum shower width (mm)
    double STK_BGO_delta_position   {0};    // Linear distance between STK and BGO projections
    double STK_BGO_delta_track      {0};    // Angular distance between BGO/STK tracks (deg)
    int track_X_clusters            {0};    // Number of requested clusters on X tracks
    int track_Y_clusters            {0};    // Number of requested clusters on Y tracks
    int track_X_holes               {0};    // Number of holes on X view
    int track_Y_holes               {0};    // Number of holes on Y view
    double STK_PSD_delta_position   {0};    // Linear distance between STK progection and PSD seed strip
    double psd_min_energy           {0};    // Minimum energy per PSD bar
    double PSD_sharge_sum           {0};    // Sum of PSD charges on X and Y views
    double PSD_single_charge        {0};    // PSD charge cut on single view
    double STK_single_charge        {0};    // STK charge cut on single view
};

struct event_display {
    double evt_min_energy {0};
    double evt_max_energy {0};
    int min_track_X_clusters {0};
    int min_track_Y_clusters {0};
    int track_X_clusters {0};
    int track_Y_clusters {0};
};

class config {
    public:
        config(const std::string energy_config_path);
        ~config(){};
        const cuts_config GetCutsConfig();
        const event_display GetEventDisplayConfig();

    private:
        std::string parse_config_file(std::string wd,std::string config_file);
        void get_config_info(std::string parsed_config);

        cuts_config cuts;
        event_display event_display_config;
};

#endif