#include "binning.h"

inline std::vector<float> createLogBinning(
	const double eMin,
	const double eMax,
	const std::size_t n_bins)
{
	std::vector<float> binning(n_bins + 1, 0);
	double log_interval = (log10(eMax) - log10(eMin)) / n_bins;
	for (unsigned int bIdx = 0; bIdx <= n_bins; ++bIdx)
		binning[bIdx] = pow(10, log10(eMin) + bIdx * log_interval);

	return binning;
}

std::vector<float> read_binning_from_config(
    const std::string wd)
{
    //std::string cwd = GetCurrentWorkingDir();
    std::string configPath = wd;
    configPath += "/fluxConfig.txt";
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
    double min_event_energy, max_event_energy = 0;

    while (input_stream >> tmp_str)
    {
        // Load flux params
        if (!strcmp(tmp_str.c_str(), "n_energy_bins"))
            input_stream >> n_bins;

        // Load cuts variables
        if (!strcmp(tmp_str.c_str(), "min_event_energy"))
        {
            input_stream >> tmp_str;
            min_event_energy = stod(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "max_event_energy"))
        {
            input_stream >> tmp_str;
            max_event_energy = stod(tmp_str, &sz);
        }
    }

    // Build energy binning vector
    return createLogBinning(
        min_event_energy, 
        max_event_energy, 
        n_bins);
}
