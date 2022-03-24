#ifndef VARS_H
#define VARS_H

#include <vector>
#include "Dmp/DmpGeoStruct.h"

#include "TF1.h"
#include "TVector3.h"

struct cosine_correction_functions
{
	std::vector<TF1> sumrms_fitfunc;
    std::vector<TF1> sumrms_fitfunc_err;
    std::vector<TF1> flast_fitfunc;
    std::vector<TF1> flast_fitfunc_err;
};

struct best_lambda
{
    std::vector<std::vector<double>> rms;
    std::vector<double> sumrms;
    std::vector<std::vector<double>> fraclayer;
    std::vector<double> fraclast;
    std::vector<double> xtrl;

    std::vector<std::vector<double>> rms_norm_mean;
    std::vector<std::vector<double>> rms_norm_rms;
    std::vector<double> sumrms_norm_mean;
    std::vector<double> sumrms_norm_rms;
    std::vector<std::vector<double>> fraclayer_norm_mean;
    std::vector<std::vector<double>> fraclayer_norm_rms;
    std::vector<double> fraclast_norm_mean;
    std::vector<double> fraclast_norm_rms;
    std::vector<double> xtrl_norm_mean;
    std::vector<double> xtrl_norm_rms;

    void initSize(const int size)
    {
        rms							.resize(size);
        sumrms						.resize(size);
        fraclayer					.resize(size);
        fraclast					.resize(size);
        xtrl						.resize(size);

        rms_norm_mean				.resize(size);
        rms_norm_rms				.resize(size);
        sumrms_norm_mean			.resize(size);
        sumrms_norm_rms				.resize(size);
        fraclayer_norm_mean			.resize(size);
        fraclayer_norm_rms			.resize(size);
        fraclast_norm_mean			.resize(size);
        fraclast_norm_rms			.resize(size);
        xtrl_norm_mean				.resize(size);
        xtrl_norm_rms				.resize(size);

        for (int idx=0; idx<size; ++idx)
        {
            rms[idx] 					= std::vector<double> (DAMPE_bgo_nLayers, 0);
            sumrms[idx] 				= 0;
            fraclayer[idx] 				= std::vector<double> (DAMPE_bgo_nLayers, 0);
            fraclast[idx] 				= 0;
            xtrl[idx] 					= 0;

            rms_norm_mean[idx] 			= std::vector<double> (DAMPE_bgo_nLayers, 0);
            rms_norm_rms[idx] 			= std::vector<double> (DAMPE_bgo_nLayers, 0);
            sumrms_norm_mean[idx] 		= 0;
            sumrms_norm_rms[idx] 		= 0;
            fraclayer_norm_mean[idx] 	= std::vector<double> (DAMPE_bgo_nLayers, 0);
            fraclayer_norm_rms[idx] 	= std::vector<double> (DAMPE_bgo_nLayers, 0);
            fraclast_norm_mean[idx] 	= 0;
            fraclast_norm_rms[idx] 		= 0;
            xtrl_norm_mean[idx] 		= 0;
            xtrl_norm_rms[idx] 			= 0;
        }

    };
};

struct bdt_vars
{
	std::vector<double> rms;
	double sumrms;
	std::vector<double> fraclayer;
	double fraclastlayer;
	double xtrl;
	double xtrl_spectator;
	double corrected_energy_gev;

	void Reset()
	{
		rms 						= std::vector<double> (DAMPE_bgo_nLayers, 0);
		sumrms 						= -999;
		fraclayer 					= std::vector<double> (DAMPE_bgo_nLayers, 0);
		fraclastlayer 				= -999;
		xtrl 						= -999;
		xtrl_spectator				= -999;
		corrected_energy_gev		= -999;
	}
};

class vars
{
public:
	vars(
        const std::string cosine_regularize_path,
		const std::string box_cox_regularize_path,
		const bool verbose);
	~vars(){};
	
	void PrintWeights();

	const bdt_vars GetVars(
		const std::vector<double> &rms,
		const double sumrms,
		const std::vector<double> &fraclayer,
		const double fraclastlayer,
		const double corrected_energy_gev,
		const std::vector<float> &energy_binning,
		const TVector3& bgo_direction);

private:
	void load_cosine_corrections(std::string cosine_regularize_path, const bool verbose);
	void load_box_cox_corrections(std::string box_cox_regularize_path, const bool verbose);
    
    int energy_nbins;
    
	cosine_correction_functions cosine_corrections;
	best_lambda box_cox_correction_parameters;
	bdt_vars my_vars;
};

#endif