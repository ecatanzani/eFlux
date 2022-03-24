#ifndef BDT_H
#define BDT_H

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>

#include "TF1.h"
#include "TVector3.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include "Dmp/DmpGeoStruct.h"

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

struct classifier_vars {
    float rmslayer_norm_1 {0};
    float rmslayer_norm_2 {0};
    float rmslayer_norm_3 {0};
    float rmslayer_norm_4 {0};
    float rmslayer_norm_5 {0};
    float rmslayer_norm_6 {0};
    float rmslayer_norm_7 {0};
    float rmslayer_norm_8 {0};
    float rmslayer_norm_9 {0};
    float rmslayer_norm_10 {0};
    float rmslayer_norm_11 {0};
    float rmslayer_norm_12 {0};
    float rmslayer_norm_13 {0};
    float rmslayer_norm_14 {0};
    float fraclayer_norm_1 {0};
    float fraclayer_norm_2 {0};
    float fraclayer_norm_3 {0};
    float fraclayer_norm_4 {0};
    float fraclayer_norm_5 {0};
    float fraclayer_norm_6 {0};
    float fraclayer_norm_7 {0};
    float fraclayer_norm_8 {0};
    float fraclayer_norm_9 {0};
    float fraclayer_norm_10 {0};
    float fraclayer_norm_11 {0};
    float fraclayer_norm_12 {0};
    float fraclayer_norm_13 {0};
    float fraclayer_norm_14 {0};
    float sumrms_norm {0};
    float fraclastlayer_norm {0};
    float xtrl_norm {0};
    float xtrl {0};
};

class bdt
{
public:
	bdt(
		const std::string bdt_config_file, 
		const std::string bdt_learning_method,
		const std::string cosine_regularize_path,
		const std::string box_cox_regularize_path,
		const bool verbose);
	~bdt(){};
	
	void PrintWeights();

	const double ComputeMVA(
		const std::vector<double> &rms,
		const double sumrms,
		const std::vector<double> &fraclayer,
		const double fraclastlayer,
		const double corrected_energy_gev,
		const std::vector<float> &energy_binning,
		const TVector3& bgo_direction);

    const double GetLowEnergyBDTCut();
    const double GetMidEnergyBDTCut();
    const double GetHighEnergyBDTCut();

private:
	std::string parse_config_file(std::string bdt_config_file);
	void get_config_info(std::string parsed_config);
	void get_methods();
	void link_reader_vars(std::shared_ptr<TMVA::Reader> reader);
    void sync_vars();
    void bookMVA(std::shared_ptr<TMVA::Reader> reader, const std::string weights);
	void load_cosine_corrections(std::string cosine_regularize_path, const bool verbose);
	void load_box_cox_corrections(std::string box_cox_regularize_path, const bool verbose);
	
	std::string le_weights;
	std::string me_weights;
	std::string he_weights;

    double le_c_cut {0};
    double me_c_cut {0};
    double he_c_cut {0};
	
	std::string method;
	int energy_nbins;
	std::map<std::string, int> methods_map;

	std::shared_ptr<TMVA::Reader> LE_reader;
    std::shared_ptr<TMVA::Reader> ME_reader;
    std::shared_ptr<TMVA::Reader> HE_reader;

	cosine_correction_functions cosine_corrections;
	best_lambda box_cox_correction_parameters;
	bdt_vars vars;
    classifier_vars tmva_vars;
};

#endif