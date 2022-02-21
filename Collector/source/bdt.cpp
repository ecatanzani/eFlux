#include "bdt.h"

#include "TFile.h"
#include "TTree.h"

bdt::bdt(
    const std::string bdt_config_file, 
    const std::string bdt_learning_method,
    const std::string cosine_regularize_path,
    const std::string box_cox_regularize_path,
    const bool verbose)
{
    // Parse bdt config file
	get_config_info(parse_config_file(bdt_config_file));
    // Initialize BDT method
	method = bdt_learning_method;
    // Check methods with the provided one
    get_methods();
    // Initialize TMVA instance
    TMVA::Tools::Instance();
    // Initialize TMVA readers
    LE_reader = std::make_shared<TMVA::Reader>();
    ME_reader = std::make_shared<TMVA::Reader>();
    HE_reader = std::make_shared<TMVA::Reader>();
    // Link readers with the variables
    link_reader_vars(LE_reader);
    link_reader_vars(ME_reader);
    link_reader_vars(HE_reader);
    // Book MVA
    bookMVA(LE_reader, le_weights);
    bookMVA(ME_reader, me_weights);
    bookMVA(HE_reader, he_weights);
    // Load cosine angular corrections
    load_cosine_corrections(cosine_regularize_path, verbose);
    // Load box-cox lambda corrections
    load_box_cox_corrections(box_cox_regularize_path, verbose);
}

std::string bdt::parse_config_file(std::string bdt_config_file)
{
	std::ifstream input_file(bdt_config_file.c_str());
	if (!input_file.is_open())
	{
		std::cerr << "\nInput config file not found [" << bdt_config_file << "]\n\n";
		exit(100);
	}
	std::string input_string(
		(std::istreambuf_iterator<char>(input_file)),
		(std::istreambuf_iterator<char>()));
	input_file.close();
	return input_string;
}

void bdt::get_config_info(std::string parsed_config)
{
	std::string tmp_str;
	std::istringstream input_stream(parsed_config);
	std::string::size_type sz;

	while (input_stream >> tmp_str)
	{
		if (!strcmp(tmp_str.c_str(), "low_energy_weights")) input_stream >> le_weights;
		if (!strcmp(tmp_str.c_str(), "medium_energy_weights")) input_stream >> me_weights;
		if (!strcmp(tmp_str.c_str(), "high_energy_weights")) input_stream >> he_weights;
	}
}

void bdt::PrintWeights() {

	std::cout << "\n\n**** Training Weights ****\n";
	std::cout << "***********************\n\n";
	std::cout << "Low Energy  (10 GeV - 100 GeV) : " << le_weights << std::endl;
	std::cout << "Mean energy (100 GeV - 1 TeV)  : " << me_weights << std::endl;
	std::cout << "High Energy (1 TeV - 10 TeV)   : " << he_weights << std::endl;
	std::cout << "\n***********************\n";
}

void bdt::get_methods() {
   
    methods_map["Cuts"] = 0;
    methods_map["CutsD"] = 0;
    methods_map["CutsPCA"] = 0;
    methods_map["CutsGA"] = 0;
    methods_map["CutsSA"] = 0;

    methods_map["Likelihood"] = 0;
    methods_map["LikelihoodD"] = 0;
    methods_map["LikelihoodPCA"] = 0;
    methods_map["LikelihoodKDE"] = 0;
    methods_map["LikelihoodMIX"] = 0;

    methods_map["PDERS"] = 0;
    methods_map["PDERSD"] = 0;
    methods_map["PDERSPCA"] = 0;
    methods_map["PDEFoam"] = 0;
    methods_map["PDEFoamBoost"] = 0;
    methods_map["KNN"] = 0;

    methods_map["LD"] = 0;
    methods_map["Fisher"] = 0;
    methods_map["FisherG"] = 0;
    methods_map["BoostedFisher"] = 0;
    methods_map["HMatrix"] = 0;

    methods_map["FDA_GA"] = 0;
    methods_map["FDA_SA"] = 0;
    methods_map["FDA_MC"] = 0;
    methods_map["FDA_MT"] = 0;
    methods_map["FDA_GAMT"] = 0;
    methods_map["FDA_MCMT"] = 0;

    methods_map["MLP"] = 0;
    methods_map["MLPBFGS"] = 0;
    methods_map["MLPBNN"] = 0;
    methods_map["CFMlpANN"] = 0;
    methods_map["TMlpANN"] = 0;

    methods_map["SVM"] = 0;

    methods_map["BDT"] = 0;
    methods_map["BDTG"] = 0;
    methods_map["BDTB"] = 0;
    methods_map["BDTD"] = 0;
    methods_map["BDTF"] = 0;

    methods_map["RuleFit"] = 0;

    auto linked_method = false;
    for (auto &&dmethod : methods_map)
        if (!strcmp(method.c_str(), dmethod.first.c_str())) {
            dmethod.second = 1;
            linked_method = true;
            break;
        }

    if (!linked_method) {
        std::cerr << "\nERROR: No match found in TMVA default methods...\n\n";
        exit(100);
    }
}

void bdt::link_reader_vars(std::shared_ptr<TMVA::Reader> reader)
{
    reader->AddVariable("rmslayer_norm_1", &vars.rms[0]);
    reader->AddVariable("rmslayer_norm_2", &vars.rms[1]);
    reader->AddVariable("rmslayer_norm_3", &vars.rms[2]);
    reader->AddVariable("rmslayer_norm_4", &vars.rms[3]);
    reader->AddVariable("rmslayer_norm_5", &vars.rms[4]);
    reader->AddVariable("rmslayer_norm_6", &vars.rms[5]);
    reader->AddVariable("rmslayer_norm_7", &vars.rms[6]);
    reader->AddVariable("rmslayer_norm_8", &vars.rms[7]);
    reader->AddVariable("rmslayer_norm_9", &vars.rms[8]);
    reader->AddVariable("rmslayer_norm_10", &vars.rms[9]);
    reader->AddVariable("rmslayer_norm_11", &vars.rms[10]);
    reader->AddVariable("rmslayer_norm_12", &vars.rms[11]);
    reader->AddVariable("rmslayer_norm_13", &vars.rms[12]);
    reader->AddVariable("rmslayer_norm_14", &vars.rms[13]);

    reader->AddVariable("fraclayer_norm_1", &vars.fraclayer[0]);
    reader->AddVariable("fraclayer_norm_2", &vars.fraclayer[1]);
    reader->AddVariable("fraclayer_norm_3", &vars.fraclayer[2]);
    reader->AddVariable("fraclayer_norm_4", &vars.fraclayer[3]);
    reader->AddVariable("fraclayer_norm_5", &vars.fraclayer[4]);
    reader->AddVariable("fraclayer_norm_6", &vars.fraclayer[5]);
    reader->AddVariable("fraclayer_norm_7", &vars.fraclayer[6]);
    reader->AddVariable("fraclayer_norm_8", &vars.fraclayer[7]);
    reader->AddVariable("fraclayer_norm_9", &vars.fraclayer[8]);
    reader->AddVariable("fraclayer_norm_10", &vars.fraclayer[9]);
    reader->AddVariable("fraclayer_norm_11", &vars.fraclayer[10]);
    reader->AddVariable("fraclayer_norm_12", &vars.fraclayer[11]);
    reader->AddVariable("fraclayer_norm_13", &vars.fraclayer[12]);
    reader->AddVariable("fraclayer_norm_14", &vars.fraclayer[13]);

    reader->AddVariable("sumrms_norm", &vars.sumrms);
    reader->AddVariable("fraclastlayer_norm", &vars.fraclastlayer);
    reader->AddVariable("xtrl_norm", &vars.xtrl);
    reader->AddSpectator("xtrl", &vars.xtrl_spectator);
}

void bdt::bookMVA(std::shared_ptr<TMVA::Reader> reader, const std::string weights)
{
    reader->BookMVA(method.c_str(), weights.c_str());
}

void bdt::load_cosine_corrections(std::string cosine_regularize_path, const bool verbose)
{
    TFile *fitfile = TFile::Open(cosine_regularize_path.c_str(), "READ");
    if (!fitfile->IsOpen())
    {
        std::cerr << "\n\nError reading summary fit TTree [" << cosine_regularize_path << "]\n\n";
        exit(100);
    }
    else
        if (verbose)
            std::cout << "\nReading fitting summary [" << cosine_regularize_path << "]\n";

    // **** Fit function parameters
    const char* tree_name = "fit_summary";
    const char* func_form = "pol3";
    const int npars = 4;
    const double lvalue = 0;
    const double rvalue = 1; 
    // *****

    auto my_corrections_tree = static_cast<TTree*>(fitfile->Get(tree_name));
    std::vector<double>* flast_pars         {nullptr};
    std::vector<double>* flast_err_pars     {nullptr};
    std::vector<double>* sumrms_pars        {nullptr};
    std::vector<double>* sumrms_err_pars    {nullptr};
    my_corrections_tree->SetBranchAddress("flast_pars", &flast_pars);
    my_corrections_tree->SetBranchAddress("flast_err_pars", &flast_err_pars);
    my_corrections_tree->SetBranchAddress("sumrms_pars", &sumrms_pars);
    my_corrections_tree->SetBranchAddress("sumrms_err_pars", &sumrms_err_pars);

    energy_nbins = my_corrections_tree->GetEntries();
    cosine_corrections.sumrms_fitfunc.resize(energy_nbins);
    cosine_corrections.sumrms_fitfunc_err.resize(energy_nbins);
    cosine_corrections.flast_fitfunc.resize(energy_nbins);
    cosine_corrections.flast_fitfunc_err.resize(energy_nbins);

    for (int bin_idx=0; bin_idx<energy_nbins; ++bin_idx)
    {
        my_corrections_tree->GetEntry(bin_idx);

        cosine_corrections.sumrms_fitfunc[bin_idx]      = TF1(
            (std::string(tree_name) + std::string("_sumrms_fitfunc_") + std::to_string(bin_idx)).c_str(), 
            func_form, 
            lvalue, 
            rvalue);
        cosine_corrections.sumrms_fitfunc_err[bin_idx]  = TF1(
            (std::string(tree_name) + std::string("_sumrms_fitfunc_err_") + std::to_string(bin_idx)).c_str(), 
            func_form, 
            lvalue, 
            rvalue);
        cosine_corrections.flast_fitfunc[bin_idx]       = TF1(
            (std::string(tree_name) + std::string("_flast_fitfunc_") + std::to_string(bin_idx)).c_str(), 
            func_form, 
            lvalue, 
            rvalue);
        cosine_corrections.flast_fitfunc_err[bin_idx]   = TF1(
            (std::string(tree_name) + std::string("_flast_fitfunc_err_") + std::to_string(bin_idx)).c_str(), 
            func_form, 
            lvalue, 
            rvalue);

        for (int par_idx=0; par_idx<npars; ++par_idx)
        {
            cosine_corrections.sumrms_fitfunc[bin_idx]      .SetParameter(par_idx, sumrms_pars->at(par_idx));
            cosine_corrections.sumrms_fitfunc_err[bin_idx]  .SetParameter(par_idx, sumrms_err_pars->at(par_idx));
            cosine_corrections.flast_fitfunc[bin_idx]       .SetParameter(par_idx, flast_pars->at(par_idx));
            cosine_corrections.flast_fitfunc_err[bin_idx]   .SetParameter(par_idx, flast_err_pars->at(par_idx));
        }
        
    }

    fitfile->Close();
}

void bdt::load_box_cox_corrections(std::string box_cox_regularize_path, const bool verbose)
{
    TFile *lambda_tree_file = TFile::Open(box_cox_regularize_path.c_str(), "READ");
    if (!lambda_tree_file->IsOpen())
    {
        std::cerr << "\n\nError reading best lambda TTree [" << box_cox_regularize_path << "]\n\n";
        exit(100);
    }
    else
    {
        if (verbose)
            std::cout << "\n\nReading best lambda TTree [" << box_cox_regularize_path << "]";
    }

    box_cox_correction_parameters.initSize(energy_nbins);
    auto my_corrections_tree = static_cast<TTree*>(lambda_tree_file->Get("corrections_tree"));
    
    unsigned int energy_bin;
    std::vector<double>* best_rms_lambda             {nullptr};
    std::vector<double>* best_fraclayer_lambda       {nullptr};
    std::vector<double>* rms_norm_mean               {nullptr};
    std::vector<double>* rms_norm_rms                {nullptr};
    std::vector<double>* fraclayer_norm_mean         {nullptr};
    std::vector<double>* fraclayer_norm_rms          {nullptr};
    
    double best_sumrms_lambda;
    double best_fraclast_lambda;
    double best_xtrl_lambda;
    double sumrms_norm_mean;
    double sumrms_norm_rms;
    double fraclast_norm_mean;
    double fraclast_norm_rms;
    double xtrl_norm_mean;
    double xtrl_norm_rms;

    my_corrections_tree->SetBranchAddress("energy_bin", &energy_bin);
    my_corrections_tree->SetBranchAddress("best_rms_lambda", &best_rms_lambda);
    my_corrections_tree->SetBranchAddress("best_fraclayer_lambda", &best_fraclayer_lambda);
    my_corrections_tree->SetBranchAddress("rms_norm_mean", &rms_norm_mean);
    my_corrections_tree->SetBranchAddress("rms_norm_rms", &rms_norm_rms);
    my_corrections_tree->SetBranchAddress("fraclayer_norm_mean", &fraclayer_norm_mean);
    my_corrections_tree->SetBranchAddress("fraclayer_norm_rms", &fraclayer_norm_rms);
    my_corrections_tree->SetBranchAddress("best_sumrms_lambda", &best_sumrms_lambda);
    my_corrections_tree->SetBranchAddress("best_fraclast_lambda", &best_fraclast_lambda);
    my_corrections_tree->SetBranchAddress("best_xtrl_lambda", &best_xtrl_lambda);
    my_corrections_tree->SetBranchAddress("sumrms_norm_mean", &sumrms_norm_mean);
    my_corrections_tree->SetBranchAddress("sumrms_norm_rms", &sumrms_norm_rms);
    my_corrections_tree->SetBranchAddress("fraclast_norm_mean", &fraclast_norm_mean);
    my_corrections_tree->SetBranchAddress("fraclast_norm_rms", &fraclast_norm_rms);
    my_corrections_tree->SetBranchAddress("xtrl_norm_mean", &xtrl_norm_mean);
    my_corrections_tree->SetBranchAddress("xtrl_norm_rms", &xtrl_norm_rms);

    for (int bin_idx=0; bin_idx<energy_nbins; ++bin_idx)
    {
        my_corrections_tree->GetEntry(bin_idx);

        box_cox_correction_parameters.rms[energy_bin-1]                     = *best_rms_lambda;
        box_cox_correction_parameters.sumrms[energy_bin-1]                  = best_sumrms_lambda;
        box_cox_correction_parameters.fraclayer[energy_bin-1]               = *best_fraclayer_lambda;
        box_cox_correction_parameters.fraclast[energy_bin-1]                = best_fraclast_lambda;
        box_cox_correction_parameters.xtrl[energy_bin-1]                    = best_xtrl_lambda;

        box_cox_correction_parameters.rms_norm_mean[energy_bin-1]           = *rms_norm_mean;
        box_cox_correction_parameters.rms_norm_rms[energy_bin-1]            = *rms_norm_rms;

        box_cox_correction_parameters.sumrms_norm_mean[energy_bin-1]        = sumrms_norm_mean;
        box_cox_correction_parameters.sumrms_norm_rms[energy_bin-1]         = sumrms_norm_rms;

        box_cox_correction_parameters.fraclayer_norm_mean[energy_bin-1]     = *fraclayer_norm_mean;
        box_cox_correction_parameters.fraclayer_norm_rms[energy_bin-1]      = *fraclayer_norm_rms;

        box_cox_correction_parameters.fraclast_norm_mean[energy_bin-1]      = fraclast_norm_mean;
        box_cox_correction_parameters.fraclast_norm_rms[energy_bin-1]       = fraclast_norm_rms;

        box_cox_correction_parameters.xtrl_norm_mean[energy_bin-1]          = xtrl_norm_mean;
        box_cox_correction_parameters.xtrl_norm_rms[energy_bin-1]           = xtrl_norm_rms;
    }

    lambda_tree_file->Close();
}

inline double xtrl_computation(const double sumRms, const double lastFracLayer)
{
	return lastFracLayer != -1 ? 0.125e-6 * pow(sumRms, 4) * lastFracLayer : -999;
}

const double bdt::ComputeMVA(
		const std::vector<double> &rms,
		const double sumrms,
		const std::vector<double> &fraclayer,
		const double fraclastlayer,
		const double corrected_energy_gev,
        const std::vector<float> &energy_binning,
        const TVector3& bgo_direction)
        {
            // Initialize mva value
            double mva_result {-999};
            // Reset the variables struct
            vars.Reset();
            // Initialize with the new event values
            double xtrl = xtrl_computation(sumrms, fraclastlayer);
            vars.corrected_energy_gev = static_cast<float>(corrected_energy_gev);
            vars.xtrl_spectator = static_cast<float>(xtrl);
            // Apply cosine regularization
            auto get_energy_bin = [&energy_binning] (const double corrected_energy) -> int
            {
                int energybin {1};
                for (size_t idx=0; idx<energy_binning.size()-1; ++idx)
                    if (corrected_energy>=energy_binning[idx] && corrected_energy<energy_binning[idx+1])
                        energybin = idx+1;
                return energybin;
            };

            auto regularize_sumrms = [=](double sumrms, int energy_bin, TVector3 bgodir) -> double 
            {
                // Initialize regularized sumrms variable
                double reg_sumrms = sumrms;
                // Initialize BGO cosine from directrion
                double bgocosine = bgodir.CosTheta();
                // Regularize sumrms
                reg_sumrms -= cosine_corrections.sumrms_fitfunc[energy_bin - 1].Eval(bgocosine);
                reg_sumrms /= cosine_corrections.sumrms_fitfunc_err[energy_bin - 1].Eval(bgocosine);
                return reg_sumrms;
            };

            auto regularize_flast = [=](double flast, int energy_bin, TVector3 bgodir) -> double 
            {
                // Initialize regularized sumrms variable
                double reg_flast = flast;
                // Initialize BGO cosine from directrion
                double bgocosine = bgodir.CosTheta();
                // Regularize sumrms
                reg_flast -= cosine_corrections.flast_fitfunc[energy_bin - 1].Eval(bgocosine);
                reg_flast /= cosine_corrections.flast_fitfunc_err[energy_bin - 1].Eval(bgocosine);
                return reg_flast;
            };

            vars.sumrms = static_cast<float>(regularize_sumrms(sumrms, get_energy_bin(corrected_energy_gev), bgo_direction));
            vars.fraclastlayer = static_cast<float>(regularize_flast(fraclastlayer, get_energy_bin(corrected_energy_gev), bgo_direction));
            // Apply box-cox regularization
            auto gaussianize_rmslayer = [=](const std::vector<double> input_rmslayer, const int energy_bin) -> std::vector<float>
            {
                auto gaussianize_elm = [](const std::vector<double> elm, const std::vector<double> lambda) -> std::vector<float> 
                {
                    std::vector<float> elm_cp (elm.size(), 0);
                    for (unsigned int idx=0; idx<elm.size(); ++idx)
                        elm_cp[idx] = static_cast<float>(lambda[idx] ? (exp(lambda[idx]*elm[idx])-1)/lambda[idx] : elm[idx]);
                    return elm_cp;
                };
                
                return gaussianize_elm(input_rmslayer, box_cox_correction_parameters.rms[energy_bin-1]);
            };

            auto gaussianize_sumrms = [=](const double input_sumrmslayer, const int energy_bin) -> float
            {
                auto gaussianize_elm = [](const double elm, const double lambda) -> float
                {
                    float elm_cp = static_cast<float>(lambda ? (exp(lambda*elm)-1)/lambda : elm);
                    return elm_cp;
                };
                
                return gaussianize_elm(input_sumrmslayer, box_cox_correction_parameters.sumrms[energy_bin-1]);
            };

            auto gaussianize_fraclayer = [=](const std::vector<double> input_fraclayer, const int energy_bin) -> std::vector<float>
            {
                auto gaussianize_elm = [](const std::vector<double> elm, const std::vector<double> lambda) -> std::vector<float> 
                {
                    std::vector<float> elm_cp (elm.size(), 0);
                    for (unsigned int idx=0; idx<elm.size(); ++idx)
                        elm_cp[idx] = static_cast<float>(lambda[idx] ? (exp(lambda[idx]*elm[idx])-1)/lambda[idx] : elm[idx]);
                    return elm_cp;
                };
                
                return gaussianize_elm(input_fraclayer, box_cox_correction_parameters.fraclayer[energy_bin-1]);
            };

            auto gaussianize_fraclastlayer = [=](const double input_fraclayer, const int energy_bin) -> float
            {
                auto gaussianize_elm = [](const double elm, const double lambda) -> float
                {
                    float elm_cp = static_cast<float>(lambda ? (exp(lambda*elm)-1)/lambda : elm);
                    return elm_cp;
                };

                return gaussianize_elm(input_fraclayer, box_cox_correction_parameters.fraclast[energy_bin-1]);
            };
            
            auto gaussianize_xtrl = [=](const double input_xtrl, const int energy_bin) -> float
            {
                auto gaussianize_elm = [](const double elm, const double lambda) -> float
                {
                    float elm_cp = static_cast<float>(lambda ? (exp(lambda*elm)-1)/lambda : elm);
                    return elm_cp;
                };

                return gaussianize_elm(input_xtrl, box_cox_correction_parameters.xtrl[energy_bin-1]);
            };

            vars.rms = gaussianize_rmslayer(rms, get_energy_bin(corrected_energy_gev));
            vars.sumrms = gaussianize_sumrms(vars.sumrms, get_energy_bin(corrected_energy_gev));
            vars.fraclayer = gaussianize_fraclayer(fraclayer, get_energy_bin(corrected_energy_gev));
            vars.fraclastlayer = gaussianize_fraclastlayer(vars.fraclastlayer, get_energy_bin(corrected_energy_gev));
            vars.xtrl = gaussianize_xtrl(xtrl, get_energy_bin(corrected_energy_gev));
            // Compute the MVA
            if (vars.corrected_energy_gev>=10 && vars.corrected_energy_gev<100) {
                mva_result = LE_reader->EvaluateMVA(method.c_str());
            }
            else if (vars.corrected_energy_gev>=100 && vars.corrected_energy_gev<1000) {
                mva_result = ME_reader->EvaluateMVA(method.c_str());
            }
            else if (vars.corrected_energy_gev>=1000 && vars.corrected_energy_gev<=10000) {
                mva_result = HE_reader->EvaluateMVA(method.c_str());
            }

            return mva_result;
        }