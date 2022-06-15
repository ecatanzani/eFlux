//#include "vars.h"
//#include "bdt_support.h"
#include "bdt.h"
#include "efficiency.h"
#include "list_parser.h"
#include "energy_config.h"

#include <memory>
#include <tuple>

#include "TF1.h"
#include "TFile.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TGraphErrors.h"
#include <ROOT/RDataFrame.hxx>

inline std::tuple<std::tuple<std::shared_ptr<TF1>, std::shared_ptr<TF1>>, std::tuple<std::shared_ptr<TGraphErrors>, std::shared_ptr<TGraphErrors>>> extractCorrectionFunctions(
    const char* file_name, 
    const char* tf1_shift_name="shift_fit_function",
    const char* tf1_sigma_name="sigma_ratio_fit_function",
    const char* gr_shift_name="gr_bin_shift",
    const char* gr_sigma_ratio_name="gr_bin_ratio") {

        TFile *infile = TFile::Open(file_name, "READ");
        if (infile->IsZombie()) {
            std::cout << "Error opening input TF1 file: [" << file_name << "]\n\n";
            exit(100);
        }

        std::shared_ptr<TF1> tf1_shift = std::shared_ptr<TF1>(static_cast<TF1*>(infile->Get(tf1_shift_name)));
        std::shared_ptr<TF1> tf1_sigma = std::shared_ptr<TF1>(static_cast<TF1*>(infile->Get(tf1_sigma_name)));
        std::shared_ptr<TGraphErrors> gr_shift = std::shared_ptr<TGraphErrors>(static_cast<TGraphErrors*>(infile->Get(gr_shift_name)));
        std::shared_ptr<TGraphErrors> gr_sigma = std::shared_ptr<TGraphErrors>(static_cast<TGraphErrors*>(infile->Get(gr_sigma_ratio_name)));

        return std::make_tuple(std::make_tuple(tf1_shift, tf1_sigma), std::make_tuple(gr_shift, gr_sigma));
    }

void buildEfficiency(const in_args input_args)
{     
    // Parse input file list
    std::shared_ptr<parser> evt_parser = std::make_unique<parser>(input_args.input_list, input_args.verbose);
    auto isMC = [](const char* tree_name) -> bool {
        std::size_t found = std::string(tree_name).find(std::string("MC"));
        if (found!=std::string::npos)
            return true;
        else
            return false;
    };

    // Check if file is MC or not
    bool mc_file = isMC(evt_parser->GetEvtTree()->GetName());
    
    // Parse local config file
    std::shared_ptr<energy_config> config = std::make_shared<energy_config>(input_args.energy_config_file);
    
    // Extract energy binning from config file
    auto energy_binning = config->GetEnergyBinning();
    auto energy_nbins = (int)energy_binning.size() - 1;
    
    // Load BDT config class
    std::shared_ptr<bdt> bdt_evt = std::make_shared<bdt>(
		input_args.bdt_config_file,
		input_args.bdt_learning_method,
		input_args.cosine_regularize_path,
		input_args.box_cox_regularize_path,
		input_args.verbose);
    
    // Extract the efficiency correction functions
    auto corrections = extractCorrectionFunctions(input_args.eff_corr_function.c_str());

    auto shift_function = std::get<0>(std::get<0>(corrections));
    auto sigma_function = std::get<1>(std::get<0>(corrections));
    auto gr_shift = std::get<0>(std::get<1>(corrections));
    auto gr_sigma = std::get<1>(std::get<1>(corrections));

    /*
    TMVA::Tools::Instance();
    check_bdt_learnign_method(input_args.bdt_learning_method);
    auto weights = get_config_info(parse_config_file(input_args.bdt_config_file));

    classifier_vars tmva_vars;

    std::shared_ptr<TMVA::Reader> LE_reader = std::make_shared<TMVA::Reader>();
    std::shared_ptr<TMVA::Reader> ME_reader = std::make_shared<TMVA::Reader>();
    std::shared_ptr<TMVA::Reader> HE_reader = std::make_shared<TMVA::Reader>();

    link_reader_vars(LE_reader, tmva_vars);
    link_reader_vars(ME_reader, tmva_vars);
    link_reader_vars(HE_reader, tmva_vars);

    bookMVA(LE_reader, weights.le_weights, input_args.bdt_learning_method);
    bookMVA(ME_reader, weights.me_weights, input_args.bdt_learning_method);
    bookMVA(HE_reader, weights.he_weights, input_args.bdt_learning_method);

    std::shared_ptr<vars> trans_vars = std::make_shared<vars>(
        input_args.cosine_regularize_path, 
        input_args.box_cox_regularize_path, 
        input_args.verbose);
    */

    // Get chain entries
    const double _entries = evt_parser->GetEvtTree()->GetEntries();
    if (input_args.verbose)
    {
        config->PrintActiveFilters();
        std::cout << "Total number of events: " << _entries;
    }
    // Initialize RDF
    ROOT::EnableImplicitMT(input_args.threads);
    ROOT::RDataFrame fr(*(evt_parser->GetEvtTree()));

    // Compute XTRL
    auto compute_xtrl = [](const double last_layer_energy_fraction, const double sumrms) -> double {
        return last_layer_energy_fraction != -999 ? 0.125e-6 * pow(sumrms, 4) * last_layer_energy_fraction : -999;
    };

    // Compute BDT
    /*
    auto compute_bdt = [&trans_vars, &tmva_vars, &energy_binning, &input_args, LE_reader, ME_reader, HE_reader](
        const std::vector<double> rms_layer, 
        const double sumrms, 
        const std::vector<double> energy_frac_layer,
        const double energy_fraction_last_layer,
        const double energy_corr_gev,
        const TVector3 bgo_trajectory) -> double 
        {
            auto transformed_vars = 
            trans_vars->GetVars(
                rms_layer,
                sumrms,
                energy_frac_layer,
                energy_fraction_last_layer,
                energy_corr_gev,
                energy_binning,
                bgo_trajectory);

            double tmva_classifier {-999};

            sync_vars(transformed_vars, tmva_vars);
            
            if (transformed_vars.corrected_energy_gev>=10 && transformed_vars.corrected_energy_gev<100) {
                tmva_classifier = LE_reader->EvaluateMVA(input_args.bdt_learning_method.c_str());
            }
            else if (transformed_vars.corrected_energy_gev>=100 && transformed_vars.corrected_energy_gev<1000) {
                tmva_classifier = ME_reader->EvaluateMVA(input_args.bdt_learning_method.c_str());
            }
            else if (transformed_vars.corrected_energy_gev>=1000 && transformed_vars.corrected_energy_gev<=10000) {
                tmva_classifier = HE_reader->EvaluateMVA(input_args.bdt_learning_method.c_str());
            }

            return tmva_classifier;
        };
    */
    auto compute_bdt = [bdt_evt, energy_binning](
        const std::vector<double> rms_layer, 
        const double sumrms, 
        const std::vector<double> energy_frac_layer,
        const double energy_fraction_last_layer,
        const double energy_corr_gev,
        const TVector3 bgo_trajectory) -> double 
        {
            return bdt_evt->ComputeMVA(
                rms_layer,
                sumrms,
                energy_frac_layer,
                energy_fraction_last_layer,
                energy_corr_gev,
                energy_binning,
                bgo_trajectory);
        };    

    // Compute XTRL cut
    auto xtrl_loose_cut = [](const double bgo_total_energy_gev, const double xtrl) -> bool {
        auto xtrl_cut_value {3*log10(bgo_total_energy_gev/200) + 12};
        if (xtrl<xtrl_cut_value && xtrl!=-999)
            return true;
        else
            return false;
    };

    auto xtrl_loose_reject_cut = [](const double bgo_total_energy_gev, const double xtrl) -> bool {
        auto xtrl_cut_value {3*log10(bgo_total_energy_gev/200) + 12};
        if (xtrl>xtrl_cut_value && xtrl!=-999)
            return true;
        else
            return false;
    };

    // Compute the energy binning giving the corrected energy
    auto get_energy_bin = [&energy_binning] (const double corrected_energy) -> int
    {
        int energybin {1};
        for (size_t idx=0; idx<energy_binning.size()-1; ++idx)
            if (corrected_energy>=energy_binning[idx] && corrected_energy<energy_binning[idx+1])
                energybin = idx+1;
        return energybin;
    };

    // Compute BDT cut
    auto get_bdt_cut = [bdt_evt, shift_function, sigma_function, gr_shift, gr_sigma] (const double energy_gev, const int energy_bin) -> double {
        /*
        double cut {0};
        double cut_corr {0};

        if (energy_gev>=10 && energy_gev<100)
            cut = bdt_evt->GetLowEnergyBDTCut();
        else if (energy_gev>=100 && energy_gev<1000)
            cut = energy_bin != 32 ? bdt_evt->GetMidEnergyBDTCut() : bdt_evt->GetMidEnergyBDTCut() - 0.07;
        else if (energy_gev>=1000 && energy_gev<=10000)     
            cut = bdt_evt->GetHighEnergyBDTCut();
        
        if (energy_bin<=33)
            cut_corr = (cut - gr_shift->GetPointY(energy_bin-1))/gr_sigma->GetPointY(energy_bin-1); 
        else
            cut_corr = (cut - shift_function->Eval(energy_gev))/sigma_function->Eval(energy_gev);
        return cut_corr;
        */

        double cut {0};
        if (energy_gev>=10 && energy_gev<100)
            cut = bdt_evt->GetBDTCut_10_100();
        else if (energy_gev>=100 && energy_gev<250)
            cut = bdt_evt->GetBDTCut_100_250();
        else if (energy_gev>=250 && energy_gev<500)     
            cut = bdt_evt->GetBDTCut_250_500();
        else if (energy_gev>=500 && energy_gev<1000)     
            cut = bdt_evt->GetBDTCut_500_1000();
        else if (energy_gev>=1000 && energy_gev<3000)     
            cut = bdt_evt->GetBDTCut_1000_3000();
        else if (energy_gev>=3000)     
            cut = bdt_evt->GetBDTCut_3000();

        return cut;
    };

    if (input_args.verbose)
        std::cout << "\n\nAnalysis running\n\n";
    
    auto get_let_prescale = [](const double geographic_latitude) -> double {
        return geographic_latitude>=-20 && geographic_latitude<=20 ? 1/8. : 1/64.;
    };
    
    auto get_unb_prescale = [](const double geographic_latitude) -> double {
        return geographic_latitude>=-20 && geographic_latitude<=20 ? 1/512. : 1/2048.;
    };

    // Reverse BDT cut
    const bool reverse_bdt_cut {false};

    // Trigger histos
    auto h_trigger_efficiency_accepted_het_over_let_tight_xtrl = mc_file ? 
                                            fr.Filter("trigger_efficiency_preselection==1 && trigger_efficiency_preselection_is_het==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_trigger_efficiency_accepted_het_over_let_tight_xtrl", "HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev")

                                            :

                                            fr.Filter("trigger_efficiency_preselection==1 && trigger_efficiency_preselection_is_het==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Define("trigger_w", get_let_prescale, {"geo_lat"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_over_let_tight_xtrl", "HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev", "trigger_w");

    auto h_trigger_efficiency_accepted_het_over_unb_tight_xtrl = mc_file ?
                                            fr.Filter("trigger_efficiency_preselection==1 && trigger_efficiency_preselection_is_het==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_trigger_efficiency_accepted_het_over_unb_tight_xtrl", "HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev")

                                            :

                                            fr.Filter("trigger_efficiency_preselection==1 && trigger_efficiency_preselection_is_het==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Define("trigger_w", get_unb_prescale, {"geo_lat"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_over_unb_tight_xtrl", "HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev", "trigger_w");

    auto h_trigger_efficiency_accepted_het_let_tight_xtrl = mc_file ? 
                                            fr.Filter("trigger_efficiency_preselection==1 && (trigger_efficiency_preselection_is_het==1 || trigger_efficiency_preselection_is_let==1)")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_trigger_efficiency_accepted_het_let_tight_xtrl", "LET + HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev")

                                            :

                                            fr.Filter("trigger_efficiency_preselection==1 && (trigger_efficiency_preselection_is_het==1 || trigger_efficiency_preselection_is_let==1)")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Define("trigger_w", [&get_let_prescale](const bool is_het, const double geo_lat) -> double { if (is_het) return get_let_prescale(geo_lat); else return 1;}, {"trigger_efficiency_preselection_is_het", "geo_lat"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_trigger_efficiency_accepted_het_let_tight_xtrl", "LET + HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev", "trigger_w");


    auto h_trigger_efficiency_accepted_het_unb_tight_xtrl = mc_file ? 
                                            fr.Filter("trigger_efficiency_preselection==1 && (trigger_efficiency_preselection_is_het==1 || trigger_efficiency_preselection_is_unb==1)")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_trigger_efficiency_accepted_het_unb_tight_xtrl", "UNB + HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev")

                                            :

                                            fr.Filter("trigger_efficiency_preselection==1 && (trigger_efficiency_preselection_is_het==1 || trigger_efficiency_preselection_is_unb==1)")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Define("trigger_w", [&get_unb_prescale](const bool is_het, const double geo_lat) -> double { if (is_het) return get_unb_prescale(geo_lat); else return 1;}, {"trigger_efficiency_preselection_is_het", "geo_lat"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_trigger_efficiency_accepted_het_unb_tight_xtrl", "UNB + HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev", "trigger_w");
    
    auto h_trigger_efficiency_accepted_het_over_let_loose_xtrl = mc_file ? 
                                            fr.Filter("trigger_efficiency_preselection==1 && trigger_efficiency_preselection_is_het==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_over_let_loose_xtrl", "HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev")

                                            :

                                            fr.Filter("trigger_efficiency_preselection==1 && trigger_efficiency_preselection_is_het==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Define("trigger_w", get_let_prescale, {"geo_lat"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_over_let_loose_xtrl", "HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev", "trigger_w");

    auto h_trigger_efficiency_accepted_het_over_unb_loose_xtrl = mc_file ? 
                                            fr.Filter("trigger_efficiency_preselection==1 && trigger_efficiency_preselection_is_het==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_over_unb_loose_xtrl", "HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev")

                                            :

                                            fr.Filter("trigger_efficiency_preselection==1 && trigger_efficiency_preselection_is_het==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Define("trigger_w", get_unb_prescale, {"geo_lat"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_over_unb_loose_xtrl", "HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev", "trigger_w");


    auto h_trigger_efficiency_accepted_het_let_loose_xtrl = mc_file ? 
                                            fr.Filter("trigger_efficiency_preselection==1 && (trigger_efficiency_preselection_is_het==1 || trigger_efficiency_preselection_is_let==1)")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_let_loose_xtrl", "LET + HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev")

                                            :

                                            fr.Filter("trigger_efficiency_preselection==1 && (trigger_efficiency_preselection_is_het==1 || trigger_efficiency_preselection_is_let==1)")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Define("trigger_w", [&get_let_prescale](const bool is_het, const double geo_lat) -> double { if (is_het) return get_let_prescale(geo_lat); else return 1;}, {"trigger_efficiency_preselection_is_het", "geo_lat"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_let_loose_xtrl", "LET + HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev", "trigger_w");

    auto h_trigger_efficiency_accepted_het_unb_loose_xtrl = mc_file ? 
                                            fr.Filter("trigger_efficiency_preselection==1 && (trigger_efficiency_preselection_is_het==1 || trigger_efficiency_preselection_is_unb==1)")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_unb_loose_xtrl", "LET + HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev")

                                            :

                                            fr.Filter("trigger_efficiency_preselection==1 && (trigger_efficiency_preselection_is_het==1 || trigger_efficiency_preselection_is_unb==1)")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Define("trigger_w", [&get_unb_prescale](const bool is_het, const double geo_lat) -> double { if (is_het) return get_unb_prescale(geo_lat); else return 1;}, {"trigger_efficiency_preselection_is_het", "geo_lat"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_unb_loose_xtrl", "LET + HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev", "trigger_w");

    auto h_trigger_efficiency_accepted_het_over_let_bdt = mc_file ? 
                                            fr.Filter("trigger_efficiency_preselection==1 && trigger_efficiency_preselection_is_het==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_over_let_bdt", "HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev")

                                            :

                                            fr.Filter("trigger_efficiency_preselection==1 && trigger_efficiency_preselection_is_het==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Define("trigger_w", get_let_prescale, {"geo_lat"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_over_let_bdt", "HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev", "trigger_w");

     auto h_trigger_efficiency_accepted_het_over_unb_bdt = mc_file ? 
                                            fr.Filter("trigger_efficiency_preselection==1 && trigger_efficiency_preselection_is_het==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_over_unb_bdt", "HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev")

                                            :

                                            fr.Filter("trigger_efficiency_preselection==1 && trigger_efficiency_preselection_is_het==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Define("trigger_w", get_unb_prescale, {"geo_lat"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_over_unb_bdt", "HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev", "trigger_w");

    auto h_trigger_efficiency_accepted_het_let_bdt = mc_file ?
                                            fr.Filter("trigger_efficiency_preselection==1 && (trigger_efficiency_preselection_is_het==1 || trigger_efficiency_preselection_is_let==1)")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_let_bdt", "LET + HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev")

                                            :

                                            fr.Filter("trigger_efficiency_preselection==1 && (trigger_efficiency_preselection_is_het==1 || trigger_efficiency_preselection_is_let==1)")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("trigger_w", [&get_let_prescale](const bool is_het, const double geo_lat) -> double { if (is_het) return get_let_prescale(geo_lat); else return 1;}, {"trigger_efficiency_preselection_is_het", "geo_lat"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_let_bdt", "LET + HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev", "trigger_w");

    auto h_trigger_efficiency_accepted_het_unb_bdt = mc_file ? 
                                            fr.Filter("trigger_efficiency_preselection==1 && (trigger_efficiency_preselection_is_het==1 || trigger_efficiency_preselection_is_unb==1)")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_unb_bdt", "LET + HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev")

                                            :

                                            fr.Filter("trigger_efficiency_preselection==1 && (trigger_efficiency_preselection_is_het==1 || trigger_efficiency_preselection_is_unb==1)")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("trigger_w", [&get_unb_prescale](const bool is_het, const double geo_lat) -> double { if (is_het) return get_unb_prescale(geo_lat); else return 1;}, {"trigger_efficiency_preselection_is_het", "geo_lat"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_unb_bdt", "LET + HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev", "trigger_w");

    // MaxRMS histos
    auto h_maxrms_efficiency_accepted_tight_xtrl = fr.Filter("maxrms_efficiency_preselection==1 && maxrms_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_maxrms_efficiency_accepted_tight_xtrl", "HET Trigger + MaxRMS Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_maxrms_efficiency_total_tight_xtrl = fr.Filter("maxrms_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_maxrms_efficiency_total_tight_xtrl", "HET Trigger + MaxRMS", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
    
    auto h_maxrms_efficiency_accepted_loose_xtrl = fr.Filter("maxrms_efficiency_preselection==1 && maxrms_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_maxrms_efficiency_accepted_loose_xtrl", "HET Trigger + MaxRMS Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_maxrms_efficiency_total_loose_xtrl = fr.Filter("maxrms_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_maxrms_efficiency_total_loose_xtrl", "HET Trigger + MaxRMS", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_maxrms_efficiency_accepted_bdt = fr.Filter("maxrms_efficiency_preselection==1 && maxrms_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_maxrms_efficiency_accepted_bdt", "HET Trigger + MaxRMS Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_maxrms_efficiency_total_bdt = fr.Filter("maxrms_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_maxrms_efficiency_total_bdt", "HET Trigger + MaxRMS", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    // nbarlayer13 histos
    auto h_nbarlayer13_efficiency_accepted_tight_xtrl = fr.Filter("nbarlayer13_efficiency_preselection==1 && nbarlayer13_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_nbarlayer13_efficiency_accepted_tight_xtrl", "HET Trigger + nbarlayer13 Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_nbarlayer13_efficiency_total_tight_xtrl = fr.Filter("nbarlayer13_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_nbarlayer13_efficiency_total_tight_xtrl", "HET Trigger + nbarlayer13", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_nbarlayer13_efficiency_accepted_loose_xtrl = fr.Filter("nbarlayer13_efficiency_preselection==1 && nbarlayer13_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_nbarlayer13_efficiency_accepted_loose_xtrl", "HET Trigger + nbarlayer13 Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_nbarlayer13_efficiency_total_loose_xtrl = fr.Filter("nbarlayer13_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_nbarlayer13_efficiency_total_loose_xtrl", "HET Trigger + nbarlayer13", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_nbarlayer13_efficiency_accepted_bdt = fr.Filter("nbarlayer13_efficiency_preselection==1 && nbarlayer13_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_nbarlayer13_efficiency_accepted_bdt", "HET Trigger + nbarlayer13 Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_nbarlayer13_efficiency_total_bdt = fr.Filter("nbarlayer13_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_nbarlayer13_efficiency_total_bdt", "HET Trigger + nbarlayer13", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    // MaxRMS and nbarlayer13 histos
    auto h_maxrms_and_nbarlayer13_efficiency_accepted_tight_xtrl = fr.Filter("maxrms_and_nbarlayer13_efficiency_preselection==1 && maxrms_and_nbarlayer13_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_maxrms_and_nbarlayer13_efficiency_accepted_tight_xtrl", "HET Trigger + maxrms & nbarlayer13 Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
    auto h_maxrms_and_nbarlayer13_efficiency_total_tight_xtrl = fr.Filter("maxrms_and_nbarlayer13_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_maxrms_and_nbarlayer13_efficiency_total_tight_xtrl", "HET Trigger + maxrms & nbarlayer13", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_maxrms_and_nbarlayer13_efficiency_accepted_loose_xtrl = fr.Filter("maxrms_and_nbarlayer13_efficiency_preselection==1 && maxrms_and_nbarlayer13_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_maxrms_and_nbarlayer13_efficiency_accepted_loose_xtrl", "HET Trigger + maxrms & nbarlayer13 Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_maxrms_and_nbarlayer13_efficiency_total_loose_xtrl = fr.Filter("maxrms_and_nbarlayer13_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_maxrms_and_nbarlayer13_efficiency_total_loose_xtrl", "HET Trigger + maxrms & nbarlayer13", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_maxrms_and_nbarlayer13_efficiency_accepted_bdt = fr.Filter("maxrms_and_nbarlayer13_efficiency_preselection==1 && maxrms_and_nbarlayer13_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_maxrms_and_nbarlayer13_efficiency_accepted_bdt", "HET Trigger + maxrms & nbarlayer13 Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_maxrms_and_nbarlayer13_efficiency_total_bdt = fr.Filter("maxrms_and_nbarlayer13_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_maxrms_and_nbarlayer13_efficiency_total_bdt", "HET Trigger + maxrms & nbarlayer13", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    // sumRMS low energy histos
    auto h_sumrms_low_energy_efficiency_accepted_tight_xtrl = fr.Filter("sumrms_low_energy_cut_efficiency_preselection==1 && sumrms_low_energy_cut_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_sumrms_low_energy_efficiency_accepted_tight_xtrl", "HET Trigger + sumrms low energy cut Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
    auto h_sumrms_low_energy_efficiency_total_tight_xtrl = fr.Filter("sumrms_low_energy_cut_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_sumrms_low_energy_efficiency_total_tight_xtrl", "HET Trigger + sumrms low energy cut", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_sumrms_low_energy_efficiency_accepted_loose_xtrl = fr.Filter("sumrms_low_energy_cut_efficiency_preselection==1 && sumrms_low_energy_cut_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_sumrms_low_energy_efficiency_accepted_loose_xtrl", "HET Trigger + sumrms low energy cut Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_sumrms_low_energy_efficiency_total_loose_xtrl = fr.Filter("sumrms_low_energy_cut_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_sumrms_low_energy_efficiency_total_loose_xtrl", "HET Trigger + sumrms low energy cut", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_sumrms_low_energy_efficiency_accepted_bdt = fr.Filter("sumrms_low_energy_cut_efficiency_preselection==1 && sumrms_low_energy_cut_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_sumrms_low_energy_efficiency_accepted_bdt", "HET Trigger + sumrms low energy cut Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_sumrms_low_energy_efficiency_total_bdt = fr.Filter("sumrms_low_energy_cut_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_sumrms_low_energy_efficiency_total_bdt", "HET Trigger + sumrms low energy cut", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    // Track Selection histos
    auto h_track_efficiency_accepted_tight_xtrl = fr.Filter("track_efficiency_preselection==1 && track_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_track_efficiency_accepted_tight_xtrl", "HET Trigger + Track Selection Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_track_efficiency_total_tight_xtrl = fr.Filter("track_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_track_efficiency_total_tight_xtrl", "HET Trigger + Track Selection", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
    
    auto h_track_efficiency_accepted_loose_xtrl = fr.Filter("track_efficiency_preselection==1 && track_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_track_efficiency_accepted_loose_xtrl", "HET Trigger + Track Selection Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_track_efficiency_total_loose_xtrl = fr.Filter("track_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_track_efficiency_total_loose_xtrl", "HET Trigger + Track Selection", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_track_efficiency_accepted_bdt = fr.Filter("track_efficiency_preselection==1 && track_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_track_efficiency_accepted_bdt", "HET Trigger + Track Selection Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_track_efficiency_total_bdt = fr.Filter("track_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_track_efficiency_total_bdt", "HET Trigger + Track Selection", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    // STK 1 RM cut
    auto h_stk_1rm_efficiency_accepted_tight_xtrl = fr.Filter("stk_1rm_cut_efficiency_preselection==1 && stk_1rm_cut_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_stk_1rm_efficiency_accepted_tight_xtrl", "HET Trigger + STK 1 RM Selection Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_stk_1rm_efficiency_total_tight_xtrl = fr.Filter("stk_1rm_cut_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_stk_1rm_efficiency_total_tight_xtrl", "HET Trigger + STK 1 RM Selection", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
    
    auto h_stk_1rm_efficiency_accepted_loose_xtrl = fr.Filter("stk_1rm_cut_efficiency_preselection==1 && stk_1rm_cut_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_stk_1rm_efficiency_accepted_loose_xtrl", "HET Trigger + STK 1 RM Selection Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_stk_1rm_efficiency_total_loose_xtrl = fr.Filter("stk_1rm_cut_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_stk_1rm_efficiency_total_loose_xtrl", "HET Trigger + STK 1 RM Selection", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_stk_1rm_efficiency_accepted_bdt = fr.Filter("stk_1rm_cut_efficiency_preselection==1 && stk_1rm_cut_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_stk_1rm_efficiency_accepted_bdt", "HET Trigger + STK 1 RM Selection Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_stk_1rm_efficiency_total_bdt = fr.Filter("stk_1rm_cut_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_stk_1rm_efficiency_total_bdt", "HET Trigger + STK 1 RM Selection", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    // Clusters on first STK layer
    auto h_clusters_on_first_STK_layer_within_psd_fvolume_accepted_tight_xtrl =  fr.Filter("track_efficiency_preselection==1 && track_efficiency_preselection_accepted==1")
                                                        .Filter("evtfilter_psd_fiducial_volume==1")
                                                        .Filter("STK_chargeX!=-999 && STK_chargeY!=-999")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                                        .Histo1D({"h_clusters_on_first_STK_layer_within_psd_fvolume_accepted_tight_xtrl", "Track Selection + PSD fvolume + STK clusters on first layer", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_clusters_on_first_STK_layer_within_psd_fvolume_total_tight_xtrl =  fr.Filter("track_efficiency_preselection==1 && track_efficiency_preselection_accepted==1")
                                                        .Filter("evtfilter_psd_fiducial_volume==1")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                                        .Histo1D({"h_clusters_on_first_STK_layer_within_psd_fvolume_total_tight_xtrl", "Track Selection + PSD fvolume", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_clusters_on_first_STK_layer_outside_psd_fvolume_accepted_tight_xtrl =  fr.Filter("track_efficiency_preselection==1 && track_efficiency_preselection_accepted==1")
                                                        .Filter("evtfilter_psd_fiducial_volume==0")
                                                        .Filter("STK_chargeX!=-999 && STK_chargeY!=-999")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                                        .Histo1D({"h_clusters_on_first_STK_layer_outside_psd_fvolume_accepted_tight_xtrl", "Track Selection + PSD fvolume + STK clusters on first layer", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_clusters_on_first_STK_layer_outside_psd_fvolume_total_tight_xtrl =  fr.Filter("track_efficiency_preselection==1 && track_efficiency_preselection_accepted==1")
                                                        .Filter("evtfilter_psd_fiducial_volume==0")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                                        .Histo1D({"h_clusters_on_first_STK_layer_outside_psd_fvolume_total_tight_xtrl", "Track Selection + PSD fvolume", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_clusters_on_first_STK_layer_within_psd_fvolume_accepted_loose_xtrl =  fr.Filter("track_efficiency_preselection==1 && track_efficiency_preselection_accepted==1")
                                                        .Filter("evtfilter_psd_fiducial_volume==1")
                                                        .Filter("STK_chargeX!=-999 && STK_chargeY!=-999")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Define("energy_gev", "energy * 0.001")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                                        .Histo1D({"h_clusters_on_first_STK_layer_within_psd_fvolume_accepted_loose_xtrl", "Track Selection + PSD fvolume + STK clusters on first layer", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_clusters_on_first_STK_layer_within_psd_fvolume_total_loose_xtrl =  fr.Filter("track_efficiency_preselection==1 && track_efficiency_preselection_accepted==1")
                                                        .Filter("evtfilter_psd_fiducial_volume==1")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Define("energy_gev", "energy * 0.001")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                                        .Histo1D({"h_clusters_on_first_STK_layer_within_psd_fvolume_total_loose_xtrl", "Track Selection + PSD fvolume", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_clusters_on_first_STK_layer_outside_psd_fvolume_accepted_loose_xtrl =  fr.Filter("track_efficiency_preselection==1 && track_efficiency_preselection_accepted==1")
                                                        .Filter("evtfilter_psd_fiducial_volume==0")
                                                        .Filter("STK_chargeX!=-999 && STK_chargeY!=-999")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Define("energy_gev", "energy * 0.001")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                                        .Histo1D({"h_clusters_on_first_STK_layer_outside_psd_fvolume_accepted_loose_xtrl", "Track Selection + PSD fvolume + STK clusters on first layer", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_clusters_on_first_STK_layer_outside_psd_fvolume_total_loose_xtrl =  fr.Filter("track_efficiency_preselection==1 && track_efficiency_preselection_accepted==1")
                                                        .Filter("evtfilter_psd_fiducial_volume==0")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Define("energy_gev", "energy * 0.001")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                                        .Histo1D({"h_clusters_on_first_STK_layer_outside_psd_fvolume_total_loose_xtrl", "Track Selection + PSD fvolume", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_clusters_on_first_STK_layer_within_psd_fvolume_accepted_bdt =  fr.Filter("track_efficiency_preselection==1 && track_efficiency_preselection_accepted==1")
                                                        .Filter("evtfilter_psd_fiducial_volume==1")
                                                        .Filter("STK_chargeX!=-999 && STK_chargeY!=-999")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                                        .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                                        .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                                        .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                                        .Histo1D({"h_clusters_on_first_STK_layer_within_psd_fvolume_accepted_bdt", "Track Selection + PSD fvolume + STK clusters on first layer", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_clusters_on_first_STK_layer_within_psd_fvolume_total_bdt =  fr.Filter("track_efficiency_preselection==1 && track_efficiency_preselection_accepted==1")
                                                        .Filter("evtfilter_psd_fiducial_volume==1")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                                        .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                                        .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                                        .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                                        .Histo1D({"h_clusters_on_first_STK_layer_within_psd_fvolume_total_bdt", "Track Selection + PSD fvolume", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_clusters_on_first_STK_layer_outside_psd_fvolume_accepted_bdt =  fr.Filter("track_efficiency_preselection==1 && track_efficiency_preselection_accepted==1")
                                                        .Filter("evtfilter_psd_fiducial_volume==0")
                                                        .Filter("STK_chargeX!=-999 && STK_chargeY!=-999")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                                        .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                                        .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                                        .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                                        .Histo1D({"h_clusters_on_first_STK_layer_outside_psd_fvolume_accepted_bdt", "Track Selection + PSD fvolume + STK clusters on first layer", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_clusters_on_first_STK_layer_outside_psd_fvolume_total_bdt =  fr.Filter("track_efficiency_preselection==1 && track_efficiency_preselection_accepted==1")
                                                        .Filter("evtfilter_psd_fiducial_volume==0")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                                        .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                                        .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                                        .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                                        .Histo1D({"h_clusters_on_first_STK_layer_outside_psd_fvolume_total_bdt", "Track Selection + PSD fvolume", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    // Track Selection histos within STK fiducial volume
    auto h_track_efficiency_stk_fvolume_accepted_tight_xtrl = fr.Filter("track_efficiency_preselection==1 && evtfilter_stk_fiducial_volume==1 && track_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_track_efficiency_stk_fvolume_accepted_tight_xtrl", "HET Trigger + Track Selection Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
    auto h_track_efficiency_stk_fvolume_total_tight_xtrl = fr.Filter("track_efficiency_preselection==1 && evtfilter_stk_fiducial_volume==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_track_efficiency_stk_fvolume_total_tight_xtrl", "HET Trigger + Track Selection", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
    
    auto h_track_efficiency_stk_fvolume_accepted_loose_xtrl = fr.Filter("track_efficiency_preselection==1 && evtfilter_stk_fiducial_volume==1 && track_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_track_efficiency_stk_fvolume_accepted_loose_xtrl", "HET Trigger + Track Selection Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
    auto h_track_efficiency_stk_fvolume_total_loose_xtrl = fr.Filter("track_efficiency_preselection==1 && evtfilter_stk_fiducial_volume==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_track_efficiency_stk_fvolume_total_loose_xtrl", "HET Trigger + Track Selection", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_track_efficiency_stk_fvolume_accepted_bdt = fr.Filter("track_efficiency_preselection==1 && evtfilter_stk_fiducial_volume==1 && track_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_track_efficiency_stk_fvolume_accepted_bdt", "HET Trigger + Track Selection Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_track_efficiency_stk_fvolume_total_bdt = fr.Filter("track_efficiency_preselection==1 && evtfilter_stk_fiducial_volume==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_track_efficiency_stk_fvolume_total_bdt", "HET Trigger + Track Selection", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    // PSD-STK match histos
    auto h_psdstkmatch_efficiency_accepted_tight_xtrl = fr.Filter("psdstkmatch_efficiency_preselection==1 && psdstkmatch_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_psdstkmatch_efficiency_accepted_tight_xtrl", "HET Trigger + PSD-STK Match Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
    auto h_psdstkmatch_efficiency_total_tight_xtrl = fr.Filter("psdstkmatch_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_psdstkmatch_efficiency_total_tight_xtrl", "HET Trigger + PSD-STK Match", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_psdstkmatch_efficiency_accepted_loose_xtrl = fr.Filter("psdstkmatch_efficiency_preselection==1 && psdstkmatch_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_psdstkmatch_efficiency_accepted_loose_xtrl", "HET Trigger + PSD-STK Match Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_psdstkmatch_efficiency_total_loose_xtrl = fr.Filter("psdstkmatch_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_psdstkmatch_efficiency_total_loose_xtrl", "HET Trigger + PSD-STK Match", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_psdstkmatch_efficiency_accepted_bdt = fr.Filter("psdstkmatch_efficiency_preselection==1 && psdstkmatch_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_psdstkmatch_efficiency_accepted_bdt", "HET Trigger + PSD-STK Match Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_psdstkmatch_efficiency_total_bdt = fr.Filter("psdstkmatch_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_psdstkmatch_efficiency_total_bdt", "HET Trigger + PSD-STK Match", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    // PSD charge histos
    auto h_psdcharge_efficiency_accepted_tight_xtrl = fr.Filter("psdcharge_efficiency_preselection==1 && psdcharge_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_psdcharge_efficiency_accepted_tight_xtrl", "HET Trigger + PSD Charge Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
    auto h_psdcharge_efficiency_total_tight_xtrl = fr.Filter("psdcharge_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_psdcharge_efficiency_total_tight_xtrl", "HET Trigger + PSD Charge Match", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
    
    auto h_psdcharge_efficiency_accepted_loose_xtrl = fr.Filter("psdcharge_efficiency_preselection==1 && psdcharge_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_psdcharge_efficiency_accepted_loose_xtrl", "HET Trigger + PSD Charge Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_psdcharge_efficiency_total_loose_xtrl = fr.Filter("psdcharge_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_psdcharge_efficiency_total_loose_xtrl", "HET Trigger + PSD Charge Match", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_psdcharge_efficiency_accepted_bdt = fr.Filter("psdcharge_efficiency_preselection==1 && psdcharge_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_psdcharge_efficiency_accepted_bdt", "HET Trigger + PSD Charge Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_psdcharge_efficiency_total_bdt = fr.Filter("psdcharge_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_psdcharge_efficiency_total_bdt", "HET Trigger + PSD Charge Match", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    // STK charge histos
    auto h_stkcharge_efficiency_accepted_tight_xtrl = fr.Filter("stkcharge_efficiency_preselection==1 && stkcharge_efficiency_preselection_accepted==1 && STK_chargeX!=-999 && STK_chargeY!=-999")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_stkcharge_efficiency_accepted_tight_xtrl", "HET Trigger + STK Charge Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
    auto h_stkcharge_efficiency_total_tight_xtrl = fr.Filter("stkcharge_efficiency_preselection==1 && STK_chargeX!=-999 && STK_chargeY!=-999")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_stkcharge_efficiency_total_tight_xtrl", "HET Trigger + STK Charge Match", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
    
    auto h_stkcharge_efficiency_accepted_loose_xtrl = fr.Filter("stkcharge_efficiency_preselection==1 && stkcharge_efficiency_preselection_accepted==1 && STK_chargeX!=-999 && STK_chargeY!=-999")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_stkcharge_efficiency_accepted_loose_xtrl", "HET Trigger + STK Charge Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_stkcharge_efficiency_total_loose_xtrl = fr.Filter("stkcharge_efficiency_preselection==1 && STK_chargeX!=-999 && STK_chargeY!=-999")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_stkcharge_efficiency_total_loose_xtrl", "HET Trigger + STK Charge", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_stkcharge_efficiency_accepted_bdt = fr.Filter("stkcharge_efficiency_preselection==1 && stkcharge_efficiency_preselection_accepted==1 && STK_chargeX!=-999 && STK_chargeY!=-999")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_stkcharge_efficiency_accepted_bdt", "HET Trigger + STK Charge Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_stkcharge_efficiency_total_bdt = fr.Filter("stkcharge_efficiency_preselection==1 && STK_chargeX!=-999 && STK_chargeY!=-999")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                            .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                            .Filter([reverse_bdt_cut] (const double tmva_value, const double tmva_cut) {return !reverse_bdt_cut ? tmva_value > tmva_cut : tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                            .Histo1D({"h_stkcharge_efficiency_total_bdt", "HET Trigger + STK Charge", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    // XTRL vs STK cosine histo
    auto h_xtrl_stk_cosine = fr.Filter("HET_trigger==1 && evtfilter_correct_bgo_reco==1")
                            .Define("corr_energy_gev", "energy_corr * 0.001")
                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                            .Histo2D({"h_xtrl_stk_cosine", "XTRL vs STK direction; cosine STK direction #cos(#theta); xtrl", 100, 0, 1, 300, 0, 300}, "STK_bestTrack_costheta", "xtrl_evt");

    auto h_xtrl_stk_cosine_zoom = fr.Filter("HET_trigger==1 && evtfilter_correct_bgo_reco==1")
                            .Define("corr_energy_gev", "energy_corr * 0.001")
                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                            .Histo2D({"h_xtrl_stk_cosine_zoom", "XTRL vs STK direction; cosine STK direction #cos(#theta); xtrl", 100, 0, 1, 100, 0, 10}, "STK_bestTrack_costheta", "xtrl_evt");
    
    auto h_xtrl_bgo_cosine = fr.Filter("HET_trigger==1 && evtfilter_correct_bgo_reco==1")
                            .Define("corr_energy_gev", "energy_corr * 0.001")
                            .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                            .Histo2D({"h_xtrl_bgo_cosine", "XTRL vs BGO direction; cosine BGO direction #cos(#theta); xtrl", 100, 0, 1, 300, 0, 300}, "bgorec_cosine", "xtrl_evt");

    auto h_xtrl_bgo_cosine_zoom = fr.Filter("HET_trigger==1 && evtfilter_correct_bgo_reco==1")
                            .Define("corr_energy_gev", "energy_corr * 0.001")
                            .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                            .Histo2D({"h_xtrl_bgo_cosine_zoom", "XTRL vs BGO direction; cosine BGO direction #cos(#theta); xtrl", 100, 0, 1, 100, 0, 10}, "bgorec_cosine", "xtrl_evt");

    // Basic PSD and STK charge plots
    auto h_psd_charge_after_stk_cut = fr.Filter("psdcharge_efficiency_preselection==1")
                            .Histo1D({"h_psd_charge_after_stk_cut", "PSD charge after STK cut; PSD Charge; entries", 100, 0, 40}, "PSD_charge");
    auto h_stk_charge_after_stk_cut = fr.Filter("psdcharge_efficiency_preselection==1")
                            .Histo1D({"h_stk_charge_after_stk_cut", "STK charge after STK cut; STK Charge; entries", 100, 0, 40}, "STK_charge");
    auto h_stk_charge = fr.Filter("evtfilter_stk_charge_measurement==1")
                            .Histo1D({"h_stk_charge", "STK charge; STK Charge; entries", 100, 0, 40}, "STK_charge");
    auto h_stk_charge_after_psd_charge_cut = fr.Filter("evtfilter_all_cut==1")
                            .Histo1D({"h_stk_charge_after_psd_charge_cut", "STK charge; STK Charge; entries", 100, 0, 40}, "STK_charge");

    // PSD fiducial volume
    auto n_events_presel_psd_fvolume    = *(fr.Filter("evtfilter_track_selection_cut==true").Count());
    auto n_events_within_psd_fvolume    = *(fr.Filter("evtfilter_track_selection_cut==true").Filter("evtfilter_psd_fiducial_volume==true").Count());
    auto n_events_outside_psd_fvolume   = *(fr.Filter("evtfilter_track_selection_cut==true").Filter("evtfilter_psd_fiducial_volume==false").Count());

    std::cout << "\n*** PSD Fiducial Volume ***\n\n";
    std::cout << "Number of events after BGO and STK selection: " << n_events_presel_psd_fvolume;
    std::cout << "\nNumber of events within the PSD fiducial volume: " << n_events_within_psd_fvolume;
    std::cout << "\nNumber of events outside the PSD fiducial volume: " << n_events_outside_psd_fvolume;
    std::cout << "\n\n***************************\n\n";

    // PSD fiducial volume distance plots
    auto h_psd_stk_match_distance_x_20_100 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                .Histo1D({"h_psd_stk_match_distance_x_20_100", "PSD/STK match - X distance", 200, -300, 300}, "SPD_STK_match_X_distance");

    auto h_psd_stk_match_distance_x_100_250 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                .Histo1D({"h_psd_stk_match_distance_x_100_250", "PSD/STK match - X distance", 200, -300, 300}, "SPD_STK_match_X_distance");

    auto h_psd_stk_match_distance_x_250_500 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                .Histo1D({"h_psd_stk_match_distance_x_250_500", "PSD/STK match - X distance",  200, -300, 300}, "SPD_STK_match_X_distance");

    auto h_psd_stk_match_distance_x_500_1000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                .Histo1D({"h_psd_stk_match_distance_x_500_1000", "PSD/STK match - X distance",  200, -300, 300}, "SPD_STK_match_X_distance");
    
    auto h_psd_stk_match_distance_x_1000_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                .Histo1D({"h_psd_stk_match_distance_x_1000_3000", "PSD/STK match - X distance", 200, -300, 300}, "SPD_STK_match_X_distance");

    auto h_psd_stk_match_distance_x_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=3000")
                                                .Histo1D({"h_psd_stk_match_distance_x_3000", "PSD/STK match - X distance", 200, -300, 300}, "SPD_STK_match_X_distance");

    auto h_psd_stk_match_distance_y_20_100 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                .Histo1D({"h_psd_stk_match_distance_y_20_100", "PSD/STK match - X distance", 200, -300, 300}, "SPD_STK_match_Y_distance");

    auto h_psd_stk_match_distance_y_100_250 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                .Histo1D({"h_psd_stk_match_distance_y_100_250", "PSD/STK match - Y distance", 200, -300, 300}, "SPD_STK_match_Y_distance");

    auto h_psd_stk_match_distance_y_250_500 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                .Histo1D({"h_psd_stk_match_distance_y_250_500", "PSD/STK match - Y distance",  200, -300, 300}, "SPD_STK_match_Y_distance");

    auto h_psd_stk_match_distance_y_500_1000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                .Histo1D({"h_psd_stk_match_distance_y_500_1000", "PSD/STK match - Y distance",  200, -300, 300}, "SPD_STK_match_Y_distance");
    
    auto h_psd_stk_match_distance_y_1000_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                .Histo1D({"h_psd_stk_match_distance_y_1000_3000", "PSD/STK match - Y distance", 200, -300, 300}, "SPD_STK_match_Y_distance");

    auto h_psd_stk_match_distance_y_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=3000")
                                                .Histo1D({"h_psd_stk_match_distance_y_3000", "PSD/STK match - Y distance", 200, -300, 300}, "SPD_STK_match_Y_distance");

    // PSD fiducial volume distance plots within PSD fiducial volume
    auto h_psd_stk_match_distance_x_within_psd_fvolume_20_100 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                .Histo1D({"h_psd_stk_match_distance_x_within_psd_fvolume_20_100", "PSD/STK match - X distance", 200, -300, 300}, "SPD_STK_match_X_distance_fiducial_volume");

    auto h_psd_stk_match_distance_x_within_psd_fvolume_100_250 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                .Histo1D({"h_psd_stk_match_distance_x_within_psd_fvolume_100_250", "PSD/STK match - X distance", 200, -300, 300}, "SPD_STK_match_X_distance_fiducial_volume");

    auto h_psd_stk_match_distance_x_within_psd_fvolume_250_500 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                .Histo1D({"h_psd_stk_match_distance_x_within_psd_fvolume_250_500", "PSD/STK match - X distance",  200, -300, 300}, "SPD_STK_match_X_distance_fiducial_volume");

    auto h_psd_stk_match_distance_x_within_psd_fvolume_500_1000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                .Histo1D({"h_psd_stk_match_distance_x_within_psd_fvolume_500_1000", "PSD/STK match - X distance",  200, -300, 300}, "SPD_STK_match_X_distance_fiducial_volume");
    
    auto h_psd_stk_match_distance_x_within_psd_fvolume_1000_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                .Histo1D({"h_psd_stk_match_distance_x_within_psd_fvolume_1000_3000", "PSD/STK match - X distance", 200, -300, 300}, "SPD_STK_match_X_distance_fiducial_volume");

    auto h_psd_stk_match_distance_x_within_psd_fvolume_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=3000")
                                                .Histo1D({"h_psd_stk_match_distance_x_within_psd_fvolume_3000", "PSD/STK match - X distance", 200, -300, 300}, "SPD_STK_match_X_distance_fiducial_volume");

    auto h_psd_stk_match_distance_y_within_psd_fvolume_20_100 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                .Histo1D({"h_psd_stk_match_distance_y_within_psd_fvolume_20_100", "PSD/STK match - X distance", 200, -300, 300}, "SPD_STK_match_Y_distance_fiducial_volume");

    auto h_psd_stk_match_distance_y_within_psd_fvolume_100_250 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                .Histo1D({"h_psd_stk_match_distance_y_within_psd_fvolume_100_250", "PSD/STK match - Y distance", 200, -300, 300}, "SPD_STK_match_Y_distance_fiducial_volume");

    auto h_psd_stk_match_distance_y_within_psd_fvolume_250_500 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                .Histo1D({"h_psd_stk_match_distance_y_within_psd_fvolume_250_500", "PSD/STK match - Y distance",  200, -300, 300}, "SPD_STK_match_Y_distance_fiducial_volume");

    auto h_psd_stk_match_distance_y_within_psd_fvolume_500_1000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                .Histo1D({"h_psd_stk_match_distance_y_within_psd_fvolume_500_1000", "PSD/STK match - Y distance",  200, -300, 300}, "SPD_STK_match_Y_distance_fiducial_volume");
    
    auto h_psd_stk_match_distance_y_within_psd_fvolume_1000_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                .Histo1D({"h_psd_stk_match_distance_y_within_psd_fvolume_1000_3000", "PSD/STK match - Y distance", 200, -300, 300}, "SPD_STK_match_Y_distance_fiducial_volume");

    auto h_psd_stk_match_distance_y_within_psd_fvolume_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=3000")
                                                .Histo1D({"h_psd_stk_match_distance_y_within_psd_fvolume_3000", "PSD/STK match - Y distance", 200, -300, 300}, "SPD_STK_match_Y_distance_fiducial_volume");

    // PSD fiducial volume distance plots outside PSD fiducial volume
    auto h_psd_stk_match_distance_x_outside_psd_fvolume_20_100 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                .Histo1D({"h_psd_stk_match_distance_x_outside_psd_fvolume_20_100", "PSD/STK match - X distance", 200, -300, 300}, "SPD_STK_match_X_distance");

    auto h_psd_stk_match_distance_x_outside_psd_fvolume_100_250 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                .Histo1D({"h_psd_stk_match_distance_x_outside_psd_fvolume_100_250", "PSD/STK match - X distance", 200, -300, 300}, "SPD_STK_match_X_distance");

    auto h_psd_stk_match_distance_x_outside_psd_fvolume_250_500 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                .Histo1D({"h_psd_stk_match_distance_x_outside_psd_fvolume_250_500", "PSD/STK match - X distance",  200, -300, 300}, "SPD_STK_match_X_distance");

    auto h_psd_stk_match_distance_x_outside_psd_fvolume_500_1000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                .Histo1D({"h_psd_stk_match_distance_x_outside_psd_fvolume_500_1000", "PSD/STK match - X distance",  200, -300, 300}, "SPD_STK_match_X_distance");
    
    auto h_psd_stk_match_distance_x_outside_psd_fvolume_1000_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                .Histo1D({"h_psd_stk_match_distance_x_outside_psd_fvolume_1000_3000", "PSD/STK match - X distance", 200, -300, 300}, "SPD_STK_match_X_distance");

    auto h_psd_stk_match_distance_x_outside_psd_fvolume_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=3000")
                                                .Histo1D({"h_psd_stk_match_distance_x_outside_psd_fvolume_3000", "PSD/STK match - X distance", 200, -300, 300}, "SPD_STK_match_X_distance");

    auto h_psd_stk_match_distance_y_outside_psd_fvolume_20_100 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                .Histo1D({"h_psd_stk_match_distance_y_outside_psd_fvolume_20_100", "PSD/STK match - X distance", 200, -300, 300}, "SPD_STK_match_Y_distance");

    auto h_psd_stk_match_distance_y_outside_psd_fvolume_100_250 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                .Histo1D({"h_psd_stk_match_distance_y_outside_psd_fvolume_100_250", "PSD/STK match - Y distance", 200, -300, 300}, "SPD_STK_match_Y_distance");

    auto h_psd_stk_match_distance_y_outside_psd_fvolume_250_500 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                .Histo1D({"h_psd_stk_match_distance_y_outside_psd_fvolume_250_500", "PSD/STK match - Y distance",  200, -300, 300}, "SPD_STK_match_Y_distance");

    auto h_psd_stk_match_distance_y_outside_psd_fvolume_500_1000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                .Histo1D({"h_psd_stk_match_distance_y_outside_psd_fvolume_500_1000", "PSD/STK match - Y distance",  200, -300, 300}, "SPD_STK_match_Y_distance");
    
    auto h_psd_stk_match_distance_y_outside_psd_fvolume_1000_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                .Histo1D({"h_psd_stk_match_distance_y_outside_psd_fvolume_1000_3000", "PSD/STK match - Y distance", 200, -300, 300}, "SPD_STK_match_Y_distance");

    auto h_psd_stk_match_distance_y_outside_psd_fvolume_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=3000")
                                                .Histo1D({"h_psd_stk_match_distance_y_outside_psd_fvolume_3000", "PSD/STK match - Y distance", 200, -300, 300}, "SPD_STK_match_Y_distance");

    // STK charge plots within PSD fiducial volume and withoud PSD charge cut
    auto h_stk_charge_psd_fvolume_no_psd_cut_20_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_20_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_psd_fvolume_no_psd_cut_100_250 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_100_250", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_psd_fvolume_no_psd_cut_250_500 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_250_500", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_no_psd_cut_500_1000 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_500_1000", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_no_psd_cut_1000_3000 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_1000_3000", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_no_psd_cut_3000 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=3000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_3000", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    // STK charge plots within PSD fiducial volume and withoud PSD charge cut (xtrl signal cut)
    auto h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=3000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    // STK charge plots within PSD fiducial volume and withoud PSD charge cut (xtrl background cut)
    auto h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=3000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    // STK charge plots within PSD fiducial volume and with PSD charge cut
    auto h_stk_charge_psd_fvolume_psd_cut_20_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_20_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_psd_fvolume_psd_cut_100_250 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_100_250", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_psd_fvolume_psd_cut_250_500 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_250_500", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_psd_cut_500_1000 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_500_1000", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_psd_cut_1000_3000 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_1000_3000", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_psd_cut_3000 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=3000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_3000", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    // STK charge plots within PSD fiducial volume and with PSD charge cut (xtrl signal cut)
    auto h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=3000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    // STK charge plots within PSD fiducial volume and with PSD charge cut (xtrl background cut)
    auto h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=3000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    // STK charge plots outside PSD fiducial volume
    auto h_stk_charge_20_100 = fr.Filter("evtfilter_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                .Histo2D({"h_stk_charge_20_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_100_250 = fr.Filter("evtfilter_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                .Histo2D({"h_stk_charge_100_250", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_250_500 = fr.Filter("evtfilter_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                .Histo2D({"h_stk_charge_250_500", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_500_1000 = fr.Filter("evtfilter_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                .Histo2D({"h_stk_charge_500_1000", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_1000_3000 = fr.Filter("evtfilter_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                .Histo2D({"h_stk_charge_1000_3000", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_3000 = fr.Filter("evtfilter_psd_fiducial_volume==false")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Filter("corr_energy_gev>=3000")
                                            .Histo2D({"h_stk_charge_cut_3000", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    // STK charge plots outside PSD fiducial volume (xtrl signal cut)
    auto h_stk_charge_20_100_xtrl_12 = fr.Filter("evtfilter_psd_fiducial_volume==false")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                        .Histo2D({"h_stk_charge_20_100_xtrl_12", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_100_250_xtrl_12 = fr.Filter("evtfilter_psd_fiducial_volume==false")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                        .Histo2D({"h_stk_charge_100_250_xtrl_12", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_250_500_xtrl_12 = fr.Filter("evtfilter_psd_fiducial_volume==false")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                        .Histo2D({"h_stk_charge_250_500_xtrl_12", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_500_1000_xtrl_12 = fr.Filter("evtfilter_psd_fiducial_volume==false")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                        .Histo2D({"h_stk_charge_500_1000_xtrl_12", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_1000_3000_xtrl_12 = fr.Filter("evtfilter_psd_fiducial_volume==false")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                        .Histo2D({"h_stk_charge_1000_3000_xtrl_12", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_3000_xtrl_12 = fr.Filter("evtfilter_psd_fiducial_volume==false")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=3000")
                                                        .Histo2D({"h_stk_charge_3000_xtrl_12", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    // STK charge plots outside PSD fiducial volume and withoud PSD charge cut (xtrl background cut)
    auto h_stk_charge_20_100_xtrl_12_100 = fr.Filter("evtfilter_psd_fiducial_volume==false")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                        .Histo2D({"h_stk_charge_20_100_xtrl_12_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_100_250_xtrl_12_100 = fr.Filter("evtfilter_psd_fiducial_volume==false")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                        .Histo2D({"h_stk_charge_100_250_xtrl_12_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_250_500_xtrl_12_100 = fr.Filter("evtfilter_psd_fiducial_volume==false")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                        .Histo2D({"h_stk_charge_250_500_xtrl_12_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_500_1000_xtrl_12_100 = fr.Filter("evtfilter_psd_fiducial_volume==false")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                        .Histo2D({"h_stk_charge_500_1000_xtrl_12_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_1000_3000_xtrl_12_100 = fr.Filter("evtfilter_psd_fiducial_volume==false")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                        .Histo2D({"h_stk_charge_1000_3000_xtrl_12_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_3000_xtrl_12_100 = fr.Filter("evtfilter_psd_fiducial_volume==false")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=3000")
                                                        .Histo2D({"h_stk_charge_3000_xtrl_12_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    // PSD charge plots within PSD fiducial volume
    auto h_psd_charge_20_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                        .Histo2D({"h_psd_charge_20_100", "PSD charge; PSD X Charge; PSD Y Charge", 200, 0, 50, 200, 0, 50}, "PSD_chargeX", "PSD_chargeY");
    
    auto h_psd_charge_100_250 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                        .Histo2D({"h_psd_charge_100_250", "PSD charge; PSD X Charge; PSD Y Charge", 200, 0, 50, 200, 0, 50}, "PSD_chargeX", "PSD_chargeY");
    
    auto h_psd_charge_250_500 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                        .Histo2D({"h_psd_charge_250_500", "PSD charge; PSD X Charge; PSD Y Charge", 200, 0, 50, 200, 0, 50}, "PSD_chargeX", "PSD_chargeY");

    auto h_psd_charge_500_1000 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                        .Histo2D({"h_psd_charge_500_1000", "PSD charge; PSD X Charge; PSD Y Charge", 200, 0, 50, 200, 0, 50}, "PSD_chargeX", "PSD_chargeY");

    auto h_psd_charge_1000_3000 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                        .Histo2D({"h_psd_charge_1000_3000", "PSD charge; PSD X Charge; PSD Y Charge", 200, 0, 50, 200, 0, 50}, "PSD_chargeX", "PSD_chargeY");

    auto h_psd_charge_3000 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=3000")
                                                        .Histo2D({"h_psd_charge_3000", "PSD charge; PSD X Charge; PSD Y Charge", 200, 0, 50, 200, 0, 50}, "PSD_chargeX", "PSD_chargeY");

    // PSD charge plots within PSD fiducial volume (xtrl signal cut)
    auto h_psd_charge_20_100_xtrl_12 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                        .Histo2D({"h_psd_charge_20_100_xtrl_12", "PSD charge; PSD X Charge; PSD Y Charge", 200, 0, 50, 200, 0, 50}, "PSD_chargeX", "PSD_chargeY");
    
    auto h_psd_charge_100_250_xtrl_12 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                        .Histo2D({"h_psd_charge_100_250_xtrl_12", "PSD charge; PSD X Charge; PSD Y Charge", 200, 0, 50, 200, 0, 50}, "PSD_chargeX", "PSD_chargeY");
    
    auto h_psd_charge_250_500_xtrl_12 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                        .Histo2D({"h_psd_charge_250_500_xtrl_12", "PSD charge; PSD X Charge; PSD Y Charge", 200, 0, 50, 200, 0, 50}, "PSD_chargeX", "PSD_chargeY");

    auto h_psd_charge_500_1000_xtrl_12 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                        .Histo2D({"h_psd_charge_500_1000_xtrl_12", "PSD charge; PSD X Charge; PSD Y Charge", 200, 0, 50, 200, 0, 50}, "PSD_chargeX", "PSD_chargeY");

    auto h_psd_charge_1000_3000_xtrl_12 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                        .Histo2D({"h_psd_charge_1000_3000_xtrl_12", "PSD charge; PSD X Charge; PSD Y Charge", 200, 0, 50, 200, 0, 50}, "PSD_chargeX", "PSD_chargeY");

    auto h_psd_charge_3000_xtrl_12 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt<12")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=3000")
                                                        .Histo2D({"h_psd_charge_3000_xtrl_12", "PSD charge; PSD X Charge; PSD Y Charge", 200, 0, 50, 200, 0, 50}, "PSD_chargeX", "PSD_chargeY");

    // PSD charge plots within PSD fiducial volume (xtrl background cut)
    auto h_psd_charge_20_100_xtrl_12_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                        .Histo2D({"h_psd_charge_20_100_xtrl_12_100", "PSD charge; PSD X Charge; PSD Y Charge", 200, 0, 50, 200, 0, 50}, "PSD_chargeX", "PSD_chargeY");
    
    auto h_psd_charge_100_250_xtrl_12_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                        .Histo2D({"h_psd_charge_100_250_xtrl_12_100", "PSD charge; PSD X Charge; PSD Y Charge", 200, 0, 50, 200, 0, 50}, "PSD_chargeX", "PSD_chargeY");
    
    auto h_psd_charge_250_500_xtrl_12_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                        .Histo2D({"h_psd_charge_250_500_xtrl_12_100", "PSD charge; PSD X Charge; PSD Y Charge", 200, 0, 50, 200, 0, 50}, "PSD_chargeX", "PSD_chargeY");

    auto h_psd_charge_500_1000_xtrl_12_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                        .Histo2D({"h_psd_charge_500_1000_xtrl_12_100", "PSD charge; PSD X Charge; PSD Y Charge", 200, 0, 50, 200, 0, 50}, "PSD_chargeX", "PSD_chargeY");

    auto h_psd_charge_1000_3000_xtrl_12_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                        .Histo2D({"h_psd_charge_1000_3000_xtrl_12_100", "PSD charge; PSD X Charge; PSD Y Charge", 200, 0, 50, 200, 0, 50}, "PSD_chargeX", "PSD_chargeY");

    auto h_psd_charge_3000_xtrl_12_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                        .Filter("xtrl_evt>12 && xtrl<100")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=3000")
                                                        .Histo2D({"h_psd_charge_3000_xtrl_12_100", "PSD charge; PSD X Charge; PSD Y Charge", 200, 0, 50, 200, 0, 50}, "PSD_chargeX", "PSD_chargeY");     

    // STK cleaning cuts after all cuts
    auto createLogBinning = [](const double min, const double max, const std::size_t n_bins) -> std::vector<double> {
        std::vector<double> binning (n_bins + 1, 0);
        double log_interval = (log10(max) - log10(min)) / n_bins;
        for (unsigned int bIdx = 0; bIdx <= n_bins; ++bIdx)
            binning[bIdx] = pow(10, log10(min) + bIdx * log_interval);

        return binning;
    };

    auto nStkClu1Rm_log_binning = createLogBinning(1, 1e+3, 80);
    auto StkEcore1rm_log_binning = createLogBinning(1e+2, 1e+6, 500);

    auto h_stk_cleaning_20_100 = fr.Filter("evtfilter_all_cut==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                    .Histo2D({"h_stk_cleaning_20_100", "STK cleaning cut; nStkClu1Rm; StkEcore1rm", (int)nStkClu1Rm_log_binning.size()-1, &nStkClu1Rm_log_binning[0], (int)StkEcore1rm_log_binning.size()-1, &StkEcore1rm_log_binning[0]}, "nStkClu1Rm_stk", "stkEcore1Rm_stk");

    auto h_stk_cleaning_100_250 = fr.Filter("evtfilter_all_cut==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                    .Histo2D({"h_stk_cleaning_100_250", "STK cleaning cut; nStkClu1Rm; StkEcore1rm", (int)nStkClu1Rm_log_binning.size()-1, &nStkClu1Rm_log_binning[0], (int)StkEcore1rm_log_binning.size()-1, &StkEcore1rm_log_binning[0]}, "nStkClu1Rm_stk", "stkEcore1Rm_stk");

    auto h_stk_cleaning_250_500 = fr.Filter("evtfilter_all_cut==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                    .Histo2D({"h_stk_cleaning_250_500", "STK cleaning cut; nStkClu1Rm; StkEcore1rm", (int)nStkClu1Rm_log_binning.size()-1, &nStkClu1Rm_log_binning[0], (int)StkEcore1rm_log_binning.size()-1, &StkEcore1rm_log_binning[0]}, "nStkClu1Rm_stk", "stkEcore1Rm_stk");

    auto h_stk_cleaning_500_1000 = fr.Filter("evtfilter_all_cut==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                    .Histo2D({"h_stk_cleaning_500_1000", "STK cleaning cut; nStkClu1Rm; StkEcore1rm", (int)nStkClu1Rm_log_binning.size()-1, &nStkClu1Rm_log_binning[0], (int)StkEcore1rm_log_binning.size()-1, &StkEcore1rm_log_binning[0]}, "nStkClu1Rm_stk", "stkEcore1Rm_stk");

    auto h_stk_cleaning_1000_3000 = fr.Filter("evtfilter_all_cut==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                    .Histo2D({"h_stk_cleaning_1000_3000", "STK cleaning cut; nStkClu1Rm; StkEcore1rm", (int)nStkClu1Rm_log_binning.size()-1, &nStkClu1Rm_log_binning[0], (int)StkEcore1rm_log_binning.size()-1, &StkEcore1rm_log_binning[0]}, "nStkClu1Rm_stk", "stkEcore1Rm_stk");

    auto h_stk_cleaning_3000 = fr.Filter("evtfilter_all_cut==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=3000")
                                    .Histo2D({"h_stk_cleaning_3000", "STK cleaning cut; nStkClu1Rm; StkEcore1rm", (int)nStkClu1Rm_log_binning.size()-1, &nStkClu1Rm_log_binning[0], (int)StkEcore1rm_log_binning.size()-1, &StkEcore1rm_log_binning[0]}, "nStkClu1Rm_stk", "stkEcore1Rm_stk");

    // STK cleaning cuts after all cuts (xtrl signal cut)
    auto h_stk_cleaning_20_100_xtrl_12 = fr.Filter("evtfilter_all_cut==true")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter("xtrl_evt<12")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                    .Histo2D({"h_stk_cleaning_20_100_xtrl_12", "STK cleaning cut; nStkClu1Rm; StkEcore1rm", (int)nStkClu1Rm_log_binning.size()-1, &nStkClu1Rm_log_binning[0], (int)StkEcore1rm_log_binning.size()-1, &StkEcore1rm_log_binning[0]}, "nStkClu1Rm_stk", "stkEcore1Rm_stk");

    auto h_stk_cleaning_100_250_xtrl_12 = fr.Filter("evtfilter_all_cut==true")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter("xtrl_evt<12")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                    .Histo2D({"h_stk_cleaning_100_250_xtrl_12", "STK cleaning cut; nStkClu1Rm; StkEcore1rm", (int)nStkClu1Rm_log_binning.size()-1, &nStkClu1Rm_log_binning[0], (int)StkEcore1rm_log_binning.size()-1, &StkEcore1rm_log_binning[0]}, "nStkClu1Rm_stk", "stkEcore1Rm_stk");

    auto h_stk_cleaning_250_500_xtrl_12 = fr.Filter("evtfilter_all_cut==true")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter("xtrl_evt<12")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                    .Histo2D({"h_stk_cleaning_250_500_xtrl_12", "STK cleaning cut; nStkClu1Rm; StkEcore1rm", (int)nStkClu1Rm_log_binning.size()-1, &nStkClu1Rm_log_binning[0], (int)StkEcore1rm_log_binning.size()-1, &StkEcore1rm_log_binning[0]}, "nStkClu1Rm_stk", "stkEcore1Rm_stk");

    auto h_stk_cleaning_500_1000_xtrl_12 = fr.Filter("evtfilter_all_cut==true")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter("xtrl_evt<12")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                    .Histo2D({"h_stk_cleaning_500_1000_xtrl_12", "STK cleaning cut; nStkClu1Rm; StkEcore1rm", (int)nStkClu1Rm_log_binning.size()-1, &nStkClu1Rm_log_binning[0], (int)StkEcore1rm_log_binning.size()-1, &StkEcore1rm_log_binning[0]}, "nStkClu1Rm_stk", "stkEcore1Rm_stk");

    auto h_stk_cleaning_1000_3000_xtrl_12 = fr.Filter("evtfilter_all_cut==true")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter("xtrl_evt<12")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                    .Histo2D({"h_stk_cleaning_1000_3000_xtrl_12", "STK cleaning cut; nStkClu1Rm; StkEcore1rm", (int)nStkClu1Rm_log_binning.size()-1, &nStkClu1Rm_log_binning[0], (int)StkEcore1rm_log_binning.size()-1, &StkEcore1rm_log_binning[0]}, "nStkClu1Rm_stk", "stkEcore1Rm_stk");

    auto h_stk_cleaning_3000_xtrl_12 = fr.Filter("evtfilter_all_cut==true")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter("xtrl_evt<12")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=3000")
                                    .Histo2D({"h_stk_cleaning_3000_xtrl_12", "STK cleaning cut; nStkClu1Rm; StkEcore1rm", (int)nStkClu1Rm_log_binning.size()-1, &nStkClu1Rm_log_binning[0], (int)StkEcore1rm_log_binning.size()-1, &StkEcore1rm_log_binning[0]}, "nStkClu1Rm_stk", "stkEcore1Rm_stk");

    // STK cleaning cuts after all cuts (xtrl background cut)
    auto h_stk_cleaning_20_100_xtrl_12_100 = fr.Filter("evtfilter_all_cut==true")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter("xtrl_evt>12 && xtrl<100")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                    .Histo2D({"h_stk_cleaning_20_100_xtrl_12_100", "STK cleaning cut; nStkClu1Rm; StkEcore1rm", (int)nStkClu1Rm_log_binning.size()-1, &nStkClu1Rm_log_binning[0], (int)StkEcore1rm_log_binning.size()-1, &StkEcore1rm_log_binning[0]}, "nStkClu1Rm_stk", "stkEcore1Rm_stk");

    auto h_stk_cleaning_100_250_xtrl_12_100 = fr.Filter("evtfilter_all_cut==true")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter("xtrl_evt>12 && xtrl<100")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                    .Histo2D({"h_stk_cleaning_100_250_xtrl_12_100", "STK cleaning cut; nStkClu1Rm; StkEcore1rm", (int)nStkClu1Rm_log_binning.size()-1, &nStkClu1Rm_log_binning[0], (int)StkEcore1rm_log_binning.size()-1, &StkEcore1rm_log_binning[0]}, "nStkClu1Rm_stk", "stkEcore1Rm_stk");

    auto h_stk_cleaning_250_500_xtrl_12_100 = fr.Filter("evtfilter_all_cut==true")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter("xtrl_evt>12 && xtrl<100")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                    .Histo2D({"h_stk_cleaning_250_500_xtrl_12_100", "STK cleaning cut; nStkClu1Rm; StkEcore1rm", (int)nStkClu1Rm_log_binning.size()-1, &nStkClu1Rm_log_binning[0], (int)StkEcore1rm_log_binning.size()-1, &StkEcore1rm_log_binning[0]}, "nStkClu1Rm_stk", "stkEcore1Rm_stk");

    auto h_stk_cleaning_500_1000_xtrl_12_100 = fr.Filter("evtfilter_all_cut==true")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter("xtrl_evt>12 && xtrl<100")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                    .Histo2D({"h_stk_cleaning_500_1000_xtrl_12_100", "STK cleaning cut; nStkClu1Rm; StkEcore1rm", (int)nStkClu1Rm_log_binning.size()-1, &nStkClu1Rm_log_binning[0], (int)StkEcore1rm_log_binning.size()-1, &StkEcore1rm_log_binning[0]}, "nStkClu1Rm_stk", "stkEcore1Rm_stk");

    auto h_stk_cleaning_1000_3000_xtrl_12_100 = fr.Filter("evtfilter_all_cut==true")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter("xtrl_evt>12 && xtrl<100")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                    .Histo2D({"h_stk_cleaning_1000_3000_xtrl_12_100", "STK cleaning cut; nStkClu1Rm; StkEcore1rm", (int)nStkClu1Rm_log_binning.size()-1, &nStkClu1Rm_log_binning[0], (int)StkEcore1rm_log_binning.size()-1, &StkEcore1rm_log_binning[0]}, "nStkClu1Rm_stk", "stkEcore1Rm_stk");

    auto h_stk_cleaning_3000_xtrl_12_100 = fr.Filter("evtfilter_all_cut==true")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter("xtrl_evt>12 && xtrl<100")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=3000")
                                    .Histo2D({"h_stk_cleaning_3000_xtrl_12_100", "STK cleaning cut; nStkClu1Rm; StkEcore1rm", (int)nStkClu1Rm_log_binning.size()-1, &nStkClu1Rm_log_binning[0], (int)StkEcore1rm_log_binning.size()-1, &StkEcore1rm_log_binning[0]}, "nStkClu1Rm_stk", "stkEcore1Rm_stk");

    
    // Rvalue variable
    auto createLinearBinning = [](const double min, const double max, const std::size_t n_bins) -> std::vector<double>
    {
        double h = (max - min) / n_bins;
        std::vector<double> xs(n_bins + 1);
        std::vector<double>::iterator x;
        double val;
        for (x = xs.begin(), val = min; x != xs.end(); ++x, val += h)
            *x = val;

        return xs;
    };
    auto rvalue_binning = createLinearBinning(0, 10, 1e+2);
    auto lvalue_binning = createLinearBinning(0, 1, 1e+2);

    auto h_rvalue = fr.Filter("evtfilter_all_cut==true")
                        .Define("corr_energy_gev", "energy_corr * 0.001")
                        .Histo2D({"h_rvalue", "R value; Corrected Energy [GeV]; RValue", energy_nbins, &energy_binning[0], (int)rvalue_binning.size()-1, &rvalue_binning[0]}, "corr_energy_gev", "rvalue");
    
    auto h_rvalue_xtrl_12 = fr.Filter("evtfilter_all_cut==true")
                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                .Filter("xtrl_evt<12")
                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                .Histo2D({"h_rvalue_xtrl_12", "R value; Corrected Energy [GeV]; RValue", energy_nbins, &energy_binning[0], (int)rvalue_binning.size()-1, &rvalue_binning[0]}, "corr_energy_gev", "rvalue");

    auto h_rvalue_xtrl_12_100 = fr.Filter("evtfilter_all_cut==true")
                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                .Filter("xtrl_evt>12 && xtrl<100")
                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                .Histo2D({"h_rvalue_xtrl_12_100", "R value; Corrected Energy [GeV]; RValue", energy_nbins, &energy_binning[0], (int)rvalue_binning.size()-1, &rvalue_binning[0]}, "corr_energy_gev", "rvalue");
    
    // Lvalue variable
    auto h_lvalue = fr.Filter("evtfilter_all_cut==true")
                        .Define("corr_energy_gev", "energy_corr * 0.001")
                        .Histo2D({"h_lvalue", "L value; Corrected Energy [GeV]; LValue", energy_nbins, &energy_binning[0], (int)lvalue_binning.size()-1, &lvalue_binning[0]}, "corr_energy_gev", "lvalue");
    
    auto h_lvalue_xtrl_12 = fr.Filter("evtfilter_all_cut==true")
                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                .Filter("xtrl_evt<12")
                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                .Histo2D({"h_lvalue_xtrl_12", "L value; Corrected Energy [GeV]; LValue", energy_nbins, &energy_binning[0], (int)lvalue_binning.size()-1, &lvalue_binning[0]}, "corr_energy_gev", "lvalue");

    auto h_lvalue_xtrl_12_100 = fr.Filter("evtfilter_all_cut==true")
                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                .Filter("xtrl_evt>12 && xtrl<100")
                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                .Histo2D({"h_lvalue_xtrl_12_100", "L value; Corrected Energy [GeV]; LValue", energy_nbins, &energy_binning[0], (int)lvalue_binning.size()-1, &lvalue_binning[0]}, "corr_energy_gev", "lvalue");

    // sumRMS vs BGO polar angle cosine
    auto sumRms_cosine_bins = createLogBinning(10, 3e+3, 1e+2);
    auto cosine_bins = createLinearBinning(0, 1, 1e+2);

    auto h_sumrms_cosine_20_100 = fr.Filter("evtfilter_maxRms_cut==true")
                        .Define("corr_energy_gev", "energy_corr * 0.001")
                        .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                        .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                        .Histo2D({"h_sumrms_cosine_20_100", "sumRMS vs BGO polar angle cosine; #cos(#theta); sumRMS [mm]", (int)cosine_bins.size() -1, &cosine_bins[0], (int)sumRms_cosine_bins.size()-1, &sumRms_cosine_bins[0]}, "bgorec_cosine", "sumRms");

    auto h_sumrms_cosine_100_250 = fr.Filter("evtfilter_maxRms_cut==true")
                        .Define("corr_energy_gev", "energy_corr * 0.001")
                        .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                        .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                        .Histo2D({"h_sumrms_cosine_100_250", "sumRMS vs BGO polar angle cosine; #cos(#theta); sumRMS [mm]", (int)cosine_bins.size() -1, &cosine_bins[0], (int)sumRms_cosine_bins.size()-1, &sumRms_cosine_bins[0]}, "bgorec_cosine", "sumRms");

    // STK track selection distance variables (events rejected by the classifiers)
    auto h_trackselection_distance_vars_20_100_xtrl_tight = fr.Filter("track_efficiency_preselection_accepted==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter("xtrl_evt>8.5 && xtrl_evt!= -999")
                                    .Define("STK_bestTrack_linear_distance_STK_BGO", "sqrt(pow(STK_bestTrack_STK_BGO_topX_distance, 2) + pow(STK_bestTrack_STK_BGO_topY_distance, 2))")
                                    .Histo2D({"h_trackselection_distance_vars_20_100_xtrl_tight", "STK Track Selection distancest; linear distance [mm]; angular distance [deg]", 100, 0, 100, 100, 0, 50}, "STK_bestTrack_linear_distance_STK_BGO", "STK_bestTrack_angular_distance_STK_BGO");

    auto h_trackselection_distance_vars_100_250_xtrl_tight = fr.Filter("track_efficiency_preselection_accepted==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter("xtrl_evt>8.5 && xtrl_evt!= -999")
                                    .Define("STK_bestTrack_linear_distance_STK_BGO", "sqrt(pow(STK_bestTrack_STK_BGO_topX_distance, 2) + pow(STK_bestTrack_STK_BGO_topY_distance, 2))")
                                    .Histo2D({"h_trackselection_distance_vars_100_250_xtrl_tight", "STK Track Selection distancest; linear distance [mm]; angular distance [deg]", 100, 0, 100, 100, 0, 50}, "STK_bestTrack_linear_distance_STK_BGO", "STK_bestTrack_angular_distance_STK_BGO");

    auto h_trackselection_distance_vars_250_500_xtrl_tight = fr.Filter("track_efficiency_preselection_accepted==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter("xtrl_evt>8.5 && xtrl_evt!= -999")
                                    .Define("STK_bestTrack_linear_distance_STK_BGO", "sqrt(pow(STK_bestTrack_STK_BGO_topX_distance, 2) + pow(STK_bestTrack_STK_BGO_topY_distance, 2))")
                                    .Histo2D({"h_trackselection_distance_vars_250_500_xtrl_tight", "STK Track Selection distancest; linear distance [mm]; angular distance [deg]", 100, 0, 100, 100, 0, 50}, "STK_bestTrack_linear_distance_STK_BGO", "STK_bestTrack_angular_distance_STK_BGO");

    auto h_trackselection_distance_vars_500_1000_xtrl_tight = fr.Filter("track_efficiency_preselection_accepted==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter("xtrl_evt>8.5 && xtrl_evt!= -999")
                                    .Define("STK_bestTrack_linear_distance_STK_BGO", "sqrt(pow(STK_bestTrack_STK_BGO_topX_distance, 2) + pow(STK_bestTrack_STK_BGO_topY_distance, 2))")
                                    .Histo2D({"h_trackselection_distance_vars_500_1000_xtrl_tight", "STK Track Selection distancest; linear distance [mm]; angular distance [deg]", 100, 0, 100, 100, 0, 50}, "STK_bestTrack_linear_distance_STK_BGO", "STK_bestTrack_angular_distance_STK_BGO");

    auto h_trackselection_distance_vars_1000_3000_xtrl_tight = fr.Filter("track_efficiency_preselection_accepted==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter("xtrl_evt>8.5 && xtrl_evt!= -999")
                                    .Define("STK_bestTrack_linear_distance_STK_BGO", "sqrt(pow(STK_bestTrack_STK_BGO_topX_distance, 2) + pow(STK_bestTrack_STK_BGO_topY_distance, 2))")
                                    .Histo2D({"h_trackselection_distance_vars_1000_3000_xtrl_tight", "STK Track Selection distancest; linear distance [mm]; angular distance [deg]", 100, 0, 100, 100, 0, 50}, "STK_bestTrack_linear_distance_STK_BGO", "STK_bestTrack_angular_distance_STK_BGO");

    auto h_trackselection_distance_vars_3000_xtrl_tight = fr.Filter("track_efficiency_preselection_accepted==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("corr_energy_gev>=3000")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter("xtrl_evt>8.5 && xtrl_evt!= -999")
                                    .Define("STK_bestTrack_linear_distance_STK_BGO", "sqrt(pow(STK_bestTrack_STK_BGO_topX_distance, 2) + pow(STK_bestTrack_STK_BGO_topY_distance, 2))")
                                    .Histo2D({"h_trackselection_distance_vars_3000_xtrl_tight", "STK Track Selection distancest; linear distance [mm]; angular distance [deg]", 100, 0, 100, 100, 0, 50}, "STK_bestTrack_linear_distance_STK_BGO", "STK_bestTrack_angular_distance_STK_BGO");

    auto h_trackselection_distance_vars_20_100_xtrl_loose = fr.Filter("track_efficiency_preselection_accepted==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Define("energy_gev", "energy * 0.001")
                                    .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter(xtrl_loose_reject_cut, {"energy_gev", "xtrl_evt"})
                                    .Define("STK_bestTrack_linear_distance_STK_BGO", "sqrt(pow(STK_bestTrack_STK_BGO_topX_distance, 2) + pow(STK_bestTrack_STK_BGO_topY_distance, 2))")
                                    .Histo2D({"h_trackselection_distance_vars_20_100_xtrl_loose", "STK Track Selection distancest; linear distance [mm]; angular distance [deg]", 100, 0, 100, 100, 0, 50}, "STK_bestTrack_linear_distance_STK_BGO", "STK_bestTrack_angular_distance_STK_BGO");

    auto h_trackselection_distance_vars_100_250_xtrl_loose = fr.Filter("track_efficiency_preselection_accepted==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Define("energy_gev", "energy * 0.001")
                                    .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter(xtrl_loose_reject_cut, {"energy_gev", "xtrl_evt"})
                                    .Define("STK_bestTrack_linear_distance_STK_BGO", "sqrt(pow(STK_bestTrack_STK_BGO_topX_distance, 2) + pow(STK_bestTrack_STK_BGO_topY_distance, 2))")
                                    .Histo2D({"h_trackselection_distance_vars_100_250_xtrl_loose", "STK Track Selection distancest; linear distance [mm]; angular distance [deg]", 100, 0, 100, 100, 0, 50}, "STK_bestTrack_linear_distance_STK_BGO", "STK_bestTrack_angular_distance_STK_BGO");

    auto h_trackselection_distance_vars_250_500_xtrl_loose = fr.Filter("track_efficiency_preselection_accepted==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Define("energy_gev", "energy * 0.001")
                                    .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter(xtrl_loose_reject_cut, {"energy_gev", "xtrl_evt"})
                                    .Define("STK_bestTrack_linear_distance_STK_BGO", "sqrt(pow(STK_bestTrack_STK_BGO_topX_distance, 2) + pow(STK_bestTrack_STK_BGO_topY_distance, 2))")
                                    .Histo2D({"h_trackselection_distance_vars_250_500_xtrl_loose", "STK Track Selection distancest; linear distance [mm]; angular distance [deg]", 100, 0, 100, 100, 0, 50}, "STK_bestTrack_linear_distance_STK_BGO", "STK_bestTrack_angular_distance_STK_BGO");

    auto h_trackselection_distance_vars_500_1000_xtrl_loose = fr.Filter("track_efficiency_preselection_accepted==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Define("energy_gev", "energy * 0.001")
                                    .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter(xtrl_loose_reject_cut, {"energy_gev", "xtrl_evt"})
                                    .Define("STK_bestTrack_linear_distance_STK_BGO", "sqrt(pow(STK_bestTrack_STK_BGO_topX_distance, 2) + pow(STK_bestTrack_STK_BGO_topY_distance, 2))")
                                    .Histo2D({"h_trackselection_distance_vars_500_1000_xtrl_loose", "STK Track Selection distancest; linear distance [mm]; angular distance [deg]", 100, 0, 100, 100, 0, 50}, "STK_bestTrack_linear_distance_STK_BGO", "STK_bestTrack_angular_distance_STK_BGO");

    auto h_trackselection_distance_vars_1000_3000_xtrl_loose = fr.Filter("track_efficiency_preselection_accepted==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Define("energy_gev", "energy * 0.001")
                                    .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter(xtrl_loose_reject_cut, {"energy_gev", "xtrl_evt"})
                                    .Define("STK_bestTrack_linear_distance_STK_BGO", "sqrt(pow(STK_bestTrack_STK_BGO_topX_distance, 2) + pow(STK_bestTrack_STK_BGO_topY_distance, 2))")
                                    .Histo2D({"h_trackselection_distance_vars_1000_3000_xtrl_loose", "STK Track Selection distancest; linear distance [mm]; angular distance [deg]", 100, 0, 100, 100, 0, 50}, "STK_bestTrack_linear_distance_STK_BGO", "STK_bestTrack_angular_distance_STK_BGO");

    auto h_trackselection_distance_vars_3000_xtrl_loose = fr.Filter("track_efficiency_preselection_accepted==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Define("energy_gev", "energy * 0.001")
                                    .Filter("corr_energy_gev>=3000")
                                    .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                    .Filter(xtrl_loose_reject_cut, {"energy_gev", "xtrl_evt"})
                                    .Define("STK_bestTrack_linear_distance_STK_BGO", "sqrt(pow(STK_bestTrack_STK_BGO_topX_distance, 2) + pow(STK_bestTrack_STK_BGO_topY_distance, 2))")
                                    .Histo2D({"h_trackselection_distance_vars_3000_xtrl_loose", "STK Track Selection distancest; linear distance [mm]; angular distance [deg]", 100, 0, 100, 100, 0, 50}, "STK_bestTrack_linear_distance_STK_BGO", "STK_bestTrack_angular_distance_STK_BGO");

    auto h_trackselection_distance_vars_20_100_bdt = fr.Filter("track_efficiency_preselection_accepted==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Define("energy_gev", "energy * 0.001")
                                    .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                    .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                    .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                    .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                    .Filter([] (const double tmva_value, const double tmva_cut) {return tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                    .Define("STK_bestTrack_linear_distance_STK_BGO", "sqrt(pow(STK_bestTrack_STK_BGO_topX_distance, 2) + pow(STK_bestTrack_STK_BGO_topY_distance, 2))")
                                    .Histo2D({"h_trackselection_distance_vars_20_100_bdt", "STK Track Selection distancest; linear distance [mm]; angular distance [deg]", 100, 0, 100, 100, 0, 50}, "STK_bestTrack_linear_distance_STK_BGO", "STK_bestTrack_angular_distance_STK_BGO");

    auto h_trackselection_distance_vars_100_250_bdt = fr.Filter("track_efficiency_preselection_accepted==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Define("energy_gev", "energy * 0.001")
                                    .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                    .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                    .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                    .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                    .Filter([] (const double tmva_value, const double tmva_cut) {return tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                    .Define("STK_bestTrack_linear_distance_STK_BGO", "sqrt(pow(STK_bestTrack_STK_BGO_topX_distance, 2) + pow(STK_bestTrack_STK_BGO_topY_distance, 2))")
                                    .Histo2D({"h_trackselection_distance_vars_100_250_bdt", "STK Track Selection distancest; linear distance [mm]; angular distance [deg]", 100, 0, 100, 100, 0, 50}, "STK_bestTrack_linear_distance_STK_BGO", "STK_bestTrack_angular_distance_STK_BGO");

    auto h_trackselection_distance_vars_250_500_bdt = fr.Filter("track_efficiency_preselection_accepted==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Define("energy_gev", "energy * 0.001")
                                    .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                    .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                    .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                    .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                    .Filter([] (const double tmva_value, const double tmva_cut) {return tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                    .Define("STK_bestTrack_linear_distance_STK_BGO", "sqrt(pow(STK_bestTrack_STK_BGO_topX_distance, 2) + pow(STK_bestTrack_STK_BGO_topY_distance, 2))")
                                    .Histo2D({"h_trackselection_distance_vars_250_500_bdt", "STK Track Selection distancest; linear distance [mm]; angular distance [deg]", 100, 0, 100, 100, 0, 50}, "STK_bestTrack_linear_distance_STK_BGO", "STK_bestTrack_angular_distance_STK_BGO");

    auto h_trackselection_distance_vars_500_1000_bdt = fr.Filter("track_efficiency_preselection_accepted==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Define("energy_gev", "energy * 0.001")
                                    .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                    .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                    .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                    .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                    .Filter([] (const double tmva_value, const double tmva_cut) {return tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                    .Define("STK_bestTrack_linear_distance_STK_BGO", "sqrt(pow(STK_bestTrack_STK_BGO_topX_distance, 2) + pow(STK_bestTrack_STK_BGO_topY_distance, 2))")
                                    .Histo2D({"h_trackselection_distance_vars_500_1000_bdt", "STK Track Selection distancest; linear distance [mm]; angular distance [deg]", 100, 0, 100, 100, 0, 50}, "STK_bestTrack_linear_distance_STK_BGO", "STK_bestTrack_angular_distance_STK_BGO");

    auto h_trackselection_distance_vars_1000_3000_bdt = fr.Filter("track_efficiency_preselection_accepted==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Define("energy_gev", "energy * 0.001")
                                    .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                    .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                    .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                    .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                    .Filter([] (const double tmva_value, const double tmva_cut) {return tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                    .Define("STK_bestTrack_linear_distance_STK_BGO", "sqrt(pow(STK_bestTrack_STK_BGO_topX_distance, 2) + pow(STK_bestTrack_STK_BGO_topY_distance, 2))")
                                    .Histo2D({"h_trackselection_distance_vars_1000_3000_bdt", "STK Track Selection distancest; linear distance [mm]; angular distance [deg]", 100, 0, 100, 100, 0, 50}, "STK_bestTrack_linear_distance_STK_BGO", "STK_bestTrack_angular_distance_STK_BGO");

    auto h_trackselection_distance_vars_3000_bdt = fr.Filter("track_efficiency_preselection_accepted==true")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Define("energy_gev", "energy * 0.001")
                                    .Filter("corr_energy_gev>=3000")
                                    .Define("energy_bin", get_energy_bin, {"corr_energy_gev"})
                                    .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                    .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                    .Filter([] (const double tmva_value, const double tmva_cut) {return tmva_value < tmva_cut; }, {"bdt_evt", "bdt_cut"})
                                    .Define("STK_bestTrack_linear_distance_STK_BGO", "sqrt(pow(STK_bestTrack_STK_BGO_topX_distance, 2) + pow(STK_bestTrack_STK_BGO_topY_distance, 2))")
                                    .Histo2D({"h_trackselection_distance_vars_3000_bdt", "STK Track Selection distancest; linear distance [mm]; angular distance [deg]", 100, 0, 100, 100, 0, 50}, "STK_bestTrack_linear_distance_STK_BGO", "STK_bestTrack_angular_distance_STK_BGO");

    TFile *output_file = TFile::Open(input_args.output_path.c_str(), "RECREATE");
    if (output_file->IsZombie()) {
        std::cerr << "Error writing output ROOT file [" << input_args.output_path << "]\n\n";
        exit(100);
    }

    h_trigger_efficiency_accepted_het_over_let_tight_xtrl               ->Write();
    h_trigger_efficiency_accepted_het_over_unb_tight_xtrl               ->Write();
    h_trigger_efficiency_accepted_het_let_tight_xtrl                    ->Write();
    h_trigger_efficiency_accepted_het_unb_tight_xtrl                    ->Write();
    h_trigger_efficiency_accepted_het_over_let_bdt                      ->Write();
    h_trigger_efficiency_accepted_het_over_unb_bdt                      ->Write();
    h_trigger_efficiency_accepted_het_let_bdt                           ->Write();
    h_trigger_efficiency_accepted_het_unb_bdt                           ->Write();
    h_maxrms_efficiency_accepted_tight_xtrl                             ->Write();
    h_maxrms_efficiency_total_tight_xtrl                                ->Write();
    h_maxrms_efficiency_accepted_bdt                                    ->Write();
    h_maxrms_efficiency_total_bdt                                       ->Write();
    h_nbarlayer13_efficiency_accepted_tight_xtrl                        ->Write();
    h_nbarlayer13_efficiency_total_tight_xtrl                           ->Write();
    h_nbarlayer13_efficiency_accepted_bdt                               ->Write();
    h_nbarlayer13_efficiency_total_bdt                                  ->Write();
    h_maxrms_and_nbarlayer13_efficiency_accepted_tight_xtrl             ->Write();
    h_maxrms_and_nbarlayer13_efficiency_total_tight_xtrl                ->Write();
    h_maxrms_and_nbarlayer13_efficiency_accepted_bdt                    ->Write();
    h_maxrms_and_nbarlayer13_efficiency_total_bdt                       ->Write();
    h_sumrms_low_energy_efficiency_accepted_tight_xtrl                  ->Write();
    h_sumrms_low_energy_efficiency_total_tight_xtrl                     ->Write();
    h_sumrms_low_energy_efficiency_accepted_loose_xtrl                  ->Write();
    h_sumrms_low_energy_efficiency_total_loose_xtrl                     ->Write();
    h_sumrms_low_energy_efficiency_accepted_bdt                         ->Write();
    h_sumrms_low_energy_efficiency_total_bdt                            ->Write();
    h_track_efficiency_accepted_tight_xtrl                              ->Write();
    h_track_efficiency_total_tight_xtrl                                 ->Write();
    h_track_efficiency_accepted_bdt                                     ->Write();
    h_track_efficiency_total_bdt                                        ->Write();
    h_stk_1rm_efficiency_accepted_tight_xtrl                            ->Write();
    h_stk_1rm_efficiency_total_tight_xtrl                               ->Write();
    h_stk_1rm_efficiency_accepted_loose_xtrl                            ->Write();
    h_stk_1rm_efficiency_total_loose_xtrl                               ->Write();
    h_stk_1rm_efficiency_accepted_bdt                                   ->Write();
    h_stk_1rm_efficiency_total_bdt                                      ->Write();

    h_clusters_on_first_STK_layer_within_psd_fvolume_accepted_tight_xtrl           ->Write();
    h_clusters_on_first_STK_layer_within_psd_fvolume_total_tight_xtrl              ->Write();
    h_clusters_on_first_STK_layer_outside_psd_fvolume_accepted_tight_xtrl          ->Write();
    h_clusters_on_first_STK_layer_outside_psd_fvolume_total_tight_xtrl             ->Write();

    h_clusters_on_first_STK_layer_within_psd_fvolume_accepted_loose_xtrl           ->Write();
    h_clusters_on_first_STK_layer_within_psd_fvolume_total_loose_xtrl              ->Write();
    h_clusters_on_first_STK_layer_outside_psd_fvolume_accepted_loose_xtrl          ->Write();
    h_clusters_on_first_STK_layer_outside_psd_fvolume_total_loose_xtrl             ->Write();

    h_clusters_on_first_STK_layer_within_psd_fvolume_accepted_bdt                   ->Write();
    h_clusters_on_first_STK_layer_within_psd_fvolume_total_bdt                      ->Write();
    h_clusters_on_first_STK_layer_outside_psd_fvolume_accepted_bdt                  ->Write();
    h_clusters_on_first_STK_layer_outside_psd_fvolume_total_bdt                     ->Write();

    h_track_efficiency_stk_fvolume_accepted_tight_xtrl                  ->Write();
    h_track_efficiency_stk_fvolume_total_tight_xtrl                     ->Write();
    h_track_efficiency_stk_fvolume_accepted_bdt                         ->Write();
    h_track_efficiency_stk_fvolume_total_bdt                            ->Write();
    h_psdstkmatch_efficiency_accepted_tight_xtrl                        ->Write();
    h_psdstkmatch_efficiency_total_tight_xtrl                           ->Write();
    h_psdstkmatch_efficiency_accepted_bdt                               ->Write();
    h_psdstkmatch_efficiency_total_bdt                                  ->Write();
    h_psdcharge_efficiency_accepted_tight_xtrl                          ->Write();
    h_psdcharge_efficiency_total_tight_xtrl                             ->Write();
    h_psdcharge_efficiency_accepted_bdt                                 ->Write();
    h_psdcharge_efficiency_total_bdt                                    ->Write();
    h_stkcharge_efficiency_accepted_tight_xtrl                          ->Write();
    h_stkcharge_efficiency_total_tight_xtrl                             ->Write();
    h_stkcharge_efficiency_accepted_bdt                                 ->Write();
    h_stkcharge_efficiency_total_bdt                                    ->Write();
    h_trigger_efficiency_accepted_het_over_let_loose_xtrl               ->Write();
    h_trigger_efficiency_accepted_het_over_unb_loose_xtrl               ->Write();
    h_trigger_efficiency_accepted_het_let_loose_xtrl                    ->Write();
    h_trigger_efficiency_accepted_het_unb_loose_xtrl                    ->Write();
    h_maxrms_efficiency_accepted_loose_xtrl                             ->Write();
    h_maxrms_efficiency_total_loose_xtrl                                ->Write();
    h_nbarlayer13_efficiency_accepted_loose_xtrl                        ->Write();
    h_nbarlayer13_efficiency_total_loose_xtrl                           ->Write();
    h_maxrms_and_nbarlayer13_efficiency_accepted_loose_xtrl             ->Write();
    h_maxrms_and_nbarlayer13_efficiency_total_loose_xtrl                ->Write();
    h_track_efficiency_accepted_loose_xtrl                              ->Write();
    h_track_efficiency_total_loose_xtrl                                 ->Write();
    h_track_efficiency_stk_fvolume_accepted_loose_xtrl                  ->Write();
    h_track_efficiency_stk_fvolume_total_loose_xtrl                     ->Write();
    h_psdstkmatch_efficiency_accepted_loose_xtrl                        ->Write();
    h_psdstkmatch_efficiency_total_loose_xtrl                           ->Write();
    h_psdcharge_efficiency_accepted_loose_xtrl                          ->Write();
    h_psdcharge_efficiency_total_loose_xtrl                             ->Write();
    h_stkcharge_efficiency_accepted_loose_xtrl                          ->Write();
    h_stkcharge_efficiency_total_loose_xtrl                             ->Write();           

    output_file->mkdir("xtrl");
    output_file->cd("xtrl");

    h_xtrl_stk_cosine                                                   ->Write();
    h_xtrl_bgo_cosine                                                   ->Write();
    h_xtrl_stk_cosine_zoom                                              ->Write();
    h_xtrl_bgo_cosine_zoom                                              ->Write();
    h_psd_charge_after_stk_cut                                          ->Write();
    h_stk_charge_after_stk_cut                                          ->Write();
    h_stk_charge                                                        ->Write();
    h_stk_charge_after_psd_charge_cut                                   ->Write();

    output_file->mkdir("psd_matching_distance");
    output_file->cd("psd_matching_distance");

    h_psd_stk_match_distance_x_20_100                                   ->Write();
    h_psd_stk_match_distance_x_100_250                                  ->Write();
    h_psd_stk_match_distance_x_250_500                                  ->Write();
    h_psd_stk_match_distance_x_500_1000                                 ->Write();
    h_psd_stk_match_distance_x_1000_3000                                ->Write();
    h_psd_stk_match_distance_x_3000                                     ->Write();
    h_psd_stk_match_distance_y_20_100                                   ->Write();
    h_psd_stk_match_distance_y_100_250                                  ->Write();
    h_psd_stk_match_distance_y_250_500                                  ->Write();
    h_psd_stk_match_distance_y_500_1000                                 ->Write();
    h_psd_stk_match_distance_y_1000_3000                                ->Write();
    h_psd_stk_match_distance_y_3000                                     ->Write();

    h_psd_stk_match_distance_x_within_psd_fvolume_20_100                ->Write();
    h_psd_stk_match_distance_x_within_psd_fvolume_100_250               ->Write();
    h_psd_stk_match_distance_x_within_psd_fvolume_250_500               ->Write();
    h_psd_stk_match_distance_x_within_psd_fvolume_500_1000              ->Write();
    h_psd_stk_match_distance_x_within_psd_fvolume_1000_3000             ->Write();
    h_psd_stk_match_distance_x_within_psd_fvolume_3000                  ->Write();
    h_psd_stk_match_distance_y_within_psd_fvolume_20_100                ->Write();
    h_psd_stk_match_distance_y_within_psd_fvolume_100_250               ->Write();
    h_psd_stk_match_distance_y_within_psd_fvolume_250_500               ->Write();
    h_psd_stk_match_distance_y_within_psd_fvolume_500_1000              ->Write();
    h_psd_stk_match_distance_y_within_psd_fvolume_1000_3000             ->Write();
    h_psd_stk_match_distance_y_within_psd_fvolume_3000                  ->Write();

    h_psd_stk_match_distance_x_outside_psd_fvolume_20_100               ->Write();
    h_psd_stk_match_distance_x_outside_psd_fvolume_100_250              ->Write();
    h_psd_stk_match_distance_x_outside_psd_fvolume_250_500              ->Write();
    h_psd_stk_match_distance_x_outside_psd_fvolume_500_1000             ->Write();
    h_psd_stk_match_distance_x_outside_psd_fvolume_1000_3000            ->Write();
    h_psd_stk_match_distance_x_outside_psd_fvolume_3000                 ->Write();
    h_psd_stk_match_distance_y_outside_psd_fvolume_20_100               ->Write();
    h_psd_stk_match_distance_y_outside_psd_fvolume_100_250              ->Write();
    h_psd_stk_match_distance_y_outside_psd_fvolume_250_500              ->Write();
    h_psd_stk_match_distance_y_outside_psd_fvolume_500_1000             ->Write();
    h_psd_stk_match_distance_y_outside_psd_fvolume_1000_3000            ->Write();
    h_psd_stk_match_distance_y_outside_psd_fvolume_3000                 ->Write();

    output_file->mkdir("stk_cleaning_cut");
    output_file->cd("stk_cleaning_cut");

    h_stk_cleaning_20_100                                               ->Write();
    h_stk_cleaning_100_250                                              ->Write();
    h_stk_cleaning_250_500                                              ->Write();
    h_stk_cleaning_500_1000                                             ->Write();
    h_stk_cleaning_1000_3000                                            ->Write();
    h_stk_cleaning_3000                                                 ->Write();

    h_stk_cleaning_20_100_xtrl_12                                       ->Write();
    h_stk_cleaning_100_250_xtrl_12                                      ->Write();
    h_stk_cleaning_250_500_xtrl_12                                      ->Write();
    h_stk_cleaning_500_1000_xtrl_12                                     ->Write();
    h_stk_cleaning_1000_3000_xtrl_12                                    ->Write();
    h_stk_cleaning_3000_xtrl_12                                         ->Write();

    h_stk_cleaning_20_100_xtrl_12_100                                   ->Write();
    h_stk_cleaning_100_250_xtrl_12_100                                  ->Write();
    h_stk_cleaning_250_500_xtrl_12_100                                  ->Write();
    h_stk_cleaning_500_1000_xtrl_12_100                                 ->Write();
    h_stk_cleaning_1000_3000_xtrl_12_100                                ->Write();
    h_stk_cleaning_3000_xtrl_12_100                                     ->Write();

    output_file->mkdir("rvalue");
    output_file->cd("rvalue");

    h_rvalue                                                            ->Write();
    h_rvalue_xtrl_12                                                    ->Write();
    h_rvalue_xtrl_12_100                                                ->Write();

    output_file->mkdir("lvalue");
    output_file->cd("lvalue");

    h_lvalue                                                            ->Write();
    h_lvalue_xtrl_12                                                    ->Write();
    h_lvalue_xtrl_12_100                                                ->Write();

    output_file->mkdir("stk_charge");
    output_file->cd("stk_charge");

    h_stk_charge_psd_fvolume_no_psd_cut_20_100                          ->Write();
    h_stk_charge_psd_fvolume_no_psd_cut_100_250                         ->Write();
    h_stk_charge_psd_fvolume_no_psd_cut_250_500                         ->Write();
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000                        ->Write();
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000                       ->Write();
    h_stk_charge_psd_fvolume_no_psd_cut_3000                            ->Write();
     
    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12                  ->Write();
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12                 ->Write();
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12                 ->Write();
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12                ->Write();
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12               ->Write();
    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12                    ->Write();

    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_100              ->Write();
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_100             ->Write();
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_100             ->Write();
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_100            ->Write();
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_100           ->Write();
    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_100                ->Write();

    h_stk_charge_psd_fvolume_psd_cut_20_100                             ->Write();
    h_stk_charge_psd_fvolume_psd_cut_100_250                            ->Write();
    h_stk_charge_psd_fvolume_psd_cut_250_500                            ->Write();
    h_stk_charge_psd_fvolume_psd_cut_500_1000                           ->Write();
    h_stk_charge_psd_fvolume_psd_cut_1000_3000                          ->Write();
    h_stk_charge_psd_fvolume_psd_cut_3000                               ->Write();

    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12                     ->Write();
    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12                    ->Write();
    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12                    ->Write();
    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12                   ->Write();
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12                  ->Write();
    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12                       ->Write();

    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_100                 ->Write();
    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_100                ->Write();
    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_100                ->Write();
    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_100               ->Write();
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_100              ->Write();
    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_100                   ->Write();

    h_stk_charge_20_100                                                 ->Write();
    h_stk_charge_100_250                                                ->Write();
    h_stk_charge_250_500                                                ->Write();
    h_stk_charge_500_1000                                               ->Write();
    h_stk_charge_1000_3000                                              ->Write();
    h_stk_charge_3000                                                   ->Write();

    h_stk_charge_20_100_xtrl_12                                         ->Write();
    h_stk_charge_100_250_xtrl_12                                        ->Write();
    h_stk_charge_250_500_xtrl_12                                        ->Write();
    h_stk_charge_500_1000_xtrl_12                                       ->Write();
    h_stk_charge_1000_3000_xtrl_12                                      ->Write();
    h_stk_charge_3000_xtrl_12                                           ->Write();

    h_stk_charge_20_100_xtrl_12_100                                     ->Write();
    h_stk_charge_100_250_xtrl_12_100                                    ->Write();
    h_stk_charge_250_500_xtrl_12_100                                    ->Write();
    h_stk_charge_500_1000_xtrl_12_100                                   ->Write();
    h_stk_charge_1000_3000_xtrl_12_100                                  ->Write();
    h_stk_charge_3000_xtrl_12_100                                       ->Write();

    output_file->mkdir("psd_charge");
    output_file->cd("psd_charge");

    h_psd_charge_20_100                                                 ->Write();
    h_psd_charge_100_250                                                ->Write();
    h_psd_charge_250_500                                                ->Write();
    h_psd_charge_500_1000                                               ->Write();
    h_psd_charge_1000_3000                                              ->Write();
    h_psd_charge_3000                                                   ->Write();

    h_psd_charge_20_100_xtrl_12                                         ->Write();
    h_psd_charge_100_250_xtrl_12                                        ->Write();
    h_psd_charge_250_500_xtrl_12                                        ->Write();
    h_psd_charge_500_1000_xtrl_12                                       ->Write();
    h_psd_charge_1000_3000_xtrl_12                                      ->Write();
    h_psd_charge_3000_xtrl_12                                           ->Write();

    h_psd_charge_20_100_xtrl_12_100                                     ->Write();
    h_psd_charge_100_250_xtrl_12_100                                    ->Write();
    h_psd_charge_250_500_xtrl_12_100                                    ->Write();
    h_psd_charge_500_1000_xtrl_12_100                                   ->Write();
    h_psd_charge_1000_3000_xtrl_12_100                                  ->Write();
    h_psd_charge_3000_xtrl_12_100                                       ->Write();

    output_file->mkdir("sumrms");
    output_file->cd("sumrms");

    h_sumrms_cosine_20_100                                              ->Write();
    h_sumrms_cosine_100_250                                             ->Write();

    output_file->mkdir("trackselection_distance");
    output_file->cd("trackselection_distance");

    h_trackselection_distance_vars_20_100_xtrl_tight                    ->Write();
    h_trackselection_distance_vars_100_250_xtrl_tight                   ->Write();
    h_trackselection_distance_vars_250_500_xtrl_tight                   ->Write();
    h_trackselection_distance_vars_500_1000_xtrl_tight                  ->Write();
    h_trackselection_distance_vars_1000_3000_xtrl_tight                 ->Write();
    h_trackselection_distance_vars_3000_xtrl_tight                      ->Write();
    
    h_trackselection_distance_vars_20_100_xtrl_loose                    ->Write();
    h_trackselection_distance_vars_100_250_xtrl_loose                   ->Write();
    h_trackselection_distance_vars_250_500_xtrl_loose                   ->Write();
    h_trackselection_distance_vars_500_1000_xtrl_loose                  ->Write();
    h_trackselection_distance_vars_1000_3000_xtrl_loose                 ->Write();
    h_trackselection_distance_vars_3000_xtrl_loose                      ->Write();

    h_trackselection_distance_vars_20_100_bdt                           ->Write();
    h_trackselection_distance_vars_100_250_bdt                          ->Write();
    h_trackselection_distance_vars_250_500_bdt                          ->Write();
    h_trackselection_distance_vars_500_1000_bdt                         ->Write();
    h_trackselection_distance_vars_1000_3000_bdt                        ->Write();
    h_trackselection_distance_vars_3000_bdt                             ->Write();

    output_file->Close();
}   