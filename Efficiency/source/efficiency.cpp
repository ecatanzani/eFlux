#include "bdt.h"
#include "efficiency.h"
#include "list_parser.h"
#include "energy_config.h"

#include <memory>

#include "TFile.h"
#include <ROOT/RDataFrame.hxx>

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
    auto compute_bdt = [bdt_evt, energy_binning](
        const std::vector<double> rms_layer, 
        const double sumrms, 
        const std::vector<double> energy_frac_layer,
        const double energy_fraction_last_layer,
        const double energy_corr_gev,
        const TVector3 bgo_trajectory) -> double 
        {
            /*
            return bdt_evt->ComputeMVA(
                rms_layer,
                sumrms,
                energy_frac_layer,
                energy_fraction_last_layer,
                energy_corr_gev,
                energy_binning,
                bgo_trajectory);
            */

           return 10;
        };

    // Compute XTRL cut
    auto xtrl_loose_cut = [](const double bgo_total_energy_gev, const double xtrl) -> bool {
        auto xtrl_cut_value {3*log10(bgo_total_energy_gev/200) + 12};
        if (xtrl<xtrl_cut_value && xtrl!=-999)
            return true;
        else
            return false;
    };

    // Compute BDT cut
    const double bdt_cut_value=0.25;
    auto bdt_cut = [bdt_cut_value](const double bdt) -> bool {
        return bdt > bdt_cut_value;
    };

    if (input_args.verbose)
        std::cout << "\n\nAnalysis running\n\n";
    
    auto get_let_prescale = [](const double geographic_latitude) -> double {
        return geographic_latitude>=-20 && geographic_latitude<=20 ? 1/8. : 1/64.;
    };
    
    auto get_unb_prescale = [](const double geographic_latitude) -> double {
        return geographic_latitude>=-20 && geographic_latitude<=20 ? 1/512. : 1/2048.;
    };

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
    
    auto h_trigger_efficiency_accepted_het_over_let_bdt = mc_file ? 
                                            fr.Filter("trigger_efficiency_preselection==1 && trigger_efficiency_preselection_is_het==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Filter(bdt_cut, {"bdt_evt"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_over_let_bdt", "HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev")

                                            :

                                            fr.Filter("trigger_efficiency_preselection==1 && trigger_efficiency_preselection_is_het==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Filter(bdt_cut, {"bdt_evt"})
                                            .Define("trigger_w", get_let_prescale, {"geo_lat"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_over_let_bdt", "HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev", "trigger_w");

     auto h_trigger_efficiency_accepted_het_over_unb_bdt = mc_file ? 
                                            fr.Filter("trigger_efficiency_preselection==1 && trigger_efficiency_preselection_is_het==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Filter(bdt_cut, {"bdt_evt"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_over_unb_bdt", "HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev")

                                            :

                                            fr.Filter("trigger_efficiency_preselection==1 && trigger_efficiency_preselection_is_het==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Filter(bdt_cut, {"bdt_evt"})
                                            .Define("trigger_w", get_unb_prescale, {"geo_lat"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_over_unb_bdt", "HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev", "trigger_w");

    auto h_trigger_efficiency_accepted_het_let_bdt = mc_file ?
                                            fr.Filter("trigger_efficiency_preselection==1 && (trigger_efficiency_preselection_is_het==1 || trigger_efficiency_preselection_is_let==1)")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Filter(bdt_cut, {"bdt_evt"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_let_bdt", "LET + HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev")

                                            :

                                            fr.Filter("trigger_efficiency_preselection==1 && (trigger_efficiency_preselection_is_het==1 || trigger_efficiency_preselection_is_let==1)")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("trigger_w", [&get_let_prescale](const bool is_het, const double geo_lat) -> double { if (is_het) return get_let_prescale(geo_lat); else return 1;}, {"trigger_efficiency_preselection_is_het", "geo_lat"})
                                            .Filter(bdt_cut, {"bdt_evt"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_let_bdt", "LET + HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev", "trigger_w");

    auto h_trigger_efficiency_accepted_het_unb_bdt = mc_file ? 
                                            fr.Filter("trigger_efficiency_preselection==1 && (trigger_efficiency_preselection_is_het==1 || trigger_efficiency_preselection_is_unb==1)")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Filter(bdt_cut, {"bdt_evt"})
                                            .Histo1D({"h_trigger_efficiency_accepted_het_unb_bdt", "LET + HET Trigger", energy_nbins, &energy_binning[0]}, "corr_energy_gev")

                                            :

                                            fr.Filter("trigger_efficiency_preselection==1 && (trigger_efficiency_preselection_is_het==1 || trigger_efficiency_preselection_is_unb==1)")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Define("trigger_w", [&get_unb_prescale](const bool is_het, const double geo_lat) -> double { if (is_het) return get_unb_prescale(geo_lat); else return 1;}, {"trigger_efficiency_preselection_is_het", "geo_lat"})
                                            .Filter(bdt_cut, {"bdt_evt"})
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
    
    auto h_maxrms_efficiency_accepted_bdt = fr.Filter("maxrms_efficiency_preselection==1 && maxrms_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Filter(bdt_cut, {"bdt_evt"})
                                            .Histo1D({"h_maxrms_efficiency_accepted_bdt", "HET Trigger + MaxRMS Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
    auto h_maxrms_efficiency_total_bdt = fr.Filter("maxrms_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Filter(bdt_cut, {"bdt_evt"})
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

    auto h_nbarlayer13_efficiency_accepted_bdt = fr.Filter("nbarlayer13_efficiency_preselection==1 && nbarlayer13_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Filter(bdt_cut, {"bdt_evt"})
                                            .Histo1D({"h_nbarlayer13_efficiency_accepted_bdt", "HET Trigger + nbarlayer13 Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
    auto h_nbarlayer13_efficiency_total_bdt = fr.Filter("nbarlayer13_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Filter(bdt_cut, {"bdt_evt"})
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

    auto h_maxrms_and_nbarlayer13_efficiency_accepted_bdt = fr.Filter("maxrms_and_nbarlayer13_efficiency_preselection==1 && maxrms_and_nbarlayer13_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Filter(bdt_cut, {"bdt_evt"})
                                            .Histo1D({"h_maxrms_and_nbarlayer13_efficiency_accepted_bdt", "HET Trigger + maxrms & nbarlayer13 Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
    auto h_maxrms_and_nbarlayer13_efficiency_total_bdt = fr.Filter("maxrms_and_nbarlayer13_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Filter(bdt_cut, {"bdt_evt"})
                                            .Histo1D({"h_maxrms_and_nbarlayer13_efficiency_total_bdt", "HET Trigger + maxrms & nbarlayer13", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

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
    
    auto h_track_efficiency_accepted_bdt = fr.Filter("track_efficiency_preselection==1 && track_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Filter(bdt_cut, {"bdt_evt"})
                                            .Histo1D({"h_track_efficiency_accepted_bdt", "HET Trigger + Track Selection Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
    auto h_track_efficiency_total_bdt = fr.Filter("track_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Filter(bdt_cut, {"bdt_evt"})
                                            .Histo1D({"h_track_efficiency_total_bdt", "HET Trigger + Track Selection", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

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

    auto h_psdstkmatch_efficiency_accepted_bdt = fr.Filter("psdstkmatch_efficiency_preselection==1 && psdstkmatch_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Filter(bdt_cut, {"bdt_evt"})
                                            .Histo1D({"h_psdstkmatch_efficiency_accepted_bdt", "HET Trigger + PSD-STK Match Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
    auto h_psdstkmatch_efficiency_total_bdt = fr.Filter("psdstkmatch_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Filter(bdt_cut, {"bdt_evt"})
                                            .Histo1D({"h_psdstkmatch_efficiency_total_bdt", "HET Trigger + PSD-STK Match", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    // PSD-STK match histos within PSD fiducial volume
    auto h_psdstkmatch_efficiency_psd_fvolume_accepted_tight_xtrl = fr.Filter("psdstkmatch_efficiency_preselection==1 && evtfilter_psd_stk_match_cut_psd_fiducial_volume==1 && psdstkmatch_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_psdstkmatch_efficiency_psd_fvolume_accepted_tight_xtrl", "HET Trigger + PSD-STK Match Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
    auto h_psdstkmatch_efficiency_psd_fvolume_total_tight_xtrl = fr.Filter("psdstkmatch_efficiency_preselection==1 && evtfilter_psd_stk_match_cut_psd_fiducial_volume==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                            .Histo1D({"h_psdstkmatch_efficiency_psd_fvolume_total_tight_xtrl", "HET Trigger + PSD-STK Match", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    auto h_psdstkmatch_efficiency_psd_fvolume_accepted_bdt = fr.Filter("psdstkmatch_efficiency_preselection==1 && evtfilter_psd_stk_match_cut_psd_fiducial_volume==1 && psdstkmatch_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Filter(bdt_cut, {"bdt_evt"})
                                            .Histo1D({"h_psdstkmatch_efficiency_psd_fvolume_accepted_bdt", "HET Trigger + PSD-STK Match Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
    auto h_psdstkmatch_efficiency_psd_fvolume_total_bdt = fr.Filter("psdstkmatch_efficiency_preselection==1 && evtfilter_psd_stk_match_cut_psd_fiducial_volume==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Filter(bdt_cut, {"bdt_evt"})
                                            .Histo1D({"h_psdstkmatch_efficiency_psd_fvolume_total_bdt", "HET Trigger + PSD-STK Match", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

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
    
    auto h_psdcharge_efficiency_accepted_bdt = fr.Filter("psdcharge_efficiency_preselection==1 && psdcharge_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Filter(bdt_cut, {"bdt_evt"})
                                            .Histo1D({"h_psdcharge_efficiency_accepted_bdt", "HET Trigger + PSD Charge Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
    auto h_psdcharge_efficiency_total_bdt = fr.Filter("psdcharge_efficiency_preselection==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("bdt_evt", compute_bdt, {"rmsLayer", "sumRms", "fracLayer", "fracLast_13", "corr_energy_gev", "BGOrec_trajectoryDirection2D"})
                                            .Filter(bdt_cut, {"bdt_evt"})
                                            .Histo1D({"h_psdcharge_efficiency_total_bdt", "HET Trigger + PSD Charge Match", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    // Trigger histos
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
    // MaxRMS histos
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

    // nbarlayer13 histos
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

    // MaxRMS and nbarlayer13 histos
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

    // Track Selection histos
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

    // PSD-STK match histos
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

    // PSD-STK match histos within PSD fiducial volume
    auto h_psdstkmatch_efficiency_psd_fvolume_accepted_loose_xtrl = fr.Filter("psdstkmatch_efficiency_preselection==1 && evtfilter_psd_stk_match_cut_psd_fiducial_volume==1 && psdstkmatch_efficiency_preselection_accepted==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_psdstkmatch_efficiency_psd_fvolume_accepted_loose_xtrl", "HET Trigger + PSD-STK Match Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
    auto h_psdstkmatch_efficiency_psd_fvolume_total_loose_xtrl = fr.Filter("psdstkmatch_efficiency_preselection==1 && evtfilter_psd_stk_match_cut_psd_fiducial_volume==1")
                                            .Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Define("energy_gev", "energy * 0.001")
                                            .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                            .Filter(xtrl_loose_cut, {"energy_gev", "xtrl_evt"})
                                            .Histo1D({"h_psdstkmatch_efficiency_psd_fvolume_total_loose_xtrl", "HET Trigger + PSD-STK Match", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

    // PSD charge histos
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
    auto n_events_within_psd_fvolume    = *(fr.Filter("evtfilter_track_selection_cut==true").Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true").Count());
    auto n_events_outside_psd_fvolume   = *(fr.Filter("evtfilter_track_selection_cut==true").Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false").Count());

    std::cout << "\n\n*** PSD Fiducial Volume ***\n\n";
    std::cout << "Number of events after BGO and STK selection: " << n_events_presel_psd_fvolume;
    std::cout << "\nNumber of events within the PSD fiducial volume: " << n_events_within_psd_fvolume;
    std::cout << "\nNumber of events outside the PSD fiducial volume: " << n_events_outside_psd_fvolume;
    std::cout << "\n\n***************************\n\n";

    // PSD fiducial volume distance plots
    auto h_psd_stk_match_distance_x_20_100 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                .Histo1D({"h_psd_stk_match_distance_x_20_100", "PSD/STK match - X distance", 600, -300, 300}, "SPD_STK_match_X_distance");

    auto h_psd_stk_match_distance_x_100_250 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                .Histo1D({"h_psd_stk_match_distance_x_100_250", "PSD/STK match - X distance", 600, -300, 300}, "SPD_STK_match_X_distance");

    auto h_psd_stk_match_distance_x_250_500 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                .Histo1D({"h_psd_stk_match_distance_x_250_500", "PSD/STK match - X distance",  600, -300, 300}, "SPD_STK_match_X_distance");

    auto h_psd_stk_match_distance_x_500_1000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                .Histo1D({"h_psd_stk_match_distance_x_500_1000", "PSD/STK match - X distance",  600, -300, 300}, "SPD_STK_match_X_distance");
    
    auto h_psd_stk_match_distance_x_1000_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                .Histo1D({"h_psd_stk_match_distance_x_1000_3000", "PSD/STK match - X distance", 600, -300, 300}, "SPD_STK_match_X_distance");

    auto h_psd_stk_match_distance_x_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=3000")
                                                .Histo1D({"h_psd_stk_match_distance_x_3000", "PSD/STK match - X distance", 600, -300, 300}, "SPD_STK_match_X_distance");

    auto h_psd_stk_match_distance_y_20_100 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                .Histo1D({"h_psd_stk_match_distance_y_20_100", "PSD/STK match - X distance", 600, -300, 300}, "SPD_STK_match_Y_distance");

    auto h_psd_stk_match_distance_y_100_250 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                .Histo1D({"h_psd_stk_match_distance_y_100_250", "PSD/STK match - Y distance", 600, -300, 300}, "SPD_STK_match_Y_distance");

    auto h_psd_stk_match_distance_y_250_500 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                .Histo1D({"h_psd_stk_match_distance_y_250_500", "PSD/STK match - Y distance",  600, -300, 300}, "SPD_STK_match_Y_distance");

    auto h_psd_stk_match_distance_y_500_1000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                .Histo1D({"h_psd_stk_match_distance_y_500_1000", "PSD/STK match - Y distance",  600, -300, 300}, "SPD_STK_match_Y_distance");
    
    auto h_psd_stk_match_distance_y_1000_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                .Histo1D({"h_psd_stk_match_distance_y_1000_3000", "PSD/STK match - Y distance", 600, -300, 300}, "SPD_STK_match_Y_distance");

    auto h_psd_stk_match_distance_y_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=3000")
                                                .Histo1D({"h_psd_stk_match_distance_y_3000", "PSD/STK match - Y distance", 600, -300, 300}, "SPD_STK_match_Y_distance");

    // PSD fiducial volume distance plots within PSD fiducial volume
    auto h_psd_stk_match_distance_x_within_psd_fvolume_20_100 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                .Histo1D({"h_psd_stk_match_distance_x_within_psd_fvolume_20_100", "PSD/STK match - X distance", 600, -300, 300}, "SPD_STK_match_X_distance_fiducial_volume");

    auto h_psd_stk_match_distance_x_within_psd_fvolume_100_250 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                .Histo1D({"h_psd_stk_match_distance_x_within_psd_fvolume_100_250", "PSD/STK match - X distance", 600, -300, 300}, "SPD_STK_match_X_distance_fiducial_volume");

    auto h_psd_stk_match_distance_x_within_psd_fvolume_250_500 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                .Histo1D({"h_psd_stk_match_distance_x_within_psd_fvolume_250_500", "PSD/STK match - X distance",  600, -300, 300}, "SPD_STK_match_X_distance_fiducial_volume");

    auto h_psd_stk_match_distance_x_within_psd_fvolume_500_1000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                .Histo1D({"h_psd_stk_match_distance_x_within_psd_fvolume_500_1000", "PSD/STK match - X distance",  600, -300, 300}, "SPD_STK_match_X_distance_fiducial_volume");
    
    auto h_psd_stk_match_distance_x_within_psd_fvolume_1000_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                .Histo1D({"h_psd_stk_match_distance_x_within_psd_fvolume_1000_3000", "PSD/STK match - X distance", 600, -300, 300}, "SPD_STK_match_X_distance_fiducial_volume");

    auto h_psd_stk_match_distance_x_within_psd_fvolume_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=3000")
                                                .Histo1D({"h_psd_stk_match_distance_x_within_psd_fvolume_3000", "PSD/STK match - X distance", 600, -300, 300}, "SPD_STK_match_X_distance_fiducial_volume");

    auto h_psd_stk_match_distance_y_within_psd_fvolume_20_100 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                .Histo1D({"h_psd_stk_match_distance_y_within_psd_fvolume_20_100", "PSD/STK match - X distance", 600, -300, 300}, "SPD_STK_match_Y_distance_fiducial_volume");

    auto h_psd_stk_match_distance_y_within_psd_fvolume_100_250 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                .Histo1D({"h_psd_stk_match_distance_y_within_psd_fvolume_100_250", "PSD/STK match - Y distance", 600, -300, 300}, "SPD_STK_match_Y_distance_fiducial_volume");

    auto h_psd_stk_match_distance_y_within_psd_fvolume_250_500 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                .Histo1D({"h_psd_stk_match_distance_y_within_psd_fvolume_250_500", "PSD/STK match - Y distance",  600, -300, 300}, "SPD_STK_match_Y_distance_fiducial_volume");

    auto h_psd_stk_match_distance_y_within_psd_fvolume_500_1000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                .Histo1D({"h_psd_stk_match_distance_y_within_psd_fvolume_500_1000", "PSD/STK match - Y distance",  600, -300, 300}, "SPD_STK_match_Y_distance_fiducial_volume");
    
    auto h_psd_stk_match_distance_y_within_psd_fvolume_1000_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                .Histo1D({"h_psd_stk_match_distance_y_within_psd_fvolume_1000_3000", "PSD/STK match - Y distance", 600, -300, 300}, "SPD_STK_match_Y_distance_fiducial_volume");

    auto h_psd_stk_match_distance_y_within_psd_fvolume_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=3000")
                                                .Histo1D({"h_psd_stk_match_distance_y_within_psd_fvolume_3000", "PSD/STK match - Y distance", 600, -300, 300}, "SPD_STK_match_Y_distance_fiducial_volume");

    // PSD fiducial volume distance plots outside PSD fiducial volume
    auto h_psd_stk_match_distance_x_outside_psd_fvolume_20_100 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                .Histo1D({"h_psd_stk_match_distance_x_outside_psd_fvolume_20_100", "PSD/STK match - X distance", 600, -300, 300}, "SPD_STK_match_X_distance");

    auto h_psd_stk_match_distance_x_outside_psd_fvolume_100_250 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                .Histo1D({"h_psd_stk_match_distance_x_outside_psd_fvolume_100_250", "PSD/STK match - X distance", 600, -300, 300}, "SPD_STK_match_X_distance");

    auto h_psd_stk_match_distance_x_outside_psd_fvolume_250_500 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                .Histo1D({"h_psd_stk_match_distance_x_outside_psd_fvolume_250_500", "PSD/STK match - X distance",  600, -300, 300}, "SPD_STK_match_X_distance");

    auto h_psd_stk_match_distance_x_outside_psd_fvolume_500_1000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                .Histo1D({"h_psd_stk_match_distance_x_outside_psd_fvolume_500_1000", "PSD/STK match - X distance",  600, -300, 300}, "SPD_STK_match_X_distance");
    
    auto h_psd_stk_match_distance_x_outside_psd_fvolume_1000_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                .Histo1D({"h_psd_stk_match_distance_x_outside_psd_fvolume_1000_3000", "PSD/STK match - X distance", 600, -300, 300}, "SPD_STK_match_X_distance");

    auto h_psd_stk_match_distance_x_outside_psd_fvolume_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=3000")
                                                .Histo1D({"h_psd_stk_match_distance_x_outside_psd_fvolume_3000", "PSD/STK match - X distance", 600, -300, 300}, "SPD_STK_match_X_distance");

    auto h_psd_stk_match_distance_y_outside_psd_fvolume_20_100 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                .Histo1D({"h_psd_stk_match_distance_y_outside_psd_fvolume_20_100", "PSD/STK match - X distance", 600, -300, 300}, "SPD_STK_match_Y_distance");

    auto h_psd_stk_match_distance_y_outside_psd_fvolume_100_250 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                .Histo1D({"h_psd_stk_match_distance_y_outside_psd_fvolume_100_250", "PSD/STK match - Y distance", 600, -300, 300}, "SPD_STK_match_Y_distance");

    auto h_psd_stk_match_distance_y_outside_psd_fvolume_250_500 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                .Histo1D({"h_psd_stk_match_distance_y_outside_psd_fvolume_250_500", "PSD/STK match - Y distance",  600, -300, 300}, "SPD_STK_match_Y_distance");

    auto h_psd_stk_match_distance_y_outside_psd_fvolume_500_1000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                .Histo1D({"h_psd_stk_match_distance_y_outside_psd_fvolume_500_1000", "PSD/STK match - Y distance",  600, -300, 300}, "SPD_STK_match_Y_distance");
    
    auto h_psd_stk_match_distance_y_outside_psd_fvolume_1000_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                .Histo1D({"h_psd_stk_match_distance_y_outside_psd_fvolume_1000_3000", "PSD/STK match - Y distance", 600, -300, 300}, "SPD_STK_match_Y_distance");

    auto h_psd_stk_match_distance_y_outside_psd_fvolume_3000 = fr.Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev>=3000")
                                                .Histo1D({"h_psd_stk_match_distance_y_outside_psd_fvolume_3000", "PSD/STK match - Y distance", 600, -300, 300}, "SPD_STK_match_Y_distance");

    // STK charge plots within PSD fiducial volume and withoud PSD charge cut
    auto h_stk_charge_psd_fvolume_no_psd_cut_20_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_20_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_psd_fvolume_no_psd_cut_100_250 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_100_250", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_psd_fvolume_no_psd_cut_250_500 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_250_500", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_no_psd_cut_500_1000 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_500_1000", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_no_psd_cut_1000_3000 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_1000_3000", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_no_psd_cut_3000 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=3000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_no_psd_cut_3000", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    // STK charge plots within PSD fiducial volume and with PSD charge cut
    auto h_stk_charge_psd_fvolume_psd_cut_20_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_20_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_psd_fvolume_psd_cut_100_250 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_100_250", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_psd_fvolume_psd_cut_250_500 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_250_500", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_psd_cut_500_1000 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_500_1000", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_psd_cut_1000_3000 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_1000_3000", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_fvolume_psd_cut_3000 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==true")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=3000")
                                                        .Histo2D({"h_stk_charge_psd_fvolume_psd_cut_3000", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    // STK charge plots outside PSD fiducial volume and withoud PSD charge cut
    auto h_stk_charge_no_psd_cut_20_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                        .Histo2D({"h_stk_charge_no_psd_cut_20_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_no_psd_cut_100_250 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                        .Histo2D({"h_stk_charge_no_psd_cut_100_250", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_no_psd_cut_250_500 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                        .Histo2D({"h_stk_charge_no_psd_cut_250_500", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_no_psd_cut_500_1000 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                        .Histo2D({"h_stk_charge_no_psd_cut_500_1000", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_no_psd_cut_1000_3000 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                        .Histo2D({"h_stk_charge_no_psd_cut_1000_3000", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_no_psd_cut_3000 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=3000")
                                                        .Histo2D({"h_stk_charge_no_psd_cut_3000", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    // STK charge plots outside PSD fiducial volume and with PSD charge cut
    auto h_stk_charge_psd_cut_20_100 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                                        .Histo2D({"h_stk_charge_psd_cut_20_100", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_psd_cut_100_250 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                                        .Histo2D({"h_stk_charge_psd_cut_100_250", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");
    
    auto h_stk_charge_psd_cut_250_500 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                                        .Histo2D({"h_stk_charge_psd_cut_250_500", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_cut_500_1000 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                                        .Histo2D({"h_stk_charge_psd_cut_500_1000", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_cut_1000_3000 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                                        .Histo2D({"h_stk_charge_psd_cut_1000_3000", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

    auto h_stk_charge_psd_cut_3000 = fr.Filter("evtfilter_psd_stk_match_cut==true")
                                                        .Filter("evtfilter_psd_stk_match_cut_psd_fiducial_volume==false")
                                                        .Filter("evtfilter_psd_charge_cut==true")
                                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev>=3000")
                                                        .Histo2D({"h_stk_charge_psd_cut_3000", "STK charge; STK X Charge; STK Y Charge", 500, 0, 150, 500, 0, 150}, "STK_chargeX", "STK_chargeY");

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
    h_track_efficiency_accepted_tight_xtrl                              ->Write();
    h_track_efficiency_total_tight_xtrl                                 ->Write();
    h_track_efficiency_accepted_bdt                                     ->Write();
    h_track_efficiency_total_bdt                                        ->Write();
    h_psdstkmatch_efficiency_accepted_tight_xtrl                        ->Write();
    h_psdstkmatch_efficiency_total_tight_xtrl                           ->Write();
    h_psdstkmatch_efficiency_accepted_bdt                               ->Write();
    h_psdstkmatch_efficiency_total_bdt                                  ->Write();
    h_psdstkmatch_efficiency_psd_fvolume_accepted_tight_xtrl            ->Write();
    h_psdstkmatch_efficiency_psd_fvolume_total_tight_xtrl               ->Write();
    h_psdstkmatch_efficiency_psd_fvolume_accepted_bdt                   ->Write();
    h_psdstkmatch_efficiency_psd_fvolume_total_bdt                      ->Write();
    h_psdcharge_efficiency_accepted_tight_xtrl                          ->Write();
    h_psdcharge_efficiency_total_tight_xtrl                             ->Write();
    h_psdcharge_efficiency_accepted_bdt                                 ->Write();
    h_psdcharge_efficiency_total_bdt                                    ->Write();
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
    h_psdstkmatch_efficiency_accepted_loose_xtrl                        ->Write();
    h_psdstkmatch_efficiency_total_loose_xtrl                           ->Write();
    h_psdstkmatch_efficiency_psd_fvolume_accepted_loose_xtrl            ->Write();
    h_psdstkmatch_efficiency_psd_fvolume_total_loose_xtrl               ->Write();
    h_psdcharge_efficiency_accepted_loose_xtrl                          ->Write();
    h_psdcharge_efficiency_total_loose_xtrl                             ->Write();

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

    output_file->mkdir("stk_charge");
    output_file->cd("stk_charge");

    h_stk_charge_psd_fvolume_no_psd_cut_20_100                          ->Write();
    h_stk_charge_psd_fvolume_no_psd_cut_100_250                         ->Write();
    h_stk_charge_psd_fvolume_no_psd_cut_250_500                         ->Write();
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000                        ->Write();
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000                       ->Write();
    h_stk_charge_psd_fvolume_no_psd_cut_3000                            ->Write();

    h_stk_charge_psd_fvolume_psd_cut_20_100                             ->Write();
    h_stk_charge_psd_fvolume_psd_cut_100_250                            ->Write();
    h_stk_charge_psd_fvolume_psd_cut_250_500                            ->Write();
    h_stk_charge_psd_fvolume_psd_cut_500_1000                           ->Write();
    h_stk_charge_psd_fvolume_psd_cut_1000_3000                          ->Write();
    h_stk_charge_psd_fvolume_psd_cut_3000                               ->Write();

    h_stk_charge_no_psd_cut_20_100                                      ->Write();
    h_stk_charge_no_psd_cut_100_250                                     ->Write();
    h_stk_charge_no_psd_cut_250_500                                     ->Write();
    h_stk_charge_no_psd_cut_500_1000                                    ->Write();
    h_stk_charge_no_psd_cut_1000_3000                                   ->Write();
    h_stk_charge_no_psd_cut_3000                                        ->Write();

    h_stk_charge_psd_cut_20_100                                         ->Write();
    h_stk_charge_psd_cut_100_250                                        ->Write();
    h_stk_charge_psd_cut_250_500                                        ->Write();
    h_stk_charge_psd_cut_500_1000                                       ->Write();
    h_stk_charge_psd_cut_1000_3000                                      ->Write();
    h_stk_charge_psd_cut_3000                                           ->Write();

    output_file->Close();
}   