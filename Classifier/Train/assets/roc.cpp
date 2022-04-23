#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include <memory>
#include <tuple>
#include <vector>
#include <iostream>

struct bdt_vars {
    /*
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
    */

    float rmslayer_1 {0};
    float rmslayer_2 {0};
    float rmslayer_3 {0};
    float rmslayer_4 {0};
    float rmslayer_5 {0};
    float rmslayer_6 {0};
    float rmslayer_7 {0};
    float rmslayer_8 {0};
    float rmslayer_9 {0};
    float rmslayer_10 {0};
    float rmslayer_11 {0};
    float rmslayer_12 {0};
    float rmslayer_13 {0};
    float rmslayer_14 {0};
    float fraclayer_1 {0};
    float fraclayer_2 {0};
    float fraclayer_3 {0};
    float fraclayer_4 {0};
    float fraclayer_5 {0};
    float fraclayer_6 {0};
    float fraclayer_7 {0};
    float fraclayer_8 {0};
    float fraclayer_9 {0};
    float fraclayer_10 {0};
    float fraclayer_11 {0};
    float fraclayer_12 {0};
    float fraclayer_13 {0};
    float fraclayer_14 {0};
    float sumRms {0};
    float fracLast {0};
    float xtrl {0};

};

struct roc_vars {
    double proton_rejection_xtrl {0};
    double proton_rejection_bdt {0};
    double electron_efficiency_xtrl {0};
    double proton_efficiency_xtrl {0};
    double electron_efficiency_bdt {0};
    double proton_efficiency_bdt {0};
};

inline const char* get_tree_name(const char* file) {
    TFile* input_file = TFile::Open(file, "READ");
    if (!input_file->IsOpen()) {
        std::cerr << "\n\nError reading input file [" << file << "]\n\n";
        exit(100);
    }
    std::string tree_name;
    for (TObject* keyAsObject : *input_file->GetListOfKeys()) {
        auto key = dynamic_cast<TKey*>(keyAsObject);
        if (!strcmp(key->GetClassName(), "TTree"))
            tree_name = static_cast<std::string>(key->GetName());
    }
    input_file->Close();
    return tree_name.c_str();
}

inline void load_bdt_vars(std::shared_ptr<TMVA::Reader> tmva_reader,  bdt_vars &my_bdt_vars) {

    /*
    tmva_reader->AddVariable("rmslayer_norm_1", &my_bdt_vars.rmslayer_norm_1);
    tmva_reader->AddVariable("rmslayer_norm_2", &my_bdt_vars.rmslayer_norm_2);
    tmva_reader->AddVariable("rmslayer_norm_3", &my_bdt_vars.rmslayer_norm_3);
    tmva_reader->AddVariable("rmslayer_norm_4", &my_bdt_vars.rmslayer_norm_4);
    tmva_reader->AddVariable("rmslayer_norm_5", &my_bdt_vars.rmslayer_norm_5);
    tmva_reader->AddVariable("rmslayer_norm_6", &my_bdt_vars.rmslayer_norm_6);
    tmva_reader->AddVariable("rmslayer_norm_7", &my_bdt_vars.rmslayer_norm_7);
    tmva_reader->AddVariable("rmslayer_norm_8", &my_bdt_vars.rmslayer_norm_8);
    tmva_reader->AddVariable("rmslayer_norm_9", &my_bdt_vars.rmslayer_norm_9);
    tmva_reader->AddVariable("rmslayer_norm_10", &my_bdt_vars.rmslayer_norm_10);
    tmva_reader->AddVariable("rmslayer_norm_11", &my_bdt_vars.rmslayer_norm_11);
    tmva_reader->AddVariable("rmslayer_norm_12", &my_bdt_vars.rmslayer_norm_12);
    tmva_reader->AddVariable("rmslayer_norm_13", &my_bdt_vars.rmslayer_norm_13);
    tmva_reader->AddVariable("rmslayer_norm_14", &my_bdt_vars.rmslayer_norm_14);
    tmva_reader->AddVariable("fraclayer_norm_1", &my_bdt_vars.fraclayer_norm_1);
    tmva_reader->AddVariable("fraclayer_norm_2", &my_bdt_vars.fraclayer_norm_2);
    tmva_reader->AddVariable("fraclayer_norm_3", &my_bdt_vars.fraclayer_norm_3);
    tmva_reader->AddVariable("fraclayer_norm_4", &my_bdt_vars.fraclayer_norm_4);
    tmva_reader->AddVariable("fraclayer_norm_5", &my_bdt_vars.fraclayer_norm_5);
    tmva_reader->AddVariable("fraclayer_norm_6", &my_bdt_vars.fraclayer_norm_6);
    tmva_reader->AddVariable("fraclayer_norm_7", &my_bdt_vars.fraclayer_norm_7);
    tmva_reader->AddVariable("fraclayer_norm_8", &my_bdt_vars.fraclayer_norm_8);
    tmva_reader->AddVariable("fraclayer_norm_9", &my_bdt_vars.fraclayer_norm_9);
    tmva_reader->AddVariable("fraclayer_norm_10", &my_bdt_vars.fraclayer_norm_10);
    tmva_reader->AddVariable("fraclayer_norm_11", &my_bdt_vars.fraclayer_norm_11);
    tmva_reader->AddVariable("fraclayer_norm_12", &my_bdt_vars.fraclayer_norm_12);
    tmva_reader->AddVariable("fraclayer_norm_13", &my_bdt_vars.fraclayer_norm_13);
    tmva_reader->AddVariable("fraclayer_norm_14", &my_bdt_vars.fraclayer_norm_14);
    tmva_reader->AddVariable("sumrms_norm", &my_bdt_vars.sumrms_norm);
    tmva_reader->AddVariable("fraclastlayer_norm", &my_bdt_vars.fraclastlayer_norm);
    tmva_reader->AddVariable("xtrl_norm", &my_bdt_vars.xtrl_norm);
    tmva_reader->AddSpectator("xtrl", &my_bdt_vars.xtrl);
    */

    tmva_reader->AddVariable("rmslayer_1", &my_bdt_vars.rmslayer_1);
    tmva_reader->AddVariable("rmslayer_2", &my_bdt_vars.rmslayer_2);
    tmva_reader->AddVariable("rmslayer_3", &my_bdt_vars.rmslayer_3);
    tmva_reader->AddVariable("rmslayer_4", &my_bdt_vars.rmslayer_4);
    tmva_reader->AddVariable("rmslayer_5", &my_bdt_vars.rmslayer_5);
    tmva_reader->AddVariable("rmslayer_6", &my_bdt_vars.rmslayer_6);
    tmva_reader->AddVariable("rmslayer_7", &my_bdt_vars.rmslayer_7);
    tmva_reader->AddVariable("rmslayer_8", &my_bdt_vars.rmslayer_8);
    tmva_reader->AddVariable("rmslayer_9", &my_bdt_vars.rmslayer_9);
    tmva_reader->AddVariable("rmslayer_10", &my_bdt_vars.rmslayer_10);
    tmva_reader->AddVariable("rmslayer_11", &my_bdt_vars.rmslayer_11);
    tmva_reader->AddVariable("rmslayer_12", &my_bdt_vars.rmslayer_12);
    tmva_reader->AddVariable("rmslayer_13", &my_bdt_vars.rmslayer_13);
    tmva_reader->AddVariable("rmslayer_14", &my_bdt_vars.rmslayer_14);
    tmva_reader->AddVariable("fraclayer_1", &my_bdt_vars.fraclayer_1);
    tmva_reader->AddVariable("fraclayer_2", &my_bdt_vars.fraclayer_2);
    tmva_reader->AddVariable("fraclayer_3", &my_bdt_vars.fraclayer_3);
    tmva_reader->AddVariable("fraclayer_4", &my_bdt_vars.fraclayer_4);
    tmva_reader->AddVariable("fraclayer_5", &my_bdt_vars.fraclayer_5);
    tmva_reader->AddVariable("fraclayer_6", &my_bdt_vars.fraclayer_6);
    tmva_reader->AddVariable("fraclayer_7", &my_bdt_vars.fraclayer_7);
    tmva_reader->AddVariable("fraclayer_8", &my_bdt_vars.fraclayer_8);
    tmva_reader->AddVariable("fraclayer_9", &my_bdt_vars.fraclayer_9);
    tmva_reader->AddVariable("fraclayer_10", &my_bdt_vars.fraclayer_10);
    tmva_reader->AddVariable("fraclayer_11", &my_bdt_vars.fraclayer_11);
    tmva_reader->AddVariable("fraclayer_12", &my_bdt_vars.fraclayer_12);
    tmva_reader->AddVariable("fraclayer_13", &my_bdt_vars.fraclayer_13);
    tmva_reader->AddVariable("fraclayer_14", &my_bdt_vars.fraclayer_14);
    tmva_reader->AddVariable("sumRms", &my_bdt_vars.sumRms);
    tmva_reader->AddVariable("fracLast", &my_bdt_vars.fracLast);
    tmva_reader->AddVariable("xtrl", &my_bdt_vars.xtrl);
}

inline void branch_roc_tree(std::shared_ptr<TTree> roc_tree, roc_vars &my_roc_vars) {
    roc_tree->Branch("proton_rejection_xtrl", &my_roc_vars.proton_rejection_xtrl, "proton_rejection_xtrl/D");
    roc_tree->Branch("proton_rejection_bdt", &my_roc_vars.proton_rejection_bdt, "proton_rejection_bdt/D");
    roc_tree->Branch("electron_efficiency_xtrl", &my_roc_vars.electron_efficiency_xtrl, "electron_efficiency_xtrl/D");
    roc_tree->Branch("proton_efficiency_xtrl", &my_roc_vars.proton_efficiency_xtrl, "proton_efficiency_xtrl/D");
    roc_tree->Branch("electron_efficiency_bdt", &my_roc_vars.electron_efficiency_bdt, "electron_efficiency_bdt/D");
    roc_tree->Branch("proton_efficiency_bdt", &my_roc_vars.proton_efficiency_bdt, "proton_efficiency_bdt/D");
}

inline bool in_energy_range (const double energy, const double emin, const double emax) {
    return energy >= emin && energy <= emax ? true : false;
}

std::tuple<double, double> get_xtrl_efficiencies(
    std::shared_ptr<TTreeReader> electron_reader, 
    std::shared_ptr<TTreeReader> proton_reader, 
    const double & xtrl_cut_value, 
    const double &emin, 
    const double &emax);

std::tuple<double, double> get_bdt_efficiencies(
    std::shared_ptr<TTreeReader> electron_reader, 
    std::shared_ptr<TTreeReader> proton_reader, 
    const double &bdt_cut_value, 
    const double &emin, 
    const double &emax,
    const char* bdt_weights_file,
    bdt_vars &my_bdt_vars,
    std::shared_ptr<TMVA::Reader> tmva_reader);

void roc(
    const char* electron_mc_file,
    const char* proton_mc_file,
    const char* bdt_weights_file,
    const char* outfile,
    const double emin,
    const double emax,
    const unsigned int samples,
    const double xtrl_min,
    const double xtrl_max,
    const double bdt_min,
    const double bdt_max) {

        // Read electron MC file and extract the reader
        TFile *electron_file = TFile::Open(electron_mc_file, "READ");
        if (!electron_file->IsOpen()) {
            std::cerr << "\n\nError reading input file [" << electron_mc_file << "\n\n";
            exit(100);
        }
        std::shared_ptr<TTreeReader> electron_reader = std::make_shared<TTreeReader>(get_tree_name(electron_mc_file), electron_file);

        // Read proton MC file and extract the reader
        TFile *proton_file = TFile::Open(proton_mc_file, "READ");
        if (!proton_file->IsOpen()) {
            std::cerr << "\n\nError reading input file [" << proton_mc_file << "\n\n";
            exit(100);
        }
        std::shared_ptr<TTreeReader> proton_reader = std::make_shared<TTreeReader>(get_tree_name(proton_mc_file), proton_file); 
        
        // Initialize TMVA reader
        TMVA::Tools::Instance();
        std::shared_ptr<TMVA::Reader> tmva_reader = std::make_shared<TMVA::Reader>();
        bdt_vars my_bdt_vars;
        load_bdt_vars(tmva_reader, my_bdt_vars);
        tmva_reader->BookMVA("BDT method", bdt_weights_file);

        // Write output file and branch the final Tree
        TFile *output = TFile::Open(outfile, "RECREATE");
        if (!output->IsOpen()) {
            std::cerr << "\n\nError writing output file [" << outfile << "\n\n";
            exit(100);
        }
        std::shared_ptr<TTree> roc_tree = std::make_shared<TTree>("roc_tree", "Proton Rejection vs Signal Selection Tree values");
        roc_vars my_roc_vars;
        branch_roc_tree(roc_tree, my_roc_vars);

        // Define vectors for TGraphs
        std::vector<double> xtrl_signal_eff (samples, 0);
        std::vector<double> bdt_signal_eff (samples, 0);
        std::vector<double> xtrl_proton_rejection (samples, 0);
        std::vector<double> bdt_proton_rejection (samples, 0);
        std::vector<double> xtrl_proton_rejection_power (samples, 0);
        std::vector<double> bdt_proton_rejection_power (samples, 0);
        std::vector<double> xtrl_signal_eff_over_proton_eff (samples, 0);
        std::vector<double> bdt_signal_eff_over_proton_eff (samples, 0);

        // Loop over samples
        for (unsigned int ev_sample {0}; ev_sample<samples; ++ev_sample) {
            std::cout << "\nNumber of processed steps [" << ev_sample+1 << "]";

            auto xtrl_cut_value = xtrl_min + ev_sample * (xtrl_max-xtrl_min)/samples;
            auto bdt_cut_value = bdt_min + ev_sample * (bdt_max-bdt_min)/samples;
            
            std::tie(my_roc_vars.electron_efficiency_xtrl, my_roc_vars.proton_efficiency_xtrl) = get_xtrl_efficiencies(electron_reader, proton_reader, xtrl_cut_value, emin, emax);
            std::tie(my_roc_vars.electron_efficiency_bdt, my_roc_vars.proton_efficiency_bdt) = get_bdt_efficiencies(electron_reader, proton_reader, bdt_cut_value, emin, emax, bdt_weights_file, my_bdt_vars, tmva_reader);
           
            my_roc_vars.proton_rejection_xtrl = 1 - my_roc_vars.proton_efficiency_xtrl;
            my_roc_vars.proton_rejection_bdt = 1 - my_roc_vars.proton_efficiency_bdt;

            xtrl_signal_eff[ev_sample] = my_roc_vars.electron_efficiency_xtrl;
            bdt_signal_eff[ev_sample] = my_roc_vars.electron_efficiency_bdt;
            xtrl_proton_rejection[ev_sample] = my_roc_vars.proton_rejection_xtrl;
            bdt_proton_rejection[ev_sample] = my_roc_vars.proton_rejection_bdt;
            xtrl_proton_rejection_power[ev_sample] = my_roc_vars.proton_efficiency_xtrl ? 1./my_roc_vars.proton_efficiency_xtrl : proton_reader->GetEntries();
            bdt_proton_rejection_power[ev_sample] = my_roc_vars.proton_efficiency_bdt ? 1./my_roc_vars.proton_efficiency_bdt : proton_reader->GetEntries();
            xtrl_signal_eff_over_proton_eff[ev_sample] = xtrl_proton_rejection_power[ev_sample]*my_roc_vars.electron_efficiency_xtrl;
            bdt_signal_eff_over_proton_eff[ev_sample] = bdt_proton_rejection_power[ev_sample]*my_roc_vars.electron_efficiency_bdt;
            
            roc_tree->Fill();
        }

        // Build TGraphs
        TGraph gr_xtrl_roc(samples, &xtrl_signal_eff[0], &xtrl_proton_rejection[0]);
        TGraph gr_bdt_roc(samples, &bdt_signal_eff[0], &bdt_proton_rejection[0]);
        TGraph gr_xtrl_proton_rejection_power(samples, &xtrl_signal_eff[0], &xtrl_proton_rejection_power[0]);
        TGraph gr_bdt_proton_rejection_power(samples, &bdt_signal_eff[0], &bdt_proton_rejection_power[0]);
        TGraph gr_xtrl_signal_eff_over_proton_eff(samples, &xtrl_signal_eff[0], &xtrl_signal_eff_over_proton_eff[0]);
        TGraph gr_bdt_signal_eff_over_proton_eff(samples, &bdt_signal_eff[0], &bdt_signal_eff_over_proton_eff[0]);

        gr_xtrl_roc.SetName("gr_xtrl_roc");
        gr_xtrl_roc.GetXaxis()->SetTitle("Signal efficiency");
        gr_xtrl_roc.GetYaxis()->SetTitle("Proton rejection");

        gr_bdt_roc.SetName("gr_bdt_roc");
        gr_bdt_roc.GetXaxis()->SetTitle("Signal efficiency");
        gr_bdt_roc.GetYaxis()->SetTitle("Proton rejection");

        gr_xtrl_proton_rejection_power.SetName("gr_xtrl_proton_rejection_power");
        gr_xtrl_proton_rejection_power.GetXaxis()->SetTitle("Signal efficiency");
        gr_xtrl_proton_rejection_power.GetYaxis()->SetTitle("1/Proton efficiency = Proton rejection power");

        gr_bdt_proton_rejection_power.SetName("gr_bdt_proton_rejection_power");
        gr_bdt_proton_rejection_power.GetXaxis()->SetTitle("Signal efficiency");
        gr_bdt_proton_rejection_power.GetYaxis()->SetTitle("1/Proton efficiency = Proton rejection power");

        gr_xtrl_signal_eff_over_proton_eff.SetName("gr_xtrl_signal_eff_over_proton_eff");
        gr_xtrl_signal_eff_over_proton_eff.GetXaxis()->SetTitle("Signal efficiency");
        gr_xtrl_signal_eff_over_proton_eff.GetYaxis()->SetTitle("Signal efficiency/Proton efficiency");

        gr_bdt_signal_eff_over_proton_eff.SetName("gr_bdt_signal_eff_over_proton_eff");
        gr_bdt_signal_eff_over_proton_eff.GetXaxis()->SetTitle("Signal efficiency");
        gr_bdt_signal_eff_over_proton_eff.GetYaxis()->SetTitle("Signal efficiency/Proton efficiency");

        gr_xtrl_roc.Write();
        gr_bdt_roc.Write();
        gr_xtrl_proton_rejection_power.Write();
        gr_bdt_proton_rejection_power.Write();
        gr_xtrl_signal_eff_over_proton_eff.Write();
        gr_bdt_signal_eff_over_proton_eff.Write();

        output->Write();
        output->Close();

    }

std::tuple<double, double> get_xtrl_efficiencies(
    std::shared_ptr<TTreeReader> electron_reader, 
    std::shared_ptr<TTreeReader> proton_reader, 
    const double & xtrl_cut_value, 
    const double &emin, 
    const double &emax) {

        TTreeReaderValue<double> electron_xtrl(*electron_reader, "xtrl");
        TTreeReaderValue<double> proton_xtrl(*proton_reader, "xtrl");
        TTreeReaderValue<double> corrected_energy_electron(*electron_reader, "energy_corr");
        TTreeReaderValue<double> corrected_energy_proton(*proton_reader, "energy_corr");

        double gev {0.001};

        unsigned int electron_events {0};
        unsigned int electron_survived_xtrl_cut {0};
        unsigned int electron_rejected_xtrl_cut {0};

        unsigned int proton_events {0};
        unsigned int proton_survived_xtrl_cut {0};
        unsigned int proton_rejected_xtrl_cut {0};

        while (electron_reader->Next()) {
            if (in_energy_range(*(corrected_energy_electron)*gev, emin, emax)) {
                electron_events += 1;
                if (*electron_xtrl<xtrl_cut_value)
                    electron_survived_xtrl_cut += 1;
                else
                    electron_rejected_xtrl_cut += 1;
            }
        }

        while (proton_reader->Next()) {
            if (in_energy_range(*(corrected_energy_proton)*gev, emin, emax)) {
                proton_events += 1;
                if (*proton_xtrl<xtrl_cut_value)
                    proton_survived_xtrl_cut += 1;
                else
                    proton_rejected_xtrl_cut += 1;
            }
        }
        
        electron_reader->Restart();
        proton_reader->Restart();
        
        double signal_selection {-999};
        double background_efficiency {-999};
        double proton_rejection {-999};

        if (electron_events)
            signal_selection = (double)electron_survived_xtrl_cut/electron_events;
        
        if (proton_events) {
            background_efficiency = (double)proton_survived_xtrl_cut/proton_events;
            proton_rejection = 1 - background_efficiency;
        }
       
        //return std::tuple<double, double> (signal_selection, proton_rejection);
        return std::tuple<double, double> (signal_selection, background_efficiency);
        //return std::tuple<double, double> ((double)electron_survived_xtrl_cut/electron_events, (double)proton_survived_xtrl_cut/proton_events);
    }

std::tuple<double, double> get_bdt_efficiencies(
    std::shared_ptr<TTreeReader> electron_reader, 
    std::shared_ptr<TTreeReader> proton_reader, 
    const double &bdt_cut_value, 
    const double &emin, 
    const double &emax,
    const char* bdt_weights_file,
    bdt_vars &my_bdt_vars,
    std::shared_ptr<TMVA::Reader> tmva_reader) {
        
        TTreeReaderValue<double> corrected_energy_electron(*electron_reader, "energy_corr");
        TTreeReaderValue<double> corrected_energy_proton(*proton_reader, "energy_corr");

        /*
        TTreeReaderValue<double> rmslayer_norm_1_electron(*electron_reader, "rmslayer_norm_1");
        TTreeReaderValue<double> rmslayer_norm_2_electron(*electron_reader, "rmslayer_norm_2");
        TTreeReaderValue<double> rmslayer_norm_3_electron(*electron_reader, "rmslayer_norm_3");
        TTreeReaderValue<double> rmslayer_norm_4_electron(*electron_reader, "rmslayer_norm_4");
        TTreeReaderValue<double> rmslayer_norm_5_electron(*electron_reader, "rmslayer_norm_5");
        TTreeReaderValue<double> rmslayer_norm_6_electron(*electron_reader, "rmslayer_norm_6");
        TTreeReaderValue<double> rmslayer_norm_7_electron(*electron_reader, "rmslayer_norm_7");
        TTreeReaderValue<double> rmslayer_norm_8_electron(*electron_reader, "rmslayer_norm_8");
        TTreeReaderValue<double> rmslayer_norm_9_electron(*electron_reader, "rmslayer_norm_9");
        TTreeReaderValue<double> rmslayer_norm_10_electron(*electron_reader, "rmslayer_norm_10");
        TTreeReaderValue<double> rmslayer_norm_11_electron(*electron_reader, "rmslayer_norm_11");
        TTreeReaderValue<double> rmslayer_norm_12_electron(*electron_reader, "rmslayer_norm_12");
        TTreeReaderValue<double> rmslayer_norm_13_electron(*electron_reader, "rmslayer_norm_13");
        TTreeReaderValue<double> rmslayer_norm_14_electron(*electron_reader, "rmslayer_norm_14");

        TTreeReaderValue<double> fraclayer_norm_1_electron(*electron_reader, "fraclayer_norm_1");
        TTreeReaderValue<double> fraclayer_norm_2_electron(*electron_reader, "fraclayer_norm_2");
        TTreeReaderValue<double> fraclayer_norm_3_electron(*electron_reader, "fraclayer_norm_3");
        TTreeReaderValue<double> fraclayer_norm_4_electron(*electron_reader, "fraclayer_norm_4");
        TTreeReaderValue<double> fraclayer_norm_5_electron(*electron_reader, "fraclayer_norm_5");
        TTreeReaderValue<double> fraclayer_norm_6_electron(*electron_reader, "fraclayer_norm_6");
        TTreeReaderValue<double> fraclayer_norm_7_electron(*electron_reader, "fraclayer_norm_7");
        TTreeReaderValue<double> fraclayer_norm_8_electron(*electron_reader, "fraclayer_norm_8");
        TTreeReaderValue<double> fraclayer_norm_9_electron(*electron_reader, "fraclayer_norm_9");
        TTreeReaderValue<double> fraclayer_norm_10_electron(*electron_reader, "fraclayer_norm_10");
        TTreeReaderValue<double> fraclayer_norm_11_electron(*electron_reader, "fraclayer_norm_11");
        TTreeReaderValue<double> fraclayer_norm_12_electron(*electron_reader, "fraclayer_norm_12");
        TTreeReaderValue<double> fraclayer_norm_13_electron(*electron_reader, "fraclayer_norm_13");
        TTreeReaderValue<double> fraclayer_norm_14_electron(*electron_reader, "fraclayer_norm_14");

        TTreeReaderValue<double> sumrms_norm_electron(*electron_reader, "sumrms_norm");
        TTreeReaderValue<double> fraclastlayer_norm_electron(*electron_reader, "fraclastlayer_norm");
        TTreeReaderValue<double> xtrl_norm_electron(*electron_reader, "xtrl_norm");
        TTreeReaderValue<double> xtrl_electron(*electron_reader, "xtrl");

        TTreeReaderValue<double> rmslayer_norm_1_proton(*proton_reader, "rmslayer_norm_1");
        TTreeReaderValue<double> rmslayer_norm_2_proton(*proton_reader, "rmslayer_norm_2");
        TTreeReaderValue<double> rmslayer_norm_3_proton(*proton_reader, "rmslayer_norm_3");
        TTreeReaderValue<double> rmslayer_norm_4_proton(*proton_reader, "rmslayer_norm_4");
        TTreeReaderValue<double> rmslayer_norm_5_proton(*proton_reader, "rmslayer_norm_5");
        TTreeReaderValue<double> rmslayer_norm_6_proton(*proton_reader, "rmslayer_norm_6");
        TTreeReaderValue<double> rmslayer_norm_7_proton(*proton_reader, "rmslayer_norm_7");
        TTreeReaderValue<double> rmslayer_norm_8_proton(*proton_reader, "rmslayer_norm_8");
        TTreeReaderValue<double> rmslayer_norm_9_proton(*proton_reader, "rmslayer_norm_9");
        TTreeReaderValue<double> rmslayer_norm_10_proton(*proton_reader, "rmslayer_norm_10");
        TTreeReaderValue<double> rmslayer_norm_11_proton(*proton_reader, "rmslayer_norm_11");
        TTreeReaderValue<double> rmslayer_norm_12_proton(*proton_reader, "rmslayer_norm_12");
        TTreeReaderValue<double> rmslayer_norm_13_proton(*proton_reader, "rmslayer_norm_13");
        TTreeReaderValue<double> rmslayer_norm_14_proton(*proton_reader, "rmslayer_norm_14");

        TTreeReaderValue<double> fraclayer_norm_1_proton(*proton_reader, "fraclayer_norm_1");
        TTreeReaderValue<double> fraclayer_norm_2_proton(*proton_reader, "fraclayer_norm_2");
        TTreeReaderValue<double> fraclayer_norm_3_proton(*proton_reader, "fraclayer_norm_3");
        TTreeReaderValue<double> fraclayer_norm_4_proton(*proton_reader, "fraclayer_norm_4");
        TTreeReaderValue<double> fraclayer_norm_5_proton(*proton_reader, "fraclayer_norm_5");
        TTreeReaderValue<double> fraclayer_norm_6_proton(*proton_reader, "fraclayer_norm_6");
        TTreeReaderValue<double> fraclayer_norm_7_proton(*proton_reader, "fraclayer_norm_7");
        TTreeReaderValue<double> fraclayer_norm_8_proton(*proton_reader, "fraclayer_norm_8");
        TTreeReaderValue<double> fraclayer_norm_9_proton(*proton_reader, "fraclayer_norm_9");
        TTreeReaderValue<double> fraclayer_norm_10_proton(*proton_reader, "fraclayer_norm_10");
        TTreeReaderValue<double> fraclayer_norm_11_proton(*proton_reader, "fraclayer_norm_11");
        TTreeReaderValue<double> fraclayer_norm_12_proton(*proton_reader, "fraclayer_norm_12");
        TTreeReaderValue<double> fraclayer_norm_13_proton(*proton_reader, "fraclayer_norm_13");
        TTreeReaderValue<double> fraclayer_norm_14_proton(*proton_reader, "fraclayer_norm_14");

        TTreeReaderValue<double> sumrms_norm_proton(*proton_reader, "sumrms_norm");
        TTreeReaderValue<double> fraclastlayer_norm_proton(*proton_reader, "fraclastlayer_norm");
        TTreeReaderValue<double> xtrl_norm_proton(*proton_reader, "xtrl_norm");
        TTreeReaderValue<double> xtrl_proton(*proton_reader, "xtrl");
        */

        TTreeReaderValue<double> rmslayer_norm_1_electron(*electron_reader, "rmslayer_1");
        TTreeReaderValue<double> rmslayer_norm_2_electron(*electron_reader, "rmslayer_2");
        TTreeReaderValue<double> rmslayer_norm_3_electron(*electron_reader, "rmslayer_3");
        TTreeReaderValue<double> rmslayer_norm_4_electron(*electron_reader, "rmslayer_4");
        TTreeReaderValue<double> rmslayer_norm_5_electron(*electron_reader, "rmslayer_5");
        TTreeReaderValue<double> rmslayer_norm_6_electron(*electron_reader, "rmslayer_6");
        TTreeReaderValue<double> rmslayer_norm_7_electron(*electron_reader, "rmslayer_7");
        TTreeReaderValue<double> rmslayer_norm_8_electron(*electron_reader, "rmslayer_8");
        TTreeReaderValue<double> rmslayer_norm_9_electron(*electron_reader, "rmslayer_9");
        TTreeReaderValue<double> rmslayer_norm_10_electron(*electron_reader, "rmslayer_10");
        TTreeReaderValue<double> rmslayer_norm_11_electron(*electron_reader, "rmslayer_11");
        TTreeReaderValue<double> rmslayer_norm_12_electron(*electron_reader, "rmslayer_12");
        TTreeReaderValue<double> rmslayer_norm_13_electron(*electron_reader, "rmslayer_13");
        TTreeReaderValue<double> rmslayer_norm_14_electron(*electron_reader, "rmslayer_14");

        TTreeReaderValue<double> fraclayer_norm_1_electron(*electron_reader, "fraclayer_1");
        TTreeReaderValue<double> fraclayer_norm_2_electron(*electron_reader, "fraclayer_2");
        TTreeReaderValue<double> fraclayer_norm_3_electron(*electron_reader, "fraclayer_3");
        TTreeReaderValue<double> fraclayer_norm_4_electron(*electron_reader, "fraclayer_4");
        TTreeReaderValue<double> fraclayer_norm_5_electron(*electron_reader, "fraclayer_5");
        TTreeReaderValue<double> fraclayer_norm_6_electron(*electron_reader, "fraclayer_6");
        TTreeReaderValue<double> fraclayer_norm_7_electron(*electron_reader, "fraclayer_7");
        TTreeReaderValue<double> fraclayer_norm_8_electron(*electron_reader, "fraclayer_8");
        TTreeReaderValue<double> fraclayer_norm_9_electron(*electron_reader, "fraclayer_9");
        TTreeReaderValue<double> fraclayer_norm_10_electron(*electron_reader, "fraclayer_10");
        TTreeReaderValue<double> fraclayer_norm_11_electron(*electron_reader, "fraclayer_11");
        TTreeReaderValue<double> fraclayer_norm_12_electron(*electron_reader, "fraclayer_12");
        TTreeReaderValue<double> fraclayer_norm_13_electron(*electron_reader, "fraclayer_13");
        TTreeReaderValue<double> fraclayer_norm_14_electron(*electron_reader, "fraclayer_14");

        TTreeReaderValue<double> sumrms_norm_electron(*electron_reader, "sumRms");
        TTreeReaderValue<double> fraclastlayer_norm_electron(*electron_reader, "fracLast");
        TTreeReaderValue<double> xtrl_norm_electron(*electron_reader, "xtrl");

        TTreeReaderValue<double> rmslayer_norm_1_proton(*proton_reader, "rmslayer_1");
        TTreeReaderValue<double> rmslayer_norm_2_proton(*proton_reader, "rmslayer_2");
        TTreeReaderValue<double> rmslayer_norm_3_proton(*proton_reader, "rmslayer_3");
        TTreeReaderValue<double> rmslayer_norm_4_proton(*proton_reader, "rmslayer_4");
        TTreeReaderValue<double> rmslayer_norm_5_proton(*proton_reader, "rmslayer_5");
        TTreeReaderValue<double> rmslayer_norm_6_proton(*proton_reader, "rmslayer_6");
        TTreeReaderValue<double> rmslayer_norm_7_proton(*proton_reader, "rmslayer_7");
        TTreeReaderValue<double> rmslayer_norm_8_proton(*proton_reader, "rmslayer_8");
        TTreeReaderValue<double> rmslayer_norm_9_proton(*proton_reader, "rmslayer_9");
        TTreeReaderValue<double> rmslayer_norm_10_proton(*proton_reader, "rmslayer_10");
        TTreeReaderValue<double> rmslayer_norm_11_proton(*proton_reader, "rmslayer_11");
        TTreeReaderValue<double> rmslayer_norm_12_proton(*proton_reader, "rmslayer_12");
        TTreeReaderValue<double> rmslayer_norm_13_proton(*proton_reader, "rmslayer_13");
        TTreeReaderValue<double> rmslayer_norm_14_proton(*proton_reader, "rmslayer_14");

        TTreeReaderValue<double> fraclayer_norm_1_proton(*proton_reader, "fraclayer_1");
        TTreeReaderValue<double> fraclayer_norm_2_proton(*proton_reader, "fraclayer_2");
        TTreeReaderValue<double> fraclayer_norm_3_proton(*proton_reader, "fraclayer_3");
        TTreeReaderValue<double> fraclayer_norm_4_proton(*proton_reader, "fraclayer_4");
        TTreeReaderValue<double> fraclayer_norm_5_proton(*proton_reader, "fraclayer_5");
        TTreeReaderValue<double> fraclayer_norm_6_proton(*proton_reader, "fraclayer_6");
        TTreeReaderValue<double> fraclayer_norm_7_proton(*proton_reader, "fraclayer_7");
        TTreeReaderValue<double> fraclayer_norm_8_proton(*proton_reader, "fraclayer_8");
        TTreeReaderValue<double> fraclayer_norm_9_proton(*proton_reader, "fraclayer_9");
        TTreeReaderValue<double> fraclayer_norm_10_proton(*proton_reader, "fraclayer_10");
        TTreeReaderValue<double> fraclayer_norm_11_proton(*proton_reader, "fraclayer_11");
        TTreeReaderValue<double> fraclayer_norm_12_proton(*proton_reader, "fraclayer_12");
        TTreeReaderValue<double> fraclayer_norm_13_proton(*proton_reader, "fraclayer_13");
        TTreeReaderValue<double> fraclayer_norm_14_proton(*proton_reader, "fraclayer_14");

        TTreeReaderValue<double> sumrms_norm_proton(*proton_reader, "sumRms");
        TTreeReaderValue<double> fraclastlayer_norm_proton(*proton_reader, "fracLast");
        TTreeReaderValue<double> xtrl_norm_proton(*proton_reader, "xtrl");

        double gev {0.001};

        unsigned int electron_events {0};
        unsigned int electron_survived_xtrl_cut {0};
        unsigned int electron_rejected_xtrl_cut {0};

        unsigned int proton_events {0};
        unsigned int proton_survived_xtrl_cut {0};
        unsigned int proton_rejected_xtrl_cut {0};

        while (electron_reader->Next()) {
            if (in_energy_range(*(corrected_energy_electron)*gev, emin, emax)) {
                electron_events += 1;
                
                /*
                my_bdt_vars.rmslayer_norm_1 = *rmslayer_norm_1_electron;
                my_bdt_vars.rmslayer_norm_2 = *rmslayer_norm_2_electron;
                my_bdt_vars.rmslayer_norm_3 = *rmslayer_norm_3_electron;
                my_bdt_vars.rmslayer_norm_4 = *rmslayer_norm_4_electron;
                my_bdt_vars.rmslayer_norm_5 = *rmslayer_norm_5_electron;
                my_bdt_vars.rmslayer_norm_6 = *rmslayer_norm_6_electron;
                my_bdt_vars.rmslayer_norm_7 = *rmslayer_norm_7_electron;
                my_bdt_vars.rmslayer_norm_8 = *rmslayer_norm_8_electron;
                my_bdt_vars.rmslayer_norm_9 = *rmslayer_norm_9_electron;
                my_bdt_vars.rmslayer_norm_10 = *rmslayer_norm_10_electron;
                my_bdt_vars.rmslayer_norm_11 = *rmslayer_norm_11_electron;
                my_bdt_vars.rmslayer_norm_12 = *rmslayer_norm_12_electron;
                my_bdt_vars.rmslayer_norm_13 = *rmslayer_norm_13_electron;
                my_bdt_vars.rmslayer_norm_14 = *rmslayer_norm_14_electron;

                my_bdt_vars.fraclayer_norm_1 = *fraclayer_norm_1_electron;
                my_bdt_vars.fraclayer_norm_2 = *fraclayer_norm_2_electron;
                my_bdt_vars.fraclayer_norm_3 = *fraclayer_norm_3_electron;
                my_bdt_vars.fraclayer_norm_4 = *fraclayer_norm_4_electron;
                my_bdt_vars.fraclayer_norm_5 = *fraclayer_norm_5_electron;
                my_bdt_vars.fraclayer_norm_6 = *fraclayer_norm_6_electron;
                my_bdt_vars.fraclayer_norm_7 = *fraclayer_norm_7_electron;
                my_bdt_vars.fraclayer_norm_8 = *fraclayer_norm_8_electron;
                my_bdt_vars.fraclayer_norm_9 = *fraclayer_norm_9_electron;
                my_bdt_vars.fraclayer_norm_10 = *fraclayer_norm_10_electron;
                my_bdt_vars.fraclayer_norm_11 = *fraclayer_norm_11_electron;
                my_bdt_vars.fraclayer_norm_12 = *fraclayer_norm_12_electron;
                my_bdt_vars.fraclayer_norm_13 = *fraclayer_norm_13_electron;
                my_bdt_vars.fraclayer_norm_14 = *fraclayer_norm_14_electron;
                
                my_bdt_vars.sumrms_norm = *sumrms_norm_electron;
                my_bdt_vars.fraclastlayer_norm = *fraclastlayer_norm_electron;
                my_bdt_vars.xtrl_norm = *xtrl_norm_electron;
                my_bdt_vars.xtrl = *xtrl_electron;
                */

                my_bdt_vars.rmslayer_1 = *rmslayer_norm_1_electron;
                my_bdt_vars.rmslayer_2 = *rmslayer_norm_2_electron;
                my_bdt_vars.rmslayer_3 = *rmslayer_norm_3_electron;
                my_bdt_vars.rmslayer_4 = *rmslayer_norm_4_electron;
                my_bdt_vars.rmslayer_5 = *rmslayer_norm_5_electron;
                my_bdt_vars.rmslayer_6 = *rmslayer_norm_6_electron;
                my_bdt_vars.rmslayer_7 = *rmslayer_norm_7_electron;
                my_bdt_vars.rmslayer_8 = *rmslayer_norm_8_electron;
                my_bdt_vars.rmslayer_9 = *rmslayer_norm_9_electron;
                my_bdt_vars.rmslayer_10 = *rmslayer_norm_10_electron;
                my_bdt_vars.rmslayer_11 = *rmslayer_norm_11_electron;
                my_bdt_vars.rmslayer_12 = *rmslayer_norm_12_electron;
                my_bdt_vars.rmslayer_13 = *rmslayer_norm_13_electron;
                my_bdt_vars.rmslayer_14 = *rmslayer_norm_14_electron;

                my_bdt_vars.fraclayer_1 = *fraclayer_norm_1_electron;
                my_bdt_vars.fraclayer_2 = *fraclayer_norm_2_electron;
                my_bdt_vars.fraclayer_3 = *fraclayer_norm_3_electron;
                my_bdt_vars.fraclayer_4 = *fraclayer_norm_4_electron;
                my_bdt_vars.fraclayer_5 = *fraclayer_norm_5_electron;
                my_bdt_vars.fraclayer_6 = *fraclayer_norm_6_electron;
                my_bdt_vars.fraclayer_7 = *fraclayer_norm_7_electron;
                my_bdt_vars.fraclayer_8 = *fraclayer_norm_8_electron;
                my_bdt_vars.fraclayer_9 = *fraclayer_norm_9_electron;
                my_bdt_vars.fraclayer_10 = *fraclayer_norm_10_electron;
                my_bdt_vars.fraclayer_11 = *fraclayer_norm_11_electron;
                my_bdt_vars.fraclayer_12 = *fraclayer_norm_12_electron;
                my_bdt_vars.fraclayer_13 = *fraclayer_norm_13_electron;
                my_bdt_vars.fraclayer_14 = *fraclayer_norm_14_electron;
                
                my_bdt_vars.sumRms = *sumrms_norm_electron;
                my_bdt_vars.fracLast = *fraclastlayer_norm_electron;
                my_bdt_vars.xtrl = *xtrl_norm_electron;

                double bdt_est = tmva_reader->EvaluateMVA( "BDT method" );

                if (bdt_est>bdt_cut_value)
                    electron_survived_xtrl_cut += 1;
                else
                    electron_rejected_xtrl_cut += 1;
            }
        }

        while (proton_reader->Next()) {
            if (in_energy_range(*(corrected_energy_proton)*gev, emin, emax)) {
                proton_events += 1;
                
                /*
                my_bdt_vars.rmslayer_norm_1 = *rmslayer_norm_1_proton;
                my_bdt_vars.rmslayer_norm_2 = *rmslayer_norm_2_proton;
                my_bdt_vars.rmslayer_norm_3 = *rmslayer_norm_3_proton;
                my_bdt_vars.rmslayer_norm_4 = *rmslayer_norm_4_proton;
                my_bdt_vars.rmslayer_norm_5 = *rmslayer_norm_5_proton;
                my_bdt_vars.rmslayer_norm_6 = *rmslayer_norm_6_proton;
                my_bdt_vars.rmslayer_norm_7 = *rmslayer_norm_7_proton;
                my_bdt_vars.rmslayer_norm_8 = *rmslayer_norm_8_proton;
                my_bdt_vars.rmslayer_norm_9 = *rmslayer_norm_9_proton;
                my_bdt_vars.rmslayer_norm_10 = *rmslayer_norm_10_proton;
                my_bdt_vars.rmslayer_norm_11 = *rmslayer_norm_11_proton;
                my_bdt_vars.rmslayer_norm_12 = *rmslayer_norm_12_proton;
                my_bdt_vars.rmslayer_norm_13 = *rmslayer_norm_13_proton;
                my_bdt_vars.rmslayer_norm_14 = *rmslayer_norm_14_proton;

                my_bdt_vars.fraclayer_norm_1 = *fraclayer_norm_1_proton;
                my_bdt_vars.fraclayer_norm_2 = *fraclayer_norm_2_proton;
                my_bdt_vars.fraclayer_norm_3 = *fraclayer_norm_3_proton;
                my_bdt_vars.fraclayer_norm_4 = *fraclayer_norm_4_proton;
                my_bdt_vars.fraclayer_norm_5 = *fraclayer_norm_5_proton;
                my_bdt_vars.fraclayer_norm_6 = *fraclayer_norm_6_proton;
                my_bdt_vars.fraclayer_norm_7 = *fraclayer_norm_7_proton;
                my_bdt_vars.fraclayer_norm_8 = *fraclayer_norm_8_proton;
                my_bdt_vars.fraclayer_norm_9 = *fraclayer_norm_9_proton;
                my_bdt_vars.fraclayer_norm_10 = *fraclayer_norm_10_proton;
                my_bdt_vars.fraclayer_norm_11 = *fraclayer_norm_11_proton;
                my_bdt_vars.fraclayer_norm_12 = *fraclayer_norm_12_proton;
                my_bdt_vars.fraclayer_norm_13 = *fraclayer_norm_13_proton;
                my_bdt_vars.fraclayer_norm_14 = *fraclayer_norm_14_proton;
                
                my_bdt_vars.sumrms_norm = *sumrms_norm_proton;
                my_bdt_vars.fraclastlayer_norm = *fraclastlayer_norm_proton;
                my_bdt_vars.xtrl_norm = *xtrl_norm_proton;
                my_bdt_vars.xtrl = *xtrl_proton;
                */

                my_bdt_vars.rmslayer_1 = *rmslayer_norm_1_proton;
                my_bdt_vars.rmslayer_2 = *rmslayer_norm_2_proton;
                my_bdt_vars.rmslayer_3 = *rmslayer_norm_3_proton;
                my_bdt_vars.rmslayer_4 = *rmslayer_norm_4_proton;
                my_bdt_vars.rmslayer_5 = *rmslayer_norm_5_proton;
                my_bdt_vars.rmslayer_6 = *rmslayer_norm_6_proton;
                my_bdt_vars.rmslayer_7 = *rmslayer_norm_7_proton;
                my_bdt_vars.rmslayer_8 = *rmslayer_norm_8_proton;
                my_bdt_vars.rmslayer_9 = *rmslayer_norm_9_proton;
                my_bdt_vars.rmslayer_10 = *rmslayer_norm_10_proton;
                my_bdt_vars.rmslayer_11 = *rmslayer_norm_11_proton;
                my_bdt_vars.rmslayer_12 = *rmslayer_norm_12_proton;
                my_bdt_vars.rmslayer_13 = *rmslayer_norm_13_proton;
                my_bdt_vars.rmslayer_14 = *rmslayer_norm_14_proton;

                my_bdt_vars.fraclayer_1 = *fraclayer_norm_1_proton;
                my_bdt_vars.fraclayer_2 = *fraclayer_norm_2_proton;
                my_bdt_vars.fraclayer_3 = *fraclayer_norm_3_proton;
                my_bdt_vars.fraclayer_4 = *fraclayer_norm_4_proton;
                my_bdt_vars.fraclayer_5 = *fraclayer_norm_5_proton;
                my_bdt_vars.fraclayer_6 = *fraclayer_norm_6_proton;
                my_bdt_vars.fraclayer_7 = *fraclayer_norm_7_proton;
                my_bdt_vars.fraclayer_8 = *fraclayer_norm_8_proton;
                my_bdt_vars.fraclayer_9 = *fraclayer_norm_9_proton;
                my_bdt_vars.fraclayer_10 = *fraclayer_norm_10_proton;
                my_bdt_vars.fraclayer_11 = *fraclayer_norm_11_proton;
                my_bdt_vars.fraclayer_12 = *fraclayer_norm_12_proton;
                my_bdt_vars.fraclayer_13 = *fraclayer_norm_13_proton;
                my_bdt_vars.fraclayer_14 = *fraclayer_norm_14_proton;
                
                my_bdt_vars.sumRms = *sumrms_norm_proton;
                my_bdt_vars.fracLast = *fraclastlayer_norm_proton;
                my_bdt_vars.xtrl = *xtrl_norm_proton;

                double bdt_est = tmva_reader->EvaluateMVA( "BDT method" );

                if (bdt_est>bdt_cut_value)
                    proton_survived_xtrl_cut += 1;
                else
                    proton_rejected_xtrl_cut += 1;
            }
        }

        electron_reader->Restart();
        proton_reader->Restart();
        
        double signal_selection {-999};
        double background_efficiency {-999};
        double proton_rejection {-999};

        if (electron_events)
            signal_selection = (double)electron_survived_xtrl_cut/electron_events;
        
        if (proton_events) {
            background_efficiency = (double)proton_survived_xtrl_cut/proton_events;
            proton_rejection = 1 - background_efficiency;
        }
       
        //return std::tuple<double, double> (signal_selection, proton_rejection);
        return std::tuple<double, double> (signal_selection, background_efficiency);
        //return std::tuple<double, double> ((double)electron_survived_xtrl_cut/electron_events, (double)proton_survived_xtrl_cut/proton_events);
    }