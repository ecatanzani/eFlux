#include <vector>
#include <memory>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include <ROOT/RDataFrame.hxx>

#define dampe_bgo_layers 14
struct energy_config
{
    std::size_t n_bins;
    double min_event_energy = -999;
    double max_event_energy = -999;
    std::vector<float> energy_binning;

    void createLogBinning()
    {
        if (min_event_energy!= -999 && max_event_energy != -999)
        {
            energy_binning.resize(n_bins + 1);
            double log_interval = (log10(max_event_energy) - log10(min_event_energy)) / n_bins;
            for (unsigned int bIdx = 0; bIdx <= n_bins; ++bIdx)
                energy_binning[bIdx] = pow(10, log10(min_event_energy) + bIdx * log_interval);
        }
    }
};

std::vector<float> getEnergyBins(const std::string config_file);
void fit(
    std::vector<ROOT::RDF::RResultPtr<TGraph>> &gr_rms,
    std::vector<ROOT::RDF::RResultPtr<TGraph>> &gr_rms_energylog,
    std::vector<ROOT::RDF::RResultPtr<TGraph>> &gr_elf,
    std::vector<ROOT::RDF::RResultPtr<TGraph>> &gr_elf_energylog,
    ROOT::RDF::RResultPtr<TGraph> gr_sumrms,
    ROOT::RDF::RResultPtr<TGraph> gr_sumrms_energylog,
    ROOT::RDF::RResultPtr<TGraph> gr_ell,
    ROOT::RDF::RResultPtr<TGraph> gr_ell_energylog,
    ROOT::RDF::RResultPtr<TGraph> gr_xtrl,
    ROOT::RDF::RResultPtr<TGraph> gr_xtrl_energylog);


void getLambdaProfiles(
    const char* input_lambda_tree,
    const char* config_file,
    const double spectrum_index = -2.7,
    const char* output_file = "lambda_profiles.root")
{
    auto energy_bins = getEnergyBins(config_file);

    TFile* infile = TFile::Open(input_lambda_tree, "READ");
    if (infile->IsZombie())
    {
        std::cerr << "\n\nError reading input file [" << input_lambda_tree << "]\n\n";
        exit(100);
    }

    std::shared_ptr<TTree> lambda_tree = std::shared_ptr<TTree>(static_cast<TTree*>(infile->Get("corrections_tree")));

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame _lambda_fr(*lambda_tree);

    std::vector<ROOT::RDF::RResultPtr<TGraph>> gr_rms (dampe_bgo_layers);
    std::vector<ROOT::RDF::RResultPtr<TGraph>> gr_elf (dampe_bgo_layers);
    std::vector<ROOT::RDF::RResultPtr<TGraph>> gr_rms_energylog (dampe_bgo_layers);
    std::vector<ROOT::RDF::RResultPtr<TGraph>> gr_elf_energylog (dampe_bgo_layers);
    
    auto getEnergyValue = [&energy_bins, &spectrum_index] (const int bin) -> double
    {
        auto dene = energy_bins[bin] - energy_bins[bin-1];
        if (spectrum_index != -1)
            return pow(fabs((pow(energy_bins[bin], spectrum_index + 1) - pow(energy_bins[bin-1], spectrum_index + 1)) / ((spectrum_index + 1) * dene)), 1. / spectrum_index);
        else
            return dene / log(energy_bins[bin] / energy_bins[bin-1]);
    };

    for (int l_idx=0; l_idx<dampe_bgo_layers; ++l_idx)
    {
        auto get_lambda_layer = [l_idx](const std::vector<double> lambdas) -> double {return lambdas[l_idx]; };
        gr_rms[l_idx] = _lambda_fr.Define("lambda_layer", get_lambda_layer, {"best_rms_lambda"}).Graph("energy_bin", "lambda_layer");
        gr_elf[l_idx] = _lambda_fr.Define("lambda_layer", get_lambda_layer, {"best_fraclayer_lambda"}).Graph("energy_bin", "lambda_layer");

        gr_rms_energylog[l_idx] = _lambda_fr.Define("lambda_layer", get_lambda_layer, {"best_rms_lambda"})
                                                .Define("energy_value", getEnergyValue, {"energy_bin"})
                                                .Graph("energy_value", "lambda_layer");

        gr_elf_energylog[l_idx] = _lambda_fr.Define("lambda_layer", get_lambda_layer, {"best_fraclayer_lambda"})
                                                .Define("energy_value", getEnergyValue, {"energy_bin"})
                                                .Graph("energy_value", "lambda_layer");

        gr_rms[l_idx]->SetName((std::string("gr_rms_layer_") + std::to_string(l_idx)).c_str());
        gr_rms[l_idx]->GetXaxis()->SetTitle("energy bin");
        gr_rms[l_idx]->GetYaxis()->SetTitle("#lambda");

        gr_rms_energylog[l_idx]->SetName((std::string("gr_rms_energylog_layer_") + std::to_string(l_idx)).c_str());
        gr_rms_energylog[l_idx]->GetXaxis()->SetTitle("Energy [GeV]");
        gr_rms_energylog[l_idx]->GetYaxis()->SetTitle("#lambda");
        
        gr_elf[l_idx]->SetName((std::string("gr_elf_layer_") + std::to_string(l_idx)).c_str());
        gr_elf[l_idx]->GetXaxis()->SetTitle("energy bin");
        gr_elf[l_idx]->GetYaxis()->SetTitle("#lambda");

        gr_elf_energylog[l_idx]->SetName((std::string("gr_elf_energylog_layer_") + std::to_string(l_idx)).c_str());
        gr_elf_energylog[l_idx]->GetXaxis()->SetTitle("Energy [GeV]");
        gr_elf_energylog[l_idx]->GetYaxis()->SetTitle("#lambda");
    }

    auto gr_sumrms = _lambda_fr.Graph("energy_bin", "best_sumrms_lambda");
    auto gr_ell = _lambda_fr.Graph("energy_bin", "best_fraclast_lambda");
    auto gr_xtrl = _lambda_fr.Graph("energy_bin", "best_xtrl_lambda");

    auto gr_sumrms_energylog = _lambda_fr.Define("energy_value", getEnergyValue, {"energy_bin"}).Graph("energy_value", "best_sumrms_lambda");
    auto gr_ell_energylog = _lambda_fr.Define("energy_value", getEnergyValue, {"energy_bin"}).Graph("energy_value", "best_fraclast_lambda");
    auto gr_xtrl_energylog = _lambda_fr.Define("energy_value", getEnergyValue, {"energy_bin"}).Graph("energy_value", "best_xtrl_lambda");

    gr_sumrms->SetName("gr_sumrms");
    gr_sumrms->GetXaxis()->SetTitle("energy bin");
    gr_sumrms->GetYaxis()->SetTitle("#lambda");

    gr_ell->SetName("gr_ell");
    gr_ell->GetXaxis()->SetTitle("energy bin");
    gr_ell->GetYaxis()->SetTitle("#lambda");

    gr_xtrl->SetName("gr_xtrl");
    gr_xtrl->GetXaxis()->SetTitle("energy bin");
    gr_xtrl->GetYaxis()->SetTitle("#lambda");

    gr_sumrms_energylog->SetName("gr_sumrms_energylog");
    gr_sumrms_energylog->GetXaxis()->SetTitle("Energy [GeV]");
    gr_sumrms_energylog->GetYaxis()->SetTitle("#lambda");

    gr_ell_energylog->SetName("gr_ell_energylog");
    gr_ell_energylog->GetXaxis()->SetTitle("Energy [GeV]");
    gr_ell_energylog->GetYaxis()->SetTitle("#lambda");

    gr_xtrl_energylog->SetName("gr_xtrl_energylog");
    gr_xtrl_energylog->GetXaxis()->SetTitle("Energy [GeV]");
    gr_xtrl_energylog->GetYaxis()->SetTitle("#lambda");

    TFile* outfile = TFile::Open(output_file, "RECREATE");
    if (outfile->IsZombie())
    {
        std::cerr << "\n\nError reading input file [" << output_file << "]\n\n";
        exit(100);
    }

    outfile->mkdir("RMS");
    outfile->cd("RMS");
    for (int l_idx=0; l_idx<dampe_bgo_layers; ++l_idx)
    {
        gr_rms[l_idx]->Write();
        gr_rms_energylog[l_idx]->Write();
    }

    outfile->mkdir("sumRMS");
    outfile->cd("sumRMS");
    gr_sumrms->Write();
    gr_sumrms_energylog->Write();

    outfile->mkdir("ELF");
    outfile->cd("ELF");
    for (int l_idx=0; l_idx<dampe_bgo_layers; ++l_idx)
    {
        gr_elf[l_idx]->Write();
        gr_elf_energylog[l_idx]->Write();
    }

    outfile->mkdir("ELL");
    outfile->cd("ELL");
    gr_ell->Write();
    gr_ell_energylog->Write();

    outfile->mkdir("XTRL");
    outfile->cd("XTRL");
    gr_xtrl->Write();
    gr_xtrl_energylog->Write();

    outfile->Close();
}

std::vector<float> getEnergyBins(const std::string config_file)
{
    std::ifstream input_file(config_file.c_str());
    if (!input_file.is_open())
	{
		std::cerr << "\nInput config file not found [" << config_file << "]\n\n";
		exit(100);
	}
    std::string input_string((std::istreambuf_iterator<char>(input_file)), (std::istreambuf_iterator<char>()));
	input_file.close();
    std::string tmp_str;
	std::istringstream input_stream(input_string);
	std::string::size_type sz;
    auto energy = energy_config();
    while (input_stream >> tmp_str)
	{
		if (!strcmp(tmp_str.c_str(), "n_energy_bins"))
			input_stream >> energy.n_bins;
		if (!strcmp(tmp_str.c_str(), "min_event_energy"))
		{
			input_stream >> tmp_str;
			energy.min_event_energy = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "max_event_energy"))
		{
			input_stream >> tmp_str;
			energy.max_event_energy = stod(tmp_str, &sz);
		}
    }

    energy.createLogBinning();
    return energy.energy_binning;
}

void fit(
    std::vector<ROOT::RDF::RResultPtr<TGraph>> &gr_rms,
    std::vector<ROOT::RDF::RResultPtr<TGraph>> &gr_elf,
    ROOT::RDF::RResultPtr<TGraph> gr_sumrms,
    ROOT::RDF::RResultPtr<TGraph> gr_ell,
    ROOT::RDF::RResultPtr<TGraph> gr_xtrl)
{
    gr_sumrms->Fit("pol1");
    gr_sumrms->Fit("pol1");
    gr_sumrms->Fit("pol1");

}