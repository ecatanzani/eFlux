#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iostream>

#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TMath.h"
#include "TGraphErrors.h"

const double generation_vertex_radius = 1.381976597885342;

std::vector<double> getEnergyBinning(const char* config_file);
std::vector<double> createLogBinning(const double eMin, const double eMax, const std::size_t n_bins);
double wtsydp(const double minene, const double maxene, const double index);

void buildAcceptance(
    const char* acc_config_file,
    const char* full_acceptance_histos, 
    const char* output_file = "acceptance.root")
{
    auto energy_binning = getEnergyBinning(acc_config_file);

    TFile* input_file = TFile::Open(full_acceptance_histos, "READ");
    if (input_file->IsZombie())
    {
        std::cerr << "\n\nError opening input file: [" << full_acceptance_histos << "]\n\n";
        exit(100);
    }

    TH1D* h_gen                     = static_cast<TH1D*>(input_file->Get("h_gen"));
    TH1D* h_geometric               = static_cast<TH1D*>(input_file->Get("h_geometric"));
    TH1D* h_geometric_trigger       = static_cast<TH1D*>(input_file->Get("h_geometric_trigger"));
    TH1D* h_bgo_fiducial            = static_cast<TH1D*>(input_file->Get("h_bgo_fiducial"));
    TH1D* h_all_cut                 = static_cast<TH1D*>(input_file->Get("h_all_cut"));

    h_gen                           ->SetDirectory(0);
    h_geometric                     ->SetDirectory(0);
    h_geometric_trigger             ->SetDirectory(0);
    h_bgo_fiducial                  ->SetDirectory(0);
    h_all_cut                       ->SetDirectory(0);

    input_file->Close();

    TH1D* h_geometric_factor        = static_cast<TH1D*>(h_geometric->Clone("h_geometric_factor"));
    TH1D* h_acc_geometric_trigger   = static_cast<TH1D*>(h_geometric_trigger->Clone("h_acc_geometric_trigger"));
    TH1D* h_acc_bgo_fiducial        = static_cast<TH1D*>(h_bgo_fiducial->Clone("h_acc_bgo_fiducial"));
    TH1D* h_acc_all_cut             = static_cast<TH1D*>(h_all_cut->Clone("h_acc_all_cut"));

    h_geometric_factor              ->Divide(h_gen);
    h_acc_geometric_trigger         ->Divide(h_gen);
    h_acc_bgo_fiducial              ->Divide(h_gen);
    h_acc_all_cut                   ->Divide(h_gen);

    double genSurface               = 4 * TMath::Pi() * pow(generation_vertex_radius, 2) / 2;
    double scaleFactor              = TMath::Pi() * genSurface;

    h_geometric_factor              ->Scale(scaleFactor);
    h_acc_geometric_trigger         ->Scale(scaleFactor);
    h_acc_bgo_fiducial              ->Scale(scaleFactor);
    h_acc_all_cut                   ->Scale(scaleFactor);

    h_geometric_factor              ->GetYaxis()->SetTitle("geometric factor [m^{2} sr]");
    h_acc_geometric_trigger         ->GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    h_acc_bgo_fiducial              ->GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    h_acc_all_cut                   ->GetYaxis()->SetTitle("acceptance [m^{2} sr]");

    std::vector<double> geometric_factor                (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_geometric_trigger           (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_bgo_fiducial                (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_all_cut                     (h_geometric_factor->GetNbinsX(), 0);

    std::vector<double> geometric_factor_err            (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_geometric_trigger_err       (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_bgo_fiducial_err            (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_all_cut_err                 (h_geometric_factor->GetNbinsX(), 0);

    std::vector<double> energy                          (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> energy_err                      (h_geometric_factor->GetNbinsX(), 0);

    for (unsigned int idx=0; idx<energy.size(); ++idx)
    {
        energy[idx] =                       wtsydp(energy_binning[idx], energy_binning[idx+1], -2.7);
        geometric_factor[idx] =             h_geometric_factor->GetBinContent(idx+1);
        geometric_factor_err[idx] =         h_geometric_factor->GetBinError(idx+1);
        acc_geometric_trigger[idx] =        h_acc_geometric_trigger->GetBinContent(idx+1);
        acc_geometric_trigger_err[idx] =    h_acc_geometric_trigger->GetBinError(idx+1);
        acc_bgo_fiducial[idx] =             h_acc_bgo_fiducial->GetBinContent(idx+1);
        acc_bgo_fiducial_err[idx] =         h_acc_bgo_fiducial->GetBinError(idx+1);
        acc_all_cut[idx] =                  h_acc_all_cut->GetBinContent(idx+1);
        acc_all_cut_err[idx] =              h_acc_all_cut->GetBinError(idx+1);
    }
    
    TGraphErrors gr_geometric_factor            ((int)geometric_factor.size(), &energy[0], &geometric_factor[0], &energy_err[0], &geometric_factor_err[0]);
    TGraphErrors gr_acc_geometric_trigger       ((int)geometric_factor.size(), &energy[0], &acc_geometric_trigger[0], &energy_err[0], &acc_geometric_trigger_err[0]);
    TGraphErrors gr_acc_bgo_fiducial            ((int)geometric_factor.size(), &energy[0], &acc_bgo_fiducial[0], &energy_err[0], &acc_bgo_fiducial_err[0]);
    TGraphErrors gr_acc_all_cut                 ((int)geometric_factor.size(), &energy[0], &acc_all_cut[0], &energy_err[0], &acc_all_cut_err[0]);

    gr_geometric_factor.SetName("gr_geometric_factor");
    gr_acc_geometric_trigger.SetName("gr_acc_geometric_trigger");
    gr_acc_bgo_fiducial.SetName("gr_acc_bgo_fiducial");
    gr_acc_all_cut.SetName("gr_acc_all_cut");

    gr_geometric_factor.SetTitle("Geometric Factor");
    gr_acc_geometric_trigger.SetTitle("Acceptance - geometric + trigger");
    gr_acc_bgo_fiducial.SetTitle("Acceptance - BGO fiducial");
    gr_acc_all_cut.SetTitle("Acceptance - all cuts");

    gr_geometric_factor.GetXaxis()->SetTitle("Energy [GeV]");
    gr_acc_geometric_trigger.GetXaxis()->SetTitle("Energy [GeV]");
    gr_acc_bgo_fiducial.GetXaxis()->SetTitle("Energy [GeV]");
    gr_acc_all_cut.GetXaxis()->SetTitle("Energy [GeV]");

    gr_geometric_factor.GetYaxis()->SetTitle("geometric factor [m^{2} sr]");
    gr_acc_geometric_trigger.GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    gr_acc_bgo_fiducial.GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    gr_acc_all_cut.GetYaxis()->SetTitle("acceptance [m^{2} sr]");

    TFile* outfile = TFile::Open(output_file, "RECREATE");
    if (outfile->IsZombie())
    {
        std::cerr << "\n\nError writing output acceptance file: [" << output_file << "]\n\n";
        exit(100);
    }

    gr_geometric_factor.Write();
    gr_acc_geometric_trigger.Write();
    gr_acc_bgo_fiducial.Write();
    gr_acc_all_cut.Write();

    outfile->Close();

    std::cout << "\n\nOutput file has been written: [" << output_file << "]\n\n";
}

std::vector<double> getEnergyBinning(const char* config_file)
{
    std::size_t n_bins = 0;
    double min_event_energy = -999;
    double max_event_energy = -999;
    std::vector<double> energy_binning;

    std::ifstream input_file(config_file);
    if (!input_file.is_open())
	{
		std::cerr << "\nInput acceptance config file not found [" << config_file << "]\n\n";
		exit(100);
	}
    std::string input_string(
		(std::istreambuf_iterator<char>(input_file)),
		(std::istreambuf_iterator<char>()));
	input_file.close();
    std::string tmp_str;
	std::istringstream input_stream(input_string);
	std::string::size_type sz;

    while (input_stream >> tmp_str)
	{
		if (!strcmp(tmp_str.c_str(), "n_energy_bins"))
			input_stream >> n_bins;
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

    return createLogBinning(min_event_energy, max_event_energy, n_bins);
}   

std::vector<double> createLogBinning(
    const double eMin,
	const double eMax,
	const std::size_t n_bins)
{
    std::vector<double> binning(n_bins+1, 0);
	double log_interval = (log10(eMax) - log10(eMin)) / n_bins;
	for (unsigned int bIdx = 0; bIdx <= n_bins; ++bIdx)
		binning[bIdx] = pow(10, log10(eMin) + bIdx * log_interval);

	return binning;
}

double wtsydp(const double minene, const double maxene, const double index)
{
    auto dene = maxene - minene;
    if (index != -1)
        return pow(fabs((pow(maxene, index + 1) - pow(minene, index + 1)) / ((index + 1) * dene)), 1. / index);
    else
        return dene / log(maxene / minene);
}