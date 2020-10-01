#include "extractor.h"

extractor::extractor(const deps_paths paths, const bool verbose)
{
	if (verbose)
		std::cout << "\n**** Extracting deps...\n\n";
	extract_electron_info(paths.electron_path, verbose);
	extract_proton_info(paths.proton_path, verbose);
	extract_data_info(paths.data_path, verbose);
	if (verbose)
		std::cout << "\n******************\n\n";
}

void extractor::extract_acceptance(std::shared_ptr<TH1D> &acc_ptr, TFile *input_file)
{
	acc_ptr = std::shared_ptr<TH1D>(static_cast<TH1D *>(input_file->Get("Acceptance_histos/h_acceptance_all_cut")));
	acc_ptr->SetDirectory(0);
}

void extractor::bin_consistency_check(const std::shared_ptr<TH1D> histo)
{
	if (histo->GetNbinsX() != n_bins)
	{
		std::cerr << "\nBin error mismatch: " << histo->GetName() << " != class bin content (" << n_bins << "\n\n";
		exit(156);
	}
}

void extractor::extract_xtrl_info(std::vector<std::shared_ptr<TH1D>> &xtrl_vec, TFile *input_file)
{
	// Set XTRL vector dimension
	if (n_bins)
		xtrl_vec.resize(n_bins);
	else
	{
		std::cerr << "\nBin error: binning set to 0 !\n\n";
		exit(256);
	}
	// Reading XTRL distributions
	for (int bIdx = 0; bIdx < n_bins; ++bIdx)
	{
		std::string xtrl_bin_name = "xtrl/h_xtrl_bin_" + std::to_string(bIdx);
		xtrl_vec[bIdx] = std::shared_ptr<TH1D>(static_cast<TH1D *>(input_file->Get(xtrl_bin_name.c_str())));
		xtrl_vec[bIdx]->SetDirectory(0);
	}
}

void extractor::extract_electron_info(const std::string path, const bool verbose)
{
	TFile *input_file = TFile::Open(path.c_str(), "READ");
	if (input_file->IsZombie())
	{
		std::cerr << "\nError reading input electron MC file [" << path << "]\n\n";
		exit(100);
	}
	if (verbose)
		std::cout << "Reading input electron MC file [" << path << "]" << std::endl;

	// Read acceptance
	extract_acceptance(electron_acceptance, input_file);
	// Set binning
	n_bins = electron_acceptance->GetNbinsX();
	// Extract XTRL info
	extract_xtrl_info(electron_mc_xtrl, input_file);
	// Close input file
	input_file->Close();
}

void extractor::extract_proton_info(const std::string path, const bool verbose)
{
	TFile *input_file = TFile::Open(path.c_str(), "READ");
	if (input_file->IsZombie())
	{
		std::cerr << "\nError reading input proton MC file [" << path << "]\n\n";
		exit(100);
	}
	if (verbose)
		std::cout << "Reading input proton MC file [" << path << "]" << std::endl;

	// Read acceptance
	extract_acceptance(proton_acceptance, input_file);
	// Bin consistency check
	bin_consistency_check(proton_acceptance);
	// Extract XTRL info
	extract_xtrl_info(proton_mc_xtrl, input_file);
	// Close input file
	input_file->Close();
}

void extractor::extract_data_info(const std::string path, const bool verbose)
{
	TFile *input_file = TFile::Open(path.c_str(), "READ");
	if (input_file->IsZombie())
	{
		std::cerr << "\nError reading input DATA file [" << path << "]\n\n";
		exit(100);
	}
	if (verbose)
		std::cout << "Reading input DATA file [" << path << "]" << std::endl;
	// Read data electron counts
	e_counts = std::shared_ptr<TH1D>(static_cast<TH1D *>(input_file->Get("h_all_cut_ce")));
	e_counts->SetDirectory(0);
	// Extract XTRL info
	extract_xtrl_info(data_xtrl, input_file);
	// Close input file
	input_file->Close();
}

std::shared_ptr<TH1D> extractor::GetElectronAcceptance()
{
	if (electron_acceptance)
		return electron_acceptance;
	else
		return std::shared_ptr<TH1D>(nullptr);
}

std::shared_ptr<TH1D> extractor::GetProtonAcceptance()
{
	if (proton_acceptance)
		return proton_acceptance;
	else
		return std::shared_ptr<TH1D>(nullptr);
}

std::vector<std::shared_ptr<TH1D>> extractor::GetElectronTemplates()
{
	if (!electron_mc_xtrl.empty())
		return electron_mc_xtrl;
	else
		return std::vector<std::shared_ptr<TH1D>>();
}

std::vector<std::shared_ptr<TH1D>> extractor::GetProtonTemplates()
{
	if (!proton_mc_xtrl.empty())
		return proton_mc_xtrl;
	else
		return std::vector<std::shared_ptr<TH1D>>();
}

std::vector<std::shared_ptr<TH1D>> extractor::GetData()
{
	if (!data_xtrl.empty())
		return data_xtrl;
	else
		return std::vector<std::shared_ptr<TH1D>>();
}

std::shared_ptr<TH1D> extractor::GetDataCounts()
{
	if (e_counts)
		return e_counts;
	else
		return std::shared_ptr<TH1D>(nullptr);
}

const int extractor::GetBins()
{
	return n_bins;
}
