#ifndef EXTRACTOR_H
#define EXTRACTOR_H

#include <memory>
#include <vector>

#include "TH1D.h"
#include "TFile.h"
#include "TDirectory.h"

#include "main.h"

class extractor
{
public:
	extractor(const deps_paths paths, const bool verbose = false);
	~extractor(){};
	std::shared_ptr<TH1D> GetElectronAcceptance();
	std::shared_ptr<TH1D> GetProtonAcceptance();
	std::vector<std::shared_ptr<TH1D>> GetElectronTemplates();
	std::vector<std::shared_ptr<TH1D>> GetProtonTemplates();
	std::vector<std::shared_ptr<TH1D>> GetData();
	std::shared_ptr<TH1D> GetDataCounts();
	const int GetBins();
	void SaveResults();

private:
	void extract_electron_info(const std::string path, const bool verbose);
	void extract_proton_info(const std::string path, const bool verbose);
	void extract_data_info(const std::string path, const bool verbose);
	void extract_acceptance(std::shared_ptr<TH1D> &acc_ptr, TFile *input_file);
	void bin_consistency_check(const std::shared_ptr<TH1D> histo);
	void extract_xtrl_info(std::vector<std::shared_ptr<TH1D>> &xtrl_vec, TFile *input_file);

	std::shared_ptr<TH1D> electron_acceptance;
	std::shared_ptr<TH1D> proton_acceptance;
	std::vector<std::shared_ptr<TH1D>> electron_mc_xtrl;
	std::vector<std::shared_ptr<TH1D>> proton_mc_xtrl;
	std::vector<std::shared_ptr<TH1D>> data_xtrl;
	std::shared_ptr<TH1D> e_counts;
	int n_bins = 0;
};

#endif