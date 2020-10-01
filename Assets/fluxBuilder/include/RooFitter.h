#ifndef ROOFITTER_H
#define ROOFITTER_H

#include <memory>
#include <vector>
#include <string>

#include "TH1D.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"

class RooFitter
{
public:
	RooFitter(
		const std::vector<std::shared_ptr<TH1D>> &in_data,
		const std::vector<std::shared_ptr<TH1D>> &in_electron_templates,
		const std::vector<std::shared_ptr<TH1D>> &in_proton_templates,
		const unsigned int n_bins,
		const bool verbose,
		const bool clamping);
	RooFitter(
		const std::vector<std::shared_ptr<TH1D>> &in_data,
		const std::vector<std::shared_ptr<TH1D>> &in_electron_templates,
		const std::vector<std::shared_ptr<TH1D>> &in_proton_templates,
		const std::vector<std::vector<double>> &guess,
		const std::vector<std::vector<bool>> &fix_guess,
		const unsigned int n_bins,
		const bool verbose,
		const bool clamping);
	~RooFitter(){};
	void norm_templates();
	bool GetVerboseStatus();
	void Fit(const std::vector<std::shared_ptr<TH1D>> &in_data);

private:
	// Initialization
	void init();
	void init_data(const std::vector<std::shared_ptr<TH1D>> &in_data);
	void init_template(
		const std::vector<std::shared_ptr<TH1D>> &in_electron_templates,
		const std::vector<std::shared_ptr<TH1D>> &in_proton_templates);
	void init_guess(
		const std::vector<std::vector<double>> &guess,
		const std::vector<std::vector<bool>> &fix_guess);
	void set_verbose_status();
	void set_zero_roofit_verbosity();
	void reset_roofit_verbosity();
	// Fit facility
	void SetRooVars();

	// Data
	std::vector<unsigned int> data_events;
	std::vector<double> data_xmin;
	std::vector<double> data_xmax;
	// Templates
	std::vector<std::vector<std::shared_ptr<TH1D>>> norm_template;
	// Result
	std::vector<std::vector<double>> res;
	std::vector<std::vector<double>> res_err;
	// Bins
	unsigned int bins = 0;
	// Options
	std::vector<std::vector<double>> initial_guess;
	std::vector<std::vector<bool>> fix_to_initial_guess;
	bool verbosity = false;
	bool kClamping = false;
	// RooFit variables
	std::vector<std::shared_ptr<RooRealVar>> roo_data_var;
	std::vector<std::vector<std::shared_ptr<RooRealVar>>> roo_comp_var;

	unsigned int _s_default = 2;
	double _d_default = -1;
	bool _b_default = false;
};

#endif