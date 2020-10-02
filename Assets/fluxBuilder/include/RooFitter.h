#ifndef ROOFITTER_H
#define ROOFITTER_H

#include <memory>
#include <vector>
#include <string>

#include "TH1D.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TCanvas.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooCBShape.h"
#include "RooBifurGauss.h"
#include "RooMinuit.h"
#include "RooArgList.h"
#include "RooMsgService.h"

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
	bool GetVerboseStatus();
	void Fit();
	void SaveResults(const std::string out_path, const bool release_flag=true);

private:
	// Initialization
	void init();
	void init_data(const std::vector<std::shared_ptr<TH1D>> &in_data);
	void init_template(const std::vector<std::vector<std::shared_ptr<TH1D>>> &in_templates);
	void init_guess(
		const std::vector<std::vector<double>> &guess,
		const std::vector<std::vector<bool>> &fix_guess);
	void set_verbose_status();
	void set_zero_roofit_verbosity();
	void reset_roofit_verbosity();
	// Fit facility
	void normalize_templates();
	void SetRooVars();
	void SetRooTemplates();
	void SetRooModel();
	void SetRooData();
	void PerformFit();
	void SetResult();
	void GetFitResult();
	// Plot facility
	void SuperimposeResults(TFile &outfile);

	// Data
	std::vector<unsigned int> data_events;
	std::vector<double> data_xmin;
	std::vector<double> data_xmax;
	std::vector<std::shared_ptr<TH1D>> data;
	// Templates
	std::vector<std::vector<std::shared_ptr<TH1D>>> norm_template;
	// Result
	std::vector<std::vector<double>> res;
	std::vector<std::vector<double>> res_err;
	std::vector<double> tot_err;
	std::vector<double> tot_res;
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
	std::vector<std::vector<std::shared_ptr<RooDataHist>>> roo_datahist_pdf;
	std::vector<std::vector<std::shared_ptr<RooHistPdf>>> roo_pdf;
	std::vector<std::shared_ptr<RooArgList>> roo_list_pdf;
	std::vector<std::shared_ptr<RooArgList>> roo_list_comp_var;
	std::vector<std::shared_ptr<RooAddPdf>> roo_model;
	std::vector<std::shared_ptr<RooDataHist>> roo_dataset;
	// Result histos
	std::vector<std::shared_ptr<TH1D>> roo_result;
	std::vector<std::vector<std::shared_ptr<TH1D>>> roo_result_comp;

	unsigned int _s_default = 2;
	double _d_default = -1;
	double _r_default = 0;
	bool _b_default = false;
};

#endif