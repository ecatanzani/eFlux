#include "RooFitter.h"

RooFitter::RooFitter(
	const std::vector<std::shared_ptr<TH1D>> &in_data,
	const std::vector<std::shared_ptr<TH1D>> &in_electron_templates,
	const std::vector<std::shared_ptr<TH1D>> &in_proton_templates,
	const unsigned int n_bins,
	const bool verbose,
	const bool clamping)
{
	if (!n_bins)
	{
		std::cout << "\nBinning could not be 0\n";
		exit(356);
	}
	bins = n_bins;
	verbosity = verbose;
	kClamping = clamping;
	init();
	init_data(in_data);
	init_template(in_electron_templates, in_proton_templates);
	set_verbose_status();
}

RooFitter::RooFitter(
	const std::vector<std::shared_ptr<TH1D>> &in_data,
	const std::vector<std::shared_ptr<TH1D>> &in_electron_templates,
	const std::vector<std::shared_ptr<TH1D>> &in_proton_templates,
	const std::vector<std::vector<double>> &guess,
	const std::vector<std::vector<bool>> &fix_guess,
	const unsigned int n_bins,
	const bool verbose,
	const bool clamping)
{
	if (!n_bins)
	{
		std::cout << "\nBinning could not be 0\n";
		exit(356);
	}
	bins = n_bins;
	verbosity = verbose;
	kClamping = clamping;
	init();
	init_data(in_data);
	init_template(in_electron_templates, in_proton_templates);
	set_verbose_status();
}

void RooFitter::init()
{
	data_events.resize(bins);
	data_xmin.resize(bins);
	data_xmax.resize(bins);
	norm_template.resize(bins);
	res.resize(bins);
	res_err.resize(bins);
	initial_guess.resize(bins);
	fix_to_initial_guess.resize(bins);
	roo_data_var.resize(bins);
	roo_comp_var.resize(bins);
	roo_datahist_pdf.resize(bins);
	roo_pdf.resize(bins);
	for (unsigned int idx=0; idx<bins; ++idx)
	{
		data_events[idx] = _d_default;
		data_xmin[idx] = _d_default;
		data_xmax[idx] = _d_default;
		norm_template[idx] = std::vector<std::shared_ptr<TH1D>>(_s_default, std::shared_ptr<TH1D>(nullptr));
		res[idx] = std::vector<double>(_s_default, _d_default);
		res_err[idx] = std::vector<double>(_s_default, _d_default);
		initial_guess[idx] = std::vector<double>(_s_default, _d_default);
		fix_to_initial_guess[idx] = std::vector<bool>(_s_default, _b_default);
		roo_comp_var[idx] = std::vector<std::shared_ptr<RooRealVar>> (_s_default, std::make_shared<RooRealVar>());
		roo_datahist_pdf[idx] = std::vector<std::shared_ptr<RooDataHist>> (_s_default, std::shared_ptr<RooDataHist>());
		roo_pdf[idx] = std::vector<std::shared_ptr<RooHistPdf>> (_s_default, std::shared_ptr<RooHistPdf>());
	}
}

void RooFitter::init_data(const std::vector<std::shared_ptr<TH1D>> &in_data)
{
	
	for (unsigned int idx=0; idx<bins; ++idx)
	{
		data_events[idx] = in_data[idx]->GetEntries();
		data_xmin[idx] = in_data[idx]->GetXaxis()->GetXmin();
		data_xmax[idx] = in_data[idx]->GetXaxis()->GetXmax();
	}
}

void RooFitter::init_template(
	const std::vector<std::shared_ptr<TH1D>> &in_electron_templates,
	const std::vector<std::shared_ptr<TH1D>> &in_proton_templates)
{
	std::vector<std::string> _temp_name = {"TemplateFit_e", "TemplateFit_p"};
	for (unsigned int idx=0; idx<bins; ++idx)
	{
		norm_template[idx][0] = std::shared_ptr<TH1D>(static_cast<TH1D *>(in_electron_templates[0]->Clone(_temp_name[0].c_str())));
		norm_template[idx][1] = std::shared_ptr<TH1D>(static_cast<TH1D *>(in_proton_templates[1]->Clone(_temp_name[1].c_str())));
	}
}

void RooFitter::init_guess(
	const std::vector<std::vector<double>> &guess,
	const std::vector<std::vector<bool>> &fix_guess)
{
	for (unsigned int idx=0; idx<bins; ++idx)
		for (unsigned int comp=0; comp<_s_default; ++comp)
		{
			initial_guess[idx][comp] = guess[idx][comp];
			fix_to_initial_guess[idx][comp] = fix_guess[idx][comp];
		}
}

void RooFitter::normalize_templates()
{
	for (auto&& _norm_v : norm_template)
		for (auto&& _elm : _norm_v)
			_elm->Scale(static_cast<double>(1) / _elm->Integral("width"));
}

bool RooFitter::GetVerboseStatus()
{
	return verbosity;
}

void RooFitter::SetRooVars()
{	
	double _initial_guess = 0;
	double _low_limit = 0;
	double _high_limit = 0;
	std::string roo_real_var_name;

	for (unsigned int idx=0; idx<bins; ++idx)
	{
		// Set X RooFit variable
		roo_data_var[idx] = std::make_shared<RooRealVar>("x", "x", data_xmin[idx], data_xmax[idx]);
		// Set component RooFit variables
		for (unsigned int comp=0; comp<_s_default; ++comp)
		{
			_initial_guess = data_events[idx]/static_cast<double>(_s_default);
			_low_limit = 0;
			_high_limit = data_events[idx];
			if (initial_guess[idx][comp]!=_d_default)
			{
				_initial_guess = initial_guess[idx][comp];
				if (fix_to_initial_guess[idx][comp])
					_low_limit = _high_limit = initial_guess[idx][comp];
				else
				{
					if (!kClamping)
					{
						_low_limit = -7*sqrt(data_events[idx]);
						_high_limit = data_events[idx] + 7*sqrt(data_events[idx]);
					}	
				}	
			}
			roo_real_var_name = "roo_comp_var_" + std::to_string(comp);
			roo_comp_var[idx][comp] = std::make_shared<RooRealVar>(
				roo_real_var_name.c_str(), 
				roo_real_var_name.c_str(), 
				_initial_guess,
				_low_limit,
				_high_limit);
		}
	}
}

void RooFitter::SetRooTemplates()
{
	std::string roo_pdfdh_name, roo_pdfh_name;

	for (unsigned int idx=0; idx<bins; ++idx)
	{	
		roo_pdfdh_name = "pdfdh_";
		roo_pdfh_name = "pdfh_";
		for(unsigned int comp=0; comp<_s_default; ++comp)
		{
			roo_pdfdh_name += std::to_string(comp);
			roo_pdfh_name += std::to_string(comp);
			roo_datahist_pdf[idx][comp] = std::make_shared<RooDataHist>(
				roo_pdfdh_name.c_str(),
				roo_pdfdh_name.c_str(),
				RooArgList(*roo_data_var[idx]),
				RooFit::Import(*norm_template[idx][comp]));
			roo_pdf[idx][comp] = std::make_shared<RooHistPdf>(
				roo_pdfh_name.c_str(),
				roo_pdfh_name.c_str(),
				RooArgSet(*roo_data_var[idx]),
				*roo_datahist_pdf[idx][comp]);
		}
	}
}

void RooFitter::Fit(const std::vector<std::shared_ptr<TH1D>> &in_data)
{
	// Normalize templates
	normalize_templates();
	// Set RooFit var
	SetRooVars();
	// Set RooFit templates
	SetRooTemplates();
}

void RooFitter::set_verbose_status()
{
	if (verbosity)
		set_zero_roofit_verbosity();
}

void RooFitter::set_zero_roofit_verbosity()
{
	RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization);
	RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
	RooMsgService::instance().getStream(1).removeTopic(RooFit::Fitting);
	RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval);
	RooMsgService::instance().getStream(1).removeTopic(RooFit::Caching);
	RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling);
	RooMsgService::instance().getStream(1).removeTopic(RooFit::InputArguments);
	RooMsgService::instance().getStream(1).removeTopic(RooFit::DataHandling);
	RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);
	RooMsgService::instance().getStream(0).removeTopic(RooFit::Minimization);
	RooMsgService::instance().getStream(0).removeTopic(RooFit::Plotting);
	RooMsgService::instance().getStream(0).removeTopic(RooFit::Fitting);
	RooMsgService::instance().getStream(0).removeTopic(RooFit::Eval);
	RooMsgService::instance().getStream(0).removeTopic(RooFit::Caching);
	RooMsgService::instance().getStream(0).removeTopic(RooFit::ObjectHandling);
	RooMsgService::instance().getStream(0).removeTopic(RooFit::InputArguments);
	RooMsgService::instance().getStream(0).removeTopic(RooFit::DataHandling);
	RooMsgService::instance().getStream(0).removeTopic(RooFit::Generation);
	RooMsgService::instance().getStream(0).removeTopic(RooFit::Integration);
	RooMsgService::instance().getStream(0).removeTopic(RooFit::LinkStateMgmt);
	RooMsgService::instance().getStream(0).removeTopic(RooFit::Optimization);
	RooMsgService::instance().getStream(0).removeTopic(RooFit::Tracing);
	RooMsgService::instance().getStream(0).removeTopic(RooFit::Contents);
	RooMsgService::instance().getStream(0).removeTopic(RooFit::NumIntegration);
	RooMsgService::instance().setSilentMode(true);
}

void RooFitter::reset_roofit_verbosity()
{
	RooMsgService::instance().reset();
}