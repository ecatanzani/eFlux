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
	roo_var.resize(bins);
	for (auto it = begin(data_events); it != end(data_events); ++it)
	{
		static auto idx = std::distance(begin(data_events), it);
		data_events[idx] = -1;
		data_xmin[idx] = -1;
		data_xmax[idx] = -1;
		norm_template[idx] = std::vector<std::shared_ptr<TH1D>>(2, std::shared_ptr<TH1D>(nullptr));
		res[idx] = std::vector<double>(2, -1);
		res_err[idx] = std::vector<double>(2, -1);
		initial_guess[idx] = std::vector<double>(2, -1);
		fix_to_initial_guess[idx] = std::vector<bool>(2, false);
	}
}

void RooFitter::init_data(const std::vector<std::shared_ptr<TH1D>> &in_data)
{
	for (auto it = begin(in_data); it != end(in_data); ++it)
	{
		static auto idx = std::distance(begin(in_data), it);
		data_events[idx] = in_data[idx]->GetEntries();
		data_xmin[idx] = in_data[idx]->GetXaxis()->GetXmin();
		data_xmax[idx] = in_data[idx]->GetXaxis()->GetXmax();
	}
}

void RooFitter::init_template(
	const std::vector<std::shared_ptr<TH1D>> &in_electron_templates,
	const std::vector<std::shared_ptr<TH1D>> &in_proton_templates)
{
	for (auto it = begin(in_electron_templates); it != end(in_electron_templates); ++it)
	{
		static auto idx = std::distance(begin(in_electron_templates), it);
		static std::vector<std::string> _temp_name = {"TemplateFit_e", "TemplateFit_p"};
		for (int comp = 0; comp < 2; ++comp)
			norm_template[idx][comp] = std::shared_ptr<TH1D>(static_cast<TH1D *>(in_electron_templates[comp]->Clone(_temp_name[comp].c_str())));
	}
}

void RooFitter::init_guess(
	const std::vector<std::vector<double>> &guess,
	const std::vector<std::vector<bool>> &fix_guess)
{
	for (auto it = begin(initial_guess); it != end(initial_guess); ++it)
	{
		static auto idx = std::distance(begin(initial_guess), it);
		(*it)[0] = guess[idx][0];
		(*it)[1] = guess[idx][1];
		fix_to_initial_guess[idx][0] = fix_guess[idx][0];
		fix_to_initial_guess[idx][1] = fix_guess[idx][1];
	}
}

void RooFitter::norm_templates()
{
	for (auto &&elm : norm_template)
		for (auto it = begin(elm); it != end(elm); ++it)
			(*it)->Scale(1. / elm[0]->Integral("width"));
}

bool RooFitter::GetVerboseStatus()
{
	return verbosity;
}

void RooFitter::SetRooVar()
{
	auto it = begin(roo_var);
	for (auto &&var : roo_var)
	{
		static auto idx = std::distance(begin(roo_var), it);
		var = RooRealVar("x", "x", data_xmin[idx], data_xmax[idx]);
		++it;
	}
}

void RooFitter::Fit(const std::vector<std::shared_ptr<TH1D>> &in_data)
{
	// Normalize templates
	norm_templates();
	// Set RooFit var
	SetRooVar();
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