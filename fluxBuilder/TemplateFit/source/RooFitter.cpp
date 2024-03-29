#include "RooFitter.h"
#include "TLegend.h"

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
	std::vector<std::vector<std::shared_ptr<TH1D>>> in_templates {in_electron_templates, in_proton_templates};
	
	if (verbosity)
		std::cout << "***** RooFitter Data Analysis Class *****";
	if (verbosity)
		std::cout << std::endl << std::endl << "Initializing class...";
	init();
	if (verbosity)
		std::cout << std::endl << "Initializing data variables...";
	init_data(in_data);
	if (verbosity)
		std::cout << std::endl << "Initializing templates...";
	init_template(in_templates);
	if (verbosity)
		std::cout << std::endl << "Initializing RooFit verbosity...\n\n";
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
	std::vector<std::vector<std::shared_ptr<TH1D>>> in_templates {in_electron_templates, in_proton_templates};

	if (verbosity)
		std::cout << "***** RooFitter Data Analysis Class *****";
	if (verbosity)
		std::cout << std::endl << "Initializing class...";
	init();
	if (verbosity)
		std::cout << std::endl << "Initializing starting guess...";
	init_guess(guess, fix_guess);
	if (verbosity)
		std::cout << std::endl << "Initializing data variables...";
	init_data(in_data);
	if (verbosity)
		std::cout << std::endl << "Initializing templates...";
	init_template(in_templates);
	if (verbosity)
		std::cout << std::endl << "Initializing RooFit verbosity...\n\n";
	set_verbose_status();
}

void RooFitter::init()
{
	// Data
	data_events.resize(bins);
	data_xmin.resize(bins);
	data_xmax.resize(bins);
	data.resize(bins);
	// Templates
	norm_template.resize(bins);
	// Results
	res.resize(bins);
	res_err.resize(bins);
	tot_err.resize(bins);
	tot_res.resize(bins);
	// Options
	initial_guess.resize(bins);
	fix_to_initial_guess.resize(bins);
	// RooFit variables
	roo_data_var.resize(bins);
	roo_comp_var.resize(bins);
	roo_datahist_pdf.resize(bins);
	roo_pdf.resize(bins);
	roo_list_pdf.resize(bins);
	roo_list_comp_var.resize(bins);
	roo_model.resize(bins);
	roo_dataset.resize(bins);
	// Result histos
	roo_result.resize(bins);
	roo_result_comp.resize(bins);
	roo_proton_sample.resize(bins);

	for (unsigned int idx=0; idx<bins; ++idx)
	{
		// Data
		data_events[idx] = _d_default;
		data_xmin[idx] = _d_default;
		data_xmax[idx] = _d_default;
		data[idx] = std::make_shared<TH1D>();
		// Templates
		norm_template[idx] = std::vector<std::shared_ptr<TH1D>>(_s_default, std::make_shared<TH1D>());
		// Results
		res[idx] = std::vector<double>(_s_default, _d_default);
		res_err[idx] = std::vector<double>(_s_default, _d_default);
		tot_err[idx] = _r_default;
		tot_res[idx] = _r_default;
		// options
		initial_guess[idx] = std::vector<double>(_s_default, _d_default);
		fix_to_initial_guess[idx] = std::vector<bool>(_s_default, _b_default);
		// RooFit variables
		roo_data_var[idx] = std::make_shared<RooRealVar>();
		roo_comp_var[idx] = std::vector<std::shared_ptr<RooRealVar>> (_s_default, std::make_shared<RooRealVar>());
		roo_datahist_pdf[idx] = std::vector<std::shared_ptr<RooDataHist>> (_s_default, std::make_shared<RooDataHist>());
		roo_pdf[idx] = std::vector<std::shared_ptr<RooHistPdf>> (_s_default, std::make_shared<RooHistPdf>());
		roo_list_pdf[idx] = std::make_shared<RooArgList>();
		roo_list_comp_var[idx] = std::make_shared<RooArgList>();
		roo_model[idx] = std::make_shared<RooAddPdf>();
		roo_dataset[idx] = std::make_shared<RooDataHist>();
		// Result histos
		roo_result[idx] = std::make_shared<TH1D>();
		roo_result_comp[idx] = std::vector<std::shared_ptr<TH1D>> (_s_default, std::make_shared<TH1D>());
		roo_proton_sample[idx] = _d_default;
	}
}

void RooFitter::init_data(const std::vector<std::shared_ptr<TH1D>> &in_data)
{
	for (unsigned int idx=0; idx<bins; ++idx)
	{
		data_events[idx] = in_data[idx]->GetEntries();
		data_xmin[idx] = in_data[idx]->GetXaxis()->GetXmin();
		data_xmax[idx] = in_data[idx]->GetXaxis()->GetXmax();
		data[idx] = in_data[idx];
	}
}

void RooFitter::init_template(const std::vector<std::vector<std::shared_ptr<TH1D>>> &in_templates)
{
	std::vector<std::string> _temp_name;
	for (unsigned int idx=0; idx<bins; ++idx)
	{
		_temp_name = {"TemplateFit_e_", "TemplateFit_p_"};
		for (auto&& _elm : _temp_name)
			_elm += std::to_string(idx);

		for (unsigned int comp=0; comp<_s_default; ++comp)
			norm_template[idx][comp] = std::shared_ptr<TH1D>(static_cast<TH1D *>(in_templates[comp][idx]->Clone(_temp_name[comp].c_str())));
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
		roo_data_var[idx] = std::make_shared<RooRealVar>(
			"x", 
			"x", 
			data_xmin[idx], 
			data_xmax[idx]);
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
		for(unsigned int comp=0; comp<_s_default; ++comp)
		{
			roo_pdfdh_name = "pdfdh_" + std::to_string(comp);
			roo_pdfh_name = "pdfh_" + std::to_string(comp);
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

void RooFitter::SetRooModel()
{
	for (unsigned int idx=0; idx<bins; ++idx)
	{
		for (unsigned int comp=0; comp<_s_default; ++comp)
		{
			roo_list_pdf[idx]->add(*roo_pdf[idx][comp]);
			roo_list_comp_var[idx]->add(*roo_comp_var[idx][comp]);
		}
		roo_model[idx] = std::make_shared<RooAddPdf>(
			"model",
			"model",
			*roo_list_pdf[idx],
			*roo_list_comp_var[idx]);
	}
}

void RooFitter::SetRooData()
{
	for (unsigned int idx=0; idx<bins; ++idx)
	{
		roo_dataset[idx] = std::make_shared<RooDataHist>(
			"sumhist",
			"sumhsit",
			*roo_data_var[idx],
			RooFit::Import(*data[idx]));
#ifdef _DEBUG 
			std::cout << "\n\nRooDataHist [" << idx << "]:";
			roo_dataset[idx]->Print("v");
#endif	
	}
}

void RooFitter::PerformFit()
{
	for (unsigned int idx=0; idx<bins; ++idx)
		roo_model[idx]->fitTo(
			*roo_dataset[idx], 
			RooFit::Extended(true),
			RooFit::Save(true),
			RooFit::SumW2Error(true),
			RooFit::Verbose(verbosity),
			RooFit::PrintLevel(verbosity?-1:3),
			RooFit::Warnings(verbosity));
}

void RooFitter::SetResult()
{
	std::vector<std::string> res_name;
	for (unsigned int idx=0; idx<bins; ++idx)
	{
		res_name =  {"TemplateFitResult_", "TemplateFitResult_e_", "TemplateFitResult_p_"};
		for (auto&& _elm : res_name)
			_elm += std::to_string(idx);
		roo_result[idx] = std::shared_ptr<TH1D>(static_cast<TH1D*>(data[idx]->Clone(res_name[0].c_str())));
		roo_result[idx]->SetTitle(res_name[0].c_str());
		roo_result[idx]->Reset();
		for (unsigned int comp=0; comp<_s_default; ++comp)
		{
			roo_result_comp[idx][comp] = std::shared_ptr<TH1D>(static_cast<TH1D*>(data[idx]->Clone(res_name[comp+1].c_str())));
			roo_result_comp[idx][comp]->SetTitle(res_name[comp+1].c_str());
			roo_result_comp[idx][comp]->Reset();
			comp ? roo_result_comp[idx][comp]->SetLineColor(38) : roo_result_comp[idx][comp]->SetLineColor(46);
		}
	}
}

void RooFitter::GetFitResult()
{
	for (unsigned int idx=0; idx<bins; ++idx)
	{
		// Extract fit results and pars errors
		for(unsigned int comp=0; comp<_s_default; ++comp)
		{
			res[idx][comp] = roo_comp_var[idx][comp]->getVal();
			res_err[idx][comp] = roo_comp_var[idx][comp]->getError();
			tot_res[idx] += res[idx][comp];
			tot_err[idx] += pow(res_err[idx][comp], 2); //Fit using covariance matrix
			// build final histos
			roo_result[idx]->Add(
				norm_template[idx][comp].get(),
				res[idx][comp]/norm_template[idx][comp]->Integral());
			roo_result_comp[idx][comp]->Add(
				norm_template[idx][comp].get(),
				res[idx][comp]/norm_template[idx][comp]->Integral());
			if (comp)
				roo_proton_sample[idx] = roo_result_comp[idx][1]->Integral();
		}
		tot_err[idx] = sqrt(tot_err[idx]);
	}
}

void RooFitter::remove_zeros()
{
	for (auto&& _elm : norm_template)
		for (auto&& _ptr_i : _elm)
			for (int bIdx=1; bIdx<=_ptr_i->GetXaxis()->GetNbins(); ++bIdx)
				if (_ptr_i->GetBinContent(bIdx) < 1e-32)
					_ptr_i->SetBinContent(bIdx, 1e-32);
}

void RooFitter::Fit()
{
#ifdef _DEBUG 
	SaveResults(
		"debug_roofitter_init.root", 
		false,
		false);
#endif	
	remove_zeros();
	// Normalize templates
	normalize_templates();
#ifdef _DEBUG 
	SaveResults(
		"debug_roofitter_normtemplates.root", 
		false,
		false);
#endif		
	// Set RooFit var
	SetRooVars();
	// Set RooFit templates
	SetRooTemplates();
	// Create the model
	SetRooModel();
	// Create the data-set
	SetRooData();
#ifdef _DEBUG 
	SaveResults(
		"debug_roofitter_modelbeforefit.root", 
		false,
		true);
#endif
	// Fit
	PerformFit();
	// Get result
	SetResult();
	GetFitResult();
}

void RooFitter::SaveResults(
	const std::string out_path, 
	const bool release_flag,
	const bool model_ready)
{
	TFile outfile(out_path.c_str(), "RECREATE");
	if (outfile.IsZombie())
	{
		std::cerr << "\n\nError writing output ROOT file: [" << out_path << "]\n\n";
		exit(100);
	}
	
	// Create data folder
	auto data_folder = outfile.mkdir("data");
	data_folder->cd();
	for (auto& _elm : data)
		if (_elm)
			_elm->Write();

	// Create templates folder
	auto e_template_folder = outfile.mkdir("electron_templates");
	auto p_template_folder = outfile.mkdir("proton_templates");
	
	for (auto& _elm : norm_template)
		for (unsigned int comp=0; comp<_s_default; ++comp)
		{
			!comp ? e_template_folder->cd() : p_template_folder->cd();
			if (_elm[comp])
				_elm[comp]->Write();
		}

	if (release_flag)
	{
		// Create final histo dir
		auto result_folder = outfile.mkdir("results");
		result_folder->cd();
		for (auto& _elm : roo_result)
			if (_elm)
				_elm->Write();
			
		// Create final histo dir
		auto e_result_folder = outfile.mkdir("electron_results");
		auto p_result_folder = outfile.mkdir("proton_results");
		for (auto& _elm : roo_result_comp)
			for (unsigned int comp=0; comp<_s_default; ++comp)
			{
				!comp ? e_result_folder->cd() : p_result_folder->cd();
				if (_elm[comp])
					_elm[comp]->Write();
			}
	}
	
	if (model_ready)
	{
		// Create models dir
		auto model_dir = outfile.mkdir("models");
		std::vector<std::unique_ptr<RooPlot>> roo_plots (bins);
		std::string plot_name;
		model_dir->cd();
		for (unsigned int idx=0; idx<bins; ++idx)
		{
			plot_name = "TemplateFitModel_" + std::to_string(idx);
			roo_plots[idx] = std::unique_ptr<RooPlot>(static_cast<RooPlot*>(roo_data_var[idx]->frame()));
			roo_model[idx]->plotOn(roo_plots[idx].get());
			roo_plots[idx]->SetName(plot_name.c_str());
			roo_plots[idx]->SetTitle(plot_name.c_str());
			roo_plots[idx]->Write();
		}
	}

	outfile.cd();
	
	// Save component results into TTree
	double _e_val, _e_err, _p_val, _p_err;
	double _data, _data_min, _data_max;
	double _init_guess_e, _init_guess_p;
	bool _fix_guess_e, _fix_guess_p;
	
	TTree component_tree("result_tree", "result_tree");
	component_tree.Branch("e_component_val", &_e_val, "_e_val/D");
	component_tree.Branch("e_component_err", &_e_err, "_e_err/D");
	component_tree.Branch("p_component_val", &_p_val, "_p_val/D");
	component_tree.Branch("p_component_err", &_p_err, "_p_err/D");
	component_tree.Branch("data_xtrl_events", &_data, "_data/D");
	component_tree.Branch("data_min", &_data_min, "_data_min/D");
	component_tree.Branch("data_max", &_data_max, "_data_max/D");
	component_tree.Branch("init_guess_e", &_init_guess_e, "_init_guess_e/D");
	component_tree.Branch("init_guess_p", &_init_guess_p, "_init_guess_p/D");
	component_tree.Branch("fix_init_guess_e", &_fix_guess_e, "_fix_guess_e/O");
	component_tree.Branch("fix_init_guess_p", &_fix_guess_p, "_fix_guess_p/O");

	for (unsigned int idx=0; idx<bins; ++idx)
	{
		_e_val = res[idx][0];
		_p_val = res[idx][1];
		_e_err = res_err[idx][0];
		_p_err = res_err[idx][0];
		_data = data_events[idx];
		_data_min = data_xmin[idx];
		_data_max = data_xmax[idx];
		_init_guess_e = initial_guess[idx][0];
		_init_guess_p = initial_guess[idx][1];
		_fix_guess_e = fix_to_initial_guess[idx][0];
		_fix_guess_e = fix_to_initial_guess[idx][1];
		component_tree.Fill();
	}
	component_tree.Write();
		
	// Superimpose results
	if (release_flag)
		SuperimposeResults(outfile);

	outfile.Close();
}

void RooFitter::SuperimposeResults(TFile &outfile)
{
	std::vector<std::unique_ptr<TCanvas>> canvas (bins);
	std::vector<std::unique_ptr<TLegend>> _legend (bins);
	auto roo_result_cp = roo_result;
	auto roo_result_comp_cp = roo_result_comp;
	auto data_cp = data;
	std::string cvs_name; 
	int _line_width = 3;
	
	for (unsigned int idx=0; idx<bins; ++idx)
	{
		// Setting tmp canvas
		cvs_name = "TemplateFit_bin_" + std::to_string(idx);
		_legend[idx] = std::make_unique<TLegend>(.1,.7,.3,.9, cvs_name.c_str());
		_legend[idx]->AddEntry(data_cp[idx].get(), data_cp[idx]->GetTitle(), "l");
		_legend[idx]->AddEntry(roo_result_cp[idx].get(), roo_result_cp[idx]->GetTitle(), "l");
		canvas[idx] = std::make_unique<TCanvas>(cvs_name.c_str(), cvs_name.c_str());
		canvas[idx]->cd();
		// Setting roo_result_cp layout
		canvas[idx]->SetLogy();
		roo_result_cp[idx]->SetLineWidth(_line_width);
		roo_result_cp[idx]->SetLineColor(6);
		data_cp[idx]->SetLineWidth(_line_width);
		// Draw
		data_cp[idx]->Draw();
		data_cp[idx]->SetStats(0);	// Remove TStats
		roo_result_cp[idx]->Draw("histo same");
		for (auto& _elm : roo_result_comp_cp[idx])
		{
			// Setting component layout
			_elm->SetLineWidth(_line_width);
			_legend[idx]->AddEntry(_elm.get(), _elm->GetTitle(), "l");
			// Draw
			_elm->Draw("histo same");
		}
		_legend[idx]->Draw("same");
	}
	auto canvas_dir = outfile.mkdir("superimpose_cvs");
	canvas_dir->cd();
	for(auto& _elm : canvas)
		_elm->Write();
}

void RooFitter::set_verbose_status()
{
	if (!verbosity)
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