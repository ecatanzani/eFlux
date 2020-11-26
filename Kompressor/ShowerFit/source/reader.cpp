#include "reader.h"
#include "DAMPE_geo_structure.h"

#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TMath.h"
#include "TGraph.h"
#include "TVector3.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include <vector>

void reader(
    const std::string wd,
    const std::string inputList,
    const std::string outputPath,
    const bool _VERBOSE,
    const bool mc)
{
    std::unique_ptr<parser> evt_parser = std::make_unique<parser>(inputList, mc, _VERBOSE);
    std::shared_ptr<config> _config = std::make_shared<config>(wd, mc);
    const double _entries = evt_parser->GetEvtTree()->GetEntries();
    if (_VERBOSE)
    {
        _config->PrintActiveFilters();
        std::cout << "Total number of events: " << _entries;
    }
    if (mc)
        mc_reader(evt_parser->GetEvtTree(), _config, _entries, outputPath, _VERBOSE);
}

double function(double *x, double *par)
{
    double xx = x[0];
    double func = par[0] * par[1] * (TMath::Power(par[1] * xx, par[2] - 1) * TMath::Exp(-par[1] * xx)) / TMath::Gamma(par[2]);
    return func;
}

void mc_reader(
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<config> _config,
    const double _entries,
    const std::string outputPath,
    const bool _VERBOSE)
{
    TTreeReader _reader(evtch.get());
    TTreeReaderValue<TVector3> bgo_trajectory(_reader, "BGOrec_trajectoryDirection2D");
    TTreeReaderValue<double> raw_energy(_reader, "energy");
    TTreeReaderValue<double> corr_energy(_reader, "energy_corr");
    TTreeReaderValue<std::vector<double>> layer_energy(_reader, "eLayer");
    TTreeReaderValue<bool> bgo_good_event(_reader, "evtfilter_good_event");

    const double bgoX0 = 11.2;       // BGO X0 in mm
    const double bgoEc = 10.50;      // BGO critical energy for electrons in MeV
    const double b_shower_par = 0.5; // b parameter electromagnetic shower in BGO
    const double _gev = 0.001;
    auto energy_binning = _config->GetEnergyBinning();
    std::vector<std::unique_ptr<TGraph>> gr_profile;

    // Histos
    std::unique_ptr<TH1D> h_a_shower_par = std::make_unique<TH1D>("h_a_shower_par", "a Shower Parameter", 100, 0, 100);
    std::unique_ptr<TH1D> h_b_shower_par = std::make_unique<TH1D>("h_b_shower_par", "b Shower Parameter", 100, 0, 1);
    std::unique_ptr<TH1D> h_corr_energy = std::make_unique<TH1D>("h_corr_energy", "DAMPE Corrected Energy", (int)energy_binning.size() - 1, &energy_binning[0]);
    std::unique_ptr<TH1D> h_energy_fit = std::make_unique<TH1D>("h_energy_fit", "Energy Fit", (int)energy_binning.size() - 1, &energy_binning[0]);
    std::unique_ptr<TH1D> h_energy_ratio = std::make_unique<TH1D>("h_energy_ratio", "Energy ratio; E_{corr,dampe}/E_{fit}; counts", 100, 0, 10);
    std::unique_ptr<TH1D> h_chi2_ndof = std::make_unique<TH1D>("h_chi2_ndof", "#chi^{2}/d.o.f.", 1000, 0, 1e+5);

    while (_reader.Next())
    {
        if (*bgo_good_event)
        {
            const double a_shower_par = 1 + b_shower_par * (TMath::Log(*raw_energy / bgoEc) - 0.5); // a parameter electromagnetic shower in BGO
            auto bgo_cosine = bgo_trajectory->CosTheta();
            auto raw_energy_gev = *raw_energy * _gev;
            auto corr_energy_gev = *corr_energy * _gev;
            std::vector<double> t_bgo(DAMPE_bgo_nLayers, -999);

            for (int idx = 0; idx < DAMPE_bgo_nLayers; ++idx)
                t_bgo[idx] = (BGO_bar_lateral * (idx + 1)) / (bgoX0 * bgo_cosine);

#if _DEBUG
            std::cout << "\nBGO costheta: " << bgo_cosine;
            std::cout << "\nShower profile scheme:\n";
            for (int idx = 0; idx < DAMPE_bgo_nLayers; ++idx)
                std::cout << "\nLayer " << idx << ": Energy: " << (*layer_energy)[idx] << " t: " << t_bgo[idx];
            std::cout << std::endl;
            std::cout << "\nDAMPE raw energy (GeV): " << raw_energy_gev;
            std::cout << "\nDAMPE corrected energy (GeV): " << corr_energy_gev << std::endl;
#endif
            gr_profile.push_back(std::make_unique<TGraph>(DAMPE_bgo_nLayers, &t_bgo.at(0), &layer_energy->at(0)));
            gr_profile.back()->SetName((std::string("shower_profile_evt_") + std::to_string(_reader.GetCurrentEntry())).c_str());
            gr_profile.back()->SetTitle((std::string("Shower Profile - Event ") + std::to_string(_reader.GetCurrentEntry())).c_str());
            TF1 fitfunc("fitfunc", function, t_bgo[0], t_bgo[DAMPE_bgo_nLayers - 1], 3);
            fitfunc.SetParameter(0, *raw_energy);
            //fitfunc.SetParLimits(0, 0, 1e+10);
            //fitfunc.FixParameter(1, b_shower_par);
            //fitfunc.FixParameter(2, a_shower_par);
            fitfunc.SetParameter(1, b_shower_par);
            fitfunc.SetParameter(2, a_shower_par);
            fitfunc.SetParNames("bgo_energy", "b", "a");
#if _DEBUG
            gr_profile.back()->Fit(&fitfunc, "IR");
#else
            gr_profile.back()->Fit(&fitfunc, "qIR");
#endif

            // Fill histos
            h_a_shower_par->Fill(fitfunc.GetParameter(2));
            h_b_shower_par->Fill(fitfunc.GetParameter(1));
            h_energy_fit->Fill(fitfunc.GetParameter(0) * _gev);
            h_energy_ratio->Fill((*corr_energy/fitfunc.GetParameter(0)));
            h_chi2_ndof->Fill(fitfunc.GetChisquare() / fitfunc.GetNDF());
            h_corr_energy->Fill(*corr_energy);
        }
    }

    TFile *_outfile = TFile::Open(outputPath.c_str(), "RECREATE");
    if (_outfile->IsZombie())
    {
        std::cerr << "\n\nError writing output ROOT file: [" << outputPath << "]\n\n";
        exit(100);
    }

    h_a_shower_par->Write();
    h_b_shower_par->Write();
    h_energy_fit->Write();
    h_energy_ratio->Write();
    h_chi2_ndof->Write();
    h_corr_energy->Write();

    for (auto &_elm : gr_profile)
        _elm->Write();

    _outfile->Close();
}