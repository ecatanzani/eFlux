#include "reader.h"
#include "DAMPE_geo_structure.h"

#include "TF1.h"
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

    const double bgoX0 = 11.2;                                                              // BGO X0 in mm
    const double bgoEc = 10.50;                                                             // BGO critical energy for electrons in MeV
    const double b_shower_par = 0.5;                                                        // b parameter electromagnetic shower in BGO
    const double a_shower_par = 1 + b_shower_par * (TMath::Log(*raw_energy / bgoEc) - 0.5); // a parameter electromagnetic shower in BGO

    double _gev = 0.001;
    auto bgo_cosine = bgo_trajectory->CosTheta();
    auto raw_energy_gev = *raw_energy * _gev;
    auto corr_energy_gev = *corr_energy * _gev;

    std::vector<double> t_bgo(DAMPE_bgo_nLayers, -999);
    for (int idx = 0; idx < DAMPE_bgo_nLayers; ++idx)
        t_bgo[idx] = (BGO_bar_lateral * (idx + 1)) / (bgoX0 * bgo_cosine);

    while (_reader.Next())
    {
        #if _DEBUG
            std::cout << "\nBGO costheta: " << bgo_cosine;
            std::cout << "\nShower profile scheme:\n";
            for (int idx = 0; idx < DAMPE_bgo_nLayers; ++idx)
                std::cout << "\nLayer " << idx << ": Energy: " << layer_energy->at(idx) << " t: " << t_bgo.at(idx);
            std::cout << std::endl;
            std::cout << "\nDAMPE raw energy (GeV): " << raw_energy_gev;
            std::cout << "\nDAMPE corrected energy (GeV): " << corr_energy_gev << std::endl;
        #endif
        TGraph gr_profile(DAMPE_bgo_nLayers, &t_bgo.at(0), &layer_energy->at(0));
        TF1 fitfunc("fitfunc", function, t_bgo[0], t_bgo[DAMPE_bgo_nLayers-1], 3);
        fitfunc.SetParameter(0, *raw_energy);
        //fitfunc.SetParLimits(0, 0, 1e+10);
        //fitfunc.FixParameter(1, b_shower_par);
        //fitfunc.FixParameter(2, a_shower_par);
        fitfunc.SetParameter(1, b_shower_par);
        fitfunc.SetParameter(2, a_shower_par);
        fitfunc.SetParNames("bgo_energy", "b", "a");
        #if _DEBUG    
            gr_profile.Fit(&fitfunc,"IR");
        #else
            gr_profile.Fit(&fitfunc,"qIR");
        #endif
    }
}