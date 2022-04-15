#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "TMath.h"
#include "TAxis.h"
#include "TF1.h"
#include "TKey.h"
#include "TPad.h"
#include "TFile.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

const int DAMPE_bgo_nLayers {14};
const double chi2_limit {100};

inline const char* get_tree_name(const char* file) {
    TFile* input_file = TFile::Open(file, "READ");
    if (!input_file->IsOpen()) {
        std::cerr << "\n\nError reading input file [" << file << "]\n\n";
        exit(100);
    }
    std::string tree_name;
    for (TObject* keyAsObject : *input_file->GetListOfKeys()) {
        auto key = dynamic_cast<TKey*>(keyAsObject);
        if (!strcmp(key->GetClassName(), "TTree"))
            tree_name = static_cast<std::string>(key->GetName());
    }
    input_file->Close();
    return tree_name.c_str();
}

double shower_profile_fit_function(double *x, double *par)
{
    double xx = x[0];
    double func = par[0] * par[1] * (TMath::Power(par[1] * xx, par[2] - 1) * TMath::Exp(-par[1] * xx)) / TMath::Gamma(par[2]);
    return func;
}

void produceShowerProfilePlots(
    const char* input_file_name,
    const double emin_gev,
    const double emax_gev,
    int number_of_plots = -1,
    const bool apply_chi2_limit = false)
    {
        TFile *input_file = TFile::Open(input_file_name, "READ");
        if (input_file->IsZombie()) {
            std::cerr << "\n\nError opening input file [" << input_file << "]\n\n";
            exit(100);
        }

        TTreeReader myReader(get_tree_name(input_file_name), input_file);

        TTreeReaderValue<std::vector<double>> layer_energy (myReader, "eLayer");
        TTreeReaderValue<std::vector<double>> t_bgo (myReader, "t_bgo");
        
        TTreeReaderValue<double> energy (myReader, "energy");
        TTreeReaderValue<double> energy_corr (myReader, "energy_corr");

        TTreeReaderValue<bool> good_bgo_reco (myReader, "evtfilter_correct_bgo_reco");
        TTreeReaderValue<bool> evtfilter_maxRms_cut (myReader, "evtfilter_maxRms_cut");

        std::vector<TGraph> shower_profiles;
    
        // Shower variables
        const double bgoX0 				{11.2};     // BGO X0 in mm
        const double bgoEc 				{10.50};    // BGO critical energy for electrons in MeV
        const double b_shower_par 		{0.5}; 		// b parameter electromagnetic shower in BGO
        double a_shower_par;
        
        unsigned int evt_number {0};

        while (myReader.Next())
        {
            /*
            if (!*good_bgo_reco)
                continue;
            */
            if (!*evtfilter_maxRms_cut)
                continue;

            if (((*energy_corr*0.001) < emin_gev) || ((*energy_corr*0.001) > emax_gev))
                continue;

            // Define the shower profile TGraph
            TGraph shower_profile_gr(DAMPE_bgo_nLayers, &t_bgo->at(0), &layer_energy->at(0));
            auto gr_title = std::string("Shower Profile - event ") + std::to_string(evt_number) + std::string(" - corrected energy: ") + std::to_string((*energy_corr)*0.001) + std::string(" GeV; X_0; layer energy [MeV]");
            auto gr_name = std::string("shower_profile_evt_") + std::to_string(evt_number) + std::string("_corrected_energy_") + std::to_string((*energy_corr)*0.001) + std::string("_gev");
            shower_profile_gr.SetTitle(gr_title.c_str());
            shower_profile_gr.SetName(gr_name.c_str());
            shower_profile_gr.SetMarkerStyle(105);

            // Fit the shower profile
            a_shower_par = 1 + b_shower_par * (TMath::Log(*energy/bgoEc) - 0.5);

            TF1 fitfunc((std::string("fitfunc_") + std::to_string(evt_number)).c_str(), shower_profile_fit_function, t_bgo->at(0), t_bgo->at(DAMPE_bgo_nLayers-1), 3);
            fitfunc.SetParameter(0, *energy);
            fitfunc.SetParameter(1, b_shower_par);
            fitfunc.SetParameter(2, a_shower_par);
            fitfunc.SetNpx(10000);
            
            shower_profile_gr.Fit(&fitfunc, "qR");

            // Add the shower to the collection
            if (apply_chi2_limit)
            {
                if(fitfunc.GetChisquare()/fitfunc.GetNDF() < chi2_limit)
                    shower_profiles.push_back(shower_profile_gr);
            }
            else
                shower_profiles.push_back(shower_profile_gr);

            // Update event counter
            ++evt_number;
        }

        // Build the final canvas
        TCanvas print_canvas("print_canvas", "print_canvas");

        if (number_of_plots == -1)
            number_of_plots = static_cast<int>(shower_profiles.size());
        else
        {
            if (number_of_plots > static_cast<int>(shower_profiles.size()))
                number_of_plots = static_cast<int>(shower_profiles.size());
        }

        for (int elm = 0; elm < number_of_plots; ++elm)
        {
            
            // Set Range
            double ymin = TMath::MinElement(shower_profiles[elm].GetN(), shower_profiles[elm].GetY());
            double ymax = TMath::MaxElement(shower_profiles[elm].GetN(), shower_profiles[elm].GetY());
            if ((ymin - 10) < 0)
                ymin = 0;
            ymax += 20;

            shower_profiles[elm].GetYaxis()->SetRangeUser(ymin, ymax);
            shower_profiles[elm].Draw("AP");
            //gPad->SetLogy();
            gPad->SetGrid(1,1);

            if (!elm)
                print_canvas.Print("shower_profiles.pdf(","Title:Shower Profile Elm");
            else if (elm < number_of_plots -1)
                print_canvas.Print("shower_profiles.pdf","Title:Shower Profile Elm");
            else
                print_canvas.Print("shower_profiles.pdf)","Title:Shower Profile Elm");
        }

    }
