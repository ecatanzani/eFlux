#include <iostream>
#include <memory>
#include <vector>

#include "TAxis.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"

void buildEnergyResolutionPlot(
    const char* file_path, 
    const char* gr_name="energy_resolution_graph",
    const double energy_th = 20,
    const char* output_file_path = "energy_resolution_plot.root") {

        std::unique_ptr<TFile> file = std::make_unique<TFile>(file_path, "READ");
        if (file->IsZombie()) {
            std::cout << "Error opening input file [" << file_path << "]\n\n";
            exit(100);
        }

        TGraph* gr_resolution = static_cast<TGraph*>(file->Get(gr_name));

        std::vector<double> energy;
        std::vector<double> resolution;

        for (int pidx {0}; pidx<gr_resolution->GetN(); ++pidx) {
            auto tmp_energy = gr_resolution->GetPointX(pidx);
            auto tmp_resolution = gr_resolution->GetPointY(pidx)*100;

            if (tmp_energy>energy_th) {
                energy.push_back(tmp_energy);
                resolution.push_back(tmp_resolution);
            }
        }

        std::unique_ptr<TGraph> gr_edit_resolution = std::make_unique<TGraph>(static_cast<int>(energy.size()), energy.data(), resolution.data());
        gr_edit_resolution->GetXaxis()->SetTitle("Energy [GeV]");
        gr_edit_resolution->GetYaxis()->SetTitle("Energy Resolution [%]");
        gr_edit_resolution->SetTitle("Energy Resolution");
        gr_edit_resolution->SetMarkerStyle(20);

        std::unique_ptr<TCanvas> canvas_energy_resolution = std::make_unique<TCanvas>("canvas_energy_resolution", "canvas_energy_resolution", 500, 500);
        canvas_energy_resolution->SetTicks();

        gr_edit_resolution->Draw("ALP");

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        std::unique_ptr<TFile> ofile = std::make_unique<TFile>(output_file_path, "RECREATE");
        if (ofile->IsZombie()) {
            std::cout << "Error writing output file [" << output_file_path << "]\n\n";
            exit(100);
        }

        gr_edit_resolution->Write();
        canvas_energy_resolution->Write();

        ofile->Close();

    }