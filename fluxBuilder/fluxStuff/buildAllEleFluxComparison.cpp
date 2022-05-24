#include <iostream>
#include <vector>
#include <memory>

#include "TFile.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

void buildAllEleFluxComparison(
    const char* input_flux_path,
    const char* input_reference_path,
    const char* output_path,
    double input_flux_systematics = 0.025) {

        auto apply_sistematic_correction = [input_flux_systematics](std::shared_ptr<TGraphErrors> graph) -> std::shared_ptr<TGraphErrors> {
            for (int pidx=0; pidx<graph->GetN(); ++pidx)
                graph->SetPointError(pidx, 0, graph->GetErrorY(pidx) + graph->GetPointY(pidx)*input_flux_systematics);
            return graph;
        };

        auto build_flux_ratio = [](std::shared_ptr<TGraphErrors> flux_input, std::shared_ptr<TGraphAsymmErrors> flux_reference) -> std::shared_ptr<TGraphErrors> {
            std::vector<double> diff;
            std::vector<double> energy;
            std::vector<double> err_x, err_y;

            for (int pidx=0; pidx<flux_reference->GetN(); ++pidx) {
                double reference_energy_val = flux_reference->GetX()[pidx];
                double diff_val = (flux_input->Eval(reference_energy_val) - flux_reference->GetY()[pidx]) / flux_reference->GetErrorY(pidx);
                diff.push_back(diff_val);
                energy.push_back(reference_energy_val);
                err_x.push_back(0);
                err_y.push_back(1.0);
            }

            std::shared_ptr<TGraphErrors> graph = std::make_shared<TGraphErrors>(energy.size(), &energy[0], &diff[0], &err_x[0], &err_y[0]);

            graph->SetName("gr_flux_ratio");
            graph->SetTitle("Flux ratio");
            graph->GetXaxis()->SetTitle("Energy [GeV]");
            graph->GetYaxis()->SetTitle("(flux - DAMPE Nature)/err [#sigma]");

            return graph;
        };

        // Open the input files and extract plots
        TFile input_flux(input_flux_path, "READ");
        if (input_flux.IsZombie()) {
            std::cerr << "\n\nError opening input flux [" << input_flux_path << "]\n\n";
            exit(100);
        }

        std::shared_ptr<TGraphErrors> gr_input_flux = apply_sistematic_correction(std::shared_ptr<TGraphErrors>(static_cast<TGraphErrors*>(input_flux.Get("gr_flux_E3"))));
        gr_input_flux->SetLineWidth(2);
        gr_input_flux->SetLineColor(kRed);
        gr_input_flux->SetMarkerColor(kRed);
        gr_input_flux->SetMarkerStyle(20);

        gr_input_flux->SetTitle("DAMPE (2016/01 - 2021/11)");

        input_flux.Close();

        TFile input_reference(input_reference_path, "READ");
        if (input_reference.IsZombie()) {
            std::cerr << "\n\nError opening input reference [" << input_reference_path << "]\n\n";
            exit(100);
        }

        std::shared_ptr<TGraphAsymmErrors> gr_input_reference = std::shared_ptr<TGraphAsymmErrors>(static_cast<TGraphAsymmErrors*>(input_reference.Get("graph1")));
        gr_input_reference->SetLineColor(kBlue);
        gr_input_reference->SetMarkerColor(kBlue);
        gr_input_reference->SetMarkerStyle(20);

        gr_input_reference->SetTitle("DAMPE Nature (2017)");

        input_reference.Close();

        // Build the flux ratio
        auto gr_ratio_sigma = build_flux_ratio(gr_input_flux, gr_input_reference);

        // Create output file
        TFile output_file(output_path, "RECREATE");
        if (output_file.IsZombie()) {
            std::cerr << "\n\nError writing output ROOT file [" << output_path << "]\n\n";
            exit(100);
        }

        // Create output canvas
        TCanvas output_canvas("output_canvas", "output_canvas");
        output_canvas.SetTicks();

        gr_input_flux->GetXaxis()->SetLimits(20, 1e+4);
        gr_input_flux->SetMinimum(0);
        gr_input_flux->SetMaximum(300);
        
        gr_input_reference->Draw("AP");
        gr_input_flux->Draw("P");

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        auto legend = output_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("lp");
        }

        output_canvas.Write();

        TCanvas output_canvas_ratio("output_canvas_ratio", "output_canvas_ratio");
        output_canvas_ratio.SetTicks();

        gr_ratio_sigma->SetMinimum(-5);
        gr_ratio_sigma->SetMaximum(5);

        gr_ratio_sigma->SetLineWidth(2);
        gr_ratio_sigma->SetLineColor(kBlack);
        gr_ratio_sigma->SetMarkerColor(kBlack);
        gr_ratio_sigma->SetMarkerStyle(20);

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        gr_ratio_sigma->Draw("AP");

        output_canvas_ratio.Write();

        output_file.Close();

    }
