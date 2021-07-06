#include <vector>
#include <memory>
#include <string>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include <ROOT/RDataFrame.hxx>

#define dampe_bgo_layers 14

void getLambdaProfiles(
    const char* input_lambda_tree,
    const char* output_file = "lambda_profiles.root")
{
    TFile* infile = TFile::Open(input_lambda_tree, "READ");
    if (infile->IsZombie())
    {
        std::cerr << "\n\nError reading input file [" << input_lambda_tree << "]\n\n";
        exit(100);
    }

    std::shared_ptr<TTree> lambda_tree = std::shared_ptr<TTree>(static_cast<TTree*>(infile->Get("corrections_tree")));

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame _lambda_fr(*lambda_tree);

    std::vector<ROOT::RDF::RResultPtr<TGraph>> gr_rms (dampe_bgo_layers);
    std::vector<ROOT::RDF::RResultPtr<TGraph>> gr_elf (dampe_bgo_layers);

    for (int l_idx=0; l_idx<dampe_bgo_layers; ++l_idx)
    {
        auto get_lambda_layer = [l_idx](const std::vector<double> lambdas) -> double {return lambdas[l_idx]; };
        gr_rms[l_idx] = _lambda_fr.Define("lambda_layer", get_lambda_layer, {"best_rms_lambda"}).Graph("energy_bin", "lambda_layer");
        gr_rms[l_idx]->SetName((std::string("gr_rms_layer_") + std::to_string(l_idx)).c_str());
        gr_rms[l_idx]->GetXaxis()->SetTitle("energy bin");
        gr_rms[l_idx]->GetYaxis()->SetTitle("#lambda");
        gr_elf[l_idx] = _lambda_fr.Define("lambda_layer", get_lambda_layer, {"best_fraclayer_lambda"}).Graph("energy_bin", "lambda_layer");
        gr_elf[l_idx]->SetName((std::string("gr_elf_layer_") + std::to_string(l_idx)).c_str());
        gr_elf[l_idx]->GetXaxis()->SetTitle("energy bin");
        gr_elf[l_idx]->GetYaxis()->SetTitle("#lambda");
    }

    auto gr_sumrms = _lambda_fr.Graph("energy_bin", "best_sumrms_lambda");
    auto gr_ell = _lambda_fr.Graph("energy_bin", "best_fraclast_lambda");
    auto gr_xtrl = _lambda_fr.Graph("energy_bin", "best_xtrl_lambda");

    gr_sumrms->SetName("gr_sumrms");
    gr_sumrms->GetXaxis()->SetTitle("energy bin");
    gr_sumrms->GetYaxis()->SetTitle("#lambda");

    gr_ell->SetName("gr_ell");
    gr_ell->GetXaxis()->SetTitle("energy bin");
    gr_ell->GetYaxis()->SetTitle("#lambda");

    gr_xtrl->SetName("gr_xtrl");
    gr_xtrl->GetXaxis()->SetTitle("energy bin");
    gr_xtrl->GetYaxis()->SetTitle("#lambda");

    TFile* outfile = TFile::Open(output_file, "RECREATE");
    if (outfile->IsZombie())
    {
        std::cerr << "\n\nError reading input file [" << output_file << "]\n\n";
        exit(100);
    }


    outfile->mkdir("RMS");
    outfile->cd("RMS");
    for (int l_idx=0; l_idx<dampe_bgo_layers; ++l_idx)
        gr_rms[l_idx]->Write();

    outfile->mkdir("sumRMS");
    outfile->cd("sumRMS");
    gr_sumrms->Write();

    outfile->mkdir("ELF");
    outfile->cd("ELF");
    for (int l_idx=0; l_idx<dampe_bgo_layers; ++l_idx)
        gr_elf[l_idx]->Write();
    
    outfile->mkdir("ELL");
    outfile->cd("ELL");
    gr_ell->Write();

    outfile->mkdir("XTRL");
    outfile->cd("XTRL");
    gr_xtrl->Write();

    outfile->Close();
}