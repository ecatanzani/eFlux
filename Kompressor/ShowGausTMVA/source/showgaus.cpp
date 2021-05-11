#include "main.h"
#include "gaus.h"
#include "parser.h"
#include "histos.h"

void showgaus(in_args input_args)
{
    auto evtch = getchain(
        input_args.input_list, 
        input_args.mc_flag, 
        input_args.verbose);
    ROOT::EnableImplicitMT(input_args.threads);
    ROOT::RDataFrame _data_fr(*evtch);
    if (input_args.verbose) std::cout << "\n\nTotal number of events: " << *(_data_fr.Count());
    auto hrmslayer = getrmslayerhistos();
    auto hefraclayer = getenergyfractionlayerhistos();

    for (int bidx=0; bidx<nenergybin; ++bidx)
    {
        auto bin_filter = [&bidx](const int energy_bin) -> bool { return energy_bin == bidx; };
        _data_fr.Filter(bin_filter, {"energy_bin"})
                .Foreach([&hrmslayer, &hefraclayer, &bidx](const std::vector<double> rms, const std::vector<double> efrac, const double energyw)
                {
                    for (unsigned int lidx=0; lidx<bgolayers; ++lidx)
                    {
                        hrmslayer[bidx][lidx]->Fill(rms[lidx], energyw);
                        hefraclayer[bidx][lidx]->Fill(efrac[lidx], energyw);
                    }
                }, {"rmsLayer_gauss", "fracLayer_gauss", "simu_energy_w_corr"});
    }

    TFile* out = TFile::Open(input_args.output_path.c_str(), "RECREATE");
    if (out->IsZombie())
    {
        std::cerr << "\n\nError writing output ROOT file: [" << input_args.output_path << "]\n\n";
        exit(100);
    }

    for (unsigned int bidx=0; bidx<nenergybin; ++bidx)
    {
        out->mkdir((std::string("RMS/energybin_") + std::to_string(bidx+1)).c_str());
        out->cd((std::string("RMS/energybin_") + std::to_string(bidx+1)).c_str());
        for (unsigned int lidx=0; lidx<bgolayers; ++lidx)
            hrmslayer[bidx][lidx]->Write();
    
        out->mkdir((std::string("ELF/energybin_") + std::to_string(bidx+1)).c_str());
        out->cd((std::string("ELF/energybin_") + std::to_string(bidx+1)).c_str());
        for (unsigned int lidx=0; lidx<bgolayers; ++lidx)
            hefraclayer[bidx][lidx]->Write();
    }
    
    out->Close();
}

