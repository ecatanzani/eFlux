#include "main.h"
#include "gaus.h"
#include "parser.h"
#include "histos.h"

#include "TCanvas.h"

void resumegaus(in_args input_args)
{
    auto hrmslayer = getrmslayerhistos_ff(input_args.input_list);
    auto hefraclayer = getenergyfractionlayerhistos_ff(input_args.input_list);

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
    
    std::vector<std::shared_ptr<TCanvas>> c_rms (nenergybin);
    std::vector<std::shared_ptr<TCanvas>> c_efraction (nenergybin);

    for (unsigned int bidx=0; bidx<nenergybin; ++bidx)
    {
        c_rms[bidx] = std::make_shared<TCanvas> ((std::string("rms_energybin_") + std::to_string(bidx+1)).c_str(), (std::string("RMS distributions - Energy bin ") + std::to_string(bidx+1)).c_str());
        c_rms[bidx]->Divide(7,2);
        for (unsigned int lidx=0; lidx<bgolayers; ++lidx)
        {
            c_rms[bidx]->cd(lidx+1);
            hrmslayer[bidx][lidx]->Draw();
        }
        out->cd((std::string("RMS/energybin_") + std::to_string(bidx+1)).c_str());
        c_rms[bidx]->cd(0);
        c_rms[bidx]->Write();

        c_efraction[bidx] = std::make_shared<TCanvas> ((std::string("elf_energybin_") + std::to_string(bidx+1)).c_str(), (std::string("Energy fraction distributions - Energy bin ") + std::to_string(bidx+1)).c_str());
        c_efraction[bidx]->Divide(7,2);
        for (unsigned int lidx=0; lidx<bgolayers; ++lidx)
        {
            c_efraction[bidx]->cd(lidx+1);
            hefraclayer[bidx][lidx]->Draw();
        }
        out->cd((std::string("ELF/energybin_") + std::to_string(bidx+1)).c_str());
        c_efraction[bidx]->cd(0);
        c_efraction[bidx]->Write();
    }

    out->Close();
}

