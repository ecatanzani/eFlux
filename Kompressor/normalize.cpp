#include <string>
#include <vector>
#include <memory>
#include <iostream>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TProfile.h"
#include "TDirectory.h"

inline std::vector<std::shared_ptr<TH1D>> projection(TH2D* histo)
{
    std::vector<std::shared_ptr<TH1D>> projections (histo->GetNbinsX());
    for (int ibin=1; ibin<=histo->GetNbinsX(); ++ibin)
        projections[ibin-1] = std::shared_ptr<TH1D>(histo->ProjectionY((std::string("projectionY_") + std::string(histo->GetName()) + std::string("_bin_") + std::to_string(ibin)).c_str(), ibin, ibin));
    return projections;
}

void SliceNormalization(
    const char* full_histo_path, 
    const char* normalized_histo_path)
{
    auto slicenorm = [](TH2D* histo) { 
        for (int bidx=1; bidx<=histo->GetNbinsX(); ++bidx)
	{
		auto scale_factor = 1/histo->Integral(bidx, bidx, 1, histo->GetNbinsY());
		for (int bidy=1; bidy<=histo->GetNbinsY(); ++bidy)
		{
            		histo->SetBinContent(bidx, bidy, histo->GetBinContent(bidx, bidy)*scale_factor);
			histo->SetBinError(bidx, bidy, histo->GetBinError(bidx, bidy)*scale_factor);
        	}
	}
	return histo; };

    TFile *_input = TFile::Open(full_histo_path, "READ");
    if (_input->IsZombie())
    {
        std::cerr << "\n\nError opening input file: [" << full_histo_path << "]\n\n";
        exit(100);
    }

    // Simu histos
    auto h_energy_diff2D = slicenorm(static_cast<TH2D*>(_input->Get("Simu/h_energy_diff2D")));
    auto h_energy_diff2D_corr = slicenorm(static_cast<TH2D*>(_input->Get("Simu/h_energy_diff2D_corr")));
    auto h_energy_unfold = slicenorm(static_cast<TH2D*>(_input->Get("Simu/h_energy_unfold")));
    auto h_energy_unfold_corr = slicenorm(static_cast<TH2D*>(_input->Get("Simu/h_energy_unfold_corr")));

    auto projectionY_energy_diff2D = projection(h_energy_diff2D);
    auto projectionY_energy_diff2D_corr = projection(h_energy_diff2D_corr);   

    h_energy_diff2D->SetDirectory(0);
    h_energy_diff2D_corr->SetDirectory(0);
    h_energy_unfold->SetDirectory(0);
    h_energy_unfold_corr->SetDirectory(0);

    // XTRL histo
    auto h_xtrl = slicenorm(static_cast<TH2D*>(_input->Get("BGO/h_xtrl")));

    h_xtrl->SetDirectory(0);
    
    TFile *_output = TFile::Open(normalized_histo_path, "RECREATE");
    if (_output->IsZombie())
    {
        std::cerr << "\n\nError writing output file: [" << normalized_histo_path << "]\n\n";
        exit(100);
    }

    auto simu_dir = _output->mkdir("Simu");
    simu_dir->cd();

    h_energy_diff2D->Write();
    h_energy_diff2D_corr->Write();
    h_energy_unfold->Write();
    h_energy_unfold_corr->Write();

    auto profile_dir = _output->mkdir("Simu/Projections");
    _output->cd("Simu/Projections");
    for (auto& _elm : projectionY_energy_diff2D)
        _elm->Write();
    for (auto& _elm : projectionY_energy_diff2D_corr)
        _elm->Write();

    auto xtrl_dir = _output->mkdir("xtrl");
    xtrl_dir->cd();

    h_xtrl->Write();

    _output->Close();

}
