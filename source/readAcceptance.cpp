#include "myHeader.h"
#include "acceptanceFit.h"

TF1 readAcceptance(
    TFile &outFile, 
    const bool verbose,
    const bool myAcceptance)
{
    // Reading acceptance
    if (myAcceptance)
    {
        auto accDir = outFile.GetDirectory("acceptance");
        auto acceptanceGr = (static_cast<TH1D*> (accDir->Get("acceptance")));
        auto fitter = fitAcceptance<TH1D>(acceptanceGr,outFile,verbose);
        return fitter;
    }
    else
    {
        const char *geoFilePath = "effHistos/Acc_vs_Ekin_fid_ele_v600.root";
        const char *selFilePath = "effHistos/hist_tight_DNN_histograms_allElectron-v6r0p0_1GeV_10TeV.root";
        // Reading acceptance from external file
        TFile geoFile(geoFilePath, "READ");
        if (!geoFile.IsOpen())
        {
            std::cerr << "\n\nError opening geometric factor TFile: " << geoFilePath << std::endl;
            exit(123);
        }
        auto gFactor = (static_cast<TH1F*> (geoFile.Get("Acc")));
        gFactor->SetDirectory(0);
        geoFile.Close();
        // Fit input Acceptance function
        auto fitter = fitAcceptance<TH1F>(gFactor, outFile, verbose);
        // Reading efficiencies from external file
        TFile selFile(selFilePath, "READ");
        if (!selFile.IsOpen())
        {
            std::cerr << "\n\nError opening efficiency selection TFile: " << selFilePath << std::endl;
            exit(123);
        }
        auto selEfficiency_beforeSelection = (static_cast<TH1D*> (selFile.Get("BgoTotalEcorr_fidBgo_log_L")));
        auto selEfficiency_alfterSelection = (static_cast<TH1D*> (selFile.Get("BgoTotalEcorr_sel_allCases_xtrLoose_log_L")));
        selEfficiency_beforeSelection->SetDirectory(0);
        selEfficiency_alfterSelection->SetDirectory(0);
        selFile.Close();
        /*
        Sumw2() just applied to each single histo when they have been computed. WARNING !!! This changes if other files will be used

        gFactor->Sumw2();
        selEfficiency_beforeSelection->Sumw2();
        selEfficiency_alfterSelection->Sumw2();
        */
        // Compute the final efficiency
        auto selEff = (static_cast<TH1D*> (selEfficiency_alfterSelection->Clone("selEff")));
        selEff->Divide(selEfficiency_beforeSelection);
        selEff->Multiply(gFactor->GetFunction(fitter.GetName()));
        auto acceptanceGr = (static_cast<TH1D*> (selEff->Clone("acceptance")));
        fitter = fitAcceptance<TH1D>(acceptanceGr,outFile,verbose);
         // Creating a TDirectory for the acceptance histo
        TDirectory *aDir = outFile.mkdir("Acceptance");
        aDir->cd();
        selEfficiency_beforeSelection->Write();
        selEfficiency_alfterSelection->Write();
        selEff->Write();
        acceptanceGr->Write();
        // Returning to main dir on the output TFile
        outFile.cd();
        return fitter;
    }
}