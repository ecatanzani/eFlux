#include "myHeader.h"
#include "acceptanceFit.h"

void readAcceptance(TH1D &acceptance, TFile &outFile, const bool verbose)
{

#if 1

    const char *accPath = "effHistos/acceptance.root";
    TFile accFile(accPath, "READ");
    if (!accFile.IsOpen())
    {
        std::cerr << "\n\nError opening geometric factor TFile: " << accPath << std::endl;
        exit(123);
    }

    new (&acceptance)(TH1D)(*(TH1D *)accFile.Get("Acc"));
    accFile.Close();

    // Creating a TDirectory for the acceptance histo
    TDirectory *aDir = outFile.mkdir("Acceptance");
    aDir->cd();

    acceptance.Write();

    // Returning to main dir on the output TFile
    outFile.cd();

#else

    TH1D *gFactor = nullptr;
    TH1D *selEfficiency_beforeSelection = nullptr;
    TH1D *selEfficiency_alfterSelection = nullptr;
    TH1D *selEff = nullptr;

    const char *geoFilePath = "effHistos/Acc_vs_Ekin_fid_ele_v600.root";
    const char *selFilePath = "effHistos/hist_tight_DNN_histograms_allElectron-v6r0p0_1GeV_10TeV.root";

    TFile geoFile(geoFilePath, "READ");
    if (!geoFile.IsOpen())
    {
        std::cerr << "\n\nError opening geometric factor TFile: " << geoFilePath << std::endl;
        exit(123);
    }

    // Getting MC Geometrical Factor from file
    //geoFile.GetObject("Acc",gFactor);
    gFactor = (TH1D *)geoFile.Get("Acc");
    gFactor->SetDirectory(0);
    geoFile.Close();
    auto fitter = fitGFactor(gFactor, outFile, verbose);

    TFile selFile(selFilePath, "READ");
    if (!selFile.IsOpen())
    {
        std::cerr << "\n\nError opening efficiency selection TFile: " << selFilePath << std::endl;
        exit(123);
    }
    //selFile.GetObject("BgoTotalEcorr_fidBgo_log_L",selEfficiency_beforeSelection);
    //selFile.GetObject("BgoTotalEcorr_sel_allCases_xtrLoose_log_L",selEfficiency_alfterSelection);
    selEfficiency_beforeSelection = (TH1D *)selFile.Get("BgoTotalEcorr_fidBgo_log_L");
    selEfficiency_alfterSelection = (TH1D *)selFile.Get("BgoTotalEcorr_sel_allCases_xtrLoose_log_L");
    selEfficiency_beforeSelection->SetDirectory(0);
    selEfficiency_alfterSelection->SetDirectory(0);
    selFile.Close();

    /*
        Sumw2() just applied to each single histo

    gFactor->Sumw2();
    selEfficiency_beforeSelection->Sumw2();
    selEfficiency_alfterSelection->Sumw2();
    */

    selEfficiency_alfterSelection->Divide(selEfficiency_beforeSelection);
    selEff = (TH1D *)selEfficiency_alfterSelection->Clone("selEff");
    selEfficiency_alfterSelection->Multiply(gFactor->GetFunction(fitter.GetName()));

    new (&acceptance)(TH1D)(*(TH1D *)selEfficiency_alfterSelection->Clone("Acceptance"));

    // Creating a TDirectory for the acceptance histo
    TDirectory *aDir = outFile.mkdir("Acceptance");
    aDir->cd();

    selEfficiency_beforeSelection->Write();
    selEfficiency_alfterSelection->Write();
    selEff->Write();
    acceptance.Write();

    // Returning to main dir on the output TFile
    outFile.cd();

#endif
}