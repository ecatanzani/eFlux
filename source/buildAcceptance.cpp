#include "myHeader.h"

void buildAcceptance(TH1D &acceptance, TFile &outFile)
{
    TH1D* gFactor = nullptr;
    TH1D* selEfficiency_beforeSelection = nullptr;
    TH1D* selEfficiency_alfterSelection = nullptr;

    const char* geoFilePath = "effHistos/Acc_vs_Ekin_fid_ele_v600.root";
    const char* selFilePath = "effHistos/hist_tight_DNN_histograms_allElectron-v6r0p0_1GeV_10TeV.root";

    TFile geoFile(geoFilePath,"READ");
    if(!geoFile.IsOpen())
    {
        std::cerr << "\n\nError opening geometric factor TFile: " << geoFilePath << std::endl;
        exit(123);
    }
    
    // Getting MC Geometrical Factor from file
    //geoFile.GetObject("Acc",gFactor);
    gFactor = (TH1D*) geoFile.Get("Acc");
    gFactor->SetDirectory(0);
    geoFile.Close();
    fitGFactor(gFactor,outFile);
    
    TFile selFile(selFilePath,"READ");
    if(!selFile.IsOpen())
    {
        std::cerr << "\n\nError opening efficiency selection TFile: " << selFilePath << std::endl;
        exit(123);
    }
    //selFile.GetObject("BgoTotalEcorr_fidBgo_log_L",selEfficiency_beforeSelection);
    //selFile.GetObject("BgoTotalEcorr_sel_allCases_xtrLoose_log_L",selEfficiency_alfterSelection);
    selEfficiency_beforeSelection = (TH1D*) selFile.Get("BgoTotalEcorr_fidBgo_log_L");
    selEfficiency_alfterSelection = (TH1D*) selFile.Get("BgoTotalEcorr_sel_allCases_xtrLoose_log_L");
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
    //gFactor->Multiply(selEfficiency_alfterSelection);

    new (&acceptance) (TH1D) (*(TH1D*)gFactor->Clone("Acceptance"));
    

}