#include "train_utils.h"

std::map<std::string, int> GetTMVAMethods(const std::vector<std::string> mymethod) {
    std::map<std::string, int> Use;

    Use["Cuts"] = 0;
    Use["CutsD"] = 0;
    Use["CutsPCA"] = 0;
    Use["CutsGA"] = 0;
    Use["CutsSA"] = 0;

    Use["Likelihood"] = 0;
    Use["LikelihoodD"] = 0;
    Use["LikelihoodPCA"] = 0;
    Use["LikelihoodKDE"] = 0;
    Use["LikelihoodMIX"] = 0;

    Use["PDERS"] = 0;
    Use["PDERSD"] = 0;
    Use["PDERSPCA"] = 0;
    Use["PDEFoam"] = 0;
    Use["PDEFoamBoost"] = 0;
    Use["KNN"] = 0;

    Use["LD"] = 0;
    Use["Fisher"] = 0;
    Use["FisherG"] = 0;
    Use["BoostedFisher"] = 0;
    Use["HMatrix"] = 0;

    Use["FDA_GA"] = 0;
    Use["FDA_SA"] = 0;
    Use["FDA_MC"] = 0;
    Use["FDA_MT"] = 0;
    Use["FDA_GAMT"] = 0;
    Use["FDA_MCMT"] = 0;

    Use["MLP"] = 0;
    Use["MLPBFGS"] = 0;
    Use["MLPBNN"] = 0;
    Use["CFMlpANN"] = 0;
    Use["TMlpANN"] = 0;

    Use["SVM"] = 0;

    Use["BDT"] = 0;
    Use["BDTG"] = 0;
    Use["BDTB"] = 0;
    Use["BDTD"] = 0;
    Use["BDTF"] = 0;

    Use["RuleFit"] = 0;

    auto linked_method = false;
    for (auto &&dmethod : Use)
        for (const auto& i_mtd: mymethod)
            if (!strcmp(i_mtd.c_str(), dmethod.first.c_str())) {
                dmethod.second = 1;
                linked_method = true;
                break;
            }

    if (!linked_method) {
        std::cerr << "\nERROR: No match found in TMVA default methods...\n\n";
        exit(100);
    }

    return Use;
}

void SetTMVAVariables(std::shared_ptr<TMVA::DataLoader> dataloader, const train_vars vars) {
    SetAllTMVAVariables(dataloader);
#if 0    
    if (vars.all_vars) SetAllTMVAVariables(dataloader);
    else if (vars.no_nud) SetNoNUDTMVAVariables(dataloader);
    else if (vars.nud_only) SetNUDTMVAVariables(dataloader);
#endif
}

void SetAllTMVAVariables(std::shared_ptr<TMVA::DataLoader> dataloader)
{
    dataloader->AddVariable("rmslayer_norm_1", "rmslayer_norm_1", "units", 'D');
    dataloader->AddVariable("rmslayer_norm_2", "rmslayer_norm_2", "units", 'D');
    dataloader->AddVariable("rmslayer_norm_3", "rmslayer_norm_3", "units", 'D');
    dataloader->AddVariable("rmslayer_norm_4", "rmslayer_norm_4", "units", 'D');
    dataloader->AddVariable("rmslayer_norm_5", "rmslayer_norm_5", "units", 'D');
    dataloader->AddVariable("rmslayer_norm_6", "rmslayer_norm_6", "units", 'D');
    dataloader->AddVariable("rmslayer_norm_7", "rmslayer_norm_7", "units", 'D');
    dataloader->AddVariable("rmslayer_norm_8", "rmslayer_norm_8", "units", 'D');
    dataloader->AddVariable("rmslayer_norm_9", "rmslayer_norm_9", "units", 'D');
    dataloader->AddVariable("rmslayer_norm_10", "rmslayer_norm_10", "units", 'D');
    dataloader->AddVariable("rmslayer_norm_11", "rmslayer_norm_11", "units", 'D');
    dataloader->AddVariable("rmslayer_norm_12", "rmslayer_norm_12", "units", 'D');
    dataloader->AddVariable("rmslayer_norm_13", "rmslayer_norm_13", "units", 'D');
    dataloader->AddVariable("rmslayer_norm_14", "rmslayer_norm_14", "units", 'D');
    
    dataloader->AddVariable("fraclayer_norm_1", "fraclayer_norm_1", "units", 'D');
    dataloader->AddVariable("fraclayer_norm_2", "fraclayer_norm_2", "units", 'D');
    dataloader->AddVariable("fraclayer_norm_3", "fraclayer_norm_3", "units", 'D');
    dataloader->AddVariable("fraclayer_norm_4", "fraclayer_norm_4", "units", 'D');
    dataloader->AddVariable("fraclayer_norm_5", "fraclayer_norm_5", "units", 'D');
    dataloader->AddVariable("fraclayer_norm_6", "fraclayer_norm_6", "units", 'D');
    dataloader->AddVariable("fraclayer_norm_7", "fraclayer_norm_7", "units", 'D');
    dataloader->AddVariable("fraclayer_norm_8", "fraclayer_norm_8", "units", 'D');
    dataloader->AddVariable("fraclayer_norm_9", "fraclayer_norm_9", "units", 'D');
    dataloader->AddVariable("fraclayer_norm_10", "fraclayer_norm_10", "units", 'D');
    dataloader->AddVariable("fraclayer_norm_11", "fraclayer_norm_11", "units", 'D');
    dataloader->AddVariable("fraclayer_norm_12", "fraclayer_norm_12", "units", 'D');
    dataloader->AddVariable("fraclayer_norm_13", "fraclayer_norm_13", "units", 'D');
    dataloader->AddVariable("fraclayer_norm_14", "fraclayer_norm_14", "units", 'D');
    
    dataloader->AddVariable("sumrms_norm", "sumrms_norm", "units", 'D');
    dataloader->AddVariable("fraclastlayer_norm", "fraclastlayer_norm", "units", 'D');
    dataloader->AddVariable("xtrl_norm", "xtrl_norm", "units", 'D');

    dataloader->AddSpectator("xtrl", "xtrl", "units", 'D');
}

extern void SetNoNUDTMVAVariables(std::shared_ptr<TMVA::DataLoader> dataloader)
{
    dataloader->AddVariable("STK_bestTrack_npoints", "STK_bestTrack_npoints", "units", 'I');
    //dataloader->AddVariable("sumRms_reg", "sumRms_reg", "units", 'D');
    //dataloader->AddVariable("fracLast_reg", "fracLast_reg", "units", 'D');

    dataloader->AddVariable("rmsLayer_1", "rmsLayer_1", "units", 'D');
    dataloader->AddVariable("rmsLayer_2", "rmsLayer_2", "units", 'D');
    dataloader->AddVariable("rmsLayer_3", "rmsLayer_3", "units", 'D');
    dataloader->AddVariable("rmsLayer_4", "rmsLayer_4", "units", 'D');
    dataloader->AddVariable("rmsLayer_5", "rmsLayer_5", "units", 'D');
    dataloader->AddVariable("rmsLayer_6", "rmsLayer_6", "units", 'D');
    dataloader->AddVariable("rmsLayer_7", "rmsLayer_7", "units", 'D');
    dataloader->AddVariable("rmsLayer_8", "rmsLayer_8", "units", 'D');
    dataloader->AddVariable("rmsLayer_9", "rmsLayer_9", "units", 'D');
    dataloader->AddVariable("rmsLayer_10", "rmsLayer_10", "units", 'D');
    dataloader->AddVariable("rmsLayer_11", "rmsLayer_11", "units", 'D');
    dataloader->AddVariable("rmsLayer_12", "rmsLayer_12", "units", 'D');
    dataloader->AddVariable("rmsLayer_13", "rmsLayer_13", "units", 'D');
    dataloader->AddVariable("rmsLayer_14", "rmsLayer_14", "units", 'D');
    dataloader->AddVariable("fracLayer_1", "fracLayer_1", "units", 'D');
    dataloader->AddVariable("fracLayer_2", "fracLayer_2", "units", 'D');
    dataloader->AddVariable("fracLayer_3", "fracLayer_3", "units", 'D');
    dataloader->AddVariable("fracLayer_4", "fracLayer_4", "units", 'D');
    dataloader->AddVariable("fracLayer_5", "fracLayer_5", "units", 'D');
    dataloader->AddVariable("fracLayer_6", "fracLayer_6", "units", 'D');
    dataloader->AddVariable("fracLayer_7", "fracLayer_7", "units", 'D');
    dataloader->AddVariable("fracLayer_8", "fracLayer_8", "units", 'D');
    dataloader->AddVariable("fracLayer_9", "fracLayer_9", "units", 'D');
    dataloader->AddVariable("fracLayer_10", "fracLayer_10", "units", 'D');
    dataloader->AddVariable("fracLayer_11", "fracLayer_11", "units", 'D');
    dataloader->AddVariable("fracLayer_12", "fracLayer_12", "units", 'D');
    dataloader->AddVariable("fracLayer_13", "fracLayer_13", "units", 'D');
    dataloader->AddVariable("fracLayer_14", "fracLayer_14", "units", 'D');

    //dataloader->AddVariable("lastBGOLayer", "lastBGOLayer", "units", 'I');
    dataloader->AddVariable("nBGOentries", "nBGOentries", "units", 'I');

    dataloader->AddVariable("energy_1R_radius_1", "energy_1R_radius_1", "units", 'D');
    dataloader->AddVariable("energy_1R_radius_2", "energy_1R_radius_2", "units", 'D');
    dataloader->AddVariable("energy_1R_radius_3", "energy_1R_radius_3", "units", 'D');
    dataloader->AddVariable("energy_1R_radius_4", "energy_1R_radius_4", "units", 'D');
    dataloader->AddVariable("energy_1R_radius_5", "energy_1R_radius_5", "units", 'D');
    dataloader->AddVariable("energy_1R_radius_6", "energy_1R_radius_6", "units", 'D');
    dataloader->AddVariable("energy_1R_radius_7", "energy_1R_radius_7", "units", 'D');
    dataloader->AddVariable("energy_1R_radius_8", "energy_1R_radius_8", "units", 'D');
    dataloader->AddVariable("energy_1R_radius_9", "energy_1R_radius_9", "units", 'D');
    dataloader->AddVariable("energy_1R_radius_10", "energy_1R_radius_10", "units", 'D');
    dataloader->AddVariable("energy_1R_radius_11", "energy_1R_radius_11", "units", 'D');
    dataloader->AddVariable("energy_1R_radius_12", "energy_1R_radius_12", "units", 'D');
    dataloader->AddVariable("energy_1R_radius_13", "energy_1R_radius_13", "units", 'D');
    dataloader->AddVariable("energy_1R_radius_14", "energy_1R_radius_14", "units", 'D');

    dataloader->AddVariable("energy_2R_radius_1", "energy_2R_radius_1", "units", 'D');
    dataloader->AddVariable("energy_2R_radius_2", "energy_2R_radius_2", "units", 'D');
    dataloader->AddVariable("energy_2R_radius_3", "energy_2R_radius_3", "units", 'D');
    dataloader->AddVariable("energy_2R_radius_4", "energy_2R_radius_4", "units", 'D');
    dataloader->AddVariable("energy_2R_radius_5", "energy_2R_radius_5", "units", 'D');
    dataloader->AddVariable("energy_2R_radius_6", "energy_2R_radius_6", "units", 'D');
    dataloader->AddVariable("energy_2R_radius_7", "energy_2R_radius_7", "units", 'D');
    dataloader->AddVariable("energy_2R_radius_8", "energy_2R_radius_8", "units", 'D');
    dataloader->AddVariable("energy_2R_radius_9", "energy_2R_radius_9", "units", 'D');
    dataloader->AddVariable("energy_2R_radius_10", "energy_2R_radius_10", "units", 'D');
    dataloader->AddVariable("energy_2R_radius_11", "energy_2R_radius_11", "units", 'D');
    dataloader->AddVariable("energy_2R_radius_12", "energy_2R_radius_12", "units", 'D');
    dataloader->AddVariable("energy_2R_radius_13", "energy_2R_radius_13", "units", 'D');
    dataloader->AddVariable("energy_2R_radius_14", "energy_2R_radius_14", "units", 'D');

    dataloader->AddVariable("energy_3R_radius_1", "energy_3R_radius_1", "units", 'D');
    dataloader->AddVariable("energy_3R_radius_2", "energy_3R_radius_2", "units", 'D');
    dataloader->AddVariable("energy_3R_radius_3", "energy_3R_radius_3", "units", 'D');
    dataloader->AddVariable("energy_3R_radius_4", "energy_3R_radius_4", "units", 'D');
    dataloader->AddVariable("energy_3R_radius_5", "energy_3R_radius_5", "units", 'D');
    dataloader->AddVariable("energy_3R_radius_6", "energy_3R_radius_6", "units", 'D');
    dataloader->AddVariable("energy_3R_radius_7", "energy_3R_radius_7", "units", 'D');
    dataloader->AddVariable("energy_3R_radius_8", "energy_3R_radius_8", "units", 'D');
    dataloader->AddVariable("energy_3R_radius_9", "energy_3R_radius_9", "units", 'D');
    dataloader->AddVariable("energy_3R_radius_10", "energy_3R_radius_10", "units", 'D');
    dataloader->AddVariable("energy_3R_radius_11", "energy_3R_radius_11", "units", 'D');
    dataloader->AddVariable("energy_3R_radius_12", "energy_3R_radius_12", "units", 'D');
    dataloader->AddVariable("energy_3R_radius_13", "energy_3R_radius_13", "units", 'D');
    dataloader->AddVariable("energy_3R_radius_14", "energy_3R_radius_14", "units", 'D');

    dataloader->AddVariable("energy_5R_radius_1", "energy_5R_radius_1", "units", 'D');
    dataloader->AddVariable("energy_5R_radius_2", "energy_5R_radius_2", "units", 'D');
    dataloader->AddVariable("energy_5R_radius_3", "energy_5R_radius_3", "units", 'D');
    dataloader->AddVariable("energy_5R_radius_4", "energy_5R_radius_4", "units", 'D');
    dataloader->AddVariable("energy_5R_radius_5", "energy_5R_radius_5", "units", 'D');
    dataloader->AddVariable("energy_5R_radius_6", "energy_5R_radius_6", "units", 'D');
    dataloader->AddVariable("energy_5R_radius_7", "energy_5R_radius_7", "units", 'D');
    dataloader->AddVariable("energy_5R_radius_8", "energy_5R_radius_8", "units", 'D');
    dataloader->AddVariable("energy_5R_radius_9", "energy_5R_radius_9", "units", 'D');
    dataloader->AddVariable("energy_5R_radius_10", "energy_5R_radius_10", "units", 'D');
    dataloader->AddVariable("energy_5R_radius_11", "energy_5R_radius_11", "units", 'D');
    dataloader->AddVariable("energy_5R_radius_12", "energy_5R_radius_12", "units", 'D');
    dataloader->AddVariable("energy_5R_radius_13", "energy_5R_radius_13", "units", 'D');
    dataloader->AddVariable("energy_5R_radius_14", "energy_5R_radius_14", "units", 'D');

    dataloader->AddVariable("xtrl", "xtrl", "units", 'D');
}

void SetNUDTMVAVariables(std::shared_ptr<TMVA::DataLoader> dataloader)
{
    dataloader->AddVariable("NUD_total_ADC_nud_total_adc", "NUD_total_ADC_nud_total_adc", "units", 'I');
    dataloader->AddVariable("NUD_ADC_min", "NUD_ADC_min", "units", 'I');
    dataloader->AddVariable("NUD_ADC_max", "NUD_ADC_max", "units", 'I');
    dataloader->AddVariable("NUD_ADC_rms", "NUD_ADC_rms", "units", 'D');
}

void SetTMVACuts(
    TCut &signal_cuts, 
    TCut &background_cuts, 
    const bool verbose)
{
#if 0    
    TCut energyrange = "";
    signal_cuts += energyrange;
    background_cuts += energyrange;
#endif

    // Remove bad events values on fraclastlayer_norm valiable (maybe due to gaussianization process)
    TCut fraclastlayer_norm_errs = "!(fraclastlayer_norm<-100 || fraclastlayer_norm>100)";
    signal_cuts+=fraclastlayer_norm_errs;
    background_cuts+=fraclastlayer_norm_errs;

    // Remove infinite values and nans
    TCut removenans_rmsl1 = "!(TMath::IsNaN(rmslayer_norm_1) || !(TMath::Finite(rmslayer_norm_1)))";
    signal_cuts+=removenans_rmsl1;
    background_cuts+=removenans_rmsl1;

    TCut removenans_rmsl2 = "!(TMath::IsNaN(rmslayer_norm_2) || !(TMath::Finite(rmslayer_norm_2)))";
    signal_cuts+=removenans_rmsl2;
    background_cuts+=removenans_rmsl2;

    TCut removenans_rmsl3 = "!(TMath::IsNaN(rmslayer_norm_3) || !(TMath::Finite(rmslayer_norm_3)))";
    signal_cuts+=removenans_rmsl3;
    background_cuts+=removenans_rmsl3;

    TCut removenans_rmsl4 = "!(TMath::IsNaN(rmslayer_norm_4) || !(TMath::Finite(rmslayer_norm_4)))";
    signal_cuts+=removenans_rmsl4;
    background_cuts+=removenans_rmsl4;

    TCut removenans_rmsl5 = "!(TMath::IsNaN(rmslayer_norm_5) || !(TMath::Finite(rmslayer_norm_5)))";
    signal_cuts+=removenans_rmsl5;
    background_cuts+=removenans_rmsl5;

    TCut removenans_rmsl6 = "!(TMath::IsNaN(rmslayer_norm_6) || !(TMath::Finite(rmslayer_norm_6)))";
    signal_cuts+=removenans_rmsl6;
    background_cuts+=removenans_rmsl6;

    TCut removenans_rmsl7 = "!(TMath::IsNaN(rmslayer_norm_7) || !(TMath::Finite(rmslayer_norm_7)))";
    signal_cuts+=removenans_rmsl7;
    background_cuts+=removenans_rmsl7;

    TCut removenans_rmsl8 = "!(TMath::IsNaN(rmslayer_norm_8) || !(TMath::Finite(rmslayer_norm_8)))";
    signal_cuts+=removenans_rmsl8;
    background_cuts+=removenans_rmsl8;

    TCut removenans_rmsl9 = "!(TMath::IsNaN(rmslayer_norm_9) || !(TMath::Finite(rmslayer_norm_9)))";
    signal_cuts+=removenans_rmsl9;
    background_cuts+=removenans_rmsl9;

    TCut removenans_rmsl10 = "!(TMath::IsNaN(rmslayer_norm_10) || !(TMath::Finite(rmslayer_norm_10)))";
    signal_cuts+=removenans_rmsl10;
    background_cuts+=removenans_rmsl10;

    TCut removenans_rmsl11 = "!(TMath::IsNaN(rmslayer_norm_11) || !(TMath::Finite(rmslayer_norm_11)))";
    signal_cuts+=removenans_rmsl11;
    background_cuts+=removenans_rmsl11;

    TCut removenans_rmsl12 = "!(TMath::IsNaN(rmslayer_norm_12) || !(TMath::Finite(rmslayer_norm_12)))";
    signal_cuts+=removenans_rmsl12;
    background_cuts+=removenans_rmsl12;

    TCut removenans_rmsl13 = "!(TMath::IsNaN(rmslayer_norm_13) || !(TMath::Finite(rmslayer_norm_13)))";
    signal_cuts+=removenans_rmsl13;
    background_cuts+=removenans_rmsl13;

    TCut removenans_rmsl14 = "!(TMath::IsNaN(rmslayer_norm_14) || !(TMath::Finite(rmslayer_norm_14)))";
    signal_cuts+=removenans_rmsl14;
    background_cuts+=removenans_rmsl14;

    TCut removenans_fracl1 = "!(TMath::IsNaN(fraclayer_norm_1) || !(TMath::Finite(fraclayer_norm_1)))";
    signal_cuts+=removenans_fracl1;
    background_cuts+=removenans_fracl1;

    TCut removenans_fracl2 = "!(TMath::IsNaN(fraclayer_norm_2) || !(TMath::Finite(fraclayer_norm_2)))";
    signal_cuts+=removenans_fracl2;
    background_cuts+=removenans_fracl2;

    TCut removenans_fracl3 = "!(TMath::IsNaN(fraclayer_norm_3) || !(TMath::Finite(fraclayer_norm_3)))";
    signal_cuts+=removenans_fracl3;
    background_cuts+=removenans_fracl3;

    TCut removenans_fracl4 = "!(TMath::IsNaN(fraclayer_norm_4) || !(TMath::Finite(fraclayer_norm_4)))";
    signal_cuts+=removenans_fracl4;
    background_cuts+=removenans_fracl4;

    TCut removenans_fracl5 = "!(TMath::IsNaN(fraclayer_norm_5) || !(TMath::Finite(fraclayer_norm_5)))";
    signal_cuts+=removenans_fracl5;
    background_cuts+=removenans_fracl5;

    TCut removenans_fracl6 = "!(TMath::IsNaN(fraclayer_norm_6) || !(TMath::Finite(fraclayer_norm_6)))";
    signal_cuts+=removenans_fracl6;
    background_cuts+=removenans_fracl6;

    TCut removenans_fracl7 = "!(TMath::IsNaN(fraclayer_norm_7) || !(TMath::Finite(fraclayer_norm_7)))";
    signal_cuts+=removenans_fracl7;
    background_cuts+=removenans_fracl7;

    TCut removenans_fracl8 = "!(TMath::IsNaN(fraclayer_norm_8) || !(TMath::Finite(fraclayer_norm_8)))";
    signal_cuts+=removenans_fracl8;
    background_cuts+=removenans_fracl8;

    TCut removenans_fracl9 = "!(TMath::IsNaN(fraclayer_norm_9) || !(TMath::Finite(fraclayer_norm_9)))";
    signal_cuts+=removenans_fracl9;
    background_cuts+=removenans_fracl9;

    TCut removenans_fracl10 = "!(TMath::IsNaN(fraclayer_norm_10) || !(TMath::Finite(fraclayer_norm_10)))";
    signal_cuts+=removenans_fracl10;
    background_cuts+=removenans_fracl10;

    TCut removenans_fracl11 = "!(TMath::IsNaN(fraclayer_norm_11) || !(TMath::Finite(fraclayer_norm_11)))";
    signal_cuts+=removenans_fracl11;
    background_cuts+=removenans_fracl11;

    TCut removenans_fracl12 = "!(TMath::IsNaN(fraclayer_norm_12) || !(TMath::Finite(fraclayer_norm_12)))";
    signal_cuts+=removenans_fracl12;
    background_cuts+=removenans_fracl12;

    TCut removenans_fracl13 = "!(TMath::IsNaN(fraclayer_norm_13) || !(TMath::Finite(fraclayer_norm_13)))";
    signal_cuts+=removenans_fracl13;
    background_cuts+=removenans_fracl13;

    TCut removenans_fracl14 = "!(TMath::IsNaN(fraclayer_norm_14) || !(TMath::Finite(fraclayer_norm_14)))";
    signal_cuts+=removenans_fracl14;
    background_cuts+=removenans_fracl14;

    TCut removenans_sumrms = "!(TMath::IsNaN(sumrms_norm) || !(TMath::Finite(sumrms_norm)))";
    signal_cuts+=removenans_sumrms;
    background_cuts+=removenans_sumrms;

    TCut removenans_fracll = "!(TMath::IsNaN(fraclastlayer_norm) || !(TMath::Finite(fraclastlayer_norm)))";
    signal_cuts+=removenans_fracll;
    background_cuts+=removenans_fracll;

    TCut removenans_xtrl = "!(TMath::IsNaN(xtrl_norm) || !(TMath::Finite(xtrl_norm)))";
    signal_cuts+=removenans_xtrl;
    background_cuts+=removenans_xtrl;
    
    if (verbose)
    {
        std::cout << "\n\nSignal cuts....\n\n";
        signal_cuts.Print();
        std::cout << "\n\nBackground cuts....\n\n";
        background_cuts.Print();
        std::cout << "\n\n";
    }
}

void BookMethods(
    std::shared_ptr<TMVA::Factory> factory,
    std::shared_ptr<TMVA::DataLoader> dataloader,
    std::map<std::string, int> Use)
{
    // Cut optimisation
    if (Use["Cuts"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kCuts, "Cuts",
                            "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart");

    if (Use["CutsD"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kCuts, "CutsD",
                            "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate");

    if (Use["CutsPCA"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kCuts, "CutsPCA",
                            "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA");

    if (Use["CutsGA"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kCuts, "CutsGA",
                            "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95");

    if (Use["CutsSA"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kCuts, "CutsSA",
                            "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale");

    // Likelihood ("naive Bayes estimator")
    if (Use["Likelihood"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kLikelihood, "Likelihood",
                            "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50");

    // Decorrelated likelihood
    if (Use["LikelihoodD"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kLikelihood, "LikelihoodD",
                            "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate");

    // PCA-transformed likelihood
    if (Use["LikelihoodPCA"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kLikelihood, "LikelihoodPCA",
                            "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA");

    // Use a kernel density estimator to approximate the PDFs
    if (Use["LikelihoodKDE"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kLikelihood, "LikelihoodKDE",
                            "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50");

    // Use a variable-dependent mix of splines and kernel density estimator
    if (Use["LikelihoodMIX"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kLikelihood, "LikelihoodMIX",
                            "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50");

    // Test the multi-dimensional probability density estimator
    // here are the options strings for the MinMax and RMS methods, respectively:
    //
    //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
    //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
    if (Use["PDERS"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kPDERS, "PDERS",
                            "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600");

    if (Use["PDERSD"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kPDERS, "PDERSD",
                            "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate");

    if (Use["PDERSPCA"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kPDERS, "PDERSPCA",
                            "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA");

    // Multi-dimensional likelihood estimator using self-adapting phase-space binning
    if (Use["PDEFoam"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kPDEFoam, "PDEFoam",
                            "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T");

    if (Use["PDEFoamBoost"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kPDEFoam, "PDEFoamBoost",
                            "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T");

    // K-Nearest Neighbour classifier (KNN)
    if (Use["KNN"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kKNN, "KNN",
                            "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim");

    // H-Matrix (chi2-squared) method
    if (Use["HMatrix"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kHMatrix, "HMatrix", "!H:!V:VarTransform=None");

    // Linear discriminant (same as Fisher discriminant)
    if (Use["LD"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10");

    // Fisher discriminant (same as LD)
    if (Use["Fisher"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10");

    // Fisher with Gauss-transformed input variables
    if (Use["FisherG"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss");

    // Composite classifier: ensemble (tree) of boosted Fisher classifiers
    if (Use["BoostedFisher"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kFisher, "BoostedFisher",
                            "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring");

    // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
    if (Use["FDA_MC"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kFDA, "FDA_MC",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1");

    if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
        factory->BookMethod(dataloader.get(), TMVA::Types::kFDA, "FDA_GA",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=100:Cycles=2:Steps=5:Trim=True:SaveBestGen=1");

    if (Use["FDA_SA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
        factory->BookMethod(dataloader.get(), TMVA::Types::kFDA, "FDA_SA",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale");

    if (Use["FDA_MT"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kFDA, "FDA_MT",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch");

    if (Use["FDA_GAMT"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kFDA, "FDA_GAMT",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim");

    if (Use["FDA_MCMT"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kFDA, "FDA_MCMT",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20");

    // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
    if (Use["MLP"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator");

    if (Use["MLPBFGS"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator");

    if (Use["MLPBNN"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=60:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator"); // BFGS training with bayesian regulators

    // Multi-architecture DNN implementation.
    if (Use["DNN_CPU"] or Use["DNN_GPU"])
    {
        // General layout.
        TString layoutString("Layout=TANH|128,TANH|128,TANH|128,LINEAR");

        // Define Training strategy. One could define multiple stratgey string separated by the "|" delimiter

        TString trainingStrategyString = ("TrainingStrategy=LearningRate=1e-2,Momentum=0.9,"
                                          "ConvergenceSteps=20,BatchSize=100,TestRepetitions=1,"
                                          "WeightDecay=1e-4,Regularization=None,"
                                          "DropConfig=0.0+0.5+0.5+0.5");

        // General Options.
        TString dnnOptions("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
                           "WeightInitialization=XAVIERUNIFORM");
        dnnOptions.Append(":");
        dnnOptions.Append(layoutString);
        dnnOptions.Append(":");
        dnnOptions.Append(trainingStrategyString);

        // Cuda implementation.
        if (Use["DNN_GPU"])
        {
            TString gpuOptions = dnnOptions + ":Architecture=GPU";
            factory->BookMethod(dataloader.get(), TMVA::Types::kDL, "DNN_GPU", gpuOptions);
        }
        // Multi-core CPU implementation.
        if (Use["DNN_CPU"])
        {
            TString cpuOptions = dnnOptions + ":Architecture=CPU";
            factory->BookMethod(dataloader.get(), TMVA::Types::kDL, "DNN_CPU", cpuOptions);
        }
    }

    // CF(Clermont-Ferrand)ANN
    if (Use["CFMlpANN"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N"); // n_cycles:#nodes:#nodes:...

    // Tmlp(Root)ANN
    if (Use["TMlpANN"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"); // n_cycles:#nodes:#nodes:...

    // Support Vector Machine
    if (Use["SVM"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm");

    // Boosted Decision Trees
    if (Use["BDTG"]) // Gradient Boost
        factory->BookMethod(dataloader.get(), TMVA::Types::kBDT, "BDTG",
                            "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2");

    if (Use["BDT"]) // Adaptive Boost
        factory->BookMethod(dataloader.get(), TMVA::Types::kBDT, "BDT",
                            "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");

    if (Use["BDTB"]) // Bagging
        factory->BookMethod(dataloader.get(), TMVA::Types::kBDT, "BDTB",
                            "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20");

    if (Use["BDTD"]) // Decorrelation + Adaptive Boost
        factory->BookMethod(dataloader.get(), TMVA::Types::kBDT, "BDTD",
                            "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate");

    if (Use["BDTF"]) // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
        factory->BookMethod(dataloader.get(), TMVA::Types::kBDT, "BDTF",
                            "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20");

    // RuleFit -- TMVA implementation of Friedman's method
    if (Use["RuleFit"])
        factory->BookMethod(dataloader.get(), TMVA::Types::kRuleFit, "RuleFit",
                            "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02");

    // For an example of the category classifier usage, see: TMVAClassificationCategory
    //
    // --------------------------------------------------------------------------------------------------
    //  Now you can optimize the setting (configuration) of the MVAs using the set of training events
    // STILL EXPERIMENTAL and only implemented for BDT's !
    //
    //     factory->OptimizeAllMethods("SigEffAtBkg0.01","Scan");
    //     factory->OptimizeAllMethods("ROCIntegral","FitGA");
    //
    // --------------------------------------------------------------------------------------------------
}