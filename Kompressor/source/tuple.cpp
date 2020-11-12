#include "tuple.h"

#include "TF1.h"
#include "TROOT.h"
#include "TMath.h"
#include "TGraph.h"
#include "TMinuit.h"

void ntuple::set_core_address()
{
    // Trigger
    evtch->SetBranchAddress("mip1_trigger", &mip1_trigger);
    evtch->SetBranchAddress("mip2_trigger", &mip2_trigger);
    evtch->SetBranchAddress("HET_trigger", &HET_trigger);
    evtch->SetBranchAddress("LET_trigger", &LET_trigger);
    evtch->SetBranchAddress("MIP_trigger", &MIP_trigger);
    evtch->SetBranchAddress("general_trigger", &general_trigger);
    // STK
    evtch->SetBranchAddress("STK_bestTrack_npoints", &STK_bestTrack_npoints);
    evtch->SetBranchAddress("STK_bestTrack_nholesX", &STK_bestTrack_nholesX);
    evtch->SetBranchAddress("STK_bestTrack_nholesY", &STK_bestTrack_nholesY);
    evtch->SetBranchAddress("STK_bestTrack_slopeX", &STK_bestTrack_slopeX);
    evtch->SetBranchAddress("STK_bestTrack_slopeY", &STK_bestTrack_slopeY);
    evtch->SetBranchAddress("STK_bestTrack_interceptX", &STK_bestTrack_interceptX);
    evtch->SetBranchAddress("STK_bestTrack_interceptY", &STK_bestTrack_interceptY);
    evtch->SetBranchAddress("STK_bestTrack_costheta", &STK_bestTrack_costheta);
    evtch->SetBranchAddress("STK_bestTrack_phi", &STK_bestTrack_phi);
    evtch->SetBranchAddress("STK_bestTrack_extr_BGO_topX", &STK_bestTrack_extr_BGO_topX);
    evtch->SetBranchAddress("STK_bestTrack_extr_BGO_topY", &STK_bestTrack_extr_BGO_topY);
    evtch->SetBranchAddress("STK_bestTrack_STK_BGO_topX_distance", &STK_bestTrack_STK_BGO_topX_distance);
    evtch->SetBranchAddress("STK_bestTrack_STK_BGO_topY_distance", &STK_bestTrack_STK_BGO_topY_distance);
    evtch->SetBranchAddress("STK_bestTrack_angular_distance_STK_BGO", &STK_bestTrack_angular_distance_STK_BGO);
    evtch->SetBranchAddress("STK_chargeX", &STK_chargeX);
    evtch->SetBranchAddress("STK_chargeY", &STK_chargeY);
    evtch->SetBranchAddress("STK_charge", &STK_charge);
    // BGO
    evtch->SetBranchAddress("energy", &energy);
    evtch->SetBranchAddress("energy_corr", &energy_corr);
    evtch->SetBranchAddress("eLayer", &eLayer);
    evtch->SetBranchAddress("BGOrec_slopeX", &BGOrec_slopeX);
    evtch->SetBranchAddress("BGOrec_slopeY", &BGOrec_slopeY);
    evtch->SetBranchAddress("BGOrec_interceptX", &BGOrec_interceptX);
    evtch->SetBranchAddress("BGOrec_interceptY", &BGOrec_interceptY);
    evtch->SetBranchAddress("sumRms", &sumRms);
    evtch->SetBranchAddress("rmsLayer", &rmsLayer);
    evtch->SetBranchAddress("fracLayer", &fracLayer);
    evtch->SetBranchAddress("fracLast", &fracLast);
    evtch->SetBranchAddress("fracLast_13", &fracLast_13);
    evtch->SetBranchAddress("lastBGOLayer", &lastBGOLayer);
    evtch->SetBranchAddress("nBGOentries", &nBGOentries);
    evtch->SetBranchAddress("energy_1R_radius", &energy_1R_radius);
    evtch->SetBranchAddress("energy_2R_radius", &energy_2R_radius);
    evtch->SetBranchAddress("energy_3R_radius", &energy_3R_radius);
    evtch->SetBranchAddress("energy_5R_radius", &energy_5R_radius);
    // PSD
    evtch->SetBranchAddress("PSD_chargeX", &PSD_chargeX);
    evtch->SetBranchAddress("PSD_chargeY", &PSD_chargeY);
    evtch->SetBranchAddress("PSD_charge", &PSD_charge);
    // NUD
    evtch->SetBranchAddress("NUD_ADC", &nud_adc);
    evtch->SetBranchAddress("NUD_total_ADC", &nud_total_adc);
    evtch->SetBranchAddress("NUD_max_ADC", &nud_max_adc);
    evtch->SetBranchAddress("NUD_max_channel_ID", &nud_max_channel_id);
    // Classifiers
    evtch->SetBranchAddress("xtr", &xtr);
    evtch->SetBranchAddress("xtrl", &xtrl);
    // Filters
    evtch->SetBranchAddress("evtfilter_out_energy_range", &evtfilter_out_energy_range);
    evtch->SetBranchAddress("evtfilter_evt_triggered", &evtfilter_evt_triggered);
    evtch->SetBranchAddress("evtfilter_correct_bgo_reco", &evtfilter_correct_bgo_reco);
    evtch->SetBranchAddress("evtfilter_good_event", &evtfilter_good_event);
    evtch->SetBranchAddress("evtfilter_geometric", &evtfilter_geometric);
    evtch->SetBranchAddress("evtfilter_BGO_fiducial", &evtfilter_BGO_fiducial);
    evtch->SetBranchAddress("evtfilter_BGO_fiducial_maxElayer_cut", &evtfilter_BGO_fiducial_maxElayer_cut);
    evtch->SetBranchAddress("evtfilter_BGO_fiducial_maxBarLayer_cut", &evtfilter_BGO_fiducial_maxBarLayer_cut);
    evtch->SetBranchAddress("evtfilter_BGO_fiducial_BGOTrackContainment_cut", &evtfilter_BGO_fiducial_BGOTrackContainment_cut);
    evtch->SetBranchAddress("evtfilter_nBarLayer13_cut", &evtfilter_nBarLayer13_cut);
    evtch->SetBranchAddress("evtfilter_maxRms_cut", &evtfilter_maxRms_cut);
    evtch->SetBranchAddress("evtfilter_track_selection_cut", &evtfilter_track_selection_cut);
    evtch->SetBranchAddress("evtfilter_psd_stk_match_cut", &evtfilter_psd_stk_match_cut);
    evtch->SetBranchAddress("evtfilter_psd_charge_cut", &evtfilter_psd_charge_cut);
    evtch->SetBranchAddress("evtfilter_stk_charge_cut", &evtfilter_stk_charge_cut);
    evtch->SetBranchAddress("evtfilter_psd_charge_measurement", &evtfilter_psd_charge_measurement);
    evtch->SetBranchAddress("evtfilter_stk_charge_measurement", &evtfilter_stk_charge_measurement);
    evtch->SetBranchAddress("evtfilter_xtrl_tight_cut", &evtfilter_xtrl_tight_cut);
    evtch->SetBranchAddress("evtfilter_xtrl_loose_cut", &evtfilter_xtrl_loose_cut);
    evtch->SetBranchAddress("evtfilter_all_cut", &evtfilter_all_cut);
    evtch->SetBranchAddress("cut_nBarLayer13", &cut_nBarLayer13);
    evtch->SetBranchAddress("cut_maxRms", &cut_maxRms);
    evtch->SetBranchAddress("cut_track_selection", &cut_track_selection);
    evtch->SetBranchAddress("cut_psd_stk_match", &cut_psd_stk_match);
    evtch->SetBranchAddress("cut_psd_charge", &cut_psd_charge);
    evtch->SetBranchAddress("cut_stk_charge", &cut_stk_charge);
    evtch->SetBranchAddress("nActiveCuts", &nActiveCuts);
}

double function(double *x, double *par)
{
    double xx = x[0];
    double func = par[0] * par[1] * (TMath::Power(par[1] * xx, par[2] - 1) * TMath::Exp(-par[1] * xx)) / TMath::Gamma(par[2]);
    return func;
}

const std::vector<double> ntuple::fit_shower_profile(const double costheta)
{
    const double bgoX0 = 11.2;                                                         // BGO X0 in mm
    const double bgoEc = 10.50;                                                        // BGO critical energy for electrons in MeV
    const double b_shower_par = 0.5;                                                   // b parameter electromagnetic shower in BGO
    const double a_shower_par = 1 + b_shower_par * (TMath::Log(energy / bgoEc) - 0.5); // a parameter electromagnetic shower in BGO

    std::vector<double> t_bgo(DAMPE_bgo_nLayers, -999);
    for (int idx = 0; idx < DAMPE_bgo_nLayers; ++idx)
        t_bgo[idx] = (BGO_bar_lateral * (idx + 1)) / (bgoX0 * costheta);

#if _DEBUG
    std::cout << "\nBGO layer displacement: " << BGO_layer_displacement;
    std::cout << "\nSTK best track costheta: " << STK_bestTrack_costheta;
    //std::cout << "\nBGO costheta: " << bgo_costheta;
    std::cout << "\nShower profile scheme:\n";
    for (int idx = 0; idx < DAMPE_bgo_nLayers; ++idx)
        std::cout << "\nLayer " << idx << ": Energy: " << eLayer->at(idx) << " t: " << t_bgo.at(idx);
    std::cout << std::endl;
    std::cout << "\nDAMPE raw energy: " << energy;
    std::cout << "\nDAMPE corrected energy: " << energy_corr << std::endl;
#endif

    TGraph gr_profile(DAMPE_bgo_nLayers, &t_bgo.at(0), &eLayer->at(0));
    TF1 fitfunc("fitfunc", function, t_bgo[0], t_bgo[DAMPE_bgo_nLayers-1], 3);
    fitfunc.SetParameter(0, energy);
    //fitfunc.SetParLimits(0, 0, 1e+10);
    //fitfunc.FixParameter(1, b_shower_par);
    //fitfunc.FixParameter(2, a_shower_par);
    fitfunc.SetParameter(1, b_shower_par);
    fitfunc.SetParameter(2, a_shower_par);
    fitfunc.SetParNames("bgo_energy", "b", "a");
#if _DEBUG    
    gr_profile.Fit(&fitfunc,"IR");
#else
    gr_profile.Fit(&fitfunc,"qIR");
#endif
    std::vector<double> fit_res(3, -999);
    for (int idx=0; idx<3; ++idx)
        fit_res[idx] = fitfunc.GetParameter(idx);

    return fit_res;
}

const double ntuple::compute_bgoreco_costheta()
{
    /*
    const double bgo_deltaX = fabs(BGOrec_slopeX * BGO_TopZ + BGOrec_interceptX - BGOrec_slopeX * BGO_BottomZ + BGOrec_interceptX);
    const double bgo_deltaY = fabs(BGOrec_slopeY * BGO_TopZ + BGOrec_interceptY - BGOrec_slopeY * BGO_BottomZ + BGOrec_interceptY);
    const double bgo_costheta = TMath::Cos(TMath::ATan2(bgo_deltaY, bgo_deltaX));
    */
    const double rise = BGO_TopZ - BGO_BottomZ;
    const double bgo_deltaX = fabs(BGOrec_slopeX * BGO_TopZ + BGOrec_interceptX - BGOrec_slopeX * BGO_BottomZ + BGOrec_interceptX);
    const double bgo_deltaY = fabs(BGOrec_slopeY * BGO_TopZ + BGOrec_interceptY - BGOrec_slopeY * BGO_BottomZ + BGOrec_interceptY);
    const double run = TMath::Sqrt(TMath::Power(bgo_deltaX, 2) + TMath::Power(bgo_deltaY, 2));
    const double bgo_costheta = TMath::Cos(fabs(rise/run));
    return bgo_costheta;
}