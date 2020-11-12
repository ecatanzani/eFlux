#include "mc_tuple.h"

#include <cmath>
#include <iostream>

#include "TFile.h"

mc_tuple::mc_tuple(const std::shared_ptr<TChain> chain)
{
    evtch = chain;
    set_core_address();
    set_simu_address();
}

void mc_tuple::set_simu_address()
{
    evtch->SetBranchAddress("simu_energy", &simu_energy);
    evtch->SetBranchAddress("simu_energy_w", &simu_energy_w);
    evtch->SetBranchAddress("energy_corr_w", &corr_energy_w);
    evtch->SetBranchAddress("simu_position", &simuPosition);
    evtch->SetBranchAddress("simu_momentum", &simuMomentum);
    evtch->SetBranchAddress("simuSlopeX", &simuSlopeX);
    evtch->SetBranchAddress("simuSlopeY", &simuSlopeY);
    evtch->SetBranchAddress("simuInterceptX", &simuInterceptX);
    evtch->SetBranchAddress("simuInterceptY", &simuInterceptY);
    evtch->SetBranchAddress("evtfilter_geometric_before_trigger", &evtfilter_geometric_before_trigger);
    evtch->SetBranchAddress("evtfilter_trigger_check", &evtfilter_trigger_check);
}

void mc_tuple::InitHistos(const std::vector<float> binning)
{
    _histos = std::make_shared<mc_histos>(binning);
    weight_shift = pow(binning[0], 2);
}

void mc_tuple::GetEntry(int idx)
{
    evtch->GetEntry(idx);
    // Set correct energy weight
    simu_energy_w *= weight_shift;
    fill_histos();
}

void mc_tuple::fill_histos()
{
    const double _GeV = 0.001;
    if (!evtfilter_out_energy_range)
    {
        if (evtfilter_geometric_before_trigger)
        {
            _histos->FillGeoFactor(simu_energy * _GeV, simu_energy_w);
            if (evtfilter_evt_triggered)
            {
                _histos->FillTrigger(simu_energy * _GeV, simu_energy_w);
                if (evtfilter_correct_bgo_reco)
                {
                    auto reco_bidx = _histos->GetEnergyBin(energy_corr * _GeV);
                    if (reco_bidx != -999)
                    {
                        auto BGOrec_costheta = compute_bgoreco_costheta();
                        _histos->FillIncoming(simu_energy * _GeV, simu_energy_w);
                        _histos->FillSimu(
                            simu_energy,
                            simu_energy_w,
                            energy,
                            energy_corr,
                            simuSlopeX,
                            simuSlopeY,
                            simuInterceptX,
                            simuInterceptY,
                            BGOrec_slopeX,
                            BGOrec_slopeY,
                            BGOrec_interceptX,
                            BGOrec_interceptY);
                        _histos->FillBGO(
                            energy,
                            energy_corr,
                            simu_energy_w,
                            reco_bidx,
                            fracLayer,
                            eLayer,
                            rmsLayer,
                            energy_1R_radius,
                            energy_2R_radius,
                            energy_3R_radius,
                            energy_5R_radius,
                            sumRms,
                            fracLast,
                            fracLast_13,
                            lastBGOLayer,
                            nBGOentries,
                            BGOrec_slopeX,
                            BGOrec_slopeY,
                            BGOrec_interceptX,
                            BGOrec_interceptY,
                            BGOrec_costheta);
                        _histos->FillSumRmsCosine(
                            energy_corr * _GeV,
                            simu_energy_w,
                            reco_bidx,
                            sumRms,
                            BGOrec_costheta);
                        _histos->FillSumRmsFLast(
                            energy_corr * _GeV,
                            simu_energy_w,
                            reco_bidx,
                            sumRms,
                            fracLast,
                            fracLast_13);
                        _histos->FillNUD(
                            simu_energy_w,
                            nud_adc,
                            nud_total_adc,
                            nud_max_adc,
                            nud_max_channel_id);
                        /*    
                        // Shower profile fit
                        auto fit_res = fit_shower_profile();
                        _histos->FillBGOShowerFit(
                            fit_res,
                            energy_corr,
                            simu_energy_w);
                        */
                        if (evtfilter_good_event)
                        {
                            if (evtfilter_geometric)
                            {
                                _histos->FillGeometric(simu_energy * _GeV, simu_energy_w);
                                if (evtfilter_BGO_fiducial_maxElayer_cut)
                                    _histos->FillGeometricMaxLayer(simu_energy * _GeV, simu_energy_w);
                                if (evtfilter_BGO_fiducial_maxBarLayer_cut)
                                    _histos->FillGeometricMaxBar(simu_energy * _GeV, simu_energy_w);
                                if (evtfilter_BGO_fiducial_BGOTrackContainment_cut)
                                    _histos->FillGeometricBGOTrack(simu_energy * _GeV, simu_energy_w);
                                if (evtfilter_BGO_fiducial)
                                    _histos->FillGeometricBGOFiducial(simu_energy * _GeV, simu_energy_w);
                                if (evtfilter_all_cut)
                                    _histos->FillGeometricAll(simu_energy * _GeV, simu_energy_w);
                            }

                            if (evtfilter_BGO_fiducial_maxElayer_cut)
                                _histos->FillMaxLayer(simu_energy * _GeV, simu_energy_w);
                            if (evtfilter_BGO_fiducial_maxBarLayer_cut)
                                _histos->FillMaxBar(simu_energy * _GeV, simu_energy_w);
                            if (evtfilter_BGO_fiducial_BGOTrackContainment_cut)
                                _histos->FillBGOTrack(simu_energy * _GeV, simu_energy_w);

                            if (evtfilter_BGO_fiducial)
                            {
                                _histos->FillBGOFiducial(simu_energy * _GeV, simu_energy_w);
                                if (evtfilter_nBarLayer13_cut)
                                    _histos->FillBGOFiducialBarL13(simu_energy * _GeV, simu_energy_w);
                                if (evtfilter_maxRms_cut)
                                    _histos->FillBGOFiducialMaxRms(simu_energy * _GeV, simu_energy_w);
                                if (evtfilter_track_selection_cut)
                                    _histos->FillBGOTrack(simu_energy * _GeV, simu_energy_w);
                                if (evtfilter_psd_stk_match_cut)
                                    _histos->FillBGOFiducialPsdStk(simu_energy * _GeV, simu_energy_w);
                                if (evtfilter_psd_charge_cut)
                                    _histos->FillBGOFiducialPsdCharge(simu_energy * _GeV, simu_energy_w);
                                if (evtfilter_stk_charge_cut)
                                    _histos->FillBGOFiducialStkCharge(simu_energy * _GeV, simu_energy_w);
                                if (evtfilter_all_cut)
                                    _histos->FillBGOFiducialAll(simu_energy * _GeV, simu_energy_w);
                            }

                            if (evtfilter_psd_charge_measurement)
                                _histos->FillPsdCharge(simu_energy_w, PSD_chargeX, PSD_chargeY, PSD_charge);
                            if (evtfilter_psd_charge_cut)
                                _histos->FillPsdCharge(simu_energy_w, PSD_chargeX, PSD_chargeY, PSD_charge, true);
                            if (evtfilter_stk_charge_measurement)
                            {
                                _histos->FillStkCharge(simu_energy_w, STK_chargeX, STK_chargeY, STK_charge);
                                _histos->FillStkCosine(simu_energy_w, STK_bestTrack_costheta, BGOrec_costheta);
                                _histos->FillBGOCosine(simu_energy_w, BGOrec_costheta);
                                // Shower profile fit
                                auto fit_res = fit_shower_profile(STK_bestTrack_costheta);
                                _histos->FillBGOShowerFit(
                                    fit_res,
                                    energy_corr,
                                    simu_energy_w);
                            }
                            if (evtfilter_stk_charge_cut)
                                _histos->FillStkCharge(simu_energy_w, STK_chargeX, STK_chargeY, STK_charge, true);

                            if (evtfilter_all_cut)
                            {
                                _histos->FillAllCut(simu_energy * _GeV, simu_energy_w);
                                _histos->FillClassifier(energy_corr * _GeV, simu_energy_w, xtrl);
                            }
                        }
                    }
                }
            }
            else
                _histos->FillIncoming(simu_energy * _GeV, simu_energy_w);
        }
    }
}

void mc_tuple::WriteHistos(const std::string output_path)
{
    TFile *outfile = TFile::Open(output_path.c_str(), "RECREATE");
    if (outfile->IsZombie())
    {
        std::cerr << "\n\nError writing TFile [" << output_path << "]\n\n";
        exit(100);
    }
    _histos->WriteCore(outfile);
    _histos->Write(outfile);
    outfile->Close();
}