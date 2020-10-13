#include "environment.h"
#include "read_parallel.h"
#include "aggregate_events.h"
#include "binning.h"

#include "TTreeReader.h"
#include "TTreeReaderValue.h"

void read_parallel(
    const std::string inputPath,
    TFile &outFile,
    const std::string wd,
    const bool verbose)
{
    auto dmpch = aggregateDataEventsTChain(
        inputPath,
        verbose);
    const auto nevents = dmpch->GetEntries();
    const auto logEBins = read_binning_from_config(wd);

    // Linking TChain
    // Create the TTreeReader object
    TTreeReader reader(dmpch.get());
    // Create the TTreeReaderValue objects

    // Trigger
    TTreeReaderValue<bool> unbiased_trigger(reader, "unbiased_trigger");
    TTreeReaderValue<bool> mip1_trigger(reader, "mip1_trigger");
    TTreeReaderValue<bool> mip2_trigger(reader, "mip2_trigger");
    TTreeReaderValue<bool> HET_trigger(reader, "HET_trigger");
    TTreeReaderValue<bool> LET_trigger(reader, "LET_trigger");
    TTreeReaderValue<bool> MIP_trigger(reader, "MIP_trigger");
    TTreeReaderValue<bool> general_trigger(reader, "general_trigger");
    // Time
    TTreeReaderValue<unsigned int> second(reader, "second");
    TTreeReaderValue<unsigned int> msecond(reader, "msecond");
    // STK
    TTreeReaderValue<unsigned int> STK_bestTrack_npoints(reader, "STK_bestTrack_npoints");
    TTreeReaderValue<unsigned int> STK_bestTrack_nholesX(reader, "STK_bestTrack_nholesX");
    TTreeReaderValue<unsigned int> STK_bestTrack_nholesY(reader, "STK_bestTrack_nholesY");
    TTreeReaderValue<double> STK_bestTrack_slopeX(reader, "STK_bestTrack_slopeX");
    TTreeReaderValue<double> STK_bestTrack_slopeY(reader, "STK_bestTrack_slopeY");
    TTreeReaderValue<double> STK_bestTrack_interceptX(reader, "STK_bestTrack_interceptX");
    TTreeReaderValue<double> STK_bestTrack_interceptY(reader, "STK_bestTrack_interceptY");
    TTreeReaderValue<double> STK_bestTrack_costheta(reader, "STK_bestTrack_costheta");
    TTreeReaderValue<double> STK_bestTrack_phi(reader, "STK_bestTrack_phi");
    TTreeReaderValue<double> STK_bestTrack_extr_BGO_topX(reader, "STK_bestTrack_extr_BGO_topX");
    TTreeReaderValue<double> STK_bestTrack_extr_BGO_topY(reader, "STK_bestTrack_extr_BGO_topY");
    TTreeReaderValue<double> STK_bestTrack_STK_BGO_topX_distance(reader, "STK_bestTrack_STK_BGO_topX_distance");
    TTreeReaderValue<double> STK_bestTrack_STK_BGO_topY_distance(reader, "STK_bestTrack_STK_BGO_topY_distance");
    TTreeReaderValue<double> STK_bestTrack_angular_distance_STK_BGO(reader, "STK_bestTrack_angular_distance_STK_BGO");
    TTreeReaderValue<double> STK_chargeX(reader, "STK_chargeX");
    TTreeReaderValue<double> STK_chargeY(reader, "STK_chargeY");
    TTreeReaderValue<double> STK_charge(reader, "STK_charge");
    // BGO
    TTreeReaderValue<double> energy(reader, "energy");
    TTreeReaderValue<double> energy_corr(reader, "energy_corr");
    TTreeReaderValue<double> BGOrec_slopeX(reader, "BGOrec_slopeX");
    TTreeReaderValue<double> BGOrec_slopeY(reader, "BGOrec_slopeY");
    TTreeReaderValue<double> BGOrec_interceptX(reader, "BGOrec_interceptX");
    TTreeReaderValue<double> BGOrec_interceptY(reader, "BGOrec_interceptY");
    TTreeReaderValue<double> sumRms(reader, "sumRms");
    TTreeReaderValue<std::vector<double>> fracLayer(reader, "fracLayer");
    TTreeReaderValue<double> fracLast(reader, "fracLast");
    TTreeReaderValue<double> fracLast_13(reader, "fracLast_13");
    TTreeReaderValue<unsigned int> lastBGOLayer(reader, "lastBGOLayer");
    TTreeReaderValue<unsigned int> nBGOentries(reader, "nBGOentries");
    // PSD
    TTreeReaderValue<double> PSD_chargeX(reader, "PSD_chargeX");
    TTreeReaderValue<double> PSD_chargeY(reader, "PSD_chargeY");
    TTreeReaderValue<double> PSD_charge(reader, "PSD_charge");
    // Classifier
    TTreeReaderValue<double> xtr(reader, "xtr");
    TTreeReaderValue<double> xtrl(reader, "xtrl");
    // Attitude
    TTreeReaderValue<double> glat(reader, "glat");
    TTreeReaderValue<double> glon(reader, "glon");
    // Cuts
    TTreeReaderValue<unsigned int> nActiveCuts(reader, "nActiveCuts");
    TTreeReaderValue<bool> evtfilter_geometric(reader, "evtfilter_geometric");
    TTreeReaderValue<bool> evtfilter_BGO_fiducial(reader, "evtfilter_BGO_fiducial");
    TTreeReaderValue<bool> evtfilter_all_cut(reader, "evtfilter_all_cut");
    TTreeReaderValue<bool> evtfilter_all_cut_no_xtrl(reader, "evtfilter_all_cut_no_xtrl");
    TTreeReaderValue<bool> evtfilter_BGO_fiducial_maxElayer_cut(reader, "evtfilter_BGO_fiducial_maxElayer_cut");
    TTreeReaderValue<bool> evtfilter_BGO_fiducial_maxBarLayer_cut(reader, "evtfilter_BGO_fiducial_maxBarLayer_cut");
    TTreeReaderValue<bool> evtfilter_BGO_fiducial_BGOTrackContainment_cut(reader, "evtfilter_BGO_fiducial_BGOTrackContainment_cut");
    TTreeReaderValue<bool> evtfilter_nBarLayer13_cut(reader, "evtfilter_nBarLayer13_cut");
    TTreeReaderValue<bool> evtfilter_maxRms_cut(reader, "evtfilter_maxRms_cut");
    TTreeReaderValue<bool> evtfilter_track_selection_cut(reader, "evtfilter_track_selection_cut");
    TTreeReaderValue<bool> evtfilter_psd_stk_match_cut(reader, "evtfilter_psd_stk_match_cut");
    TTreeReaderValue<bool> evtfilter_psd_charge_cut(reader, "evtfilter_psd_charge_cut");
    TTreeReaderValue<bool> evtfilter_stk_charge_cut(reader, "evtfilter_stk_charge_cut");
    TTreeReaderValue<bool> evtfilter_xtrl_cut(reader, "evtfilter_xtrl_cut");
    TTreeReaderValue<bool> evtfilter_psd_charge_measurement(reader, "evtfilter_psd_charge_measurement");
    TTreeReaderValue<bool> evtfilter_stk_charge_measurement(reader, "evtfilter_stk_charge_measurement");

    for (auto i : ROOT::TSeqUL(nevents))
    {
        dmpch->GetEntry(i);
    }
}