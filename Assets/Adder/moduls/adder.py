from ROOT import TFile, TH1D, TH2D, TGraph
import os


def compute_final_histos_mc(condor_dir_list, opts):
    
    # Cut histos
    h_geo_factor = TH1D()
    h_incoming = TH1D()
    h_trigger = TH1D()
    h_geometric_cut = TH1D()
    h_maxElayer_cut = TH1D()
    h_maxBarLayer_cut = TH1D()
    h_BGOTrackContainment_cut = TH1D()
    h_BGO_fiducial_cut = TH1D()
    h_all_cut = TH1D()

    h_geo_factor_nw = TH1D()
    h_incoming_nw = TH1D()
    h_trigger_nw = TH1D()
    h_geometric_cut_nw = TH1D()
    h_maxElayer_cut_nw = TH1D()
    h_maxBarLayer_cut_nw = TH1D()
    h_BGOTrackContainment_cut_nw = TH1D()
    h_BGO_fiducial_cut_nw = TH1D()
    h_all_cut_nw = TH1D()

    # Cuts && Geometric Cut
    h_geometric_maxElayer_cut = TH1D()
    h_geometric_maxBarLayer_cut = TH1D()
    h_geometric_BGOTrackContainment_cut = TH1D()
    h_geometric_BGO_fiducial_cut = TH1D()
    h_geometric_all_cut = TH1D()

    h_geometric_maxElayer_cut_nw = TH1D()
    h_geometric_maxBarLayer_cut_nw = TH1D()
    h_geometric_BGOTrackContainment_cut_nw = TH1D()
    h_geometric_BGO_fiducial_cut_nw = TH1D()
    h_geometric_all_cut_nw = TH1D()

    # Cuts && BGO fiducial volume cut
    h_BGOfiducial_nBarLayer13_cut = TH1D()
    h_BGOfiducial_maxRms_cut = TH1D()
    h_BGOfiducial_track_selection_cut = TH1D()
    h_BGOfiducial_psd_stk_match_cut = TH1D()
    h_BGOfiducial_psd_charge_cut = TH1D()
    h_BGOfiducial_stk_charge_cut = TH1D()
    h_BGOfiducial_xtrl_cut = TH1D()
    h_BGOfiducial_all_cut = TH1D()

    h_BGOfiducial_nBarLayer13_cut_nw = TH1D()
    h_BGOfiducial_maxRms_cut_nw = TH1D()
    h_BGOfiducial_track_selection_cut_nw = TH1D()
    h_BGOfiducial_psd_stk_match_cut_nw = TH1D()
    h_BGOfiducial_psd_charge_cut_nw = TH1D()
    h_BGOfiducial_stk_charge_cut_nw = TH1D()
    h_BGOfiducial_xtrl_cut_nw = TH1D()
    h_BGOfiducial_all_cut_nw = TH1D()

    # Analysis histos - simu and reco energy of incoming events
    h_BGOrec_energy = TH1D()
    h_simu_energy = TH1D()
    h_energy_diff = TH1D()
    h_energy_diff2D = TH2D()
    h_energy_unfold = TH2D()

    # Ratio of layer energy respect to total BGO energy
    h_layer_max_energy_ratio = TH1D()
    
    h_layer_energy_ratio = []
    for idx in range(0,14):
        h_layer_energy_ratio.append(TH1D())

    # Pre Geometric Cut
    # Top X and Y
    h_preGeo_BGOrec_topX_vs_realX = TH1D()
    h_preGeo_BGOrec_topY_vs_realY = TH1D()
    # Slope X and Y
    h_preGeo_real_slopeX = TH1D()
    h_preGeo_real_slopeY = TH1D()
    h_preGeo_BGOrec_slopeX = TH1D()
    h_preGeo_BGOrec_slopeY = TH1D()
    # Intercept X and Y
    h_preGeo_real_interceptX = TH1D()
    h_preGeo_real_interceptY = TH1D()
    h_preGeo_BGOrec_interceptX = TH1D()
    h_preGeo_BGOrec_interceptY = TH1D()
    # Top Maps
    h_preGeo_real_topMap = TH2D()
    h_preGeo_BGOreco_topMap = TH2D()
    # Bottom Maps
    h_preGeo_real_bottomMap = TH2D()
    h_preGeo_BGOreco_bottomMap = TH2D()
    # Map of events outside the "real" first BGO layer
    h_noBGOenergy_real_topMap = TH2D()

    # After Geometric Cut
    # Top X and Y
    h_geo_BGOrec_topX_vs_realX = TH1D()
    h_geo_BGOrec_topY_vs_realY = TH1D()
    # Slope X and Y
    h_geo_real_slopeX = TH1D()
    h_geo_real_slopeY = TH1D()
    h_geo_BGOrec_slopeX = TH1D()
    h_geo_BGOrec_slopeY = TH1D()
    # Intercept X and Y
    h_geo_real_interceptX = TH1D()
    h_geo_real_interceptY = TH1D()
    h_geo_BGOrec_interceptX = TH1D()
    h_geo_BGOrec_interceptY = TH1D()
    # Top Maps
    h_geo_real_topMap = TH2D()
    h_geo_BGOreco_topMap = TH2D()
    # Bottom Maps
    h_geo_real_bottomMap = TH2D()
    h_geo_BGOreco_bottomMap = TH2D()

    # sumRms - cosine correlation
    h_sumRms_cosine = []
    sumRms_cosine_20_100 = TH2D()
    sumRms_cosine_100_250 = TH2D()
    sumRms_cosine_250_500 = TH2D()
    sumRms_cosine_500_1000 = TH2D()
    sumRms_cosine_1000_3000 = TH2D()
    sumRms_cosine_3000_10000 = TH2D()

    # XTRL histos
    h_xtrl_energy_int = TH1D()
    h_xtrl = TH2D()
    e_discrimination_last = TH2D()
    e_discrimination_last_20_100 = TH2D()
    e_discrimination_last_100_250 = TH2D()
    e_discrimination_last_250_500 = TH2D()
    e_discrimination_last_500_1000 = TH2D()
    e_discrimination_last_1000_3000 = TH2D()
    e_discrimination_last_3000_10000 = TH2D()
    e_discrimination = TH2D()
    e_discrimination_20_100 = TH2D()
    e_discrimination_100_250 = TH2D()
    e_discrimination_250_500 = TH2D()
    e_discrimination_500_1000 = TH2D()
    e_discrimination_1000_3000 = TH2D()
    e_discrimination_3000_10000 = TH2D()
    h_bin_xtrl = []

    # PSD charge histos
    h_psd_chargeX = TH1D()
    h_psd_chargeY = TH1D()
    h_psd_charge = TH1D()
    h_psd_charge2D = TH2D()

    h_psd_selected_chargeX = TH1D()
    h_psd_selected_chargeY = TH1D()
    h_psd_selected_charge = TH1D()
    h_psd_selected_charge2D = TH2D()

    # STK charge histos
    h_stk_chargeX = TH1D()
    h_stk_chargeY = TH1D()
    h_stk_charge = TH1D()
    h_stk_charge2D = TH2D()

    h_stk_selected_chargeX = TH1D()
    h_stk_selected_chargeY = TH1D()
    h_stk_selected_charge = TH1D()
    h_stk_selected_charge2D = TH2D()

    # Proton background
    h_background_under_xtrl_cut = TH1D()
    h_background_over_xtrl_cut = TH1D()

    for dIdx, tmp_dir in enumerate(condor_dir_list):
        tmp_dir += "/outFiles"
        tmp_dir_list = os.listdir(tmp_dir)
        for elm in tmp_dir_list:
            if elm.startswith("analysisOutFile_"):
                rFile_path = tmp_dir + "/" + elm

        # Open ROOT output file
        rFile = TFile.Open(rFile_path, "READ")
        if rFile.IsOpen():
            if opts.verbose:
                if dIdx == 0:
                    print('\nReading file {}: {}'.format((dIdx+1), rFile_path))
                else:
                    print('Reading file {}: {}'.format((dIdx+1), rFile_path))
        else:
            print('Error reading file {}: {}'.format((dIdx+1), rFile_path))
            sys.exit()
        
        # Check the keys
        if rFile.GetNkeys() == 0:
            continue

        # Reading histos
        h_geo_factor_tmp = rFile.Get("h_geo_factor")
        h_incoming_tmp = rFile.Get("h_incoming")
        h_trigger_tmp = rFile.Get("h_trigger")
        h_geometric_cut_tmp = rFile.Get("h_geometric_cut")
        h_maxElayer_cut_tmp = rFile.Get("h_maxElayer_cut")
        h_maxBarLayer_cut_tmp = rFile.Get("h_maxBarLayer_cut")
        h_BGOTrackContainment_cut_tmp = rFile.Get("h_BGOTrackContainment_cut")
        h_BGO_fiducial_cut_tmp = rFile.Get("h_BGO_fiducial_cut")
        h_all_cut_tmp = rFile.Get("h_all_cut")

        n_energy_bins = h_geo_factor_tmp.GetNbinsX()
        # Init bin_xtrl and sumRms
        for idx in range(n_energy_bins):
            h_bin_xtrl.append(TH1D())
            h_sumRms_cosine.append(TH2D())

        h_geo_factor_tmp_nw = rFile.Get("nw_histos/h_geo_factor_w")
        h_incoming_tmp_nw = rFile.Get("nw_histos/h_incoming_w")
        h_trigger_tmp_nw = rFile.Get("nw_histos/h_trigger_w")
        h_geometric_cut_tmp_nw = rFile.Get("nw_histos/h_geometric_cut_w")
        h_maxElayer_cut_tmp_nw = rFile.Get("nw_histos/h_maxElayer_cut_w")
        h_maxBarLayer_cut_tmp_nw = rFile.Get("nw_histos/h_maxBarLayer_cut_w")
        h_BGOTrackContainment_cut_tmp_nw = rFile.Get("nw_histos/h_BGOTrackContainment_cut_w")
        h_BGO_fiducial_cut_tmp_nw = rFile.Get("nw_histos/h_BGO_fiducial_cut_w")
        h_all_cut_tmp_nw = rFile.Get("nw_histos/h_all_cut_w")

        h_geometric_maxElayer_cut_tmp = rFile.Get("h_geometric_maxElayer_cut")
        h_geometric_maxBarLayer_cut_tmp = rFile.Get("h_geometric_maxBarLayer_cut")
        h_geometric_BGOTrackContainment_cut_tmp = rFile.Get("h_geometric_BGOTrackContainment_cut")
        h_geometric_BGO_fiducial_cut_tmp = rFile.Get("h_geometric_BGO_fiducial_cut")
        h_geometric_all_cut_tmp = rFile.Get("h_geometric_all_cut")

        h_geometric_maxElayer_cut_tmp_nw = rFile.Get("nw_histos/h_geometric_maxElayer_cut_w")
        h_geometric_maxBarLayer_cut_tmp_nw = rFile.Get("nw_histos/h_geometric_maxBarLayer_cut_w")
        h_geometric_BGOTrackContainment_cut_tmp_nw = rFile.Get("nw_histos/h_geometric_BGOTrackContainment_cut_w")
        h_geometric_BGO_fiducial_cut_tmp_nw = rFile.Get("nw_histos/h_geometric_BGO_fiducial_cut_w")
        h_geometric_all_cut_tmp_nw = rFile.Get("nw_histos/h_geometric_all_cut_w")

        h_BGOfiducial_nBarLayer13_cut_tmp = rFile.Get("h_BGOfiducial_nBarLayer13_cut")
        h_BGOfiducial_maxRms_cut_tmp = rFile.Get("h_BGOfiducial_maxRms_cut")
        h_BGOfiducial_track_selection_cut_tmp = rFile.Get("h_BGOfiducial_track_selection_cut")
        h_BGOfiducial_psd_stk_match_cut_tmp = rFile.Get("h_BGOfiducial_psd_stk_match_cut")
        h_BGOfiducial_psd_charge_cut_tmp = rFile.Get("h_BGOfiducial_psd_charge_cut")
        h_BGOfiducial_stk_charge_cut_tmp = rFile.Get("h_BGOfiducial_stk_charge_cut")
        h_BGOfiducial_xtrl_cut_tmp = rFile.Get("h_BGOfiducial_xtrl_cut")
        h_BGOfiducial_all_cut_tmp = rFile.Get("h_BGOfiducial_all_cut")

        h_BGOfiducial_nBarLayer13_cut_tmp_nw = rFile.Get("nw_histos/h_BGOfiducial_nBarLayer13_cut_w")
        h_BGOfiducial_maxRms_cut_tmp_nw = rFile.Get("nw_histos/h_BGOfiducial_maxRms_cut_w")
        h_BGOfiducial_track_selection_cut_tmp_nw = rFile.Get("nw_histos/h_BGOfiducial_track_selection_cut_w")
        h_BGOfiducial_psd_stk_match_cut_tmp_nw = rFile.Get("nw_histos/h_BGOfiducial_psd_stk_match_cut_w")
        h_BGOfiducial_psd_charge_cut_tmp_nw = rFile.Get("nw_histos/h_BGOfiducial_psd_charge_cut_w")
        h_BGOfiducial_stk_charge_cut_tmp_nw = rFile.Get("nw_histos/h_BGOfiducial_stk_charge_cut_w")
        h_BGOfiducial_xtrl_cut_tmp_nw = rFile.Get("nw_histos/h_BGOfiducial_xtrl_cut_w")
        h_BGOfiducial_all_cut_tmp_nw = rFile.Get("nw_histos/h_BGOfiducial_all_cut_w")

        h_preGeo_BGOrec_topX_vs_realX_tmp = rFile.Get("Analysis_preGeoCut/h_preGeo_BGOrec_topX_vs_realX")
        h_preGeo_BGOrec_topY_vs_realY_tmp = rFile.Get("Analysis_preGeoCut/h_preGeo_BGOrec_topY_vs_realY")
        h_preGeo_real_slopeX_tmp = rFile.Get("Analysis_preGeoCut/h_preGeo_real_slopeX")
        h_preGeo_real_slopeY_tmp = rFile.Get("Analysis_preGeoCut/h_preGeo_real_slopeY")
        h_preGeo_BGOrec_slopeX_tmp = rFile.Get("Analysis_preGeoCut/h_preGeo_BGOrec_slopeX")
        h_preGeo_BGOrec_slopeY_tmp = rFile.Get("Analysis_preGeoCut/h_preGeo_BGOrec_slopeY")
        h_preGeo_real_interceptX_tmp = rFile.Get("Analysis_preGeoCut/h_preGeo_real_interceptX")
        h_preGeo_real_interceptY_tmp = rFile.Get("Analysis_preGeoCut/h_preGeo_real_interceptY")
        h_preGeo_BGOrec_interceptX_tmp = rFile.Get("Analysis_preGeoCut/h_preGeo_BGOrec_interceptX")
        h_preGeo_BGOrec_interceptY_tmp = rFile.Get("Analysis_preGeoCut/h_preGeo_BGOrec_interceptY")
        h_preGeo_real_topMap_tmp = rFile.Get("Analysis_preGeoCut/h_preGeo_real_topMap")
        h_preGeo_BGOreco_topMap_tmp = rFile.Get("Analysis_preGeoCut/h_preGeo_BGOreco_topMap")
        h_preGeo_real_bottomMap_tmp = rFile.Get("Analysis_preGeoCut/h_preGeo_real_bottomMap")
        h_preGeo_BGOreco_bottomMap_tmp = rFile.Get("Analysis_preGeoCut/h_preGeo_BGOreco_bottomMap")
        h_noBGOenergy_real_topMap_tmp = rFile.Get("Analysis_preGeoCut/h_noBGOenergy_real_topMap")

        h_geo_BGOrec_topX_vs_realX_tmp = rFile.Get("Analysis_GeoCut/h_geo_BGOrec_topX_vs_realX")
        h_geo_BGOrec_topY_vs_realY_tmp = rFile.Get("Analysis_GeoCut/h_geo_BGOrec_topY_vs_realY")
        h_geo_real_slopeX_tmp = rFile.Get("Analysis_GeoCut/h_geo_real_slopeX")
        h_geo_real_slopeY_tmp = rFile.Get("Analysis_GeoCut/h_geo_real_slopeY")
        h_geo_BGOrec_slopeX_tmp = rFile.Get("Analysis_GeoCut/h_geo_BGOrec_slopeX")
        h_geo_BGOrec_slopeY_tmp = rFile.Get("Analysis_GeoCut/h_geo_BGOrec_slopeY")
        h_geo_real_interceptX_tmp = rFile.Get("Analysis_GeoCut/h_geo_real_interceptX")
        h_geo_real_interceptY_tmp = rFile.Get("Analysis_GeoCut/h_geo_real_interceptY")
        h_geo_BGOrec_interceptX_tmp = rFile.Get("Analysis_GeoCut/h_geo_BGOrec_interceptX")
        h_geo_BGOrec_interceptY_tmp = rFile.Get("Analysis_GeoCut/h_geo_BGOrec_interceptY")
        h_geo_real_topMap_tmp = rFile.Get("Analysis_GeoCut/h_geo_real_topMap")
        h_geo_BGOreco_topMap_tmp = rFile.Get("Analysis_GeoCut/h_geo_BGOreco_topMap")
        h_geo_real_bottomMap_tmp = rFile.Get("Analysis_GeoCut/h_geo_real_bottomMap")
        h_geo_BGOreco_bottomMap_tmp = rFile.Get("Analysis_GeoCut/h_geo_BGOreco_bottomMap")
        
        h_BGOrec_energy_tmp = rFile.Get("BGO_Energy/h_BGOrec_energy")
        h_simu_energy_tmp = rFile.Get("BGO_Energy/h_simu_energy")
        h_energy_diff_tmp = rFile.Get("BGO_Energy/h_energy_diff")
        h_energy_diff2D_tmp = rFile.Get("BGO_Energy/h_energy_diff2D")
        h_energy_unfold_tmp = rFile.Get("BGO_Energy/h_energy_unfold")
        h_layer_max_energy_ratio_tmp = rFile.Get("BGO_Energy/h_layer_max_energy_ratio")
        h_layer_energy_ratio_tmp = []
        for idx in range(0,14):
            histoName = "BGO_Energy/h_layer_energy_ratio_"
            histoName += str(idx)
            h_layer_energy_ratio_tmp.append(rFile.Get(histoName))

        h_sumRms_cosine_tmp = []
        for idx in range(n_energy_bins):
            histoName = "BGO_Energy/sumRms_cosine_"
            histoName += str(idx)
            h_sumRms_cosine_tmp.append(rFile.Get(histoName))

        sumRms_cosine_20_100_tmp = rFile.Get("BGO_Energy/sumRms_cosine_20_100")
        sumRms_cosine_100_250_tmp = rFile.Get("BGO_Energy/sumRms_cosine_100_250")
        sumRms_cosine_250_500_tmp = rFile.Get("BGO_Energy/sumRms_cosine_250_500")
        sumRms_cosine_500_1000_tmp = rFile.Get("BGO_Energy/sumRms_cosine_500_1000")
        sumRms_cosine_1000_3000_tmp = rFile.Get("BGO_Energy/sumRms_cosine_1000_3000")
        sumRms_cosine_3000_10000_tmp = rFile.Get("BGO_Energy/sumRms_cosine_3000_10000")

        h_xtrl_energy_int_tmp = rFile.Get("xtrl/h_xtrl_energy_int")
        h_xtrl_tmp = rFile.Get("xtrl/h_xtrl")
        e_discrimination_last_tmp = rFile.Get("xtrl/e_discrimination_last")
        e_discrimination_last_20_100_tmp = rFile.Get("xtrl/e_discrimination_last_20_100")
        e_discrimination_last_100_250_tmp = rFile.Get("xtrl/e_discrimination_last_100_250")
        e_discrimination_last_250_500_tmp = rFile.Get("xtrl/e_discrimination_last_250_500")
        e_discrimination_last_500_1000_tmp = rFile.Get("xtrl/e_discrimination_last_500_1000")
        e_discrimination_last_1000_3000_tmp = rFile.Get("xtrl/e_discrimination_last_1000_3000")
        e_discrimination_last_3000_10000_tmp = rFile.Get("xtrl/e_discrimination_last_3000_10000")

        e_discrimination_tmp = rFile.Get("xtrl/e_discrimination")
        e_discrimination_20_100_tmp = rFile.Get("xtrl/e_discrimination_20_100")
        e_discrimination_100_250_tmp = rFile.Get("xtrl/e_discrimination_100_250")
        e_discrimination_250_500_tmp = rFile.Get("xtrl/e_discrimination_250_500")
        e_discrimination_500_1000_tmp = rFile.Get("xtrl/e_discrimination_500_1000")
        e_discrimination_1000_3000_tmp = rFile.Get("xtrl/e_discrimination_1000_3000")
        e_discrimination_3000_10000_tmp = rFile.Get("xtrl/e_discrimination_3000_10000")

        h_bin_xtrl_tmp = []
        for idx in range(n_energy_bins):
            histoName = "xtrl/h_xtrl_bin_"
            histoName += str(idx)
            h_bin_xtrl_tmp.append(rFile.Get(histoName))

        h_psd_chargeX_tmp = rFile.Get("PSDcharge/h_psd_chargeX")
        h_psd_chargeY_tmp = rFile.Get("PSDcharge/h_psd_chargeY")
        h_psd_charge_tmp = rFile.Get("PSDcharge/h_psd_charge")
        h_psd_charge2D_tmp = rFile.Get("PSDcharge/h_psd_charge2D")

        h_psd_selected_chargeX_tmp = rFile.Get("PSDcharge/h_psd_selected_chargeX")
        h_psd_selected_chargeY_tmp = rFile.Get("PSDcharge/h_psd_selected_chargeY")
        h_psd_selected_charge_tmp = rFile.Get("PSDcharge/h_psd_selected_charge")
        h_psd_selected_charge2D_tmp = rFile.Get("PSDcharge/h_psd_selected_charge2D")

        h_stk_chargeX_tmp = rFile.Get("STKcharge/h_stk_chargeX")
        h_stk_chargeY_tmp = rFile.Get("STKcharge/h_stk_chargeY")
        h_stk_charge_tmp = rFile.Get("STKcharge/h_stk_charge")
        h_stk_charge2D_tmp = rFile.Get("STKcharge/h_stk_charge2D")

        h_stk_selected_chargeX_tmp = rFile.Get("STKcharge/h_stk_selected_chargeX")
        h_stk_selected_chargeY_tmp = rFile.Get("STKcharge/h_stk_selected_chargeY")
        h_stk_selected_charge_tmp = rFile.Get("STKcharge/h_stk_selected_charge")
        h_stk_selected_charge2D_tmp = rFile.Get("STKcharge/h_stk_selected_charge2D")

        h_background_under_xtrl_cut_tmp = rFile.Get("mc_ancillary/h_background_under_xtrl_cut")
        h_background_over_xtrl_cut_tmp = rFile.Get("mc_ancillary/h_background_over_xtrl_cut")

        # Unlink histos
        h_geo_factor_tmp.SetDirectory(0)
        h_incoming_tmp.SetDirectory(0)
        h_trigger_tmp.SetDirectory(0)
        h_geometric_cut_tmp.SetDirectory(0)
        h_maxElayer_cut_tmp.SetDirectory(0)
        h_maxBarLayer_cut_tmp.SetDirectory(0)
        h_BGOTrackContainment_cut_tmp.SetDirectory(0)
        h_BGO_fiducial_cut_tmp.SetDirectory(0)
        h_all_cut_tmp.SetDirectory(0)

        h_geo_factor_tmp_nw.SetDirectory(0)
        h_incoming_tmp_nw.SetDirectory(0)
        h_trigger_tmp_w.SetDirectory(0)
        h_geometric_cut_tmp_nw.SetDirectory(0)
        h_maxElayer_cut_tmp_nw.SetDirectory(0)
        h_maxBarLayer_cut_tmp_nw.SetDirectory(0)
        h_BGOTrackContainment_cut_tmp_nw.SetDirectory(0)
        h_BGO_fiducial_cut_tmp_nw.SetDirectory(0)
        h_all_cut_tmp_nw.SetDirectory(0)
        
        h_geometric_maxElayer_cut_tmp.SetDirectory(0)
        h_geometric_maxBarLayer_cut_tmp.SetDirectory(0)
        h_geometric_BGOTrackContainment_cut_tmp.SetDirectory(0)
        h_geometric_BGO_fiducial_cut_tmp.SetDirectory(0)
        h_geometric_all_cut_tmp.SetDirectory(0)

        h_geometric_maxElayer_cut_tmp_nw.SetDirectory(0)
        h_geometric_maxBarLayer_cut_tmp_nw.SetDirectory(0)
        h_geometric_BGOTrackContainment_cut_tmp_nw.SetDirectory(0)
        h_geometric_BGO_fiducial_cut_tmp_nw.SetDirectory(0)
        h_geometric_all_cut_tmp_nw.SetDirectory(0)

        h_BGOfiducial_nBarLayer13_cut_tmp.SetDirectory(0)
        h_BGOfiducial_maxRms_cut_tmp.SetDirectory(0)
        h_BGOfiducial_track_selection_cut_tmp.SetDirectory(0)
        h_BGOfiducial_psd_stk_match_cut_tmp.SetDirectory(0)
        h_BGOfiducial_psd_charge_cut_tmp.SetDirectory(0)
        h_BGOfiducial_stk_charge_cut_tmp.SetDirectory(0)
        h_BGOfiducial_xtrl_cut_tmp.SetDirectory(0)
        h_BGOfiducial_all_cut_tmp.SetDirectory(0)

        h_BGOfiducial_nBarLayer13_cut_tmp_nw.SetDirectory(0)
        h_BGOfiducial_maxRms_cut_tmp_nw.SetDirectory(0)
        h_BGOfiducial_track_selection_cut_tmp_nw.SetDirectory(0)
        h_BGOfiducial_psd_stk_match_cut_tmp_nw.SetDirectory(0)
        h_BGOfiducial_psd_charge_cut_tmp_nw.SetDirectory(0)
        h_BGOfiducial_stk_charge_cut_tmp_nw.SetDirectory(0)
        h_BGOfiducial_xtrl_cut_tmp_nw.SetDirectory(0)
        h_BGOfiducial_all_cut_tmp_nw.SetDirectory(0)

        h_preGeo_BGOrec_topX_vs_realX_tmp.SetDirectory(0)
        h_preGeo_BGOrec_topY_vs_realY_tmp.SetDirectory(0)
        h_preGeo_real_slopeX_tmp.SetDirectory(0)
        h_preGeo_real_slopeY_tmp.SetDirectory(0)
        h_preGeo_BGOrec_slopeX_tmp.SetDirectory(0)
        h_preGeo_BGOrec_slopeY_tmp.SetDirectory(0)
        h_preGeo_real_interceptX_tmp.SetDirectory(0)
        h_preGeo_real_interceptY_tmp.SetDirectory(0)
        h_preGeo_BGOrec_interceptX_tmp.SetDirectory(0)
        h_preGeo_BGOrec_interceptY_tmp.SetDirectory(0)
        h_preGeo_real_topMap_tmp.SetDirectory(0)
        h_preGeo_BGOreco_topMap_tmp.SetDirectory(0)
        h_preGeo_real_bottomMap_tmp.SetDirectory(0)
        h_preGeo_BGOreco_bottomMap_tmp.SetDirectory(0)
        h_noBGOenergy_real_topMap_tmp.SetDirectory(0)

        h_geo_BGOrec_topX_vs_realX_tmp.SetDirectory(0)
        h_geo_BGOrec_topY_vs_realY_tmp.SetDirectory(0)
        h_geo_real_slopeX_tmp.SetDirectory(0)
        h_geo_real_slopeY_tmp.SetDirectory(0)
        h_geo_BGOrec_slopeX_tmp.SetDirectory(0)
        h_geo_BGOrec_slopeY_tmp.SetDirectory(0)
        h_geo_real_interceptX_tmp.SetDirectory(0)
        h_geo_real_interceptY_tmp.SetDirectory(0)
        h_geo_BGOrec_interceptX_tmp.SetDirectory(0)
        h_geo_BGOrec_interceptY_tmp.SetDirectory(0)
        h_geo_real_topMap_tmp.SetDirectory(0)
        h_geo_BGOreco_topMap_tmp.SetDirectory(0)
        h_geo_real_bottomMap_tmp.SetDirectory(0)
        h_geo_BGOreco_bottomMap_tmp.SetDirectory(0)

        h_BGOrec_energy_tmp.SetDirectory(0)
        h_simu_energy_tmp.SetDirectory(0)
        h_energy_diff_tmp.SetDirectory(0)
        h_energy_diff2D_tmp.SetDirectory(0)
        h_energy_unfold_tmp.SetDirectory(0)
        h_layer_max_energy_ratio_tmp.SetDirectory(0)
        for idx in range(0,14):
            h_layer_energy_ratio_tmp[idx].SetDirectory(0)
        
        for idx in range(n_energy_bins):
            h_sumRms_cosine_tmp[idx].SetDirectory(0)
        sumRms_cosine_20_100_tmp.SetDirectory(0)
        sumRms_cosine_100_250_tmp.SetDirectory(0)
        sumRms_cosine_250_500_tmp.SetDirectory(0)
        sumRms_cosine_500_1000_tmp.SetDirectory(0)
        sumRms_cosine_1000_3000_tmp.SetDirectory(0)
        sumRms_cosine_3000_10000_tmp.SetDirectory(0)

        h_xtrl_energy_int_tmp.SetDirectory(0)
        h_xtrl_tmp.SetDirectory(0)
        e_discrimination_last_tmp.SetDirectory(0)
        e_discrimination_last_20_100_tmp.SetDirectory(0)
        e_discrimination_last_100_250_tmp.SetDirectory(0)
        e_discrimination_last_250_500_tmp.SetDirectory(0)
        e_discrimination_last_500_1000_tmp.SetDirectory(0)
        e_discrimination_last_1000_3000_tmp.SetDirectory(0)
        e_discrimination_last_3000_10000_tmp.SetDirectory(0)

        e_discrimination_tmp.SetDirectory(0)
        e_discrimination_20_100_tmp.SetDirectory(0)
        e_discrimination_100_250_tmp.SetDirectory(0)
        e_discrimination_250_500_tmp.SetDirectory(0)
        e_discrimination_500_1000_tmp.SetDirectory(0)
        e_discrimination_1000_3000_tmp.SetDirectory(0)
        e_discrimination_3000_10000_tmp.SetDirectory(0)

        for idx in range(n_energy_bins):
            h_bin_xtrl_tmp[idx].SetDirectory(0)

        h_psd_chargeX_tmp.SetDirectory(0)
        h_psd_chargeY_tmp.SetDirectory(0)
        h_psd_charge_tmp.SetDirectory(0)
        h_psd_charge2D_tmp.SetDirectory(0)

        h_psd_selected_chargeX_tmp.SetDirectory(0)
        h_psd_selected_chargeY_tmp.SetDirectory(0)
        h_psd_selected_charge_tmp.SetDirectory(0)
        h_psd_selected_charge2D_tmp.SetDirectory(0)

        h_stk_chargeX_tmp.SetDirectory(0)
        h_stk_chargeY_tmp.SetDirectory(0)
        h_stk_charge_tmp.SetDirectory(0)
        h_stk_charge2D_tmp.SetDirectory(0)

        h_stk_selected_chargeX_tmp.SetDirectory(0)
        h_stk_selected_chargeY_tmp.SetDirectory(0)
        h_stk_selected_charge_tmp.SetDirectory(0)
        h_stk_selected_charge2D_tmp.SetDirectory(0)

        h_background_under_xtrl_cut_tmp.SetDirectory(0)
        h_background_over_xtrl_cut_tmp.SetDirectory(0)

        # Clone output file
        rFile.Close()

        # Add histos
        if dIdx == 0:
            
            h_geo_factor = h_geo_factor_tmp.Clone("h_geo_factor")
            h_incoming = h_incoming_tmp.Clone("h_incoming")
            h_trigger = h_trigger_tmp.Clone("h_trigger")
            h_geometric_cut = h_geometric_cut_tmp.Clone("h_geometric_cut")
            h_maxElayer_cut = h_maxElayer_cut_tmp.Clone("h_maxElayer_cut")
            h_maxBarLayer_cut = h_maxBarLayer_cut_tmp.Clone("h_maxBarLayer_cut")
            h_BGOTrackContainment_cut = h_BGOTrackContainment_cut_tmp.Clone("h_BGOTrackContainment_cut")
            h_BGO_fiducial_cut = h_BGO_fiducial_cut_tmp.Clone("h_BGO_fiducial_cut")
            h_all_cut = h_all_cut_tmp.Clone("h_all_cut")

            h_geo_factor_nw = h_geo_factor_tmp_nw.Clone("h_geo_factor_nw")
            h_incoming_nw = h_incoming_tmp_nw.Clone("h_incoming_nw")
            h_trigger_nw = h_trigger_tmp_nw.Clone("h_trigger_nw")
            h_geometric_cut_nw = h_geometric_cut_tmp_nw.Clone("h_geometric_cut_nw")
            h_maxElayer_cut_nw = h_maxElayer_cut_tmp_nw.Clone("h_maxElayer_cut_nw")
            h_maxBarLayer_cut_nw = h_maxBarLayer_cut_tmp_nw.Clone("h_maxBarLayer_cut_nw")
            h_BGOTrackContainment_cut_nw = h_BGOTrackContainment_cut_tmp_nw.Clone("h_BGOTrackContainment_cut_nw")
            h_BGO_fiducial_cut_nw = h_BGO_fiducial_cut_tmp_nw.Clone("h_BGO_fiducial_cut_nw")
            h_all_cut_nw = h_all_cut_tmp_nw.Clone("h_all_cut_nw")

            h_geometric_maxElayer_cut = h_geometric_maxElayer_cut_tmp.Clone("h_geometric_maxElayer_cut")
            h_geometric_maxBarLayer_cut = h_geometric_maxBarLayer_cut_tmp.Clone("h_geometric_maxBarLayer_cut")
            h_geometric_BGOTrackContainment_cut = h_geometric_BGOTrackContainment_cut_tmp.Clone("h_geometric_BGOTrackContainment_cut")
            h_geometric_BGO_fiducial_cut = h_geometric_BGO_fiducial_cut_tmp.Clone("h_geometric_BGO_fiducial_cut")
            h_geometric_all_cut = h_geometric_all_cut_tmp.Clone("h_geometric_all_cut")

            h_geometric_maxElayer_cut_nw = h_geometric_maxElayer_cut_tmp_nw.Clone("h_geometric_maxElayer_cut_nw")
            h_geometric_maxBarLayer_cut_nw = h_geometric_maxBarLayer_cut_tmp_nw.Clone("h_geometric_maxBarLayer_cut_nw")
            h_geometric_BGOTrackContainment_cut_nw = h_geometric_BGOTrackContainment_cut_tmp_nw.Clone("h_geometric_BGOTrackContainment_cut_nw")
            h_geometric_BGO_fiducial_cut_nw = h_geometric_BGO_fiducial_cut_tmp_nw.Clone("h_geometric_BGO_fiducial_cut_nw")
            h_geometric_all_cut_nw = h_geometric_all_cut_tmp_nw.Clone("h_geometric_all_cut_nw")

            h_BGOfiducial_nBarLayer13_cut = h_BGOfiducial_nBarLayer13_cut_tmp.Clone("h_BGOfiducial_nBarLayer13_cut")
            h_BGOfiducial_maxRms_cut = h_BGOfiducial_maxRms_cut_tmp.Clone("h_BGOfiducial_maxRms_cut")
            h_BGOfiducial_track_selection_cut = h_BGOfiducial_track_selection_cut_tmp.Clone("h_BGOfiducial_track_selection_cut")
            h_BGOfiducial_psd_stk_match_cut = h_BGOfiducial_psd_stk_match_cut_tmp.Clone("h_BGOfiducial_psd_stk_match_cut")
            h_BGOfiducial_psd_charge_cut = h_BGOfiducial_psd_charge_cut_tmp.Clone("h_BGOfiducial_psd_charge_cut")
            h_BGOfiducial_stk_charge_cut = h_BGOfiducial_stk_charge_cut_tmp.Clone("h_BGOfiducial_stk_charge_cut")
            h_BGOfiducial_xtrl_cut = h_BGOfiducial_xtrl_cut_tmp.Clone("h_BGOfiducial_xtrl_cut")
            h_BGOfiducial_all_cut = h_BGOfiducial_all_cut_tmp.Clone("h_BGOfiducial_all_cut")

            h_BGOfiducial_nBarLayer13_cut_nw = h_BGOfiducial_nBarLayer13_cut_tmp_nw.Clone("h_BGOfiducial_nBarLayer13_cut_nw")
            h_BGOfiducial_maxRms_cut_nw = h_BGOfiducial_maxRms_cut_tmp_nw.Clone("h_BGOfiducial_maxRms_cut_nw")
            h_BGOfiducial_track_selection_cut_nw = h_BGOfiducial_track_selection_cut_tmp_nw.Clone("h_BGOfiducial_track_selection_cut_nw")
            h_BGOfiducial_psd_stk_match_cut_nw = h_BGOfiducial_psd_stk_match_cut_tmp_nw.Clone("h_BGOfiducial_psd_stk_match_cut_nw")
            h_BGOfiducial_psd_charge_cut_nw = h_BGOfiducial_psd_charge_cut_tmp_nw.Clone("h_BGOfiducial_psd_charge_cut_nw")
            h_BGOfiducial_stk_charge_cut_nw = h_BGOfiducial_stk_charge_cut_tmp_nw.Clone("h_BGOfiducial_stk_charge_cut_nw")
            h_BGOfiducial_xtrl_cut_nw = h_BGOfiducial_xtrl_cut_tmp_nw.Clone("h_BGOfiducial_xtrl_cut_nw")
            h_BGOfiducial_all_cut_nw = h_BGOfiducial_all_cut_tmp_nw.Clone("h_BGOfiducial_all_cut_nw")

            h_preGeo_BGOrec_topX_vs_realX = h_preGeo_BGOrec_topX_vs_realX_tmp.Clone("h_preGeo_BGOrec_topX_vs_realX")
            h_preGeo_BGOrec_topY_vs_realY = h_preGeo_BGOrec_topY_vs_realY_tmp.Clone("h_preGeo_BGOrec_topY_vs_realY")
            h_preGeo_real_slopeX = h_preGeo_real_slopeX_tmp.Clone("h_preGeo_real_slopeX")
            h_preGeo_real_slopeY = h_preGeo_real_slopeY_tmp.Clone("h_preGeo_real_slopeY")
            h_preGeo_BGOrec_slopeX = h_preGeo_BGOrec_slopeX_tmp.Clone("h_preGeo_BGOrec_slopeX")
            h_preGeo_BGOrec_slopeY = h_preGeo_BGOrec_slopeY_tmp.Clone("h_preGeo_BGOrec_slopeY")
            h_preGeo_real_interceptX = h_preGeo_real_interceptX_tmp.Clone("h_preGeo_real_interceptX")
            h_preGeo_real_interceptY = h_preGeo_real_interceptY_tmp.Clone("h_preGeo_real_interceptY")
            h_preGeo_BGOrec_interceptX = h_preGeo_BGOrec_interceptX_tmp.Clone("h_preGeo_BGOrec_interceptX")
            h_preGeo_BGOrec_interceptY = h_preGeo_BGOrec_interceptY_tmp.Clone("h_preGeo_BGOrec_interceptY")
            h_preGeo_real_topMap = h_preGeo_real_topMap_tmp.Clone("h_preGeo_real_topMap")
            h_preGeo_BGOreco_topMap = h_preGeo_BGOreco_topMap_tmp.Clone("h_preGeo_BGOreco_topMap")
            h_preGeo_real_bottomMap = h_preGeo_real_bottomMap_tmp.Clone("h_preGeo_real_bottomMap")
            h_preGeo_BGOreco_bottomMap = h_preGeo_BGOreco_bottomMap_tmp.Clone("h_preGeo_BGOreco_bottomMap")
            h_noBGOenergy_real_topMap = h_noBGOenergy_real_topMap_tmp.Clone("h_noBGOenergy_real_topMap")

            h_geo_BGOrec_topX_vs_realX = h_geo_BGOrec_topX_vs_realX_tmp.Clone("h_geo_BGOrec_topX_vs_realX")
            h_geo_BGOrec_topY_vs_realY = h_geo_BGOrec_topY_vs_realY_tmp.Clone("h_geo_BGOrec_topY_vs_realY")
            h_geo_real_slopeX = h_geo_real_slopeX_tmp.Clone("h_geo_real_slopeX")
            h_geo_real_slopeY = h_geo_real_slopeY_tmp.Clone("h_geo_real_slopeY")
            h_geo_BGOrec_slopeX = h_geo_BGOrec_slopeX_tmp.Clone("h_geo_BGOrec_slopeX")
            h_geo_BGOrec_slopeY = h_geo_BGOrec_slopeY_tmp.Clone("h_geo_BGOrec_slopeY")
            h_geo_real_interceptX = h_geo_real_interceptX_tmp.Clone("h_geo_real_interceptX")
            h_geo_real_interceptY = h_geo_real_interceptY_tmp.Clone("h_geo_real_interceptY")
            h_geo_BGOrec_interceptX = h_geo_BGOrec_interceptX_tmp.Clone("h_geo_BGOrec_interceptX")
            h_geo_BGOrec_interceptY = h_geo_BGOrec_interceptY_tmp.Clone("h_geo_BGOrec_interceptY")
            h_geo_real_topMap = h_geo_real_topMap_tmp.Clone("h_geo_real_topMap")
            h_geo_BGOreco_topMap = h_geo_BGOreco_topMap_tmp.Clone("h_geo_BGOreco_topMap")
            h_geo_real_bottomMap = h_geo_real_bottomMap_tmp.Clone("h_geo_real_bottomMap")
            h_geo_BGOreco_bottomMap = h_geo_BGOreco_bottomMap_tmp.Clone("h_geo_BGOreco_bottomMap")

            h_BGOrec_energy = h_BGOrec_energy_tmp.Clone("h_BGOrec_energy")
            h_simu_energy = h_simu_energy_tmp.Clone("h_simu_energy")
            h_energy_diff = h_energy_diff_tmp.Clone("h_energy_diff")
            h_energy_diff2D = h_energy_diff2D_tmp.Clone("h_energy_diff2D")
            h_energy_unfold = h_energy_unfold_tmp.Clone("h_energy_unfold")
            h_layer_max_energy_ratio = h_layer_max_energy_ratio_tmp.Clone("h_layer_max_energy_ratio")
            for idx in range(0,14):
                h_ratio_name = "h_layer_energy_ratio_"
                h_ratio_name += str(idx)
                h_layer_energy_ratio[idx] = h_layer_energy_ratio_tmp[idx].Clone(h_ratio_name)

            for idx in range(n_energy_bins):
                histoName = "sumRms_cosine_"
                histoName += str(idx)
                h_sumRms_cosine[idx] = h_sumRms_cosine_tmp.Clone(histoName)
            sumRms_cosine_20_100 = sumRms_cosine_20_100_tmp.Clone("sumRms_cosine_20_100")
            sumRms_cosine_100_250 = sumRms_cosine_100_250_tmp.Clone("sumRms_cosine_100_250")
            sumRms_cosine_250_500 = sumRms_cosine_250_500_tmp.Clone("sumRms_cosine_250_500")
            sumRms_cosine_500_1000 = sumRms_cosine_500_1000_tmp.Clone("sumRms_cosine_500_1000")
            sumRms_cosine_1000_3000 = sumRms_cosine_1000_3000_tmp.Clone("sumRms_cosine_1000_3000")
            sumRms_cosine_3000_10000 = sumRms_cosine_3000_10000_tmp.Clone("sumRms_cosine_3000_10000")

            h_xtrl_energy_int = h_xtrl_energy_int_tmp.Clone("h_xtrl_energy_int")
            h_xtrl = h_xtrl_tmp.Clone("h_xtrl")
            e_discrimination_last = e_discrimination_last_tmp.Clone("e_discrimination_last")
            e_discrimination_last_20_100 = e_discrimination_last_20_100_tmp.Clone("e_discrimination_last_20_100")
            e_discrimination_last_100_250 = e_discrimination_last_100_250_tmp.Clone("e_discrimination_last_100_250")
            e_discrimination_last_250_500 = e_discrimination_last_250_500_tmp.Clone("e_discrimination_last_250_500")
            e_discrimination_last_500_1000 = e_discrimination_last_500_1000_tmp.Clone("e_discrimination_last_500_1000")
            e_discrimination_last_1000_3000 = e_discrimination_last_1000_3000_tmp.Clone("e_discrimination_last_1000_3000")
            e_discrimination_last_3000_10000 = e_discrimination_last_3000_10000_tmp.Clone("e_discrimination_last_3000_10000")

            e_discrimination = e_discrimination_tmp.Clone("e_discrimination")
            e_discrimination_20_100 = e_discrimination_20_100_tmp.Clone("e_discrimination_20_100")
            e_discrimination_100_250 = e_discrimination_100_250_tmp.Clone("e_discrimination_100_250")
            e_discrimination_250_500 = e_discrimination_250_500_tmp.Clone("e_discrimination_250_500")
            e_discrimination_500_1000 = e_discrimination_500_1000_tmp.Clone("e_discrimination_500_1000")
            e_discrimination_1000_3000 = e_discrimination_1000_3000_tmp.Clone("e_discrimination_1000_3000")
            e_discrimination_3000_10000 = e_discrimination_3000_10000_tmp.Clone("e_discrimination_3000_10000")

            for idx in range(n_energy_bins):
                histoName = "h_xtrl_bin_"
                histoName += str(idx)
                h_bin_xtrl[idx] = h_bin_xtrl_tmp[idx].Clone(histoName)

            h_psd_chargeX = h_psd_chargeX_tmp.Clone("h_psd_chargeX")
            h_psd_chargeY = h_psd_chargeY_tmp.Clone("h_psd_chargeY")
            h_psd_charge = h_psd_charge_tmp.Clone("h_psd_charge")
            h_psd_charge2D = h_psd_charge2D_tmp.Clone("h_psd_charge2D")

            h_psd_selected_chargeX = h_psd_selected_chargeX_tmp.Clone("h_psd_selected_chargeX")
            h_psd_selected_chargeY = h_psd_selected_chargeY_tmp.Clone("h_psd_selected_chargeY")
            h_psd_selected_charge = h_psd_selected_charge_tmp.Clone("h_psd_selected_charge")
            h_psd_selected_charge2D = h_psd_selected_charge2D_tmp.Clone("h_psd_selected_charge2D")

            h_stk_chargeX = h_stk_chargeX_tmp.Clone("h_stk_chargeX")
            h_stk_chargeY = h_stk_chargeY_tmp.Clone("h_stk_chargeY")
            h_stk_charge = h_stk_charge_tmp.Clone("h_stk_charge")
            h_stk_charge2D = h_stk_charge2D_tmp.Clone("h_stk_charge2D")

            h_stk_selected_chargeX = h_stk_selected_chargeX_tmp.Clone("h_stk_selected_chargeX")
            h_stk_selected_chargeY = h_stk_selected_chargeY_tmp.Clone("h_stk_selected_chargeY")
            h_stk_selected_charge = h_stk_selected_charge_tmp.Clone("h_stk_selected_charge")
            h_stk_selected_charge2D = h_stk_selected_charge2D_tmp.Clone("h_stk_selected_charge2D")
            
            h_background_under_xtrl_cut = h_background_under_xtrl_cut_tmp.Clone("h_background_under_xtrl_cut")
            h_background_over_xtrl_cut = h_background_over_xtrl_cut_tmp.Clone("h_background_over_xtrl_cut")

        else:

            h_geo_factor.Add(h_geo_factor_tmp)
            h_incoming.Add(h_incoming_tmp)
            h_trigger.Add(h_trigger_tmp)
            h_geometric_cut.Add(h_geometric_cut_tmp)
            h_maxElayer_cut.Add(h_maxElayer_cut_tmp)
            h_maxBarLayer_cut.Add(h_maxBarLayer_cut_tmp)
            h_BGOTrackContainment_cut.Add(h_BGOTrackContainment_cut_tmp)
            h_BGO_fiducial_cut.Add(h_BGO_fiducial_cut_tmp)
            h_all_cut.Add(h_all_cut_tmp)

            h_geo_factor_nw.Add(h_geo_factor_tmp_nw)
            h_incoming_nw.Add(h_incoming_tmp_nw)
            h_trigger_nw.Add(h_trigger_tmp_nw)
            h_geometric_cut_nw.Add(h_geometric_cut_tmp_nw)
            h_maxElayer_cut_nw.Add(h_maxElayer_cut_tmp_nw)
            h_maxBarLayer_cut_nw.Add(h_maxBarLayer_cut_tmp_nw)
            h_BGOTrackContainment_cut_nw.Add(h_BGOTrackContainment_cut_tmp_nw)
            h_BGO_fiducial_cut_nw.Add(h_BGO_fiducial_cut_tmp_nw)
            h_all_cut_nw.Add(h_all_cut_tmp_nw)

            h_geometric_maxElayer_cut.Add(h_geometric_maxElayer_cut_tmp)
            h_geometric_maxBarLayer_cut.Add(h_geometric_maxBarLayer_cut_tmp)
            h_geometric_BGOTrackContainment_cut.Add(h_geometric_BGOTrackContainment_cut_tmp)
            h_geometric_BGO_fiducial_cut.Add(h_geometric_BGO_fiducial_cut_tmp)
            h_geometric_all_cut.Add(h_geometric_all_cut_tmp)

            h_geometric_maxElayer_cut_nw.Add(h_geometric_maxElayer_cut_tmp_nw)
            h_geometric_maxBarLayer_cut_nw.Add(h_geometric_maxBarLayer_cut_tmp_nw)
            h_geometric_BGOTrackContainment_cut_nw.Add(h_geometric_BGOTrackContainment_cut_tmp_nw)
            h_geometric_BGO_fiducial_cut_nw.Add(h_geometric_BGO_fiducial_cut_tmp_nw)
            h_geometric_all_cut_nw.Add(h_geometric_all_cut_tmp_nw)

            h_BGOfiducial_nBarLayer13_cut.Add(h_BGOfiducial_nBarLayer13_cut_tmp)
            h_BGOfiducial_maxRms_cut.Add(h_BGOfiducial_maxRms_cut_tmp)
            h_BGOfiducial_track_selection_cut.Add(h_BGOfiducial_track_selection_cut_tmp)
            h_BGOfiducial_psd_stk_match_cut.Add(h_BGOfiducial_psd_stk_match_cut_tmp)
            h_BGOfiducial_psd_charge_cut.Add(h_BGOfiducial_psd_charge_cut_tmp)
            h_BGOfiducial_stk_charge_cut.Add(h_BGOfiducial_stk_charge_cut_tmp)
            h_BGOfiducial_xtrl_cut.Add(h_BGOfiducial_xtrl_cut_tmp)
            h_BGOfiducial_all_cut.Add(h_BGOfiducial_all_cut_tmp)

            h_BGOfiducial_nBarLayer13_cut_nw.Add(h_BGOfiducial_nBarLayer13_cut_tmp_nw)
            h_BGOfiducial_maxRms_cut_nw.Add(h_BGOfiducial_maxRms_cut_tmp_nw)
            h_BGOfiducial_track_selection_cut_nw.Add(h_BGOfiducial_track_selection_cut_tmp_nw)
            h_BGOfiducial_psd_stk_match_cut_nw.Add(h_BGOfiducial_psd_stk_match_cut_tmp_nw)
            h_BGOfiducial_psd_charge_cut_nw.Add(h_BGOfiducial_psd_charge_cut_tmp_nw)
            h_BGOfiducial_stk_charge_cut_nw.Add(h_BGOfiducial_stk_charge_cut_tmp_nw)
            h_BGOfiducial_xtrl_cut_nw.Add(h_BGOfiducial_xtrl_cut_tmp_nw)
            h_BGOfiducial_all_cut_nw.Add(h_BGOfiducial_all_cut_tmp_nw)

            h_preGeo_BGOrec_topX_vs_realX.Add(h_preGeo_BGOrec_topX_vs_realX_tmp)
            h_preGeo_BGOrec_topY_vs_realY.Add(h_preGeo_BGOrec_topY_vs_realY_tmp)
            h_preGeo_real_slopeX.Add(h_preGeo_real_slopeX_tmp)
            h_preGeo_real_slopeY.Add(h_preGeo_real_slopeY_tmp)
            h_preGeo_BGOrec_slopeX.Add(h_preGeo_BGOrec_slopeX_tmp)
            h_preGeo_BGOrec_slopeY.Add(h_preGeo_BGOrec_slopeY_tmp)
            h_preGeo_real_interceptX.Add(h_preGeo_real_interceptX_tmp)
            h_preGeo_real_interceptY.Add(h_preGeo_real_interceptY_tmp)
            h_preGeo_BGOrec_interceptX.Add(h_preGeo_BGOrec_interceptX_tmp)
            h_preGeo_BGOrec_interceptY.Add(h_preGeo_BGOrec_interceptY_tmp)
            h_preGeo_real_topMap.Add(h_preGeo_real_topMap_tmp)
            h_preGeo_BGOreco_topMap.Add(h_preGeo_BGOreco_topMap_tmp)
            h_preGeo_real_bottomMap.Add(h_preGeo_real_bottomMap_tmp)
            h_preGeo_BGOreco_bottomMap.Add(h_preGeo_BGOreco_bottomMap_tmp)
            h_noBGOenergy_real_topMap.Add(h_noBGOenergy_real_topMap_tmp)

            h_geo_BGOrec_topX_vs_realX.Add(h_geo_BGOrec_topX_vs_realX_tmp)
            h_geo_BGOrec_topY_vs_realY.Add(h_geo_BGOrec_topY_vs_realY_tmp)
            h_geo_real_slopeX.Add(h_geo_real_slopeX_tmp)
            h_geo_real_slopeY.Add(h_geo_real_slopeY_tmp)
            h_geo_BGOrec_slopeX.Add(h_geo_BGOrec_slopeX_tmp)
            h_geo_BGOrec_slopeY.Add(h_geo_BGOrec_slopeY_tmp)
            h_geo_real_interceptX.Add(h_geo_real_interceptX_tmp)
            h_geo_real_interceptY.Add(h_geo_real_interceptY_tmp)
            h_geo_BGOrec_interceptX.Add(h_geo_BGOrec_interceptX_tmp)
            h_geo_BGOrec_interceptY.Add(h_geo_BGOrec_interceptY_tmp)
            h_geo_real_topMap.Add(h_geo_real_topMap_tmp)
            h_geo_BGOreco_topMap.Add(h_geo_BGOreco_topMap_tmp)
            h_geo_real_bottomMap.Add(h_geo_real_bottomMap_tmp)
            h_geo_BGOreco_bottomMap.Add(h_geo_BGOreco_bottomMap_tmp)

            h_BGOrec_energy.Add(h_BGOrec_energy_tmp)
            h_simu_energy.Add(h_simu_energy_tmp)
            h_energy_diff.Add(h_energy_diff_tmp)
            h_energy_diff2D.Add(h_energy_diff2D_tmp)
            h_energy_unfold.Add(h_energy_unfold_tmp)
            h_layer_max_energy_ratio.Add(h_layer_max_energy_ratio_tmp)
            for idx in range(0,14):
                h_layer_energy_ratio[idx].Add(h_layer_energy_ratio_tmp[idx])

            for idx in range(n_energy_bins):
                h_sumRms_cosine[idx].Add(h_sumRms_cosine_tmp[idx])
            sumRms_cosine_20_100.Add(sumRms_cosine_20_100_tmp)
            sumRms_cosine_100_250.Add(sumRms_cosine_100_250_tmp)
            sumRms_cosine_250_500.Add(sumRms_cosine_250_500_tmp)
            sumRms_cosine_500_1000.Add(sumRms_cosine_500_1000_tmp)
            sumRms_cosine_1000_3000.Add(sumRms_cosine_1000_3000_tmp)
            sumRms_cosine_3000_10000.Add(sumRms_cosine_3000_10000_tmp)

            h_xtrl_energy_int.Add(h_xtrl_energy_int_tmp)
            h_xtrl.Add(h_xtrl_tmp)

            e_discrimination_last.Add(e_discrimination_last_tmp)
            e_discrimination_last_20_100.Add(e_discrimination_last_20_100_tmp)
            e_discrimination_last_100_250.Add(e_discrimination_last_100_250_tmp)
            e_discrimination_last_250_500.Add(e_discrimination_last_250_500_tmp)
            e_discrimination_last_500_1000.Add(e_discrimination_last_500_1000_tmp)
            e_discrimination_last_1000_3000.Add(e_discrimination_last_1000_3000_tmp)
            e_discrimination_last_3000_10000.Add(e_discrimination_last_3000_10000_tmp)

            e_discrimination.Add(e_discrimination_tmp)
            e_discrimination_20_100.Add(e_discrimination_20_100_tmp)
            e_discrimination_100_250.Add(e_discrimination_100_250_tmp)
            e_discrimination_250_500.Add(e_discrimination_250_500_tmp)
            e_discrimination_500_1000.Add(e_discrimination_500_1000_tmp)
            e_discrimination_1000_3000.Add(e_discrimination_1000_3000_tmp)
            e_discrimination_3000_10000.Add(e_discrimination_3000_10000_tmp)

            for idx in range(n_energy_bins):
                h_bin_xtrl[idx].Add(h_bin_xtrl_tmp[idx])

            h_psd_chargeX.Add(h_psd_chargeX_tmp)
            h_psd_chargeY.Add(h_psd_chargeY_tmp)
            h_psd_charge.Add(h_psd_charge_tmp)
            h_psd_charge2D.Add(h_psd_charge2D_tmp)

            h_psd_selected_chargeX.Add(h_psd_selected_chargeX_tmp)
            h_psd_selected_chargeY.Add(h_psd_selected_chargeY_tmp)
            h_psd_selected_charge.Add(h_psd_selected_charge_tmp)
            h_psd_selected_charge2D.Add(h_psd_selected_charge2D_tmp)

            h_stk_chargeX.Add(h_stk_chargeX_tmp)
            h_stk_chargeY.Add(h_stk_chargeY_tmp)
            h_stk_charge.Add(h_stk_charge_tmp)
            h_stk_charge2D.Add(h_stk_charge2D_tmp)

            h_stk_selected_chargeX.Add(h_stk_selected_chargeX_tmp)
            h_stk_selected_chargeY.Add(h_stk_selected_chargeY_tmp)
            h_stk_selected_charge.Add(h_stk_selected_charge_tmp)
            h_stk_selected_charge2D.Add(h_stk_selected_charge2D_tmp)

            h_background_under_xtrl_cut.Add(h_background_under_xtrl_cut_tmp)
            h_background_over_xtrl_cut.Add(h_background_over_xtrl_cut_tmp)

    # Create output file for full histos
    fOut = TFile.Open(opts.output, "RECREATE")
    if fOut.IsOpen():
        if opts.verbose:
            print('Output TFile has been created: {}'.format(opts.output))
    else:
        print('Error creating output TFile: {}'.format(opts.output))
        sys.exit()

    # Writing final histos to file
    h_geo_factor.Write()
    h_incoming.Write()
    h_trigger.Write()
    h_geometric_cut.Write()
    h_maxElayer_cut.Write()
    h_maxBarLayer_cut.Write()
    h_BGOTrackContainment_cut.Write()
    h_BGO_fiducial_cut.Write()
    h_all_cut.Write()
   
    h_geometric_maxElayer_cut.Write()
    h_geometric_maxBarLayer_cut.Write()
    h_geometric_BGOTrackContainment_cut.Write()
    h_geometric_BGO_fiducial_cut.Write()
    h_geometric_all_cut.Write()

    h_BGOfiducial_nBarLayer13_cut.Write()
    h_BGOfiducial_maxRms_cut.Write()
    h_BGOfiducial_track_selection_cut.Write()
    h_BGOfiducial_psd_stk_match_cut.Write()
    h_BGOfiducial_psd_charge_cut.Write()
    h_BGOfiducial_stk_charge_cut.Write()
    h_BGOfiducial_xtrl_cut.Write()
    h_BGOfiducial_all_cut.Write()
    
    fOut.mkdir("nw_histos")
    fOut.cd("nw_histos")

    h_geo_factor_nw.Write()
    h_incoming_nw.Write()
    h_trigger_nw.Write()
    h_geometric_cut_nw.Write()
    h_maxElayer_cut_nw.Write()
    h_maxBarLayer_cut_nw.Write()
    h_BGOTrackContainment_cut_nw.Write()
    h_BGO_fiducial_cut_nw.Write()
    h_all_cut_nw.Write()
   
    h_geometric_maxElayer_cut_nw.Write()
    h_geometric_maxBarLayer_cut_nw.Write()
    h_geometric_BGOTrackContainment_cut_nw.Write()
    h_geometric_BGO_fiducial_cut_nw.Write()
    h_geometric_all_cut_nw.Write()

    h_BGOfiducial_nBarLayer13_cut_nw.Write()
    h_BGOfiducial_maxRms_cut_nw.Write()
    h_BGOfiducial_track_selection_cut_nw.Write()
    h_BGOfiducial_psd_stk_match_cut_nw.Write()
    h_BGOfiducial_psd_charge_cut_nw.Write()
    h_BGOfiducial_stk_charge_cut_nw.Write()
    h_BGOfiducial_xtrl_cut_nw.Write()
    h_BGOfiducial_all_cut_nw.Write()

    fOut.mkdir("Analysis_preGeoCut")
    fOut.cd("Analysis_preGeoCut")

    h_preGeo_BGOrec_topX_vs_realX.Write()
    h_preGeo_BGOrec_topY_vs_realY.Write()
    h_preGeo_real_slopeX.Write()
    h_preGeo_real_slopeY.Write()
    h_preGeo_BGOrec_slopeX.Write()
    h_preGeo_BGOrec_slopeY.Write()
    h_preGeo_real_interceptX.Write()
    h_preGeo_real_interceptY.Write()
    h_preGeo_BGOrec_interceptX.Write()
    h_preGeo_BGOrec_interceptY.Write()
    h_preGeo_real_topMap.Write()
    h_preGeo_BGOreco_topMap.Write()
    h_preGeo_real_bottomMap.Write()
    h_preGeo_BGOreco_bottomMap.Write()
    h_noBGOenergy_real_topMap.Write()
    
    fOut.mkdir("Analysis_GeoCut")
    fOut.cd("Analysis_GeoCut")

    h_geo_BGOrec_topX_vs_realX.Write()
    h_geo_BGOrec_topY_vs_realY.Write()
    h_geo_real_slopeX.Write()
    h_geo_real_slopeY.Write()
    h_geo_BGOrec_slopeX.Write()
    h_geo_BGOrec_slopeY.Write()
    h_geo_real_interceptX.Write()
    h_geo_real_interceptY.Write()
    h_geo_BGOrec_interceptX.Write()
    h_geo_BGOrec_interceptY.Write()
    h_geo_real_topMap.Write()
    h_geo_BGOreco_topMap.Write()
    h_geo_real_bottomMap.Write()
    h_geo_BGOreco_bottomMap.Write()
    
    fOut.mkdir("BGO_Energy")
    fOut.cd("BGO_Energy")

    h_BGOrec_energy.Write()
    h_simu_energy.Write()
    h_energy_diff.Write()
    h_energy_diff2D.Write()
    h_energy_unfold.Write()
    h_layer_max_energy_ratio.Write()
    
    for idx in range(0,14):
        h_layer_energy_ratio[idx].Write()
    
    for idx in range(n_energy_bins):
        h_sumRms_cosine[idx].Write()
    
    sumRms_cosine_20_100.Write()
    sumRms_cosine_100_250.Write()
    sumRms_cosine_250_500.Write()
    sumRms_cosine_500_1000.Write()
    sumRms_cosine_1000_3000.Write()
    sumRms_cosine_3000_10000.Write()

    fOut.mkdir("xtrl")
    fOut.cd("xtrl")

    h_xtrl_energy_int.Write()
    h_xtrl.Write()   
    e_discrimination_last.Write()
    e_discrimination_last_20_100.Write()
    e_discrimination_last_100_250.Write()
    e_discrimination_last_250_500.Write()
    e_discrimination_last_500_1000.Write()
    e_discrimination_last_1000_3000.Write()
    e_discrimination_last_3000_10000.Write()
    e_discrimination.Write()
    e_discrimination_20_100.Write()
    e_discrimination_100_250.Write()
    e_discrimination_250_500.Write()
    e_discrimination_500_1000.Write()
    e_discrimination_1000_3000.Write()
    e_discrimination_3000_10000.Write()

    for idx in range(n_energy_bins):
        h_bin_xtrl[idx].Write()

    fOut.mkdir("PSDcharge")
    fOut.cd("PSDcharge")

    h_psd_chargeX.Write()
    h_psd_chargeY.Write()
    h_psd_charge.Write()
    h_psd_charge2D.Write()

    h_psd_selected_chargeX.Write()
    h_psd_selected_chargeY.Write()
    h_psd_selected_charge.Write()
    h_psd_selected_charge2D.Write()

    fOut.mkdir("STKcharge")
    fOut.cd("STKcharge")

    h_stk_chargeX.Write()
    h_stk_chargeY.Write()
    h_stk_charge.Write()
    h_stk_charge2D.Write()

    h_stk_selected_chargeX.Write()
    h_stk_selected_chargeY.Write()
    h_stk_selected_charge.Write()
    h_stk_selected_charge2D.Write()

    fOut.mkdir("mc_ancillary")
    fOut.cd("mc_ancillary")

    h_background_under_xtrl_cut.Write()
    h_background_over_xtrl_cut.Write()

    # Closing output file
    fOut.Close()


def compute_final_histos_data(condor_dir_list, opts):
    
    # Cut histos
    h_trigger = TH1D()
    h_geometric_cut = TH1D()
    h_maxElayer_cut = TH1D()
    h_maxBarLayer_cut = TH1D()
    h_BGOTrackContainment_cut = TH1D()
    h_BGO_fiducial_cut = TH1D()
    h_all_cut = TH1D()

    # Cuts && Geometric Cut
    h_geometric_maxElayer_cut = TH1D()
    h_geometric_maxBarLayer_cut = TH1D()
    h_geometric_BGOTrackContainment_cut = TH1D()
    h_geometric_BGO_fiducial_cut = TH1D()
    h_geometric_all_cut = TH1D()

    # Cuts && BGO fiducial volume cut
    h_BGOfiducial_nBarLayer13_cut = TH1D()
    h_BGOfiducial_maxRms_cut = TH1D()
    h_BGOfiducial_track_selection_cut = TH1D()
    h_BGOfiducial_psd_stk_match_cut = TH1D()
    h_BGOfiducial_psd_charge_cut = TH1D()
    h_BGOfiducial_stk_charge_cut = TH1D()
    h_BGOfiducial_xtrl_cut = TH1D()
    h_BGOfiducial_all_cut = TH1D()

    # Analysis histos - reco energy of incoming events
    h_BGOrec_energy = TH1D()
   
    # Ratio of layer energy respect to total BGO energy
    h_layer_max_energy_ratio = TH1D()
    
    h_layer_energy_ratio = []
    for idx in range(0,14):
        h_layer_energy_ratio.append(TH1D())

    # After Geometric Cut
    # Slope X and Y
    h_geo_BGOrec_slopeX = TH1D()
    h_geo_BGOrec_slopeY = TH1D()
    # Intercept X and Y
    h_geo_BGOrec_interceptX = TH1D()
    h_geo_BGOrec_interceptY = TH1D()
    # Top Maps
    h_geo_BGOreco_topMap = TH2D()
    # Bottom Maps
    h_geo_BGOreco_bottomMap = TH2D()

    # sumRms - cosine correlation
    h_sumRms_cosine = []
    sumRms_cosine_20_100 = TH2D()
    sumRms_cosine_100_250 = TH2D()
    sumRms_cosine_250_500 = TH2D()
    sumRms_cosine_500_1000 = TH2D()
    sumRms_cosine_1000_3000 = TH2D()
    sumRms_cosine_3000_10000 = TH2D()

    # XTRL histos
    h_xtrl_energy_int = TH1D()
    h_xtrl = TH2D()
    e_discrimination_last = TH2D()
    e_discrimination_last_20_100 = TH2D()
    e_discrimination_last_100_250 = TH2D()
    e_discrimination_last_250_500 = TH2D()
    e_discrimination_last_500_1000 = TH2D()
    e_discrimination_last_1000_3000 = TH2D()
    e_discrimination_last_3000_10000 = TH2D()
    e_discrimination = TH2D()
    e_discrimination_20_100 = TH2D()
    e_discrimination_100_250 = TH2D()
    e_discrimination_250_500 = TH2D()
    e_discrimination_500_1000 = TH2D()
    e_discrimination_1000_3000 = TH2D()
    e_discrimination_3000_10000 = TH2D()
    h_bin_xtrl = []

    # PSD charge histos
    h_psd_chargeX = TH1D()
    h_psd_chargeY = TH1D()
    h_psd_charge = TH1D()
    h_psd_charge2D = TH2D()

    h_psd_selected_chargeX = TH1D()
    h_psd_selected_chargeY = TH1D()
    h_psd_selected_charge = TH1D()
    h_psd_selected_charge2D = TH2D()

    # STK charge histos
    h_stk_chargeX = TH1D()
    h_stk_chargeY = TH1D()
    h_stk_charge = TH1D()
    h_stk_charge2D = TH2D()

    h_stk_selected_chargeX = TH1D()
    h_stk_selected_chargeY = TH1D()
    h_stk_selected_charge = TH1D()
    h_stk_selected_charge2D = TH2D()

    # Proton background
    h_background_under_xtrl_cut = TH1D()
    h_background_over_xtrl_cut = TH1D()

    for dIdx, tmp_dir in enumerate(condor_dir_list):
        tmp_dir += "/outFiles"
        tmp_dir_list = os.listdir(tmp_dir)
        for elm in tmp_dir_list:
            if elm.startswith("analysisOutFile_"):
                rFile_path = tmp_dir + "/" + elm

        # Open ROOT output file
        rFile = TFile.Open(rFile_path, "READ")
        if rFile.IsOpen():
            if opts.verbose:
                if dIdx == 0:
                    print('\nReading file {}: {}'.format((dIdx+1), rFile_path))
                else:
                    print('Reading file {}: {}'.format((dIdx+1), rFile_path))
        else:
            print('Error reading file {}: {}'.format((dIdx+1), rFile_path))
            sys.exit()
        
        # Check the keys
        if rFile.GetNkeys() == 0:
            continue
        
        # Reading histos
        h_trigger_tmp = rFile.Get("h_trigger")
        h_geometric_cut_tmp = rFile.Get("h_geometric_cut")
        h_maxElayer_cut_tmp = rFile.Get("h_maxElayer_cut")
        h_maxBarLayer_cut_tmp = rFile.Get("h_maxBarLayer_cut")
        h_BGOTrackContainment_cut_tmp = rFile.Get("h_BGOTrackContainment_cut")
        h_BGO_fiducial_cut_tmp = rFile.Get("h_BGO_fiducial_cut")
        h_all_cut_tmp = rFile.Get("h_all_cut")

        n_energy_bins = h_trigger_tmp.GetNbinsX()
        # Init bin_xtrl
        for idx in range(n_energy_bins):
            h_bin_xtrl.append(TH1D())
            h_sumRms_cosine.append(TH2D())

        h_geometric_maxElayer_cut_tmp = rFile.Get("h_geometric_maxElayer_cut")
        h_geometric_maxBarLayer_cut_tmp = rFile.Get("h_geometric_maxBarLayer_cut")
        h_geometric_BGOTrackContainment_cut_tmp = rFile.Get("h_geometric_BGOTrackContainment_cut")
        h_geometric_BGO_fiducial_cut_tmp = rFile.Get("h_geometric_BGO_fiducial_cut")
        h_geometric_all_cut_tmp = rFile.Get("h_geometric_all_cut")

        h_BGOfiducial_nBarLayer13_cut_tmp = rFile.Get("h_BGOfiducial_nBarLayer13_cut")
        h_BGOfiducial_maxRms_cut_tmp = rFile.Get("h_BGOfiducial_maxRms_cut")
        h_BGOfiducial_track_selection_cut_tmp = rFile.Get("h_BGOfiducial_track_selection_cut")
        h_BGOfiducial_psd_stk_match_cut_tmp = rFile.Get("h_BGOfiducial_psd_stk_match_cut")
        h_BGOfiducial_psd_charge_cut_tmp = rFile.Get("h_BGOfiducial_psd_charge_cut")
        h_BGOfiducial_stk_charge_cut_tmp = rFile.Get("h_BGOfiducial_stk_charge_cut")
        h_BGOfiducial_xtrl_cut_tmp = rFile.Get("h_BGOfiducial_xtrl_cut")
        h_BGOfiducial_all_cut_tmp = rFile.Get("h_BGOfiducial_all_cut")
        
        h_geo_BGOrec_slopeX_tmp = rFile.Get("Analysis_GeoCut/h_geo_BGOrec_slopeX")
        h_geo_BGOrec_slopeY_tmp = rFile.Get("Analysis_GeoCut/h_geo_BGOrec_slopeY")
        h_geo_BGOrec_interceptX_tmp = rFile.Get("Analysis_GeoCut/h_geo_BGOrec_interceptX")
        h_geo_BGOrec_interceptY_tmp = rFile.Get("Analysis_GeoCut/h_geo_BGOrec_interceptY")
        h_geo_BGOreco_topMap_tmp = rFile.Get("Analysis_GeoCut/h_geo_BGOreco_topMap")
        h_geo_BGOreco_bottomMap_tmp = rFile.Get("Analysis_GeoCut/h_geo_BGOreco_bottomMap")
        
        h_BGOrec_energy_tmp = rFile.Get("BGO_Energy/h_BGOrec_energy")
        h_layer_max_energy_ratio_tmp = rFile.Get("BGO_Energy/h_layer_max_energy_ratio")
        h_layer_energy_ratio_tmp = []
        for idx in range(0,14):
            histoName = "BGO_Energy/h_layer_energy_ratio_"
            histoName += str(idx)
            h_layer_energy_ratio_tmp.append(rFile.Get(histoName))

        h_sumRms_cosine_tmp = []
        for idx in range(n_energy_bins):
            histoName = "BGO_Energy/sumRms_cosine_"
            histoName += str(idx)
            h_sumRms_cosine_tmp.append(rFile.Get(histoName))
        
        sumRms_cosine_20_100_tmp = rFile.Get("BGO_Energy/sumRms_cosine_20_100")
        sumRms_cosine_100_250_tmp = rFile.Get("BGO_Energy/sumRms_cosine_100_250")
        sumRms_cosine_250_500_tmp = rFile.Get("BGO_Energy/sumRms_cosine_250_500")
        sumRms_cosine_500_1000_tmp = rFile.Get("BGO_Energy/sumRms_cosine_500_1000")
        sumRms_cosine_1000_3000_tmp = rFile.Get("BGO_Energy/sumRms_cosine_1000_3000")
        sumRms_cosine_3000_10000_tmp = rFile.Get("BGO_Energy/sumRms_cosine_3000_10000")

        h_xtrl_energy_int_tmp = rFile.Get("xtrl/h_xtrl_energy_int")
        h_xtrl_tmp = rFile.Get("xtrl/h_xtrl")
        e_discrimination_last_tmp = rFile.Get("xtrl/e_discrimination_last")
        e_discrimination_last_20_100_tmp = rFile.Get("xtrl/e_discrimination_last_20_100")
        e_discrimination_last_100_250_tmp = rFile.Get("xtrl/e_discrimination_last_100_250")
        e_discrimination_last_250_500_tmp = rFile.Get("xtrl/e_discrimination_last_250_500")
        e_discrimination_last_500_1000_tmp = rFile.Get("xtrl/e_discrimination_last_500_1000")
        e_discrimination_last_1000_3000_tmp = rFile.Get("xtrl/e_discrimination_last_1000_3000")
        e_discrimination_last_3000_10000_tmp = rFile.Get("xtrl/e_discrimination_last_3000_10000")

        e_discrimination_tmp = rFile.Get("xtrl/e_discrimination")
        e_discrimination_20_100_tmp = rFile.Get("xtrl/e_discrimination_20_100")
        e_discrimination_100_250_tmp = rFile.Get("xtrl/e_discrimination_100_250")
        e_discrimination_250_500_tmp = rFile.Get("xtrl/e_discrimination_250_500")
        e_discrimination_500_1000_tmp = rFile.Get("xtrl/e_discrimination_500_1000")
        e_discrimination_1000_3000_tmp = rFile.Get("xtrl/e_discrimination_1000_3000")
        e_discrimination_3000_10000_tmp = rFile.Get("xtrl/e_discrimination_3000_10000")

        h_bin_xtrl_tmp = []
        for idx in range(n_energy_bins):
            histoName = "xtrl/h_xtrl_bin_"
            histoName += str(idx)
            h_bin_xtrl_tmp.append(rFile.Get(histoName))

        h_psd_chargeX_tmp = rFile.Get("PSDcharge/h_psd_chargeX")
        h_psd_chargeY_tmp = rFile.Get("PSDcharge/h_psd_chargeY")
        h_psd_charge_tmp = rFile.Get("PSDcharge/h_psd_charge")
        h_psd_charge2D_tmp = rFile.Get("PSDcharge/h_psd_charge2D")

        h_psd_selected_chargeX_tmp = rFile.Get("PSDcharge/h_psd_selected_chargeX")
        h_psd_selected_chargeY_tmp = rFile.Get("PSDcharge/h_psd_selected_chargeY")
        h_psd_selected_charge_tmp = rFile.Get("PSDcharge/h_psd_selected_charge")
        h_psd_selected_charge2D_tmp = rFile.Get("PSDcharge/h_psd_selected_charge2D")

        h_stk_chargeX_tmp = rFile.Get("STKcharge/h_stk_chargeX")
        h_stk_chargeY_tmp = rFile.Get("STKcharge/h_stk_chargeY")
        h_stk_charge_tmp = rFile.Get("STKcharge/h_stk_charge")
        h_stk_charge2D_tmp = rFile.Get("STKcharge/h_stk_charge2D")

        h_stk_selected_chargeX_tmp = rFile.Get("STKcharge/h_stk_selected_chargeX")
        h_stk_selected_chargeY_tmp = rFile.Get("STKcharge/h_stk_selected_chargeY")
        h_stk_selected_charge_tmp = rFile.Get("STKcharge/h_stk_selected_charge")
        h_stk_selected_charge2D_tmp = rFile.Get("STKcharge/h_stk_selected_charge2D")

        h_background_under_xtrl_cut_tmp = rFile.Get("proton_background/h_background_under_xtrl_cut")
        h_background_over_xtrl_cut_tmp = rFile.Get("proton_background/h_background_over_xtrl_cut")

        # Unlink histos
        h_trigger_tmp.SetDirectory(0)
        h_geometric_cut_tmp.SetDirectory(0)
        h_maxElayer_cut_tmp.SetDirectory(0)
        h_maxBarLayer_cut_tmp.SetDirectory(0)
        h_BGOTrackContainment_cut_tmp.SetDirectory(0)
        h_BGO_fiducial_cut_tmp.SetDirectory(0)
        h_all_cut_tmp.SetDirectory(0)
        
        h_geometric_maxElayer_cut_tmp.SetDirectory(0)
        h_geometric_maxBarLayer_cut_tmp.SetDirectory(0)
        h_geometric_BGOTrackContainment_cut_tmp.SetDirectory(0)
        h_geometric_BGO_fiducial_cut_tmp.SetDirectory(0)
        h_geometric_all_cut_tmp.SetDirectory(0)

        h_BGOfiducial_nBarLayer13_cut_tmp.SetDirectory(0)
        h_BGOfiducial_maxRms_cut_tmp.SetDirectory(0)
        h_BGOfiducial_track_selection_cut_tmp.SetDirectory(0)
        h_BGOfiducial_psd_stk_match_cut_tmp.SetDirectory(0)
        h_BGOfiducial_psd_charge_cut_tmp.SetDirectory(0)
        h_BGOfiducial_stk_charge_cut_tmp.SetDirectory(0)
        h_BGOfiducial_xtrl_cut_tmp.SetDirectory(0)
        h_BGOfiducial_all_cut_tmp.SetDirectory(0)
        
        h_geo_BGOrec_slopeX_tmp.SetDirectory(0)
        h_geo_BGOrec_slopeY_tmp.SetDirectory(0)
        h_geo_BGOrec_interceptX_tmp.SetDirectory(0)
        h_geo_BGOrec_interceptY_tmp.SetDirectory(0)
        h_geo_BGOreco_topMap_tmp.SetDirectory(0)
        h_geo_BGOreco_bottomMap_tmp.SetDirectory(0)
        
        h_BGOrec_energy_tmp.SetDirectory(0)
        h_layer_max_energy_ratio_tmp.SetDirectory(0)
        for idx in range(0,14):
            h_layer_energy_ratio_tmp[idx].SetDirectory(0)

        for idx in range(n_energy_bins):
            h_sumRms_cosine_tmp[idx].SetDirectory(0)
        sumRms_cosine_20_100_tmp.SetDirectory(0)
        sumRms_cosine_100_250_tmp.SetDirectory(0)
        sumRms_cosine_250_500_tmp.SetDirectory(0)
        sumRms_cosine_500_1000_tmp.SetDirectory(0)
        sumRms_cosine_1000_3000_tmp.SetDirectory(0)
        sumRms_cosine_3000_10000_tmp.SetDirectory(0)

        h_xtrl_energy_int_tmp.SetDirectory(0)
        h_xtrl_tmp.SetDirectory(0)
        e_discrimination_last_tmp.SetDirectory(0)
        e_discrimination_last_20_100_tmp.SetDirectory(0)
        e_discrimination_last_100_250_tmp.SetDirectory(0)
        e_discrimination_last_250_500_tmp.SetDirectory(0)
        e_discrimination_last_500_1000_tmp.SetDirectory(0)
        e_discrimination_last_1000_3000_tmp.SetDirectory(0)
        e_discrimination_last_3000_10000_tmp.SetDirectory(0)

        e_discrimination_tmp.SetDirectory(0)
        e_discrimination_20_100_tmp.SetDirectory(0)
        e_discrimination_100_250_tmp.SetDirectory(0)
        e_discrimination_250_500_tmp.SetDirectory(0)
        e_discrimination_500_1000_tmp.SetDirectory(0)
        e_discrimination_1000_3000_tmp.SetDirectory(0)
        e_discrimination_3000_10000_tmp.SetDirectory(0)

        for idx in range(n_energy_bins):
            h_bin_xtrl_tmp[idx].SetDirectory(0)

        h_psd_chargeX_tmp.SetDirectory(0)
        h_psd_chargeY_tmp.SetDirectory(0)
        h_psd_charge_tmp.SetDirectory(0)
        h_psd_charge2D_tmp.SetDirectory(0)

        h_psd_selected_chargeX_tmp.SetDirectory(0)
        h_psd_selected_chargeY_tmp.SetDirectory(0)
        h_psd_selected_charge_tmp.SetDirectory(0)
        h_psd_selected_charge2D_tmp.SetDirectory(0)

        h_stk_chargeX_tmp.SetDirectory(0)
        h_stk_chargeY_tmp.SetDirectory(0)
        h_stk_charge_tmp.SetDirectory(0)
        h_stk_charge2D_tmp.SetDirectory(0)

        h_stk_selected_chargeX_tmp.SetDirectory(0)
        h_stk_selected_chargeY_tmp.SetDirectory(0)
        h_stk_selected_charge_tmp.SetDirectory(0)
        h_stk_selected_charge2D_tmp.SetDirectory(0)

        h_background_under_xtrl_cut_tmp.SetDirectory(0)
        h_background_over_xtrl_cut_tmp.SetDirectory(0)
        
        # Clone output file
        rFile.Close()

        # Add histos
        if dIdx == 0:
            
            h_trigger = h_trigger_tmp.Clone("h_trigger")
            h_geometric_cut = h_geometric_cut_tmp.Clone("h_geometric_cut")
            h_maxElayer_cut = h_maxElayer_cut_tmp.Clone("h_maxElayer_cut")
            h_maxBarLayer_cut = h_maxBarLayer_cut_tmp.Clone("h_maxBarLayer_cut")
            h_BGOTrackContainment_cut = h_BGOTrackContainment_cut_tmp.Clone("h_BGOTrackContainment_cut")
            h_BGO_fiducial_cut = h_BGO_fiducial_cut_tmp.Clone("h_BGO_fiducial_cut")
            h_all_cut = h_all_cut_tmp.Clone("h_all_cut")
            
            h_geometric_maxElayer_cut = h_geometric_maxElayer_cut_tmp.Clone("h_geometric_maxElayer_cut")
            h_geometric_maxBarLayer_cut = h_geometric_maxBarLayer_cut_tmp.Clone("h_geometric_maxBarLayer_cut")
            h_geometric_BGOTrackContainment_cut = h_geometric_BGOTrackContainment_cut_tmp.Clone("h_geometric_BGOTrackContainment_cut")
            h_geometric_BGO_fiducial_cut = h_geometric_BGO_fiducial_cut_tmp.Clone("h_geometric_BGO_fiducial_cut")
            h_geometric_all_cut = h_geometric_all_cut_tmp.Clone("h_geometric_all_cut")

            h_BGOfiducial_nBarLayer13_cut = h_BGOfiducial_nBarLayer13_cut_tmp.Clone("h_BGOfiducial_nBarLayer13_cut")
            h_BGOfiducial_maxRms_cut = h_BGOfiducial_maxRms_cut_tmp.Clone("h_BGOfiducial_maxRms_cut")
            h_BGOfiducial_track_selection_cut = h_BGOfiducial_track_selection_cut_tmp.Clone("h_BGOfiducial_track_selection_cut")
            h_BGOfiducial_psd_stk_match_cut = h_BGOfiducial_psd_stk_match_cut_tmp.Clone("h_BGOfiducial_psd_stk_match_cut")
            h_BGOfiducial_psd_charge_cut = h_BGOfiducial_psd_charge_cut_tmp.Clone("h_BGOfiducial_psd_charge_cut")
            h_BGOfiducial_stk_charge_cut = h_BGOfiducial_stk_charge_cut_tmp.Clone("h_BGOfiducial_stk_charge_cut")
            h_BGOfiducial_xtrl_cut = h_BGOfiducial_xtrl_cut_tmp.Clone("h_BGOfiducial_xtrl_cut")
            h_BGOfiducial_all_cut = h_BGOfiducial_all_cut_tmp.Clone("h_BGOfiducial_all_cut")
            
            h_geo_BGOrec_slopeX = h_geo_BGOrec_slopeX_tmp.Clone("h_geo_BGOrec_slopeX")
            h_geo_BGOrec_slopeY = h_geo_BGOrec_slopeY_tmp.Clone("h_geo_BGOrec_slopeY")
            h_geo_BGOrec_interceptX = h_geo_BGOrec_interceptX_tmp.Clone("h_geo_BGOrec_interceptX")
            h_geo_BGOrec_interceptY = h_geo_BGOrec_interceptY_tmp.Clone("h_geo_BGOrec_interceptY")
            h_geo_BGOreco_topMap = h_geo_BGOreco_topMap_tmp.Clone("h_geo_BGOreco_topMap")
            h_geo_BGOreco_bottomMap = h_geo_BGOreco_bottomMap_tmp.Clone("h_geo_BGOreco_bottomMap")
            
            h_BGOrec_energy = h_BGOrec_energy_tmp.Clone("h_BGOrec_energy")
            h_layer_max_energy_ratio = h_layer_max_energy_ratio_tmp.Clone("h_layer_max_energy_ratio")
            for idx in range(0,14):
                h_ratio_name = "h_layer_energy_ratio_"
                h_ratio_name += str(idx)
                h_layer_energy_ratio[idx] = h_layer_energy_ratio_tmp[idx].Clone(h_ratio_name)

            for idx in range(n_energy_bins):
                histoName = "sumRms_cosine_"
                histoName += str(idx)
                h_sumRms_cosine[idx] = h_sumRms_cosine_tmp.Clone(histoName)
            sumRms_cosine_20_100 = sumRms_cosine_20_100_tmp.Clone("sumRms_cosine_20_100")
            sumRms_cosine_100_250 = sumRms_cosine_100_250_tmp.Clone("sumRms_cosine_100_250")
            sumRms_cosine_250_500 = sumRms_cosine_250_500_tmp.Clone("sumRms_cosine_250_500")
            sumRms_cosine_500_1000 = sumRms_cosine_500_1000_tmp.Clone("sumRms_cosine_500_1000")
            sumRms_cosine_1000_3000 = sumRms_cosine_1000_3000_tmp.Clone("sumRms_cosine_1000_3000")
            sumRms_cosine_3000_10000 = sumRms_cosine_3000_10000_tmp.Clone("sumRms_cosine_3000_10000")

            h_xtrl_energy_int = h_xtrl_energy_int_tmp.Clone("h_xtrl_energy_int")
            h_xtrl = h_xtrl_tmp.Clone("h_xtrl")
            e_discrimination_last = e_discrimination_last_tmp.Clone("e_discrimination_last")
            e_discrimination_last_20_100 = e_discrimination_last_20_100_tmp.Clone("e_discrimination_last_20_100")
            e_discrimination_last_100_250 = e_discrimination_last_100_250_tmp.Clone("e_discrimination_last_100_250")
            e_discrimination_last_250_500 = e_discrimination_last_250_500_tmp.Clone("e_discrimination_last_250_500")
            e_discrimination_last_500_1000 = e_discrimination_last_500_1000_tmp.Clone("e_discrimination_last_500_1000")
            e_discrimination_last_1000_3000 = e_discrimination_last_1000_3000_tmp.Clone("e_discrimination_last_1000_3000")
            e_discrimination_last_3000_10000 = e_discrimination_last_3000_10000_tmp.Clone("e_discrimination_last_3000_10000")

            e_discrimination = e_discrimination_tmp.Clone("e_discrimination")
            e_discrimination_20_100 = e_discrimination_20_100_tmp.Clone("e_discrimination_20_100")
            e_discrimination_100_250 = e_discrimination_100_250_tmp.Clone("e_discrimination_100_250")
            e_discrimination_250_500 = e_discrimination_250_500_tmp.Clone("e_discrimination_250_500")
            e_discrimination_500_1000 = e_discrimination_500_1000_tmp.Clone("e_discrimination_500_1000")
            e_discrimination_1000_3000 = e_discrimination_1000_3000_tmp.Clone("e_discrimination_1000_3000")
            e_discrimination_3000_10000 = e_discrimination_3000_10000_tmp.Clone("e_discrimination_3000_10000")

            for idx in range(n_energy_bins):
                histoName = "h_xtrl_bin_"
                histoName += str(idx)
                h_bin_xtrl[idx] = h_bin_xtrl_tmp[idx].Clone(histoName)
            h_xtrl = h_xtrl_tmp.Clone("h_xtrl")

            h_psd_chargeX = h_psd_chargeX_tmp.Clone("h_psd_chargeX")
            h_psd_chargeY = h_psd_chargeY_tmp.Clone("h_psd_chargeY")
            h_psd_charge = h_psd_charge_tmp.Clone("h_psd_charge")
            h_psd_charge2D = h_psd_charge2D_tmp.Clone("h_psd_charge2D")

            h_psd_selected_chargeX = h_psd_selected_chargeX_tmp.Clone("h_psd_selected_chargeX")
            h_psd_selected_chargeY = h_psd_selected_chargeY_tmp.Clone("h_psd_selected_chargeY")
            h_psd_selected_charge = h_psd_selected_charge_tmp.Clone("h_psd_selected_charge")
            h_psd_selected_charge2D = h_psd_selected_charge2D_tmp.Clone("h_psd_selected_charge2D")

            h_stk_chargeX = h_stk_chargeX_tmp.Clone("h_stk_chargeX")
            h_stk_chargeY = h_stk_chargeY_tmp.Clone("h_stk_chargeY")
            h_stk_charge = h_stk_charge_tmp.Clone("h_stk_charge")
            h_stk_charge2D = h_stk_charge2D_tmp.Clone("h_stk_charge2D")

            h_stk_selected_chargeX = h_stk_selected_chargeX_tmp.Clone("h_stk_selected_chargeX")
            h_stk_selected_chargeY = h_stk_selected_chargeY_tmp.Clone("h_stk_selected_chargeY")
            h_stk_selected_charge = h_stk_selected_charge_tmp.Clone("h_stk_selected_charge")
            h_stk_selected_charge2D = h_stk_selected_charge2D_tmp.Clone("h_stk_selected_charge2D")

            h_background_under_xtrl_cut = h_background_under_xtrl_cut_tmp.Clone("h_background_under_xtrl_cut")
            h_background_over_xtrl_cut = h_background_over_xtrl_cut_tmp.Clone("h_background_over_xtrl_cut")

        else:
            
            h_trigger.Add(h_trigger_tmp)
            h_geometric_cut.Add(h_geometric_cut_tmp)
            h_maxElayer_cut.Add(h_maxElayer_cut_tmp)
            h_maxBarLayer_cut.Add(h_maxBarLayer_cut_tmp)
            h_BGOTrackContainment_cut.Add(h_BGOTrackContainment_cut_tmp)
            h_BGO_fiducial_cut.Add(h_BGO_fiducial_cut_tmp)
            h_all_cut.Add(h_all_cut_tmp)

            h_geometric_maxElayer_cut.Add(h_geometric_maxElayer_cut_tmp)
            h_geometric_maxBarLayer_cut.Add(h_geometric_maxBarLayer_cut_tmp)
            h_geometric_BGOTrackContainment_cut.Add(h_geometric_BGOTrackContainment_cut_tmp)
            h_geometric_BGO_fiducial_cut.Add(h_geometric_BGO_fiducial_cut_tmp)
            h_geometric_all_cut.Add(h_geometric_all_cut_tmp)

            h_BGOfiducial_nBarLayer13_cut.Add(h_BGOfiducial_nBarLayer13_cut_tmp)
            h_BGOfiducial_maxRms_cut.Add(h_BGOfiducial_maxRms_cut_tmp)
            h_BGOfiducial_track_selection_cut.Add(h_BGOfiducial_track_selection_cut_tmp)
            h_BGOfiducial_psd_stk_match_cut.Add(h_BGOfiducial_psd_stk_match_cut_tmp)
            h_BGOfiducial_psd_charge_cut.Add(h_BGOfiducial_psd_charge_cut_tmp)
            h_BGOfiducial_stk_charge_cut.Add(h_BGOfiducial_stk_charge_cut_tmp)
            h_BGOfiducial_xtrl_cut.Add(h_BGOfiducial_xtrl_cut_tmp)
            h_BGOfiducial_all_cut.Add(h_BGOfiducial_all_cut_tmp)
            
            h_geo_BGOrec_slopeX.Add(h_geo_BGOrec_slopeX_tmp)
            h_geo_BGOrec_slopeY.Add(h_geo_BGOrec_slopeY_tmp)
            h_geo_BGOrec_interceptX.Add(h_geo_BGOrec_interceptX_tmp)
            h_geo_BGOrec_interceptY.Add(h_geo_BGOrec_interceptY_tmp)
            h_geo_BGOreco_topMap.Add(h_geo_BGOreco_topMap_tmp)
            h_geo_BGOreco_bottomMap.Add(h_geo_BGOreco_bottomMap_tmp)
            
            h_BGOrec_energy.Add(h_BGOrec_energy_tmp)
            h_layer_max_energy_ratio.Add(h_layer_max_energy_ratio_tmp)
            for idx in range(0,14):
                h_layer_energy_ratio[idx].Add(h_layer_energy_ratio_tmp[idx])

            for idx in range(n_energy_bins):
                h_sumRms_cosine[idx].Add(h_sumRms_cosine_tmp[idx])
            sumRms_cosine_20_100.Add(sumRms_cosine_20_100_tmp)
            sumRms_cosine_100_250.Add(sumRms_cosine_100_250_tmp)
            sumRms_cosine_250_500.Add(sumRms_cosine_250_500_tmp)
            sumRms_cosine_500_1000.Add(sumRms_cosine_500_1000_tmp)
            sumRms_cosine_1000_3000.Add(sumRms_cosine_1000_3000_tmp)
            sumRms_cosine_3000_10000.Add(sumRms_cosine_3000_10000_tmp)

            h_xtrl_energy_int.Add(h_xtrl_energy_int_tmp)
            h_xtrl.Add(h_xtrl_tmp)
            e_discrimination_last.Add(e_discrimination_last_tmp)
            e_discrimination_last_20_100.Add(e_discrimination_last_20_100_tmp)
            e_discrimination_last_100_250.Add(e_discrimination_last_100_250_tmp)
            e_discrimination_last_250_500.Add(e_discrimination_last_250_500_tmp)
            e_discrimination_last_500_1000.Add(e_discrimination_last_500_1000_tmp)
            e_discrimination_last_1000_3000.Add(e_discrimination_last_1000_3000_tmp)
            e_discrimination_last_3000_10000.Add(e_discrimination_last_3000_10000_tmp)

            e_discrimination.Add(e_discrimination_tmp)
            e_discrimination_20_100.Add(e_discrimination_20_100_tmp)
            e_discrimination_100_250.Add(e_discrimination_100_250_tmp)
            e_discrimination_250_500.Add(e_discrimination_250_500_tmp)
            e_discrimination_500_1000.Add(e_discrimination_500_1000_tmp)
            e_discrimination_1000_3000.Add(e_discrimination_1000_3000_tmp)
            e_discrimination_3000_10000.Add(e_discrimination_3000_10000_tmp)

            for idx in range(n_energy_bins):
                h_bin_xtrl[idx].Add(h_bin_xtrl_tmp[idx])
    
            h_psd_chargeX.Add(h_psd_chargeX_tmp)
            h_psd_chargeY.Add(h_psd_chargeY_tmp)
            h_psd_charge.Add(h_psd_charge_tmp)
            h_psd_charge2D.Add(h_psd_charge2D_tmp)

            h_psd_selected_chargeX.Add(h_psd_selected_chargeX_tmp)
            h_psd_selected_chargeY.Add(h_psd_selected_chargeY_tmp)
            h_psd_selected_charge.Add(h_psd_selected_charge_tmp)
            h_psd_selected_charge2D.Add(h_psd_selected_charge2D_tmp)

            h_stk_chargeX.Add(h_stk_chargeX_tmp)
            h_stk_chargeY.Add(h_stk_chargeY_tmp)
            h_stk_charge.Add(h_stk_charge_tmp)
            h_stk_charge2D.Add(h_stk_charge2D_tmp)

            h_stk_selected_chargeX.Add(h_stk_selected_chargeX_tmp)
            h_stk_selected_chargeY.Add(h_stk_selected_chargeY_tmp)
            h_stk_selected_charge.Add(h_stk_selected_charge_tmp)
            h_stk_selected_charge2D.Add(h_stk_selected_charge2D_tmp)

            h_background_under_xtrl_cut.Add(h_background_under_xtrl_cut_tmp)
            h_background_over_xtrl_cut.Add(h_background_over_xtrl_cut_tmp)

    # Create output file for full histos
    fOut = TFile.Open(opts.output, "RECREATE")
    if fOut.IsOpen():
        if opts.verbose:
            print('Output TFile has been created: {}'.format(opts.output))
    else:
        print('Error creating output TFile: {}'.format(opts.output))
        sys.exit()

    # Writing final histos to file
    h_trigger.Write()
    h_geometric_cut.Write()
    h_maxElayer_cut.Write()
    h_maxBarLayer_cut.Write()
    h_BGOTrackContainment_cut.Write()
    h_BGO_fiducial_cut.Write()
    h_all_cut.Write()
   
    h_geometric_maxElayer_cut.Write()
    h_geometric_maxBarLayer_cut.Write()
    h_geometric_BGOTrackContainment_cut.Write()
    h_geometric_BGO_fiducial_cut.Write()
    h_geometric_all_cut.Write()

    h_BGOfiducial_nBarLayer13_cut.Write()
    h_BGOfiducial_maxRms_cut.Write()
    h_BGOfiducial_track_selection_cut.Write()
    h_BGOfiducial_psd_stk_match_cut.Write()
    h_BGOfiducial_psd_charge_cut.Write()
    h_BGOfiducial_stk_charge_cut.Write()
    h_BGOfiducial_xtrl_cut.Write()
    h_BGOfiducial_all_cut.Write()
    
    fOut.mkdir("Analysis_GeoCut")
    fOut.cd("Analysis_GeoCut")
    
    h_geo_BGOrec_slopeX.Write()
    h_geo_BGOrec_slopeY.Write()
    h_geo_BGOrec_interceptX.Write()
    h_geo_BGOrec_interceptY.Write()
    h_geo_BGOreco_topMap.Write()
    h_geo_BGOreco_bottomMap.Write()
    
    fOut.mkdir("BGO_Energy")
    fOut.cd("BGO_Energy")

    h_BGOrec_energy.Write()
    h_layer_max_energy_ratio.Write()
    for idx in range(0,14):
        h_layer_energy_ratio[idx].Write()
   
    fOut.mkdir("xtrl")
    fOut.cd("xtrl")

    h_xtrl_energy_int.Write()
    h_xtrl.Write()   
    e_discrimination_last.Write()
    e_discrimination_last_20_100.Write()
    e_discrimination_last_100_250.Write()
    e_discrimination_last_250_500.Write()
    e_discrimination_last_500_1000.Write()
    e_discrimination_last_1000_3000.Write()
    e_discrimination_last_3000_10000.Write()
    e_discrimination.Write()
    e_discrimination_20_100.Write()
    e_discrimination_100_250.Write()
    e_discrimination_250_500.Write()
    e_discrimination_500_1000.Write()
    e_discrimination_1000_3000.Write()
    e_discrimination_3000_10000.Write()

    for idx in range(n_energy_bins):
        h_bin_xtrl[idx].Write()

    for idx in range(n_energy_bins):
        h_sumRms_cosine[idx].Write()
    
    sumRms_cosine_20_100.Write()
    sumRms_cosine_100_250.Write()
    sumRms_cosine_250_500.Write()
    sumRms_cosine_500_1000.Write()
    sumRms_cosine_1000_3000.Write()
    sumRms_cosine_3000_10000.Write()

    fOut.mkdir("PSDcharge")
    fOut.cd("PSDcharge")

    h_psd_chargeX.Write()
    h_psd_chargeY.Write()
    h_psd_charge.Write()
    h_psd_charge2D.Write()

    h_psd_selected_chargeX.Write()
    h_psd_selected_chargeY.Write()
    h_psd_selected_charge.Write()
    h_psd_selected_charge2D.Write()

    fOut.mkdir("STKcharge")
    fOut.cd("STKcharge")

    h_stk_chargeX.Write()
    h_stk_chargeY.Write()
    h_stk_charge.Write()
    h_stk_charge2D.Write()

    h_stk_selected_chargeX.Write()
    h_stk_selected_chargeY.Write()
    h_stk_selected_charge.Write()
    h_stk_selected_charge2D.Write()

    fOut.mkdir("mc_ancillary")
    fOut.cd("mc_ancillary")

    h_background_under_xtrl_cut.Write()
    h_background_over_xtrl_cut.Write()
    
    # Closing output file
    fOut.Close()