from ROOT import TFile, TH1D, TH2D
import os
import sys


class mc_histos():

    def __init__(self):

        # Cut histos
        self.h_geo_factor = TH1D()
        self.h_incoming = TH1D()
        self.h_trigger = TH1D()
        self.h_geometric_cut = TH1D()
        self.h_maxElayer_cut = TH1D()
        self.h_maxBarLayer_cut = TH1D()
        self.h_BGOTrackContainment_cut = TH1D()
        self.h_BGO_fiducial_cut = TH1D()
        self.h_all_cut = TH1D()

        self.h_geo_factor_nw = TH1D()
        self.h_incoming_nw = TH1D()
        self.h_trigger_nw = TH1D()
        self.h_geometric_cut_nw = TH1D()
        self.h_maxElayer_cut_nw = TH1D()
        self.h_maxBarLayer_cut_nw = TH1D()
        self.h_BGOTrackContainment_cut_nw = TH1D()
        self.h_BGO_fiducial_cut_nw = TH1D()
        self.h_all_cut_nw = TH1D()

        # Cuts && Geometric Cut
        self.h_geometric_maxElayer_cut = TH1D()
        self.h_geometric_maxBarLayer_cut = TH1D()
        self.h_geometric_BGOTrackContainment_cut = TH1D()
        self.h_geometric_BGO_fiducial_cut = TH1D()
        self.h_geometric_all_cut = TH1D()

        self.h_geometric_maxElayer_cut_nw = TH1D()
        self.h_geometric_maxBarLayer_cut_nw = TH1D()
        self.h_geometric_BGOTrackContainment_cut_nw = TH1D()
        self.h_geometric_BGO_fiducial_cut_nw = TH1D()
        self.h_geometric_all_cut_nw = TH1D()

        # Cuts && BGO fiducial volume cut
        self.h_BGOfiducial_nBarLayer13_cut = TH1D()
        self.h_BGOfiducial_maxRms_cut = TH1D()
        self.h_BGOfiducial_track_selection_cut = TH1D()
        self.h_BGOfiducial_psd_stk_match_cut = TH1D()
        self.h_BGOfiducial_psd_charge_cut = TH1D()
        self.h_BGOfiducial_stk_charge_cut = TH1D()
        self.h_BGOfiducial_xtrl_cut = TH1D()
        self.h_BGOfiducial_all_cut = TH1D()

        self.h_BGOfiducial_nBarLayer13_cut_nw = TH1D()
        self.h_BGOfiducial_maxRms_cut_nw = TH1D()
        self.h_BGOfiducial_track_selection_cut_nw = TH1D()
        self.h_BGOfiducial_psd_stk_match_cut_nw = TH1D()
        self.h_BGOfiducial_psd_charge_cut_nw = TH1D()
        self.h_BGOfiducial_stk_charge_cut_nw = TH1D()
        self.h_BGOfiducial_xtrl_cut_nw = TH1D()
        self.h_BGOfiducial_all_cut_nw = TH1D()

        # Analysis histos - simu and reco energy of incoming events
        self.h_BGOrec_energy = TH1D()
        self.h_simu_energy = TH1D()
        self.h_energy_diff = TH1D()
        self.h_energy_diff2D = TH2D()
        self.h_energy_unfold = TH2D()

        # Ratio of layer energy respect to total BGO energy
        self.h_layer_max_energy_ratio = TH1D()

        self.h_layer_energy_ratio = []
        for idx in range(0, 14):
            self.h_layer_energy_ratio.append(TH1D())

        # Pre Geometric Cut
        # Top X and Y
        self.h_preGeo_BGOrec_topX_vs_realX = TH1D()
        self.h_preGeo_BGOrec_topY_vs_realY = TH1D()
        # Slope X and Y
        self.h_preGeo_real_slopeX = TH1D()
        self.h_preGeo_real_slopeY = TH1D()
        self.h_preGeo_BGOrec_slopeX = TH1D()
        self.h_preGeo_BGOrec_slopeY = TH1D()
        # Intercept X and Y
        self.h_preGeo_real_interceptX = TH1D()
        self.h_preGeo_real_interceptY = TH1D()
        self.h_preGeo_BGOrec_interceptX = TH1D()
        self.h_preGeo_BGOrec_interceptY = TH1D()
        # Top Maps
        self.h_preGeo_real_topMap = TH2D()
        self.h_preGeo_BGOreco_topMap = TH2D()
        # Bottom Maps
        self.h_preGeo_real_bottomMap = TH2D()
        self.h_preGeo_BGOreco_bottomMap = TH2D()
        # Map of events outside the "real" first BGO layer
        self.h_noBGOenergy_real_topMap = TH2D()

        # After Geometric Cut
        # Top X and Y
        self.h_geo_BGOrec_topX_vs_realX = TH1D()
        self.h_geo_BGOrec_topY_vs_realY = TH1D()
        # Slope X and Y
        self.h_geo_real_slopeX = TH1D()
        self.h_geo_real_slopeY = TH1D()
        self.h_geo_BGOrec_slopeX = TH1D()
        self.h_geo_BGOrec_slopeY = TH1D()
        # Intercept X and Y
        self.h_geo_real_interceptX = TH1D()
        self.h_geo_real_interceptY = TH1D()
        self.h_geo_BGOrec_interceptX = TH1D()
        self.h_geo_BGOrec_interceptY = TH1D()
        # Top Maps
        self.h_geo_real_topMap = TH2D()
        self.h_geo_BGOreco_topMap = TH2D()
        # Bottom Maps
        self.h_geo_real_bottomMap = TH2D()
        self.h_geo_BGOreco_bottomMap = TH2D()

        # sumRms - cosine correlation
        self.h_sumRms_cosine = []
        self.sumRms_cosine_20_100 = TH2D()
        self.sumRms_cosine_100_250 = TH2D()
        self.sumRms_cosine_250_500 = TH2D()
        self.sumRms_cosine_500_1000 = TH2D()
        self.sumRms_cosine_1000_3000 = TH2D()
        self.sumRms_cosine_3000_10000 = TH2D()

        # XTRL histos
        self.h_xtrl_energy_int = TH1D()
        self.h_xtrl = TH2D()
        self.e_discrimination_last = TH2D()
        self.e_discrimination_last_20_100 = TH2D()
        self.e_discrimination_last_100_250 = TH2D()
        self.e_discrimination_last_250_500 = TH2D()
        self.e_discrimination_last_500_1000 = TH2D()
        self.e_discrimination_last_1000_3000 = TH2D()
        self.e_discrimination_last_3000_10000 = TH2D()
        self.e_discrimination = TH2D()
        self.e_discrimination_20_100 = TH2D()
        self.e_discrimination_100_250 = TH2D()
        self.e_discrimination_250_500 = TH2D()
        self.e_discrimination_500_1000 = TH2D()
        self.e_discrimination_1000_3000 = TH2D()
        self.e_discrimination_3000_10000 = TH2D()
        self.h_bin_xtrl = []

        # PSD charge histos
        self.h_psd_chargeX = TH1D()
        self.h_psd_chargeY = TH1D()
        self.h_psd_charge = TH1D()
        self.h_psd_charge2D = TH2D()

        self.h_psd_selected_chargeX = TH1D()
        self.h_psd_selected_chargeY = TH1D()
        self.h_psd_selected_charge = TH1D()
        self.h_psd_selected_charge2D = TH2D()

        # STK charge histos
        self.h_stk_chargeX = TH1D()
        self.h_stk_chargeY = TH1D()
        self.h_stk_charge = TH1D()
        self.h_stk_charge2D = TH2D()

        self.h_stk_selected_chargeX = TH1D()
        self.h_stk_selected_chargeY = TH1D()
        self.h_stk_selected_charge = TH1D()
        self.h_stk_selected_charge2D = TH2D()

        self.first_file_read = True
        self.n_energy_bins = 0

    def add_file(self, rFile):
        
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

        if self.first_file_read:
            self.n_energy_bins = h_geo_factor_tmp.GetNbinsX()
            # Init bin_xtrl and sumRms
            for idx in range(self.n_energy_bins):
                self.h_bin_xtrl.append(TH1D())
                self.h_sumRms_cosine.append(TH2D())

        h_geo_factor_tmp_nw = rFile.Get("nw_histos/h_geo_factor_nw")
        h_incoming_tmp_nw = rFile.Get("nw_histos/h_incoming_nw")
        h_trigger_tmp_nw = rFile.Get("nw_histos/h_trigger_nw")
        h_geometric_cut_tmp_nw = rFile.Get("nw_histos/h_geometric_cut_nw")
        h_maxElayer_cut_tmp_nw = rFile.Get("nw_histos/h_maxElayer_cut_nw")
        h_maxBarLayer_cut_tmp_nw = rFile.Get("nw_histos/h_maxBarLayer_cut_nw")
        h_BGOTrackContainment_cut_tmp_nw = rFile.Get(
            "nw_histos/h_BGOTrackContainment_cut_nw")
        h_BGO_fiducial_cut_tmp_nw = rFile.Get(
            "nw_histos/h_BGO_fiducial_cut_nw")
        h_all_cut_tmp_nw = rFile.Get("nw_histos/h_all_cut_nw")

        h_geometric_maxElayer_cut_tmp = rFile.Get("h_geometric_maxElayer_cut")
        h_geometric_maxBarLayer_cut_tmp = rFile.Get(
            "h_geometric_maxBarLayer_cut")
        h_geometric_BGOTrackContainment_cut_tmp = rFile.Get(
            "h_geometric_BGOTrackContainment_cut")
        h_geometric_BGO_fiducial_cut_tmp = rFile.Get(
            "h_geometric_BGO_fiducial_cut")
        h_geometric_all_cut_tmp = rFile.Get("h_geometric_all_cut")

        h_geometric_maxElayer_cut_tmp_nw = rFile.Get(
            "nw_histos/h_geometric_maxElayer_cut_nw")
        h_geometric_maxBarLayer_cut_tmp_nw = rFile.Get(
            "nw_histos/h_geometric_maxBarLayer_cut_nw")
        h_geometric_BGOTrackContainment_cut_tmp_nw = rFile.Get(
            "nw_histos/h_geometric_BGOTrackContainment_cut_nw")
        h_geometric_BGO_fiducial_cut_tmp_nw = rFile.Get(
            "nw_histos/h_geometric_BGO_fiducial_cut_nw")
        h_geometric_all_cut_tmp_nw = rFile.Get(
            "nw_histos/h_geometric_all_cut_nw")

        h_BGOfiducial_nBarLayer13_cut_tmp = rFile.Get(
            "h_BGOfiducial_nBarLayer13_cut")
        h_BGOfiducial_maxRms_cut_tmp = rFile.Get("h_BGOfiducial_maxRms_cut")
        h_BGOfiducial_track_selection_cut_tmp = rFile.Get(
            "h_BGOfiducial_track_selection_cut")
        h_BGOfiducial_psd_stk_match_cut_tmp = rFile.Get(
            "h_BGOfiducial_psd_stk_match_cut")
        h_BGOfiducial_psd_charge_cut_tmp = rFile.Get(
            "h_BGOfiducial_psd_charge_cut")
        h_BGOfiducial_stk_charge_cut_tmp = rFile.Get(
            "h_BGOfiducial_stk_charge_cut")
        h_BGOfiducial_xtrl_cut_tmp = rFile.Get("h_BGOfiducial_xtrl_cut")
        h_BGOfiducial_all_cut_tmp = rFile.Get("h_BGOfiducial_all_cut")

        h_BGOfiducial_nBarLayer13_cut_tmp_nw = rFile.Get(
            "nw_histos/h_BGOfiducial_nBarLayer13_cut_nw")
        h_BGOfiducial_maxRms_cut_tmp_nw = rFile.Get(
            "nw_histos/h_BGOfiducial_maxRms_cut_nw")
        h_BGOfiducial_track_selection_cut_tmp_nw = rFile.Get(
            "nw_histos/h_BGOfiducial_track_selection_cut_nw")
        h_BGOfiducial_psd_stk_match_cut_tmp_nw = rFile.Get(
            "nw_histos/h_BGOfiducial_psd_stk_match_cut_nw")
        h_BGOfiducial_psd_charge_cut_tmp_nw = rFile.Get(
            "nw_histos/h_BGOfiducial_psd_charge_cut_nw")
        h_BGOfiducial_stk_charge_cut_tmp_nw = rFile.Get(
            "nw_histos/h_BGOfiducial_stk_charge_cut_nw")
        h_BGOfiducial_xtrl_cut_tmp_nw = rFile.Get(
            "nw_histos/h_BGOfiducial_xtrl_cut_nw")
        h_BGOfiducial_all_cut_tmp_nw = rFile.Get(
            "nw_histos/h_BGOfiducial_all_cut_nw")

        h_preGeo_BGOrec_topX_vs_realX_tmp = rFile.Get(
            "Analysis_preGeoCut/h_preGeo_BGOrec_topX_vs_realX")
        h_preGeo_BGOrec_topY_vs_realY_tmp = rFile.Get(
            "Analysis_preGeoCut/h_preGeo_BGOrec_topY_vs_realY")
        h_preGeo_real_slopeX_tmp = rFile.Get(
            "Analysis_preGeoCut/h_preGeo_real_slopeX")
        h_preGeo_real_slopeY_tmp = rFile.Get(
            "Analysis_preGeoCut/h_preGeo_real_slopeY")
        h_preGeo_BGOrec_slopeX_tmp = rFile.Get(
            "Analysis_preGeoCut/h_preGeo_BGOrec_slopeX")
        h_preGeo_BGOrec_slopeY_tmp = rFile.Get(
            "Analysis_preGeoCut/h_preGeo_BGOrec_slopeY")
        h_preGeo_real_interceptX_tmp = rFile.Get(
            "Analysis_preGeoCut/h_preGeo_real_interceptX")
        h_preGeo_real_interceptY_tmp = rFile.Get(
            "Analysis_preGeoCut/h_preGeo_real_interceptY")
        h_preGeo_BGOrec_interceptX_tmp = rFile.Get(
            "Analysis_preGeoCut/h_preGeo_BGOrec_interceptX")
        h_preGeo_BGOrec_interceptY_tmp = rFile.Get(
            "Analysis_preGeoCut/h_preGeo_BGOrec_interceptY")
        h_preGeo_real_topMap_tmp = rFile.Get(
            "Analysis_preGeoCut/h_preGeo_real_topMap")
        h_preGeo_BGOreco_topMap_tmp = rFile.Get(
            "Analysis_preGeoCut/h_preGeo_BGOreco_topMap")
        h_preGeo_real_bottomMap_tmp = rFile.Get(
            "Analysis_preGeoCut/h_preGeo_real_bottomMap")
        h_preGeo_BGOreco_bottomMap_tmp = rFile.Get(
            "Analysis_preGeoCut/h_preGeo_BGOreco_bottomMap")
        h_noBGOenergy_real_topMap_tmp = rFile.Get(
            "Analysis_preGeoCut/h_noBGOenergy_real_topMap")

        h_geo_BGOrec_topX_vs_realX_tmp = rFile.Get(
            "Analysis_GeoCut/h_geo_BGOrec_topX_vs_realX")
        h_geo_BGOrec_topY_vs_realY_tmp = rFile.Get(
            "Analysis_GeoCut/h_geo_BGOrec_topY_vs_realY")
        h_geo_real_slopeX_tmp = rFile.Get("Analysis_GeoCut/h_geo_real_slopeX")
        h_geo_real_slopeY_tmp = rFile.Get("Analysis_GeoCut/h_geo_real_slopeY")
        h_geo_BGOrec_slopeX_tmp = rFile.Get(
            "Analysis_GeoCut/h_geo_BGOrec_slopeX")
        h_geo_BGOrec_slopeY_tmp = rFile.Get(
            "Analysis_GeoCut/h_geo_BGOrec_slopeY")
        h_geo_real_interceptX_tmp = rFile.Get(
            "Analysis_GeoCut/h_geo_real_interceptX")
        h_geo_real_interceptY_tmp = rFile.Get(
            "Analysis_GeoCut/h_geo_real_interceptY")
        h_geo_BGOrec_interceptX_tmp = rFile.Get(
            "Analysis_GeoCut/h_geo_BGOrec_interceptX")
        h_geo_BGOrec_interceptY_tmp = rFile.Get(
            "Analysis_GeoCut/h_geo_BGOrec_interceptY")
        h_geo_real_topMap_tmp = rFile.Get("Analysis_GeoCut/h_geo_real_topMap")
        h_geo_BGOreco_topMap_tmp = rFile.Get(
            "Analysis_GeoCut/h_geo_BGOreco_topMap")
        h_geo_real_bottomMap_tmp = rFile.Get(
            "Analysis_GeoCut/h_geo_real_bottomMap")
        h_geo_BGOreco_bottomMap_tmp = rFile.Get(
            "Analysis_GeoCut/h_geo_BGOreco_bottomMap")

        h_BGOrec_energy_tmp = rFile.Get("BGO_Energy/h_BGOrec_energy")
        h_simu_energy_tmp = rFile.Get("BGO_Energy/h_simu_energy")
        h_energy_diff_tmp = rFile.Get("BGO_Energy/h_energy_diff")
        h_energy_diff2D_tmp = rFile.Get("BGO_Energy/h_energy_diff2D")
        h_energy_unfold_tmp = rFile.Get("BGO_Energy/h_energy_unfold")
        h_layer_max_energy_ratio_tmp = rFile.Get(
            "BGO_Energy/h_layer_max_energy_ratio")
        h_layer_energy_ratio_tmp = []
        for idx in range(0, 14):
            histoName = "BGO_Energy/h_layer_energy_ratio_"
            histoName += str(idx)
            h_layer_energy_ratio_tmp.append(rFile.Get(histoName))

        h_sumRms_cosine_tmp = []
        for idx in range(self.n_energy_bins):
            histoName = "BGO_Energy/sumRms_cosine_"
            histoName += str(idx)
            h_sumRms_cosine_tmp.append(rFile.Get(histoName))

        sumRms_cosine_20_100_tmp = rFile.Get("BGO_Energy/sumRms_cosine_20_100")
        sumRms_cosine_100_250_tmp = rFile.Get(
            "BGO_Energy/sumRms_cosine_100_250")
        sumRms_cosine_250_500_tmp = rFile.Get(
            "BGO_Energy/sumRms_cosine_250_500")
        sumRms_cosine_500_1000_tmp = rFile.Get(
            "BGO_Energy/sumRms_cosine_500_1000")
        sumRms_cosine_1000_3000_tmp = rFile.Get(
            "BGO_Energy/sumRms_cosine_1000_3000")
        sumRms_cosine_3000_10000_tmp = rFile.Get(
            "BGO_Energy/sumRms_cosine_3000_10000")

        h_xtrl_energy_int_tmp = rFile.Get("xtrl/h_xtrl_energy_int")
        h_xtrl_tmp = rFile.Get("xtrl/h_xtrl")
        e_discrimination_last_tmp = rFile.Get("xtrl/e_discrimination_last")
        '''
        e_discrimination_last_20_100_tmp = rFile.Get("xtrl/e_discrimination_last_20_100")
        e_discrimination_last_100_250_tmp = rFile.Get("xtrl/e_discrimination_last_100_250")
        e_discrimination_last_250_500_tmp = rFile.Get("xtrl/e_discrimination_last_250_500")
        e_discrimination_last_500_1000_tmp = rFile.Get("xtrl/e_discrimination_last_500_1000")
        e_discrimination_last_1000_3000_tmp = rFile.Get("xtrl/e_discrimination_last_1000_3000")
        e_discrimination_last_3000_10000_tmp = rFile.Get("xtrl/e_discrimination_last_3000_10000")
        '''

        e_discrimination_tmp = rFile.Get("xtrl/e_discrimination")
        '''
        e_discrimination_20_100_tmp = rFile.Get("xtrl/e_discrimination_20_100")
        e_discrimination_100_250_tmp = rFile.Get("xtrl/e_discrimination_100_250")
        e_discrimination_250_500_tmp = rFile.Get("xtrl/e_discrimination_250_500")
        e_discrimination_500_1000_tmp = rFile.Get("xtrl/e_discrimination_500_1000")
        e_discrimination_1000_3000_tmp = rFile.Get("xtrl/e_discrimination_1000_3000")
        e_discrimination_3000_10000_tmp = rFile.Get("xtrl/e_discrimination_3000_10000")
        '''

        h_bin_xtrl_tmp = []
        for idx in range(self.n_energy_bins):
            histoName = "xtrl/h_xtrl_bin_"
            histoName += str(idx)
            h_bin_xtrl_tmp.append(rFile.Get(histoName))

        h_psd_chargeX_tmp = rFile.Get("PSDcharge/h_psd_chargeX")
        h_psd_chargeY_tmp = rFile.Get("PSDcharge/h_psd_chargeY")
        h_psd_charge_tmp = rFile.Get("PSDcharge/h_psd_charge")
        h_psd_charge2D_tmp = rFile.Get("PSDcharge/h_psd_charge2D")

        h_psd_selected_chargeX_tmp = rFile.Get(
            "PSDcharge/h_psd_selected_chargeX")
        h_psd_selected_chargeY_tmp = rFile.Get(
            "PSDcharge/h_psd_selected_chargeY")
        h_psd_selected_charge_tmp = rFile.Get(
            "PSDcharge/h_psd_selected_charge")
        h_psd_selected_charge2D_tmp = rFile.Get(
            "PSDcharge/h_psd_selected_charge2D")

        h_stk_chargeX_tmp = rFile.Get("STKcharge/h_stk_chargeX")
        h_stk_chargeY_tmp = rFile.Get("STKcharge/h_stk_chargeY")
        h_stk_charge_tmp = rFile.Get("STKcharge/h_stk_charge")
        h_stk_charge2D_tmp = rFile.Get("STKcharge/h_stk_charge2D")

        h_stk_selected_chargeX_tmp = rFile.Get(
            "STKcharge/h_stk_selected_chargeX")
        h_stk_selected_chargeY_tmp = rFile.Get(
            "STKcharge/h_stk_selected_chargeY")
        h_stk_selected_charge_tmp = rFile.Get(
            "STKcharge/h_stk_selected_charge")
        h_stk_selected_charge2D_tmp = rFile.Get(
            "STKcharge/h_stk_selected_charge2D")

        h_background_under_xtrl_cut_tmp = rFile.Get(
            "mc_ancillary/h_background_under_xtrl_cut")
        h_background_over_xtrl_cut_tmp = rFile.Get(
            "mc_ancillary/h_background_over_xtrl_cut")

        if self.first_file_read:

            self.h_geo_factor = h_geo_factor_tmp.Clone("h_geo_factor")
            self.h_incoming = h_incoming_tmp.Clone("h_incoming")
            self.h_trigger = h_trigger_tmp.Clone("h_trigger")
            self.h_geometric_cut = h_geometric_cut_tmp.Clone("h_geometric_cut")
            self.h_maxElayer_cut = h_maxElayer_cut_tmp.Clone("h_maxElayer_cut")
            self.h_maxBarLayer_cut = h_maxBarLayer_cut_tmp.Clone(
                "h_maxBarLayer_cut")
            self.h_BGOTrackContainment_cut = h_BGOTrackContainment_cut_tmp.Clone(
                "h_BGOTrackContainment_cut")
            self.h_BGO_fiducial_cut = h_BGO_fiducial_cut_tmp.Clone(
                "h_BGO_fiducial_cut")
            self.h_all_cut = h_all_cut_tmp.Clone("h_all_cut")

            self.h_geo_factor_nw = h_geo_factor_tmp_nw.Clone("h_geo_factor_nw")
            self.h_incoming_nw = h_incoming_tmp_nw.Clone("h_incoming_nw")
            self.h_trigger_nw = h_trigger_tmp_nw.Clone("h_trigger_nw")
            self.h_geometric_cut_nw = h_geometric_cut_tmp_nw.Clone(
                "h_geometric_cut_nw")
            self.h_maxElayer_cut_nw = h_maxElayer_cut_tmp_nw.Clone(
                "h_maxElayer_cut_nw")
            self.h_maxBarLayer_cut_nw = h_maxBarLayer_cut_tmp_nw.Clone(
                "h_maxBarLayer_cut_nw")
            self.h_BGOTrackContainment_cut_nw = h_BGOTrackContainment_cut_tmp_nw.Clone(
                "h_BGOTrackContainment_cut_nw")
            self.h_BGO_fiducial_cut_nw = h_BGO_fiducial_cut_tmp_nw.Clone(
                "h_BGO_fiducial_cut_nw")
            self.h_all_cut_nw = h_all_cut_tmp_nw.Clone("h_all_cut_nw")

            self.h_geometric_maxElayer_cut = h_geometric_maxElayer_cut_tmp.Clone(
                "h_geometric_maxElayer_cut")
            self.h_geometric_maxBarLayer_cut = h_geometric_maxBarLayer_cut_tmp.Clone(
                "h_geometric_maxBarLayer_cut")
            self.h_geometric_BGOTrackContainment_cut = h_geometric_BGOTrackContainment_cut_tmp.Clone(
                "h_geometric_BGOTrackContainment_cut")
            self.h_geometric_BGO_fiducial_cut = h_geometric_BGO_fiducial_cut_tmp.Clone(
                "h_geometric_BGO_fiducial_cut")
            self.h_geometric_all_cut = h_geometric_all_cut_tmp.Clone(
                "h_geometric_all_cut")

            self.h_geometric_maxElayer_cut_nw = h_geometric_maxElayer_cut_tmp_nw.Clone(
                "h_geometric_maxElayer_cut_nw")
            self.h_geometric_maxBarLayer_cut_nw = h_geometric_maxBarLayer_cut_tmp_nw.Clone(
                "h_geometric_maxBarLayer_cut_nw")
            self.h_geometric_BGOTrackContainment_cut_nw = h_geometric_BGOTrackContainment_cut_tmp_nw.Clone(
                "h_geometric_BGOTrackContainment_cut_nw")
            self.h_geometric_BGO_fiducial_cut_nw = h_geometric_BGO_fiducial_cut_tmp_nw.Clone(
                "h_geometric_BGO_fiducial_cut_nw")
            self.h_geometric_all_cut_nw = h_geometric_all_cut_tmp_nw.Clone(
                "h_geometric_all_cut_nw")

            self.h_BGOfiducial_nBarLayer13_cut = h_BGOfiducial_nBarLayer13_cut_tmp.Clone(
                "h_BGOfiducial_nBarLayer13_cut")
            self.h_BGOfiducial_maxRms_cut = h_BGOfiducial_maxRms_cut_tmp.Clone(
                "h_BGOfiducial_maxRms_cut")
            self.h_BGOfiducial_track_selection_cut = h_BGOfiducial_track_selection_cut_tmp.Clone(
                "h_BGOfiducial_track_selection_cut")
            self.h_BGOfiducial_psd_stk_match_cut = h_BGOfiducial_psd_stk_match_cut_tmp.Clone(
                "h_BGOfiducial_psd_stk_match_cut")
            self.h_BGOfiducial_psd_charge_cut = h_BGOfiducial_psd_charge_cut_tmp.Clone(
                "h_BGOfiducial_psd_charge_cut")
            self.h_BGOfiducial_stk_charge_cut = h_BGOfiducial_stk_charge_cut_tmp.Clone(
                "h_BGOfiducial_stk_charge_cut")
            self.h_BGOfiducial_xtrl_cut = h_BGOfiducial_xtrl_cut_tmp.Clone(
                "h_BGOfiducial_xtrl_cut")
            self.h_BGOfiducial_all_cut = h_BGOfiducial_all_cut_tmp.Clone(
                "h_BGOfiducial_all_cut")

            self.h_BGOfiducial_nBarLayer13_cut_nw = h_BGOfiducial_nBarLayer13_cut_tmp_nw.Clone(
                "h_BGOfiducial_nBarLayer13_cut_nw")
            self.h_BGOfiducial_maxRms_cut_nw = h_BGOfiducial_maxRms_cut_tmp_nw.Clone(
                "h_BGOfiducial_maxRms_cut_nw")
            self.h_BGOfiducial_track_selection_cut_nw = h_BGOfiducial_track_selection_cut_tmp_nw.Clone(
                "h_BGOfiducial_track_selection_cut_nw")
            self.h_BGOfiducial_psd_stk_match_cut_nw = h_BGOfiducial_psd_stk_match_cut_tmp_nw.Clone(
                "h_BGOfiducial_psd_stk_match_cut_nw")
            self.h_BGOfiducial_psd_charge_cut_nw = h_BGOfiducial_psd_charge_cut_tmp_nw.Clone(
                "h_BGOfiducial_psd_charge_cut_nw")
            self.h_BGOfiducial_stk_charge_cut_nw = h_BGOfiducial_stk_charge_cut_tmp_nw.Clone(
                "h_BGOfiducial_stk_charge_cut_nw")
            self.h_BGOfiducial_xtrl_cut_nw = h_BGOfiducial_xtrl_cut_tmp_nw.Clone(
                "h_BGOfiducial_xtrl_cut_nw")
            self.h_BGOfiducial_all_cut_nw = h_BGOfiducial_all_cut_tmp_nw.Clone(
                "h_BGOfiducial_all_cut_nw")

            self.h_preGeo_BGOrec_topX_vs_realX = h_preGeo_BGOrec_topX_vs_realX_tmp.Clone(
                "h_preGeo_BGOrec_topX_vs_realX")
            self.h_preGeo_BGOrec_topY_vs_realY = h_preGeo_BGOrec_topY_vs_realY_tmp.Clone(
                "h_preGeo_BGOrec_topY_vs_realY")
            self.h_preGeo_real_slopeX = h_preGeo_real_slopeX_tmp.Clone(
                "h_preGeo_real_slopeX")
            self.h_preGeo_real_slopeY = h_preGeo_real_slopeY_tmp.Clone(
                "h_preGeo_real_slopeY")
            self.h_preGeo_BGOrec_slopeX = h_preGeo_BGOrec_slopeX_tmp.Clone(
                "h_preGeo_BGOrec_slopeX")
            self.h_preGeo_BGOrec_slopeY = h_preGeo_BGOrec_slopeY_tmp.Clone(
                "h_preGeo_BGOrec_slopeY")
            self.h_preGeo_real_interceptX = h_preGeo_real_interceptX_tmp.Clone(
                "h_preGeo_real_interceptX")
            self.h_preGeo_real_interceptY = h_preGeo_real_interceptY_tmp.Clone(
                "h_preGeo_real_interceptY")
            self.h_preGeo_BGOrec_interceptX = h_preGeo_BGOrec_interceptX_tmp.Clone(
                "h_preGeo_BGOrec_interceptX")
            self.h_preGeo_BGOrec_interceptY = h_preGeo_BGOrec_interceptY_tmp.Clone(
                "h_preGeo_BGOrec_interceptY")
            self.h_preGeo_real_topMap = h_preGeo_real_topMap_tmp.Clone(
                "h_preGeo_real_topMap")
            self.h_preGeo_BGOreco_topMap = h_preGeo_BGOreco_topMap_tmp.Clone(
                "h_preGeo_BGOreco_topMap")
            self.h_preGeo_real_bottomMap = h_preGeo_real_bottomMap_tmp.Clone(
                "h_preGeo_real_bottomMap")
            self.h_preGeo_BGOreco_bottomMap = h_preGeo_BGOreco_bottomMap_tmp.Clone(
                "h_preGeo_BGOreco_bottomMap")
            self.h_noBGOenergy_real_topMap = h_noBGOenergy_real_topMap_tmp.Clone(
                "h_noBGOenergy_real_topMap")

            self.h_geo_BGOrec_topX_vs_realX = h_geo_BGOrec_topX_vs_realX_tmp.Clone(
                "h_geo_BGOrec_topX_vs_realX")
            self.h_geo_BGOrec_topY_vs_realY = h_geo_BGOrec_topY_vs_realY_tmp.Clone(
                "h_geo_BGOrec_topY_vs_realY")
            self.h_geo_real_slopeX = h_geo_real_slopeX_tmp.Clone(
                "h_geo_real_slopeX")
            self.h_geo_real_slopeY = h_geo_real_slopeY_tmp.Clone(
                "h_geo_real_slopeY")
            self.h_geo_BGOrec_slopeX = h_geo_BGOrec_slopeX_tmp.Clone(
                "h_geo_BGOrec_slopeX")
            self.h_geo_BGOrec_slopeY = h_geo_BGOrec_slopeY_tmp.Clone(
                "h_geo_BGOrec_slopeY")
            self.h_geo_real_interceptX = h_geo_real_interceptX_tmp.Clone(
                "h_geo_real_interceptX")
            self.h_geo_real_interceptY = h_geo_real_interceptY_tmp.Clone(
                "h_geo_real_interceptY")
            self.h_geo_BGOrec_interceptX = h_geo_BGOrec_interceptX_tmp.Clone(
                "h_geo_BGOrec_interceptX")
            self.h_geo_BGOrec_interceptY = h_geo_BGOrec_interceptY_tmp.Clone(
                "h_geo_BGOrec_interceptY")
            self.h_geo_real_topMap = h_geo_real_topMap_tmp.Clone(
                "h_geo_real_topMap")
            self.h_geo_BGOreco_topMap = h_geo_BGOreco_topMap_tmp.Clone(
                "h_geo_BGOreco_topMap")
            self.h_geo_real_bottomMap = h_geo_real_bottomMap_tmp.Clone(
                "h_geo_real_bottomMap")
            self.h_geo_BGOreco_bottomMap = h_geo_BGOreco_bottomMap_tmp.Clone(
                "h_geo_BGOreco_bottomMap")

            self.h_BGOrec_energy = h_BGOrec_energy_tmp.Clone("h_BGOrec_energy")
            self.h_simu_energy = h_simu_energy_tmp.Clone("h_simu_energy")
            self.h_energy_diff = h_energy_diff_tmp.Clone("h_energy_diff")
            self.h_energy_diff2D = h_energy_diff2D_tmp.Clone("h_energy_diff2D")
            self.h_energy_unfold = h_energy_unfold_tmp.Clone("h_energy_unfold")
            self.h_layer_max_energy_ratio = h_layer_max_energy_ratio_tmp.Clone(
                "h_layer_max_energy_ratio")
            for idx in range(0, 14):
                h_ratio_name = "h_layer_energy_ratio_"
                h_ratio_name += str(idx)
                self.h_layer_energy_ratio[idx] = h_layer_energy_ratio_tmp[idx].Clone(
                    h_ratio_name)

            for idx in range(self.n_energy_bins):
                histoName = "sumRms_cosine_"
                histoName += str(idx)
                self.h_sumRms_cosine[idx] = h_sumRms_cosine_tmp[idx].Clone(
                    histoName)
            self.sumRms_cosine_20_100 = sumRms_cosine_20_100_tmp.Clone(
                "sumRms_cosine_20_100")
            self.sumRms_cosine_100_250 = sumRms_cosine_100_250_tmp.Clone(
                "sumRms_cosine_100_250")
            self.sumRms_cosine_250_500 = sumRms_cosine_250_500_tmp.Clone(
                "sumRms_cosine_250_500")
            self.sumRms_cosine_500_1000 = sumRms_cosine_500_1000_tmp.Clone(
                "sumRms_cosine_500_1000")
            self.sumRms_cosine_1000_3000 = sumRms_cosine_1000_3000_tmp.Clone(
                "sumRms_cosine_1000_3000")
            self.sumRms_cosine_3000_10000 = sumRms_cosine_3000_10000_tmp.Clone(
                "sumRms_cosine_3000_10000")

            self.h_xtrl_energy_int = h_xtrl_energy_int_tmp.Clone(
                "h_xtrl_energy_int")
            self.h_xtrl = h_xtrl_tmp.Clone("h_xtrl")
            self.e_discrimination_last = e_discrimination_last_tmp.Clone(
                "e_discrimination_last")

            '''
            self.e_discrimination_last_20_100 = e_discrimination_last_20_100_tmp.Clone("e_discrimination_last_20_100")
            self.e_discrimination_last_100_250 = e_discrimination_last_100_250_tmp.Clone("e_discrimination_last_100_250")
            self.e_discrimination_last_250_500 = e_discrimination_last_250_500_tmp.Clone("e_discrimination_last_250_500")
            self.e_discrimination_last_500_1000 = e_discrimination_last_500_1000_tmp.Clone("e_discrimination_last_500_1000")
            self.e_discrimination_last_1000_3000 = e_discrimination_last_1000_3000_tmp.Clone("e_discrimination_last_1000_3000")
            self.e_discrimination_last_3000_10000 = e_discrimination_last_3000_10000_tmp.Clone("e_discrimination_last_3000_10000")
            '''

            self.e_discrimination = e_discrimination_tmp.Clone(
                "e_discrimination")
            '''
            self.e_discrimination_20_100 = e_discrimination_20_100_tmp.Clone("e_discrimination_20_100")
            self.e_discrimination_100_250 = e_discrimination_100_250_tmp.Clone("e_discrimination_100_250")
            self.e_discrimination_250_500 = e_discrimination_250_500_tmp.Clone("e_discrimination_250_500")
            self.e_discrimination_500_1000 = e_discrimination_500_1000_tmp.Clone("e_discrimination_500_1000")
            self.e_discrimination_1000_3000 = e_discrimination_1000_3000_tmp.Clone("e_discrimination_1000_3000")
            self.e_discrimination_3000_10000 = e_discrimination_3000_10000_tmp.Clone("e_discrimination_3000_10000")
            '''

            for idx in range(self.n_energy_bins):
                histoName = "h_xtrl_bin_"
                histoName += str(idx)
                self.h_bin_xtrl[idx] = h_bin_xtrl_tmp[idx].Clone(histoName)

            self.h_psd_chargeX = h_psd_chargeX_tmp.Clone("h_psd_chargeX")
            self.h_psd_chargeY = h_psd_chargeY_tmp.Clone("h_psd_chargeY")
            self.h_psd_charge = h_psd_charge_tmp.Clone("h_psd_charge")
            self.h_psd_charge2D = h_psd_charge2D_tmp.Clone("h_psd_charge2D")

            self.h_psd_selected_chargeX = h_psd_selected_chargeX_tmp.Clone(
                "h_psd_selected_chargeX")
            self.h_psd_selected_chargeY = h_psd_selected_chargeY_tmp.Clone(
                "h_psd_selected_chargeY")
            self.h_psd_selected_charge = h_psd_selected_charge_tmp.Clone(
                "h_psd_selected_charge")
            self.h_psd_selected_charge2D = h_psd_selected_charge2D_tmp.Clone(
                "h_psd_selected_charge2D")

            self.h_stk_chargeX = h_stk_chargeX_tmp.Clone("h_stk_chargeX")
            self.h_stk_chargeY = h_stk_chargeY_tmp.Clone("h_stk_chargeY")
            self.h_stk_charge = h_stk_charge_tmp.Clone("h_stk_charge")
            self.h_stk_charge2D = h_stk_charge2D_tmp.Clone("h_stk_charge2D")

            self.h_stk_selected_chargeX = h_stk_selected_chargeX_tmp.Clone(
                "h_stk_selected_chargeX")
            self.h_stk_selected_chargeY = h_stk_selected_chargeY_tmp.Clone(
                "h_stk_selected_chargeY")
            self.h_stk_selected_charge = h_stk_selected_charge_tmp.Clone(
                "h_stk_selected_charge")
            self.h_stk_selected_charge2D = h_stk_selected_charge2D_tmp.Clone(
                "h_stk_selected_charge2D")

            self.h_background_under_xtrl_cut = h_background_under_xtrl_cut_tmp.Clone(
                "h_background_under_xtrl_cut")
            self.h_background_over_xtrl_cut = h_background_over_xtrl_cut_tmp.Clone(
                "h_background_over_xtrl_cut")

            self.first_file_read = False

            # Unlink histos
            self.h_geo_factor.SetDirectory(0)
            self.h_incoming.SetDirectory(0)
            self.h_trigger.SetDirectory(0)
            self.h_geometric_cut.SetDirectory(0)
            self.h_maxElayer_cut.SetDirectory(0)
            self.h_maxBarLayer_cut.SetDirectory(0)
            self.h_BGOTrackContainment_cut.SetDirectory(0)
            self.h_BGO_fiducial_cut.SetDirectory(0)
            self.h_all_cut.SetDirectory(0)

            self.h_geo_factor_nw.SetDirectory(0)
            self.h_incoming_nw.SetDirectory(0)
            self.h_trigger_nw.SetDirectory(0)
            self.h_geometric_cut_nw.SetDirectory(0)
            self.h_maxElayer_cut_nw.SetDirectory(0)
            self.h_maxBarLayer_cut_nw.SetDirectory(0)
            self.h_BGOTrackContainment_cut_nw.SetDirectory(0)
            self.h_BGO_fiducial_cut_nw.SetDirectory(0)
            self.h_all_cut_nw.SetDirectory(0)
            
            self.h_geometric_maxElayer_cut.SetDirectory(0)
            self.h_geometric_maxBarLayer_cut.SetDirectory(0)
            self.h_geometric_BGOTrackContainment_cut.SetDirectory(0)
            self.h_geometric_BGO_fiducial_cut.SetDirectory(0)
            self.h_geometric_all_cut.SetDirectory(0)

            self.h_geometric_maxElayer_cut_nw.SetDirectory(0)
            self.h_geometric_maxBarLayer_cut_nw.SetDirectory(0)
            self.h_geometric_BGOTrackContainment_cut_nw.SetDirectory(0)
            self.h_geometric_BGO_fiducial_cut_nw.SetDirectory(0)
            self.h_geometric_all_cut_nw.SetDirectory(0)

            self.h_BGOfiducial_nBarLayer13_cut.SetDirectory(0)
            self.h_BGOfiducial_maxRms_cut.SetDirectory(0)
            self.h_BGOfiducial_track_selection_cut.SetDirectory(0)
            self.h_BGOfiducial_psd_stk_match_cut.SetDirectory(0)
            self.h_BGOfiducial_psd_charge_cut.SetDirectory(0)
            self.h_BGOfiducial_stk_charge_cut.SetDirectory(0)
            self.h_BGOfiducial_xtrl_cut.SetDirectory(0)
            self.h_BGOfiducial_all_cut.SetDirectory(0)

            self.h_BGOfiducial_nBarLayer13_cut_nw.SetDirectory(0)
            self.h_BGOfiducial_maxRms_cut_nw.SetDirectory(0)
            self.h_BGOfiducial_track_selection_cut_nw.SetDirectory(0)
            self.h_BGOfiducial_psd_stk_match_cut_nw.SetDirectory(0)
            self.h_BGOfiducial_psd_charge_cut_nw.SetDirectory(0)
            self.h_BGOfiducial_stk_charge_cut_nw.SetDirectory(0)
            self.h_BGOfiducial_xtrl_cut_nw.SetDirectory(0)
            self.h_BGOfiducial_all_cut_nw.SetDirectory(0)

            self.h_preGeo_BGOrec_topX_vs_realX.SetDirectory(0)
            self.h_preGeo_BGOrec_topY_vs_realY.SetDirectory(0)
            self.h_preGeo_real_slopeX.SetDirectory(0)
            self.h_preGeo_real_slopeY.SetDirectory(0)
            self.h_preGeo_BGOrec_slopeX.SetDirectory(0)
            self.h_preGeo_BGOrec_slopeY.SetDirectory(0)
            self.h_preGeo_real_interceptX.SetDirectory(0)
            self.h_preGeo_real_interceptY.SetDirectory(0)
            self.h_preGeo_BGOrec_interceptX.SetDirectory(0)
            self.h_preGeo_BGOrec_interceptY.SetDirectory(0)
            self.h_preGeo_real_topMap.SetDirectory(0)
            self.h_preGeo_BGOreco_topMap.SetDirectory(0)
            self.h_preGeo_real_bottomMap.SetDirectory(0)
            self.h_preGeo_BGOreco_bottomMap.SetDirectory(0)
            self.h_noBGOenergy_real_topMap.SetDirectory(0)

            self.h_geo_BGOrec_topX_vs_realX.SetDirectory(0)
            self.h_geo_BGOrec_topY_vs_realY.SetDirectory(0)
            self.h_geo_real_slopeX.SetDirectory(0)
            self.h_geo_real_slopeY.SetDirectory(0)
            self.h_geo_BGOrec_slopeX.SetDirectory(0)
            self.h_geo_BGOrec_slopeY.SetDirectory(0)
            self.h_geo_real_interceptX.SetDirectory(0)
            self.h_geo_real_interceptY.SetDirectory(0)
            self.h_geo_BGOrec_interceptX.SetDirectory(0)
            self.h_geo_BGOrec_interceptY.SetDirectory(0)
            self.h_geo_real_topMap.SetDirectory(0)
            self.h_geo_BGOreco_topMap.SetDirectory(0)
            self.h_geo_real_bottomMap.SetDirectory(0)
            self.h_geo_BGOreco_bottomMap.SetDirectory(0)

            self.h_BGOrec_energy.SetDirectory(0)
            self.h_simu_energy.SetDirectory(0)
            self.h_energy_diff.SetDirectory(0)
            self.h_energy_diff2D.SetDirectory(0)
            self.h_energy_unfold.SetDirectory(0)
            self.h_layer_max_energy_ratio.SetDirectory(0)
            for idx in range(0,14):
                self.h_layer_energy_ratio[idx].SetDirectory(0)
            
            for idx in range(self.n_energy_bins):
                self.h_sumRms_cosine[idx].SetDirectory(0)
            self.sumRms_cosine_20_100.SetDirectory(0)
            self.sumRms_cosine_100_250.SetDirectory(0)
            self.sumRms_cosine_250_500.SetDirectory(0)
            self.sumRms_cosine_500_1000.SetDirectory(0)
            self.sumRms_cosine_1000_3000.SetDirectory(0)
            self.sumRms_cosine_3000_10000.SetDirectory(0)

            self.h_xtrl_energy_int.SetDirectory(0)
            self.h_xtrl.SetDirectory(0)
            self.e_discrimination_last.SetDirectory(0)
            self.e_discrimination_last_20_100.SetDirectory(0)
            self.e_discrimination_last_100_250.SetDirectory(0)
            self.e_discrimination_last_250_500.SetDirectory(0)
            self.e_discrimination_last_500_1000.SetDirectory(0)
            self.e_discrimination_last_1000_3000.SetDirectory(0)
            self.e_discrimination_last_3000_10000.SetDirectory(0)

            self.e_discrimination.SetDirectory(0)
            self.e_discrimination_20_100.SetDirectory(0)
            self.e_discrimination_100_250.SetDirectory(0)
            self.e_discrimination_250_500.SetDirectory(0)
            self.e_discrimination_500_1000.SetDirectory(0)
            self.e_discrimination_1000_3000.SetDirectory(0)
            self.e_discrimination_3000_10000.SetDirectory(0)

            for idx in range(self.n_energy_bins):
                self.h_bin_xtrl[idx].SetDirectory(0)

            self.h_psd_chargeX.SetDirectory(0)
            self.h_psd_chargeY.SetDirectory(0)
            self.h_psd_charge.SetDirectory(0)
            self.h_psd_charge2D.SetDirectory(0)

            self.h_psd_selected_chargeX.SetDirectory(0)
            self.h_psd_selected_chargeY.SetDirectory(0)
            self.h_psd_selected_charge.SetDirectory(0)
            self.h_psd_selected_charge2D.SetDirectory(0)

            self.h_stk_chargeX.SetDirectory(0)
            self.h_stk_chargeY.SetDirectory(0)
            self.h_stk_charge.SetDirectory(0)
            self.h_stk_charge2D.SetDirectory(0)

            self.h_stk_selected_chargeX.SetDirectory(0)
            self.h_stk_selected_chargeY.SetDirectory(0)
            self.h_stk_selected_charge.SetDirectory(0)
            self.h_stk_selected_charge2D.SetDirectory(0)

            self.h_background_under_xtrl_cut.SetDirectory(0)
            self.h_background_over_xtrl_cut.SetDirectory(0)

        else:

            self.h_geo_factor.Add(h_geo_factor_tmp)
            self.h_incoming.Add(h_incoming_tmp)
            self.h_trigger.Add(h_trigger_tmp)
            self.h_geometric_cut.Add(h_geometric_cut_tmp)
            self.h_maxElayer_cut.Add(h_maxElayer_cut_tmp)
            self.h_maxBarLayer_cut.Add(h_maxBarLayer_cut_tmp)
            self.h_BGOTrackContainment_cut.Add(h_BGOTrackContainment_cut_tmp)
            self.h_BGO_fiducial_cut.Add(h_BGO_fiducial_cut_tmp)
            self.h_all_cut.Add(h_all_cut_tmp)

            self.h_geo_factor_nw.Add(h_geo_factor_tmp_nw)
            self.h_incoming_nw.Add(h_incoming_tmp_nw)
            self.h_trigger_nw.Add(h_trigger_tmp_nw)
            self.h_geometric_cut_nw.Add(h_geometric_cut_tmp_nw)
            self.h_maxElayer_cut_nw.Add(h_maxElayer_cut_tmp_nw)
            self.h_maxBarLayer_cut_nw.Add(h_maxBarLayer_cut_tmp_nw)
            self.h_BGOTrackContainment_cut_nw.Add(
                h_BGOTrackContainment_cut_tmp_nw)
            self.h_BGO_fiducial_cut_nw.Add(h_BGO_fiducial_cut_tmp_nw)
            self.h_all_cut_nw.Add(h_all_cut_tmp_nw)

            self.h_geometric_maxElayer_cut.Add(h_geometric_maxElayer_cut_tmp)
            self.h_geometric_maxBarLayer_cut.Add(
                h_geometric_maxBarLayer_cut_tmp)
            self.h_geometric_BGOTrackContainment_cut.Add(
                h_geometric_BGOTrackContainment_cut_tmp)
            self.h_geometric_BGO_fiducial_cut.Add(
                h_geometric_BGO_fiducial_cut_tmp)
            self.h_geometric_all_cut.Add(h_geometric_all_cut_tmp)

            self.h_geometric_maxElayer_cut_nw.Add(
                h_geometric_maxElayer_cut_tmp_nw)
            self.h_geometric_maxBarLayer_cut_nw.Add(
                h_geometric_maxBarLayer_cut_tmp_nw)
            self.h_geometric_BGOTrackContainment_cut_nw.Add(
                h_geometric_BGOTrackContainment_cut_tmp_nw)
            self.h_geometric_BGO_fiducial_cut_nw.Add(
                h_geometric_BGO_fiducial_cut_tmp_nw)
            self.h_geometric_all_cut_nw.Add(h_geometric_all_cut_tmp_nw)

            self.h_BGOfiducial_nBarLayer13_cut.Add(
                h_BGOfiducial_nBarLayer13_cut_tmp)
            self.h_BGOfiducial_maxRms_cut.Add(h_BGOfiducial_maxRms_cut_tmp)
            self.h_BGOfiducial_track_selection_cut.Add(
                h_BGOfiducial_track_selection_cut_tmp)
            self.h_BGOfiducial_psd_stk_match_cut.Add(
                h_BGOfiducial_psd_stk_match_cut_tmp)
            self.h_BGOfiducial_psd_charge_cut.Add(
                h_BGOfiducial_psd_charge_cut_tmp)
            self.h_BGOfiducial_stk_charge_cut.Add(
                h_BGOfiducial_stk_charge_cut_tmp)
            self.h_BGOfiducial_xtrl_cut.Add(h_BGOfiducial_xtrl_cut_tmp)
            self.h_BGOfiducial_all_cut.Add(h_BGOfiducial_all_cut_tmp)

            self.h_BGOfiducial_nBarLayer13_cut_nw.Add(
                h_BGOfiducial_nBarLayer13_cut_tmp_nw)
            self.h_BGOfiducial_maxRms_cut_nw.Add(
                h_BGOfiducial_maxRms_cut_tmp_nw)
            self.h_BGOfiducial_track_selection_cut_nw.Add(
                h_BGOfiducial_track_selection_cut_tmp_nw)
            self.h_BGOfiducial_psd_stk_match_cut_nw.Add(
                h_BGOfiducial_psd_stk_match_cut_tmp_nw)
            self.h_BGOfiducial_psd_charge_cut_nw.Add(
                h_BGOfiducial_psd_charge_cut_tmp_nw)
            self.h_BGOfiducial_stk_charge_cut_nw.Add(
                h_BGOfiducial_stk_charge_cut_tmp_nw)
            self.h_BGOfiducial_xtrl_cut_nw.Add(h_BGOfiducial_xtrl_cut_tmp_nw)
            self.h_BGOfiducial_all_cut_nw.Add(h_BGOfiducial_all_cut_tmp_nw)

            self.h_preGeo_BGOrec_topX_vs_realX.Add(
                h_preGeo_BGOrec_topX_vs_realX_tmp)
            self.h_preGeo_BGOrec_topY_vs_realY.Add(
                h_preGeo_BGOrec_topY_vs_realY_tmp)
            self.h_preGeo_real_slopeX.Add(h_preGeo_real_slopeX_tmp)
            self.h_preGeo_real_slopeY.Add(h_preGeo_real_slopeY_tmp)
            self.h_preGeo_BGOrec_slopeX.Add(h_preGeo_BGOrec_slopeX_tmp)
            self.h_preGeo_BGOrec_slopeY.Add(h_preGeo_BGOrec_slopeY_tmp)
            self.h_preGeo_real_interceptX.Add(h_preGeo_real_interceptX_tmp)
            self.h_preGeo_real_interceptY.Add(h_preGeo_real_interceptY_tmp)
            self.h_preGeo_BGOrec_interceptX.Add(h_preGeo_BGOrec_interceptX_tmp)
            self.h_preGeo_BGOrec_interceptY.Add(h_preGeo_BGOrec_interceptY_tmp)
            self.h_preGeo_real_topMap.Add(h_preGeo_real_topMap_tmp)
            self.h_preGeo_BGOreco_topMap.Add(h_preGeo_BGOreco_topMap_tmp)
            self.h_preGeo_real_bottomMap.Add(h_preGeo_real_bottomMap_tmp)
            self.h_preGeo_BGOreco_bottomMap.Add(h_preGeo_BGOreco_bottomMap_tmp)
            self.h_noBGOenergy_real_topMap.Add(h_noBGOenergy_real_topMap_tmp)

            self.h_geo_BGOrec_topX_vs_realX.Add(h_geo_BGOrec_topX_vs_realX_tmp)
            self.h_geo_BGOrec_topY_vs_realY.Add(h_geo_BGOrec_topY_vs_realY_tmp)
            self.h_geo_real_slopeX.Add(h_geo_real_slopeX_tmp)
            self.h_geo_real_slopeY.Add(h_geo_real_slopeY_tmp)
            self.h_geo_BGOrec_slopeX.Add(h_geo_BGOrec_slopeX_tmp)
            self.h_geo_BGOrec_slopeY.Add(h_geo_BGOrec_slopeY_tmp)
            self.h_geo_real_interceptX.Add(h_geo_real_interceptX_tmp)
            self.h_geo_real_interceptY.Add(h_geo_real_interceptY_tmp)
            self.h_geo_BGOrec_interceptX.Add(h_geo_BGOrec_interceptX_tmp)
            self.h_geo_BGOrec_interceptY.Add(h_geo_BGOrec_interceptY_tmp)
            self.h_geo_real_topMap.Add(h_geo_real_topMap_tmp)
            self.h_geo_BGOreco_topMap.Add(h_geo_BGOreco_topMap_tmp)
            self.h_geo_real_bottomMap.Add(h_geo_real_bottomMap_tmp)
            self.h_geo_BGOreco_bottomMap.Add(h_geo_BGOreco_bottomMap_tmp)

            self.h_BGOrec_energy.Add(h_BGOrec_energy_tmp)
            self.h_simu_energy.Add(h_simu_energy_tmp)
            self.h_energy_diff.Add(h_energy_diff_tmp)
            self.h_energy_diff2D.Add(h_energy_diff2D_tmp)
            self.h_energy_unfold.Add(h_energy_unfold_tmp)
            self.h_layer_max_energy_ratio.Add(h_layer_max_energy_ratio_tmp)
            for idx in range(0, 14):
                self.h_layer_energy_ratio[idx].Add(
                    h_layer_energy_ratio_tmp[idx])

            for idx in range(self.n_energy_bins):
                self.h_sumRms_cosine[idx].Add(h_sumRms_cosine_tmp[idx])
            self.sumRms_cosine_20_100.Add(sumRms_cosine_20_100_tmp)
            self.sumRms_cosine_100_250.Add(sumRms_cosine_100_250_tmp)
            self.sumRms_cosine_250_500.Add(sumRms_cosine_250_500_tmp)
            self.sumRms_cosine_500_1000.Add(sumRms_cosine_500_1000_tmp)
            self.sumRms_cosine_1000_3000.Add(sumRms_cosine_1000_3000_tmp)
            self.sumRms_cosine_3000_10000.Add(sumRms_cosine_3000_10000_tmp)

            self.h_xtrl_energy_int.Add(h_xtrl_energy_int_tmp)
            self.h_xtrl.Add(h_xtrl_tmp)

            self.e_discrimination_last.Add(e_discrimination_last_tmp)
            '''
            self.e_discrimination_last_20_100.Add(e_discrimination_last_20_100_tmp)
            self.e_discrimination_last_100_250.Add(e_discrimination_last_100_250_tmp)
            self.e_discrimination_last_250_500.Add(e_discrimination_last_250_500_tmp)
            self.e_discrimination_last_500_1000.Add(e_discrimination_last_500_1000_tmp)
            self.e_discrimination_last_1000_3000.Add(e_discrimination_last_1000_3000_tmp)
            self.e_discrimination_last_3000_10000.Add(e_discrimination_last_3000_10000_tmp)
            '''

            self.e_discrimination.Add(e_discrimination_tmp)
            '''
            self.e_discrimination_20_100.Add(e_discrimination_20_100_tmp)
            self.e_discrimination_100_250.Add(e_discrimination_100_250_tmp)
            self.e_discrimination_250_500.Add(e_discrimination_250_500_tmp)
            self.e_discrimination_500_1000.Add(e_discrimination_500_1000_tmp)
            self.e_discrimination_1000_3000.Add(e_discrimination_1000_3000_tmp)
            self.e_discrimination_3000_10000.Add(e_discrimination_3000_10000_tmp)
            '''

            for idx in range(self.n_energy_bins):
                self.h_bin_xtrl[idx].Add(h_bin_xtrl_tmp[idx])

            self.h_psd_chargeX.Add(h_psd_chargeX_tmp)
            self.h_psd_chargeY.Add(h_psd_chargeY_tmp)
            self.h_psd_charge.Add(h_psd_charge_tmp)
            self.h_psd_charge2D.Add(h_psd_charge2D_tmp)

            self.h_psd_selected_chargeX.Add(h_psd_selected_chargeX_tmp)
            self.h_psd_selected_chargeY.Add(h_psd_selected_chargeY_tmp)
            self.h_psd_selected_charge.Add(h_psd_selected_charge_tmp)
            self.h_psd_selected_charge2D.Add(h_psd_selected_charge2D_tmp)

            self.h_stk_chargeX.Add(h_stk_chargeX_tmp)
            self.h_stk_chargeY.Add(h_stk_chargeY_tmp)
            self.h_stk_charge.Add(h_stk_charge_tmp)
            self.h_stk_charge2D.Add(h_stk_charge2D_tmp)

            self.h_stk_selected_chargeX.Add(h_stk_selected_chargeX_tmp)
            self.h_stk_selected_chargeY.Add(h_stk_selected_chargeY_tmp)
            self.h_stk_selected_charge.Add(h_stk_selected_charge_tmp)
            self.h_stk_selected_charge2D.Add(h_stk_selected_charge2D_tmp)

            self.h_background_under_xtrl_cut.Add(
                h_background_under_xtrl_cut_tmp)
            self.h_background_over_xtrl_cut.Add(h_background_over_xtrl_cut_tmp)
    
    def write_histos(self, out_file_path, verbose):
        
        # Create output file for full histos
        fOut = TFile.Open(out_file_path, "RECREATE")
        if fOut.IsOpen():
            if verbose:
                print('Output TFile has been created: {}'.format(out_file_path))
        else:
            print('Error creating output TFile: {}'.format(out_file_path))
            sys.exit()

        # Writing final histos to file
        self.h_geo_factor.Write()
        self.h_incoming.Write()
        self.h_trigger.Write()
        self.h_geometric_cut.Write()
        self.h_maxElayer_cut.Write()
        self.h_maxBarLayer_cut.Write()
        self.h_BGOTrackContainment_cut.Write()
        self.h_BGO_fiducial_cut.Write()
        self.h_all_cut.Write()
    
        self.h_geometric_maxElayer_cut.Write()
        self.h_geometric_maxBarLayer_cut.Write()
        self.h_geometric_BGOTrackContainment_cut.Write()
        self.h_geometric_BGO_fiducial_cut.Write()
        self.h_geometric_all_cut.Write()

        self.h_BGOfiducial_nBarLayer13_cut = Write()
        self.h_BGOfiducial_maxRms_cut.Write()
        self.h_BGOfiducial_track_selection_cut.Write()
        self.h_BGOfiducial_psd_stk_match_cut.Write()
        self.h_BGOfiducial_psd_charge_cut.Write()
        self.h_BGOfiducial_stk_charge_cut.Write()
        self.h_BGOfiducial_xtrl_cut.Write()
        self.h_BGOfiducial_all_cut.Write()
        
        fOut.mkdir("nw_histos")
        fOut.cd("nw_histos")

        self.h_geo_factor_nw.Write()
        self.h_incoming_nw.Write()
        self.h_trigger_nw.Write()
        self.h_geometric_cut_nw.Write()
        self.h_maxElayer_cut_nw.Write()
        self.h_maxBarLayer_cut_nw.Write()
        self.h_BGOTrackContainment_cut_nw.Write()
        self.h_BGO_fiducial_cut_nw.Write()
        self.h_all_cut_nw.Write()
    
        self.h_geometric_maxElayer_cut_nw.Write()
        self.h_geometric_maxBarLayer_cut_nw.Write()
        self.h_geometric_BGOTrackContainment_cut_nw.Write()
        self.h_geometric_BGO_fiducial_cut_nw.Write()
        self.h_geometric_all_cut_nw.Write()

        self.h_BGOfiducial_nBarLayer13_cut_nw.Write()
        self.h_BGOfiducial_maxRms_cut_nw.Write()
        self.h_BGOfiducial_track_selection_cut_nw.Write()
        self.h_BGOfiducial_psd_stk_match_cut_nw.Write()
        self.h_BGOfiducial_psd_charge_cut_nw.Write()
        self.h_BGOfiducial_stk_charge_cut_nw.Write()
        self.h_BGOfiducial_xtrl_cut_nw.Write()
        self.h_BGOfiducial_all_cut_nw.Write()

        fOut.mkdir("Analysis_preGeoCut")
        fOut.cd("Analysis_preGeoCut")

        self.h_preGeo_BGOrec_topX_vs_realX.Write()
        self.h_preGeo_BGOrec_topY_vs_realY.Write()
        self.h_preGeo_real_slopeX.Write()
        self.h_preGeo_real_slopeY.Write()
        self.h_preGeo_BGOrec_slopeX.Write()
        self.h_preGeo_BGOrec_slopeY.Write()
        self.h_preGeo_real_interceptX.Write()
        self.h_preGeo_real_interceptY.Write()
        self.h_preGeo_BGOrec_interceptX.Write()
        self.h_preGeo_BGOrec_interceptY.Write()
        self.h_preGeo_real_topMap.Write()
        self.h_preGeo_BGOreco_topMap.Write()
        self.h_preGeo_real_bottomMap.Write()
        self.h_preGeo_BGOreco_bottomMap.Write()
        self.h_noBGOenergy_real_topMap.Write()
        
        fOut.mkdir("Analysis_GeoCut")
        fOut.cd("Analysis_GeoCut")

        self.h_geo_BGOrec_topX_vs_realX.Write()
        self.h_geo_BGOrec_topY_vs_realY.Write()
        self.h_geo_real_slopeX.Write()
        self.h_geo_real_slopeY.Write()
        self.h_geo_BGOrec_slopeX.Write()
        self.h_geo_BGOrec_slopeY.Write()
        self.h_geo_real_interceptX.Write()
        self.h_geo_real_interceptY.Write()
        self.h_geo_BGOrec_interceptX.Write()
        self.h_geo_BGOrec_interceptY.Write()
        self.h_geo_real_topMap.Write()
        self.h_geo_BGOreco_topMap.Write()
        self.h_geo_real_bottomMap.Write()
        self.h_geo_BGOreco_bottomMap.Write()
        
        fOut.mkdir("BGO_Energy")
        fOut.cd("BGO_Energy")

        self.h_BGOrec_energy.Write()
        self.h_simu_energy.Write()
        self.h_energy_diff.Write()
        self.h_energy_diff2D.Write()
        self.h_energy_unfold.Write()
        self.h_layer_max_energy_ratio.Write()
        
        for idx in range(0,14):
            self.h_layer_energy_ratio[idx].Write()
        
        for idx in range(self.n_energy_bins):
            self.h_sumRms_cosine[idx].Write()
        
        self.sumRms_cosine_20_100.Write()
        self.sumRms_cosine_100_250.Write()
        self.sumRms_cosine_250_500.Write()
        self.sumRms_cosine_500_1000.Write()
        self.sumRms_cosine_1000_3000.Write()
        self.sumRms_cosine_3000_10000.Write()

        fOut.mkdir("xtrl")
        fOut.cd("xtrl")

        self.h_xtrl_energy_int.Write()
        self.h_xtrl.Write()   
        self.e_discrimination_last.Write()
        self.e_discrimination_last_20_100.Write()
        self.e_discrimination_last_100_250.Write()
        self.e_discrimination_last_250_500.Write()
        self.e_discrimination_last_500_1000.Write()
        self.e_discrimination_last_1000_3000.Write()
        self.e_discrimination_last_3000_10000.Write()
        self.e_discrimination.Write()
        self.e_discrimination_20_100.Write()
        self.e_discrimination_100_250.Write()
        self.e_discrimination_250_500.Write()
        self.e_discrimination_500_1000.Write()
        self.e_discrimination_1000_3000.Write()
        self.e_discrimination_3000_10000.Write()

        for idx in range(self.n_energy_bins):
            self.h_bin_xtrl[idx].Write()

        fOut.mkdir("PSDcharge")
        fOut.cd("PSDcharge")

        self.h_psd_chargeX.Write()
        self.h_psd_chargeY.Write()
        self.h_psd_charge.Write()
        self.h_psd_charge2D.Write()

        self.h_psd_selected_chargeX.Write()
        self.h_psd_selected_chargeY.Write()
        self.h_psd_selected_charge.Write()
        self.h_psd_selected_charge2D.Write()

        fOut.mkdir("STKcharge")
        fOut.cd("STKcharge")

        self.h_stk_chargeX.Write()
        self.h_stk_chargeY.Write()
        self.h_stk_charge.Write()
        self.h_stk_charge2D.Write()

        self.h_stk_selected_chargeX.Write()
        self.h_stk_selected_chargeY.Write()
        self.h_stk_selected_charge.Write()
        self.h_stk_selected_charge2D.Write()

        fOut.mkdir("mc_ancillary")
        fOut.cd("mc_ancillary")

        self.h_background_under_xtrl_cut.Write()
        self.h_background_over_xtrl_cut.Write()

        # Closing output file
        fOut.Close()

def compute_final_histos_mc(condor_dir_list, opts):

    my_mc_histos = mc_histos()

    for dIdx, tmp_dir in enumerate(condor_dir_list):
        tmp_dir += "/outFiles"
        tmp_dir_list = os.listdir(tmp_dir)
        for elm in tmp_dir_list:
            if elm.startswith("analysisOutFile_"):
                rFile_path = tmp_dir + "/" + elm

        # Open ROOT input file
        rFile = TFile.Open(rFile_path, "READ")
        if rFile.IsOpen():
            if opts.verbose:
                print('Reading file {}: {}'.format((dIdx+1), rFile_path))
        else:
            print('Error reading file {}: {}'.format((dIdx+1), rFile_path))
            sys.exit()

        # Check the keys
        if rFile.GetNkeys() == 0:
            continue

        my_mc_histos.add_file(rFile)
        rFile.Close()

    my_mc_histos.write_histos(opts.output, opts.verbose)