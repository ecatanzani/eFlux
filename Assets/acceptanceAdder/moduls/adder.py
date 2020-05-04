from ROOT import TFile, TH1D, TGraph


def compute_final_histos(condor_dir_list, opts):
    # Full statistic histos
    h_incoming = TH1D()
    h_gometric_cut = TH1D()
    h_maxElayer_cut = TH1D()
    h_maxBarLayer_cut = TH1D()
    h_BGOTrackContainment_cut = TH1D()
    h_BGO_fiducial_cut = TH1D()
    h_nBarLayer13_cut = TH1D()
    h_maxRms_cut = TH1D()
    h_track_selection_cut = TH1D()
    h_xtrl_cut = TH1D()
    h_psd_charge_cut = TH1D()
    h_all_cut = TH1D()

    #Ratio histos
    h_ratio_gometric_cut = TH1D()
    h_ratio_maxElayer_cut = TH1D()
    h_ratio_maxBarLayer_cut = TH1D()
    h_ratio_BGOTrackContainment_cut = TH1D()
    h_ratio_BGO_fiducial = TH1D()
    h_ratio_nBarLayer13_cut = TH1D()
    h_ratio_maxRms_cut = TH1D()
    h_ratio_track_selection_cut = TH1D()
    h_ratio_xtrl_cut = TH1D()
    h_ratio_psd_charge_cut = TH1D()
    h_ratio_all_cut = TH1D()

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
        
        # Reading histos
        h_incoming_tmp = rFile.Get("h_incoming")
        h_gometric_cut_tmp = rFile.Get("h_gometric_cut")
        h_maxElayer_cut_tmp = rFile.Get("h_maxElayer_cut")
        h_maxBarLayer_cut_tmp = rFile.Get("h_maxBarLayer_cut")
        h_BGOTrackContainment_cut_tmp = rFile.Get("h_BGOTrackContainment_cut")
        h_BGO_fiducial_cut_tmp = rFile.Get("h_BGO_fiducial")
        h_nBarLayer13_cut_tmp = rFile.Get("h_nBarLayer13_cut")
        h_maxRms_cut_tmp = rFile.Get("h_maxRms_cut")
        h_track_selection_cut_tmp = rFile.Get("h_track_selection_cut")
        h_xtrl_cut_tmp = rFile.Get("h_xtrl_cut")
        h_psd_charge_cut_tmp = rFile.Get("h_psd_charge_cut")
        h_all_cut_tmp = rFile.Get("h_all_cut")

        h_incoming_tmp.SetDirectory(0)
        h_gometric_cut_tmp.SetDirectory(0)
        h_maxElayer_cut_tmp.SetDirectory(0)
        h_maxBarLayer_cut_tmp.SetDirectory(0)
        h_BGOTrackContainment_cut_tmp.SetDirectory(0)
        h_BGO_fiducial_cut_tmp.SetDirectory(0)
        h_nBarLayer13_cut_tmp.SetDirectory(0)
        h_maxRms_cut_tmp.SetDirectory(0)
        h_track_selection_cut_tmp.SetDirectory(0)
        h_xtrl_cut_tmp.SetDirectory(0)
        h_psd_charge_cut_tmp.SetDirectory(0)
        h_all_cut_tmp.SetDirectory(0)

        # Clone output file
        rFile.Close()

        # Add histos
        if dIdx == 0:
            h_incoming = h_incoming_tmp.Clone("h_incoming")
            h_gometric_cut = h_gometric_cut_tmp.Clone("h_gometric_cut")
            h_maxElayer_cut = h_maxElayer_cut_tmp.Clone("h_maxElayer_cut")
            h_maxBarLayer_cut = h_maxBarLayer_cut_tmp.Clone("h_maxBarLayer_cut")
            h_BGOTrackContainment_cut = h_BGOTrackContainment_cut_tmp.Clone("h_BGOTrackContainment_cut")
            h_BGO_fiducial_cut = h_BGO_fiducial_cut_tmp.Clone("h_BGO_fiducial_cut")
            h_nBarLayer13_cut = h_nBarLayer13_cut_tmp.Clone("h_nBarLayer13_cut")
            h_maxRms_cut = h_maxRms_cut_tmp.Clone("h_maxRms_cut")
            h_track_selection_cut = h_track_selection_cut_tmp.Clone("h_track_selection_cut")
            h_xtrl_cut = h_xtrl_cut_tmp.Clone("h_xtrl_cut")
            h_psd_charge_cut = h_psd_charge_cut_tmp.Clone("h_psd_charge_cut")
            h_all_cut = h_all_cut_tmp.Clone("h_all_cut")

        else:
            h_incoming.Add(h_incoming_tmp)
            h_gometric_cut.Add(h_gometric_cut_tmp)
            h_maxElayer_cut.Add(h_maxElayer_cut_tmp)
            h_maxBarLayer_cut.Add(h_maxBarLayer_cut_tmp)
            h_BGOTrackContainment_cut.Add(h_BGOTrackContainment_cut_tmp)
            h_BGO_fiducial_cut.Add(h_BGO_fiducial_cut_tmp)
            h_nBarLayer13_cut.Add(h_nBarLayer13_cut_tmp)
            h_maxRms_cut.Add(h_maxRms_cut_tmp)
            h_track_selection_cut.Add(h_track_selection_cut_tmp)
            h_xtrl_cut.Add(h_xtrl_cut_tmp)
            h_psd_charge_cut.Add(h_psd_charge_cut_tmp)
            h_all_cut.Add(h_all_cut_tmp)

    # Create output file for full histos
    fOut = TFile.Open(opts.output, "RECREATE")
    if fOut.IsOpen():
        if opts.verbose:
            print('Output TFile has been created: {}'.format(opts.output))
    else:
        print('Error creating output TFile: {}'.format(opts.output))
        sys.exit()

    # Writing final histos to file
    h_incoming.Write()
    h_gometric_cut.Write()
    h_maxElayer_cut.Write()
    h_maxBarLayer_cut.Write()
    h_BGOTrackContainment_cut.Write()
    h_BGO_fiducial_cut.Write()
    h_nBarLayer13_cut.Write()
    h_maxRms_cut.Write()
    h_track_selection_cut.Write()
    h_xtrl_cut.Write()
    h_psd_charge_cut.Write()
    h_all_cut.Write()
   
    # Closing output file
    fOut.Close()