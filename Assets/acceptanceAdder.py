import os
import sys
from datetime import datetime
from argparse import ArgumentParser
from ROOT import TFile, TH1D, TGraph
import math

start = datetime.now()

def getListOfFiles(condor_wd):
    list_dir = []
    skipped_dirs = []
    for tmp_dir in os.listdir(condor_wd):
        if tmp_dir.startswith('job_'):
            full_dir_path = condor_wd + "/" + tmp_dir
            expected_condor_outDir = full_dir_path + "/outFiles"
            if os.path.isdir(expected_condor_outDir):
                list_dir.append(full_dir_path)
            else:
                skipped_dirs.append(full_dir_path)
    return list_dir, skipped_dirs


def compute_final_histos(condor_dir_list, opts):
    h_incoming = TH1D()
    h_maxElateral_cut = TH1D()
    h_maxBarLayer_cut = TH1D()
    h_BGOTrackContainment_cut = TH1D()
    h_nBarLayer13_cut = TH1D()
    h_maxRms_cut = TH1D()
    h_all_cut = TH1D()

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
        h_maxElateral_cut_tmp = rFile.Get("h_maxElateral_cut")
        h_maxBarLayer_cut_tmp = rFile.Get("h_maxBarLayer_cut")
        h_BGOTrackContainment_cut_tmp = rFile.Get("h_BGOTrackContainment_cut")
        h_nBarLayer13_cut_tmp = rFile.Get("h_nBarLayer13_cut")
        h_maxRms_cut_tmp = rFile.Get("h_maxRms_cut")
        h_all_cut_tmp = rFile.Get("h_all_cut")
        
        h_incoming_tmp.SetDirectory(0)
        h_maxElateral_cut_tmp.SetDirectory(0)
        h_maxBarLayer_cut_tmp.SetDirectory(0)
        h_BGOTrackContainment_cut_tmp.SetDirectory(0)
        h_nBarLayer13_cut_tmp.SetDirectory(0)
        h_maxRms_cut_tmp.SetDirectory(0)
        h_all_cut_tmp.SetDirectory(0)

        # Clone output file
        rFile.Close()

        # Add histos
        if dIdx == 0:
            h_incoming = h_incoming_tmp.Clone("h_incoming")
            h_maxElateral_cut = h_maxElateral_cut_tmp.Clone("h_maxElateral_cut")
            h_maxBarLayer_cut = h_maxBarLayer_cut_tmp.Clone("h_maxBarLayer_cut")
            h_BGOTrackContainment_cut = h_BGOTrackContainment_cut_tmp.Clone("h_BGOTrackContainment_cut")
            h_nBarLayer13_cut = h_nBarLayer13_cut_tmp.Clone("h_nBarLayer13_cut")
            h_maxRms_cut = h_maxRms_cut_tmp.Clone("h_maxRms_cut")
            h_all_cut = h_all_cut_tmp.Clone("h_all_cut")
        else:
            h_incoming.Add(h_incoming_tmp)
            h_maxElateral_cut.Add(h_maxElateral_cut_tmp)
            h_maxBarLayer_cut.Add(h_maxBarLayer_cut_tmp)
            h_BGOTrackContainment_cut.Add(h_BGOTrackContainment_cut_tmp)
            h_nBarLayer13_cut.Add(h_nBarLayer13_cut_tmp)
            h_maxRms_cut.Add(h_maxRms_cut_tmp)
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
    h_maxElateral_cut.Write()
    h_maxBarLayer_cut.Write()
    h_BGOTrackContainment_cut.Write()
    h_nBarLayer13_cut.Write()
    h_maxRms_cut.Write()
    h_all_cut.Write()

    # Closing output file
    fOut.Close()


def main(args=None):
    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="Acceptance Adder")

    parser.add_argument("-i", "--input", type=str,
                        dest='input', help='Input condor jobs WD')
    parser.add_argument("-o", "--output", type=str,
                        dest='output', help='Output ROOT TFile')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False,
                        action='store_true', help='run in high verbosity mode')

    opts = parser.parse_args(args)

    good_dirs, skipped_dirs = getListOfFiles(opts.input)
    if opts.verbose:
        print('Found {} GOOD condor directories'.format(len(good_dirs)))
    
    print('Found {} BAD condor directories...'.format(len(skipped_dirs)))
    for idx, elm in enumerate(skipped_dirs):
        print('Skipped {} directory: {}'.format(idx, elm))
    
    compute_final_histos(good_dirs, opts)


if __name__ == "__main__":
    main()
    print(datetime.now()-start)
