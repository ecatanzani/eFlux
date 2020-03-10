import os
import sys
from datetime import datetime
from argparse import ArgumentParser
from ROOT import TFile, TH1D, TGraph
import math

start = datetime.now()


def getGenRadius(opts):
    config_params = []
    radius = 0

    with open(opts.config, "r") as _config:
        for line in _config:
            for word in line.split():
                config_params.append(word)
    for idx, word in enumerate(config_params):
        if word == "generation_vertex_radius":
            radius = float(config_params[idx+1])

    return radius


def getListOfFiles(condor_wd):
    list_dir = []
    for tmp_dir in os.listdir(condor_wd):
        if tmp_dir.startswith('job_'):
            full_dir_path = condor_wd + "/" + tmp_dir
            list_dir.append(full_dir_path)
    return list_dir


def compute_final_histo(condor_dir_list, opts):
    accepted_histo = TH1D()
    incoming_histo = TH1D()

    for dIdx, tmp_dir in enumerate(condor_dir_list):
        tmp_dir += "/outFiles"
        tmp_dir_list = os.listdir(tmp_dir)
        for elm in tmp_dir_list:
            if elm.startswith("analysisOutFile_"):
                rFile_path = tmp_dir + "/" + elm

        rFile = TFile.Open(rFile_path, "READ")
        if rFile.IsOpen():
            if opts.verbose:
                print('Reading file {}: {}'.format((dIdx+1), rFile_path))
        else:
            print('Error reading file {}: {}'.format((dIdx+1), rFile_path))
            sys.exit()
        h_all_cut_tmp = rFile.Get("h_all_cut")
        h_incoming_tmp = rFile.Get("h_incoming")
        h_all_cut_tmp.SetDirectory(0)
        h_incoming_tmp.SetDirectory(0)
        rFile.Close()
        if dIdx == 0:
            accepted_histo = h_all_cut_tmp.Clone("acceptance")
            incoming_histo = h_incoming_tmp.Clone("incoming")
        else:
            accepted_histo.Add(h_all_cut_tmp)
            incoming_histo.Add(h_incoming_tmp)

    sph_radius = getGenRadius(opts)
    sph_surfice = 4*math.pi*math.pow(sph_radius, 2)
    accepted_histo.Divide(incoming_histo)
    accepted_histo.Scale(sph_surfice)

    return accepted_histo


def main(args=None):
    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="Acceptance Adder")

    parser.add_argument("-i", "--input", type=str,
                        dest='input', help='Input condor jobs WD')
    parser.add_argument("-o", "--output", type=str,
                        dest='output', help='Output ROOT TFile')
    parser.add_argument("-c", "--config", type=str,
                        dest='config', help='Acceptance config file')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False,
                        action='store_true', help='run in high verbosity mode')

    opts = parser.parse_args(args)

    file_list = getListOfFiles(opts.input)
    if opts.verbose:
        print('Found {} condor directories...'.format(len(file_list)))

    full_acc_histo = compute_final_histo(file_list, opts)

    fOut = TFile.Open(opts.output, "RECREATE")
    if fOut.IsOpen():
        print('Output TFile has been created: {}'.format(opts.output))
    else:
        print('Error creating output TFile: {}'.format(opts.output))

    full_acc_histo.Write()

    fOut.Close()


if __name__ == "__main__":
    main()
    print(datetime.now()-start)
