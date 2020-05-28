import os
from ROOT import TFile


def getListOfFiles(condor_wd):
    list_dir = []
    skipped_dirs = []
    for tmp_dir in os.listdir(condor_wd):
        if tmp_dir.startswith('job_'):
            full_dir_path = condor_wd + "/" + tmp_dir
            expected_condor_outDir = full_dir_path + "/outFiles"
            if os.path.isdir(expected_condor_outDir):
                _list_dir = os.listdir(expected_condor_outDir)
                tmp_acc_full_path = ""
                for file in _list_dir:
                    if file.endswith(".root"):
                        tmp_acc_full_path = expected_condor_outDir + "/" + file
                        break
                if os.path.isfile(tmp_acc_full_path):
                    tmp_acc_file = TFile.Open(tmp_acc_full_path, "READ")
                    if tmp_acc_file.IsOpen():
                        list_dir.append(full_dir_path)
                    else:
                        skipped_dirs.append(full_dir_path)
                else:
                    skipped_dirs.append(full_dir_path)
            else:
                skipped_dirs.append(full_dir_path)
    return list_dir, skipped_dirs
