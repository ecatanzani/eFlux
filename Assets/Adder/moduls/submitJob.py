import os
import subprocess


def submit_add_job(opts):

    if opts.verbose:
        print("\n **** Submitting adder HTCondor job...\n")

    job_wd = opts.job
    # Find out paths
    outputPath = job_wd + str("/output.log")
    logPath = job_wd + str("/output.clog")
    errPath = job_wd + str("/output.err")
    bashScriptPath = job_wd + str("/script.sh")
    py_exe = opts.location + str("/addTmpFiles.py")

    # Write sub file
    subFilePath = job_wd + str("/crawler.sub")
    try:
        with open(subFilePath, 'w') as outSub:
            outSub.write("universe = vanilla\n")
            outSub.write('executable = {}\n'.format(bashScriptPath))
            outSub.write('output = {}\n'.format(outputPath))
            outSub.write('error = {}\n'.format(errPath))
            outSub.write('log = {}\n'.format(logPath))
            outSub.write("ShouldTransferFiles = YES\n")
            outSub.write("WhenToTransferOutput = ON_EXIT\n")
            outSub.write("queue 1")
    except OSError:
        print('ERROR creating HTCondor sub file in: {}'.format(job_wd))
        raise
    else:
        if opts.verbose:
            print('HTCondor sub file created in: {}'.format(job_wd))

    # Write bash script
    try:
        with open(bashScriptPath, "w") as outScript:
            outScript.write("#!/usr/bin/env bash\n")
            outScript.write(
                "source /cvmfs/dampe.cern.ch/centos7/etc/setup.sh\n")
            outScript.write("dampe_init trunk\n")
            if opts.simulation:
                outScript.write('python {} -l {} -i {} -s -o {} -v'.format(py_exe,
                                                                           opts.location, opts.input, opts.output))
            if opts.data:
                outScript.write('python {} -l {} -i {} -d -o {} -v'.format(py_exe,
                                                                           opts.location, opts.input, opts.output))
    except OSError:
        print('ERROR creating HTCondor bash script file in: {}'.format(job_wd))
        raise
    else:
        if opts.verbose:
            print('HTCondor bash script file created in: {}'.format(job_wd))

    # Make bash script executable
    subprocess.run('chmod +x {}'.format(bashScriptPath),
                   shell=True, check=True)

    # Submit HTCondor job
    subprocess.run(
        ['condor_submit -name sn-01.cr.cnaf.infn.it -spool {}'.format(subFilePath)], shell=True, check=True)
