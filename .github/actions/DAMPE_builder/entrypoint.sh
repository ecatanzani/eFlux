#!/usr/bin/env bash

# Set DAMPE env
source /root/setup-externals.sh
unset DMPSWSYS
cd /dampe/releases/DmpSoftware-6-0-10/
source bin/thisdmpsw.sh

# Clone repo
_CLONE="git clone --recurse-submodules -b flux_computation https://github.com/ecatanzani/eFlux.git"
${_CLONE}

# Building
cd eFlux
make