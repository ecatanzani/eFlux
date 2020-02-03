#!/usr/bin/env bash

# Load DevToolSet-7 for c++14
source /opt/rh/devtoolset-7/enable

# Create 'outFiles' dir if does not exist
_DIRECTORY="/storage/gpfs_data/dampe/users/ecatanzani/myRepos/DAMPE/eFlux/condorJob/outFiles"
if [[ ! -d "$_DIRECTORY" ]]
then
    mkdir $_DIRECTORY
fi

# Run Script
./storage/gpfs_data/dampe/users/ecatanzani/myRepos/DAMPE/eFlux/Release/eFLux -i /storage/gpfs_data/dampe/users/ecatanzani/Data/DAMPE/nTuples/Trees/myTree.root -d $_DIRECTORY -t 34900000 -v