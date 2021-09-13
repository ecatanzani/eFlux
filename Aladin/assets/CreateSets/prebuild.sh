#!/usr/bin/env bash

if [ $# -eq 0 ]; then
    echo "ERROR: No arguments supplied"
    echo >&2 '
    ***************
    *** ABORTED ***
    ***************
    '
else
    LOCAL_ROOT6=/ROOT6/bin/thisroot.sh
    CNAF_LOCAL_ROOT6=/storage/gpfs_data/dampe/users/ecatanzani/deps/root-6.22/bin/thisroot.sh
    LOCAL_DEVTOOLSET7="scl_source enable devtoolset-7"
    CNAF_LOCAL_DEVTOOLSET7=/opt/rh/devtoolset-7/enable
    if [ $1 == "cnaf" ]; then
        ROOT=$CNAF_LOCAL_ROOT6
        DEVTOOLSET=$CNAF_LOCAL_DEVTOOLSET7
    else
        ROOT=$LOCAL_ROOT6
        DEVTOOLSET=$LOCAL_DEVTOOLSET7
    fi
    source $DEVTOOLSET
    source $ROOT
    if [ $? -eq 0 ]; then
        DYPATH="dylib"
        echo "*** Creating dynamic library folder ***"
        mkdir $DYPATH
        echo "*** Generating dictionary for std::map<double, std::vector<double>> class ***"
        rootcling -f $DYPATH/mapdict.cpp  -rmf $DYPATH/libmapdict.rootmap -rml $DYPATH/libmapdict.so  include/LinkDef.h
        if [ $? -eq 0 ]; then
            echo "*** Compile the dictionary as a shared library ***"
            g++ -shared -fPIC -o $DYPATH/libmapdict.so $DYPATH/mapdict.cpp `root-config --cflags --libs`
            if [ $? -eq 0 ]; then
                export LD_LIBRARY_PATH=$(pwd)/$DYPATH:$LD_LIBRARY_PATH
                if [ $? -eq 0 ]; then
                    echo >&2 '
    ************
    *** DONE *** 
    ************
    '
                else
                    echo >&2 '
    ***************
    *** ABORTED ***
    ***************
    '
                fi
            else
                echo >&2 '
    ***************
    *** ABORTED ***
    ***************
    '
            fi  
        else
            echo >&2 '
    ***************
    *** ABORTED ***
    ***************
    '
        fi
    else
        echo >&2 '
    ***************
    *** ABORTED ***
    ***************
    '
    fi
fi