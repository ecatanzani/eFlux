#!/usr/bin/env bash

if [ $# -eq 0 ]
  then
    echo "ERROR: No arguments supplied"
    echo >&2 '
***************
*** ABORTED ***
***************
'
    exit 1
fi

source scl_source enable devtoolset-7

LOCAL_ROOT6=/ROOT6/bin/thisroot.sh
CNAF_LOCAL_ROOT6=/storage/gpfs_data/dampe/users/ecatanzani/deps/root-6.22/bin/thisroot.sh

if [ $1 == "cnaf" ]; then
    ROOT=$CNAF_LOCAL_ROOT6
else
    ROOT=$LOCAL_ROOT6
fi

source $ROOT

DYPATH="dylib"
echo "*** Creating dynamic library folder ***"
mkdir $DYPATH
echo "*** Generating dictionary for std::map<double, std::vector<double>> class ***"
rootcling -f $DYPATH/mapdict.cpp  -rmf $DYPATH/libmapdict.rootmap -rml $DYPATH/libmapdict.so  include/LinkDef.h
echo "*** Compile the dictionary as a shared library ***"
g++ -shared -fPIC -o $DYPATH/libmapdict.so $DYPATH/mapdict.cpp `root-config --cflags --libs`
export LD_LIBRARY_PATH=$(pwd)/$DYPATH:$LD_LIBRARY_PATH

echo >&2 '
************
*** DONE *** 
************
'