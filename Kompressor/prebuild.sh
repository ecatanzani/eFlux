#!/usr/bin/env bash

source scl_source enable devtoolset-7
source /ROOT6/bin/thisroot.sh

DYPATH="dylib"
echo "Creating dynamic library folder"
mkdir $DYPATH
echo "Generating dictionary for std::map<std::vector<double>, double> class"
rootcling -f $DYPATH/mapdict.cpp  -rmf $DYPATH/libmapdict.rootmap -rml $DYPATH/libmapdict.so  include/LinkDef.h
echo "Compile the dictionary as a shared library"
g++ -shared -fPIC -o $DYPATH/libmapdict.so $DYPATH/mapdict.cpp `root-config --cflags --libs`
export LD_LIBRARY_PATH=$(pwd)/$DYPATH:$LD_LIBRARY_PATH