#ifndef REGULARIZE_H
#define REGULARIZE_H

#include <string>
#include <ROOT/RDataFrame.hxx>

extern void load_tf1s(
    std::string file_path,
    std::vector<TF1> &sumrms_fitfunc,
    std::vector<TF1> &sumrms_fitfunc_err,
    std::vector<TF1> &flast_fitfunc,
    std::vector<TF1> &flast_fitfunc_err,
    const bool verbose);

#endif