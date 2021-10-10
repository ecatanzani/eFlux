#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "TKey.h"
#include "TTree.h"
#include "TFile.h"
#include <ROOT/RDataFrame.hxx>

inline const std::string get_tree_name(const std::string file) {
    TFile* input_file = TFile::Open(file.c_str(), "READ");
    if (!input_file->IsOpen()) {
        std::cerr << "\n\nError reading input file [" << file << "]\n\n";
        exit(100);
    }
    std::string tree_name{""};
    for (TObject* keyAsObject : *input_file->GetListOfKeys()) {
        auto key = dynamic_cast<TKey*>(keyAsObject);
        if (!strcmp(key->GetClassName(), "TTree"))
            tree_name = static_cast<std::string>(key->GetName());
    }
    return tree_name;
}

void xtrl_cut(
    const char* input_file, 
    const char* output_file,
    const double xtrl_low_edge,
    const double xtrl_high_edge,
    const int threads = 1,
    const bool verbose = true) {

        if (verbose) std::cout << "\nReading input file [" << input_file << "]";
        TFile* infile = TFile::Open(input_file, "READ");
        if (!infile->IsOpen()) {
            std::cerr << "\n\nError opening input file [" << input_file << "]\n\n";
            exit(100);
        }

        std::unique_ptr<TTree> tree = std::unique_ptr<TTree>(static_cast<TTree*>(infile->Get(get_tree_name(std::string(input_file)).c_str())));
        tree->SetDirectory(0);
        infile->Close();

        ROOT::EnableImplicitMT(threads);
        ROOT::RDataFrame _fr(*tree);

        auto xtrl_cut = [&xtrl_low_edge, &xtrl_high_edge] (double xtrl) -> bool {return xtrl>=xtrl_low_edge && xtrl<=xtrl_high_edge ? true : false;};
        auto fr_xtrl_cut = _fr.Filter(xtrl_cut, {"xtrl"});

        fr_xtrl_cut.Snapshot(tree->GetName(), output_file);
        if (verbose) std::cout << "\nOutput file has been written [" << output_file << "]\n\n";
    }