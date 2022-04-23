#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <limits>

#include "TKey.h"
#include "TTree.h"
#include "TMath.h"
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

void clean_dset(
    const char* input_file, 
    const char* output_file,
    const int threads = 1,
    const bool verbose = true) {

    if (verbose) std::cout << "\nReading input file [" << input_file << "]";
    TFile* infile = TFile::Open(input_file, "READ");
    if (!infile->IsOpen()) {
        std::cerr << "\n\nError opening input file [" << input_file << "]\n\n";
        exit(100);
    }

    std::unique_ptr<TTree> tree = std::unique_ptr<TTree>(static_cast<TTree*>(infile->Get(get_tree_name(std::string(input_file)).c_str())));
    
    ROOT::EnableImplicitMT(threads);
    ROOT::RDataFrame _fr(*tree);

    auto fr_tmva = _fr.Define("rmslayer_1", "rmsLayer[0]")
                    .Define("rmslayer_2", "rmsLayer[1]")
                    .Define("rmslayer_3", "rmsLayer[2]")
                    .Define("rmslayer_4", "rmsLayer[3]")
                    .Define("rmslayer_5", "rmsLayer[4]")
                    .Define("rmslayer_6", "rmsLayer[5]")
                    .Define("rmslayer_7", "rmsLayer[6]")
                    .Define("rmslayer_8", "rmsLayer[7]")
                    .Define("rmslayer_9", "rmsLayer[8]")
                    .Define("rmslayer_10", "rmsLayer[9]")
                    .Define("rmslayer_11", "rmsLayer[10]")
                    .Define("rmslayer_12", "rmsLayer[11]")
                    .Define("rmslayer_13", "rmsLayer[12]")
                    .Define("rmslayer_14", "rmsLayer[13]")

                    .Define("fraclayer_1", "fracLayer[0]")
                    .Define("fraclayer_2", "fracLayer[1]")
                    .Define("fraclayer_3", "fracLayer[2]")
                    .Define("fraclayer_4", "fracLayer[3]")
                    .Define("fraclayer_5", "fracLayer[4]")
                    .Define("fraclayer_6", "fracLayer[5]")
                    .Define("fraclayer_7", "fracLayer[6]")
                    .Define("fraclayer_8", "fracLayer[7]")
                    .Define("fraclayer_9", "fracLayer[8]")
                    .Define("fraclayer_10", "fracLayer[9]")
                    .Define("fraclayer_11", "fracLayer[10]")
                    .Define("fraclayer_12", "fracLayer[11]")
                    .Define("fraclayer_13", "fracLayer[12]")
                    .Define("fraclayer_14", "fracLayer[13]");

    /*
    auto clean_var = [](
        const double rmslayer_norm_1,
        const double rmslayer_norm_2,
        const double rmslayer_norm_3,
        const double rmslayer_norm_4,
        const double rmslayer_norm_5,
        const double rmslayer_norm_6,
        const double rmslayer_norm_7,
        const double rmslayer_norm_8,
        const double rmslayer_norm_9,
        const double rmslayer_norm_10,
        const double rmslayer_norm_11,
        const double rmslayer_norm_12,
        const double rmslayer_norm_13,
        const double rmslayer_norm_14,
        const double fraclayer_norm_1,
        const double fraclayer_norm_2,
        const double fraclayer_norm_3,
        const double fraclayer_norm_4,
        const double fraclayer_norm_5,
        const double fraclayer_norm_6,
        const double fraclayer_norm_7,
        const double fraclayer_norm_8,
        const double fraclayer_norm_9,
        const double fraclayer_norm_10,
        const double fraclayer_norm_11,
        const double fraclayer_norm_12,
        const double fraclayer_norm_13,
        const double fraclayer_norm_14,
        const double sumrms_norm,
        const double fraclastlayer_norm,
        const double xtrl_norm,
        const double xtrl) -> bool {
        
        float float_max = std::numeric_limits<float>::max();
        float float_min = std::numeric_limits<float>::min();

        return  !(TMath::IsNaN(rmslayer_norm_1) || !(TMath::Finite(rmslayer_norm_1)) || abs(rmslayer_norm_1)<float_min || abs(rmslayer_norm_1)>float_max) &&
                !(TMath::IsNaN(rmslayer_norm_2) || !(TMath::Finite(rmslayer_norm_2)) || abs(rmslayer_norm_2)<float_min || abs(rmslayer_norm_2)>float_max) &&
                !(TMath::IsNaN(rmslayer_norm_3) || !(TMath::Finite(rmslayer_norm_3)) || abs(rmslayer_norm_3)<float_min || abs(rmslayer_norm_3)>float_max) &&
                !(TMath::IsNaN(rmslayer_norm_4) || !(TMath::Finite(rmslayer_norm_4)) || abs(rmslayer_norm_4)<float_min || abs(rmslayer_norm_4)>float_max) &&
                !(TMath::IsNaN(rmslayer_norm_5) || !(TMath::Finite(rmslayer_norm_5)) || abs(rmslayer_norm_5)<float_min || abs(rmslayer_norm_5)>float_max) &&
                !(TMath::IsNaN(rmslayer_norm_6) || !(TMath::Finite(rmslayer_norm_6)) || abs(rmslayer_norm_6)<float_min || abs(rmslayer_norm_6)>float_max) &&
                !(TMath::IsNaN(rmslayer_norm_7) || !(TMath::Finite(rmslayer_norm_7)) || abs(rmslayer_norm_7)<float_min || abs(rmslayer_norm_7)>float_max) &&
                !(TMath::IsNaN(rmslayer_norm_8) || !(TMath::Finite(rmslayer_norm_8)) || abs(rmslayer_norm_8)<float_min || abs(rmslayer_norm_8)>float_max) &&
                !(TMath::IsNaN(rmslayer_norm_9) || !(TMath::Finite(rmslayer_norm_9)) || abs(rmslayer_norm_9)<float_min || abs(rmslayer_norm_9)>float_max) &&
                !(TMath::IsNaN(rmslayer_norm_10) || !(TMath::Finite(rmslayer_norm_10)) || abs(rmslayer_norm_10)<float_min || abs(rmslayer_norm_10)>float_max) &&
                !(TMath::IsNaN(rmslayer_norm_11) || !(TMath::Finite(rmslayer_norm_11)) || abs(rmslayer_norm_11)<float_min || abs(rmslayer_norm_11)>float_max) &&
                !(TMath::IsNaN(rmslayer_norm_12) || !(TMath::Finite(rmslayer_norm_12)) || abs(rmslayer_norm_12)<float_min || abs(rmslayer_norm_12)>float_max) &&
                !(TMath::IsNaN(rmslayer_norm_13) || !(TMath::Finite(rmslayer_norm_13)) || abs(rmslayer_norm_13)<float_min || abs(rmslayer_norm_13)>float_max) &&
                !(TMath::IsNaN(rmslayer_norm_14) || !(TMath::Finite(rmslayer_norm_14)) || abs(rmslayer_norm_14)<float_min || abs(rmslayer_norm_14)>float_max) &&

                !(TMath::IsNaN(fraclayer_norm_1) || !(TMath::Finite(fraclayer_norm_1)) || abs(fraclayer_norm_1)<float_min || abs(fraclayer_norm_1)>float_max) &&
                !(TMath::IsNaN(fraclayer_norm_2) || !(TMath::Finite(fraclayer_norm_2)) || abs(fraclayer_norm_2)<float_min || abs(fraclayer_norm_2)>float_max) &&
                !(TMath::IsNaN(fraclayer_norm_3) || !(TMath::Finite(fraclayer_norm_3)) || abs(fraclayer_norm_3)<float_min || abs(fraclayer_norm_3)>float_max) &&
                !(TMath::IsNaN(fraclayer_norm_4) || !(TMath::Finite(fraclayer_norm_4)) || abs(fraclayer_norm_4)<float_min || abs(fraclayer_norm_4)>float_max) &&
                !(TMath::IsNaN(fraclayer_norm_5) || !(TMath::Finite(fraclayer_norm_5)) || abs(fraclayer_norm_5)<float_min || abs(fraclayer_norm_5)>float_max) &&
                !(TMath::IsNaN(fraclayer_norm_6) || !(TMath::Finite(fraclayer_norm_6)) || abs(fraclayer_norm_6)<float_min || abs(fraclayer_norm_6)>float_max) &&
                !(TMath::IsNaN(fraclayer_norm_7) || !(TMath::Finite(fraclayer_norm_7)) || abs(fraclayer_norm_7)<float_min || abs(fraclayer_norm_7)>float_max) &&
                !(TMath::IsNaN(fraclayer_norm_8) || !(TMath::Finite(fraclayer_norm_8)) || abs(fraclayer_norm_8)<float_min || abs(fraclayer_norm_8)>float_max) &&
                !(TMath::IsNaN(fraclayer_norm_9) || !(TMath::Finite(fraclayer_norm_9)) || abs(fraclayer_norm_9)<float_min || abs(fraclayer_norm_9)>float_max) &&
                !(TMath::IsNaN(fraclayer_norm_10) || !(TMath::Finite(fraclayer_norm_10)) || abs(fraclayer_norm_10)<float_min || abs(fraclayer_norm_10)>float_max) &&
                !(TMath::IsNaN(fraclayer_norm_11) || !(TMath::Finite(fraclayer_norm_11)) || abs(fraclayer_norm_11)<float_min || abs(fraclayer_norm_11)>float_max) &&
                !(TMath::IsNaN(fraclayer_norm_12) || !(TMath::Finite(fraclayer_norm_12)) || abs(fraclayer_norm_12)<float_min || abs(fraclayer_norm_12)>float_max) &&
                !(TMath::IsNaN(fraclayer_norm_13) || !(TMath::Finite(fraclayer_norm_13)) || abs(fraclayer_norm_13)<float_min || abs(fraclayer_norm_13)>float_max) &&
                !(TMath::IsNaN(fraclayer_norm_14) || !(TMath::Finite(fraclayer_norm_14)) || abs(fraclayer_norm_14)<float_min || abs(fraclayer_norm_14)>float_max) &&

                !(TMath::IsNaN(sumrms_norm) || !(TMath::Finite(sumrms_norm)) || abs(sumrms_norm)<float_min || abs(sumrms_norm)>float_max) &&
                !(TMath::IsNaN(fraclastlayer_norm) || !(TMath::Finite(fraclastlayer_norm)) || abs(fraclastlayer_norm)<float_min || abs(fraclastlayer_norm)>float_max) &&
                !(TMath::IsNaN(xtrl_norm) || !(TMath::Finite(xtrl_norm)) || abs(xtrl_norm)<float_min || abs(xtrl_norm)>float_max) &&
                !(TMath::IsNaN(xtrl) || !(TMath::Finite(xtrl)) || abs(xtrl)<float_min || abs(xtrl)>float_max);
    };
    */

    auto clean_var = [](
            const double rmslayer_norm_1,
            const double rmslayer_norm_2,
            const double rmslayer_norm_3,
            const double rmslayer_norm_4,
            const double rmslayer_norm_5,
            const double rmslayer_norm_6,
            const double rmslayer_norm_7,
            const double rmslayer_norm_8,
            const double rmslayer_norm_9,
            const double rmslayer_norm_10,
            const double rmslayer_norm_11,
            const double rmslayer_norm_12,
            const double rmslayer_norm_13,
            const double rmslayer_norm_14,
            const double fraclayer_norm_1,
            const double fraclayer_norm_2,
            const double fraclayer_norm_3,
            const double fraclayer_norm_4,
            const double fraclayer_norm_5,
            const double fraclayer_norm_6,
            const double fraclayer_norm_7,
            const double fraclayer_norm_8,
            const double fraclayer_norm_9,
            const double fraclayer_norm_10,
            const double fraclayer_norm_11,
            const double fraclayer_norm_12,
            const double fraclayer_norm_13,
            const double fraclayer_norm_14,
            const double sumrms_norm,
            const double fraclastlayer_norm,
            const double xtrl) -> bool {
            
            float float_max = std::numeric_limits<float>::max();
            float float_min = std::numeric_limits<float>::min();

            return  !(TMath::IsNaN(rmslayer_norm_1) || !(TMath::Finite(rmslayer_norm_1)) || abs(rmslayer_norm_1)<float_min || abs(rmslayer_norm_1)>float_max) &&
                    !(TMath::IsNaN(rmslayer_norm_2) || !(TMath::Finite(rmslayer_norm_2)) || abs(rmslayer_norm_2)<float_min || abs(rmslayer_norm_2)>float_max) &&
                    !(TMath::IsNaN(rmslayer_norm_3) || !(TMath::Finite(rmslayer_norm_3)) || abs(rmslayer_norm_3)<float_min || abs(rmslayer_norm_3)>float_max) &&
                    !(TMath::IsNaN(rmslayer_norm_4) || !(TMath::Finite(rmslayer_norm_4)) || abs(rmslayer_norm_4)<float_min || abs(rmslayer_norm_4)>float_max) &&
                    !(TMath::IsNaN(rmslayer_norm_5) || !(TMath::Finite(rmslayer_norm_5)) || abs(rmslayer_norm_5)<float_min || abs(rmslayer_norm_5)>float_max) &&
                    !(TMath::IsNaN(rmslayer_norm_6) || !(TMath::Finite(rmslayer_norm_6)) || abs(rmslayer_norm_6)<float_min || abs(rmslayer_norm_6)>float_max) &&
                    !(TMath::IsNaN(rmslayer_norm_7) || !(TMath::Finite(rmslayer_norm_7)) || abs(rmslayer_norm_7)<float_min || abs(rmslayer_norm_7)>float_max) &&
                    !(TMath::IsNaN(rmslayer_norm_8) || !(TMath::Finite(rmslayer_norm_8)) || abs(rmslayer_norm_8)<float_min || abs(rmslayer_norm_8)>float_max) &&
                    !(TMath::IsNaN(rmslayer_norm_9) || !(TMath::Finite(rmslayer_norm_9)) || abs(rmslayer_norm_9)<float_min || abs(rmslayer_norm_9)>float_max) &&
                    !(TMath::IsNaN(rmslayer_norm_10) || !(TMath::Finite(rmslayer_norm_10)) || abs(rmslayer_norm_10)<float_min || abs(rmslayer_norm_10)>float_max) &&
                    !(TMath::IsNaN(rmslayer_norm_11) || !(TMath::Finite(rmslayer_norm_11)) || abs(rmslayer_norm_11)<float_min || abs(rmslayer_norm_11)>float_max) &&
                    !(TMath::IsNaN(rmslayer_norm_12) || !(TMath::Finite(rmslayer_norm_12)) || abs(rmslayer_norm_12)<float_min || abs(rmslayer_norm_12)>float_max) &&
                    !(TMath::IsNaN(rmslayer_norm_13) || !(TMath::Finite(rmslayer_norm_13)) || abs(rmslayer_norm_13)<float_min || abs(rmslayer_norm_13)>float_max) &&
                    !(TMath::IsNaN(rmslayer_norm_14) || !(TMath::Finite(rmslayer_norm_14)) || abs(rmslayer_norm_14)<float_min || abs(rmslayer_norm_14)>float_max) &&

                    !(TMath::IsNaN(fraclayer_norm_1) || !(TMath::Finite(fraclayer_norm_1)) || abs(fraclayer_norm_1)<float_min || abs(fraclayer_norm_1)>float_max) &&
                    !(TMath::IsNaN(fraclayer_norm_2) || !(TMath::Finite(fraclayer_norm_2)) || abs(fraclayer_norm_2)<float_min || abs(fraclayer_norm_2)>float_max) &&
                    !(TMath::IsNaN(fraclayer_norm_3) || !(TMath::Finite(fraclayer_norm_3)) || abs(fraclayer_norm_3)<float_min || abs(fraclayer_norm_3)>float_max) &&
                    !(TMath::IsNaN(fraclayer_norm_4) || !(TMath::Finite(fraclayer_norm_4)) || abs(fraclayer_norm_4)<float_min || abs(fraclayer_norm_4)>float_max) &&
                    !(TMath::IsNaN(fraclayer_norm_5) || !(TMath::Finite(fraclayer_norm_5)) || abs(fraclayer_norm_5)<float_min || abs(fraclayer_norm_5)>float_max) &&
                    !(TMath::IsNaN(fraclayer_norm_6) || !(TMath::Finite(fraclayer_norm_6)) || abs(fraclayer_norm_6)<float_min || abs(fraclayer_norm_6)>float_max) &&
                    !(TMath::IsNaN(fraclayer_norm_7) || !(TMath::Finite(fraclayer_norm_7)) || abs(fraclayer_norm_7)<float_min || abs(fraclayer_norm_7)>float_max) &&
                    !(TMath::IsNaN(fraclayer_norm_8) || !(TMath::Finite(fraclayer_norm_8)) || abs(fraclayer_norm_8)<float_min || abs(fraclayer_norm_8)>float_max) &&
                    !(TMath::IsNaN(fraclayer_norm_9) || !(TMath::Finite(fraclayer_norm_9)) || abs(fraclayer_norm_9)<float_min || abs(fraclayer_norm_9)>float_max) &&
                    !(TMath::IsNaN(fraclayer_norm_10) || !(TMath::Finite(fraclayer_norm_10)) || abs(fraclayer_norm_10)<float_min || abs(fraclayer_norm_10)>float_max) &&
                    !(TMath::IsNaN(fraclayer_norm_11) || !(TMath::Finite(fraclayer_norm_11)) || abs(fraclayer_norm_11)<float_min || abs(fraclayer_norm_11)>float_max) &&
                    !(TMath::IsNaN(fraclayer_norm_12) || !(TMath::Finite(fraclayer_norm_12)) || abs(fraclayer_norm_12)<float_min || abs(fraclayer_norm_12)>float_max) &&
                    !(TMath::IsNaN(fraclayer_norm_13) || !(TMath::Finite(fraclayer_norm_13)) || abs(fraclayer_norm_13)<float_min || abs(fraclayer_norm_13)>float_max) &&
                    !(TMath::IsNaN(fraclayer_norm_14) || !(TMath::Finite(fraclayer_norm_14)) || abs(fraclayer_norm_14)<float_min || abs(fraclayer_norm_14)>float_max) &&

                    !(TMath::IsNaN(sumrms_norm) || !(TMath::Finite(sumrms_norm)) || abs(sumrms_norm)<float_min || abs(sumrms_norm)>float_max) &&
                    !(TMath::IsNaN(fraclastlayer_norm) || !(TMath::Finite(fraclastlayer_norm)) || abs(fraclastlayer_norm)<float_min || abs(fraclastlayer_norm)>float_max) &&
                    !(TMath::IsNaN(xtrl) || !(TMath::Finite(xtrl)) || abs(xtrl)<float_min || abs(xtrl)>float_max);
        };

    /*
    auto _fr_cleaned = _fr.Filter(clean_var, 
        {
            "rmslayer_norm_1", 
            "rmslayer_norm_2", 
            "rmslayer_norm_3", 
            "rmslayer_norm_4",
            "rmslayer_norm_5",
            "rmslayer_norm_6",
            "rmslayer_norm_7",
            "rmslayer_norm_8",
            "rmslayer_norm_9",
            "rmslayer_norm_10",
            "rmslayer_norm_11",
            "rmslayer_norm_12",
            "rmslayer_norm_13",
            "rmslayer_norm_14",

            "fraclayer_norm_1",
            "fraclayer_norm_2",
            "fraclayer_norm_3",
            "fraclayer_norm_4",
            "fraclayer_norm_5",
            "fraclayer_norm_6",
            "fraclayer_norm_7",
            "fraclayer_norm_8",
            "fraclayer_norm_9",
            "fraclayer_norm_10",
            "fraclayer_norm_11",
            "fraclayer_norm_12",
            "fraclayer_norm_13",
            "fraclayer_norm_14",
            
            "sumrms_norm",
            "fraclastlayer_norm",
            "xtrl_norm",
            "xtrl"});
        */

       auto _fr_cleaned = fr_tmva.Filter(clean_var, 
        {
            "rmslayer_1", 
            "rmslayer_2", 
            "rmslayer_3", 
            "rmslayer_4",
            "rmslayer_5",
            "rmslayer_6",
            "rmslayer_7",
            "rmslayer_8",
            "rmslayer_9",
            "rmslayer_10",
            "rmslayer_11",
            "rmslayer_12",
            "rmslayer_13",
            "rmslayer_14",

            "fraclayer_1",
            "fraclayer_2",
            "fraclayer_3",
            "fraclayer_4",
            "fraclayer_5",
            "fraclayer_6",
            "fraclayer_7",
            "fraclayer_8",
            "fraclayer_9",
            "fraclayer_10",
            "fraclayer_11",
            "fraclayer_12",
            "fraclayer_13",
            "fraclayer_14",
            
            "sumRms",
            "fracLast",
            "xtrl"});

    _fr_cleaned.Snapshot(tree->GetName(), output_file);
    if (verbose) {
        std::cout << "\n\n[" << *(_fr.Count())-*(_fr_cleaned.Count()) << "] events have been removed from dataset";
        std::cout << "\nOutput file has been written [" << output_file << "]\n\n";
    }
}