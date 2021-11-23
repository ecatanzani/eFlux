#include <iostream>
#include <memory>

#include "TFile.h"
#include "TH2D.h"

void countMaxRMS(const char* input_path, const char* histo_name) {
    // Read input file
    TFile *input_file = TFile::Open(input_path, "READ");
    if (!input_file->IsOpen()) {
        std::cerr << "\n\nError opening input file [" << input_path << "]\n\n";
        exit(100);
    }

    std::unique_ptr<TH2D> histo = std::unique_ptr<TH2D>(static_cast<TH2D*>(input_file->Get(histo_name)));

    double up_events {0};
    double tot_events {0};
    int bin_th = histo->GetYaxis()->FindBin(100);

    for (int bX=1; bX<=histo->GetNbinsX(); ++bX) {
        for (int bY=1; bY<=histo->GetNbinsY()+1; ++bY) {
            if (bY > bin_th)
                up_events += histo->GetBinContent(bX, bY);
            tot_events += histo->GetBinContent(bX, bY);
        }
    }

    std::cout << "\n\nEvents above threshold: " << up_events;
    std::cout << "\nEvents: " << tot_events;
    std::cout << "\n\nFraction of rejected events: " << up_events/tot_events << std::endl;    

}