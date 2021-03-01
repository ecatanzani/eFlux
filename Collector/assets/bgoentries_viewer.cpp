#include <iostream>
#include <vector>

#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include <ROOT/RDataFrame.hxx>

std::vector<double> createLogBinning(
	const double eMin,
	const double eMax,
	const std::size_t n_bins)
{
	std::vector<double> binning(n_bins + 1, 0);
	double log_interval = (log10(eMax) - log10(eMin)) / n_bins;
	for (unsigned int bIdx = 0; bIdx <= n_bins; ++bIdx)
		binning[bIdx] = pow(10, log10(eMin) + bIdx * log_interval);

	return binning;
}

void scanBGOhits(
    const char* input_file,
    const char* output_file,
    const bool verbose = true,
    const double min_bar_energy = 10,   // minimum bar energy (in MeV)
    const double min_energy = 1000,     // minimum event enrgy (in GeV)
    const double max_energy = 2000,     // maximum event energy (in GeV) 
    const std::size_t nbins=100)
{
    TFile* infile = TFile::Open(input_file, "READ");
    if (infile->IsZombie())
    {
        std::cerr << "\n\nError opening input file [" << input_file << "]" << std::endl;
        exit(100);
    }

    TIter nextkey(infile->GetListOfKeys());
    TKey *key = nullptr;
    std::shared_ptr<TTree> mytree;
    while ((key=static_cast<TKey*>(nextkey()))) 
    {
        TObject *obj = key->ReadObj();
        if (obj->IsA()->InheritsFrom(TTree::Class()))
        {
            mytree = std::shared_ptr<TTree>(static_cast<TTree*>(obj));
            break;
        }
    }

    if (verbose)
        std::cout << "\nFound TTree in input file [" << mytree->GetName() << "]\n\n";

    const double gev = 0.001;
    auto binning = createLogBinning(min_energy*gev, max_energy*gev, nbins);
    // Enable multithreading
    ROOT::EnableImplicitMT();
    // Create RDF
    ROOT::RDataFrame _data_fr(*mytree);
    // Create the RDF for BGO analysis
    auto energy_filter = [&min_energy, &max_energy, &gev] (const double energy) -> bool {
        bool status = false;
        if (energy*gev>= min_energy && energy*gev <= max_energy)
            status = true;
        return status;
    };

    auto _fr_bgo_analysis = _data_fr.Filter("evtfilter_good_event==true")
                                    .Filter(energy_filter, {"energy_corr"});

    auto count_bgo_hits = [&min_bar_energy] (const std::vector<std::vector<double>> bar_energy) -> unsigned int {
        unsigned int hits=0;
        for (unsigned int lidx=0; lidx<bar_energy.size(); ++lidx)
            for (unsigned int bidx=0; bidx<bar_energy[lidx].size(); ++bidx)
                if (bar_energy[lidx][bidx]>min_bar_energy)
                    ++hits;
        return hits;
    };

    auto h_hits_filter = _fr_bgo_analysis.Define("hits", count_bgo_hits, {"layerBarEnergy"})
                                        .Histo1D<unsigned int, double>({"h_hits_filter", "h_hits_filter", 100, 0, 400}, "hits", "simu_energy_w");

    auto hspectrum = _data_fr.Define("energy_corr_gev", [&gev](double energy) -> double {return energy*gev;}, {"energy_corr"})
                            .Histo1D({"hspectrum", "hspectrum", (int)binning.size()-1, &binning[0]}, "energy_corr_gev", "simu_energy_w");

    TFile *outfile = TFile::Open(output_file, "RECREATE");
    if (outfile->IsZombie())
    {
        std::cerr << "\n\nError writing output file [" << output_file << "]" << std::endl;
        exit(100);
    }

    h_hits_filter->Write();
    hspectrum->Write();

    outfile->Close();


}