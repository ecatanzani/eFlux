#include <iostream>
#include <vector>

#include "TF1.h"
#include "TKey.h"
#include "TFile.h"
#include "TMath.h"

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/RInterface.hxx>

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

void SpectralFit(
    const char *inputfile,
    const char *outputfile="fitresult.root",
    const bool verbose=true,
    const double emin_gev=1000,
    const double emax_gev=2000,
    const std::size_t nbins=100)
{
    TFile* infile = TFile::Open(inputfile, "READ");
    if (infile->IsZombie())
    {
        std::cerr << "\n\nError opening input file [" << inputfile << "]" << std::endl;
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

    auto binning = createLogBinning(emin_gev, emax_gev, nbins);

    // Enable multithreading
    ROOT::EnableImplicitMT();
    // Create RDF
    ROOT::RDataFrame _data_fr(*mytree);
    const double gev = 0.001;

    auto hspectrum = _data_fr.Define("energy_corr_gev", [&gev](double energy) -> double {return energy*gev;}, {"energy_corr"})
                            .Histo1D({"hspectrum", "hspectrum", (int)binning.size()-1, &binning[0]}, "energy_corr_gev", "simu_energy_w");
   
    /*
    auto hspectrum_norm = static_cast<TH1D*>(hspectrum->Clone("hspectrum_norm"));
    hspectrum_norm->Scale(1./hspectrum->GetEntries());
    */

   auto hspectrum_norm = static_cast<TH1D*>(hspectrum->DrawNormalized());
   hspectrum_norm->SetName("hspectrum_norm");

    std::unique_ptr<TF1> fitfunc = std::make_unique<TF1>("fitfunc", "TMath::Power(10, [0] + [1]*log10(x) + [2]*TMath::Power(log10(x),2))", emin_gev, emax_gev);
    
    fitfunc->SetNpx(10000);

    hspectrum_norm->Fit("fitfunc", "IR");

    TFile *outfile = TFile::Open(outputfile, "RECREATE");
    if (outfile->IsZombie())
    {
        std::cerr << "\n\nError writing output file [" << outputfile << "]" << std::endl;
        exit(100);
    }

    hspectrum->Write();
    hspectrum_norm->Write();
    fitfunc->Write();

    outfile->Close();

}
