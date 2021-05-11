#include "main.h"
#include "showgaus.h"

std::string parse_input_file(const std::string input_list)
{
	std::ifstream input_file(input_list.c_str());
	if (!input_file.is_open())
	{
		std::cerr << "\n\nError (100) reading input file list...[" << input_list << "]" << std::endl;
		exit(100);
	}
	std::string input_string((std::istreambuf_iterator<char>(input_file)), (std::istreambuf_iterator<char>()));
	input_file.close();
	return input_string;
}

std::shared_ptr<TChain> getchain(
    const std::string filelist, 
    const bool mc,
    const bool verbose)
{
    std::string tree_name = !mc ? "DmpEvtNtup_gauss" : "DmpMCEvtNtup_gauss";
    std::shared_ptr<TChain> evtch = std::make_shared<TChain>(tree_name.c_str(), "DAMPE event tree");
    std::istringstream input_stream(parse_input_file(filelist));
    std::string tmp_str;
    while (input_stream >> tmp_str)
    {
        evtch->Add(tmp_str.c_str());
        if (verbose)
            std::cout << "\nAdding " << tmp_str << " to the chain ...";
    }
    return evtch;
}

std::vector<std::vector<std::shared_ptr<TH1D>>> getrmslayerhistos()
{
    std::vector<std::vector<std::shared_ptr<TH1D>>> rmslayer (nenergybin);
    for (unsigned int bidx=0; bidx<nenergybin; ++bidx)
    {
        rmslayer[bidx] = std::vector<std::shared_ptr<TH1D>> (bgolayers);
        for (unsigned int lidx=0; lidx<bgolayers; ++lidx)
        {
            std::string tmphistoname = std::string("h_rms_energybin_") + std::to_string(bidx+1) + std::string("_layer_") + std::to_string(lidx+1);
            std::string tmphistotitle = std::string("RMS - energybin ") + std::to_string(bidx+1) + std::string(" - BGO layer ") + std::to_string(lidx+1);
            rmslayer[bidx][lidx] = std::make_shared<TH1D>(tmphistoname.c_str(), tmphistotitle.c_str(), 100, -5, 5);
        }
    }
    return rmslayer;
}

std::vector<std::vector<std::shared_ptr<TH1D>>> getenergyfractionlayerhistos()
{
    std::vector<std::vector<std::shared_ptr<TH1D>>> energyfractionlayer (nenergybin);
    for (unsigned int bidx=0; bidx<nenergybin; ++bidx)
    {
        energyfractionlayer[bidx] = std::vector<std::shared_ptr<TH1D>> (bgolayers);
        for (unsigned int lidx=0; lidx<bgolayers; ++lidx)
        {
            std::string tmphistoname = std::string("h_energyfraction_energybin_") + std::to_string(bidx+1) + std::string("_layer_") + std::to_string(lidx+1);
            std::string tmphistotitle = std::string("Energy Fraction - energybin ") + std::to_string(bidx+1) + std::string(" - BGO layer ") + std::to_string(lidx+1);
            energyfractionlayer[bidx][lidx] = std::make_shared<TH1D>(tmphistoname.c_str(), tmphistotitle.c_str(), 100, -5, 5);
        }
    }
    return energyfractionlayer;
}

void showgaus(in_args input_args)
{
    auto evtch = getchain(
        input_args.input_list, 
        input_args.mc_flag, 
        input_args.verbose);
    ROOT::EnableImplicitMT(input_args.threads);
    ROOT::RDataFrame _data_fr(*evtch);
    if (input_args.verbose) std::cout << "\n\nTotal number of events: " << *(_data_fr.Count());
    auto hrmslayer = getrmslayerhistos();
    auto hefraclayer = getenergyfractionlayerhistos();

    for (int bidx=0; bidx<nenergybin; ++bidx)
    {
        auto bin_filter = [&bidx](const int energy_bin) -> bool { return energy_bin == bidx; };
        _data_fr.Filter(bin_filter, {"energy_bin"})
                .Foreach([&hrmslayer, &hefraclayer, &bidx](const std::vector<double> rms, const std::vector<double> efrac, const double energyw)
                {
                    for (unsigned int lidx=0; lidx<bgolayers; ++lidx)
                    {
                        hrmslayer[bidx][lidx]->Fill(rms[lidx], energyw);
                        hefraclayer[bidx][lidx]->Fill(efrac[lidx], energyw);
                    }
                }, {"rmsLayer_gauss", "fracLayer_gauss", "simu_energy_w_corr"});
    }

    TFile* out = TFile::Open(input_args.output_path.c_str(), "RECREATE");
    if (out->IsZombie())
    {
        std::cerr << "\n\nError writing output ROOT file: [" << input_args.output_path << "]\n\n";
        exit(100);
    }

    for (unsigned int bidx=0; bidx<nenergybin; ++bidx)
    {
        out->mkdir((std::string("RMS/energybin_") + std::to_string(bidx+1)).c_str());
        out->cd((std::string("RMS/energybin_") + std::to_string(bidx+1)).c_str());
        for (unsigned int lidx=0; lidx<bgolayers; ++lidx)
            hrmslayer[bidx][lidx]->Write();
    
        out->mkdir((std::string("ELF/energybin_") + std::to_string(bidx+1)).c_str());
        out->cd((std::string("ELF/energybin_") + std::to_string(bidx+1)).c_str());
        for (unsigned int lidx=0; lidx<bgolayers; ++lidx)
            hefraclayer[bidx][lidx]->Write();
    }
    
    out->Close();
}
