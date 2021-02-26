#include "TFile.h"
#include "TKey.h"

#include <memory>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/RInterface.hxx>

void ReweightedHistos(
    const char *inputfile,
    const char *outputfile,
    const bool verbose=true)
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

    // Enable multithreading
    ROOT::EnableImplicitMT();
    // Create RDF
    ROOT::RDataFrame _data_fr(*mytree);

    auto h_energy = _data_fr.Histo1D("energy", "simu_energy_w");
    h_energy->SetName("h_energy");
    auto h_energy_corr = _data_fr.Histo1D("energy_corr", "simu_energy_w");
    h_energy_corr->SetName("h_energy_corr");
    auto h_sumRms = _data_fr.Histo1D("sumRms", "simu_energy_w");
    h_sumRms->SetName("h_sumRms");
    auto h_fracLast = _data_fr.Histo1D("fracLast", "simu_energy_w");
    h_fracLast->SetName("h_fracLast");

    auto h_rmsLayer_1 = _data_fr.Histo1D("rmsLayer_1", "simu_energy_w");
    h_rmsLayer_1->SetName("h_rmsLayer_1");
    auto h_rmsLayer_2 = _data_fr.Histo1D("rmsLayer_2", "simu_energy_w");
    h_rmsLayer_2->SetName("h_rmsLayer_2");
    auto h_rmsLayer_3 = _data_fr.Histo1D("rmsLayer_3", "simu_energy_w");
    h_rmsLayer_3->SetName("h_rmsLayer_3");
    auto h_rmsLayer_4 = _data_fr.Histo1D("rmsLayer_4", "simu_energy_w");
    h_rmsLayer_4->SetName("h_rmsLayer_4");
    auto h_rmsLayer_5 = _data_fr.Histo1D("rmsLayer_5", "simu_energy_w");
    h_rmsLayer_5->SetName("h_rmsLayer_5");
    auto h_rmsLayer_6 = _data_fr.Histo1D("rmsLayer_6", "simu_energy_w");
    h_rmsLayer_6->SetName("h_rmsLayer_6");
    auto h_rmsLayer_7 = _data_fr.Histo1D("rmsLayer_7", "simu_energy_w");
    h_rmsLayer_7->SetName("h_rmsLayer_7");
    auto h_rmsLayer_8 = _data_fr.Histo1D("rmsLayer_8", "simu_energy_w");
    h_rmsLayer_8->SetName("h_rmsLayer_8");
    auto h_rmsLayer_9 = _data_fr.Histo1D("rmsLayer_9", "simu_energy_w");
    h_rmsLayer_9->SetName("h_rmsLayer_9");
    auto h_rmsLayer_10 = _data_fr.Histo1D("rmsLayer_10", "simu_energy_w");
    h_rmsLayer_10->SetName("h_rmsLayer_10");
    auto h_rmsLayer_11 = _data_fr.Histo1D("rmsLayer_11", "simu_energy_w");
    h_rmsLayer_11->SetName("h_rmsLayer_11");
    auto h_rmsLayer_12 = _data_fr.Histo1D("rmsLayer_12", "simu_energy_w");
    h_rmsLayer_12->SetName("h_rmsLayer_12");
    auto h_rmsLayer_13 = _data_fr.Histo1D("rmsLayer_13", "simu_energy_w");
    h_rmsLayer_13->SetName("h_rmsLayer_13");
    auto h_rmsLayer_14 = _data_fr.Histo1D("rmsLayer_14", "simu_energy_w");
    h_rmsLayer_14->SetName("h_rmsLayer_14");

    auto h_fracLayer_1 = _data_fr.Histo1D("fracLayer_1", "simu_energy_w");
    h_fracLayer_1->SetName("h_fracLayer_1");
    auto h_fracLayer_2 = _data_fr.Histo1D("fracLayer_2", "simu_energy_w");
    h_fracLayer_2->SetName("h_fracLayer_2");
    auto h_fracLayer_3 = _data_fr.Histo1D("fracLayer_3", "simu_energy_w");
    h_fracLayer_3->SetName("h_fracLayer_3");
    auto h_fracLayer_4 = _data_fr.Histo1D("fracLayer_4", "simu_energy_w");
    h_fracLayer_4->SetName("h_fracLayer_4");
    auto h_fracLayer_5 = _data_fr.Histo1D("fracLayer_5", "simu_energy_w");
    h_fracLayer_5->SetName("h_fracLayer_5");
    auto h_fracLayer_6 = _data_fr.Histo1D("fracLayer_6", "simu_energy_w");
    h_fracLayer_6->SetName("h_fracLayer_6");
    auto h_fracLayer_7 = _data_fr.Histo1D("fracLayer_7", "simu_energy_w");
    h_fracLayer_7->SetName("h_fracLayer_7");
    auto h_fracLayer_8 = _data_fr.Histo1D("fracLayer_8", "simu_energy_w");
    h_fracLayer_8->SetName("h_fracLayer_8");
    auto h_fracLayer_9 = _data_fr.Histo1D("fracLayer_9", "simu_energy_w");
    h_fracLayer_9->SetName("h_fracLayer_9");
    auto h_fracLayer_10 = _data_fr.Histo1D("fracLayer_10", "simu_energy_w");
    h_fracLayer_10->SetName("h_fracLayer_10");
    auto h_fracLayer_11 = _data_fr.Histo1D("fracLayer_11", "simu_energy_w");
    h_fracLayer_11->SetName("h_fracLayer_11");
    auto h_fracLayer_12 = _data_fr.Histo1D("fracLayer_12", "simu_energy_w");
    h_fracLayer_12->SetName("h_fracLayer_12");
    auto h_fracLayer_13 = _data_fr.Histo1D("fracLayer_13", "simu_energy_w");
    h_fracLayer_13->SetName("h_fracLayer_13");
    auto h_fracLayer_14 = _data_fr.Histo1D("fracLayer_14", "simu_energy_w");
    h_fracLayer_14->SetName("h_fracLayer_14");

    auto h_lastbgolayer = _data_fr.Histo1D("lastBGOLayer", "simu_energy_w");
    h_lastbgolayer->SetName("h_lastbgolayer");

    auto h_nbgoentries = _data_fr.Histo1D("nBGOentries", "simu_energy_w");
    h_nbgoentries->SetName("h_nbgoentries");

    auto h_xtrl = _data_fr.Histo1D("xtrl", "simu_energy_w");
    h_xtrl->SetName("h_xtrl");

    auto h_nudadc_1 = _data_fr.Histo1D("NUD_ADC_1", "simu_energy_w");
    h_nudadc_1->SetName("h_nudadc_1");
    auto h_nudadc_2 = _data_fr.Histo1D("NUD_ADC_2", "simu_energy_w");
    h_nudadc_2->SetName("h_nudadc_2");
    auto h_nudadc_3 = _data_fr.Histo1D("NUD_ADC_3", "simu_energy_w");
    h_nudadc_3->SetName("h_nudadc_3");
    auto h_nudadc_4 = _data_fr.Histo1D("NUD_ADC_4", "simu_energy_w");
    h_nudadc_4->SetName("h_nudadc_4");

    auto h_nud_totaladc = _data_fr.Histo1D("NUD_total_ADC_nud_total_adc", "simu_energy_w");
    h_nud_totaladc->SetName("h_nud_totaladc");
    auto h_nud_maxadc = _data_fr.Histo1D("NUD_max_ADC_nud_max_adc", "simu_energy_w");
    h_nud_maxadc->SetName("h_nud_maxadc");

    TFile *outfile = TFile::Open(outputfile, "RECREATE");
    if (outfile->IsZombie())
    {
        std::cerr << "\n\nError writing output file [" << outputfile << "]" << std::endl;
        exit(100);
    }

    h_energy->Write();
    h_energy_corr->Write();
    h_sumRms->Write();
    h_fracLast->Write();

    h_rmsLayer_1->Write();
    h_rmsLayer_2->Write();
    h_rmsLayer_3->Write();
    h_rmsLayer_4->Write();
    h_rmsLayer_5->Write();
    h_rmsLayer_6->Write();
    h_rmsLayer_7->Write();
    h_rmsLayer_8->Write();
    h_rmsLayer_9->Write();
    h_rmsLayer_10->Write();
    h_rmsLayer_11->Write();
    h_rmsLayer_12->Write();
    h_rmsLayer_13->Write();
    h_rmsLayer_14->Write();

    h_fracLayer_1->Write();
    h_fracLayer_2->Write();
    h_fracLayer_3->Write();
    h_fracLayer_4->Write();
    h_fracLayer_5->Write();
    h_fracLayer_6->Write();
    h_fracLayer_7->Write();
    h_fracLayer_8->Write();
    h_fracLayer_9->Write();
    h_fracLayer_10->Write();
    h_fracLayer_11->Write();
    h_fracLayer_12->Write();
    h_fracLayer_13->Write();
    h_fracLayer_14->Write();

    h_lastbgolayer->Write();

    h_nbgoentries->Write();

    h_xtrl->Write();

    h_nudadc_1->Write();
    h_nudadc_2->Write();
    h_nudadc_3->Write();
    h_nudadc_4->Write();

    h_nud_totaladc->Write();
    h_nud_maxadc->Write();

    outfile->Close();
}