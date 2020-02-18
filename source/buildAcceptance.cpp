#include "myHeader.h"

#include <fstream>
#include <sstream>

#include "TSystem.h"

/**
 * @brief 
 * 
 * @param accInputPath 
 * @param verbose 
 */

std::string getListPath(const std::string accInputPath,const bool MC)
{   
    std::string MClist = accInputPath;
    std::string relMClist;
    if(MC)
    {
        std::string relMClist = "/MC.txt";
        MClist.append(relMClist);    
    }
    else
    {
        std::string relMClist = "/Data.txt";
        MClist.append(relMClist);
    }

    return relMClist;
}   

DmpChain* aggregateEventsDmpChain(const std::string accInputPath,const bool verbose)
{
    // ****** Access data using DAMPE Chain ******

    // Create DmpChain object
    DmpChain* dmpch = new DmpChain("CollectionTree");

    // Add MC file list to DmpChain
    dmpch->AddFromList(getListPath(accInputPath,true).c_str());
    if(verbose)
        dmpch->GetListOfFiles()->Print();
    
    return dmpch;
}

TChain* aggregateEventsTChain(const std::string accInputPath,const bool verbose)
{   
    // ****** Access data using ROOT TChain ******

    // Create TChain object
    TChain* dmpch = new TChain("CollectionTree");

    // Reading list of MC files
    std::ifstream input_file(getListPath(accInputPath,true).c_str());
    if(!input_file.is_open()) {
        std::cerr << "\nERROR 100! File not open " << getListPath(accInputPath,true) << "\n\n";
        exit(123);
    }
    std::string input_string((std::istreambuf_iterator< char >(input_file)), (std::istreambuf_iterator< char >()));
    input_file.close();
    std::string tmp_str;
    std::istringstream input_stream(input_string);
    while(input_stream>>tmp_str)
        dmpch->Add(tmp_str.c_str());
        if(verbose)
            std::cout << "\nAdding " << tmp_str << " to the chain ...";
    return dmpch;
}

void buildAcceptance(const std::string accInputPath,const bool verbose)
{
    gSystem->Load("libDmpEvent.so");

    //auto dmpch = aggregateEventsDmpChain(accInputPath,verbose);
    auto dmpch = aggregateEventsTChain(accInputPath,verbose);


    delete dmpch;
}