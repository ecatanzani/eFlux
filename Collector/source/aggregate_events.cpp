#include "aggregate_events.h"

std::shared_ptr<DmpChain> aggregateEventsDmpChain(
    const std::string listInputPath,
    const bool verbose)
{
    // ****** Access data using DAMPE Chain ******

    // Create DmpChain object
    std::shared_ptr<DmpChain> dmpch = std::make_shared<DmpChain>("CollectionTree");

    // Add MC file list to DmpChain
    //dmpch->AddFromList(getListPath(accInputPath, true).c_str());
    dmpch->AddFromList(listInputPath.c_str());
    if (verbose)
        dmpch->GetListOfFiles()->Print();

    return dmpch;
}

std::shared_ptr<TChain> aggregateEventsTChain(
    const std::string listInputPath,
    const bool verbose)
{
    // ****** Access data using ROOT TChain ******

    // Create TChain object
    std::shared_ptr<TChain> dmpch = std::make_shared<TChain>("CollectionTree");

    // Reading list of MC files
    std::ifstream input_file(listInputPath.c_str());
    if (!input_file.is_open())
    {
        std::cerr << "\n\nError (100) reading input MC file list...[" << listInputPath << "]" << std::endl;
        exit(100);
    }
    std::string input_string((std::istreambuf_iterator<char>(input_file)), (std::istreambuf_iterator<char>()));
    input_file.close();
    std::string tmp_str;
    std::istringstream input_stream(input_string);
    while (input_stream >> tmp_str)
    {
        dmpch->Add(tmp_str.c_str());
        if (verbose)
            std::cout << "\nAdding " << tmp_str << " to the chain ...";
    }

    return dmpch;
}

std::shared_ptr<TChain> aggregateDataEventsTChain(
    const std::string listInputPath,
    const bool verbose,
    const bool skimmed)
{
    // ****** Access data using ROOT TChain ******
    // Create TChain object
    std::shared_ptr<TChain> dmpch;

    // Read DATA file list
    std::ifstream input_file(listInputPath.c_str());
    if (!input_file.is_open())
    {
        std::cerr << "\n\nError (100) reading input DATA file list...[" << listInputPath << "]" << std::endl;
        exit(100);
    }
    std::string input_string((std::istreambuf_iterator<char>(input_file)), (std::istreambuf_iterator<char>()));
    input_file.close();
    std::istringstream input_stream(input_string);

    if (skimmed)
    {
        dmpch = std::make_shared<TChain>("DmpEvtNtup");
        std::string tmp_str;
        while (input_stream >> tmp_str)
        {
            // Add file to the TChain
            dmpch->Add(tmp_str.c_str());
            if (verbose)
                std::cout << "\nAdding " << tmp_str << " to the chain ...";
        }
    }
    else
    {
        dmpch = std::make_shared<TChain>("CollectionTree");
        std::string tmp_str;
        while (input_stream >> tmp_str)
        {
            // Setting up gIOSvc
            gIOSvc->Set("InData/Read", tmp_str.c_str());

            // Add file to the TChain
            dmpch->Add(tmp_str.c_str());
            if (verbose)
                std::cout << "\nAdding " << tmp_str << " to the chain ...";
        }
    }
    
    return dmpch;
}

std::shared_ptr<TChain> aggregateTupleDataEventsTChain(
    const std::string listInputPath,
    std::string &year,
    std::string &month,
    const bool verbose)
{
    // ****** Access data using ROOT TChain ******

    // Create TChain object
    std::shared_ptr<TChain> dmpch = std::make_shared<TChain>("CollectionTree");
    
    // Reading list of nTuples files
    std::ifstream input_file(listInputPath.c_str());
    if (!input_file.is_open())
    {
        std::cerr << "\n\nError (100) reading input nTuples file list...[" << listInputPath << "]" << std::endl;
        exit(100);
    }
    std::string input_string((std::istreambuf_iterator<char>(input_file)), (std::istreambuf_iterator<char>()));
    input_file.close();
    std::string tmp_str;
    std::istringstream input_stream(input_string);
    std::string::size_type sz;
    bool first_file = true;
    while (input_stream >> tmp_str)
    {
        // Setting up gIOSvc
        gIOSvc->Set("InData/Read", tmp_str.c_str());

        if (first_file)
        {   
            auto yIdx = tmp_str.find("/2A/") + 4;
            auto mIdx = yIdx + 4;
            year = tmp_str.substr(yIdx, 4);
            month = tmp_str.substr(mIdx, 2);
            /*
            auto syear = tmp_str.substr(yIdx, 4);
            auto smonth = tmp_str.substr(mIdx, 2);
            year = stoi(syear, &sz);
            month = stoi(smonth, &sz);
            */
            first_file = false;
        }

        // Add file to the TChain
        dmpch->Add(tmp_str.c_str());
        if (verbose)
            std::cout << "\nAdding " << tmp_str << " to the chain ...";
    }

    return dmpch;
}