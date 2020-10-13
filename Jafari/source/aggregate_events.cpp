#include "aggregate_events.h"

std::shared_ptr<TChain> aggregateDataEventsTChain(
    const std::string listInputPath,
    const bool verbose)
{
    // ****** Access data using ROOT TChain ******
    // Create TChain object
    std::shared_ptr<TChain> dmpch = std::make_shared<TChain>("DmpEvtNtup");

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
    std::string tmp_str;
    while (input_stream >> tmp_str)
    {
        // Add file to the TChain
        dmpch->Add(tmp_str.c_str());
        if (verbose)
            std::cout << "\nAdding element to the chain ... [" << tmp_str << "]";
    }
    
    return dmpch;
}