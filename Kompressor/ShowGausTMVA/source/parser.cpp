#include "parser.h"

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