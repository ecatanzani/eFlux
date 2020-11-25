#include "list_parser.h"

parser::parser(
    const std::string input_list,
    const bool mc,
    const bool _VERBOSE)
{
    std::string tree_name;
    mc ? tree_name = mc_tree_name : tree_name = data_tree_name;
    evtch = std::make_shared<TChain> (tree_name.c_str(), "DAMPE event tree");
    std::istringstream input_stream(parse_input_file(input_list));
    std::string tmp_str;
    while (input_stream >> tmp_str)
    {
        evtch->Add(tmp_str.c_str());
        if (_VERBOSE)
            std::cout << "\nAdding " << tmp_str << " to the chain ...";
    }
}

std::string parser::parse_input_file(const std::string input_list)
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

std::shared_ptr<TChain> parser::GetEvtTree()
{
    if (evtch)
        return evtch;
    else
        return std::shared_ptr<TChain> (nullptr);
    
}