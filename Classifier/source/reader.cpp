#include "reader.h"

#include <iostream>

std::shared_ptr<TChain> ReadTreeFromFile(const std::string input_list, const char *tree_name, const bool verbose)
{
    std::shared_ptr<TChain> chain = std::make_shared<TChain>(tree_name);
    std::istringstream input_stream(parse_input_file(input_list));
    std::string tmp_str;
    while (input_stream >> tmp_str)
    {
        chain->Add(tmp_str.c_str());
        if (verbose)
			std::cout << "\nAdding " << tmp_str << " to the chain ...";
    }
    return chain;
}

std::string parse_input_file(std::string input_list)
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