#include "reader.h"

#include <fstream>
#include <sstream>
#include <iostream>

#include "TKey.h"
#include "TFile.h"

inline const std::string get_tree_name(const std::string file) {
    TFile* input_file = TFile::Open(file.c_str(), "READ");
    if (!input_file->IsOpen()) {
        std::cerr << "\n\nError reading input file [" << file << "]\n\n";
        exit(100);
    }
    std::string tree_name{""};
    for (TObject* keyAsObject : *input_file->GetListOfKeys()) {
        auto key = dynamic_cast<TKey*>(keyAsObject);
        if (!strcmp(key->GetClassName(), "TTree"))
            tree_name = static_cast<std::string>(key->GetName());
    }
    return tree_name;
}

std::shared_ptr<TChain> BuildChain(const std::string input_list, const bool verbose) {
    std::istringstream input_stream(parse_input_file(input_list));
    std::string tmp_str;
	std::shared_ptr<TChain> chain;
	unsigned int nfiles {0};
    while (input_stream >> tmp_str) {	
		if (!nfiles) chain = std::make_shared<TChain>(get_tree_name(tmp_str).c_str());
        chain->Add(tmp_str.c_str());
        if (verbose) std::cout << "\nAdding " << tmp_str << " to the chain ...";
		++nfiles;
    }
	if (verbose) std::cout << "\n\n" << nfiles << " files have been read from list\n\n";
    return chain;
}

std::shared_ptr<TChain> ReadTreeFromFile(const std::string input_list, const char *tree_name, const bool verbose) {
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

std::string parse_input_file(std::string input_list) {
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