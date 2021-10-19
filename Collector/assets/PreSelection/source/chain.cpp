#include "chain.h"

#include <fstream>
#include <sstream>
#include <iostream>

#include "TKey.h"
#include "TFile.h"

inline const std::string get_tree_name(const std::string stream) {
    const std::string file = stream.substr(0, stream.find('\n'));
    TFile* input_file = TFile::Open(file.c_str(), "READ");
    if (!input_file->IsOpen()) {
        std::cerr << "\n\nError reading input file [" << file << "]\n\n";
        exit(100);
    }
    std::string tree_name;
    for (TObject* keyAsObject : *input_file->GetListOfKeys()) {
        auto key = dynamic_cast<TKey*>(keyAsObject);
        if (!strcmp(key->GetClassName(), "TTree"))
            tree_name = static_cast<std::string>(key->GetName());
    }
    input_file->Close();
    return tree_name;
}

inline std::string parse_input_file(const std::string input_list)
{
	std::ifstream input_file(input_list.c_str());
	if (!input_file.is_open()) {
		std::cerr << "\n\nError (100) reading input file list...[" << input_list << "]" << std::endl;
		exit(100);
	}
	std::string input_string((std::istreambuf_iterator<char>(input_file)), (std::istreambuf_iterator<char>()));
	input_file.close();
	return input_string;
}

std::shared_ptr<TChain> GetFileChain(const std::string input, const bool verbose) {
    std::istringstream input_stream(parse_input_file(input));
    std::shared_ptr<TChain> evtch = std::make_shared<TChain> (get_tree_name(input_stream.str()).c_str(), "DAMPE event tree");
    std::string tmp_str;
    while (input_stream >> tmp_str) {
        evtch->Add(tmp_str.c_str());
        if (verbose) std::cout << "\nAdding " << tmp_str << " to the chain ...";
    }
    return evtch;
}