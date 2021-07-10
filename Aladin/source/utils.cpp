#include "utils.h"

std::string expand_output_path(
    AnyOption &opt,
    std::string output)
{
    auto result = output;
    if (opt.getValue("outputdir") || opt.getValue('d'))
        result += "/aladin_output.root";
    return result;
}

std::string st1_path(std::string output)
{
    return output.substr(0, output.find_last_of('.')-1) + std::string("_st1.root");
}


std::string st2_path(std::string output)
{
    return output.substr(0, output.find_last_of('.')-1) + std::string("_st2.root");
}