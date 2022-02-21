#include "utils.h"

std::string expand_output_path(
    AnyOption &opt,
    std::string output)
{
    auto result = output;
    if (opt.getValue("outputdir") || opt.getValue('d'))
        result += "/kompressor_output.root";
    return result;
}

std::string expand_tt_output_path(
    std::string output,
    bool train)
{
    std::string result;
    if (train)
        result = output.substr(0, output.find_last_of(".")) + "_trainset.root";
    else
        result = output.substr(0, output.find_last_of(".")) + "_testset.root";
    return result;
}