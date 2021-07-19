#include "utils.h"

std::string expand_output_path(
    AnyOption &opt,
    std::string output)
{
    auto result = output;
    if (opt.getValue("outputdir") || opt.getValue('d'))
        result += "/efficiency_output.root";
    return result;
}
