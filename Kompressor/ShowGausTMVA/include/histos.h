#ifndef HISTOS_H
#define HISTOS_H

#include <memory>
#include <vector>

#include "TH1D.h"

#define bgolayers 14
#define nenergybin 50

extern std::vector<std::vector<std::shared_ptr<TH1D>>> getrmslayerhistos();
extern std::vector<std::vector<std::shared_ptr<TH1D>>> getrmslayerhistos_ff(const std::string inputfile);
extern std::vector<std::vector<std::shared_ptr<TH1D>>> getenergyfractionlayerhistos();
extern std::vector<std::vector<std::shared_ptr<TH1D>>> getenergyfractionlayerhistos_ff(const std::string inputfile);

#endif