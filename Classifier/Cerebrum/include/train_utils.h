#ifndef TRAIN_UTILS_H
#define TRAIN_UTILS_H

#include <memory>
#include <string>

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCategory.h"
#include "TMVA/IMethod.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/DataLoader.h"

#include "TCut.h"

#include "config.h"

extern std::map<std::string, int> GetTMVAMethods(const std::vector<std::string> mymethod);
extern void SetTMVAVariables(
    std::shared_ptr<TMVA::DataLoader> dataloader,
    const train_vars vars);
extern void SetAllTMVAVariables(std::shared_ptr<TMVA::DataLoader> dataloader);
extern void SetNoNUDTMVAVariables(std::shared_ptr<TMVA::DataLoader> dataloader);
extern void SetNUDTMVAVariables(std::shared_ptr<TMVA::DataLoader> dataloader);
extern void SetTMVACuts(
    TCut &signal_cuts, 
    TCut &background_cuts, 
    const bool verbose);
extern void BookMethods(
    std::shared_ptr<TMVA::Factory> factory,
    std::shared_ptr<TMVA::DataLoader> dataloader,
    std::map<std::string, int> Use);

#endif