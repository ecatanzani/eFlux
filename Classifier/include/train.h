#ifndef TRAIN_H
#define TRAIN_H

#include "main.h"

#include <memory>
#include <vector>

#include "TCut.h"
#include "TChain.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCategory.h"
#include "TMVA/IMethod.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/DataLoader.h"

extern int Train(in_args input_args);

extern std::map<std::string, int> GetTMVAMethods(const std::vector<std::string> mymethod);
extern void SetTMVAVariables(std::shared_ptr<TMVA::DataLoader> dataloader);
extern void SetTMVACuts(
    TCut &signal_cuts, 
    TCut &background_cuts, 
    const bool verbose);
    
#endif