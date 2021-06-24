#ifndef LAMBDA_CONFIG_H
#define LAMBDA_CONFIG_H

#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>

struct rms_lambdas
{
    double start = 0;
    double end = 0;
    int num = 1;
    double step = 1;
};

struct energylastfraction_lambdas
{
    double start = 0;
    double end = 0;
    int num = 1;
    double step = 1;
};

class lambda_config
{
    public:
        lambda_config(const std::string working_dir);
        ~lambda_config(){};
        void PrintLambdaSettings();
        const rms_lambdas GetRMSLambdaStruct();
        const energylastfraction_lambdas GetELFLambdaStruct();
    private:
        std::string parse_config_file(
            const std::string wd,
            const std::string config_file);
        void get_config_info(const std::string parsed_config);
        rms_lambdas rms_lambda_info;
        energylastfraction_lambdas elf_lambda_info;
};

#endif