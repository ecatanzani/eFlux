#ifndef LAMBDA_CONFIG_H
#define LAMBDA_CONFIG_H

#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>

struct lambdas
{
    double start = 0;
    double end = 0;
    int num = 1;
    int step = (end-start)/num;
};

class lambda_config
{
    public:
        lambda_config(const std::string working_dir);
        ~lambda_config(){};
        void PrintLambdaSettings();
        const lambdas GetLambdaStruct();
    private:
        std::string parse_config_file(
            const std::string wd,
            const std::string config_file);
        void get_config_info(const std::string parsed_config);
        void UpdateLambdaSteps();
        lambdas lambda_info;
};

#endif