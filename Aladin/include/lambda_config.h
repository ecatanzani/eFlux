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

struct elf_lambdas
{
    double start = 0;
    double end = 0;
    int num = 1;
    double step = 1;
};

struct sumrms_lambdas
{
    double start = 0;
    double end = 0;
    int num = 1;
    double step = 1;
};

struct ell_lambdas
{
    double start = 0;
    double end = 0;
    int num = 1;
    double step = 1;
};

struct xtrl_lambdas
{
    double start = 0;
    double end = 0;
    int num = 1;
    double step = 1;
};

struct rms_lambdas_st2
{
    double start = 0;
    double end = 0;
    int num = 1;
    double step = 1;
};

struct elf_lambdas_st2
{
    double start = 0;
    double end = 0;
    int num = 1;
    double step = 1;
};

struct sumrms_lambdas_st2
{
    double start = 0;
    double end = 0;
    int num = 1;
    double step = 1;
};

struct ell_lambdas_st2
{
    double start = 0;
    double end = 0;
    int num = 1;
    double step = 1;
};

struct xtrl_lambdas_st2
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
        const sumrms_lambdas GetSumRMSLambdaStruct();
        const elf_lambdas GetELFLambdaStruct();
        const ell_lambdas GetELLLambdaStruct();
        const xtrl_lambdas GetXTRLLambdaStruct();
        void PrintSt2LambdaSettings();
        const rms_lambdas_st2 GetSt2RMSLambdaStruct();
        const sumrms_lambdas_st2 GetSt2SumRMSLambdaStruct();
        const elf_lambdas_st2 GetSt2ELFLambdaStruct();
        const ell_lambdas_st2 GetSt2ELLLambdaStruct();
        const xtrl_lambdas_st2 GetSt2XTRLLambdaStruct();

    private:
        std::string parse_config_file(
            const std::string wd,
            const std::string config_file);
        void get_config_info(const std::string parsed_config);
        rms_lambdas rms_lambda_info;
        sumrms_lambdas sumrms_lambda_info;
        elf_lambdas elf_lambda_info;
        ell_lambdas ell_lambda_info;
        xtrl_lambdas xtrl_lambda_info;

        rms_lambdas_st2 rms_lambda_info_st2;
        sumrms_lambdas_st2 sumrms_lambda_info_st2;
        elf_lambdas_st2 elf_lambda_info_st2;
        ell_lambdas_st2 ell_lambda_info_st2;
        xtrl_lambdas_st2 xtrl_lambda_info_st2;
};

#endif