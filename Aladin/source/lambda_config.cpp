#include "lambda_config.h"

lambda_config::lambda_config(const std::string working_dir)
{
    const std::string local_path = "/Aladin/config";
    const std::string config_file_name = "lambda_config.conf";
    const auto tmp_config = working_dir.substr(0, working_dir.find("/Collector/config")) + local_path;
    get_config_info(parse_config_file(tmp_config, config_file_name));
}

std::string lambda_config::parse_config_file(
	const std::string wd,
	const std::string config_file)
{
	std::string configPath = wd + "/" + config_file;
	std::ifstream input_file(configPath.c_str());
	if (!input_file.is_open())
	{
		std::cerr << "\nInput config file not found [" << configPath << "]\n\n";
		exit(100);
	}
	std::string input_string(
		(std::istreambuf_iterator<char>(input_file)),
		(std::istreambuf_iterator<char>()));
	input_file.close();
	return input_string;
}

void lambda_config::get_config_info(const std::string parsed_config)
{
	std::string tmp_str;
	std::istringstream input_stream(parsed_config);
	std::string::size_type sz;

	while (input_stream >> tmp_str)
	{
		// First stage lambda values

		if (!strcmp(tmp_str.c_str(), "RMS_lambdas"))
		{
			while(strcmp(tmp_str.c_str(), "start_value"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			rms_lambda_info.start = stod(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "end_value"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			rms_lambda_info.end = stod(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "number_of_elements"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			rms_lambda_info.num = stoi(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "step"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			rms_lambda_info.step = stod(tmp_str, &sz);
		}

		if (!strcmp(tmp_str.c_str(), "SumRMS_lambdas"))
		{
			while(strcmp(tmp_str.c_str(), "start_value"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			sumrms_lambda_info.start = stod(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "end_value"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			sumrms_lambda_info.end = stod(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "number_of_elements"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			sumrms_lambda_info.num = stoi(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "step"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			sumrms_lambda_info.step = stod(tmp_str, &sz);
		}

		if (!strcmp(tmp_str.c_str(), "ELF_lambdas"))
		{
			while(strcmp(tmp_str.c_str(), "start_value"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			elf_lambda_info.start = stod(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "end_value"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			elf_lambda_info.end = stod(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "number_of_elements"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			elf_lambda_info.num = stoi(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "step"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			elf_lambda_info.step = stod(tmp_str, &sz);
		}

		if (!strcmp(tmp_str.c_str(), "ELL_lambdas"))
		{
			while(strcmp(tmp_str.c_str(), "start_value"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			ell_lambda_info.start = stod(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "end_value"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			ell_lambda_info.end = stod(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "number_of_elements"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			ell_lambda_info.num = stoi(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "step"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			ell_lambda_info.step = stod(tmp_str, &sz);
		}

		if (!strcmp(tmp_str.c_str(), "XTRL_lambdas"))
		{
			while(strcmp(tmp_str.c_str(), "start_value"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			xtrl_lambda_info.start = stod(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "end_value"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			xtrl_lambda_info.end = stod(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "number_of_elements"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			xtrl_lambda_info.num = stoi(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "step"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			xtrl_lambda_info.step = stod(tmp_str, &sz);
		}

		// Second stage lambda values
		if (!strcmp(tmp_str.c_str(), "RMS_lambdas_st2"))
		{
			while(strcmp(tmp_str.c_str(), "start_value"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			rms_lambda_info_st2.start = stod(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "end_value"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			rms_lambda_info_st2.end = stod(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "number_of_elements"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			rms_lambda_info_st2.num = stoi(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "step"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			rms_lambda_info_st2.step = stod(tmp_str, &sz);
		}

		if (!strcmp(tmp_str.c_str(), "SumRMS_lambdas_st2"))
		{
			while(strcmp(tmp_str.c_str(), "start_value"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			sumrms_lambda_info_st2.start = stod(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "end_value"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			sumrms_lambda_info_st2.end = stod(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "number_of_elements"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			sumrms_lambda_info_st2.num = stoi(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "step"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			sumrms_lambda_info_st2.step = stod(tmp_str, &sz);
		}

		if (!strcmp(tmp_str.c_str(), "ELF_lambdas_st2"))
		{
			while(strcmp(tmp_str.c_str(), "start_value"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			elf_lambda_info_st2.start = stod(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "end_value"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			elf_lambda_info_st2.end = stod(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "number_of_elements"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			elf_lambda_info_st2.num = stoi(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "step"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			elf_lambda_info_st2.step = stod(tmp_str, &sz);
		}

		if (!strcmp(tmp_str.c_str(), "ELL_lambdas_st2"))
		{
			while(strcmp(tmp_str.c_str(), "start_value"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			ell_lambda_info_st2.start = stod(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "end_value"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			ell_lambda_info_st2.end = stod(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "number_of_elements"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			ell_lambda_info_st2.num = stoi(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "step"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			ell_lambda_info_st2.step = stod(tmp_str, &sz);
		}

		if (!strcmp(tmp_str.c_str(), "XTRL_lambdas_st2"))
		{
			while(strcmp(tmp_str.c_str(), "start_value"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			xtrl_lambda_info_st2.start = stod(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "end_value"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			xtrl_lambda_info_st2.end = stod(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "number_of_elements"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			xtrl_lambda_info_st2.num = stoi(tmp_str, &sz);
			while(strcmp(tmp_str.c_str(), "step"))
				input_stream >> tmp_str;
			input_stream >> tmp_str;
			xtrl_lambda_info_st2.step = stod(tmp_str, &sz);
		}
	}
}

const rms_lambdas lambda_config::GetRMSLambdaStruct()
{
    return rms_lambda_info;
}

const sumrms_lambdas lambda_config::GetSumRMSLambdaStruct()
{
	return sumrms_lambda_info;
}

const elf_lambdas lambda_config::GetELFLambdaStruct()
{
    return elf_lambda_info;
}

const ell_lambdas lambda_config::GetELLLambdaStruct()
{
	return ell_lambda_info;
}

const xtrl_lambdas lambda_config::GetXTRLLambdaStruct()
{
	return xtrl_lambda_info;
}

const rms_lambdas_st2 lambda_config::GetSt2RMSLambdaStruct()
{
    return rms_lambda_info_st2;
}

const sumrms_lambdas_st2 lambda_config::GetSt2SumRMSLambdaStruct()
{
	return sumrms_lambda_info_st2;
}

const elf_lambdas_st2 lambda_config::GetSt2ELFLambdaStruct()
{
    return elf_lambda_info_st2;
}

const ell_lambdas_st2 lambda_config::GetSt2ELLLambdaStruct()
{
	return ell_lambda_info_st2;
}

const xtrl_lambdas_st2 lambda_config::GetSt2XTRLLambdaStruct()
{
	return xtrl_lambda_info_st2;
}

void lambda_config::PrintLambdaSettings()
{
    std::cout << "\n***** Input lambdas settings *****\n";

	std::cout << "\n*** RMS ***\n";
    std::cout << "\nLambda START value: [" << rms_lambda_info.start << "]";
    std::cout << "\nLambda END value: [" << rms_lambda_info.end << "]";
    std::cout << "\nLambda number of values: [" << rms_lambda_info.num << "]";
	std::cout << "\nLambda step: [" << rms_lambda_info.step << "]\n";

	std::cout << "\n*** SUM-RMS ***\n";
    std::cout << "\nLambda START value: [" << sumrms_lambda_info.start << "]";
    std::cout << "\nLambda END value: [" << sumrms_lambda_info.end << "]";
    std::cout << "\nLambda number of values: [" << sumrms_lambda_info.num << "]";
	std::cout << "\nLambda step: [" << sumrms_lambda_info.step << "]\n";

	std::cout << "\n*** ELF ***\n";
    std::cout << "\nLambda START value: [" << elf_lambda_info.start << "]";
    std::cout << "\nLambda END value: [" << elf_lambda_info.end << "]";
    std::cout << "\nLambda number of values: [" << elf_lambda_info.num << "]";
	std::cout << "\nLambda step: [" << elf_lambda_info.step << "]\n";

	std::cout << "\n*** ELL ***\n";
    std::cout << "\nLambda START value: [" << ell_lambda_info.start << "]";
    std::cout << "\nLambda END value: [" << ell_lambda_info.end << "]";
    std::cout << "\nLambda number of values: [" << ell_lambda_info.num << "]";
	std::cout << "\nLambda step: [" << ell_lambda_info.step << "]\n";

	std::cout << "\n*** XTRL ***\n";
    std::cout << "\nLambda START value: [" << xtrl_lambda_info.start << "]";
    std::cout << "\nLambda END value: [" << xtrl_lambda_info.end << "]";
    std::cout << "\nLambda number of values: [" << xtrl_lambda_info.num << "]";
	std::cout << "\nLambda step: [" << xtrl_lambda_info.step << "]";
	
    std::cout << "\n\n******************************\n\n";
}

void lambda_config::PrintSt2LambdaSettings()
{
    std::cout << "\n***** Input lambdas settings *****\n";

	std::cout << "\n*** RMS ***\n";
    std::cout << "\nLambda START value: [" << rms_lambda_info_st2.start << "]";
    std::cout << "\nLambda END value: [" << rms_lambda_info_st2.end << "]";
    std::cout << "\nLambda number of values: [" << rms_lambda_info_st2.num << "]";
	std::cout << "\nLambda step: [" << rms_lambda_info_st2.step << "]\n";

	std::cout << "\n*** SUM-RMS ***\n";
    std::cout << "\nLambda START value: [" << sumrms_lambda_info_st2.start << "]";
    std::cout << "\nLambda END value: [" << sumrms_lambda_info_st2.end << "]";
    std::cout << "\nLambda number of values: [" << sumrms_lambda_info_st2.num << "]";
	std::cout << "\nLambda step: [" << sumrms_lambda_info_st2.step << "]\n";

	std::cout << "\n*** ELF ***\n";
    std::cout << "\nLambda START value: [" << elf_lambda_info_st2.start << "]";
    std::cout << "\nLambda END value: [" << elf_lambda_info_st2.end << "]";
    std::cout << "\nLambda number of values: [" << elf_lambda_info_st2.num << "]";
	std::cout << "\nLambda step: [" << elf_lambda_info_st2.step << "]\n";

	std::cout << "\n*** ELL ***\n";
    std::cout << "\nLambda START value: [" << ell_lambda_info_st2.start << "]";
    std::cout << "\nLambda END value: [" << ell_lambda_info_st2.end << "]";
    std::cout << "\nLambda number of values: [" << ell_lambda_info_st2.num << "]";
	std::cout << "\nLambda step: [" << ell_lambda_info_st2.step << "]\n";

	std::cout << "\n*** XTRL ***\n";
    std::cout << "\nLambda START value: [" << xtrl_lambda_info_st2.start << "]";
    std::cout << "\nLambda END value: [" << xtrl_lambda_info_st2.end << "]";
    std::cout << "\nLambda number of values: [" << xtrl_lambda_info_st2.num << "]";
	std::cout << "\nLambda step: [" << xtrl_lambda_info_st2.step << "]";
	
    std::cout << "\n\n******************************\n\n";
}