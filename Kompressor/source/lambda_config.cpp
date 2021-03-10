#include "lambda_config.h"

lambda_config::lambda_config(const std::string working_dir)
{
    const std::string local_path = "/Kompressor/config";
    const std::string config_file_name = "lambda_config.conf";
    const auto tmp_config = working_dir.substr(0, working_dir.find("/Collector/config")) + local_path;
    get_config_info(parse_config_file(tmp_config, config_file_name));
    UpdateLambdaSteps();
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
		if (!strcmp(tmp_str.c_str(), "start_value"))
		{
			input_stream >> tmp_str;
			lambda_info.start = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "end_value"))
		{
			input_stream >> tmp_str;
			lambda_info.end = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "number_of_elements"))
		{
			input_stream >> tmp_str;
			lambda_info.num = stoi(tmp_str, &sz);
		}
	}
}

void lambda_config::UpdateLambdaSteps()
{
    lambda_info.step = (lambda_info.end-lambda_info.start)/lambda_info.num;
}

const lambdas lambda_config::GetLambdaStruct()
{
    return lambda_info;
}

void lambda_config::PrintLambdaSettings()
{
    std::cout << "\n***** Input lambdas settings *****\n";
    std::cout << "\nLambda START value: [" << lambda_info.start << "]";
    std::cout << "\nLambda END value: [" << lambda_info.end << "]";
    std::cout << "\nLambda number of values: [" << lambda_info.num << "]";
    std::cout << "\n\n******************************\n\n";
}