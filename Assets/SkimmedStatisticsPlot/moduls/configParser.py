def parseConfigFile():
	dConfig = {'farmAddress': "", 'data_XRDFS_skimmed_path': "", 'start_year': 0, 'end_year': 0}
	
	config_params = []
	with open("skim_xrootd.conf", "r") as _config:
		for line in _config:
			for word in line.split():
				config_params.append(word)

	for idx, word in enumerate(config_params):
		if word == "farmAddress":
			dConfig['farmAddress'] = config_params[idx+1]
		if word == "data_XRDFS_skimmed_path":
			dConfig['data_XRDFS_skimmed_path'] = config_params[idx+1]
		if word == "start_year":
			dConfig['start_year'] = int(config_params[idx+1])
		if word == "end_year":
			dConfig['end_year'] = int(config_params[idx+1])

	return dConfig