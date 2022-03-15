from datetime import datetime, date

def parseConfigFile(opts):
	dConfig = {'farmAddress': "", 'localStorage': "", 'UseLocalStorage': True, 'simu_XRDFS_path': "", 'data_XRDFS_path': "", 'data_XRDFS_skimmed_path': "", 'geometry': "", 'simu_eMin': 0, 'simu_eMax': 0, 'particle': "", 'start_date': datetime.now(), 'end_date': datetime.now(), 'jSet': ""}
	
	config_params = []
	custom_farm = False
	custom_set = False

	if opts.farm:
		farmAddress = opts.farm
		dConfig['farmAddress'] = farmAddress
		custom_farm = True
	if opts.input:
		jSet = opts.input
		dConfig['jSet'] = jSet
		custom_set = True
	
	
	with open("crawlerConfig.conf", "r") as _config:
		for line in _config:
			for word in line.split():
				config_params.append(word)

	for idx, word in enumerate(config_params):
		if word == "farmAddress":
			if not custom_farm:
				dConfig['farmAddress'] = config_params[idx+1]
		if word == "localStorage":
			dConfig['localStorage'] = config_params[idx+1]
		if word == "UseLocalStorage":
			dConfig['UseLocalStorage'] = True if config_params[idx+1] == "True" else False
		if word == "simu_XRDFS_path":
			dConfig['simu_XRDFS_path'] = config_params[idx+1]
		if word == "data_XRDFS_path":
			dConfig['data_XRDFS_path'] = config_params[idx+1]
		if word == "data_XRDFS_skimmed_path":
			dConfig['data_XRDFS_skimmed_path'] = config_params[idx+1]
		if word == "geometry":
			dConfig['geometry'] = config_params[idx+1]
		if word == "simu_eMin":
			dConfig['simu_eMin'] = float(config_params[idx+1])
		if word == "simu_eMax":
			dConfig['simu_eMax'] = float(config_params[idx+1])
		if word == "particle":
			dConfig['particle'] = config_params[idx+1]
		if word == "start_date":
			year = int(config_params[idx+1][:4])
			month = int(config_params[idx+1][5:7])
			day = int(config_params[idx+1][8:10])
			dConfig["start_date"] = date(year, month, day)
		if word == "end_date":
			year = int(config_params[idx+1][:4])
			month = int(config_params[idx+1][5:7])
			day = int(config_params[idx+1][8:10])
			dConfig["end_date"] = date(year, month, day)
		if word == "jSet":
			if not custom_set:
				dConfig['jSet'] = config_params[idx+1]
	
	return dConfig