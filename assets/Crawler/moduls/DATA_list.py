import os
import sys
import subprocess
from datetime import date, timedelta

def extract_date_from_folder(folder):
	year = int(folder[folder.rfind('/')+1:folder.rfind('/')+5])
	month = int(folder[folder.rfind('/')+5:folder.rfind('/')+7])
	day = int(folder[folder.rfind('/')+7:folder.rfind('/')+9])
	return date(year, month, day)

def extract_date_from_file(file):
	year = int(file[file.rfind('_OBS_')+5:file.rfind('_OBS_')+9])
	month = int(file[file.rfind('_OBS_')+9:file.rfind('_OBS_')+11])
	day = int(file[file.rfind('_OBS_')+11:file.rfind('_OBS_')+13])
	return date(year, month, day)

def order_by_date(file_name):
	return file_name[file_name.rfind('_OBS_')+5:file_name.rfind('_OBS_')+13]

def createDATAlist(pars, opts):
	if pars['start_date'] > pars['end_date']:
		print("ERROR: Start date could not be bigger respect to end date")
		sys.exit()

	# Crate output data file list
	dList = []
	years = []
	counters = []

	# Get stage 0 dirs --> /FM/FlightData/2A/
	if pars['UseLocalStorage']:
		dataDirs = [os.path.join(pars['localStorage'], pars['data_XRDFS_path'][1:], _dir) for _dir in os.listdir(os.path.join(pars['localStorage'], pars['data_XRDFS_path'][1:]))]
	else:
		getDataDirsCommand = 'xrdfs {} ls {}'.format(pars['farmAddress'], pars['data_XRDFS_path'])
		if opts.verbose:
			print('Executing XRDFS command: {}'.format(getDataDirsCommand))
		dataDirsOut = subprocess.run(getDataDirsCommand, shell=True, check=True, stdout=subprocess.PIPE)
		dataDirs = str.split(dataDirsOut.stdout.decode('utf-8').rstrip(), '\n')

	if opts.verbose:
		print('Collecting data from {} to {}'.format(pars['start_date'], pars['end_date']))

	# Get stage 1 dirs --> /FM/FlightData/2A/DayOfAcquisition/
	for dir_st1 in dataDirs:

		# Select interesting data dirs, within the selected time window
		if "2A/20" in dir_st1:
			# Extract date from folder
			date_folder = extract_date_from_folder(dir_st1)
			
			if date_folder < pars['start_date'] - timedelta(days=5) or date_folder > pars['end_date'] + timedelta(days=5):
				continue

			if date_folder.year not in years:
				years.append(date_folder.year)
				year_data_idx = len(years)-1
				counters.append(0)

			if pars['UseLocalStorage']:
				dataDirs_st1 = [os.path.join(pars['localStorage'], pars['data_XRDFS_path'][1:], dir_st1, _dir) for _dir in os.listdir(os.path.join(pars['localStorage'], pars['data_XRDFS_path'][1:], dir_st1))]
			else:
				getDataDirsCommand = 'xrdfs {} ls {}'.format(pars['farmAddress'], dir_st1)
				if opts.verbose:
					print('Executing XRDFS command: {}'.format(getDataDirsCommand))
				dataDirsOut = subprocess.run(getDataDirsCommand, shell=True, check=True, stdout=subprocess.PIPE)
				dataDirs_st1 = str.split(dataDirsOut.stdout.decode('utf-8').rstrip(), '\n')

			# Get stage 2 dirs --> /FM/FlightData/2A/DayOfAcquisition/DataDir
			for dir_st2 in dataDirs_st1:
				if pars['UseLocalStorage']:
					dataDirs_st2 = [os.path.join(pars['localStorage'], pars['data_XRDFS_path'][1:], dir_st2, _dir) for _dir in os.listdir(os.path.join(pars['localStorage'], pars['data_XRDFS_path'][1:], dir_st2))]
				else:
					getDataDirsCommand = 'xrdfs {} ls {}'.format(pars['farmAddress'], dir_st2)
					if opts.verbose:
						print('Executing XRDFS command: {}'.format(getDataDirsCommand))
					dataDirsOut = subprocess.run(getDataDirsCommand, shell=True, check=True, stdout=subprocess.PIPE)
					dataDirs_st2 = str.split(dataDirsOut.stdout.decode('utf-8').rstrip(), '\n')

				# Get ROOT data file
				for data_elm in dataDirs_st2:
					if data_elm.endswith('.root'):
						# Check the date of the ROOT file
						date_file = extract_date_from_file(data_elm)
						if date_file >= pars['start_date'] and date_file <= pars['end_date']:
							dList.append(data_elm)
							counters[year_data_idx] += 1
	
	# Sort data file list
	dList.sort(key=order_by_date)

	if opts.verbose:
		print('{} data files have been read...'.format(sum(counters)))
		for year_idx, year in enumerate(years):
			print('{} data files found in {} folder'.format(counters[year_idx], year))

	if opts.output:
		data_list_path = opts.output
	else:
		data_list_path =  "dataFileList.txt"

	# Write output file
	with open(data_list_path, "w") as outList:
		for elm in dList:
			if pars['UseLocalStorage']:
				outList.write(elm + "\n")
			else:
				outList.write(pars['farmAddress'] + elm + "\n")


def createSkimmedDATAlist(pars, opts):
	if pars['data_sYear'] > pars['data_eYear']:
		print("ERROR: Start year could not be bigger respect to end year")
		sys.exit()
	
	# Crate output data file list
	dList = []
	years = []
	counters = []

	# Get stage 0 dirs --> /FM/skim/6.0.0/v2/
	getDataDirsCommand = 'xrdfs {} ls {}'.format(pars['farmAddress'], pars['data_XRDFS_skimmed_path'])
	if opts.verbose:
		print('Executing XRDFS command: {}'.format(getDataDirsCommand))
	dataDirsOut = subprocess.run(getDataDirsCommand, shell=True, check=True, stdout=subprocess.PIPE)
	dataDirs = str.split(dataDirsOut.stdout.decode('utf-8').rstrip(), '\n')

	if opts.verbose:
		print('Collecting data from {} to {}'.format(pars['data_sYear'], pars['data_eYear']))

	# Get stage 1 dirs --> /FM/FlightData/2A/YearOfAcquisition/
	for dir_st1 in dataDirs:
		
		if "20" in dir_st1 and "statistics" not in dir_st1:

			# Extract year
			year = int(dir_st1[dir_st1.rfind('/')+1:])

			if year < pars['data_sYear'] or year > pars['data_eYear']:
				continue

			if year not in years:
				years.append(year)
				year_data_idx = len(years)-1
				counters.append(0)

			getDataDirsCommand = 'xrdfs {} ls {}'.format(pars['farmAddress'], dir_st1)
			if opts.verbose:
				print('Executing XRDFS command: {}'.format(getDataDirsCommand))
			dataDirsOut = subprocess.run(getDataDirsCommand, shell=True, check=True, stdout=subprocess.PIPE)
			dataDirs_st1 = str.split(dataDirsOut.stdout.decode('utf-8').rstrip(), '\n')

			# Get stage 2 dirs --> /FM/FlightData/2A/YearOfAcquisition/MonthOfAcquisition
			for dir_st2 in dataDirs_st1:
				getDataDirsCommand = 'xrdfs {} ls {}'.format(pars['farmAddress'], dir_st2)
				if opts.verbose:
					print('Executing XRDFS command: {}'.format(getDataDirsCommand))
				dataDirsOut = subprocess.run(getDataDirsCommand, shell=True, check=True, stdout=subprocess.PIPE)
				dataDirs_st2 = str.split(dataDirsOut.stdout.decode('utf-8').rstrip(), '\n')

				# Get ROOT data file
				for data_elm in dataDirs_st2:
					if data_elm.endswith('.root') and "data_photon" not in data_elm:
						dList.append(data_elm)
						counters[year_data_idx] += 1

	if opts.verbose:
		print('{} data files have been read...'.format(sum(counters)))
		for year_idx, year in enumerate(years):
			print('{} data files found in {} folder'.format(counters[year_idx], year))

	if opts.output:
		data_list_path = opts.output
	else:
		data_list_path =  "dataFileList.txt"
	with open(data_list_path, "w") as outList:
		for elm in dList:
			outList.write(pars['farmAddress'] + elm + "\n")