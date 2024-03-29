import subprocess
import json
import os
import sys
from createFinalLists import createOutSimuFile


def createMClist(pars, opts):
	# Parse masks sets
	maskInputFile = os.getcwd()
	maskInputFile = maskInputFile[0:maskInputFile.find('/Crawler')+8] + "/setMask.conf"
	try:
		dListMask = open(maskInputFile).read().splitlines()
	except OSError:
			print('ERROR reading input dataSet mask file in: {}'.format(maskInputFile))
			raise
	# Parse data json
	with open(pars['jSet'], 'rb') as set_:
		dataSets = json.load(set_)
		# Create out file list
		dList = []
		setPowerLowConfig = []
		for elm in dataSets:
			if elm["dSet"] in dListMask:
				continue
			if elm["geometry"] == pars['geometry']:
				if elm["particle"] == pars['particle']:
					if int(elm["eMin"]) >= pars['simu_eMin']:
						if int(elm["eMax"]) <= pars['simu_eMax']:
							tmpSetPath = pars['simu_XRDFS_path'] + "v" + \
								pars['geometry'] + "/" + elm["dSet"]
							dList.append(tmpSetPath)
							setPowerLowConfig.append(
								[elm['dSet'], elm['eMin'], elm['eMax'], elm['pSpectrum']])

	finalListPath = createOutSimuFile(
		opts, pars['simu_eMin'], pars['simu_eMax'], pars['particle'])
	if opts.verbose:
		print('Writing output MC file list: {}'.format(finalListPath))
	with open(finalListPath, "w") as outList:
		for set_ in dList:
			tmpCommand = 'xrdfs {} ls {}'.format(pars['farmAddress'], set_)
			if opts.verbose:
				print('Executing XRDFS command: {}'.format(tmpCommand))
			tmpOut = subprocess.run(
				tmpCommand, shell=True, check=True, stdout=subprocess.PIPE)
			setList = str.split(tmpOut.stdout.decode('utf-8').rstrip(), '\n')
			for file in setList:
				outPath = pars['farmAddress'] + file + "\n"
				outList.write(outPath)
