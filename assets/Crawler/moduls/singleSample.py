import subprocess
import json
from createFinalLists import createOutSetSimuFile


def listMCsample(pars, opts):
	dSetName = opts.spectrum
	_setPath = pars['simu_XRDFS_path'] + \
		dSetName[dSetName.find('v'):dSetName.find('_')] + "/" + dSetName
	tmpCommand = 'xrdfs {} ls {}'.format(pars['farmAddress'], _setPath)
	if opts.verbose:
		print('Getting list of files... {}'.format(tmpCommand))
	tmpCommandOut = subprocess.run(tmpCommand, shell=True,
							check=True, stdout=subprocess.PIPE)
	_setList = str.split(tmpCommandOut.stdout.decode('utf-8').rstrip(), '\n')
	if opts.verbose:
		print('Found {} reco ROOT files for the dataSet'.format(len(_setList)))
	
	listOutPath = createOutSetSimuFile(dSetName)
	with open(listOutPath, "w") as outList:
		for elm in _setList:
			tmpFinalPath = pars['farmAddress'] + elm + "\n"
			outList.write(tmpFinalPath)    

