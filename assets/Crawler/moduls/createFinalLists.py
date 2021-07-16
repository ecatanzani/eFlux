def createOutSimuFile(opts, simu_eMin, simu_eMax, particle):
	if not opts.input:
		finalPath = "simuFileList_"
		finalPath += particle
		finalPath += "_"
		finalPath += str(simu_eMin)
		finalPath += "_"
		finalPath += str(simu_eMax)
		finalPath += ".txt"
		return finalPath
	else:
		return opts.output

def createOutSetSimuFile(dSetName):
	return "simuFileList_" + dSetName