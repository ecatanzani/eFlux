def getTotalEvents(opts):
    totEvents = 0
    for file in os.listdir(opts.input):
        fNpath = str(opts.input) + "/" + str(file)
        data = np.load(fNpath,'r')
        nRows = np.size(data,0)
        totEvents += nRows
    return totEvents

def eventProcess(opts,evIdx,totEvents,kStep=1):
    kBreak = False
    if opts.verbose:
        if opts.debug:
            if (evIdx%((opts.debug)/10))==0 and evIdx!=0:
                print('\tProcessed event {} of {}'.format(evIdx,totEvents))
            if evIdx==(opts.debug-1):
                print('\tProcessed event {} of {}'.format(evIdx+1,totEvents))
                kBreak = True 
        else:
            if (evIdx%kStep)==0 and evIdx!=0:
                print('\tProcessed event {} of {}'.format(evIdx,totEvents))
    return kBreak
