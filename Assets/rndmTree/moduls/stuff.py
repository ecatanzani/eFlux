def eventProcess(opts,evIdx,totEvents,kStep=1):
    if opts.verbose:
        if (evIdx%kStep)==0 and evIdx!=0:
            print('\tProcessed event {} of {}'.format(evIdx,totEvents))
