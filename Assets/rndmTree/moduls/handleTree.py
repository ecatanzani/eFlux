import sys
from ROOT import TTree,TFile,TRandom3

def branchTree(
                opts,
                myTree,
                t_totalEnergyCorr,
                t_totalEnergy,
                t_xtrl,
                t_satPositionX,
                t_satPositionY,
                t_satPositionZ,
                t_satVelocityX,
                t_satVelocityY,
                t_satVelocityZ
                ):

    if opts.verbose:
        print("Branching output TTree...")
    myTree.Branch( 'eReco', t_totalEnergy, 'eReco/F' )
    myTree.Branch( 'eRecoCorr', t_totalEnergyCorr, 'eRecoCorr/F' )
    myTree.Branch( 'xtrl', t_xtrl, 'xtrl/F' )
    myTree.Branch( 'posizionX', t_satPositionX, 'posizionX/F' )
    myTree.Branch( 'posizionY', t_satPositionY, 'posizionY/F' )
    myTree.Branch( 'posizionZ', t_satPositionZ, 'posizionZ/F' )
    myTree.Branch( 'velocityX', t_satVelocityX, 'velocityX/F' )
    myTree.Branch( 'velocityY', t_satVelocityY, 'velocityY/F' )
    myTree.Branch( 'velocityZ', t_satVelocityZ, 'velocityZ/F' )

def shuffleTree(
                opts,
                myTree,
                t_totalEnergyCorr,
                t_totalEnergy,
                t_xtrl,
                t_satPositionX,
                t_satPositionY,
                t_satPositionZ,
                t_satVelocityX,
                t_satVelocityY,
                t_satVelocityZ,
                kStep
                ):
    
    from stuff import eventProcess

    if opts.verbose:
        print("Filling output TTree...")
    
    #Reading input TFile
    inFile = TFile(opts.input,"READ")
    
    if not inFile.IsOpen():
        print('Error opening input TFile: {}'.format(opts.input))
        sys.exit()
    
    #Extracting input TTree
    dTree = inFile.Get("collectionTree")
    
    if opts.verbose:
        print('Collected {} events from input TTree'.format(dTree.GetEntries()))
        print('{} events will be randomly extracted'.format(opts.events))

    prEvents = 0                #Number of processed events
    gRandom = TRandom3()
    #Filling output TTre
    idxList = []                #List of shuffled indexes - avoid repetitions

    while prEvents < opts.events:
        
        rndmIdx = gRandom.Integer(dTree.GetEntries())
        if len(idxList) == 0:
            idxList.append(rndmIdx)
        else:
            if rndmIdx not in idxList:
                idxList.append(rndmIdx)
            else:
                continue
        
        #Read TTree event
        dTree.GetEntry(rndmIdx)
        
        #Read branches
        t_totalEnergyCorr[0] = dTree.eRecoCorr
        t_totalEnergy[0] = dTree.eReco
        t_xtrl[0] = dTree.xtrl
        t_satPositionX[0] = dTree.posizionX
        t_satPositionY[0] = dTree.posizionY
        t_satPositionZ[0] = dTree.posizionZ
        t_satVelocityX[0] = dTree.velocityX
        t_satVelocityY[0] = dTree.velocityY
        t_satVelocityZ[0] = dTree.velocityZ

        #Fill output TTree
        myTree.Fill()

        eventProcess(opts,prEvents+1,dTree.GetEntries(),kStep)
        prEvents+=1
    
    #Close input TFile
    inFile.Close()



