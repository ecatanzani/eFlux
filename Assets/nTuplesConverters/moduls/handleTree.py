from tqdm import tqdm

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

def fillTreeChain(
                opts,
                dChain,
                nevents,
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
    
    if opts.verbose:
        print("Filling output Tree...")
        if opts.debug:
            print('Debug mode activated... max mumber of events {}'.format(opts.debug))

    for idx,event in tqdm(enumerate(dChain)):
        t_totalEnergyCorr[0] = event.tt_bgoTotalEcorr_GeV/1000.
        t_totalEnergy[0] = event.tt_bgoTotalE_GeV/1000.
        t_xtrl[0] = event.tt_Xtrl
        t_satPositionX[0] = event.tt_sat_position_x
        t_satPositionY[0] = event.tt_sat_position_y
        t_satPositionZ[0] = event.tt_sat_position_z
        t_satVelocityX[0] = event.tt_sat_velocity_x
        t_satVelocityY[0] = event.tt_sat_velocity_y
        t_satVelocityZ[0] = event.tt_sat_velocity_z
        myTree.Fill()
        if opts.verbose:
            if opts.debug:
                if (idx%((opts.debug)/10))==0 and idx!=0:
                    print('\tProcessed event {} of {}'.format(idx,nevents))
                if idx==(opts.debug-1):
                    print('\tProcessed event {} of {}'.format(idx+1,nevents))
                    break
            else:
                if (idx%kStep)==0 and idx!=0:
                    print('\tProcessed event {} of {}'.format(idx,nevents))
            

def fillTree(
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
        print("Filling output Tree...")
    
    for file in tqdm(os.listdir(opts.input)):
        fNpath = str(opts.input) + "/" + str(file)
        data = np.load(fNpath,'r')
        nRows = np.size(data,0)
        for evIdx in range(0,nRows):
            t_totalEnergyCorr[0] =  data[evIdx,44]/1000.
            t_totalEnergy[0] = data[evIdx,47]/1000.
            t_xtrl[0] = data[evIdx,49]
            ''' This still needs to be impleented into the pickle input files

            t_satPositionX[0] = 
            t_satPositionY[0] = 
            t_satPositionZ[0] =
            t_satVelocityX[0] = 
            t_satVelocityY[0] = 
            t_satVelocityZ[0] = 

            '''
            myTree.Fill()
        if opts.verbose:
            if (idx%kStep)==0:
                print('\tProcessed event {} of {}'.format(idx,nevents))
