#!/usr/bin/env python
from ROOT import * 

def main(args=None):
    gSystem.Load("libDmpService.so") 
    a=DmpSvcLiveTime.GetInstance() 
    #01/01/16-31/10/17 
    b =a.GetLiveTime(94608354,152496303); 
    
if __name__ == '__main__' :
    main()