# DAMPE - Live-Time computation
This utility is ised to compute the live-time of the DAMPE detector.
The live-time utility takes the following inputs:
- start 2A input file
- end 2A input file
- verbose (optional)

The input files can be provided both with the absolute/relative file system path or using XROOTD protocol (in this case the path must include also the host server).

The live-time utility needs the starting second and the final second of the acquisition period in order to compute the live-time. Both are automatically extracted from the 2A files provided using the DAMPEChain tool.

## Usage

```
usage: Usage: lvTimeCounter.py [options]

DAMPE live-time facility

optional arguments:
  -h, --help            show this help message and exit
  -s START_FILE, --start START_FILE
                        start 2A DAMPE data file
  -e END_FILE, --end END_FILE
                        end 2A DAMPE data file
  -v, --verbose         run in high verbosity mode
```

## Usage example

```
[user@server assets]$ python lvTimeCounter.py -s xrootd_input_start_file.root -e xrootd_input_end_file.root -v

**************************************************
      Offline software of DAMPE (DMPSW)
      version:  6.1(alpha)
**************************************************
  [Service manager]	Adding element: DmpIOSvc
  
Adding xrootd_input_start_file.root to the chain...
Adding xrootd_input_end_file.root to the chain...

DmpEvent::GetHead::Creating-DmpEvent-As-Singleton
DmpEvent::InitHead
DmpChain::GetDmpEvent::RecognizedAs::Flight

Computing LiveTime in the following interval: [94608354, 278986240]

Error in <TMySQLServer::TMySQLServer>: Code: 2003  Msg: Can't connect to MySQL server on '192.168.1.154' (110)
Error in <TMySQLServer::TMySQLServer>: Code: 2003  Msg: Can't connect to MySQL server on '172.16.0.2' (110)

LiveTime: 140143948.911
```

The example above computes the live-time relative to the period 01/01/2016 - 03/11/2021.