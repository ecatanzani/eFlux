Electron Flux
=======

This is a C++ software that computes the all-electron spectrum in a given energy range using DAMPE flight data.
 
Practical information:

* All the information contained in the flight data are incapsulated using a **TTree**, according to **ROOT** framework.
* **XROOTD** will be the default protocol to aceess remote data, both flight and MC. This protocol permits to execute the software indipendently on the target farm, avoidint to specify custom paths for each one of them or to move large amounts of data.

Software usage:

```markdown

Usage: 

 -h  --help                                                   Prints this help
 -i  --input          <path_to_input_DATA_TTree>          (*) Input data TTree - flux calculation only
 -o  --output         <path_to_output_TFile>                  Output ROOT TFile
 -d  --outputDir      <path_to_output_TFile_dir>              Output ROOT TFile directory
 -t  --lvtime         <live_time-value>                   (*) DAMPE live-time 
 -a  --acceptance     <path_to_MC_list>                   (*) Acceptance calculation
 -c  --collect        <path_to_complete_histo>            (*) Generate TGraph from final histo
 -f  --flux           <path_to_DATA_list_dir>             (*) Flux calculation
 -g  --geometry       <path_to_acceptance_ROOT_file>      (*) Acceptance file - flux calculation only
 -v  --verbose                                                Verbose output
 -p  --pedantic                                               Pedantic output
```

In order to activate the flux computation, the **--flux** flux must me used.

Aceptance calculation
---------------

The flux computation requires an acceptance, an extrnal ROOT file, as input parameter.
The software has an acceptance computing function already built-in, whose parameters can be set through the config file (below an example). This functionality is available proving the **--acceptance** flag.

The acceptance calculation requires, as input parameter, the full path of a list of **RECO** MC files, accordingly to the **XROOTD** format.

Here a default config file for the acceptance computation.

```markdown
--------- Acceptance Config File ---------

-- Cuts

---- CUT VARIABLE ----       ---- VALUE ----      ---- UNIT ----

min_event_energy                    1                   GeV
max_event_energy                    10000               GeV
energy_lRatio                       0.35
shower_axis_delta                   280                 mm
max_rms_shower_width                100                 mm
track_X_clusters                    4
track_Y_clusters                    4
track_missingHit_X                  1
track_missingHit_Y                  1
STK_BGO_delta_track                 10                  deg
STK_BGO_delta_position              40                  mm
xtrl                                8.5
STK_PSD_delta_position              40                  mm
PSD_bar_min_energy_release          0.5                 MeV

---- CUT ----                ---- ACTIVE ----

BGO_fiducial                        YES
nBarLayer13                         NO
maxRms                              NO
track_selection                     NO
xtrl_selection                      NO
psd_charge                          NO

-- Simu Params

generation_vertex_radius            1.381976597885342   m

```

This config file shows the default value for the parameters used to compute the acceptance, such as the energy range or the XTRL value. Each of the parameter can be changed in order to modify the acceptance computation. This file is parsed at the beginning of the code execution; once the program starts, future modification to the config file will NOT be considered in the current instance of the software.

Here a description of the event selection parameters:

* **min_event_energy**: minimum value for the energy range (default value set to *1 GeV*)
* **max_event_energy**: maximum value for the energy range (default value set to *10 TeV*)
* **energy_lRatio**: ratio of the energy released in a certain layer respect to the total calorimeter energy (default value set to *0.35*)
* **shower_axis_delta**: distance, expressed in mm, between the projection points (x,y) of the BGO shower on the TOP and BOTTOM layers of the calorimeter respect to its center (default value set to *280 mm*)
* **max_rms_shower_width** : maximum value of the shower width of all layers with energy more than 1% of the total calorimeter energy (default value set to *100 mm*)
* **track_X_clusters**: number of required X clusters for a track (default value set to *4*)
* **track_Y_clusters**: number of required Y clusters for a track (default value set to *4*)
* **track_missingHit_X**: maximum number of X missing points for a track (default value set to *1*)
* **track_missingHit_Y**: maximum number of Y missing points for a track (default value set to *1*)
* **STK_BGO_delta_track**: maximum mangular distance between the direction of the reconstructed BGO shower and the STK track direction (default value set to *10 deg*)
* **STK_BGO_delta_position**: maximum value of the distance between the extrapolated positions of the track and the shower to the top of the BGO (default value set to *40 mm*)
* **xtrl**: value used to discriminate between electrons(positrons) and hadrons (default value set to *8.5*)
* **STK_PSD_delta_position**: maximum value of the distance between the hit position of the PSD cluster seed strip and the extrapolated track position (default value set to *40 mm*)
* **PSD_bar_min_energy_release**: value of the minimum required energy release into the PSD bard (default value set to *0.5 MeV*)

All DAMPE simulations are based on an HALF-sphere with the detector on its center; all the simulated particles go downwards (exept for the data sets marked as **BACKENTERING**). The radius of such a sphere may vary but its default value, represented by the parameter **generation_vertex_radius** is set to 1.38 m. This value canbe found both on the **mac** file used for the simulation or studying the **DmpEventSimuPrimaries/pv_x** and **DmpTruthTrajectoriesCollection/start_x** distributions on the **RECO** files.

The last section of the config file permits to activate or deactivate some event cuts.
The geometric cut is activated by default and will be performed each time the code executes.

The list of the provided cuts is the following:

* **BGO_fiducial**
* **nBarLayer13**
* **maxRms**
* **track_selection**
* **xtrl_selection**
* **psd_charge**

**BGO_fiducial** cut ensures that all the events selected are well contained into the BGO volume (this cut is active by default).
**nBarLayer13** and **maxRms** permit to remove lateral and large showering events.
**track_selection** is used to select events with a good STK track associated.
**xtrl_selection** is used to select electrons (positrons) from a background of hadrons or low energy particles that mimic thei behaviour.

All the accessory cuts can be activated/deactivated using the flags **YES**/**NO**.

### Crawler

Crawler is an **Asset** software used to easily create the input MC **RECO** file list for the acceptance computation.

**Python 3** and **XROOTD** are required to execute Crawler.

Software usage:

```markdown
reco MC files explorer

optional arguments:
  -h, --help            		show this help message and exit
  -i INPUT, --input INPUT		json data storage
  -f FARM, --farm FARM			MC ROOT input farm
  -o OUTPUT, --output OUTPUT	output file list
  -v, --verbose         		run in high verbosity mode
```

Crawler computes the MC file list as output (the **--output** flag is optional and, if not specified, a default value will be automatically assigned accordingly to the energy range chosen). 

The software uses a config file to produce the final data list.

This is an example of the configuration file with some default values.

```markdown
---DataSets default params

farmAddress         root://farm_address//
simu_XRDFS_path     /MC/reco/
geometry            6r0p0
simu_eMin           1
simu_eMax           10000
particle            e
jSet                dataSets.json
```

* **farmAddress**: the address of the farm where the data are located (default value set to *root://farm_address//*)
* **simu_XRDFS_path**: local path for the reconstructed MC data (default value set to */MC/reco/*)
* **geometry**: release of the detector geometry (default value set to *6r0p10*)
* **simu_eMin**: minimum value of the data-set energy (default value set to *1 GeV*)
* **simu_eMax**: maximum value of the data-set energy (default value set to *10 TeV*)
* **particle**: particle ID (default value set to *e* for *electrons*)
* **jSet**: name of the local json file where all the MC data-sets are collected according to the previous parameters (default value set to *dataSets.json*, but external json files can be used, through the **-input** flag).

Crawler parses the configuration file and, accordingly to the parameters value, queries *jSet* to retieve the proper MC reco data list.

As an example, this is a real usage case of Crawler according to the previous configuration file:

```markdown
❯ python3 explorer.py -v
xrdfs root://farm_address// ls /MC/reco/v6r0p0 allElectron-v6r0p0_100GeV_10TeV
xrdfs root://farm_address// ls /MC/reco/v6r0p0/allElectron-v6r0p0_100GeV_10TeV-p2
xrdfs root://farm_address// ls /MC/reco/v6r0p0/allElectron-v6r0p0_100GeV_10TeV-p3
xrdfs root://farm_address// ls /MC/reco/v6r0p0/allElectron-v6r0p0_1GeV_100GeV
xrdfs root://farm_address// ls /MC/reco/v6r0p0/allElectron-v6r0p0_1GeV_15GeV
0:00:08.744187
```

Crawler creates a txt list with default name accordingly to the energy range:

```markdown
❯ ls -ll
total 12320
-rw-r--r--  1 enrico  staff      246  1 Mag 11:35 crawlerConfig.txt
-rw-r--r--  1 enrico  staff     1672 20 Feb 23:21 dataSets.json
-rw-r--r--  1 enrico  staff     3963  9 Mar 15:21 explorer.py
-rw-r--r--  1 enrico  staff  5786558  1 Mag 11:42 simuFileList_e_1.0_10000.0.txt
```

```markdown
❯ head simuFileList_e_1.0_10000.0.txt
root://farm_address///MC/reco/v6r0p0/allElectron-v6r0p0_100GeV_10TeV/allElectron-v6r0p0_100GeV_10TeV.noOrb.000001.reco.root
root://farm_address///MC/reco/v6r0p0/allElectron-v6r0p0_100GeV_10TeV/allElectron-v6r0p0_100GeV_10TeV.noOrb.000002.reco.root
root://farm_address///MC/reco/v6r0p0/allElectron-v6r0p0_100GeV_10TeV/allElectron-v6r0p0_100GeV_10TeV.noOrb.000003.reco.root
root://farm_address///MC/reco/v6r0p0/allElectron-v6r0p0_100GeV_10TeV/allElectron-v6r0p0_100GeV_10TeV.noOrb.000004.reco.root
root://farm_address///MC/reco/v6r0p0/allElectron-v6r0p0_100GeV_10TeV/allElectron-v6r0p0_100GeV_10TeV.noOrb.000005.reco.root
root://farm_address///MC/reco/v6r0p0/allElectron-v6r0p0_100GeV_10TeV/allElectron-v6r0p0_100GeV_10TeV.noOrb.000006.reco.root
root://farm_address///MC/reco/v6r0p0/allElectron-v6r0p0_100GeV_10TeV/allElectron-v6r0p0_100GeV_10TeV.noOrb.000007.reco.root
root://farm_address///MC/reco/v6r0p0/allElectron-v6r0p0_100GeV_10TeV/allElectron-v6r0p0_100GeV_10TeV.noOrb.000008.reco.root
root://farm_address///MC/reco/v6r0p0/allElectron-v6r0p0_100GeV_10TeV/allElectron-v6r0p0_100GeV_10TeV.noOrb.000009.reco.root
root://farm_address///MC/reco/v6r0p0/allElectron-v6r0p0_100GeV_10TeV/allElectron-v6r0p0_100GeV_10TeV.noOrb.000010.reco.root
```

Using the default *v6r0p0* geometry, the following data-sets will be available:

* allElectron-v6r0p0_1GeV_15GeV
* allElectron-v6r0p0_1GeV_100GeV
* allElectron-v6r0p0_100GeV_10TeV
* allElectron-v6r0p0_100GeV_10TeV-p2
* allElectron-v6r0p0_100GeV_10TeV-p3

### Energy binning

This is an **Asset** software to create a logaritmic energy binning for the flux and acceptance calculation. 

This software is particularly usefull if multiple data-sets are needed to cover the analysis energy range; in this case the bins should be computed in order to match the energy range of each data-set. 

This software is written in **Python** and requires the **numpy** package.

Here the usage:

```markdown
Binning Finder

optional arguments:
  -h, --help                  show this help message and exit
  -i INPUT, --input INPUT     json data storage
  -o OUTPUT, --output OUTPUT  Output binning text file
  -v, --verbose               run in high verbosity mode
```
The software reads a json configuration file (an external one can also be used, through the **--input** flag). This file stores the parameters used during the energy binning building. Here an example, showing the default configuration:

```markdown
{
    "eMin": 1,
    "eMax": 10000,
    "bins": 20,
    "junctions": []
}
```

* **eMin**: minimum energy value used in the binning (default value set to *1 GeV*)
* **eMax**: maximum energy value used in the binning (default value set to *10 TeV*)
* **bins**: number of bins (default value set to *20*)
* **junctions**: *list* of the energy edges of the different data-sets (*empty* by default)

The software may produces an output text file containing the binning, whose path needs to be specified using the **--output** flag. If the flag is not specified the binning is not written to disk.

```markdown
❯ python getEnergyBinning.py -o myBinning.txt -v

**** Binning settings:
eMin: 1
eMax: 10000
Number of bins: 20
****

Binning found!
('Number of bins: ', 20)
****

[  1.00000000e+00   1.58489319e+00   2.51188643e+00   3.98107171e+00
   6.30957344e+00   1.00000000e+01   1.58489319e+01   2.51188643e+01
   3.98107171e+01   6.30957344e+01   1.00000000e+02   1.58489319e+02
   2.51188643e+02   3.98107171e+02   6.30957344e+02   1.00000000e+03
   1.58489319e+03   2.51188643e+03   3.98107171e+03   6.30957344e+03
   1.00000000e+04]

0:00:00.004979
```

```markdown
❯ ls -ll
total 24
-rw-r--r--  1 enrico  staff    74  1 Mag 14:59 config.json
-rw-r--r--  1 enrico  staff  3688 28 Apr 10:37 getEnergyBinning.py
-rw-r--r--  1 enrico  staff   250  1 Mag 15:13 myBinning.txt
```

```markdown
❯ head myBinning.txt
1.0
1.58489319246
2.51188643151
3.98107170553
6.3095734448
10.0
15.8489319246
25.1188643151
39.8107170553
63.095734448
```

**IMPORTANT NOTE:** if the **junctions** list is empty the number of bins reported in the json config file is conserved, while in the other cases the software may modify the number of bins (at maximum *5*) in order to match the energy requirement and found the best binning.

Support
=======

Having trouble with this code ? Check out the [documentation](https://ecatanzani.github.io/eFlux/)

Prerequisites
=======

- C++14 (5.3.1 required)
- DAMPESW: Event package. The whole DAMPESW package is NOT required to build and run this code.
- ROOT (5.34.36 required)

Build Tests
=======

| CentOS |
|:--:|
| ![CentOS](https://github.com/ecatanzani/eFlux/workflows/CentOS%20-%20DAMPE%20framework/badge.svg) |
