Electron Flux
=======

This software computes the all-electron spectrum in a given energy range using DAMPE flight data.

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

### Support

Having trouble with this code ? Check out the [documentation](https://ecatanzani.github.io/eFlux/)

### Prerequisites

- C++14 (5.3.1 required)
- DAMPESW: Event package. To BUILD and RUN this code not the whole DAMPESW package is required
- ROOT (5.34.36 required)

### Build Tests

| CentOS |
|:--:|
| ![CentOS](https://github.com/ecatanzani/eFlux/workflows/CentOS%20-%20DAMPE%20framework/badge.svg) |

