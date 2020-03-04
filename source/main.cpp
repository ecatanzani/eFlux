#include "myHeader.h"
#include "acceptance.h"

#include <sstream>

int main(int argc, char **argv)
{

    AnyOption opt;

    opt.addUsage("Usage: ");
    opt.addUsage("");
    opt.addUsage(" -h  --help                                                   Prints this help");
    opt.addUsage(" -i  --input          <path_to_input_DATA_TTree>          (*) Input data TTree - flux calculation only");
    opt.addUsage(" -o  --output         <path_to_output_TFile>                  Output ROOT TFile");
    opt.addUsage(" -d  --outputDir      <path_to_output_TFile_dir>              Output ROOT TFile directory");
    opt.addUsage(" -t  --lvtime         <live_time-value>                   (*) DAMPE live-time ");
    opt.addUsage(" -a  --acceptance     <path_to_MC_list_dir>               (*) Acceptance calculation");
    opt.addUsage(" -f  --flux           <path_to_DATA_list_dir>             (*) Flux calculation");
    opt.addUsage(" -g  --geometry       <path_to_acceptance_ROOT_file>      (*) Acceptance file - flux calculation only");
    opt.addUsage(" -v  --verbose                                                Verbose output");
    opt.addUsage(" -p  --pedantic                                               Pedantic output");
    opt.addUsage("");

    opt.setFlag("help", 'h');
    opt.setOption("input", 'i');
    opt.setOption("output", 'o');
    opt.setOption("outputDir", 'd');
    opt.setOption("lvtime", 't');
    opt.setOption("acceptance", 'a');
    opt.setOption("flux", 'f');
    opt.setOption("geometry", 'g');
    opt.setFlag("verbose", 'v');
    opt.setFlag("pedantic", 'p');

    opt.processCommandArgs(argc, argv);

    /*
        Input variables

    */

    std::string inputPath;
    std::string outputPath;
    std::string accInputPath;
    std::string fluxInputPath;
    stringstream str_lvTime;

    bool verbose = false;
    bool pedantic = false;
    bool myAcceptance = false;
    bool myFlux = false;
    unsigned int lvTime = 0;

    if (!opt.hasOptions())
        opt.printUsage();

    if (opt.getFlag("help") || opt.getFlag('h'))
        opt.printUsage();
    if (opt.getValue("input") || opt.getValue('i'))
        inputPath = opt.getValue('i');
    if (opt.getValue("output") || opt.getValue('o'))
        outputPath = opt.getValue('o');
    if (opt.getValue("outputDir") || opt.getValue('d'))
        outputPath = opt.getValue('d');
    if (opt.getValue("lvtime") || opt.getValue('t'))
    {
        str_lvTime << opt.getValue('t');
        str_lvTime >> lvTime;
    }
    if (opt.getValue("acceptance") || opt.getValue('a'))
    {
        myAcceptance = true;
        accInputPath = opt.getValue('a');
    }
    if (opt.getValue("flux") || opt.getValue('f'))
    {
        myFlux = true;
        fluxInputPath = opt.getValue('f');
    }
    if (opt.getValue("geometry") || opt.getValue('g'))
        accInputPath = opt.getValue('g');
    if (opt.getFlag("verbose") || opt.getFlag('v'))
        verbose = opt.getFlag('v');
    if (opt.getFlag("pedantic") || opt.getFlag('p'))
        pedantic = opt.getFlag('p');

    if (myAcceptance)
        computeAcceptance(
            accInputPath,
            verbose,
            pedantic,
            outputPath,
            opt);

    if (myFlux)
        eCore(
            inputPath,
            outputPath,
            verbose,
            pedantic,
            lvTime,
            accInputPath,
            opt);

    return 0;
}