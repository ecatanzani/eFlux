#include "myHeader.h"
#include "anyoption.h"

int main(int argc,char* argv[])
{
    
    AnyOption opt;
    
    opt.addUsage("Usage: ");
    opt.addUsage("");
    opt.addUsage(" -h  --help                           Prints this help ");
    opt.addUsage(" -i  --input  <path_to_input_TTree>   Input data TTree ");
    opt.addUsage(" -o  --output <path_to_output_TFile>  Output ROOT TFile");
    opt.addUsage(" -v  --verbose                        Verbose output   ");
    opt.addUsage("");
    
    opt.setFlag("help",'h');
    opt.setOption("input",'i');
    opt.setOption("output",'o');
    opt.setFlag("verbose",'v');

    opt.processCommandArgs(argc,argv);
    
    /*
        Input variables

    */

    std::string inputPath;
    std::string outputPath;
    bool verbose = false;
    
   
    if(!opt.hasOptions()) 
        opt.printUsage();

    if(opt.getFlag("help") || opt.getFlag('h'))
        opt.printUsage();
    if(opt.getValue("input") || opt.getValue('i'))
        inputPath = opt.getValue('i');
    if(opt.getValue("output") || opt.getValue('o'))
        inputPath = opt.getValue('o');
    if(opt.getFlag("verbose") || opt.getFlag('v'))
        verbose = opt.getFlag('v');
    
    eCore(inputPath,outputPath,verbose);

    return 0;
}