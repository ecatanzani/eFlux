#include "myHeader.h"

#include <sstream>

int main(int argc,char* argv[])
{
    
    AnyOption opt;
    
    opt.addUsage("Usage: ");
    opt.addUsage("");
    opt.addUsage(" -h  --help                                           Prints this help ");
    opt.addUsage(" -i  --input          <path_to_input_TTree>       (*) Input data TTree ");
    opt.addUsage(" -o  --output         <path_to_output_TFile>          Output ROOT TFile");
    opt.addUsage(" -d  --outputDir      <path_to_output_TFile_dir>      Output ROOT TFile directory");
    opt.addUsage(" -t  --lvtime         <live-time-value>           (*) DAMPE live-time  ");
    opt.addUsage(" -a  --acceptance     <path_to_MC/Data_list_dir>      Acceptance calculation");
    opt.addUsage(" -v  --verbose                                        Verbose output   ");
    opt.addUsage(" -p  --pedantic                                       Pedantic output   ");
    opt.addUsage("");
    
    opt.setFlag("help",'h');
    opt.setOption("input",'i');
    opt.setOption("output",'o');
    opt.setOption("outputDir",'d');
    opt.setOption("lvtime", 't');
    opt.setOption("acceptance",'a');
    opt.setFlag("verbose",'v');
    opt.setFlag("pedantic",'p');

    opt.processCommandArgs(argc,argv);
    
    /*
        Input variables

    */

    std::string inputPath;
    std::string outputPath;
    std::string accInputPath;
    stringstream str_lvTime;
    bool verbose = false;
    bool pedantic = false;
    bool myAcceptance = false;
    unsigned int lvTime = 0;
    
   
    if(!opt.hasOptions()) 
        opt.printUsage();

    if(opt.getFlag("help") || opt.getFlag('h'))
        opt.printUsage();
    if(opt.getValue("input") || opt.getValue('i'))
        inputPath = opt.getValue('i');
    if(opt.getValue("output") || opt.getValue('o'))
        outputPath = opt.getValue('o');
    if(opt.getValue("outputDir") || opt.getValue('d'))
        outputPath = opt.getValue('d');
    if(opt.getValue("lvtime") || opt.getValue('t'))
    {
        str_lvTime << opt.getValue('t');
        str_lvTime >> lvTime;
    }
    if(opt.getValue("acceptance") || opt.getValue('a'))
    {
        myAcceptance = true;
        accInputPath = opt.getValue('a');
    }
    if(opt.getFlag("verbose") || opt.getFlag('v'))
        verbose = opt.getFlag('v');
    if(opt.getFlag("pedantic") || opt.getFlag('p'))
        pedantic = opt.getFlag('p');

#ifdef DEBUG

    if(argc>=5)
    {
        if(chechFlags(opt,inputPath,outputPath,lvTime))
        {
            eCore(inputPath,outputPath,verbose,lvTime,opt);
        }
    }
    else
    {
        std::cerr << "\n\t !!! Hey mate, insert all required fields ...\n";opt.printUsage();
    }

#else
    eCore(
            inputPath,
            outputPath,
            verbose,
            pedantic,
            lvTime,
            myAcceptance,
            accInputPath,
            opt
        );
#endif

    return 0;
}