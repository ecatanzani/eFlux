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
    opt.addUsage(" -a  --acceptance     <Mpath_to_MC/Data_dir>          Acceptance calculation");
    opt.addUsage(" -v  --verbose                                        Verbose output   ");
    opt.addUsage("");
    
    opt.setFlag("help",'h');
    opt.setOption("input",'i');
    opt.setOption("output",'o');
    opt.setOption("outputDir",'d');
    opt.setOption("lvtime", 't');
    opt.setFlag("verbose",'v');

    opt.processCommandArgs(argc,argv);
    
    /*
        Input variables

    */

    std::string inputPath;
    std::string outputPath;
    std::string accInputPath;
    stringstream str_lvTime;
    bool verbose = false;
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
    if(opt.getFlag("acceptance") || opt.getFlag('a'))
    {
        myAcceptance = true;
        accInputPath = opt.getValue('a');
    }
    if(opt.getFlag("verbose") || opt.getFlag('v'))
        verbose = opt.getFlag('v');

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

    if(argc>=5)
        eCore(
                inputPath,
                outputPath,
                verbose,
                lvTime,
                myAcceptance,
                accInputPath,
                opt
            );
    else
    {
        std::cerr << "\n\t !!! Hey mate, insert all required fields ...\n\n";opt.printUsage();
    }

#endif


    return 0;
}