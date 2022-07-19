#ifndef INPUT_H
#define INPUT_H

// for handling command-line arguments and logging for
// aaocean terminal app

#include "log.h"
#include "optionparser.h"

class input
{
public:
    unsigned int resolution;
    unsigned int spectrum;
    unsigned int seed;
    float oceanScale;
    float oceanDepth;
    float surfaceTension;
    float velocity;
    float smooth;
    float windDir;
    int windAlign;
    float reflectedWaves;
    float waveSpeed;
    float waveHeight;
    float waveChop;
    float repeatTime;
    char outputFolder[512];
    char postfix[512];

    int startFrame;
    int endFrame;
    float fps;

    int accumulateFoam;

    int logLevel;
    char helpMsg[2048];

    input();
    void help();
};

input::input()
{
    resolution      = 64;
    spectrum        = 0;
    seed            = 1;
    oceanScale      = 100.f;
    oceanDepth      = 1000.f;
    surfaceTension  = 0.0f;
    velocity        = 10.0f;
    smooth          = 0.0f;
    windDir         = 45.0f;
    windAlign       = 1;
    reflectedWaves  = 0.3f;
    waveSpeed       = 1.0f;
    waveHeight      = 1.0f;
    waveChop        = 1.0f;
    repeatTime      = 10000.f;
    sprintf(outputFolder, "/tmp");
    sprintf(postfix, "");
    startFrame      = 1;
    endFrame        = 1;
    fps             = 24;
    accumulateFoam  = 0;

    logLevel = 0;
    help();
}

void input::help()
{
    sprintf(helpMsg,"\naaOcean Terminal/Console Application. Amaan Akram -- amaan@amaanakram.com\n");
    sprintf(helpMsg, "%saaOcean is a Tessendorf vector displacement map generator\n\n", helpMsg);
    sprintf(helpMsg, "%sexample usage:\n\n", helpMsg);
    sprintf(helpMsg, "%saaOcean --resolution=512 --seed=5 --oceanscale=100 --velocity=15.1 --winddir=45 --windalign=1 \
--reflectedwaves=1.0 --speed=1.1 --waveheight=3.0 --wavechop=2.0 --startframe=101 --endframe=150 --outputfolder=C:\\temp\n\n", helpMsg);
}

// define options and usage
enum  optionIndex { UNKNOWN, HELP, RESOLUTION, SPECTRUM, SEED, OCEANSCALE, \
    OCEANDEPTH, SURFACETENSION, VELOCITY, SMOOTH, \
    WINDDIRECTION, WINDALIGN, REFLECTEDWAVES, WAVESPEED, \
    WAVEHEIGHT, WAVECHOP, STARTFRAME, ENDFRAME, FPS, ACCUMULATEFOAM, OUTPUTFOLDER, POSTFIX, LOGLEVEL};

const option::Descriptor usage[] =
{
    {UNKNOWN, 0, "", "",                        option::Arg::None,     "\nUSAGE: aaOcean --help\n" },
    {RESOLUTION, 0,"res","resolution",          option::Arg::Optional, "  --resolution=<arg>  Defines map resolution in powers of two. \n  Examples: 64, 128, 256, 512, 1024, 2048, 4096\n  Default: 64\n" },
    {SPECTRUM, 0,"spec","spectrum",             option::Arg::Optional, "  --spectrum=<arg>  Select which spectrum to run an FFT on. For Philips, use 0\n  Default: 0\n" },
    {SEED, 0,"seed","seed",                     option::Arg::Optional, "  --seed=<arg>  Seed for random number generator. \n  Different seeds produce different oceans\n  Default: 1\n" },
    {OCEANSCALE, 0,"scale","oceanscale",        option::Arg::Optional, "  --oceanscale=<arg>\n  Defines ocean size in meters.\n  Default: 100.0\n" },
    {OCEANDEPTH, 0,"depth","oceandepth",        option::Arg::Optional, "  --oceandepth=<arg>\n  Defines ocean depth in meters.\n  Default: 10000.0\n" },
    {SURFACETENSION, 0,"tens","surfacetension", option::Arg::Optional, "  --surfacetension=<arg>\n  Layers small capillary waves on top of ocean\n  Default: 0.0\n " },
    {VELOCITY, 0,"vel","velocity",              option::Arg::Optional, "  --velocity=<arg>\n  Velocity in meters/sec. Controls size of waves on ocean surface\n  Active range 1.0 to 30.0\n  Default 10\n" },
    {SMOOTH, 0,"sm","smooth",                   option::Arg::Optional, "  --smooth=<arg>\n  Removes high-frequency waves\n  Default: 0.0\n" },
    {WINDDIRECTION, 0,"dir","winddir",          option::Arg::Optional, "  --winddir=<arg>\n  Wind direction in degrees\n  Default: 45\n" },
    {WINDALIGN, 0,"align","windalign",          option::Arg::Optional, "  --windalign=<arg>\n  Stretches waves perpendicular to wind direction\n  Simulates 'stretched' waves in shallow water\n  Default: 0\n" },
    {REFLECTEDWAVES, 0,"ref","reflectedwaves",  option::Arg::Optional, "  --reflectedwaves=<arg>\n  Control waves travelling opposite wind direction\n  1.0 removes them completely, 0.0 keeps them\n  Default 0.3\n" },
    {WAVESPEED, 0,"sp","speed",                 option::Arg::Optional, "  --speed=<arg>\n  Speed multiplier\n  Default: 1.0\n" },
    {WAVEHEIGHT, 0,"ht","waveheight",           option::Arg::Optional, "  --waveheight=<arg>\n  Height multiplier\n  Default 1.0\n" },
    {WAVECHOP, 0,"chp","wavechop",              option::Arg::Optional, "  --wavechop=<arg>\n  Choppiness, makes sharp wave peaks\n  Default: 1.0\n" },
    {STARTFRAME, 0,"start","startframe",        option::Arg::Optional, "  --startframe=<arg>\n  Start Frame for multi-frame output\n  Default: 1\n" },
    {ENDFRAME, 0,"end","endframe",              option::Arg::Optional, "  --endframe=<arg>\n  End Frame for multi-frame output\n  Default: 1\n" },
    {FPS, 0,"fps","fps",                        option::Arg::Optional, "  --fps=<arg>\n  Frames Per Second. Used to compute time increment\n  Default 24.0\n" },
    {ACCUMULATEFOAM, 0,"accFoam","accFoam",     option::Arg::Optional, "  --accFoam=<arg>\n  Accumulates Foam over time\n  Default 0\n" },
    {OUTPUTFOLDER, 0,"o","outputfolder",        option::Arg::Optional, "  --outputfolder=<arg>\n  Folder to write OpenEXR sequences to\n" },
    {POSTFIX, 0,"pfix","postfix",               option::Arg::Optional, "  --postfix=<arg>\n  Postfix for output filename\n" },
    {LOGLEVEL, 0,"log","loglevel",              option::Arg::Optional, "  --loglevel=<arg>\n Set to 1 for added info about what aaOcean is doing" },
    {HELP, 0,"h", "help",                       option::Arg::None, "" },
    {0,0,0,0,0,0}
};

bool processInput(int argc, char** argv, input &oceanInput)
{
    FILELog::ReportingLevel() = logINFO;
    int level;
    int res;
    char msg[512];

    argc-=(argc>0); argv+=(argc>0); 

    option::Stats  stats(usage, argc, argv);
    option::Option* options = new option::Option[stats.options_max];
    option::Option* buffer  = new option::Option[stats.buffer_max];
    option::Parser parse(usage, argc, argv, options, buffer);

    if (parse.error())
        return 0;

    if (argc == 0)
    {
        option::printUsage(std::cout, usage);
        return 0;
    }

    for (int i = 0; i < parse.optionsCount(); ++i)
    {
        option::Option& opt = buffer[i];
        switch (opt.index())
        {
        case HELP:
            printf(oceanInput.helpMsg);
            option::printUsage(std::cout, usage);
            break;
        case LOGLEVEL:
            level = atoi(opt.arg);
            if(level > 0)
            {
                sprintf(msg, "logging DEBUG and INFO messages with log level set to %d\n", atoi(opt.arg));
                FILELog::ReportingLevel() = logDEBUG;
                LOG(logDEBUG) << msg;
            }
            break;
        case RESOLUTION:
            res = atoi(opt.arg);
            if((res & (res-1))==0 && res > 0)
            {
                oceanInput.resolution = res;
                sprintf(msg, "--resolution with argument '%s'\n", opt.arg);
                LOG(logDEBUG) << msg;
                break;
            }
            else
            {
                sprintf(msg, "invalid resolution with argument '%s'. Please use power of 2, like 64, 128..\n", opt.arg);
                LOG(logDEBUG) << msg;
                exit(0);
            }
        case SPECTRUM:
            oceanInput.spectrum = atoi(opt.arg);
            sprintf(msg, "--spectrum with argument '%s'\n", opt.arg);
            LOG(logDEBUG) << msg;
            break;
        case SEED:
            oceanInput.seed = atoi(opt.arg);
            sprintf(msg, "--seed with argument '%s'\n", opt.arg);
            LOG(logDEBUG) << msg;
            break;
        case OCEANSCALE:
            oceanInput.oceanScale = (float)atof(opt.arg);
            sprintf(msg, "--oceanscale with argument '%s'\n", opt.arg);
            LOG(logDEBUG) << msg;
            break;
        case OCEANDEPTH:
            oceanInput.oceanDepth = (float)atof(opt.arg);
            sprintf(msg, "--oceandepth with argument '%s'\n", opt.arg);
            LOG(logDEBUG) << msg;
            break;
        case SURFACETENSION:
            oceanInput.surfaceTension = (float)atof(opt.arg);
            sprintf(msg, "--surfacetension with argument '%s'\n", opt.arg);
            LOG(logDEBUG) << msg;
            break;
        case VELOCITY:
            oceanInput.velocity = (float)atof(opt.arg);
            sprintf(msg, "--velocity with argument '%s'\n", opt.arg);
            LOG(logDEBUG) << msg;
            break;
        case SMOOTH:
            oceanInput.smooth = (float)atof(opt.arg);
            sprintf(msg, "--smooth with argument '%s'\n", opt.arg);
            LOG(logDEBUG) << msg;
            break;
        case WINDDIRECTION:
            oceanInput.windDir = (float)atof(opt.arg);
            sprintf(msg, "--winddir with argument '%s'\n", opt.arg);
            LOG(logDEBUG) << msg;
            break;
        case WINDALIGN:
            oceanInput.windAlign = atoi(opt.arg);
            sprintf(msg, "--windalign with argument '%s'\n", opt.arg);
            LOG(logDEBUG) << msg;
            break;
        case REFLECTEDWAVES:
            oceanInput.reflectedWaves = (float)atof(opt.arg);
            sprintf(msg, "--reflectedwaves with argument '%s'\n", opt.arg);
            LOG(logDEBUG) << msg;
            break;
        case WAVESPEED:
            oceanInput.waveSpeed = (float)atof(opt.arg);
            sprintf(msg, "--speed with argument '%s'\n", opt.arg);
            LOG(logDEBUG) << msg;
            break;
        case WAVEHEIGHT:
            oceanInput.waveHeight = (float)atof(opt.arg);
            sprintf(msg, "--waveheight with argument '%s'\n", opt.arg);
            LOG(logDEBUG) << msg;
            break;
        case WAVECHOP:
            oceanInput.waveChop = (float)atof(opt.arg);
            sprintf(msg, "--wavechop with argument '%s'\n", opt.arg);
            LOG(logDEBUG) << msg;
            break;
        case STARTFRAME:
            oceanInput.startFrame = atoi(opt.arg);
            sprintf(msg, "--startframe with argument '%s'\n", opt.arg);
            LOG(logDEBUG) << msg;
            break;
        case ENDFRAME:
            oceanInput.endFrame = atoi(opt.arg);
            sprintf(msg, "--endframe with argument '%s'\n", opt.arg);
            LOG(logDEBUG) << msg;
            break;
        case FPS:
            oceanInput.fps = (float)atof(opt.arg);
            sprintf(msg, "--fps with argument '%s'\n", opt.arg);
            LOG(logDEBUG) << msg;
            break;
        case ACCUMULATEFOAM:
            oceanInput.accumulateFoam = atoi(opt.arg);
            sprintf(msg, "--accuFoam with argument '%s'\n", opt.arg);
            LOG(logDEBUG) << msg;
            break;
        case OUTPUTFOLDER:
            sprintf(oceanInput.outputFolder, opt.arg);
            sprintf(msg, "--outputfolder with argument '%s'\n", opt.arg);
            LOG(logDEBUG) << msg;
            break;
        case POSTFIX:
            sprintf(oceanInput.postfix, opt.arg);
            sprintf(msg, "--postfix with argument '%s'\n", opt.arg);
            LOG(logDEBUG) << msg;
            break;
        }
    }

    for (option::Option* opt = options[UNKNOWN]; opt; opt = opt->next())
        LOG(logWARNING) << "Unknown option: " << std::string(opt->name,opt->namelen) << "\n";

    for (int i = 0; i < parse.nonOptionsCount(); ++i)
        LOG(logWARNING) << "Non-option #" << i << ": " << parse.nonOption(i) << "\n";

    delete[] options;
    delete[] buffer;

    return 1;
}

#endif // INPUT_H