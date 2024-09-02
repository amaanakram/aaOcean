// aaOcean
// Author: Amaan Akram 
// https://linkedin.com/in/amaan
//
// LICENSE: 
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html
//
// A "New BSD" License for aaOcean can be obtained by contacting the author
// For more details on aaOcean and associated 3rd Party licenses, please see
// license.txt file that is part of the aaOcean repository:
// https://github.com/amaanakram/aaOcean

// usage:
// ./aaOcean --resolution=256 --seed=1 --oceanscale=100 --oceandepth=10000 
// --surfacetension=0.0 --velocity=15 --smooth=0 --winddir=45 --windalign=1 
// --reflectedwaves=0.3 --speed=1.0 --waveheight=1.0 --wavechop=1.0 -
// -startframe=101 --endframe=101 --outputfolder=/tmp --postfix=Your_Postfix

// not all of the above arguments need to be specified. See default values in 
// input.h in input constructor

#include <stdio.h>
#include <iostream>
#include <iomanip> 
#include <omp.h>
#include <string.h>
#include <math.h> 

#include "functionLib.h"
#include "input.h"
#include "aaOceanClass.cpp"
#include "openEXROutput.h"

void writeFoam(float *&foam, input &oceanInput, int currentFrame)
{
    // TODO: move this function into openExrOutput.h
    int dimension = oceanInput.resolution;

    EXR::Array2D<float> pixels(dimension, dimension);
    for (int i = 0; i < dimension; i++)
    {
        #pragma omp parallel for
        for (int j = 0; j < dimension; j++)
            pixels[j][i] = foam[i*dimension + j];
    }

    char outputFileName[512];
    char outputFolder[512];
    char postfix[512];

    strcpy(outputFolder, oceanInput.outputFolder);
    strcpy(postfix, oceanInput.postfix);
    strcpy(postfix, "_accumulatedFoam");

    genFullFilePath(&outputFileName[0], &outputFolder[0], &postfix[0], currentFrame);

    writeSingleChannelExr(&outputFileName[0], &pixels[0][0], dimension, dimension);
    LOG(logINFO) << "written accumulated foam";
}

void accumulateFoam(aaOcean *&pOcean, float *&foamCurrentFrame, float *&accumulatedFoam, input &oceanInput, int currentFrame)
{
    int size = oceanInput.resolution * oceanInput.resolution;
    pOcean->getOceanArray(foamCurrentFrame, aaOcean::eFOAM); // this overwrites foam

    if (currentFrame > oceanInput.startFrame)
    {
        #pragma omp parallel for 
        for (int i = 0; i < size; ++i)
            accumulatedFoam[i] = clamp(foamCurrentFrame[i], 0.f, 100.f) + accumulatedFoam[i] * 0.9f;
    }
    else
    {
        #pragma omp parallel for 
        for (int i = 0; i < size; ++i)
            accumulatedFoam[i] = clamp(foamCurrentFrame[i], 0.f, 100.f);
    }
}

int main(int argc, char* argv[])
{
    char msg[512];
    input oceanInput;
    LOG(logINFO) << oceanInput.helpMsg;

    if(!processInput(argc, argv, oceanInput))
        return 1;

    int dimension = oceanInput.resolution;
    LOG(logINFO) << "Starting ocean evaluation for " << dimension << "x" << dimension;

    float *foamCurrentFrame = 0;
    float *accumulatedFoam = 0;

    if (oceanInput.accumulateFoam > 0)
    {
        int size = dimension * dimension * sizeof(float);
        foamCurrentFrame = (float*)malloc(size);
        accumulatedFoam = (float*)malloc(size);
    }

    // modify dimension to account of aaOcean's
    // arbitrary resolution scaling of 4
    dimension = (int)(log((float)dimension)/ log(2.0f)) - 4;
    aaOcean *pOcean = new aaOcean;
    
    float timestep = 1.0f/oceanInput.fps;
    int currentFrame = oceanInput.startFrame;

    while(currentFrame <= oceanInput.endFrame)
    {
        int absoluteFrame = currentFrame - oceanInput.startFrame;
        float time = float(absoluteFrame) * timestep;

        sprintf(msg,"Evaluating for time %0.3f seconds and frame %d", time, currentFrame);
        LOG(logINFO) << msg;

        pOcean->input(
            dimension,                  // resolution 
            oceanInput.spectrum,        // spectrum
            oceanInput.seed,            // seed
            oceanInput.oceanScale,      // ocean scale
            oceanInput.oceanDepth,      // ocean depth
            oceanInput.surfaceTension,  // surface tension
            oceanInput.velocity,        // velocity
            oceanInput.smooth,          // cutoff/smooth
            oceanInput.windDir,         // wind dir
            oceanInput.windAlign,       // wind align
            oceanInput.reflectedWaves,  // damp
            oceanInput.waveSpeed,       // wave speed
            oceanInput.waveHeight,      // wave height
            oceanInput.waveChop,        // chop amount
            time,                       // time in seconds
            oceanInput.repeatTime,      // repeat/loop time
            1);                         // calculate foam
        
        LOG(logDEBUG) << "Logging Ocean Core messages\n" << pOcean->m_state;
        LOG(logINFO) << msg;
        
        char outputFileName[512];
        oceanDataToEXR(pOcean, 
                       &oceanInput.outputFolder[0], 
                       &oceanInput.postfix[0], 
                       currentFrame, 
                       &outputFileName[0]);

        if (oceanInput.accumulateFoam > 0)
        {
            accumulateFoam(pOcean, foamCurrentFrame, accumulatedFoam, oceanInput, currentFrame);
            writeFoam(accumulatedFoam, oceanInput, currentFrame);
        }

        sprintf(msg,"Written OpenEXR object-space vector-displacement map");
        LOG(logDEBUG) << msg;
        sprintf(msg,"OpenEXR image location: %s", &outputFileName[0]);
        LOG(logINFO) << msg;
        sprintf(msg,"OpenEXR RGB contains position, Alpha contains raw foam/spray emission data");
        LOG(logDEBUG) << msg;

        currentFrame++;
    }

    if (foamCurrentFrame)
    {
        free(foamCurrentFrame);
        free(accumulatedFoam);
    }
    delete pOcean;
    return 0;
}
