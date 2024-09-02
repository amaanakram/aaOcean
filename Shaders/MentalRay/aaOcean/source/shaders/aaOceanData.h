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

#ifndef aaOceanDataShader_H
#define aaOceanDataShader_H 

#include <shader.h>

typedef struct
{
    miVector    uv_coords;          // uv_coords
    miBoolean   use_uv_input;
    miScalar    fade;
    miInteger   resolution;
    miScalar    oceanScale;
    miScalar    oceanDepth;
    miScalar    surfaceTension;
    miInteger   seed;
    miScalar    waveHeight;
    miScalar    velocity;
    miScalar    waveSpeed;
    miScalar    chopAmount;
    miScalar    cutoff;
    miScalar    windDir;
    miScalar    damp;
    miInteger   windAlign;
    miScalar    time;
    miScalar    repeatTime;

    miScalar    gamma;
    miScalar    brightness;
    miBoolean   rawOutput;
    miBoolean   invertFoam;
    miScalar    fmin;
    miScalar    fmax;

    miScalar    layerOcean;
    miMatrix    transform;

    miBoolean   writeFile;
    miTag       outputFolder;
    miTag       postfix;
    miInteger   currentFrame;

} aaOceanDataShader_t;

#endif // aaOceanDataShader_H


