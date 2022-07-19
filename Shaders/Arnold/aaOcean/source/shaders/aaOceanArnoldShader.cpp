// aaOcean
// Author: Amaan Akram 
// www.amaanakram.com
// Outputs RGBA, with Vector Displacement in RGB, and foam in Alpha
//
// LICENSE: 
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html
//
// A "New BSD" License for aaOcean can be obtained by contacting the author
// For more details on aaOcean and associated 3rd Party licenses, please see
// license.txt file that is part of the aaOcean repository:
// https://bitbucket.org/amaanakram/aaocean

#include <limits>
#include <string>
#include "ai.h"
#include "aaOceanClass.cpp"

#ifdef WRITE_OPENEXR
#include "openEXROutput.h"
#endif /* WRITE_EXR */

AI_SHADER_NODE_EXPORT_METHODS(aaOceanMethods);

enum aaOceanParams
{
    p_uv_coords,
    p_useUVInput,
    p_fade,
    p_resolution,
    p_oceanScale,
    p_oceanDepth,
    p_surfaceTension,
    p_seed,
    p_waveHeight,
    p_velocity,
    p_waveSpeed,
    p_chopAmount,
    p_cutoff,
    p_windDir,
    p_damp,
    p_windAlign,
    p_time,
    p_repeatTime,
    p_gamma,
    p_brightness,
    p_raw,
    p_invertFoam,
    p_fMin,
    p_fMax,
    p_writeFile,
    p_outputFolder,
    p_postfix,
    p_currentFrame,
    p_rotateUV,
    p_transform
};

node_parameters
{
    AtMatrix matrix44;
    AiM4Identity(matrix44);
    AiParameterVEC ( "uv_coords"        , 1.0f, 1.0f, 1.0f);
    AiParameterBOOL( "use_uv_input"     , 0);
    AiParameterFLT ( "fade"             , 0.0f);
    AiParameterINT ( "resolution"       , 4);
    AiParameterFLT ( "oceanScale"       , 100.0f);
    AiParameterFLT ( "oceanDepth"       , 10000.0f);
    AiParameterFLT ( "surfaceTension"   , 0.0f);
    AiParameterINT ( "seed"             , 1);
    AiParameterFLT ( "waveHeight"       , 5.0f);
    AiParameterFLT ( "velocity"         , 4.5f);
    AiParameterFLT ( "waveSpeed"        , 1.0f);
    AiParameterFLT ( "chopAmount"       , 1.0f);
    AiParameterFLT ( "cutoff"           , 0.0f);
    AiParameterFLT ( "windDir"          , 45.0f);
    AiParameterFLT ( "damp"             , 0.985f);
    AiParameterINT ( "windAlign"        , 0);
    AiParameterFLT ( "time"             , 0.1f);
    AiParameterFLT ( "repeatTime"       , 1000.f);
    AiParameterFLT ( "gamma"            , 1.0f);
    AiParameterFLT ( "brightness"       , 1.0f);
    AiParameterBOOL( "raw"              , 0);
    AiParameterBOOL( "invertFoam"       , 0);
    AiParameterFLT ( "fMin"             , 0.0f);
    AiParameterFLT ( "fMax"             , 0.0f);
    AiParameterBOOL( "writeFile"        , 0);
    AiParameterStr ( "outputFolder"     , "");
    AiParameterStr ( "postfix"          , "");
    AiParameterINT ( "currentFrame"     , 1);
    AiParameterBOOL( "rotateUV"         , 0);
    AiParameterMTX ( "transform"        , matrix44);
}

node_update
{
    // retrieve ocean pointer from user-data
    aaOcean *pOcean = (aaOcean *)AiNodeGetLocalData(node);

    // main input function
    pOcean->input(
        params[p_resolution].INT,
        params[p_seed].INT,
        params[p_oceanScale].FLT,
        params[p_oceanDepth].FLT,
        params[p_surfaceTension].FLT,
        params[p_velocity].FLT,
        params[p_cutoff].FLT,
        params[p_windDir].FLT,
        params[p_windAlign].INT,
        params[p_damp].FLT,
        params[p_waveSpeed].FLT,
        params[p_waveHeight].FLT,
        params[p_chopAmount].FLT,
        params[p_time].FLT,
        params[p_repeatTime].FLT,
        TRUE);

    // see if user has requested normalized or raw foam values
    bool rawOutput = params[p_raw].BOOL;
    if(pOcean->isChoppy() && !rawOutput)
    {
        float outMin, outMax;
        float   fmin        = params[p_fMin].FLT;
        float   fmax        = params[p_fMax].FLT;
        pOcean->getFoamBounds(outMin, outMax);

        float epsilon = 1e-3f;
        if(!isfEqual(fmin, outMin, epsilon) || !isfEqual(fmax, outMax, epsilon) )
            AiMsgWarning("[aaOcean Shader] Foam Min/Max mismatch. Please set the Foam Min/Max values in foam shader to Min: %f, Max: %f", 
                        outMin, outMax);
    }
}

shader_evaluate
{
    // retrieve ocean pointer from user-data
    aaOcean *pOcean = (aaOcean*)AiNodeGetLocalData(node);

    // get our UV's
    AtPoint uvPt;
    if(AiShaderEvalParamBool(p_useUVInput))
        uvPt = AiShaderEvalParamVec(p_uv_coords);
    else
    {
        uvPt.x = sg->u;
        uvPt.y = sg->v;
    }

    // retrieve heightfield in world space
    AtPoint worldSpaceDisplacementVec;
    worldSpaceDisplacementVec.y = pOcean->getOceanData(uvPt.x, uvPt.y, aaOcean::eHEIGHTFIELD);

    if(pOcean->isChoppy())
    {
        // retrieve chop displacement
        worldSpaceDisplacementVec.x = pOcean->getOceanData(uvPt.x, uvPt.y, aaOcean::eCHOPX);
        worldSpaceDisplacementVec.z = pOcean->getOceanData(uvPt.x, uvPt.y, aaOcean::eCHOPZ);

        // retrieve foam and store it in our alpha channel
        sg->out.RGBA.a = pOcean->getOceanData(uvPt.x, uvPt.y, aaOcean::eFOAM);

        // see if user has requested normalized or raw foam values
        if(!AiShaderEvalParamBool(p_raw))
        {
            // get normalization weights
            float   fmin = AiShaderEvalParamFlt(p_fMin);
            float   fmax = AiShaderEvalParamFlt(p_fMax);
            
            // fitting to 0 - 1 range using rescale(...)
            // invert result to put foam on wave peaks
            if(AiShaderEvalParamBool(p_invertFoam))
                sg->out.RGBA.a  = 1.0f - rescale(sg->out.RGBA.a, fmin, fmax, 0.0f, 1.0f);
            else
                sg->out.RGBA.a  = rescale(sg->out.RGBA.a, fmin, fmax, 0.0f, 1.0f);

            // removing negative leftovers
            sg->out.RGBA.a  = maximum<float>(sg->out.RGBA.a, 0.0f);
            // apply gamma
            sg->out.RGBA.a  = pow(sg->out.RGBA.a, AiShaderEvalParamFlt(p_gamma));
            sg->out.RGBA.a *= AiShaderEvalParamFlt(p_brightness);
        }
    }
    else
    {
        // return only heightfield, since ocean isn't choppy, 
        // and set foam to 0.0 since it is derived from choppy oceans
        worldSpaceDisplacementVec.x = worldSpaceDisplacementVec.z = 0.0f;
        sg->out.RGBA.a = 0.f;
    }

    // get space-transform matrix
    AtMatrix matrix;
    AiM4Identity(matrix);
    AiM4Copy(matrix, *AiShaderEvalParamMtx(p_transform));

    // convert to the local space of the user-specified transform matrix
    AtPoint oceanLocalSpace;
    oceanLocalSpace.x = oceanLocalSpace.z = 0.f;
    AiM4PointByMatrixMult(&oceanLocalSpace, matrix, &worldSpaceDisplacementVec);

    // store result in output
    sg->out.RGBA.r = oceanLocalSpace.x;
    sg->out.RGBA.g = oceanLocalSpace.y;
    sg->out.RGBA.b = oceanLocalSpace.z;
}

node_initialize
{
    // store a new ocean pointer in user-data
    aaOcean *pOcean;
    pOcean = new aaOcean;
    AiNodeSetLocalData(node, pOcean);
    AiMsgInfo("[aaOcean Arnold] Created new aaOcean data");
}

node_finish
{
    // retrieve ocean pointer from user-data
    aaOcean *pOcean = (aaOcean *)AiNodeGetLocalData(node);
    
    #ifdef WRITE_OPENEXR
    if(AiNodeGetBool(node, "writeFile"))
    {
        const char *outputFolder = AiNodeGetStr(node, "outputFolder");
        if(!dirExists(outputFolder))
            AiMsgWarning("[aaOcean] Invalid folder path: %s", outputFolder);
        else
        {
            char outputFileName[255];
            sprintf(outputFileName, "none");
            const char* postfix = AiNodeGetStr(node, "postfix");
            int frame = AiNodeGetInt(node, "currentFrame");
            oceanDataToEXR(pOcean,&outputFolder[0], &postfix[0], frame, &outputFileName[0]);
            AiMsgInfo("[aaOcean Arnold] Image written to %s", outputFileName);
        }
    }
    #endif
    // cleanup ocean
    delete pOcean;
    AiMsgInfo("[aaOcean Arnold] Deleted aaOcean data");
}

node_loader
{
   if (i > 0)
      return FALSE;

   node->methods      = aaOceanMethods;
   node->output_type  = AI_TYPE_RGBA;
   node->name         = "aaOcean";
   node->node_type    = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return TRUE;
}

#ifdef WRITE_OPENEXR

#endif /* WRITE_EXR */
