// aaOcean
// Author: Amaan Akram 
// https://linkedin.com/in/amaan
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
// https://github.com/amaanakram/aaOcean

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
    p_transform,
    p_spectrum,
    p_randWeight,
    p_spectrumMult,
    p_peakSharpening,
    p_jswpfetch,
    p_swell
};

node_parameters
{
    AtMatrix matrix44 = AiM4Identity();
    AiParameterVec ( "uv_coords"        , 1.0f, 1.0f, 1.0f);
    AiParameterBool( "use_uv_input"     , 0);
    AiParameterFlt ( "fade"             , 0.0f);
    AiParameterInt ( "resolution"       , 4);
    AiParameterFlt ( "oceanScale"       , 100.0f);
    AiParameterFlt ( "oceanDepth"       , 10000.0f);
    AiParameterFlt ( "surfaceTension"   , 0.0f);
    AiParameterInt ( "seed"             , 1);
    AiParameterFlt ( "waveHeight"       , 5.0f);
    AiParameterFlt ( "velocity"         , 4.5f);
    AiParameterFlt ( "waveSpeed"        , 1.0f);
    AiParameterFlt ( "chopAmount"       , 1.0f);
    AiParameterFlt ( "cutoff"           , 0.0f);
    AiParameterFlt ( "windDir"          , 45.0f);
    AiParameterFlt ( "damp"             , 0.985f);
    AiParameterInt ( "windAlign"        , 0);
    AiParameterFlt ( "time"             , 0.1f);
    AiParameterFlt ( "repeatTime"       , 1000.f);
    AiParameterFlt ( "gamma"            , 1.0f);
    AiParameterFlt ( "brightness"       , 1.0f);
    AiParameterBool( "raw"              , 0);
    AiParameterBool( "invertFoam"       , 0);
    AiParameterFlt ( "fMin"             , 0.0f);
    AiParameterFlt ( "fMax"             , 0.0f);
    AiParameterBool( "writeFile"        , 0);
    AiParameterStr ( "outputFolder"     , "");
    AiParameterStr ( "postfix"          , "");
    AiParameterInt ( "currentFrame"     , 1);
    AiParameterBool( "rotateUV"         , 0);
    AiParameterMtx ( "transform"        , matrix44);
    AiParameterInt ( "spectrum"         , 0);
    AiParameterFlt ( "randWeight"       , 0.0f);
    AiParameterFlt ( "spectrumMult"     , 1.0f);
    AiParameterFlt ( "peakSharpening"   , 1.0f);
    AiParameterFlt ( "jswpfetch"        , 100.0f);
    AiParameterFlt ( "swell"            , 0.0f);
}

node_update
{
    Timer timer;
    // retrieve ocean pointer from user-data
    aaOcean *pOcean = reinterpret_cast<aaOcean*>(AiNodeGetLocalData(node));

    // main input function
    pOcean->input(
        AiNodeGetInt(node, "resolution"),
        AiNodeGetInt(node, "spectrum"), 
        AiNodeGetInt(node, "seed"),
        AiNodeGetFlt(node, "oceanScale"),
        AiNodeGetFlt(node, "oceanDepth"),
        AiNodeGetFlt(node, "surfaceTension"),
        AiNodeGetFlt(node, "velocity"),
        AiNodeGetFlt(node, "cutoff"),
        AiNodeGetFlt(node, "windDir"),
        AiNodeGetInt(node, "windAlign"),
        AiNodeGetFlt(node, "damp"),
        AiNodeGetFlt(node, "waveSpeed"),
        AiNodeGetFlt(node, "waveHeight"),
        AiNodeGetFlt(node, "chopAmount"),
        AiNodeGetFlt(node, "time"),
        AiNodeGetFlt(node, "repeatTime"),
        TRUE,
        AiNodeGetFlt(node, "randWeight"),
        AiNodeGetFlt(node, "spectrumMult"),
        AiNodeGetFlt(node, "peakSharpening"),
        AiNodeGetFlt(node, "jswpfetch"),
        AiNodeGetFlt(node, "swell"));

    // see if user has requested normalized or raw foam values
    bool rawOutput = AiNodeGetBool(node, "raw");
    if(pOcean->isChoppy() && !rawOutput)
    {
        float outMin, outMax;
        float   fmin = AiNodeGetFlt(node, "fMin");
        float   fmax = AiNodeGetFlt(node, "fMax");
        pOcean->getFoamBounds(outMin, outMax);

        float epsilon = 1e-3f;
        if(!isfEqual(fmin, outMin, epsilon) || !isfEqual(fmax, outMax, epsilon) )
            AiMsgWarning("[aaOcean Shader] Foam Min/Max mismatch. Please set the Foam Min/Max values in foam shader to Min: %f, Max: %f", 
                        outMin, outMax);
    }
    timer.printElapsed("[aaOcean Shader] node_update finished");
}

shader_evaluate
{
    // retrieve ocean pointer from user-data
    aaOcean* pOcean = reinterpret_cast<aaOcean*>(AiNodeGetLocalData(node));

    // get our UV's
    AtVector2 uvPt;
    if(AiShaderEvalParamBool(p_useUVInput))
        uvPt = AiShaderEvalParamVec(p_uv_coords);
    else
    {
        uvPt.x = sg->u;
        uvPt.y = sg->v;
    }

    // retrieve heightfield in world space
    AtVector worldSpaceDisplacementVec;
    worldSpaceDisplacementVec.y = pOcean->getOceanData(uvPt.x, uvPt.y, aaOcean::eHEIGHTFIELD);

    if(pOcean->isChoppy())
    {
        // retrieve chop displacement
        worldSpaceDisplacementVec.x = pOcean->getOceanData(uvPt.x, uvPt.y, aaOcean::eCHOPX);
        worldSpaceDisplacementVec.z = pOcean->getOceanData(uvPt.x, uvPt.y, aaOcean::eCHOPZ);

        // retrieve foam and store it in our alpha channel
        sg->out.RGBA().a = pOcean->getOceanData(uvPt.x, uvPt.y, aaOcean::eFOAM);

        // see if user has requested normalized or raw foam values
        if(!AiShaderEvalParamBool(p_raw))
        {
            // get normalization weights
            float   fmin = AiShaderEvalParamFlt(p_fMin);
            float   fmax = AiShaderEvalParamFlt(p_fMax);
            
            // fitting to 0 - 1 range using rescale(...)
            // invert result to put foam on wave peaks
            if(AiShaderEvalParamBool(p_invertFoam))
                sg->out.RGBA().a = 1.0f - rescale(sg->out.RGBA().a, fmin, fmax, 0.0f, 1.0f);
            else
                sg->out.RGBA().a = rescale(sg->out.RGBA().a, fmin, fmax, 0.0f, 1.0f);

            // removing negative leftovers
            sg->out.RGBA().a = maximum<float>(sg->out.RGBA().a, 0.0f);
            // apply gamma
            sg->out.RGBA().a = pow(sg->out.RGBA().a, AiShaderEvalParamFlt(p_gamma));
            sg->out.RGBA().a *= AiShaderEvalParamFlt(p_brightness);
        }
    }
    else
    {
        // return only heightfield, since ocean isn't choppy, 
        // and set foam to 0.0 since it is derived from choppy oceans
        worldSpaceDisplacementVec.x = worldSpaceDisplacementVec.z = 0.0f;
        sg->out.RGBA().a = 0.f;
    }

    // convert to the local space of the user-specified transform matrix
    AtMatrix transform_matrix = *AiShaderEvalParamMtx(p_transform);
    AtVector oceanLocalSpace = AiM4PointByMatrixMult(transform_matrix, worldSpaceDisplacementVec);

    // store result in output
    sg->out.RGBA().r = oceanLocalSpace.x;
    sg->out.RGBA().g = oceanLocalSpace.y;
    sg->out.RGBA().b = oceanLocalSpace.z;
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
     aaOcean* pOcean = reinterpret_cast<aaOcean*>(AiNodeGetLocalData(node));
    
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
