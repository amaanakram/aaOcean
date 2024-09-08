// aaOcean Spectrum Generator
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

#include <string>
#include "ai.h"
#include "spectrum.h"
#include "aaOceanClass.cpp"

struct MsgContainer
{
    AtString msgName;
    aaSpectrum *pSpectrum;
};

AI_SHADER_NODE_EXPORT_METHODS(OceanSpectrumMethods);

enum aaOceanParams
{
    p_spectrumName,
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
    p_currentFrame,
    p_spectrum,
    p_randWeight,
    p_spectrumMult,
    p_peakSharpening,
    p_jswpfetch,
    p_swell
};

node_parameters
{
    AiParameterStr ( "spectrumName"     , "OceanSpectrum");
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
    AiParameterInt ( "currentFrame"     , 1);
    AiParameterInt ( "spectrum"         , 0);
    AiParameterFlt ( "randWeight"       , 0.0f);
    AiParameterFlt ( "spectrumMult"     , 1.0f);
    AiParameterFlt ( "peakSharpening"   , 1.0f);
    AiParameterFlt ( "jswpfetch"        , 1.0f);
    AiParameterFlt ( "swell"            , 0.0f);
}

node_update
{
    Timer timer;
    // retrieve ocean pointer from user-data
    MsgContainer *msg = reinterpret_cast<MsgContainer*>(AiNodeGetLocalData(node));
    msg->msgName = AiNodeGetStr(node, "spectrumName");
    
    aaSpectrum *pSpectrum = msg->pSpectrum;

    pSpectrum->input(
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
    
    pSpectrum->BuildSpectrum();
    timer.printElapsed("[aaOcean Shader] node_update finished");
}

shader_evaluate
{
    // retrieve ocean pointer from user-data
    MsgContainer *msg = reinterpret_cast<MsgContainer*>(AiNodeGetLocalData(node));
    AiStateSetMsgPtr(msg->msgName, (void*)msg);

    // get our UV's
    AtVector2 uvPt;
    uvPt.x = sg->u;
    uvPt.y = sg->v;

    AtRGB spectrum; 
    std::tie(spectrum.r, spectrum.g) = msg->pSpectrum->GetSpectralData(uvPt.x, uvPt.y);

    // store result in output
    sg->out.RGB().r = spectrum.r;
    sg->out.RGB().g = spectrum.g;
    sg->out.RGB().b = 0.f;
}

node_initialize
{
    // store a new ocean pointer in user-data

    MsgContainer *msg = new MsgContainer;
    msg->pSpectrum = new aaSpectrum;
    AiNodeSetLocalData(node, msg);
    AiMsgInfo("[aaOcean Arnold] Created new Ocean Spectrum data");
}

node_finish
{
    // retrieve ocean pointer from user-data
    MsgContainer* msg = reinterpret_cast<MsgContainer*>(AiNodeGetLocalData(node));
    delete msg->pSpectrum;
    delete msg;
    AiMsgInfo("[Spectrum] Deleted Ocean Spectrum");
}

node_loader
{
   if (i > 0)
      return FALSE;

   node->methods      = OceanSpectrumMethods;
   node->output_type  = AI_TYPE_RGB;
   node->name         = "OceanSpectrum";
   node->node_type    = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return TRUE;
}
