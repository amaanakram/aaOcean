// aaOcean Noise Shader
// Author: Amaan Akram 
// https://linkedin.com/in/amaan
// Outputs scalar noise for use in tiling
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

#define HAVE_SSE2 1
#include <cmath>
#include "ai.h"
#include "perlin_noise.h"

AI_SHADER_NODE_EXPORT_METHODS(aaOceanNoiseMethods);

struct ShaderData {
    int     seed = 1;
    int     octaves = 3;
    float   frequency = 1;
    float   time = 1;
    float   u_mult = 1;
    float   v_mult = 1;
    float   gamma = 1;
    float   amplitude = 1;

    PerlinNoise *pNoise = nullptr;
};

enum aaOceanNoiseParams
{
    p_seed,
    p_octaves,
    p_frequency,
    p_time,
    p_uv,
    p_u_mult,
    p_v_mult,
    p_gamma,
    p_amplitude
};

node_parameters
{
    AiParameterInt ( "seed"             , 1);
    AiParameterInt ( "octaves"          , 3);
    AiParameterFlt ( "frequency"        , 10.f);
    AiParameterFlt ( "time"             , 1.f);
    AiParameterVec ( "uv"               , 0.0f, 0.0f, 0.5f);
    AiParameterFlt ( "u_mult"           , 1.0f);
    AiParameterFlt ( "v_mult"           , 1.0f);
    AiParameterFlt ( "gamma"            , 1.0f);
    AiParameterFlt ( "amplitude"        , 1.0f);
}

node_update
{
    ShaderData *data    = (ShaderData*)AiNodeGetLocalData(node);
    data->seed          = AiNodeGetInt(node, "seed");
    data->octaves       = AiNodeGetInt(node, "octaves");
    data->frequency     = std::max(AiNodeGetFlt(node, "frequency"), FLT_MIN);
    data->time          = AiNodeGetFlt(node, "time");
    data->u_mult        = AiNodeGetFlt(node, "u_mult");
    data->v_mult        = AiNodeGetFlt(node, "v_mult");
    data->gamma         = std::max(AiNodeGetFlt(node, "gamma"), FLT_MIN);
    data->amplitude     = std::max(AiNodeGetFlt(node, "amplitude"), FLT_MIN);

    if(data->pNoise)
        delete data->pNoise;
    data->pNoise = new PerlinNoise(data->seed);
}

shader_evaluate
{
    ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);
    PerlinNoise *pNoise = data->pNoise;

    AtVector uv = AiShaderEvalParamVec(p_uv);
    float u_mult = AiShaderEvalParamFlt(p_u_mult);
    float v_mult = AiShaderEvalParamFlt(p_v_mult);
    float time = AiShaderEvalParamFlt(p_time);
    float frequency = AiShaderEvalParamFlt(p_frequency);
    float amplitude = AiShaderEvalParamFlt(p_amplitude);
    float gamma = AiShaderEvalParamFlt(p_gamma);

    uv.x = uv.x * u_mult;
    uv.y = uv.y * v_mult;

    float noise = pNoise->octaveNoise(uv.x, uv.y, time, data->octaves, frequency);
    noise = (noise + 1.f) * 0.5f;
    noise = std::pow(noise, gamma) * amplitude;

    sg->out.RGB() = AtRGB(noise, noise, noise);
}

node_initialize
{
    ShaderData *data = new ShaderData();
    AiNodeSetLocalData(node, data);
}

node_finish
{
    ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);

    if (data->pNoise)
        delete data->pNoise;
    if(data)
        delete data;
}

node_loader
{
   if (i > 0)
      return false;

   node->methods      = aaOceanNoiseMethods;
   node->output_type  = AI_TYPE_RGB;
   node->name         = "aaOceanNoise";
   node->node_type    = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return true;
}
