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

#ifndef aaOceanDataShader_CPP
#define aaOceanDataShader_CPP

#include <string>
#include <shader.h>
#include "aaOceanClass.cpp"
#include "oceanStore.h"
#include "aaOceanData.h"

#ifdef WRITE_OPENEXR
#include "openEXROutput.h"
#endif

extern "C" DLLEXPORT 
miBoolean aaOceanDataShader(miColor *result, miState *state, aaOceanDataShader_t *params)
{
    // retrieve ocean pointer from user-data
    oceanStore** os;
    mi_query( miQ_FUNC_USERPTR, state, 0, (void *)&os);
    aaOcean *pOcean = (*os)->ocean;

    miScalar    layerOcean  = *mi_eval_scalar(&params->layerOcean);
    miVector    *coord      =  mi_eval_vector(&params->uv_coords);
    miScalar    fade        = *mi_eval_scalar(&params->fade);
    
    if(*mi_eval_boolean(&params->use_uv_input) == FALSE)
    {
        coord->x = state->tex_list[0].x;
        coord->y = state->tex_list[0].y;
    }

    // retrieve heightfield
    miVector oceanWorldSpace;
    oceanWorldSpace.y = pOcean->getOceanData(coord->x, coord->y, aaOcean::eHEIGHTFIELD) * (1.0f - fade);
    oceanWorldSpace.x = oceanWorldSpace.z = 0.0f;
    
    if(pOcean->isChoppy())
    {
        // retrieve chop displacement
        oceanWorldSpace.x = pOcean->getOceanData(coord->x, coord->y, aaOcean::eCHOPX)  *  (1.0f - fade);
        oceanWorldSpace.z = pOcean->getOceanData(coord->x, coord->y, aaOcean::eCHOPZ)  *  (1.0f - fade);

        // retrieve foam and store it in our alpha channel
        result->a = pOcean->getOceanData(coord->x, coord->y, aaOcean::eFOAM);

        // see if user has requested normalized or raw foam values
        if(!*mi_eval_boolean(&params->rawOutput))
        {
            miScalar    gamma       = *mi_eval_scalar(&params->gamma);
            miScalar    brightness  = *mi_eval_scalar(&params->brightness);
            miScalar    fmin        = *mi_eval_scalar(&params->fmin);
            miScalar    fmax        = *mi_eval_scalar(&params->fmax);

            // fitting to 0 - 1 range using rescale(...)
            // invert result to put foam on wave peaks
            if(*mi_eval_boolean(&params->invertFoam))
                result->a  = 1.0f - rescale(result->a, fmin, fmax, 0.0f, 1.0f);
            else
                result->a  = rescale(result->a, fmin, fmax, 0.0f, 1.0f);

            result->a  = maximum<float>(result->a, 0.0f);
            result->a  = pow(result->a, gamma);
            result->a *= brightness;
            result->a *= (1.0f - fade);
        }
    }
    else
        result->a = 0.0f;

    // convert to the local space of the user-specified transform matrix
    miVector oceanLocalSpace;
    mi_vector_transform(&oceanLocalSpace, &oceanWorldSpace, (*os)->transform);

    // return the result
    result->r = oceanLocalSpace.x;
    result->g = oceanLocalSpace.y;
    result->b = oceanLocalSpace.z;

    return miTRUE;
}

extern "C" DLLEXPORT 
void aaOceanDataShader_init(miState *state, aaOceanDataShader_t *params, miBoolean *inst_init_req)
{
    if( !params )
    {
        // request per-shader-instance initialization
        *inst_init_req = miTRUE;
    }
    else
    {
        // evaluated any previously connected ocean shaders
        miScalar layerOcean = *mi_eval_scalar(&params->layerOcean);

        // allocate memory for our ocean object and store its pointer 
        // in shader's user-data construct
        oceanStore** os = NULL;
        mi_query( miQ_FUNC_USERPTR, state, 0, (void *)&os);
        *os = (oceanStore *)mi_mem_allocate( sizeof(oceanStore));

        // create a new instance of our ocean object
        // and create a new ocean pointer for convenient syntax
        (*os)->ocean = new aaOcean;
        aaOcean *pOcean = (*os)->ocean;

        if(!mi_matrix_isnull(mi_eval_transform(&params->transform)))
            mi_matrix_copy((*os)->transform, mi_eval_transform(&params->transform));
        else
        {
            mi_matrix_ident((*os)->transform);
            mi_warning("[aaOcean Shader] Shader's Transform input is zero. Using Identity instead");    
        }
    
        // retrieve user input for shader
        pOcean->input(
        *mi_eval_integer(&params->resolution), 
        *mi_eval_integer(&params->seed),
        *mi_eval_scalar (&params->oceanScale), 
        *mi_eval_scalar (&params->oceanDepth), 
        *mi_eval_scalar (&params->surfaceTension), 
        *mi_eval_scalar (&params->velocity), 
        *mi_eval_scalar (&params->cutoff), 
        *mi_eval_scalar (&params->windDir), 
        *mi_eval_integer(&params->windAlign), 
        *mi_eval_scalar (&params->damp), 
        *mi_eval_scalar (&params->waveSpeed), 
        *mi_eval_scalar (&params->waveHeight),
        *mi_eval_scalar (&params->chopAmount), 
        *mi_eval_scalar (&params->time),
        *mi_eval_scalar (&params->repeatTime),
        miTRUE);

        // get the tag of the currently running shader so that we can call its name
        miTag shaderInst; 
        mi_query(miQ_FUNC_TAG, state, 0, &shaderInst);
        mi_info("%s. Shader ID: %d", pOcean->getState(), shaderInst);   

        miBoolean rawOutput = *mi_eval_boolean(&params->rawOutput);
        if(pOcean->isChoppy() && !rawOutput)
        {
            float outMin, outMax;
            miScalar fmin = *mi_eval_scalar(&params->fmin);
            miScalar fmax = *mi_eval_scalar(&params->fmax);
            pOcean->getFoamBounds(outMin, outMax);

            float epsilon = 1e-3f;
            if( ( !isfEqual(fmin, outMin, epsilon) || !isfEqual(fmax, outMax, epsilon) ))
                mi_warning("[aaOcean Shader] Foam Min/Max mismatch. Please set the Foam Min/Max values in foam shader to Min: %f, Max: %f", 
                            outMin, outMax);
        }

        // clear arrays that are not required during shader evaluation
        pOcean->clearResidualArrays();

        mi_info("[aaOcean Shader] Data shader initiated. Shader ID: %d, location: %p", shaderInst, pOcean); 
    }
}

char* mitag_to_string(miTag tag, char *default_value) 
{ 
    char *result = default_value; 
    if (tag != 0) 
    { 
        result = (char*)mi_db_access(tag); 
        mi_db_unpin(tag); 
    } 
    return result; 
}

extern "C" DLLEXPORT 
void aaOceanDataShader_exit(miState *state, aaOceanDataShader_t *params)
{
    if( params )
    {
        oceanStore** os = NULL;
        if(mi_query( miQ_FUNC_USERPTR, state, 0, (void *)&os))
        {
            aaOcean *pOcean = (*os)->ocean;

            #ifdef WRITE_OPENEXR
            if(params->writeFile)
            {
                char* outputFolder  = mitag_to_string(params->outputFolder,"");
                if(!dirExists(outputFolder))
                    mi_error("[aaOcean] Invalid folder path: %s", outputFolder);
                else
                {
                    char outputFileName[255];
                    sprintf(outputFileName, "none");
                    char* postfix = mitag_to_string(params->postfix,"");
                    int frame = *mi_eval_integer(&params->currentFrame);
                    oceanDataToEXR(pOcean,&outputFolder[0], &postfix[0], frame, &outputFileName[0]);
                    mi_info("[aaOcean] Image written to %s", outputFileName);
                }
            }
            #endif

            if(pOcean)
                delete pOcean;
            mi_mem_release( (void *)(*os) );

            miTag shaderInst; // tag of the currently running shader
            mi_query(miQ_FUNC_TAG, state, 0, &shaderInst);
            mi_info("[aaOcean Shader] Data Shader ID %d, terminated at location: %p", shaderInst, pOcean);
        }
    }
}

extern "C" DLLEXPORT int aaOceanDataShader_version( )
{
    return( 1 );
}

#endif /* aaOceanDataShader_CPP */