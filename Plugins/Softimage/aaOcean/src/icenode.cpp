// aaOcean
// Author: Amaan Akram 
// www.amaanakram.com
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

#include "xsi_includes.h"
#include "aaOceanClass.cpp"
#include "icenode_portDefs.h"
#include "icenode_registration.h"

void multiplyMatrix(float *InVector, float *OutVector, XSI::CDataArrayMatrix4f &transform, int transformArraySize)
{
    OutVector[0] =  InVector[0] * transform[transformArraySize].GetValue(0,0) + 
                    InVector[1] * transform[transformArraySize].GetValue(1,0) + 
                    InVector[2] * transform[transformArraySize].GetValue(2,0);

    OutVector[1] =  InVector[0] * transform[transformArraySize].GetValue(0,1) + 
                    InVector[1] * transform[transformArraySize].GetValue(1,1) + 
                    InVector[2] * transform[transformArraySize].GetValue(2,1);

    OutVector[2] =  InVector[0] * transform[transformArraySize].GetValue(0,2) + 
                    InVector[1] * transform[transformArraySize].GetValue(1,2) + 
                    InVector[2] * transform[transformArraySize].GetValue(2,2);
}

SICALLBACK aaOcean_Init( CRef& in_ctxt )
{
    // create new aaOceanClass instance
    aaOcean *pOcean = new aaOcean;

    // store ocean pointer in user-data
    Context ctxt(in_ctxt);
    ctxt.PutUserData( (CValue::siPtrType)pOcean);

    Application().LogMessage("[aaOcean ICE] : Successfully allocated memory");
    return CStatus::OK;
}

SICALLBACK aaOcean_Term( CRef& in_ctxt )
{
    Context ctxt( in_ctxt );

    // get ocean pointer from user-data and do clean up
    CValue userData = ctxt.GetUserData( );
    aaOcean *pOcean = (aaOcean *)(CValue::siPtrType)ctxt.GetUserData();
    if(pOcean)
        delete pOcean;  
    pOcean = 0;
    ctxt.PutUserData( CValue() );
    Application().LogMessage(L"[aaOcean ICE] : Successfully cleaned up memory");

    return CStatus::OK;
}

SICALLBACK aaOcean_BeginEvaluate( ICENodeContext& in_ctxt )
{
    // get ocean pointer from user-data
    aaOcean *pOcean = (aaOcean *)(CValue::siPtrType)in_ctxt.GetUserData();

    // get ICE node input port arrays
    CDataArrayLong PointID( in_ctxt, ID_IN_PointID);
    CDataArrayLong resolution( in_ctxt, ID_IN_RESOLUTION);      
    CDataArrayLong seed( in_ctxt, ID_IN_SEED);
    CDataArrayFloat waveHeight( in_ctxt, ID_IN_WAVE_HEIGHT);
    CDataArrayFloat waveSpeed( in_ctxt, ID_IN_WAVESPEED);
    CDataArrayFloat chop( in_ctxt, ID_IN_CHOP);
    CDataArrayFloat oceanScale( in_ctxt, ID_IN_OCEAN_SCALE );
    CDataArrayFloat oceanDepth( in_ctxt, ID_IN_OCEAN_DEPTH );
    CDataArrayFloat windDir( in_ctxt, ID_IN_WINDDIR );
    CDataArrayFloat cutoff( in_ctxt, ID_IN_CUTOFF);
    CDataArrayFloat velocity( in_ctxt, ID_IN_WINDVELOCITY);
    CDataArrayLong  windAlign( in_ctxt, ID_IN_WINDALIGN );
    CDataArrayFloat damp( in_ctxt, ID_IN_DAMP);
    CDataArrayBool enableFoam( in_ctxt, ID_IN_ENABLEFOAM);
    CDataArrayFloat time( in_ctxt, ID_IN_TIME);
    CDataArrayFloat loopTime( in_ctxt, ID_IN_REPEAT_TIME);
    CDataArrayFloat surfaceTension( in_ctxt, ID_IN_SURFACE_TENSION);
    CDataArrayFloat randWeight( in_ctxt, ID_IN_RAND_WEIGHT);

    pOcean->input(resolution[0], 
        seed[0],
        oceanScale[0], 
        oceanDepth[0],
        surfaceTension[0],
        velocity[0], 
        cutoff[0], 
        windDir[0], 
        windAlign[0], 
        damp[0], 
        waveSpeed[0], 
        waveHeight[0],
        chop[0], 
        time[0],
        loopTime[0],
        enableFoam[0],
        randWeight[0]);

    return CStatus::OK;
}

SICALLBACK aaOcean_Evaluate( ICENodeContext& in_ctxt )
{
    aaOcean *pOcean = (aaOcean *)(CValue::siPtrType)in_ctxt.GetUserData();

    CIndexSet indexSet(in_ctxt, ID_IN_PointID );
    CDataArrayLong PointID(in_ctxt, ID_IN_PointID );
    CDataArrayBool bEnable(in_ctxt, ID_IN_ENABLE);
    CDataArrayFloat uCoord(in_ctxt, ID_IN_U);
    CDataArrayFloat vCoord(in_ctxt, ID_IN_V);
    CDataArrayBool enableFoam( in_ctxt, ID_IN_ENABLEFOAM);
    CDataArrayMatrix4f transform(in_ctxt, ID_IN_TRANSFORM);

    const int count = PointID.GetCount();

    float worldSpaceVec[3] = {0.0f, 0.0f, 0.0f};
    float localSpaceVec[3] = {0.0f, 0.0f, 0.0f};

    int transformArraySize = 0;
    bool transformSingleton = TRUE;

    if(transform.GetCount() > 1)
        transformSingleton = FALSE;

    ULONG out_portID = in_ctxt.GetEvaluatedOutputPortID( ); 
    switch( out_portID )
    {
        case ID_OUT_OCEAN:
        {
            CDataArrayVector3f outData( in_ctxt );
            if(bEnable[0])
            {
                #pragma omp parallel for private(worldSpaceVec, localSpaceVec)
                for(int i = 0; i<count; ++i)
                {
                    // get ocean displacement vector
                    worldSpaceVec[1] = pOcean->getOceanData(uCoord[i], vCoord[i], aaOcean::eHEIGHTFIELD);
                    if(pOcean->isChoppy())
                    {
                        worldSpaceVec[0] = pOcean->getOceanData(uCoord[i], vCoord[i], aaOcean::eCHOPX);
                        worldSpaceVec[2] = pOcean->getOceanData(uCoord[i], vCoord[i], aaOcean::eCHOPZ);
                    }

                    // multiply displacement vector by input transform matrix
                    if(!transformSingleton)
                        transformArraySize = i;

                    multiplyMatrix(&worldSpaceVec[0], &localSpaceVec[0], transform, transformArraySize);

                    outData[i].PutX(localSpaceVec[0]);
                    outData[i].PutY(localSpaceVec[1]);
                    outData[i].PutZ(localSpaceVec[2]);

                    outData[i].PutX(worldSpaceVec[0]);
                    outData[i].PutY(worldSpaceVec[1]);
                    outData[i].PutZ(worldSpaceVec[2]);
                }
            }
            else
            {
                #pragma omp parallel for
                for(int i = 0; i<count; ++i)
                {
                    outData[i].PutX(0.f); 
                    outData[i].PutY(0.f); 
                    outData[i].PutZ(0.f);
                }
            }
        }
        break;

        case ID_OUT_FOAM:
        {
            if(pOcean->isChoppy() && bEnable[0] && enableFoam[0])
            {               
                CDataArrayFloat outData( in_ctxt );

                // output raw (unscaled) foam in ICE deformer
                #pragma omp parallel for
                for(int i = 0; i<count; ++i)
                    outData[i] = pOcean->getOceanData(uCoord[i], vCoord[i], aaOcean::eFOAM);
            }
        }
        break;

        case ID_OUT_EIGEN_MINUS:
        {
            CDataArrayVector3f outData( in_ctxt );
            if(pOcean->isChoppy() && bEnable[0] && enableFoam[0])
            {
                #pragma omp parallel for private(worldSpaceVec, localSpaceVec)
                for(int i = 0; i<count; ++i)
                {
                    worldSpaceVec[0] = pOcean->getOceanData(uCoord[i], vCoord[i], aaOcean::eEIGENMINUSX);
                    worldSpaceVec[2] = pOcean->getOceanData(uCoord[i], vCoord[i], aaOcean::eEIGENMINUSZ);

                    multiplyMatrix(&worldSpaceVec[0], &localSpaceVec[0], transform, transformArraySize);

                    outData[i].PutX(localSpaceVec[0]);
                    outData[i].PutY(0.0f);
                    outData[i].PutZ(localSpaceVec[2]);
                }
            }
            else
            {
                #pragma omp parallel for
                for(int i = 0; i<count; ++i){
                    outData[i].PutX(0);outData[i].PutY(0);outData[i].PutZ(0);}
            }
        }
        break;

        case ID_OUT_EIGEN_PLUS:
        {
            CDataArrayVector3f outData( in_ctxt );
            if(pOcean->isChoppy() && bEnable[0] && enableFoam[0])
            {
                #pragma omp parallel for private(worldSpaceVec, localSpaceVec)
                for(int i = 0; i<count; ++i)
                {
                    worldSpaceVec[0] = pOcean->getOceanData(uCoord[i], vCoord[i], aaOcean::eEIGENPLUSX);
                    worldSpaceVec[2] = pOcean->getOceanData(uCoord[i], vCoord[i], aaOcean::eEIGENPLUSZ);

                    multiplyMatrix(&worldSpaceVec[0], &localSpaceVec[0], transform, transformArraySize);

                    outData[i].PutX(localSpaceVec[0]);
                    outData[i].PutY(0.0f);
                    outData[i].PutZ(localSpaceVec[2]);
                }
            }
            else
            {
                #pragma omp parallel for
                for(int i = 0; i<count; ++i){
                    outData[i].PutX(0); outData[i].PutY(0); outData[i].PutZ(0);}
            }
        }
        break;
    }

    return CStatus::OK;
}