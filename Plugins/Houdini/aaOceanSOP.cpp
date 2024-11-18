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

#include "aaOceanSOP.h"
#include <iostream>
#include <SYS/SYS_Math.h>
#include <UT/UT_DSOVersion.h>
#include <UT/UT_Interrupt.h>
#include <UT/UT_Matrix3.h>
#include <UT/UT_Matrix4.h>
#include <UT/UT_ThreadedAlgorithm.h>
#include <GU/GU_Detail.h>
#include <GU/GU_PrimPoly.h>
#include <PRM/PRM_Include.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <SOP/SOP_Guide.h>
#include <UT/UT_Options.h>
#include <GEO/GEO_AttributeHandle.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Include.h>

void newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(new OP_Operator("aaOceanSOP",
        "aaOceanSOP",
        aaOceanSOP::myConstructor,
        aaOceanSOP::myTemplateList,
        1,
        1,
        0));
}

static PRM_Name spectrumNames[] =
{
    PRM_Name("Philips",  "Philips"),
    PRM_Name("Pierson-Morkowitz", "Pierson-Morkowitz"),
    PRM_Name("TMA",  "TMA"),
    PRM_Name("JONSWAP",  "JONSWAP"),
    PRM_Name(0)
};
static PRM_ChoiceList spectrumNamesMenu(PRM_CHOICELIST_SINGLE, spectrumNames);

static PRM_Name resolutionList[] =
{
    PRM_Name("2", "64x64"),
    PRM_Name("3", "128x128"),
    PRM_Name("4", "256x256"),
    PRM_Name("5", "512x512"),
    PRM_Name("6", "1024x1024"),
    PRM_Name("7", "2048x2048"),
    PRM_Name(0)
};
static PRM_ChoiceList resolutionListMenu(PRM_CHOICELIST_SINGLE, resolutionList);

static PRM_Name names[] = 
{
    PRM_Name("resolution",      "Resolution"),
    PRM_Name("seed",            "Seed"),
    PRM_Name("oceanScale",      "Ocean Scale"),
    PRM_Name("oceanDepth",      "Ocean Depth"),
    PRM_Name("surfaceTension",  "Surface Tension"),
    PRM_Name("velocity",        "Wave Size/Wind Velocity"),
    PRM_Name("cutoff",          "Wave Smooth"),
    PRM_Name("windDir",         "Wind Dir"),
    PRM_Name("windAlign",       "Wind Align"),
    PRM_Name("damp",            "Damp Reflected Waves"),
    PRM_Name("waveSpeed",       "Wave Speed"),
    PRM_Name("waveHeight",      "Wave Height"),
    PRM_Name("chop",            "Chop Amount"),
    PRM_Name("enableEigens",    "Output Eigens & Spectrum Attributes"),
    PRM_Name("timeOffset",      "Time Offset"),
    PRM_Name("loopTime",        "Loop Time"),
    PRM_Name("uvAttribute",     "UV Attribute"),
    PRM_Name("spectrum",        "Spectrum Type"),
    PRM_Name("spectrumMult",    "Spectrum Multiplier"),
    PRM_Name("peakSharpening",   "Peak Sharpening"),
    PRM_Name("fetch",           "TMA Fetch"),
    PRM_Name("swellAmount",     "Swell Amount"),
    PRM_Name("randWeight",      "Random Weight"),
    PRM_Name("time",            "Time"),
};

// defining some custom ranges and defaults

static PRM_Range        resolutionRange(PRM_RANGE_RESTRICTED, 1, PRM_RANGE_RESTRICTED, 8);
static PRM_Default      resolutionDefault(4);

static PRM_Range        oceanScaleRange(PRM_RANGE_RESTRICTED, 0.0, PRM_RANGE_UI, 200.0);
static PRM_Default      oceanScaleDefault(100.0);

static PRM_Range        oceanDepthRange(PRM_RANGE_RESTRICTED, 1.0, PRM_RANGE_RESTRICTED, 10000.0);
static PRM_Default      oceanDepthDefault(10000.0);

static PRM_Range        seedRange(PRM_RANGE_RESTRICTED, 1, PRM_RANGE_RESTRICTED, 15);
static PRM_Default      seedDefault(1);

static PRM_Range        velocityRange(PRM_RANGE_RESTRICTED, 0.0, PRM_RANGE_UI, 300.0);
static PRM_Default      velocityDefault(4.0);

static PRM_Range        waveSpeedRange(PRM_RANGE_UI, -10.0, PRM_RANGE_UI, 10.0);
static PRM_Default      waveSpeedDefault(1.0);

static PRM_Range        timeRange(PRM_RANGE_UI, 0, PRM_RANGE_UI, 1000.0);
static PRM_Default      timeRangeDefault(1000.0);

static PRM_Range        loopTimeRange(PRM_RANGE_RESTRICTED, 0.001, PRM_RANGE_UI, 1000.0);
static PRM_Default      loopTimeDefault(1000.0);

static PRM_Range        randWeightRange(PRM_RANGE_RESTRICTED, 0.0, PRM_RANGE_UI, 1.f);
static PRM_Default      randWeightDefault(0.0);

static PRM_Range        spectrumRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_RESTRICTED, 2);
static PRM_Default      spectrumDefault(0);

static PRM_Range        surftRange(PRM_RANGE_RESTRICTED, 0.001, PRM_RANGE_UI, 100.0);
static PRM_Default      surftDefault(1.0);

static PRM_Range        spectrumMultRange(PRM_RANGE_RESTRICTED, 0.001, PRM_RANGE_UI, 100.0);
static PRM_Default      spectrumMultDefault(1.0);

static PRM_Range        peakSharpeningRange(PRM_RANGE_RESTRICTED, 0.001, PRM_RANGE_UI, 6.0);
static PRM_Default      peakSharpeningDefault(1.0);

static PRM_Range        fetchRange(PRM_RANGE_RESTRICTED, 0.0001, PRM_RANGE_UI, 1000.0);
static PRM_Default      fetchDefault(100.0);

static PRM_Range        swellRange(PRM_RANGE_RESTRICTED, 0.0, PRM_RANGE_UI, 1.f);
static PRM_Default      swellDefault(0.0);


static PRM_Name		switcherName("oceanTabs");
static PRM_Default	switcher[] = {
    PRM_Default(3, "Wave Params"),	    
    PRM_Default(4, "Wind Params"),	    
    PRM_Default(3, "TMA-specific Params"),
    PRM_Default(7, "Misc"),	    
};

PRM_Template aaOceanSOP::myTemplateList[] = 
{   
    //PRM_Template(PRM_INT_E, 1, &names[0],  &resolutionDefault,  0, &resolutionRange),           // resolution   // 0
    PRM_Template(PRM_ORD, PRM_Template::PRM_EXPORT_MAX, 1, &names[0],  0, &resolutionListMenu),   // resolution   // 0
    PRM_Template(PRM_ORD, PRM_Template::PRM_EXPORT_MAX, 1, &names[17], 0, &spectrumNamesMenu),  // spectrum type
    PRM_Template(PRM_FLT_J, 1, &names[2],  &oceanScaleDefault,  0, &oceanScaleRange),           // oceanScale   // 2
    PRM_Template(PRM_INT_E, 1, &names[1],  &seedDefault,        0, &seedRange),                 // seed         // 1
    PRM_Template(PRM_FLT_J, 1, &names[23], PRMzeroDefaults,     0, &timeRange),                 // time         // 23
    PRM_Template(PRM_FLT_J, 1, &names[14], PRMzeroDefaults,     0, &PRMscaleRange),             // timeOffset   // 14
    PRM_Template(PRM_TOGGLE,1, &names[13]),                                                     // enable Foam  // 13

    PRM_Template(PRM_SWITCHER,  sizeof(switcher) / sizeof(PRM_Default), &switcherName, switcher),

    PRM_Template(PRM_FLT_J, 1, &names[11], PRMoneDefaults,      0, &PRMdivision0Range),     // waveHeight   // 11
    PRM_Template(PRM_FLT_J, 1, &names[10], &waveSpeedDefault,   0, &waveSpeedRange),        // waveSpeed    // 10
    PRM_Template(PRM_FLT_J, 1, &names[12], PRMzeroDefaults,     0, &PRMrolloffRange),       // chop         // 12

    PRM_Template(PRM_FLT_J, 1, &names[5],  &velocityDefault,    0, &velocityRange),         // velocity (Wave Size) //5
    PRM_Template(PRM_FLT_J, 1, &names[7],  PRMzeroDefaults,     0, &PRMangleRange),         // windDir      // 7
    PRM_Template(PRM_FLT_J, 1, &names[9],  PRMzeroDefaults,     0, &PRMunitRange),          // damp         // 9
    PRM_Template(PRM_INT_E, 1, &names[8],  PRMoneDefaults,      0, &PRMrolloffRange),       // windAlign    // 8

    PRM_Template(PRM_FLT_J, 1, &names[19], &peakSharpeningDefault,0, &peakSharpeningRange), // peak sharpening// 19
    PRM_Template(PRM_FLT_J, 1, &names[20], &fetchDefault,        0, &fetchRange),           // fetch // 20
    PRM_Template(PRM_FLT_J, 1, &names[21], &swellDefault,        0, &swellRange),           // swell // 21

    PRM_Template(PRM_FLT_J, 1, &names[3],  &oceanDepthDefault,  0, &oceanDepthRange),       // oceanDepth   // 3
    PRM_Template(PRM_FLT_J, 1, &names[4],  &surftDefault,       0, &surftRange),            // surfaceTension// 4
    PRM_Template(PRM_FLT_J, 1, &names[6],  PRMzeroDefaults,     0, &PRMdivision0Range),     // cutoff (Wave Smooth) // 6
    PRM_Template(PRM_STRING,1, &names[16], 0),                                              // UV Attribute // 16
    PRM_Template(PRM_FLT_J, 1, &names[15], &loopTimeDefault,    0, &loopTimeRange),         // loop time    // 15
    PRM_Template(PRM_FLT_J, 1, &names[18], &spectrumMultDefault, 0, &spectrumMultRange),    // spectrumMult // 18
    PRM_Template(PRM_FLT_J, 1, &names[22], &randWeightDefault,   0, &randWeightRange),      // randWeight //

    PRM_Template(),
};

OP_Node *aaOceanSOP::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new aaOceanSOP(net, name, op);
}

aaOceanSOP::aaOceanSOP(OP_Network *net, const char *name, OP_Operator *op)
    : SOP_Node(net, name, op)
{
    sprintf(eVecPlusName,   "eVecPlus");
    sprintf(eVecMinusName,  "eVecMinus");
    sprintf(eValuesName,    "eValues");
    sprintf(spectrumName,   "spectrum");
    enableEigens = FALSE;

    pOcean = new aaOcean;
}

aaOceanSOP::~aaOceanSOP() 
{
    if(pOcean)
        delete pOcean;
}

OP_ERROR aaOceanSOP::cookMySop(OP_Context &context)
{
    if (lockInputs(context) >= UT_ERROR_ABORT)
        return error();

    duplicateSource(0, context);
    setVariableOrder(3, 2, 0, 1);
    setCurGdh(0, myGdpHandle);
    setupLocalVars();

    // Flag the SOP as being time dependent (i.e. cook on time changes)
    flags().setTimeDep(true);

    // get the user-specified attribute that holds uv-data
    getUVAttributeName(UvAttribute);
    if(UvAttribute.length() == 0)
        UvAttribute = "uv";
    const char* UVAttribName = (const char *)UvAttribute;
    bool isVertexContext = false;

    // assumming that UVs are defined per-Vertex
    GA_ROHandleV3 UVHandle(gdp->findVertexAttribute(UVAttribName));

    if (UVHandle.isValid())
        isVertexContext = true;
    else
    {
        UVHandle = gdp->findPointAttribute(UVAttribName);
        if(!UVHandle.isValid())
        {
            // uv attribute not found
            char msg[512];
            snprintf(msg, sizeof(msg),
            "[aaOcean SOP] Specified UV Attribute '%s' not found on input geometry. "
            "Please add in either Vertex or Point class\n",
            UVAttribName);

            std::cout<<msg;
            std::cout.flush();
            addError(SOP_ERR_INVALID_SRC, msg); 
            unlockInputs();
            return error();
        }
    }

    // inputs validated. Begin writing ocean data to output handles
    // start pulling in SOP inputs and send to aaOcean 
    enableEigens = (ENABLEEIGENS() != 0);
    if(pOcean->isChoppy() && enableEigens)
        enableEigens = TRUE;
    
    // variable declarations
    float now  = context.getTime();

    pOcean->input(  RESOLUTION() + 2,
                    SPECTRUM(),
					SEED(),
                    OCEANSCALE(now),
                    OCEANDEPTH(now),
                    SURFACETENSION(now),
                    VELOCITY(now), 
                    CUTOFF(now), 
                    WINDDIR(now), 
                    WINDALIGN(), 
                    DAMP(now), 
                    WAVESPEED(now), 
                    WAVEHEIGHT(now),
                    CHOP(now), 
                    TIME(now) + TIMEOFFSET(now),
                    LOOPTIME(now),
                    enableEigens,
                    RANDWEIGHT(now),
                    SPECTRUMMULT(now),
                    PEAKSHARPENING(now),
                    FETCH(now),
                    SWELLAMOUNT(now));

    // setup local variables to output Eigens
    if(enableEigens)
    {
        eVecPlusRef  = gdp->addFloatTuple(GA_ATTRIB_POINT, eVecPlusName,    3);
        eVecMinusRef = gdp->addFloatTuple(GA_ATTRIB_POINT, eVecMinusName,   3);
        eValuesRef   = gdp->addFloatTuple(GA_ATTRIB_POINT, eValuesName,     1);

        eVecPlusHandle  = GA_RWHandleV3(eVecPlusRef.getAttribute());
        eVecMinusHandle = GA_RWHandleV3(eVecMinusRef.getAttribute());
        eValuesHandle   = GA_RWHandleF(eValuesRef.getAttribute());

        spectrumRef = gdp->addFloatTuple(GA_ATTRIB_POINT, spectrumName, 1);
        spectrumHandle = GA_RWHandleF(spectrumRef.getAttribute());
    }
    
	//#pragma omp parallel for 
    for (size_t pt_index = 0; pt_index < gdp->getNumPoints(); ++pt_index)
    {
        UT_Vector3F pos = gdp->getPos3(pt_index);
        UT_Vector3F uv;

        if(isVertexContext)
        {
            // convert from Point to Vertex context
            GA_Offset vtx_index = gdp->pointVertex(pt_index);   
            uv = UVHandle.get(vtx_index);
        }
        else
            uv = UVHandle.get(pt_index);

        const float u = uv.x();
        const float v = uv.y(); 
        
        pos.y() += pOcean->getOceanData(u, v, aaOcean::eHEIGHTFIELD);

        if(pOcean->isChoppy())
        {
            pos.x() += pOcean->getOceanData(u, v, aaOcean::eCHOPX);
            pos.z() += pOcean->getOceanData(u, v, aaOcean::eCHOPZ);
        }
        
        gdp->setPos3(pt_index, pos); // not thread-safe

        if(enableEigens)
        {
            UT_Vector3F eigenVectorPlusValue;
            UT_Vector3F eigenVectorMinusValue;
            
            eigenVectorPlusValue.x() =  pOcean->getOceanData(u, v, aaOcean::eEIGENPLUSX);
            eigenVectorPlusValue.y() =  0.0f;
            eigenVectorPlusValue.z() =  pOcean->getOceanData(u, v, aaOcean::eEIGENPLUSZ);

            eigenVectorMinusValue.x() = pOcean->getOceanData(u, v, aaOcean::eEIGENMINUSX);
            eigenVectorMinusValue.y() = 0.0f;
            eigenVectorMinusValue.z() = pOcean->getOceanData(u, v, aaOcean::eEIGENMINUSZ);

            eVecPlusHandle.set(pt_index, eigenVectorPlusValue);
            eVecMinusHandle.set(pt_index, eigenVectorMinusValue);

            float eigenValue = pOcean->getOceanData(u, v, aaOcean::eFOAM);
            eValuesHandle.set(pt_index,eigenValue);

            float spectrumValue = pOcean->getOceanData(u, v, aaOcean::eSPECTRUM);
            spectrumHandle.set(pt_index, spectrumValue);
        }
    }
    unlockInputs();
    return error();
}

const char *aaOceanSOP::inputLabel(unsigned) const
{
    return "UV'ed Geometry to simulate ocean on";
}
