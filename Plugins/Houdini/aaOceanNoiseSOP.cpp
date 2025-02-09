// aaFractalNoise
// Author: Amaan Akram 
// https://linkedin.com/in/amaan


#include "aaOceanNoiseSOP.h"
#include <iostream>
#include <SYS/SYS_Math.h>
#include <UT/UT_DSOVersion.h>
#include <UT/UT_Interrupt.h>
#include <GU/GU_Detail.h>
#include <GU/GU_PrimPoly.h>
#include <PRM/PRM_Include.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <SOP/SOP_Guide.h>
#include <UT/UT_Options.h>
#include <GEO/GEO_AttributeHandle.h>

void newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(new OP_Operator("aaOceanNoiseSOP",
        "aaOceanNoiseSOP",
        aaOceanNoiseSOP::myConstructor,
        aaOceanNoiseSOP::myTemplateList,
        1,
        1,
        0));
}

static PRM_Name names[] = 
{
    PRM_Name("seed",            "seed"),
    PRM_Name("octaves",         "octaves"),
    PRM_Name("frequency",       "Frequency"),
    PRM_Name("u_mult",          "u mult"),
    PRM_Name("v_mult",          "v mult"),
    PRM_Name("amplitude",       "amplitude"),
    PRM_Name("gamma",           "gamma"),
    PRM_Name("time",            "Time"),
    PRM_Name("uvAttribute",     "Input UV Attrib"),
    PRM_Name("outputVariable",  "Output Variable"),
};

// defining some custom ranges and defaults

static PRM_Range        seedRange(PRM_RANGE_RESTRICTED, 1, PRM_RANGE_UI, 15);
static PRM_Default      seedDefault(1);

static PRM_Range        octavesRange(PRM_RANGE_RESTRICTED, 1, PRM_RANGE_UI, 6);
static PRM_Default      octavesDefault(3);

static PRM_Range        freqRange(PRM_RANGE_RESTRICTED, 0.01, PRM_RANGE_UI, 100);
static PRM_Default      freqDefault(10);

static PRM_Range        defaultRange(PRM_RANGE_RESTRICTED, 0.0, PRM_RANGE_UI, 10.0);
static PRM_Default      defaultOne(1.0);

static PRM_Default uvDefault(0, "uv");
static PRM_Default outputDefault(0, "f_noise");

PRM_Template aaOceanNoiseSOP::myTemplateList[] = 
{   
    PRM_Template(PRM_INT_E, 1, &names[0],  &seedDefault,        0, &seedRange),         // seed        
    PRM_Template(PRM_INT_E, 1, &names[1],  &octavesDefault,     0, &octavesRange),      // octaves        
    PRM_Template(PRM_FLT_J, 1, &names[2],  &freqDefault,        0, &freqRange),         // frequency
    PRM_Template(PRM_FLT_J, 1, &names[3],  &defaultOne,         0, &defaultRange),      // u_mult
    PRM_Template(PRM_FLT_J, 1, &names[4],  &defaultOne,         0, &defaultRange),      // v_mult
    PRM_Template(PRM_FLT_J, 1, &names[5],  &defaultOne,         0, &defaultRange),      // amplitude
    PRM_Template(PRM_FLT_J, 1, &names[6],  &defaultOne,         0, &defaultRange),      // gamma
    PRM_Template(PRM_FLT_J, 1, &names[7],  &defaultOne,         0, &defaultRange),      // time
    PRM_Template(PRM_STRING,1, &names[8],  &uvDefault),
    PRM_Template(PRM_STRING,1, &names[9],  &outputDefault),

    PRM_Template(),
};

OP_Node *aaOceanNoiseSOP::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new aaOceanNoiseSOP(net, name, op);
}

aaOceanNoiseSOP::aaOceanNoiseSOP(OP_Network *net, const char *name, OP_Operator *op)
    : SOP_Node(net, name, op)
{
    sprintf(aaNoiseName, "f_noise");
    pNoise = nullptr;
}

aaOceanNoiseSOP::~aaOceanNoiseSOP() 
{
    if(pNoise)
        delete pNoise;
}

OP_ERROR aaOceanNoiseSOP::cookMySop(OP_Context &context)
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
            "[aaFractalNoise SOP] Specified UV Attribute '%s' not found on input geometry. "
            "Please add in either Vertex or Point class\n",
            UVAttribName);

            std::cout<<msg;
            std::cout.flush();
            addError(SOP_ERR_INVALID_SRC, msg); 
            unlockInputs();
            return error();
        }
    }

    UT_String output;
    evalString(output, "outputVariable", 0, context.getTime());
    if (output.isstring())
        sprintf(aaNoiseName, output.c_str());

    // variable declarations
    float now       = context.getTime();
    aaNoiseRef      = gdp->addFloatTuple(GA_ATTRIB_POINT, aaNoiseName, 1);
    aaNoiseHandle   = GA_RWHandleF(aaNoiseRef.getAttribute());

    float u_mult = U_MULT(now);
    float v_mult = V_MULT(now);
    float gamma = GAMMA(now);
    float amplitude = AMPLITUDE(now);
    int octaves = OCTAVES(now);
    int seed = SEED(now);
    float time = TIME(now);
    float frequency = FREQUENCY(now);

    Timer timer;
    if(pNoise == nullptr)
        pNoise = new PerlinNoise(seed);
    else if (pNoise->getSeed() != seed)
    {
        delete pNoise;
        pNoise = new PerlinNoise(seed);
    }

    GA_Offset pt_offset;
    GA_FOR_ALL_PTOFF(gdp, pt_offset)
    {
        UT_Vector3F uv;

        if(isVertexContext)
        {
            // convert from Point to Vertex context
            GA_Offset vtx_index = gdp->pointVertex(pt_offset);   
            if(vtx_index == GA_INVALID_OFFSET)
                continue;
            uv = UVHandle.get(vtx_index);
        }
        else
            uv = UVHandle.get(pt_offset);
        
        float u = uv.x() * u_mult;
        float v = uv.y() * v_mult; 
        float noise = pNoise->octaveNoise(u, v, time, octaves, frequency);
        noise = (noise + 1.f) * 0.5f;
        noise = std::pow(noise, gamma) * amplitude;
        aaNoiseHandle.set(pt_offset, noise);
    }

    timer.printElapsed("completed iteration");
    unlockInputs();
    return error();
}

const char *aaOceanNoiseSOP::inputLabel(unsigned) const
{
    return "UV'ed Geometry to generate aaNoise on";
}