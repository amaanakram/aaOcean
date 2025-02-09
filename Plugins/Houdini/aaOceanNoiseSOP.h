// aaOceanNoise
// Author: Amaan Akram 
// https://linkedin.com/in/amaan
//
#ifndef __aaOceanNoise_h__
#define __aaOceanNoise_h__

#include <SOP/SOP_Node.h>
#include "perlin_noise.h"
#include "timer.h"

class aaOceanNoiseSOP : public SOP_Node
{
public:
    aaOceanNoiseSOP(OP_Network *net, const char *name, OP_Operator *op);
    virtual ~aaOceanNoiseSOP();

    static PRM_Template myTemplateList[];
    static OP_Node *myConstructor(OP_Network*, const char *, OP_Operator *);

protected:
    virtual const char  *inputLabel(unsigned idx) const;

    /// Method to cook geometry for the SOP
    virtual OP_ERROR cookMySop(OP_Context &context);

private:
    // working variables
    GA_ROAttributeRef uvRef;
    const GA_Attribute *uvAttribute;
    const GA_AIFTuple *uvTuple;

    GA_RWHandleF aaNoiseHandle;
    GA_RWAttributeRef aaNoiseRef;
    char aaNoiseName[64];

    UT_String   UvAttribute;
    
    int     SEED(fpreal t)          { return evalInt("seed", 0, t); }
    int     OCTAVES(fpreal t)       { return evalInt("octaves", 0, t); }
    fpreal  FREQUENCY(fpreal t)   { return evalFloat("frequency", 0, t); }
    fpreal  U_MULT(fpreal t)        { return evalFloat("u_mult", 0, t); }
    fpreal  V_MULT(fpreal t)        { return evalFloat("v_mult", 0, t); }
    fpreal  GAMMA(fpreal t)         { return evalFloat("gamma", 0, t); }
    fpreal  AMPLITUDE(fpreal t)     { return evalFloat("amplitude", 0, t); }
    fpreal  TIME(fpreal t)          { return evalFloat("time", 0, t); }

    void    getUVAttributeName(UT_String &str){ evalString(str, "uvAttribute", 0, 0); }
    
    /// This variable is used together with the call to the "checkInputChanged"
    /// routine to notify the handles (if any) if the input has changed.
    GU_DetailGroupPair   myDetailGroupPair;

    PerlinNoise *pNoise;
};

#endif