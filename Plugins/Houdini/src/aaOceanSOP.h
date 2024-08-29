// aaOcean v2.6 Houdini SOP Deformer
// Author: Amaan Akram 
// www.amaanakram.com
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html

#ifndef __aaOceanSOP_h__
#define __aaOceanSOP_h__

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <SOP/SOP_Node.h>
#include "aaOceanClass.cpp"

class aaOceanSOP : public SOP_Node
{
public:
    aaOceanSOP(OP_Network *net, const char *name, OP_Operator *op);
    virtual ~aaOceanSOP();

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

    GA_RWAttributeRef eVecPlusRef;
    GA_RWAttributeRef eVecMinusRef;
    GA_RWAttributeRef eValuesRef;
    GA_RWAttributeRef spectrumRef;
    GA_RWHandleV3 eVecPlusHandle;
    GA_RWHandleV3 eVecMinusHandle;
    GA_RWHandleF eValuesHandle;
    GA_RWHandleF spectrumHandle;

    bool enableEigens;
    char eVecPlusName[10];
    char eVecMinusName[10];
    char eValuesName[10];
    char spectrumName[15];
    UT_String   UvAttribute;
    
    int     RESOLUTION()            { return evalInt("resolution", 0, 0); }
    int     SPECTRUM()              { return evalInt("spectrum", 0, 0); }
    int     SEED()                  { return evalInt("seed", 0, 0); }
    fpreal  OCEANSCALE(fpreal t)    { return evalFloat("oceanScale", 0, t); }
    fpreal  OCEANDEPTH(fpreal t)    { return evalFloat("oceanDepth", 0, t); }
    fpreal  SURFACETENSION(fpreal t){ return evalFloat("surfaceTension", 0, t); }

    fpreal  VELOCITY(fpreal t)      { return evalFloat("velocity", 0, t); }
    fpreal  CUTOFF(fpreal t)        { return evalFloat("cutoff", 0, t); }
    fpreal  WINDDIR(fpreal t)       { return evalFloat("windDir", 0, t); }
    int     WINDALIGN()             { return evalInt("windAlign", 0, 0); }

    fpreal  DAMP(fpreal t)          { return evalFloat("damp", 0, t); }
    fpreal  WAVESPEED(fpreal t)     { return evalFloat("waveSpeed", 0, t); }
    fpreal  WAVEHEIGHT(fpreal t)    { return evalFloat("waveHeight", 0, t); }
    fpreal  CHOP(fpreal t)          { return evalFloat("chop", 0, t); }
    int     ENABLEEIGENS()          { return evalInt("enableEigens", 0, 0); }
    fpreal  TIMEOFFSET(fpreal t)    { return evalFloat("timeOffset", 0, t); }
    fpreal  LOOPTIME(fpreal t)      { return evalFloat("loopTime", 0, t); }
    fpreal  RANDWEIGHT(fpreal t)    { return evalFloat("randWeight", 0, t); }
    fpreal  SPECTRUMMULT(fpreal t)  { return evalFloat("spectrumMult", 0, t); }
    fpreal  PEAKSHARPENING(fpreal t){ return evalFloat("peakSharpening", 0, t); }
    fpreal  FETCH(fpreal t)         { return evalFloat("fetch", 0, t); }
    fpreal  SWELLAMOUNT(fpreal t)   { return evalFloat("swellAmount", 0, t); }

    void    getUVAttributeName(UT_String &str){ evalString(str, "uvAttribute", 0, 0); }
    void    getEigenPlusAttribute(UT_String &str){ evalString(str, "eigenPlusAttribute", 0, 0); }
    void    getEigenMinusAttribute(UT_String &str){ evalString(str, "eigenMinusAttribute", 0, 0); }
    void    getEigenValueAttribute(UT_String &str){ evalString(str, "eigenValueAttribute", 0, 0); }
    void    getSpectrumAttribute(UT_String &str) { evalString(str, "spectrumAttribute", 0, 0); }
    
    /// This variable is used together with the call to the "checkInputChanged"
    /// routine to notify the handles (if any) if the input has changed.
    GU_DetailGroupPair   myDetailGroupPair;

    // pointer to our aaOcean class
    aaOcean *pOcean;
};

#endif