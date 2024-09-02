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

#include <string.h>
#include <maya/MIOStream.h>
#include <math.h>
#include <maya/MGlobal.h>

#include <maya/MPxDeformerNode.h> 
#include <maya/MItGeometry.h>
#include <maya/MPxLocatorNode.h> 

#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnUnitAttribute.h>

#include <maya/MFnPlugin.h>
#include <maya/MFnDependencyNode.h>

#include <maya/MTypeId.h> 
#include <maya/MPlug.h>

#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MArrayDataHandle.h>

#include <maya/MPoint.h>
#include <maya/MVector.h>
#include <maya/MMatrix.h>

#include <maya/MAnimControl.h>
#include <maya/MPointArray.h>
#include <maya/MFloatVectorArray.h>
#include <maya/MFnStringData.h>

#include <maya/MFloatMatrix.h>
#include <maya/MFnMatrixData.h>
#include <maya/MFloatArray.h>
#include <maya/MGeometryManager.h>
#include <maya/MGeometry.h>
#include <maya/MGeometryData.h>
#include <maya/MGeometryPrimitive.h>
#include <maya/MFnMesh.h>
#include <maya/MThreadUtils.h>

#include "aaOceanClass.cpp"

//Maya Node ID 0x20B6EF34  (548859700)  -- Randomly generated. Change if this conflicts
#define MAYA_NODE_ID 0x20B6EF34

#ifdef __GNUC__
    #if __GNUC_MINOR__ > 2 
        #warning [aaOcean] compiler version may generate code incompatible with Maya
    #endif
#endif

class aaOceanDeformer : public MPxDeformerNode
{
public:
                        aaOceanDeformer();
    virtual             ~aaOceanDeformer();
    static  void*       creator();
    static  MStatus     initialize();

    //virtual void     postConstructor();
    virtual MStatus  compute( const MPlug& plug, MDataBlock& data );
    bool getUVs(MFnMesh &mesh, MDataBlock &block);
    void getColorSets(MFnMesh &mesh, MDataBlock &block);
    void setColorSets(MFnMesh &mesh, MDataBlock &block);

    static  MObject  resolution;    
	static  MObject  spectrum;
    static  MObject  oceanSize;
    static  MObject  oceanDepth;
    static  MObject  surfaceTension;
    static  MObject  seed;  

    static  MObject  waveHeight;
    static  MObject  waveSize;
    static  MObject  waveSpeed;
    static  MObject  waveChop;
    static  MObject  waveSmooth;

    static  MObject  waveDirection;
    static  MObject  waveReflection;
    static  MObject  waveAlign;
    static  MObject  time;
    static  MObject  repeatTime;
    static  MObject  timeOffset;
    static  MObject  doFoam;
    static  MObject  invertFoam;
    static  MObject  spectrumMult;
    static  MObject  peakSharpening;
    static  MObject  jswpFetch;
    static  MObject  swell;

    static  MObject  uvMap;
    static  MObject  eigenVectorMap;
    static  MObject  eigenValueMap;

    static  MObject  inTransform;

    static  MTypeId  id;

    aaOcean* pOcean;
    
	MPointArray verts;
    MFloatArray u;
    MFloatArray v;
    MColorArray colArrayEigenVector;
    MColorArray colArrayEigenValue;
    bool foundEigenVector;
    bool foundEigenValue;
    MIntArray faceColorID;

	int numPolygons;
};

MObject     aaOceanDeformer::resolution;    
MObject     aaOceanDeformer::spectrum;
MObject     aaOceanDeformer::oceanSize; 
MObject     aaOceanDeformer::oceanDepth;    
MObject     aaOceanDeformer::surfaceTension;
MObject     aaOceanDeformer::seed;  
MObject     aaOceanDeformer::waveHeight;
MObject     aaOceanDeformer::waveSize;
MObject     aaOceanDeformer::waveSpeed;
MObject     aaOceanDeformer::waveChop;
MObject     aaOceanDeformer::waveSmooth;
MObject     aaOceanDeformer::waveDirection;
MObject     aaOceanDeformer::waveReflection;
MObject     aaOceanDeformer::waveAlign;
MObject     aaOceanDeformer::time;
MObject     aaOceanDeformer::repeatTime;
MObject     aaOceanDeformer::timeOffset;
MObject     aaOceanDeformer::doFoam;
MObject     aaOceanDeformer::invertFoam;
MObject     aaOceanDeformer::uvMap;
MObject     aaOceanDeformer::eigenVectorMap;
MObject     aaOceanDeformer::eigenValueMap;
MObject     aaOceanDeformer::inTransform;
MObject     aaOceanDeformer::spectrumMult;
MObject     aaOceanDeformer::peakSharpening;
MObject     aaOceanDeformer::jswpFetch;
MObject     aaOceanDeformer::swell;

MTypeId     aaOceanDeformer::id(MAYA_NODE_ID); //Maya Node ID 548859700

MStatus aaOceanDeformer::initialize()
{
    MFnNumericAttribute nAttrResolution;
    resolution = nAttrResolution.create( "Resolution", "resolution", MFnNumericData::kInt, 2 );
    nAttrResolution.setKeyable (true);  
    nAttrResolution.setWritable(true);
    nAttrResolution.setSoftMin(1);  
    nAttrResolution.setSoftMax(6);  
    nAttrResolution.setMin(1);  
    nAttrResolution.setMax(6);  
    addAttribute( resolution );
    attributeAffects( aaOceanDeformer::resolution, aaOceanDeformer::outputGeom);

	MFnNumericAttribute nAttrSpectrum;
	spectrum = nAttrSpectrum.create("Spectrum", "spectrum", MFnNumericData::kInt, 0);
	nAttrSpectrum.setKeyable(true);
	nAttrSpectrum.setWritable(true);
	nAttrSpectrum.setSoftMin(0);
	nAttrSpectrum.setSoftMax(3);
	nAttrSpectrum.setMin(0);
	nAttrSpectrum.setMax(3);
	addAttribute(spectrum);
	attributeAffects(aaOceanDeformer::spectrum, aaOceanDeformer::outputGeom);

    MFnNumericAttribute nAttrSpectrumMult;
    spectrumMult = nAttrSpectrumMult.create("spectrumMult", "spectrumMult", MFnNumericData::kFloat, 1.f);
    nAttrSpectrumMult.setKeyable(true);
    nAttrSpectrumMult.setWritable(true);
    nAttrSpectrumMult.setSoftMin(1.f);
    nAttrSpectrumMult.setSoftMax(100.f);
    nAttrSpectrumMult.setMin(1.f);
    addAttribute(spectrumMult);
    attributeAffects(aaOceanDeformer::spectrumMult, aaOceanDeformer::outputGeom);

    MFnNumericAttribute nAttrPeakSharpening;
    peakSharpening = nAttrPeakSharpening.create("PeakSharpening", "peakSharpening", MFnNumericData::kFloat, 1.f);
    nAttrPeakSharpening.setKeyable(true);
    nAttrPeakSharpening.setWritable(true);
    nAttrPeakSharpening.setMin(0.001f);
    nAttrPeakSharpening.setMax(6.f);
    addAttribute(peakSharpening);
    attributeAffects(aaOceanDeformer::peakSharpening, aaOceanDeformer::outputGeom);

    MFnNumericAttribute nAttrJswpFetch;
    jswpFetch = nAttrJswpFetch.create("jswpFetch", "jswpFetch", MFnNumericData::kFloat, 20.f);
    nAttrJswpFetch.setKeyable(true);
    nAttrJswpFetch.setWritable(true);
    nAttrJswpFetch.setMin(0.001f);
    nAttrJswpFetch.setMax(1000.f);
    addAttribute(jswpFetch);
    attributeAffects(aaOceanDeformer::jswpFetch, aaOceanDeformer::outputGeom);

    MFnNumericAttribute nAttrSwell;
    swell = nAttrPeakSharpening.create("Swell", "swell", MFnNumericData::kFloat, 0.f);
    nAttrSwell.setKeyable(true);
    nAttrSwell.setWritable(true);
    nAttrSwell.setMin(0.0f);
    nAttrSwell.setMax(1.0f);
    addAttribute(swell);
    attributeAffects(aaOceanDeformer::swell, aaOceanDeformer::outputGeom);

    MFnNumericAttribute nAttrOceanSize;
    oceanSize= nAttrOceanSize.create( "oceanSize", "oceanSize", MFnNumericData::kFloat, 100.f );
    nAttrOceanSize.setKeyable(  true ); 
    nAttrOceanSize.setWritable(true);
    nAttrOceanSize.setSoftMin(1.f); 
    nAttrOceanSize.setSoftMax(1000.f);      
    addAttribute( oceanSize );
    attributeAffects( aaOceanDeformer::oceanSize, aaOceanDeformer::outputGeom);

    MFnNumericAttribute nAttrOceanDepth;
    oceanDepth= nAttrOceanDepth.create( "oceanDepth", "oceanDepth", MFnNumericData::kFloat, 10000.f );
    nAttrOceanDepth.setKeyable(  true );    
    nAttrOceanDepth.setWritable(true);
    nAttrOceanDepth.setMin(1.f);
    nAttrOceanDepth.setSoftMin(1.f);    
    nAttrOceanDepth.setSoftMax(40000.f);        
    addAttribute( oceanDepth );
    attributeAffects( aaOceanDeformer::oceanDepth, aaOceanDeformer::outputGeom);

    MFnNumericAttribute nAttrSurfaceTension;
    surfaceTension= nAttrSurfaceTension.create( "surfaceTension", "surfaceTension", MFnNumericData::kFloat, 0.0f );
    nAttrSurfaceTension.setKeyable(  true );    
    nAttrSurfaceTension.setWritable(true);
    nAttrSurfaceTension.setSoftMin(0.f);    
    nAttrSurfaceTension.setSoftMax(5.f);        
    addAttribute( surfaceTension );
    attributeAffects( aaOceanDeformer::surfaceTension, aaOceanDeformer::outputGeom);

    MFnNumericAttribute nAttrSeed;
    seed = nAttrSeed.create( "Seed", "seed", MFnNumericData::kInt, 1 );
    nAttrSeed.setKeyable(  true );  
    nAttrSeed.setWritable(true);
    nAttrSeed.setSoftMin(1);    
    nAttrSeed.setSoftMax(32);   
    nAttrSeed.setMin(1);    
    nAttrSeed.setMax(64);   
    addAttribute( seed );
    attributeAffects( aaOceanDeformer::seed, aaOceanDeformer::outputGeom);

    MFnNumericAttribute nAttrWaveHeight;
    waveHeight = nAttrWaveHeight.create( "waveHeight", "waveHeight", MFnNumericData::kFloat, 2.0f );
    nAttrWaveHeight.setKeyable(  true );    
    nAttrWaveHeight.setWritable(true);
    nAttrWaveHeight.setSoftMin(0.001f); 
    nAttrWaveHeight.setSoftMax(100.f);  
    nAttrWaveHeight.setMin(0.001f);     
    addAttribute( waveHeight );
    attributeAffects( aaOceanDeformer::waveHeight, aaOceanDeformer::outputGeom);

    MFnNumericAttribute nAttrWaveSize;
    waveSize = nAttrWaveSize.create( "waveSize", "waveSize", MFnNumericData::kFloat, 4.0f );
    nAttrWaveSize.setKeyable(  true );  
    nAttrWaveSize.setWritable(true);
    nAttrWaveSize.setSoftMin(1.f);  
    nAttrWaveSize.setSoftMax(30.f);     
    addAttribute( waveSize );
    attributeAffects( aaOceanDeformer::waveSize, aaOceanDeformer::outputGeom);

    MFnNumericAttribute nAttrWaveSpeed;
    waveSpeed = nAttrWaveSpeed.create( "waveSpeed", "waveSpeed", MFnNumericData::kFloat, 1.0f );
    nAttrWaveSpeed.setKeyable(  true ); 
    nAttrWaveSpeed.setWritable(true);
    nAttrWaveSpeed.setSoftMin(1.f); 
    nAttrWaveSpeed.setSoftMax(10.f);    
    addAttribute( waveSpeed );
    attributeAffects( aaOceanDeformer::waveSpeed, aaOceanDeformer::outputGeom);

    MFnNumericAttribute nAttrWaveChop;
    waveChop = nAttrWaveChop.create( "waveChop", "waveChop", MFnNumericData::kFloat, 2.0f );
    nAttrWaveChop.setKeyable(  true );  
    nAttrWaveChop.setWritable(true);
    nAttrWaveChop.setSoftMin(0.0f); 
    nAttrWaveChop.setSoftMax(6.0f);
    addAttribute( waveChop );
    attributeAffects( aaOceanDeformer::waveChop, aaOceanDeformer::outputGeom);

    MFnNumericAttribute nAttrWaveSmooth;
    waveSmooth = nAttrWaveSmooth.create( "waveSmooth", "waveSmooth", MFnNumericData::kFloat, 0.0f);
    nAttrWaveSmooth.setKeyable(  true );    
    nAttrWaveSmooth.setWritable(true);
    nAttrWaveSmooth.setSoftMin(0.0f);   
    nAttrWaveSmooth.setSoftMax(20.f);
    addAttribute( waveSmooth );
    attributeAffects( aaOceanDeformer::waveSmooth, aaOceanDeformer::outputGeom);

    MFnNumericAttribute nAttrWaveDirection;
    waveDirection = nAttrWaveDirection.create( "waveDirection", "waveDirection", MFnNumericData::kFloat, 45.0f );
    nAttrWaveDirection.setKeyable(  true ); 
    nAttrWaveDirection.setWritable(true);
    nAttrWaveDirection.setSoftMin(0.0f);    
    nAttrWaveDirection.setSoftMax(360.0f);  
    nAttrWaveDirection.setMin(0.0f);    
    nAttrWaveDirection.setMax(360.0f);  
    addAttribute( waveDirection );
    attributeAffects( aaOceanDeformer::waveDirection, aaOceanDeformer::outputGeom);

    MFnNumericAttribute nAttrWaveReflection;
    waveReflection = nAttrWaveReflection.create( "waveReflection", "waveReflection", MFnNumericData::kFloat, 0.0f );
    nAttrWaveReflection.setKeyable(  true );    
    nAttrWaveReflection.setWritable(true);
    nAttrWaveReflection.setSoftMin(0.0f);   
    nAttrWaveReflection.setSoftMax(1.0f);   
    nAttrWaveReflection.setMin(0.0f);   
    nAttrWaveReflection.setMax(1.0f);   
    addAttribute( waveReflection );
    attributeAffects( aaOceanDeformer::waveReflection, aaOceanDeformer::outputGeom);

    MFnNumericAttribute nAttrWaveAlign;
    waveAlign = nAttrWaveAlign.create( "waveAlign", "waveAlign", MFnNumericData::kInt, 1 );
    nAttrWaveAlign.setKeyable(  true ); 
    nAttrWaveAlign.setWritable(true);
    nAttrWaveAlign.setSoftMin(0);   
    nAttrWaveAlign.setSoftMax(10);  
    nAttrWaveAlign.setMin(0);   
    nAttrWaveAlign.setMax(10);  
    addAttribute( waveAlign );
    attributeAffects( aaOceanDeformer::waveAlign, aaOceanDeformer::outputGeom);

    MFnUnitAttribute  uTime;
    time = uTime.create( "time", "time", MFnUnitAttribute::kTime);
    uTime.setHidden(true);
    uTime.setStorable(false);
    uTime.setKeyable(false);
    addAttribute( time );
    attributeAffects( aaOceanDeformer::time, aaOceanDeformer::outputGeom);

    MFnNumericAttribute nAttrRepeatTime;
    repeatTime = nAttrRepeatTime.create( "repeatTime", "repeatTime", MFnNumericData::kFloat, 1000.f );
    nAttrRepeatTime.setKeyable(  true );
    nAttrRepeatTime.setWritable(true);
    addAttribute( repeatTime );
    attributeAffects( aaOceanDeformer::repeatTime, aaOceanDeformer::outputGeom);
    
    MFnNumericAttribute nAttrTimeOffset;
    timeOffset = nAttrTimeOffset.create( "timeOffset", "timeOffset", MFnNumericData::kFloat, 0.f );
    nAttrTimeOffset.setKeyable(  true );
    nAttrTimeOffset.setWritable(true);
    addAttribute( timeOffset );
    attributeAffects( aaOceanDeformer::timeOffset, aaOceanDeformer::outputGeom);

    MFnNumericAttribute nAttrDoFoam;
    doFoam = nAttrDoFoam.create( "doFoam", "doFoam", MFnNumericData::kBoolean, 0 );
    nAttrDoFoam.setKeyable(  true );
    nAttrDoFoam.setWritable(true);
    addAttribute( doFoam );
    attributeAffects( aaOceanDeformer::doFoam, aaOceanDeformer::outputGeom);

    MFnNumericAttribute nAttrInvertFoam;
    invertFoam = nAttrInvertFoam.create( "invertFoam", "invertFoam", MFnNumericData::kBoolean, 0 );
    nAttrInvertFoam.setKeyable(  true );
    nAttrInvertFoam.setWritable(true);
    addAttribute( invertFoam );
    attributeAffects( aaOceanDeformer::invertFoam, aaOceanDeformer::outputGeom);

    MObject defaultValue;
    MFnTypedAttribute nAttrUVMap;
    defaultValue = MFnStringData().create(MString("map1"));
    uvMap = nAttrUVMap.create("uvMap", "uvMap",MFnData::kString, defaultValue);
    nAttrUVMap.setWritable(true);
    addAttribute(uvMap);
    attributeAffects( aaOceanDeformer::uvMap, aaOceanDeformer::outputGeom);

    MFnTypedAttribute nAttrEigenVectorMap;
    defaultValue = MFnStringData().create(MString("colorSetEigenVec"));
    eigenVectorMap = nAttrEigenVectorMap.create("eigenVectorMap", "eigenVectorMap",MFnData::kString, defaultValue);
    nAttrEigenVectorMap.setWritable(true);
    addAttribute(eigenVectorMap);
    attributeAffects( aaOceanDeformer::eigenVectorMap, aaOceanDeformer::outputGeom);

    MFnTypedAttribute nAttrEigenValueMap;
    defaultValue = MFnStringData().create(MString("colorSetEigenVal"));
    eigenValueMap = nAttrEigenValueMap.create("eigenValueMap", "eigenValueMap",MFnData::kString, defaultValue);
    nAttrEigenValueMap.setWritable(true);
    addAttribute(eigenValueMap);
    attributeAffects( aaOceanDeformer::eigenValueMap, aaOceanDeformer::outputGeom);

    MFnMatrixAttribute nAttrInTransform;
    inTransform = nAttrInTransform.create( "InputTransform", "InputTransform", MFnMatrixAttribute::kDouble);
    nAttrInTransform.setConnectable(true);
    nAttrInTransform.setHidden(true);
    addAttribute( inTransform );
    attributeAffects( aaOceanDeformer::inTransform, aaOceanDeformer::outputGeom);

    return MStatus::kSuccess;
}

aaOceanDeformer::aaOceanDeformer() 
{
    pOcean = new aaOcean;
    u.setLength(1);
    v.setLength(1);
    colArrayEigenVector.setLength(1);
    colArrayEigenValue.setLength(1);
    foundEigenVector = FALSE;
    foundEigenValue = FALSE;
	numPolygons = 0;
    MGlobal::displayInfo( "[aaOcean Maya] Created a new ocean patch" );
}
aaOceanDeformer::~aaOceanDeformer() 
{
    if(pOcean)
    {
        delete pOcean;
        pOcean = NULL;
        MGlobal::displayInfo( "[aaOcean Maya] Deleted ocean patch" );
    }
}

//void aaOceanDeformer::postConstructor()
//{
//    this->setDeformationDetails(MPxDeformerNode::kDeformsAll);
//}

void* aaOceanDeformer::creator()
{
    return new aaOceanDeformer();
}

// extern "C" __declspec(dllexport)
MStatus initializePlugin( MObject obj )
{
    MStatus result;
    MFnPlugin plugin( obj, "Amaan Akram", "2.6", "Any");
    result = plugin.registerNode( "aaOceanDeformer", aaOceanDeformer::id, aaOceanDeformer::creator, 
        aaOceanDeformer::initialize, MPxNode::kDeformerNode );
    return result;
}

// extern "C" __declspec(dllexport)
MStatus uninitializePlugin( MObject obj)
{
    MStatus result;
    MFnPlugin plugin( obj );
    result = plugin.deregisterNode( aaOceanDeformer::id );
    return result;
}
