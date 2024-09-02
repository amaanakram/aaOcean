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

// Defines port, group and map identifiers used for registering the ICENode
enum IDs
{
    ID_IN_PointID = 0,
    ID_IN_RESOLUTION,
    ID_IN_OCEAN_SCALE,
    ID_IN_OCEAN_DEPTH,
    ID_IN_SURFACE_TENSION,
    ID_IN_WINDDIR,
    ID_IN_CUTOFF,
    ID_IN_WINDVELOCITY,
    ID_IN_WINDALIGN,
    ID_IN_DAMP,
    ID_IN_WAVESPEED ,
    ID_IN_CHOP ,
    ID_IN_WAVE_HEIGHT ,
    ID_IN_SEED  ,
    ID_IN_ENABLE ,
    ID_IN_ENABLEFOAM ,
    ID_IN_U,
    ID_IN_V,
    ID_IN_TRANSFORM,
    ID_IN_TIME,
    ID_IN_REPEAT_TIME,
    ID_IN_RAND_WEIGHT ,

    ID_G_100 = 100,
    ID_G_101,
    ID_G_102,
    ID_G_103,
    ID_G_104,

    ID_G_200 = 200,
    ID_G_201,
    ID_G_202,
    ID_G_203,
    ID_G_204,

    ID_G_300 = 300,
    ID_OUT_OCEAN,
    ID_OUT_FOAM,
    ID_OUT_EIGEN_MINUS,
    ID_OUT_EIGEN_PLUS,
    ID_OUT_NORMALS,

    ID_UNDEF = ULONG_MAX
};

CStatus RegisteraaOcean( PluginRegistrar& in_reg )
{
    ICENodeDef nodeDef;
    nodeDef = Application().GetFactory().CreateICENodeDef(L"aaOcean");

    CStatus st;
    st = nodeDef.PutColor(133,192,192);
    st.AssertSucceeded();

    st = nodeDef.PutThreadingModel(XSI::siICENodeSingleThreading);
    st.AssertSucceeded();

    // Add input ports and groups.
    st = nodeDef.AddPortGroup(ID_G_100,1,1,L"Ocean Parameters");
    st.AssertSucceeded();
    st = nodeDef.AddPortGroup(ID_G_101,1,1,L"Wave Parameters");
    st.AssertSucceeded();
    st = nodeDef.AddPortGroup(ID_G_102,1,1,L"Wind Parameters"); 
    st.AssertSucceeded();
    st = nodeDef.AddPortGroup(ID_G_103,1,1,L"Data Inputs"); 
    st.AssertSucceeded();
    st = nodeDef.AddPortGroup(ID_G_104,1,1,L"Advanced"); 
    st.AssertSucceeded();

    st = nodeDef.AddInputPort(  ID_IN_ENABLE,
                                ID_G_100,
                                siICENodeDataBool,siICENodeStructureSingle,siICENodeContextSingleton,
                                L"Enable Ocean",L"Enable",
                                true);

    st = nodeDef.AddInputPort(  ID_IN_ENABLEFOAM,
                                ID_G_100,
                                siICENodeDataBool,siICENodeStructureSingle,siICENodeContextSingleton,
                                L"Enable Foam",L"Enable_Foam",
                                CValue(0));
    st.AssertSucceeded( ) ;

    st = nodeDef.AddInputPort(  ID_IN_RESOLUTION,
                                ID_G_100,
                                siICENodeDataLong,siICENodeStructureSingle,siICENodeContextSingleton,
                                L"Resolution",L"Resolution",
                                4,1,6,
                                ID_UNDEF, ID_UNDEF, ID_UNDEF);
    st.AssertSucceeded( ) ;

    st = nodeDef.AddInputPort(  ID_IN_OCEAN_SCALE,
                                ID_G_100,
                                siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
                                L"Ocean Size",L"Ocean_Size",
                                100.0f);
    st.AssertSucceeded( ) ;

    st = nodeDef.AddInputPort(  ID_IN_SEED,
                                ID_G_100,
                                siICENodeDataLong,siICENodeStructureSingle,siICENodeContextSingleton,
                                L"Seed",L"Seed",
                                1, 1, 20,
                                ID_UNDEF, ID_UNDEF, ID_UNDEF);
    st.AssertSucceeded( ) ;

    st = nodeDef.AddInputPort(  ID_IN_TIME,
                                ID_G_103,
                                siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
                                L"Current Time (secs)",L"CurrentTime",
                                1.0f);
    st.AssertSucceeded( ) ;

    st = nodeDef.AddInputPort(  ID_IN_WAVE_HEIGHT,
                                ID_G_101,
                                siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
                                L"Wave Height",L"Wave_Height",
                                1.0f);
    st.AssertSucceeded( ) ;

    st = nodeDef.AddInputPort(  ID_IN_WINDVELOCITY,
                                ID_G_101,
                                siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
                                L"Wave Size",L"Wave_Size",
                                5.0f, FLT_MIN, 30.f,
                                ID_UNDEF, ID_UNDEF, ID_UNDEF);
    st.AssertSucceeded( ) ;

    st = nodeDef.AddInputPort(  ID_IN_WAVESPEED,
                                ID_G_101,
                                siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
                                L"Wave Speed",L"Wave_Speed",
                                1.0f);
    st.AssertSucceeded( ) ;

    st = nodeDef.AddInputPort(  ID_IN_CHOP,
                                ID_G_101,
                                siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
                                L"Chopiness",L"Chopiness",
                                1.0f);
    st.AssertSucceeded( ) ;

    st = nodeDef.AddInputPort(  ID_IN_CUTOFF,
                                ID_G_101,
                                siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
                                L"Smooth",L"Smooth",
                                0.0f, 0.f, 100.f,
                                ID_UNDEF, ID_UNDEF, ID_UNDEF);
    st.AssertSucceeded( ) ;

    st = nodeDef.AddInputPort(  ID_IN_WINDDIR,
                                ID_G_102,
                                siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
                                L"Wind Direction",L"Wind_Direction",
                                45.0f, 0.f, 360.f,
                                ID_UNDEF, ID_UNDEF, ID_UNDEF);
    st.AssertSucceeded( ) ;

    st = nodeDef.AddInputPort(  ID_IN_DAMP,
                                ID_G_102,
                                siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
                                L"Reflected Waves",L"Reflected_Waves",
                                0.1f, 0.f, 1.f,
                                ID_UNDEF, ID_UNDEF, ID_UNDEF);
    st.AssertSucceeded( ) ;

    st = nodeDef.AddInputPort(  ID_IN_WINDALIGN,
                                ID_G_102,
                                siICENodeDataLong,siICENodeStructureSingle,siICENodeContextSingleton,
                                L"Wind Align",L"Wind_Align",
                                0);
    st.AssertSucceeded( ) ;

    st = nodeDef.AddInputPort(  ID_IN_REPEAT_TIME,
                                ID_G_102,
                                siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
                                L"Wave Loop Time (secs)",L"LoopTime",
                                1000.0f, 0.001f, 1000.f,
                                ID_UNDEF, ID_UNDEF, ID_UNDEF);
    st.AssertSucceeded( ) ;

    st = nodeDef.AddInputPort(  ID_IN_U,
                                ID_G_103,
                                siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextComponent0D,
                                L"Texture U",L"Texture_U",
                                1.0f);

    st = nodeDef.AddInputPort(  ID_IN_V,
                                ID_G_103,
                                siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextComponent0D,
                                L"Texture V",L"Texture_V",
                                1.0f);

    st = nodeDef.AddInputPort(  ID_IN_TRANSFORM,
                                ID_G_103,siICENodeDataMatrix44,siICENodeStructureSingle,siICENodeContextSingletonOrComponent0D,
                                L"Input_Transform",L"Input_Transform",
                                MATH::CMatrix4f(1.0,0.0,0.0,1.0, 
                                                0.0,1.0,0.0,1.0, 
                                                0.0,0.0,1.0,1.0, 
                                                0.0,0.0,0.0,1.0 ));

    st = nodeDef.AddInputPort(  ID_IN_PointID,
                                ID_G_103,
                                siICENodeDataLong,siICENodeStructureSingle,siICENodeContextComponent0D,
                                L"PointID",L"PointID",
                                1);
    st.AssertSucceeded( ) ;

    st = nodeDef.AddInputPort(  ID_IN_OCEAN_DEPTH,
                                ID_G_104,
                                siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
                                L"Ocean Depth",L"Ocean_Depth",
                                10000.0f);
    st.AssertSucceeded( ) ;

    st = nodeDef.AddInputPort(  ID_IN_SURFACE_TENSION,
                                ID_G_104,
                                siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
                                L"Surface Tension",L"Surface_Tension",
                                0.0f,0.0f,1.0f,
                                ID_UNDEF, ID_UNDEF, ID_UNDEF);
    st.AssertSucceeded( ) ;

    st = nodeDef.AddInputPort(  ID_IN_RAND_WEIGHT,
                                ID_G_104,
                                siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
                                L"Random Weight",L"Random Weight",
                                0.0f,0.0f,1.0f,
                                ID_UNDEF, ID_UNDEF, ID_UNDEF);
    st.AssertSucceeded( ) ;

//-----------------------------------START OUTPUT PORTS----------------

    st = nodeDef.AddPortGroup(ID_G_200);
    st.AssertSucceeded( ) ;
    st = nodeDef.AddOutputPort( ID_OUT_OCEAN,
                                ID_G_200,
                                siICENodeDataVector3,siICENodeStructureSingle,siICENodeContextComponent0D,
                                L"Ocean",L"Ocean");
    st.AssertSucceeded( ) ;

    st = nodeDef.AddPortGroup(ID_G_201);
    st.AssertSucceeded( ) ;
    st = nodeDef.AddOutputPort( ID_OUT_FOAM,
                                ID_G_201,
                                siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextComponent0D,
                                L"Foam",L"Foam");
    st.AssertSucceeded( ) ;

    st = nodeDef.AddPortGroup(ID_G_202);
    st.AssertSucceeded( ) ;
    st = nodeDef.AddOutputPort( ID_OUT_EIGEN_MINUS,
                                ID_G_202,
                                siICENodeDataVector3,siICENodeStructureSingle,siICENodeContextComponent0D,
                                L"egnVector Minus",L"egnVector_Minus");
    st.AssertSucceeded( ) ;

    st = nodeDef.AddPortGroup(ID_G_203);
    st.AssertSucceeded( ) ;

    st = nodeDef.AddOutputPort( ID_OUT_EIGEN_PLUS,
                                ID_G_203,
                                siICENodeDataVector3,siICENodeStructureSingle,siICENodeContextComponent0D,
                                L"egnVector Plus",L"egnVector_Plus");
    st.AssertSucceeded( ) ;

    /*st = nodeDef.AddPortGroup(ID_G_204);
    st.AssertSucceeded( ) ;

    st = nodeDef.AddOutputPort( ID_OUT_NORMALS,
                                ID_G_204,
                                siICENodeDataVector3,siICENodeStructureSingle,siICENodeContextComponent0D,
                                L"normals",L"normals");
    st.AssertSucceeded( ) ;*/

    PluginItem nodeItem = in_reg.RegisterICENode(nodeDef);
    nodeItem.PutCategories(L"aaOcean Suite");

    return CStatus::OK;
}
