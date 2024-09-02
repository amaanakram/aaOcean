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

#include <xsi_matrix4f.h>

SICALLBACK aaOceanDataShader_aaOceanDataShader_1_0_DefineInfo( const CRef& in_ctxt )
{
    Context ctxt(in_ctxt);

    ctxt.PutAttribute("Category", "aaOceanDataShader" );
    ctxt.PutAttribute("DisplayName", "aaOceanDataShader" );
    
    return CStatus::OK;
}

SICALLBACK aaOceanDataShader_aaOceanDataShader_1_0_Define( const CRef& in_ctxt )
{
    Application app;
    Context ctxt(in_ctxt);

    ShaderDef sdef(ctxt.GetAttribute(L"Definition"));
    sdef.AddShaderFamily(siShaderFamilyTexture, true);

    // PARAMETERS
    Factory fact = app.GetFactory();
    ShaderParamDefContainer inpdefs = sdef.GetInputParamDefs();
    ShaderParamDefContainer outpdefs = sdef.GetOutputParamDefs();

    // Output
    CRef poutopts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions outopts = ShaderParamDefOptions(poutopts);
    outpdefs.AddParamDef(L"out", siShaderDataTypeColor4, poutopts);

    // Input
    CRef uv_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions uv_opts = ShaderParamDefOptions(uv_popts);
    uv_opts.SetTexturable(true);
    uv_opts.SetInspectable(false);
    uv_opts.SetDefaultValue(0);
    uv_opts.SetLongName("UV Coords");
    inpdefs.AddParamDef(L"uv_coords", siShaderDataTypeVector3, uv_popts);

    CRef use_uv_input_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions use_uv_input_opts = ShaderParamDefOptions(use_uv_input_popts);
    use_uv_input_opts.SetTexturable(false);
    use_uv_input_opts.SetInspectable(true);
    use_uv_input_opts.SetDefaultValue(false);
    use_uv_input_opts.SetLongName("Use Input UVs");
    inpdefs.AddParamDef(L"use_uv_input", siShaderDataTypeBoolean, use_uv_input_popts);

    CRef fade_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions fade_opts = ShaderParamDefOptions(fade_popts);
    fade_opts.SetTexturable(true);
    fade_opts.SetInspectable(true);
    fade_opts.SetHardLimit(0, 1);
    fade_opts.SetDefaultValue(0);
    fade_opts.SetLongName("Fade");
    inpdefs.AddParamDef(L"fade", siShaderDataTypeScalar, fade_popts);

    CRef resolution_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions resolution_opts = ShaderParamDefOptions(resolution_popts);
    resolution_opts.SetTexturable(true);
    resolution_opts.SetInspectable(true);
    resolution_opts.SetHardLimit(3, 10);
    resolution_opts.SetSoftLimit(4, 8);
    resolution_opts.SetDefaultValue(4);
    resolution_opts.SetLongName("Resolution");
    inpdefs.AddParamDef(L"resolution", siShaderDataTypeInteger, resolution_popts);

    CRef oceanScale_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions oceanScale_opts = ShaderParamDefOptions(oceanScale_popts);
    oceanScale_opts.SetTexturable(true);
    oceanScale_opts.SetInspectable(true);
    oceanScale_opts.SetDefaultValue(100);
    oceanScale_opts.SetHardLimit(0.0, 200000);
    oceanScale_opts.SetSoftLimit(0.0001,15.0);
    oceanScale_opts.SetLongName("Ocean Size");
    inpdefs.AddParamDef(L"oceanScale", siShaderDataTypeScalar, oceanScale_popts);

    CRef oceanDepth_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions oceanDepth_opts = ShaderParamDefOptions(oceanDepth_popts);
    oceanDepth_opts.SetTexturable(true);
    oceanDepth_opts.SetInspectable(true);
    oceanDepth_opts.SetDefaultValue(10000.f);
    oceanDepth_opts.SetHardLimit(1.0e-6f, 200000.f);
    oceanDepth_opts.SetSoftLimit(100, 10000.f);
    oceanDepth_opts.SetLongName("Ocean Depth");
    inpdefs.AddParamDef(L"oceanDepth", siShaderDataTypeScalar, oceanDepth_popts);

    CRef surfaceTension_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions surfaceTension_opts = ShaderParamDefOptions(surfaceTension_popts);
    surfaceTension_opts.SetTexturable(true);
    surfaceTension_opts.SetInspectable(true);
    surfaceTension_opts.SetDefaultValue(0.0f);
    surfaceTension_opts.SetHardLimit(0.0f, 10.0f);
    surfaceTension_opts.SetSoftLimit(0.0f, 1.f);
    surfaceTension_opts.SetLongName("Surface Tension");
    inpdefs.AddParamDef(L"surfaceTension", siShaderDataTypeScalar, surfaceTension_popts);

    CRef seed_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions seed_opts = ShaderParamDefOptions(seed_popts);
    seed_opts.SetTexturable(true);
    seed_opts.SetInspectable(true);
    seed_opts.SetHardLimit(1, 20);
    seed_opts.SetDefaultValue(1);
    seed_opts.SetLongName("Seed");
    inpdefs.AddParamDef(L"seed", siShaderDataTypeInteger, seed_popts);

    CRef waveHeight_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions waveHeight_opts = ShaderParamDefOptions(waveHeight_popts);
    waveHeight_opts.SetTexturable(true);
    waveHeight_opts.SetInspectable(true);
    waveHeight_opts.SetHardLimit(0.0, 10000);
    waveHeight_opts.SetSoftLimit(0.001,15.0);
    waveHeight_opts.SetDefaultValue(1);
    waveHeight_opts.SetLongName("Wave Height");
    inpdefs.AddParamDef(L"waveHeight", siShaderDataTypeScalar, waveHeight_popts);

    CRef velocity_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions velocity_opts = ShaderParamDefOptions(velocity_popts);
    velocity_opts.SetTexturable(true);
    velocity_opts.SetInspectable(true);
    velocity_opts.SetHardLimit(0.0, 30);
    velocity_opts.SetDefaultValue(5);
    velocity_opts.SetLongName("Wave Size");
    inpdefs.AddParamDef(L"velocity", siShaderDataTypeScalar, velocity_popts);

    CRef waveSpeed_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions waveSpeed_opts = ShaderParamDefOptions(waveSpeed_popts);
    waveSpeed_opts.SetTexturable(true);
    waveSpeed_opts.SetInspectable(true);
    waveSpeed_opts.SetHardLimit(0.00001, 1000);
    waveSpeed_opts.SetSoftLimit(0.001,15.0);
    waveSpeed_opts.SetDefaultValue(1);
    waveSpeed_opts.SetLongName("Wave Speed");
    inpdefs.AddParamDef(L"waveSpeed", siShaderDataTypeScalar, waveSpeed_popts);

    CRef chopAmount_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions chopAmount_opts = ShaderParamDefOptions(chopAmount_popts);
    chopAmount_opts.SetTexturable(true);
    chopAmount_opts.SetInspectable(true);
    chopAmount_opts.SetHardLimit(0.0, 100000);
    chopAmount_opts.SetDefaultValue(1);
    chopAmount_opts.SetSoftLimit(0.001,15.0);
    chopAmount_opts.SetLongName("Chop Amount");
    inpdefs.AddParamDef(L"chopAmount", siShaderDataTypeScalar, chopAmount_popts);

    CRef cutoff_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions cutoff_opts = ShaderParamDefOptions(cutoff_popts);
    cutoff_opts.SetTexturable(true);
    cutoff_opts.SetInspectable(true);
    cutoff_opts.SetHardLimit(0.0, 1000);
    cutoff_opts.SetSoftLimit(0.0,15.0);
    cutoff_opts.SetDefaultValue(0);
    cutoff_opts.SetLongName("Smooth Amount");
    inpdefs.AddParamDef(L"cutoff", siShaderDataTypeScalar, cutoff_popts);

    CRef windDir_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions windDir_opts = ShaderParamDefOptions(windDir_popts);
    windDir_opts.SetTexturable(true);
    windDir_opts.SetInspectable(true);
    windDir_opts.SetHardLimit(0.0, 360.0);
    windDir_opts.SetDefaultValue(45);
    windDir_opts.SetLongName("Wind Direction");
    inpdefs.AddParamDef(L"windDir", siShaderDataTypeScalar, windDir_popts);

    CRef damp_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions damp_opts = ShaderParamDefOptions(damp_popts);
    damp_opts.SetTexturable(true);
    damp_opts.SetInspectable(true);
    damp_opts.SetHardLimit(0.0, 1.0);
    damp_opts.SetDefaultValue(0.985);
    damp_opts.SetLongName("Reflected Waves");
    inpdefs.AddParamDef(L"damp", siShaderDataTypeScalar, damp_popts);

    CRef windAlign_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions windAlign_opts = ShaderParamDefOptions(windAlign_popts);
    windAlign_opts.SetTexturable(true);
    windAlign_opts.SetInspectable(true);
    windAlign_opts.SetHardLimit(0, 200);
    windAlign_opts.SetSoftLimit(0,10);
    windAlign_opts.SetDefaultValue(1);
    windAlign_opts.SetLongName("Wind Align");
    inpdefs.AddParamDef(L"windAlign", siShaderDataTypeInteger, windAlign_popts);

    CRef time_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions time_opts = ShaderParamDefOptions(time_popts);
    time_opts.SetTexturable(true);
    time_opts.SetInspectable(true);
    time_opts.SetDefaultValue(0.042);
    time_opts.SetLongName("Current Time");
    inpdefs.AddParamDef(L"time", siShaderDataTypeScalar, time_popts);

    CRef repeatTime_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions repeatTime_opts = ShaderParamDefOptions(repeatTime_popts);
    repeatTime_opts.SetTexturable(true);
    repeatTime_opts.SetInspectable(true);
    repeatTime_opts.SetDefaultValue(1000.f);
    repeatTime_opts.SetLongName("Repeat Time");
    inpdefs.AddParamDef(L"repeatTime", siShaderDataTypeScalar, repeatTime_popts);

    CRef gamma_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions gamma_opts = ShaderParamDefOptions(gamma_popts);
    gamma_opts.SetTexturable(true);
    gamma_opts.SetInspectable(true);
    gamma_opts.SetHardLimit(0, 1000);
    gamma_opts.SetSoftLimit(0, 10);
    gamma_opts.SetDefaultValue(1);
    gamma_opts.SetLongName("Gamma");
    inpdefs.AddParamDef(L"gamma", siShaderDataTypeScalar, gamma_popts);

    CRef brightness_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions brightness_opts = ShaderParamDefOptions(brightness_popts);
    brightness_opts.SetTexturable(true);
    brightness_opts.SetInspectable(true);
    brightness_opts.SetHardLimit(0, 1000);
    brightness_opts.SetSoftLimit(0, 10);
    brightness_opts.SetDefaultValue(1);
    brightness_opts.SetLongName("Brightness");
    inpdefs.AddParamDef(L"brightness", siShaderDataTypeScalar, brightness_popts);

    CRef rawOutput_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions rawOutput_opts = ShaderParamDefOptions(rawOutput_popts);
    rawOutput_opts.SetTexturable(true);
    rawOutput_opts.SetInspectable(true);
    rawOutput_opts.SetDefaultValue(false);
    rawOutput_opts.SetLongName("Raw Output");
    inpdefs.AddParamDef(L"rawOutput", siShaderDataTypeBoolean, rawOutput_popts);

    CRef invertFoam_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions invertFoam_opts = ShaderParamDefOptions(invertFoam_popts);
    invertFoam_opts.SetTexturable(true);
    invertFoam_opts.SetInspectable(true);
    invertFoam_opts.SetDefaultValue(false);
    invertFoam_opts.SetLongName("Invert Foam");
    inpdefs.AddParamDef(L"invertFoam", siShaderDataTypeBoolean, invertFoam_popts);

    CRef fmin_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions fmin_opts = ShaderParamDefOptions(fmin_popts);
    fmin_opts.SetTexturable(true);
    fmin_opts.SetInspectable(true);
    fmin_opts.SetDefaultValue(0.0);
    fmin_opts.SetLongName("Min");
    inpdefs.AddParamDef(L"fmin", siShaderDataTypeScalar, fmin_popts);

    CRef fmax_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions fmax_opts = ShaderParamDefOptions(fmax_popts);
    fmax_opts.SetTexturable(true);
    fmax_opts.SetInspectable(true);
    fmax_opts.SetDefaultValue(0.0f);
    fmax_opts.SetLongName("Max");
    inpdefs.AddParamDef(L"fmax", siShaderDataTypeScalar, fmax_popts);

    CRef layerOcean_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions layerOcean_opts = ShaderParamDefOptions(layerOcean_popts);
    layerOcean_opts.SetTexturable(true);
    layerOcean_opts.SetInspectable(false);
    layerOcean_opts.SetDefaultValue(0);
    layerOcean_opts.SetLongName("Layer Ocean");
    inpdefs.AddParamDef(L"layerOcean", siShaderDataTypeScalar, layerOcean_popts);

    // Input Parameter: Matrix
    XSI::MATH::CMatrix4f mat;
    mat = mat.SetIdentity();
    ShaderParamDefOptions paramOptions;
    paramOptions = fact.CreateShaderParamDefOptions();
    paramOptions.SetLongName("transform");
    paramOptions.SetAnimatable(true);
    paramOptions.SetTexturable(true);
    paramOptions.SetDefaultValue(CValue(mat));
    paramOptions.SetInspectable(true);
    inpdefs.AddParamDef( "transform", siShaderDataTypeMatrix44, paramOptions ); 

    CRef writeFile_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions writeFile_opts = ShaderParamDefOptions(writeFile_popts);
    writeFile_opts.SetTexturable(false);
    writeFile_opts.SetInspectable(true);
    writeFile_opts.SetDefaultValue(false);
    writeFile_opts.SetLongName("Write shader data");
    inpdefs.AddParamDef(L"writeFile", siShaderDataTypeBoolean, writeFile_popts);

    CRef outputFolder_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions outputFolder_opts = ShaderParamDefOptions(outputFolder_popts);
    outputFolder_opts.SetTexturable(false);
    outputFolder_opts.SetInspectable(true);
    outputFolder_opts.SetDefaultValue("");
    outputFolder_opts.SetLongName("Output Folder");
    outputFolder_opts.SetAttribute(siUIInitialDir, "c:\\");
    inpdefs.AddParamDef(L"outputFolder", siShaderDataTypeString, outputFolder_popts);

    CRef postfix_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions postfix_opts = ShaderParamDefOptions(outputFolder_popts);
    postfix_opts.SetTexturable(false);
    postfix_opts.SetInspectable(true);
    postfix_opts.SetDefaultValue("");
    postfix_opts.SetLongName("Postfix");
    inpdefs.AddParamDef(L"postfix", siShaderDataTypeString, postfix_popts);

    CRef currentFrame_popts = fact.CreateShaderParamDefOptions();
    ShaderParamDefOptions currentFrame_opts = ShaderParamDefOptions(currentFrame_popts);
    currentFrame_opts.SetTexturable(true);
    currentFrame_opts.SetInspectable(true);
    currentFrame_opts.SetDefaultValue(1);
    currentFrame_opts.SetLongName("Current Frame");
    inpdefs.AddParamDef(L"currentFrame", siShaderDataTypeInteger, currentFrame_popts);

    PPGLayout oPPGLayout = sdef.GetPPGLayout();

    oPPGLayout.AddTab("Ocean Parameters");
    oPPGLayout.AddGroup("Ocean Parameters");
        oPPGLayout.AddItem("resolution");
        oPPGLayout.AddItem("oceanScale");
        oPPGLayout.AddItem("oceanDepth");
        oPPGLayout.AddItem("surfaceTension");
        oPPGLayout.AddItem("seed");
        oPPGLayout.AddItem("time");
        oPPGLayout.AddItem("repeatTime");
        oPPGLayout.AddItem("use_uv_input");
    oPPGLayout.EndGroup();
    oPPGLayout.AddGroup("Wave Parameters");
        oPPGLayout.AddItem("waveHeight");
        oPPGLayout.AddItem("velocity");
        oPPGLayout.AddItem("waveSpeed");
        oPPGLayout.AddItem("chopAmount");
        oPPGLayout.AddItem("cutoff");
    oPPGLayout.EndGroup();
    oPPGLayout.AddGroup("Wind Parameters");
        oPPGLayout.AddItem("windDir");
        oPPGLayout.AddItem("damp");
        oPPGLayout.AddItem("windAlign");
    oPPGLayout.EndGroup();
    oPPGLayout.AddGroup("Misc");
        oPPGLayout.AddItem("fade");
        oPPGLayout.AddItem("transform");
    oPPGLayout.EndGroup();

    oPPGLayout.AddTab("Foam Parameters");
    oPPGLayout.AddItem("rawOutput");
    oPPGLayout.AddGroup("Non raw options");
    oPPGLayout.AddItem("invertFoam");
    oPPGLayout.AddItem("gamma");
    oPPGLayout.AddItem("brightness");
    oPPGLayout.AddItem("fmin");
    oPPGLayout.AddItem("fmax");
    oPPGLayout.EndGroup();

    oPPGLayout.AddTab("File Output");
    oPPGLayout.AddItem("writeFile");
    oPPGLayout.AddGroup("File Name");
    oPPGLayout.AddItem("outputFolder", "", siControlFolder );
    oPPGLayout.AddItem("postfix");
    oPPGLayout.AddItem("currentFrame");
    oPPGLayout.EndGroup();

    // RENDERERS
    // Renderer definition
    MetaShaderRendererDef rendererDef = sdef.AddRendererDef( "Mental Ray" );
    rendererDef.PutSymbolName( "aaOceanDataShader" );
    rendererDef.PutCodePath( "aaOceanDataShader" );
    return CStatus::OK;
}