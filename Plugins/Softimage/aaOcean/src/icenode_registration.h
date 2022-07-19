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

XSI::CStatus RegisteraaOcean( XSI::PluginRegistrar& in_reg );

SICALLBACK XSILoadPlugin( PluginRegistrar& in_reg )
{
    in_reg.PutAuthor(L"Amaan Akram");
    in_reg.PutName(L"aaOceanIceDeformer");
    in_reg.PutEmail(L"amaan@amaanakram.com");
    in_reg.PutURL(L"http://www.amaanakram.com");
    in_reg.PutVersion(2,6);

    RegisteraaOcean( in_reg );

    return CStatus::OK;
}

SICALLBACK XSIUnloadPlugin( const PluginRegistrar& in_reg )
{
    CString strPluginName;
    strPluginName = in_reg.GetName();
    Application().LogMessage(strPluginName + L" has been unloaded.");
    return CStatus::OK;
}