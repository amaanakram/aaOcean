aaOcean is an implementation of Jerry Tessendorf's 2004 paper on Simulating Ocean Waves.
Author: Amaan Akram 
www.amaanakram.com

**FEATURES**
* Multi-threaded via OpenMP
* FFT spectrum sampling ensures that the resulting ocean shape is only enhanced by increasing ocean resolution, 
and not completely changed
* OpenEXR output for object-space vector displacement

Example of work done with aaOcean:
https://vimeo.com/42087457

This repository contains the following

* aaOcean core class
* aaOcean Mental Ray shader
* Softimage Shader Definitions for Mental Ray shaders
* aaOcean Arnold shader
* aaOcean Prman 19 RIS displacement shader
* aaOcean Softimage ICE deformer
* aaOcean Maya Deformer
* aaOcean Houdini SOP
* aaOcean standalong terminal/shell application
* several helper functions that I often use

*****LINUX BUILD INSTRUCTIONS*******
Please see the wiki page at 
https://bitbucket.org/amaanakram/aaocean/wiki/Build%20instructions

MAYA & LINUX:
To compile the Maya plugin, you will need gcc-4.2.4 which is ABI-compatible with gcc-4.1.x and
supports OpenMP which is used by aaOcean

Acknowledgements for help and bug fixes: Frederic Servant, Fabrice Macagno, Phil Stopford, Andrew Helmer,
The Softimage XSI Community.

LICENSE: 
aaOcean is covered by a GNU GPL v3 license, unless another license is specifically 
granted by Amaan Akram.
A "New BSD" License for aaOcean can be obtained by contacting the author
For more details on aaOcean and associated 3rd Party licenses, please see
"license.txt" file that is part of the aaOcean repository:
https://bitbucket.org/amaanakram/aaocean# aaOcean