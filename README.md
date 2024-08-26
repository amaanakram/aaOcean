aaOcean is an implementation of Jerry Tessendorf's 2004 paper on Simulating Ocean Waves.
Author: Amaan Akram 
[www.linkedin.com/in/amaan](https://www.linkedin.com/in/amaan/)
First implemented in 2008 for Softimage XSI

**FEATURES**
* Novel FFT spectrum sampling ensures that the resulting ocean shape is
  only enhanced by increasing ocean resolution, and not completely changed
* Multi-threaded via OpenMP
* OpenEXR output for object-space vector displacement

aaOcean uses, sightings, references and derivative works
- [Ships with Renderman](https://rmanwiki-26.pixar.com/space/REN26/19661569/aaOceanPrmanShader)
- [Referenced in "Empirical directional wave spectra for computer graphics" by Christopher J. Horvath](https://dl.acm.org/doi/10.1145/2791261.2791267)
- Derivates made available by 3rd parties for Cinema 4D and Foundry's Modo

Example of work done with aaOcean:
https://vimeo.com/42087457

This repository contains the following

* aaOcean core class
* aaOcean Mental Ray shader
* Softimage Shader Definitions for Mental Ray shaders
* aaOcean Arnold shader
* aaOcean Prman 19+ RIS displacement shader
* aaOcean Softimage ICE deformer
* aaOcean Maya Deformer
* aaOcean Houdini SOP
* aaOcean standalong terminal/shell application
* several helper functions that I often use

Acknowledgements for help and bug fixes: Frederic Servant, Fabrice Macagno, Andrew Helmer,
The Softimage XSI Community, Renderman community

LICENSE: 
aaOcean is covered by a GNU GPL v3 license, unless another license is specifically 
granted by Amaan Akram.
A "New BSD" License for aaOcean can be obtained by contacting the author
For more details on aaOcean and associated 3rd Party licenses, please see
"license.txt" file 
