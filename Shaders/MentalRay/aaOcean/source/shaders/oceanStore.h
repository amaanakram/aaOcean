// aaOcean Mental Ray Shaders
// simple container class for convenience -- legacy, should be removed
// Author: Amaan Akram 
// www.amaanakram.com
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html

#ifndef OCEANSTORE_H
#define OCEANSTORE_H

class oceanStore
{
public:
    aaOcean*    ocean;
    miMatrix    transform;
    oceanStore()
    {
        ocean = 0;
        mi_matrix_ident(transform);
    }   
};

#endif /* OCEANSTORE_H */