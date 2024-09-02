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