# ifndef __RESIDUE_VOLUME_H
# define __RESIDUE_VOLUME_H

# include <cmath>
# include <string>
# include <cstring>
# include <cassert>
# include <cstdlib>
# include "geometry.h"
# include "coordinate.h"
using namespace std;

# ifndef EPS
# define EPS 1e-3
# endif

const string Aminos[] = {"ALA", "ARG", "ASN", "ASP", 
                         "CYS", "GLN", "GLU", "GLY", 
                         "HIS", "ILE", "LEU", "LYS", 
                         "MET", "PHE", "PRO", "SER", 
                         "THR", "TRP", "TYR", "VAL"};

const float  Radius[] = { 1.027, 3.794, 2.178, 1.865,
                          1.640, 2.626, 2.503, 0.485,
                          2.623, 2.283, 2.318, 3.477, 
                          2.885, 3.004, 1.605, 1.579, 
                          1.686, 3.491, 3.382, 1.685};


const float  CARadius = 1.7;
typedef Sphere3D SideChainModel;


float get_radius( string const& );
int naminos();
int amino_index( string const& );
SideChainModel fix_sidechain( Coordinate const&, Coordinate const&, Coordinate const&, string const&);


float get_radius( string const& resname ){
   int n = amino_index(resname);
   return (n != -1)?Radius[n]:0;
}


int naminos() { return static_cast<int>(sizeof(Aminos)/sizeof(string)); }

int amino_index( string const& resname ){
   int n = naminos();
   for(int i=0; i < n; ++i)
      if(Aminos[i] == resname )
         return i;
   return -1;
}

SideChainModel fix_sidechain( Coordinate const& ca_0, Coordinate const& ca_1, Coordinate const& ca_2 , string const& resname ){
   float r = get_radius(resname);
   Vector3D u(ca_0.x, ca_0.y, ca_0.z), v(ca_1.x, ca_1.y, ca_1.z), w(ca_2.x, ca_2.y, ca_2.z);
   Vector3D uv = (v - u).unit_vector();
   Vector3D wv = (v - w).unit_vector();
   Vector3D ov = (uv + wv).unit_vector();
   float x = ca_1.x + ov[0] * r;
   float y = ca_1.y + ov[1] * r;
   float z = ca_1.z + ov[2] * r;
   return SideChainModel(x,y,z,r); 
}



# endif
