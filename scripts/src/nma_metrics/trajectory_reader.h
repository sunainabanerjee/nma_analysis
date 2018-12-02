# ifndef __TRAJECTORY_READER_H
# define __TRAJECTORY_READER_H

# include <set>
# include <map>
# include <string>
# include <cstdio>
# include <vector>
# include <cstring>
# include <cstdlib>
# include <cassert>
# include <fstream>
# include "utility.h"
# include "trajectory.h"
# include "coordinate.h"
using namespace std;

# define PDB_LINEWIDTH 100


Trajectory read_trajectory( const string& );

bool to_append( vector<Snapshot>& trj, Snapshot& snapshot) {
   return ( valid_snapshot(snapshot) && (trj.size() == 0 || snapshot_size(trj[0]) == snapshot_size(snapshot)) );
}

Trajectory read_trajectory( const string& trjfile ){
   char line[PDB_LINEWIDTH + 1];
   Trajectory trajectory;
   bzero( (void*)line, PDB_LINEWIDTH + 1);
   assert( is_file(trjfile.c_str()) );
   ifstream f(trjfile.c_str());

   string residueType;
   int    resid;
   float  x,y,z;
   Snapshot snapshot;
   bool   append = false;

   while( f.getline(line, PDB_LINEWIDTH) ){
     if( strncmp(line, "MODEL", 5) == 0 ){
         if( to_append(trajectory, snapshot) ) 
           trajectory.push_back(snapshot);
         clear_snapshot(snapshot);
     }
     if( (strncmp(line, "ATOM  ",6 ) == 0) && (string(line).substr(13,2) == "CA" ) ){
        residueType = string(line).substr(17,3);
        resid = atoi(string(line).substr(22,4).c_str());
        x = atof(string(line).substr(30,8).c_str());
        y = atof(string(line).substr(38,8).c_str());
        z = atof(string(line).substr(46,8).c_str());
        snapshot.residueType.push_back(residueType);
        snapshot.residueIds.push_back(resid);
        snapshot.coordinates.push_back( Coordinate(x,y,z));
     }
     bzero( (void*)line, PDB_LINEWIDTH + 1);
   }
   if( snapshot_size(snapshot) > 0 && to_append(trajectory, snapshot) ){
      trajectory.push_back(snapshot);
      clear_snapshot(snapshot);
   }
   return trajectory;
}

# endif
