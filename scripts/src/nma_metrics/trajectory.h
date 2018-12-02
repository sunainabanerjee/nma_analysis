# ifndef __TRAJECTORY_H
# define __TRAJECTORY_H

# include <cmath>
# include <vector>
# include <string>
# include <algorithm>
# include "coordinate.h"
using namespace std;

struct PDBxyz{
  vector<int>         residueIds;
  vector<string>      residueType;
  vector<Coordinate>  coordinates;
};

typedef struct PDBxyz Snapshot;
typedef vector<Snapshot> Trajectory;
typedef vector<int>  ResidueIds;

void       clear_snapshot( Snapshot&);
bool       valid_snapshot( Snapshot const& );
int        snapshot_size( Snapshot const& );
int        trajectory_size(Trajectory const& );

int        residue_position( Snapshot const&, const int );
Coordinate fetch_coordinate( Snapshot const&, const int );
string     fetch_resname(Snapshot const&, const int);

bool       snapshot_check( Snapshot const&, Snapshot const&, const bool = false );
Snapshot   snapshot_truncate( Snapshot const&, const int, const int);
float      snapshot_distance( Snapshot const&, Snapshot const&, const bool = false );
float      trajectory_distance( Trajectory const&, Trajectory const&, const bool = false, const int  = 0, const bool  = false );
int        trajectory_align( Trajectory const&, Trajectory const&, const bool = false , const bool = false);


int snapshot_size(Snapshot const& snapshot){
   return static_cast<int>(snapshot.residueType.size());
}


int trajectory_size(Trajectory const& trj ){
  return static_cast<int>(trj.size());
}


void clear_snapshot(Snapshot& snapshot ){
  snapshot.residueType.clear();
  snapshot.residueIds.clear();
  snapshot.coordinates.clear();
}


bool valid_snapshot( Snapshot& snapshot ){
   if( ( snapshot.residueType.size() == snapshot.residueIds.size() ) &&
       ( snapshot.residueType.size() == snapshot.coordinates.size() ) )
       return (snapshot_size(snapshot) > 0 );
   return false;
}


int residue_position( Snapshot const& snapshot, const int residueId ){
   typename vector<int>::const_iterator it;
   if( (it = find( snapshot.residueIds.begin(), snapshot.residueIds.end(), residueId)) == snapshot.residueIds.end() )
      return -1;
   return static_cast<int>( it - snapshot.residueIds.begin() );
}


Coordinate fetch_coordinate(Snapshot const& snapshot, const int residueId){
   int i = residue_position(snapshot, residueId );
   return ( i == -1 )?Coordinate():snapshot.coordinates[i];
}


string fetch_resname(Snapshot const& snapshot, const int residueId){
    int i = residue_position(snapshot, residueId);
    return (i == -1)?"":snapshot.residueType[i];
}


/**
 Check whether two given snapshot are comparable. It ensures the snapshot length,
 i.e., protein sizes are same, and residue Ids are comparable. 
**/
bool snapshot_check( Snapshot const& snap1, 
                     Snapshot const& snap2,
                     const bool allow_mutation )
{
  if( snapshot_size(snap1) != snapshot_size(snap2) ) return false;
  int n = snapshot_size(snap1);  
  if( allow_mutation )
    for( int i=0; i < n; ++i )
      if( (snap1.residueIds[i] != snap2.residueIds[i]) ) return false;
  else
    for( int i=0; i < n; ++i )
      if( (snap1.residueIds[i] != snap2.residueIds[i]) || (snap1.residueType[i] != snap2.residueType[i]) ) return false;
  return true;
}


Snapshot snapshot_truncate( Snapshot const& snap, const int start_resid, const int end_resid )
{
   int start, stop;
   assert( (start = residue_position(snap, start_resid)) != -1);
   assert( (stop  = residue_position(snap, end_resid)) != -1);
   assert( (start < stop) );
   Snapshot part;
   for(int i=start; i != stop; ++i ){
      part.residueIds.push_back( snap.residueIds[i] );
      part.residueType.push_back( snap.residueType[i] );
      part.coordinates.push_back( snap.coordinates[i] );
   }  
   return part;
}


float snapshot_distance( Snapshot const& snap1,
                         Snapshot const& snap2,
                         const bool allow_mutation )
{  
  assert( snapshot_check(snap1, snap2, allow_mutation) );
  int n = snapshot_size(snap1);
  float s, d;
  s = 0;
  # pragma omp parallel for private(d) reduction(+:s)
  for( int i=0; i < n; ++i ){
     d = euclidean_distance( snap1.coordinates[i], snap2.coordinates[i] );
     s += d*d; 
  }  
  return sqrt(s / n);
}


float trajectory_distance( Trajectory const& trj1, 
                           Trajectory const& trj2, 
                           const bool allow_mutation,
                           const int shift,
                           const bool flip )
{
  assert( trajectory_size(trj1) == trajectory_size(trj2) );
  int n = trajectory_size(trj1);
  float s,d;
  s = 0;
  if( flip ){
    # pragma omp parallel for private(d) reduction(+:s)
    for(int i=0; i < n; ++i ){
      d = snapshot_distance( trj1[i], trj2[(i+shift) % n], allow_mutation );
      s += d * d;
    }
  }else{
    # pragma omp parallel for private(d) reduction(+:s)
    for(int i=0; i < n; ++i ){
      d = snapshot_distance( trj1[i], trj2[ (n-i+shift) % n], allow_mutation );
      s += d * d;
    }
  }
  return sqrt( s / n);
}

int trajectory_align( Trajectory const& trj1, Trajectory const& trj2, const bool allow_mutation, const bool flip)
{
   assert( trajectory_size(trj1) == trajectory_size(trj2) );
   int n = trajectory_size(trj1);
   int shift = 0;
   float score = trajectory_distance(trj1, trj2, allow_mutation, shift);
   for( int i=1; i < n; ++i ){
     float s = trajectory_distance(trj1, trj2, allow_mutation, i, flip);
     if( s < score ){
          shift = i;
          score = s;
     }
   }
   return shift;
}

# endif
