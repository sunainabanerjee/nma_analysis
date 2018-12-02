# ifndef __TRAJECTORY_ANALYSIS_H
# define __TRAJECTORY_ANALYSIS_H

# include <set>
# include <cmath>
# include <cassert>
# include <algorithm>
# include "grid.h"
# include "octree.h"
# include "coordinate.h" 
# include "trajectory.h"
# include "residue_volume.h"
using namespace std;

typedef vector<GridSignature> TrajectorySignature;

enum Approximation {  _grid = 1, _octree = 2};


Grid3D sphere_approximation( Coordinate const& center,
                             const float radius,
                             const float precision )
{
    assert( (precision < radius) && (precision > 0) );
    float xmx = center.x + radius, ymx = center.y + radius, zmx = center.z + radius;
    float xmn = center.x - radius, ymn = center.y - radius, zmn = center.z - radius;
    Grid3D grid(xmx,ymx,zmx,xmn,ymn,zmn,precision, precision, precision);
    int n =  grid.cell_count();
    float v = grid.cell_volume();
    for( int i=0; i < n; ++i){
      Coordinate crd = grid.get_coordinate(i);
      if( euclidean_distance(center, crd) < (radius + 0.5 * precision) ){
          grid.fill_volume( crd, v);
      }
    }
    return grid;
}



void fill_sphere_in_grid( Grid3D& grid, 
                          Coordinate const& center, 
                          const float radius, 
                          const float precision,
                          const Approximation approx_type = _octree)
{
    if( approx_type == _grid )
    {
      Grid3D sphere = sphere_approximation( center, radius, precision );
      for(int j=0; j < sphere.cell_count(); ++j ){
         Coordinate crd = sphere.get_coordinate(j);
         grid.fill_volume(crd, sphere.cell_volume(j) );
      }
    }else{
       SphereIntersectionDecision decision(center, radius);
       
       Octree octree( Coordinate(center.x + radius, 
                                 center.y + radius, 
                                 center.z + radius),
                      Coordinate(center.x - radius, 
                                 center.y - radius, 
                                 center.z - radius),
                      decision);
       
       GridIndex gidx_start = grid.grid_coordinate( Coordinate(center.x - radius, 
                                                               center.y - radius, 
                                                               center.z - radius ) );
       
       GridIndex gidx_end = grid.grid_coordinate( Coordinate(center.x + radius, 
                                                             center.y + radius, 
                                                             center.z + radius ) );
       
       
       for( int i=gidx_start.xi; i <= gidx_end.xi; ++i)
          for( int j=gidx_start.yi; j <= gidx_end.yi; ++j )
             for( int k=gidx_start.zi; k <= gidx_end.zi; ++k )
             {
                AABB cell = grid.grid_cell_by_index( GridIndex(i,j,k) );
                float v = octree.intersecting_volume( cell );
                grid.fill_volume(cell.center(), v);
             }
    }
}


Grid3D fill_grid( Snapshot const& snapshot,   /* Snapshot containing all calpha positions */ 
                  Coordinate const& box_min,  /* 3D coordinate of lower corner of the bounding box */
                  Coordinate const& box_max,  /* 3D coordinate of the upper corner of the bounding box */
                  const float granularity,    /* size of each grid */
                  const float precision,      /* finer volume calculation smaller grid size */
                  const Approximation approx_type = _octree
                )
{
   assert( precision < granularity );
   assert( precision > 0. );
   int n = snapshot_size(snapshot);
   Grid3D grid( box_max.x, box_max.y, box_max.z, box_min.x, box_min.y, box_min.z, granularity, granularity, granularity );

   fill_sphere_in_grid( grid, 
                        fetch_coordinate(snapshot, snapshot.residueIds[0]), 
                        CARadius,
                        precision, 
                        approx_type );
   fill_sphere_in_grid( grid, 
                        fetch_coordinate(snapshot, snapshot.residueIds[n-1]), 
                        CARadius, 
                        precision, 
                        approx_type );

   for( int i=1; i < n-1; ++i ){
      int residue0 = snapshot.residueIds[i-1];
      int residue1 = snapshot.residueIds[i];
      int residue2 = snapshot.residueIds[i+1];

      Coordinate crd0 = fetch_coordinate(snapshot, residue0); 
      Coordinate crd1 = fetch_coordinate(snapshot, residue1);
      Coordinate crd2 = fetch_coordinate(snapshot, residue2);

      SideChainModel sm = fix_sidechain(crd0, crd1, crd2, snapshot.residueType[i]);

      fill_sphere_in_grid(grid, 
                          crd1, 
                          CARadius, 
                          precision, 
                          approx_type);
      fill_sphere_in_grid(grid, 
                          sm.center(), 
                          sm.radius(), 
                          precision, 
                          approx_type );
   }
   return grid;   
}


bool get_bounding_box( Trajectory const& trajectory, 
                       vector<int> const& residues, 
                       Coordinate& max_point, 
                       Coordinate& min_point )
{
   float mx, my, mz;  mx = my = mz = -9999.f;
   float nx, ny, nz;  nx = ny = nz =  9999.f;
   bool status = true;

   for( int i=0; i < trajectory_size(trajectory); ++i ){
     for( int j=0; j < residues.size(); ++j ){
        if( residue_position(trajectory[i], residues[j]) != -1 ){
           Coordinate crd = fetch_coordinate(trajectory[i], residues[j] );
           mx = (mx < crd.x)?crd.x:mx;
           my = (my < crd.y)?crd.y:my;
           mz = (mz < crd.z)?crd.z:mz;

           nx = (nx > crd.x)?crd.x:nx;
           ny = (ny > crd.y)?crd.y:ny;
           nz = (nz > crd.z)?crd.z:nz;
        }else{
           status = false;
        }
     }
   }
   max_point.x = mx; max_point.y = my; max_point.z = mz;
   min_point.x = nx; min_point.y = ny; min_point.z = ny;
   return status;
}


bool is_box( Coordinate const& mx_pt, Coordinate const& mn_pt ) {
   return (mx_pt.x > mn_pt.x) && (mx_pt.y > mn_pt.y) && (mx_pt.z > mn_pt.z);
}


float box_volume( Coordinate const& max_point, 
                  Coordinate const& min_point ) 
{
   assert( is_box(max_point, min_point) );
   return (max_point.x - min_point.x)*(max_point.y - min_point.y)*(max_point.z - min_point.z);
}


bool is_inside( Coordinate const& mx_point, 
                Coordinate const& mn_point, 
                Coordinate const& crd )
{
   if( ! is_box(mx_point, mn_point) ) return false;
   return ( mx_point.x > crd.x && crd.x > mn_point.x ) && 
          ( mx_point.y > crd.y && crd.y > mn_point.y ) && 
          ( mx_point.z > crd.z && crd.z > mn_point.z );
}


int box_intersect( Coordinate const& max_point1, 
                   Coordinate const& min_point1,
                   Coordinate const& max_point2, 
                   Coordinate const& min_point2 )
{
   Coordinate max_point_ref = max_point1, 
              min_point_ref = min_point1,
              max_point_tgt = max_point2,
              min_point_tgt = min_point2;

   if( box_volume(max_point1, min_point1) < box_volume(max_point2, min_point2) )
   {
      max_point_ref = max_point2;
      min_point_ref = min_point2;
      max_point_tgt = max_point1;
      min_point_tgt = min_point1;
   }
   
   int s = 0;
   if( is_inside( max_point_ref, min_point_ref, min_point_tgt ) ) s++;
   if( is_inside(max_point_ref, min_point_ref, Coordinate(min_point_tgt.x, min_point_tgt.y, max_point_tgt.z)) ) s++;
   if( is_inside(max_point_ref, min_point_ref, Coordinate(min_point_tgt.x, max_point_tgt.y, min_point_tgt.z)) ) s++;
   if( is_inside(max_point_ref, min_point_ref, Coordinate(min_point_tgt.x, max_point_tgt.y, max_point_tgt.z)) ) s++;
   if( is_inside(max_point_ref, min_point_ref, Coordinate(max_point_tgt.x, min_point_tgt.y, min_point_tgt.z)) ) s++;
   if( is_inside(max_point_ref, min_point_ref, Coordinate(max_point_tgt.x, min_point_tgt.y, max_point_tgt.z)) ) s++;
   if( is_inside(max_point_ref, min_point_ref, Coordinate(max_point_tgt.x, max_point_tgt.y, min_point_tgt.z)) ) s++;
   if( is_inside(max_point_ref, min_point_ref, max_point_tgt ) ) s++;
   
   return s;
}


bool is_intersecting( Coordinate const& max_point1, 
                      Coordinate const& min_point1, 
                      Coordinate const& max_point2, 
                      Coordinate const& min_point2 )
{
  return (box_intersect(max_point1, min_point1, max_point2, min_point2) > 0);
}


bool inclusive( Coordinate const& max_point1, 
                Coordinate const& min_point1,
                Coordinate const& max_point2, 
                Coordinate const& min_point2)
{
   return (box_intersect(max_point1, min_point1, max_point2, min_point2) == 8);
}


bool union_box( Coordinate const& max_point1, 
                Coordinate const& min_point1,
                Coordinate const& max_point2, 
                Coordinate const& min_point2,
                Coordinate& max_point,  
                Coordinate& min_point )
{
  assert( is_box(max_point1, min_point1) && is_box(max_point2, min_point2) );
  bool status = ( ! is_intersecting( max_point1, min_point1, max_point2, min_point2) );

  max_point.x = max( max_point1.x, max_point2.x );
  max_point.y = max( max_point1.y, max_point2.y );
  max_point.z = max( max_point1.z, max_point2.z );
  
  min_point.x = min( min_point1.x, min_point2.x );
  min_point.y = min( min_point1.y, min_point2.y );
  min_point.z = min( min_point1.z, min_point2.z );
  return status;
}


TrajectorySignature generate_volume_signature( Trajectory const& trajectory, 
                                                 Coordinate const& max_box,
                                                 Coordinate const& min_box, 
                                                 const float granularity, 
                                                 const float precision )
{
   TrajectorySignature vsig;
   for(int i=0; i < trajectory_size(trajectory); ++i )
   {
     Grid3D grid = fill_grid(trajectory[i], 
                             min_box, 
                             max_box, 
                             granularity, 
                             precision);
     vsig.push_back( grid.volume_signature() );
   }
   return vsig;
}


float grid_signature_difference( GridSignature const& sig1, 
                                 GridSignature const& sig2 )
{
   float score, s1, s2;
   score = s1 = s2 = 0.f;
   set<int> grid_index; 
   typename GridSignature::const_iterator it;
   for( it = sig1.begin(); it != sig1.end(); ++it )
        grid_index.insert(it->first);
   for( it = sig2.begin(); it != sig2.end(); ++it )
        grid_index.insert(it->first);
   for( typename set<int>::const_iterator i = grid_index.begin(); i != grid_index.end(); ++i )
   {
      s1 = ( (it = sig1.find(*i)) == sig1.end() )?0.f:it->second;
      s2 = ( (it = sig2.find(*i)) == sig2.end() )?0.f:it->second;
      score += fabs(s1 - s2);
   }
   return score;
}


vector<float> trajectory_signature_different( TrajectorySignature const& trj_sig1,  
                                              TrajectorySignature const& trj_sig2,
                                              const int theta = 0)
{
  assert( trj_sig1.size() == trj_sig2.size() );
  vector<float> v;
  int n = trj_sig1.size();
  for( int i=0; i < n; ++i )
     v.push_back( grid_signature_difference(trj_sig1[i], trj_sig2[ (i + theta) % n ]) );
  return v;
}


GridSignature occlusion_map( TrajectorySignature const& trj_sig ){
  set<int> grid_index;
  GridSignature omap;
  typename GridSignature::iterator iter;
  typename GridSignature::const_iterator const_iter;
  for( int i=0; i < trj_sig.size(); ++i ){
   for( const_iter = trj_sig[i].begin(); const_iter != trj_sig[i].end(); ++const_iter )
      if( (iter = omap.find(const_iter->first) ) == omap.end() )
          omap[const_iter->first] = const_iter->second;
      else if( iter->second < const_iter->second )
          iter->second = const_iter->second;
  }
  return omap;
}

# endif
