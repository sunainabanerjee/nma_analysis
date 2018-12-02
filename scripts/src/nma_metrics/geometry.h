# ifndef __GEOMETRY_H
# define __GEOMETRY_H

# include <cmath>
# include <vector>
# include <cassert>
# include <stdexcept>
# include "coordinate.h"
using namespace std;

# ifndef M_PI
# define M_PI 3.14159f
# endif

# define MINIMUM(x,y) (((x) < (y))?(x):(y))
# define MAXIMUM(x,y) (((x) > (y))?(x):(y))
# define SQ(x)        ((x)*(x))

class Sphere3D{
   private:
     Coordinate cxyz;
     float      R;
   public:
     Sphere3D(float x=0.f, float y=0.f, float z=0.f, float r=0.f):cxyz(x,y,z),R(r)
     {
       assert( this->R >= 0.f);
     }

     Sphere3D(Coordinate crd, float r=0.f):cxyz(crd),R(r)
     {
       assert( this->R >= 0.f);
     }

     Coordinate center()const {
       return this->cxyz;
     }
  
     float radius()const {
       return this->R;
     }

     bool is_inside(Coordinate const& c)const {
        return euclidean_distance(this->cxyz, c) < this->R;
     }

     float volume()const{
       return (4.f/3.f) * M_PI * this->R * this->R * this->R;
     }

     float surface_area()const{
        return 4.f * M_PI * this->R * this->R;
     }
};


class Vector3D{
  private:
     float xi;
     float yi;
     float zi;
  public:
      Vector3D( const float x = 0.f, const float y = 0.f, const float z = 0.f):xi(x),yi(y),zi(z)
      {}
     
      Vector3D( Coordinate const& c):xi(c.x),yi(c.y),zi(c.z)
      {}

      Vector3D( Coordinate const& from , Coordinate const& to):xi( to.x - from.x ), yi(to.y - from.y), zi(to.z - from.z)
      {}

      float norm()const{
         return sqrt( SQ(this->xi) + SQ(this->yi) + SQ(this->zi) ); 
      }
       
      Vector3D unit_vector()const{
         float n = this->norm();
         if( n > 0 )
           return Vector3D( this->xi/n, this->yi/n, this->zi/n );
         return *this;
      }

      size_t size()const { return 3; }

      Vector3D operator + ( Vector3D const& rhs )const {
         return Vector3D( this->xi + rhs.xi, this->yi + rhs.yi, this->zi + rhs.zi );
      }

      Vector3D operator - ( Vector3D const& rhs )const{
         return Vector3D( this->xi - rhs.xi, this->yi - rhs.yi, this->zi - rhs.zi );
      }

      float operator [] (const int index)const {
       switch( index ){
         case 0: return this->xi;
         case 1: return this->yi;
         case 2: return this->zi;
       }
       throw out_of_range("Error: Valid indices are (0,1,2)");
      }

      float& operator [] (const int index) {
       switch( index ){
         case 0: return this->xi;
         case 1: return this->yi;
         case 2: return this->zi;
       }
       throw out_of_range("Error: Valid indices are (0,1,2)");
      }

};

float dot( Vector3D const& v1, Vector3D const& v2){
  float s = 0;
  for( int i=0; i < (int)v1.size(); ++i)
    s = v1[i] * v2[i];
  return s;
}

Vector3D cross( Vector3D const& v1, Vector3D const& v2){
   return Vector3D( v1[1]*v2[2] - v1[2]*v2[1], v1[2]*v2[0] - v1[0]*v2[2], v1[0]*v2[1] - v1[1]*v2[0] );
}

Vector3D operator * ( const float x, Vector3D const& v){
   return Vector3D(x*v[0], x*v[1], x*v[2]);
}

Vector3D operator * (Vector3D const& v, const float x){
   return Vector3D(x*v[0], x*v[1], x*v[2]);
}



class Plane3D{
  private:
    Vector3D    direction;
    Coordinate  pt;   
  public:
    Plane3D( Coordinate const& c, Vector3D const& v):direction(v.unit_vector()),pt(c)
    {}

    Plane3D( float x, float y, float z, Vector3D const& v):direction(v.unit_vector()),pt(x,y,z)
    {}

    float closest_distance( Coordinate const& c)const{
      return abs((c[0] - pt[0]) * direction[0] + (c[1] - pt[1]) * direction[1] + (c[2] - pt[2]) * direction[2]);
    }
};

class AABB{
   private:
      Coordinate  min_coordinate;
      Coordinate  max_coordinate;
   public:
      AABB( const float xmx,
            const float ymx,
            const float zmx,
            const float xmn,
            const float ymn,
            const float zmn ):max_coordinate(xmx, ymx, zmx),min_coordinate(xmn, ymn, zmn)
      {
        assert( this->max_coordinate.x >= this->min_coordinate.x );
        assert( this->max_coordinate.y >= this->min_coordinate.y );
        assert( this->max_coordinate.z >= this->min_coordinate.z );
      }

      AABB( const Coordinate& max_corner,
            const Coordinate& min_corner ):max_coordinate(max_corner),min_coordinate(min_corner)
      {
        assert( this->max_coordinate.x > this->min_coordinate.x );
        assert( this->max_coordinate.y > this->min_coordinate.y );
        assert( this->max_coordinate.z > this->min_coordinate.z );
      }

      Coordinate center()const{
        return Coordinate( 0.5f * (min_coordinate.x + max_coordinate.x),
                           0.5f * (min_coordinate.y + max_coordinate.y),
                           0.5f * (min_coordinate.z + max_coordinate.z) );
      }

      float volume( )const{
        return (this->max_coordinate.x - this->min_coordinate.x) *
               (this->max_coordinate.y - this->min_coordinate.y) *
               (this->max_coordinate.z - this->min_coordinate.z);
      }

      float surface_area()const {
         float s = 0.f;
         s += 2 * (this->max_coordinate.x - this->min_coordinate.x) * (this->max_coordinate.y - this->min_coordinate.y);
         s += 2 * (this->max_coordinate.x - this->min_coordinate.x) * (this->max_coordinate.z - this->min_coordinate.z);
         s += 2 * (this->max_coordinate.z - this->min_coordinate.z) * (this->max_coordinate.y - this->min_coordinate.y);
         return s;
      }

      float min_x()const { return this->min_coordinate.x; }
      float min_y()const { return this->min_coordinate.y; }
      float min_z()const { return this->min_coordinate.z; }

      float max_x()const { return this->max_coordinate.x; }
      float max_y()const { return this->max_coordinate.y; }
      float max_z()const { return this->max_coordinate.z; }

      float dim_x()const { return this->max_coordinate.x - this->min_coordinate.x; }
      float dim_y()const { return this->max_coordinate.y - this->min_coordinate.y; }
      float dim_z()const { return this->max_coordinate.z - this->min_coordinate.z; }
      
      Coordinate max_corner()const {
          return this->max_coordinate;
      }
      
      Coordinate min_corner()const {
          return this->min_coordinate;
      }
      
      float dim()const   
      { 
          return (this->dim_x() < this->dim_y())?((this->dim_x() < this->dim_z())?this->dim_x():this->dim_z()):(this->dim_y() < this->dim_z())?this->dim_y():this->dim_z(); 
      }

      vector<Coordinate> corners()const{
         vector<Coordinate> pts;
         pts.push_back( Coordinate( this->min_coordinate.x , this->min_coordinate.y, this->min_coordinate.z) );
         pts.push_back( Coordinate( this->min_coordinate.x , this->min_coordinate.y, this->max_coordinate.z) );
         pts.push_back( Coordinate( this->min_coordinate.x , this->max_coordinate.y, this->min_coordinate.z) );
         pts.push_back( Coordinate( this->min_coordinate.x , this->max_coordinate.y, this->max_coordinate.z) );
         pts.push_back( Coordinate( this->max_coordinate.x , this->min_coordinate.y, this->min_coordinate.z) );
         pts.push_back( Coordinate( this->max_coordinate.x , this->min_coordinate.y, this->max_coordinate.z) );
         pts.push_back( Coordinate( this->max_coordinate.x , this->max_coordinate.y, this->min_coordinate.z) );
         pts.push_back( Coordinate( this->max_coordinate.x , this->max_coordinate.y, this->max_coordinate.z) );
         return pts;
      }

      bool is_inside(Coordinate const& c)const{
        return (c.x >= this->min_coordinate.x) && (c.x < this->max_coordinate.x) &&
               (c.y >= this->min_coordinate.y) && (c.y < this->max_coordinate.y) &&
               (c.z >= this->min_coordinate.z) && (c.z < this->max_coordinate.z);
      }
};

vector<AABB> split_aabb_centrally( AABB const& box ){
    Coordinate c = box.center();
    float mxx = box.max_x(), mxy = box.max_y(), mxz = box.max_z();
    float mnx = box.min_x(), mny = box.min_y(), mnz = box.min_z();
    vector<AABB> quadrants;
    quadrants.push_back( AABB( mxx, mxy, mxz, c.x, c.y, c.z ) );
    quadrants.push_back( AABB( c.x, mxy, mxz, mnx, c.y, c.z ) );
    quadrants.push_back( AABB( c.x, c.y, mxz, mnx, mny, c.z ) );
    quadrants.push_back( AABB( mxx, c.y, mxz, c.x, mny, c.z ) );
    quadrants.push_back( AABB( mxx, mxy, c.z, c.x, c.y, mnz ) );
    quadrants.push_back( AABB( c.x, mxy, c.z, mnx, c.y, mnz ) );
    quadrants.push_back( AABB( c.x, c.y, c.z, mnx, mny, mnz ) );
    quadrants.push_back( AABB( mxx, c.y, c.z, c.x, mny, mnz ) );
    return quadrants;
}


bool is_intersecting( Plane3D const& plane, Sphere3D const& sphere ){
  return plane.closest_distance( sphere.center() ) < sphere.radius();
}

bool is_intersecting( AABB const& box1, AABB const& box2){
  vector<Coordinate> vpt = box2.corners();
  for( typename vector<Coordinate>::const_iterator i = vpt.begin(); i != vpt.end(); ++i )
    if( box1.is_inside(*i) ) return true;
  return false;
}

bool is_subset( AABB const& query, AABB const& target )
{
    Coordinate qmx = query.max_corner();
    Coordinate qmn = query.min_corner();
    Coordinate tmx = target.max_corner();
    Coordinate tmn = target.min_corner();
    return ( (qmx.x <= tmx.x && qmx.x >= tmn.x) && (qmn.x <= tmx.x && qmn.x >= tmn.x) &&
             (qmx.y <= tmx.y && qmx.y >= tmn.y) && (qmn.y <= tmx.y && qmn.y >= tmn.y) &&
             (qmx.z <= tmx.z && qmx.z >= tmn.z) && (qmn.z <= tmx.z && qmn.z >= tmn.z) );
}

bool is_intersecting( AABB const& box, Sphere3D const& sph){
  Coordinate C = sph.center();
  Coordinate Bmax = box.max_corner();
  Coordinate Bmin = box.min_corner();
  float r2 = SQ(sph.radius());
  float dmin = 0.f;
  for( int i = 0; i < 3; i++ ) {
    if( C[i] < Bmin[i] ) dmin += SQ( C[i] - Bmin[i] );
    else if( C[i] > Bmax[i] ) dmin += SQ( C[i] - Bmax[i] );     
  }
  return dmin <= r2;
}

bool is_intersecting( Sphere3D const& sph, AABB const& box){
  return is_intersecting(box, sph);
}

bool is_intersecting( Sphere3D const& sph1, Sphere3D const& sph2 ){
   return ( euclidean_distance(sph1.center(), sph2.center()) < sph1.radius() + sph2.radius() );
}


float intersecting_volume( AABB const& box1, AABB const& box2 ){
   if( ! is_intersecting(box1, box2)) return 0.f;
   AABB ibox( Coordinate( MINIMUM(box1.max_x(), box2.max_x()), 
                           MINIMUM(box1.max_y(), box2.max_y()), 
                           MINIMUM(box1.max_z(), box2.max_z()) ),
               Coordinate( MAXIMUM(box1.min_x(), box2.min_x()), 
                           MAXIMUM(box1.min_y(), box2.min_y()), 
                           MAXIMUM(box1.min_z(), box2.min_z()) ) );
  return ibox.volume();
}

float intersecting_volume( Sphere3D const& sph1, Sphere3D const& sph2 ){
  if( ! is_intersecting(sph1, sph2))  return 0.f;
  float d = euclidean_distance(sph1.center(), sph2.center());
  float r = sph1.radius();
  float R = sph2.radius();
  return M_PI * SQ( R + r - d ) * (d*d + 2*d*r - 3*r*r + 2*d*R + 6 * r * R - 3*R*R)/(12.f * d);
}

AABB bounding_box( Sphere3D const& sphere ){
  Coordinate center = sphere.center();
  float r = sphere.radius();
  return AABB( center.x + r, center.y + r, center.z + r, center.x - r, center.y - r, center.z - r);
}


# endif
