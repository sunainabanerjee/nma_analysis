# ifndef __GRID_H
# define __GRID_H 1

# include <map>
# include <cmath>
# include <vector>
# include <cassert>
# include <iomanip>
# include <iostream>
# include "geometry.h"
# include "coordinate.h"
using namespace std;

# ifndef EPS
# define EPS 1e-5
# endif

# ifndef INVALID_INDEX
# define INVALID_INDEX -1
# endif

struct GridIndex{
  int xi, yi, zi;
  public:
      GridIndex(const int x, 
                const int y, 
                const int z):xi(x),
                             yi(y),
                             zi(z) 
      {}
};

bool valid_grid_index( GridIndex const& gidx ){
   return (gidx.xi != INVALID_INDEX) && (gidx.yi != INVALID_INDEX) && (gidx.zi != INVALID_INDEX);
}

typedef float VolumeFractionOccupied;
typedef map<int,float> GridSignature;

class Grid3D{

  private:
   vector<VolumeFractionOccupied> grid_cells;
   float x_max, y_max, z_max;
   float x_min, y_min, z_min;
   float dx, dy, dz;
  public:
     Grid3D( const float xmx, 
             const float ymx, 
             const float zmx, 
             const float xmn, 
             const float ymn, 
             const float zmn,
             const float hx,
             const float hy,
             const float hz):x_max(xmx),y_max(ymx),z_max(zmx),x_min(xmn),y_min(ymn),z_min(zmn),dx(hx),dy(hy),dz(hz),grid_cells()
     {
        assert( (this->x_max > this->x_min) && (this->x_max - this->x_min > this->dx) );
        assert( (this->y_max > this->y_min) && (this->y_max - this->y_min > this->dy) );
        assert( (this->z_max > this->z_min) && (this->z_max - this->z_min > this->dz) );
        assert( this->dx > EPS && this->dy > EPS && this->dz > EPS );
        this->grid_cells = vector<VolumeFractionOccupied>(this->cell_count(), 0);
     }

     inline   int          nx()const { return static_cast<int>( ceil( (this->x_max - this->x_min)/ this->dx) ); }
     inline   int          ny()const { return static_cast<int>( ceil( (this->y_max - this->y_min)/ this->dy) ); }
     inline   int          nz()const { return static_cast<int>( ceil( (this->z_max - this->z_min)/ this->dz) ); }

     inline float grid_volume()const { return static_cast<float>( (this->x_max - this->x_min) * 
                                                                  (this->y_max - this->y_min) * 
                                                                  (this->z_max - this->z_min) ); }
     inline float cell_volume(void)const { return static_cast<float>(this->dx * this->dy * this->dz); }
 
     inline float cell_volume(const int i)const {
        if(i < 0 || i > this->cell_count()) return 0;
        return static_cast<float>( this->grid_cells[i] * this->cell_volume() );
     }
       
     inline float cell_volume(GridIndex const& gidx)const {
        return this->cell_volume( this->to_linear_index( gidx ) );
     }
    
     inline int   cell_count()const  { return static_cast<int>(this->nx() * this->ny() * this->nz()); }

     inline bool inside_grid( Coordinate const& crd )const { 
         return ( (crd.x < this->x_max) && (crd.y < this->y_max) && (crd.z < this->z_max) && 
                  (crd.x > this->y_min) && (crd.y > this->y_min) && (crd.z > this->z_min) );
     }     
   
     inline bool inside_grid( const float x, const float y, const float z)const {
        return this->inside_grid( Coordinate(x,y,z) );
     }

     inline GridIndex grid_coordinate( Coordinate const& crd )const {
         GridIndex gidx(INVALID_INDEX,INVALID_INDEX,INVALID_INDEX);
         if( this->inside_grid(crd) ){
            gidx.xi = static_cast<int>( (crd.x - this->x_min)/this->dx );
            gidx.yi = static_cast<int>( (crd.y - this->y_min)/this->dy );
            gidx.zi = static_cast<int>( (crd.z - this->z_min)/this->dz ); 
         }
         return gidx;
     }
     
     inline GridIndex grid_coordinate( const float x, const float y, const float z)const {
       return this->grid_coordinate( Coordinate(x,y,z) ); 
     }

     inline AABB grid_cell( Coordinate const& crd ) const {
        return  this->grid_cell_by_index( this->grid_coordinate(crd) );
     }
     
     inline AABB grid_cell_by_index( GridIndex const& gidx ) const {
        return AABB( Coordinate( (gidx.xi+1) * this->dx + this->x_min, 
                                 (gidx.yi+1) * this->dy + this->y_min, 
                                 (gidx.zi+1) * this->dz + this->z_min ),
                     Coordinate( gidx.xi * this->dx + this->x_min, 
                                 gidx.yi * this->dy + this->y_min, 
                                 gidx.zi * this->dz + this->z_min ) );
     }
     
     inline AABB grid_cell_by_index( const int i ) const {
          return this->grid_cell_by_index( this->to_grid_index(i) );
     }

     inline int to_linear_index( GridIndex const& gidx )const {
        int index = INVALID_INDEX;
        if( valid_grid_index(gidx) )
           index = gidx.xi * this->ny() * this->nz() + gidx.yi * this->nz() + gidx.zi;
        return index;
     }

     inline GridIndex to_grid_index( const int linear_index )const {
        GridIndex gidx(INVALID_INDEX, INVALID_INDEX, INVALID_INDEX);
        if( linear_index >= 0 && linear_index < this->cell_count() ){
           gidx.xi = static_cast<int>( linear_index / (this->nz() * this->ny() ) );
           gidx.yi = static_cast<int>( (linear_index % (this->ny() * this->nz())) / this->nz()  );
           gidx.zi = static_cast<int>( (linear_index %  (this->ny() * this->nz())) % this->nz() );
        }
        return gidx;
     }

     inline Coordinate get_coordinate( GridIndex const& gidx )const {
        assert( valid_grid_index(gidx) );
        Coordinate crd( (gidx.xi + 0.5) * this->dx + this->x_min, 
                        (gidx.yi + 0.5) * this->dy + this->y_min, 
                        (gidx.zi + 0.5) * this->dz + this->z_min);
        return crd;
     }

     inline Coordinate get_coordinate( const int linear_index )const{
        return this->get_coordinate( this->to_grid_index(linear_index) );
     }

     inline bool fill_volume( Coordinate const& crd, const float vol ){
         if( ! this->inside_grid( crd ) || vol < 0.f ) return false;
            
         int idx = this->to_linear_index( this->grid_coordinate( crd ) );
         float frac = vol/(1.f * this->cell_volume());
         
         frac = (frac > 1.0)?1.:(frac < 0.)?0.:frac;
         this->grid_cells[idx] += frac;
         if (this->grid_cells[idx] > 1.0)
            this->grid_cells[idx] = 1.0;
         return true;
     }

     inline bool free_volume( Coordinate const& crd, const float vol ){
        if( ! this->inside_grid(crd) || vol < 0.f ) return false;
        int idx = this->to_linear_index( this->grid_coordinate(crd) );
        float frac = vol / (1.f * this->cell_volume() );
        this->grid_cells[idx] = (frac > this->grid_cells[idx])?0.f:this->grid_cells[idx] - frac;
        return true;
     }

     inline float filled_volume()const{
         float f = 0;
         for(int i=0; i < this->grid_cells.size(); ++i )
            f += this->grid_cells[i];
         return this->cell_volume() *  f;
     }

     inline Coordinate min_coordinate()const {
        return Coordinate( this->x_min, this->y_min, this->z_min );
     }

     inline Coordinate max_coordinate()const {
        return Coordinate( this->x_max, this->y_max, this->z_max );
     }

     inline GridSignature volume_signature()const {
        GridSignature sig;
        for(int i=0; i < this->cell_count(); ++i )
           if( this->grid_cells[i] > 0 )
              sig[i] = this->grid_cells[i];
        return sig;
     }
};

ostream& operator << (ostream& os, GridSignature const& sig){
   typename GridSignature::const_iterator it;
   for( it = sig.begin(); it != sig.end(); ++it )
       os << it->first << ":" << fixed << setprecision(3) << it->second << endl;
   return os;
}


# endif /* __GRID_H */
