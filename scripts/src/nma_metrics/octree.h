# ifndef __OCTREE_H
# define __OCTREE_H

# include <vector>
# include <cassert>
# include "geometry.h"
# include "coordinate.h"
using namespace std;

# ifndef ALLOWED_SPLIT
# define ALLOWED_SPLIT  0.1
# endif


class SplitDecision{
public:
    virtual bool   operator () (AABB const&) = 0;
    virtual bool   enclosed( AABB const& ) = 0;
    virtual bool   intersecting(AABB const& ) = 0;
    virtual float  intersecting_volume(AABB const& ) = 0;
};


class ConvexSplitDecision:public SplitDecision
{
    public:
    virtual bool operator () ( AABB const& box ) 
    {
        return (box.dim() >= ALLOWED_SPLIT) && 
               (! this->enclosed(box)) && this->intersecting(box);
    }
    
    virtual bool inside_object( Coordinate const& ) const = 0;
    
    
    virtual bool enclosed( AABB const& box ){
        int n = this->corner_touched(box);
        return ( box.corners().size() == n);
    }
    
    virtual int corner_touched( AABB const& box )const {
        vector<Coordinate> c = box.corners();
        int enclosed = 0;
        for( int i=0; i < c.size(); ++i )
            if( this->inside_object(c[i]) )
                enclosed++;
        return enclosed;
    }
    
};

class SphereIntersectionDecision:public ConvexSplitDecision
{
private:
    Sphere3D sphere;
public:
    SphereIntersectionDecision( Sphere3D const& sph):sphere(sph)
    {
    }
    
    SphereIntersectionDecision( Coordinate const& center, const float r ):sphere(center, r)
    {
    }
    
    virtual bool inside_object( Coordinate const& crd ) const
    {
        return this->sphere.is_inside(crd);
    }
    
    virtual bool intersecting(AABB const& box)
    {  
        return is_intersecting(this->sphere, box);
    }
    
    virtual float intersecting_volume(AABB const& box) 
    {
        bool isect = this->intersecting(box);
        if( ! isect ) return 0.f;
        int n = this->corner_touched(box);
        int m = box.corners().size();
        float v = (n > 0 && isect)?((n * 1.0)/m) * box.volume():this->sphere.volume();
        return v;
    }
};


class Octree{
  private:
    AABB    box;
    Octree *children[8];
    float   volume_occupied;
  public:
    Octree( Coordinate const& maxPt, 
            Coordinate const& minPt,
            SplitDecision& decider):box(maxPt,minPt),volume_occupied(0.f)
    {
        for( int i=0; i < 8; i++ )
            this->children[i] = NULL;
            
        if( decider(this->box) )
        {
           vector<AABB> qlist = split_aabb_centrally(this->box);
           assert( qlist.size() == 8 );
           for(int i=0; i < 8; ++i )
           {
              this->children[i] = new Octree(qlist[i].max_corner(), 
                                             qlist[i].min_corner(), 
                                             decider ); 
              this->volume_occupied += this->children[i]->volume_occupied;
           }
        }else{
           this->volume_occupied += decider.intersecting_volume(this->box);
        }
        if(this->volume_occupied > this->box.volume() )
             this->volume_occupied = this->box.volume();
    }
    
    Octree( Octree const& octree ):box(octree.box.max_corner(), 
                                       octree.box.min_corner()),
                                   volume_occupied(octree.volume_occupied)
    {
        for( int i=0; i < 8; ++i){
            if( octree.children[i] == NULL )
                this->children[i] = NULL;
            else
                this->children[i] = new Octree(*(octree.children[i]));
        }
    }

    float volume()const 
    {
        return this->volume_occupied;
    }

    int nchilds()const
    {
        int c = 0;
        for(int i=0; i < 8; ++i)
            if( this->children[i] != NULL )
                c++;
        return c;
    }

    Octree const& operator [] ( const int i )const 
    {
        return *( this->children[i] );
    }
    
    float intersecting_volume(AABB const& cell)const
    {
        float f = 0.f;
        if( is_subset(this->box, cell) ){
            f = this->volume();
        }else if( is_intersecting(this->box, cell) ){
            if( this->nchilds() == 0){
                f = this->volume();
            }else{
                for(int i=0; i < 8; ++i )
                    f += this->children[i]->intersecting_volume(cell);
            }
        }
        return f;
    }
   
    virtual ~Octree() 
    {
        for( int i=0; i < 8; i++)
            if( this->children[i] != NULL )
                delete this->children[i];
        this->volume_occupied = 0.f;
    }
};

# endif

