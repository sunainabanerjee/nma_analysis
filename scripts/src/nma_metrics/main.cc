# include <iostream>
# include "octree.h"
# include "geometry.h"
using namespace std;


int main(int argc, char** argv)
{
    Sphere3D sph(0.f, 0.f, 0.f, 1.f);
    Coordinate max_coordinate(  4.f, 4.f, 4.f );
    Coordinate min_coordinate( -4.f,-4.f,-4.f );
    AABB box( Coordinate(1,1,1), Coordinate(0,0,0));
    SphereIntersectionDecision decision(sph);
    Octree tree(max_coordinate, min_coordinate, decision);
    cout << "Volume estimated: " << tree.volume() << " Actual volume: " << sph.volume() << endl;
    cout << "Intersecting volume: " << tree.intersecting_volume(box) << endl;
    return 0;
}