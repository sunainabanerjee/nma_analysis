# include <iostream>
# include <fstream>
# include <string>
# include <cstdlib>
# include <cstdio>
# include <cmath>
# include <map>
# include <set>
using namespace std;

# define EPS        0.01
# define SQ(x)      ((x)*(x))


bool is_file( const char* filename ){
  ifstream f(filename);
  return f.good();
}

template<typename T>
bool compare_equality(T const& f, T const& s){
   return f == s;
}


template< typename T>
class Pair{
public:
  T first;
  T second;
public:
  Pair(T const& m , T const& r):first(m), second(r) {}
  Pair(Pair<T> const& p):first(p.first), second(p.second) {}

  bool operator == (Pair<T> const& p)const {
    return compare_equality(p.first , this->first) && compare_equality(p.second, this->second);
  }
 
  bool operator < ( Pair<T> const& p)const {
    if( compare_equality(p.first, this->first) ) return this->second < p.second;
    return this->first < p.first;
  }

};

typedef Pair<int>   iPair;
typedef Pair<float> fPair;

int main( int argc, char** argv ){

   if( argc != 2 || ! is_file( argv[1]) ){
     printf("Usage: %s  <coordinate-model file>\n", argv[0] );
     exit( 1 );
   }
 
   map<iPair, float>  x_map, y_map, z_map;
   set<int> models, residues;
   ifstream f(argv[1]);
   if( f.is_open() ){
      int modelId, residueId;
      float x, y, z;
      char sep;
      while( f >> modelId >> sep >> residueId >> sep >> x >> sep >> y >> sep >> z ){
        models.insert(modelId);
        residues.insert(residueId);
        iPair p(modelId, residueId);
        x_map[p] = x;
        y_map[p] = y;
        z_map[p] = z;
      }
      f.close();
   }
   
   for( set<int>::const_iterator ri = residues.begin(); ri != residues.end(); ri++){
     for( set<int>::const_iterator rj = residues.begin(); rj != residues.end(); rj++ ){
        if( *ri < *rj ){
          int n = 0;
          float x = 0, x2 = 0, xm = 0;
          iPair pi(1,*ri), pj(1,*rj);
          float d0 = sqrt( SQ(x_map[pi] - x_map[pj]) + SQ(y_map[pi] - y_map[pj]) + SQ(z_map[pi] - z_map[pj]) ); 
          for( set<int>::const_iterator m = models.begin(); m != models.end(); m++){
             if( *m > 1 ){
                iPair pi(*m,*ri), pj(*m,*rj);
                float d = sqrt( SQ(x_map[pi] - x_map[pj]) + SQ(y_map[pi] - y_map[pj]) + SQ(z_map[pi] - z_map[pj]) );
                x += (d - d0);
                x2 += (d - d0) * (d - d0);
                xm += abs(d - d0);
                n++;
             }
          }
          float stdDev = (x2/(n*1.0)) - SQ(x/(n*1.0));
          float umean = (xm/(n*1.0));
          printf( "%d,%d,%.3f,%.4f,%.5f\n", *ri,*rj,d0,umean,stdDev ); 
        }
     }
   }

   return 0;
}
