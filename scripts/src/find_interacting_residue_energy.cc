#include <map>
#include <unistd.h>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <cstring>
#include <cassert>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <algorithm>
using namespace std;

#define MINIMUM     0L
#define MAXIMUM     1L
#define GRID_SPACE  3.f
#define BUFFER      1.f

class Atom;
class GridIndex;
typedef map<int, vector<Atom> > PDB;
typedef map<GridIndex, vector<Atom> > ForwardLookUp;
typedef map<Atom, GridIndex> ReverseLookUp;

static float  xRange[2], yRange[2], zRange[2];

string MAIN_CHAIN[] = { "N", "CA", "C", "O" };

char ELEMENTS[] = { 'H', 'C', 'N', 'O', 'S'};

float  VDW_RADIUS[] = {1.2f, 1.7f, 1.55f, 1.52f, 1.8f };

string  HB_DONOR[] = { "N", "NE", "NH1", "NH2", "ND2", "SG", "NE2", "ND1", "NZ", "OG", "OG1", "NE1", "OH"};

string  HB_ACCEPTOR[] = { "O", "OD1", "OD2", "OE1", "OE2", "ND1", "SD", "OG", "OG1", "OH"};

float   HB_DISTANCE[] = {2.5f,  3.2f, 4.5f };

float   HB_ENERGY[] = {115.f, 40.f, 17.f};

string  VDW_POLAR_ATOM[] = { "NE1", "OE1", "ND2", "OD1" };

float   VDW_SPACING[] = {0.5f};

float   VDW_ENERGY[] = {6.f};

string  POLAR_RESIDUE[] = { "ASP", "GLU", "ARG", "LYS", "HIS" };

string  POLAR_ATOM[] = {"OD1", "OD2", "OE1", "OE2", "NH1", "NH2", "NE", "CZ", "NZ", "ND1", "NE2"};

float   POLAR_DISTANCE[] = {4.5f};

float   POLAR_ENERGY[] = {20.f};

class GridIndex{
  public:
    int x, y, z;

    GridIndex( int _x = 0, int _y = 0, int _z = 0):x(_x),y(_y),z(_z) {}

    bool operator == (GridIndex const& g)const { return (this->x == g.x) && (this->y == g.y) && (this->z == g.z); }

    bool operator != (GridIndex const& g)const { return !(*this == g); } 

    bool operator  < (GridIndex const& g)const { 
        if( this->x == g.x ){
           if( this->y == g.y )
              return this->z < g.z;
           return this->y < g.y;
        }
        return this->x < g.x;
    }

    ~GridIndex() {}
};

class Atom{
 protected:
  string resName;
  string atomName;
  int    resId, atomId;
  float  x, y, z;

 public:
  Atom(string const& rn, 
       int ri, 
       string const& an,
       int ai,
       float xx,
       float yy,
       float zz):resName(rn),atomName(an),resId(ri),atomId(ai),x(xx),y(yy),z(zz)
  {}
  
  string get_resname()const  { return this->resName;}
  string get_atomname()const { return this->atomName;}
  int get_atomid()const { return this->atomId; }
  int get_resid()const  { return this->resId; }
  float get_x()const { return this->x; }
  float get_y()const { return this->y; }
  float get_z()const { return this->z; }

  string get_restag()const {
     char name[10];
     bzero(name, sizeof(name)*sizeof(char));
     sprintf(name, "%s%d", this->resName.c_str(), this->resId);
     return string(name);
  }

  bool operator == ( Atom const& atm )const {
    return (this->resId == atm.resId) && (this->atomId == atm.atomId);
  }

  bool operator < ( Atom const& atm )const {
    return this->atomId < atm.atomId;
  }

  virtual ~Atom(){}
};

ostream& operator << ( ostream& os, Atom const& atom ){
   os << atom.get_resname() << atom.get_resid() << ":" << atom.get_atomname();
   return os;
}

inline float atom_distance(Atom const& a1, Atom const& a2){
   float dx = a1.get_x() - a2.get_x();
   float dy = a1.get_y() - a2.get_y();
   float dz = a1.get_z() - a2.get_z();
   return sqrt( dx*dx + dy*dy + dz*dz );
}

bool is_main_chain( Atom const& atom ){
   int n = sizeof(MAIN_CHAIN)/sizeof(string);
   for( int i=0; i < n; ++i )
     if( MAIN_CHAIN[i] == atom.get_atomname() )
         return true;
   return false;
}

char get_atom_type(Atom const& atom ){
  return atom.get_atomname().at(0);
}

float get_vdw_radius( Atom const& atom ){
  char s = get_atom_type(atom);
  int  n = sizeof(ELEMENTS)/sizeof(char);
  for( int i=0; i < n; ++i )
    if( s == ELEMENTS[i] )
      return VDW_RADIUS[i];
  return 1.7f; 
}

bool is_polar( Atom const& atom ){
  int n1 = sizeof(POLAR_ATOM)/sizeof(string);
  int n2 = sizeof(POLAR_RESIDUE)/sizeof(string);
  bool yes = false;
  
  for(int i=0; i < n2; ++i)
    if( POLAR_RESIDUE[i] == atom.get_resname())
      yes = true;

  if(yes)
    for(int i=0; i < n1; ++i )
      if( POLAR_ATOM[i] == atom.get_atomname() )
        return true;
  return false;
}

bool is_colliding( Atom const& a, Atom const& b){
  float d1 = get_vdw_radius(a);
  float d2 = get_vdw_radius(b);
  return (atom_distance(a,b) < (d1 + d2));
}


bool is_hb_donor(Atom const& a){
  int n = sizeof(HB_DONOR)/sizeof(string);
  for( int i=0; i < n; ++i)
    if(a.get_atomname() == HB_DONOR[i] )
       return true;
  return false;
}

bool is_hb_acceptor(Atom const& a){
  int n = sizeof(HB_ACCEPTOR)/sizeof(string);
  for( int i=0; i < n; ++i )
    if( a.get_atomname() == HB_ACCEPTOR[i] )
       return true;
  return false;
}


bool is_hbond_contact(Atom const& a, Atom const& b){
  int n = sizeof(HB_DISTANCE)/sizeof(float);
  if( a.get_resid() != b.get_resid() && atom_distance(a,b) < HB_DISTANCE[n-1] ){
     int atype = (is_hb_donor(a))?-1:(is_hb_acceptor(a))?1:0;
     int btype = (is_hb_donor(b))?-1:(is_hb_acceptor(b))?1:0;
     if( (atype + btype) == 0 && atype*btype != 0)
         return true;
  }
  return false;
}

bool is_polar_contact(Atom const& a, Atom const& b){
  return  !is_colliding(a,b) && (a.get_resid() != b.get_resid()) && ! is_hbond_contact(a,b) && is_polar(a) && is_polar(b);
}

bool is_vdw_contact(Atom const& a, Atom const& b){
   return !is_colliding(a,b) && ! is_hbond_contact(a,b) && !is_polar_contact(a,b) && a.get_resid() != b.get_resid();
}

float get_interaction_energy( Atom const& a, Atom const& b){
  float d = atom_distance(a,b);
  if( is_hbond_contact(a,b) ){
     int n = sizeof(HB_DISTANCE)/sizeof(float);
     for( int i = 0; i < n; ++i )
           if( d < HB_DISTANCE[i] )
              return HB_ENERGY[i];
  }else if( is_polar_contact(a,b) ){
     int n = sizeof(POLAR_DISTANCE)/sizeof(float);
     for( int i = 0; i < n; ++i )
        if( d <= POLAR_DISTANCE[i] )
           return POLAR_ENERGY[i];
  }else if( is_vdw_contact(a,b) ){
     float s = d - (get_vdw_radius(a) + get_vdw_radius(b));
     int n = sizeof(VDW_SPACING)/sizeof(float);
     for(int i=0; i < n; ++i)
       if( s <= VDW_SPACING[i] )
          return VDW_ENERGY[i];
  }
  return 0;  
}

class AtomPair{
 public:
  Atom first;
  Atom second;
  AtomPair( Atom const& a1, Atom const& a2):first(a1),second(a2) {}
};

template<typename K, typename V>
bool has_key( map<K,V> const& lookup, K const& key){
  typename map<K,V>::const_iterator it = lookup.find(key);
  return it != lookup.end();
}

GridIndex get_grid_index( float x, float y, float z){
  assert( ( x > xRange[MINIMUM] ) && (x < xRange[MAXIMUM]) );
  assert( ( y > yRange[MINIMUM] ) && (y < yRange[MAXIMUM]) );
  assert( ( z > zRange[MINIMUM] ) && (z < zRange[MAXIMUM]) );
  return GridIndex( int((x - xRange[MINIMUM])/GRID_SPACE),  
                    int((y - yRange[MINIMUM])/GRID_SPACE),  
                    int((z - zRange[MINIMUM])/GRID_SPACE) ); 
}


bool build_grid( PDB const& pdb, ForwardLookUp& fl, ReverseLookUp& rl){
  fl.clear();
  rl.clear();
  for( typename PDB::const_iterator i = pdb.begin(); i != pdb.end(); ++i){
    int n = i->second.size();
    for( int j=0; j < n; ++j ){
      Atom const& atm = i->second[j];
      GridIndex gidx = get_grid_index(atm.get_x(), atm.get_y(), atm.get_z());
      if( ! has_key( fl, gidx ) )
            fl[gidx] = vector<Atom>();
      fl[gidx].push_back(atm);
      rl[atm] = gidx;
    }
  }
  return true;
}


PDB parse_pdb( string const& filename ){
  ifstream f( filename.c_str() , std::ifstream::in );
  string resName, atomName;
  int    resId, atmId, counter = 0;
  float  x, y, z; 
  PDB pdb;
  while( !f.eof()){
     f >> resName >> resId >> atomName >> atmId >> x >> y >> z;
     if( ! has_key(pdb, resId) )
         pdb[resId] = vector<Atom>();
     if( counter == 0 ){
          xRange[MINIMUM] = xRange[MAXIMUM] = x;
          yRange[MINIMUM] = yRange[MAXIMUM] = y;
          zRange[MINIMUM] = zRange[MAXIMUM] = z;
     }

     xRange[MINIMUM] = ( xRange[MINIMUM] > x)?x:xRange[MINIMUM];
     xRange[MAXIMUM] = ( xRange[MAXIMUM] < x)?x:xRange[MAXIMUM];
     yRange[MINIMUM] = ( yRange[MINIMUM] > y)?y:yRange[MINIMUM];
     yRange[MAXIMUM] = ( yRange[MAXIMUM] < y)?y:yRange[MAXIMUM];
     zRange[MINIMUM] = ( zRange[MINIMUM] > z)?z:zRange[MINIMUM];
     zRange[MAXIMUM] = ( zRange[MAXIMUM] < z)?z:zRange[MAXIMUM];

     pdb[resId].push_back( Atom(resName, resId, atomName, atmId, x, y, z) );
     counter++;
  }
  f.close();
  xRange[MINIMUM] -= BUFFER; xRange[MAXIMUM] += BUFFER;
  yRange[MINIMUM] -= BUFFER; yRange[MAXIMUM] += BUFFER;
  zRange[MINIMUM] -= BUFFER; zRange[MAXIMUM] += BUFFER;
  return pdb;
}


vector<Atom> find_neighboring_residue( Atom const& atom, const float d, ForwardLookUp const& fl, ReverseLookUp const& rl){
   int nGrid = int(ceil(d/GRID_SPACE));
   GridIndex gidx = rl.find(atom)->second;
   vector<Atom> nbrs;
   for( int x=-nGrid; x <= nGrid; ++x){
     for(int y=-nGrid; y <= nGrid; ++y){
        for(int z=-nGrid; z <= nGrid; ++z){
           GridIndex gn(gidx);
           gn.x += x;
           gn.y += y;
           gn.z += z;
           if( has_key(fl, gn) ){
               vector<Atom> nbr = fl.find(gn)->second;
               int n = nbr.size();
               for(int i=0; i < n; ++i )
                  if( nbr[i].get_resid() != atom.get_resid() && atom < nbr[i] && atom_distance( nbr[i], atom ) <= d ) 
                      nbrs.push_back(nbr[i]);
           }
        }
     }
   }
   return nbrs;
}

vector<AtomPair> get_all_atom_neighbors( PDB const& pdb, const float d ){
  ForwardLookUp fl;
  ReverseLookUp rl;
  assert(build_grid(pdb, fl, rl));
  vector<AtomPair> nbrs;
  for( typename PDB::const_iterator it = pdb.begin(); it != pdb.end(); ++it ){
     int n = it->second.size();
     for( int i=0; i < n; ++i ){
         vector<Atom> vlst = find_neighboring_residue( it->second[i], d , fl, rl );
         for( typename vector<Atom>::const_iterator j = vlst.begin(); j != vlst.end(); ++j ){
           if( get_atom_type(it->second[i]) == 'H' || get_atom_type(*j) == 'H' ) continue;
           nbrs.push_back( AtomPair(it->second[i], *j) );
         }
     }
  }
  return nbrs;
}

string interaction_type(Atom const& a, Atom const& b){
  string firstAtom = (is_main_chain(a))?"M":"S";
  string secondAtom = (is_main_chain(b))?"M":"S";
  return firstAtom + secondAtom;
}

bool is_polar_interaction(Atom const& a, Atom const& b){
  return (is_polar(a) && is_polar(b));
}


void usage(char const* filename ){
  fprintf( stderr, "Usage: %s -f <xyz file> -r <cutoff radius> -m\n" , filename );
  fprintf( stderr, "m: consider mainchain mainchain interaction\n");
  exit(0);
}

int main(int argc, char** argv){
  string filename = "";
  float distance = 4.5f;
  int c;
  char const* exename = argv[0];
  bool mainChain = false;

  if( argc < 3 ) usage(exename);

  while( (c = getopt(argc, argv, "f:r:mh")) != -1){
    switch(c){
      case 'f':
        filename = optarg; break;
      case 'r':
        distance = atof(optarg); break;
      case 'm':
        mainChain = true; break;
      default:
         usage(exename);
    }
  }
  if( filename == "" || distance == 0.f) 
     usage;

  PDB pdb = parse_pdb(filename);
  vector<AtomPair> nbrs = get_all_atom_neighbors(pdb, distance);
  map< string, map<string, float> > interaction_energy;
  vector<string> allowedContacts;
  allowedContacts.push_back("MS");
  allowedContacts.push_back("SM");
  allowedContacts.push_back("SS");
  if( mainChain ){
    allowedContacts.push_back("MM");
  }
  int nallow = allowedContacts.size();

  for( typename vector<AtomPair>::iterator it = nbrs.begin(); it != nbrs.end(); ++it){
      if(is_colliding(it->first, it->second) ) continue;
      string r1 = it->first.get_restag();
      string r2 = it->second.get_restag();
     
      float e = get_interaction_energy(it->first, it->second);
      string itype = interaction_type(it->first, it->second);
      bool allowed = false;
      for( int i=0; i  < nallow; ++i )
        if(allowedContacts[i] == itype )
            allowed = true;
      if( allowed && e > 0.f ){
         if( ! has_key( interaction_energy, r1 ) )
              interaction_energy[r1] = map<string, float>();        
         if( ! has_key( interaction_energy[r1], r2) )
              interaction_energy[r1][r2] = 0.f;
         interaction_energy[r1][r2] += e;
      }
  }
  for( typename map<string, map<string, float> >::const_iterator i = interaction_energy.begin(); i != interaction_energy.end(); ++i ){
     for( typename map<string, float>::const_iterator j = i->second.begin(); j != i->second.end(); ++j )
        cout << i->first << "," << j->first << "," << setprecision(2) << j->second <<endl;
  }
  return 0; 
}
