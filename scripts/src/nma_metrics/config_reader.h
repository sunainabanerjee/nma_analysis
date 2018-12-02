# ifndef __CONFIG_READER_H
# define __CONFIG_READER_H

# include <map>
# include <cstdlib>
# include <cstring>
# include <cassert>
# include "utility.h"
# include "coordinate.h"
# include "trajectory.h"
# include "trajectory_reader.h"
using namespace std;

# define CONFIG_LINEWIDTH 200

struct BoundingBoxConfig{
  map<int, Trajectory> trajectory;
  ResidueIds    bounding_residues;
};

struct ModeCompareConfig{
   Coordinate max_point, min_point;
   Trajectory trajectory_ref, trajectory_tgt;
};

BoundingBoxConfig   read_boundingbox_configfile( char const* );
ModeCompareConfig   read_modecompare_configfile( char const* );

BoundingBoxConfig read_boundingbox_configfile( char const* filename ){
  ifstream f(filename);
  assert( f.good() );
  map<int,string>  trjfiles;
  set<int>         residues;
  char             line[CONFIG_LINEWIDTH + 1];
  bzero( line, CONFIG_LINEWIDTH + 1);
  while( f.getline(line, CONFIG_LINEWIDTH) ){
      vector<string> words = string_split(string(line),":");
      if( words.size() > 1 ){
         if(words[0] == "M" ){
           assert( words.size() == 3 );
           assert( is_integer(words[1]) );
           assert( is_file(words[2]) );
           trjfiles[atoi(words[1].c_str())] = words[2];
         }else if( words[0] == "R" ){
           assert( words.size() == 2 );
           assert( is_integer(words[1]) );
           residues.insert(atoi(words[1].c_str()));
         }
      }
      bzero( line, CONFIG_LINEWIDTH);
  }
  assert( residues.size() > 3 );
  BoundingBoxConfig config;
  for( typename map<int,string>::const_iterator i = trjfiles.begin(); i != trjfiles.end(); ++i){
    config.trajectory[i->first] = read_trajectory( i->second );
    assert( trajectory_size(config.trajectory[i->first]) > 0) ;
  }
  config.bounding_residues = ResidueIds(residues.begin(), residues.end());

  for(int i=0; i < config.bounding_residues.size(); ++i){
     for( typename map<int,Trajectory>::const_iterator j=config.trajectory.begin(); j != config.trajectory.end(); ++j ){
             assert( find( j->second[0].residueIds.begin(),
                           j->second[0].residueIds.end(),
                           config.bounding_residues[i]) != j->second[0].residueIds.end() );
     }
  }
  return config;
}


ModeCompareConfig   read_modecompare_configfile( char const* filename ){
  ifstream f(filename);
  assert( f.good() );
  map<int,string>  trjfiles;
  set<int>         residues;
  char             line[CONFIG_LINEWIDTH + 1];
  bzero( line, CONFIG_LINEWIDTH + 1);
  ModeCompareConfig config;
  while( f.getline(line, CONFIG_LINEWIDTH) ){
      vector<string> words = string_split(string(line),":");
      if( words.size() > 1 ){
         if(words[0] == "S" ){
           assert( words.size() == 2 );
           vector<string> flds = string_split(words[1], ",");
           assert(flds.size() == 3);
           config.min_point = Coordinate(atof(flds[0].c_str()),  atof(flds[1].c_str()), atof(flds[2].c_str()));
         }else if( words[0] == "E" ){
           assert( words.size() == 2 );
           vector<string> flds = string_split(words[1], ",");
           assert(flds.size() == 3);
           config.max_point = Coordinate(atof(flds[0].c_str()),  atof(flds[1].c_str()), atof(flds[2].c_str()));
         }else if( words[0] == "R"){
           assert( words.size() == 2 );
           is_file( words[1] );
           config.trajectory_ref = read_trajectory( words[1] );
         }else if( words[0] == "C"){
           assert( words.size() == 2 );
           is_file( words[1] );
           config.trajectory_tgt = read_trajectory( words[1] );
         }
      }
      bzero( line, CONFIG_LINEWIDTH);
  }
  return config;
}


# endif
