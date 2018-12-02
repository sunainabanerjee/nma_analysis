# include <map>
# include <cassert>
# include <iomanip>
# include <iostream>
# include "trajectory_reader.h"
# include "trajectory_analysis.h"
# include "config_reader.h"
using namespace std;

/**
Input Config file must contain list 
**/

int main(int argc, char** argv){
   if( argc != 2 || ! is_file(argv[1]) ){
      cerr << "Usage: " << argv[0] << " config.file " << endl;
      exit(1);
   }
   BoundingBoxConfig config = read_boundingbox_configfile( argv[1]);
   /*
   cout << "Number of modes: " << config.trajectory.size() << endl;
   for( typename map<int,Trajectory>::const_iterator it = config.trajectory.begin(); it != config.trajectory.end(); ++it ){
      cout << "Number of snapshots in mode -(" <<  it->first <<"): " << trajectory_size(it->second) << endl;
   }
   cout << "Number of residues: " << config.bounding_residues.size() << endl;
   */
  
   vector<Coordinate> max_point, min_point; 
   Coordinate maxAux, minAux, mxPt, mnPt;
   for( typename map<int,Trajectory>::const_iterator it = config.trajectory.begin(); it != config.trajectory.end(); ++it ){
      get_bounding_box( it->second, config.bounding_residues, mxPt, mnPt );
      max_point.push_back(mxPt);
      min_point.push_back(mnPt);
   }
   maxAux = max_point[0];
   minAux = min_point[0];
   for( int i=1; i < (int)max_point.size(); ++i){
      union_box(maxAux, minAux, max_point[i], min_point[i], mxPt, mnPt);
      maxAux = mxPt;
      minAux = mnPt;
   }
   cout << mxPt << endl;
   cout << mnPt << endl;    
   return 0;
}
