# include <map>
# include <cassert>
# include <iomanip>
# include <iostream>
# include "trajectory_reader.h"
# include "trajectory_analysis.h"
using namespace std;

int main(int argc, char** argv){
   if( argc != 2 || ! is_file(argv[1]) ){
      cerr << "Usage: " << argv[0] << " config.file " << endl;
      exit(1);
   }
   Config config = read_configfile( argv[1]);
   cout << "Number of snapshots: " << trajectory_size(config.trajectory) << endl;
   cout << "Number of residues: " << config.bounding_residues.size() << endl;
  
   assert( trajectory_size(config.trajectory) > 0 );

   Coordinate max_point, min_point;

   get_bounding_box( config.trajectory, config.bounding_residues, max_point, min_point );
   
   TrajectorySignature trj_sig = generate_volume_signature( config.trajectory, max_point, min_point, 1.5, .25);
 
   vector<float> trj_diff = trajectory_signature_different( trj_sig, trj_sig , 12 ); 
  
   for(int i=0; i < trj_diff.size(); ++i) cout << trj_diff[i] << endl; 

   return 0;
}
