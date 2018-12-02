# include <vector>
# include <cassert>
# include "utility.h"
# include "trajectory_analysis.h"
# include "config_reader.h"
using namespace std;

# define GRID_SIZE     1.5f
# define VOL_DISCRETE  0.3f

int  main( int argc, char** argv){

   if( argc != 2 || ! is_file(argv[1]) ){
      cerr << "Usage: " << argv[0] << " config.file " << endl;
      exit(1);
   }
   ModeCompareConfig config = read_modecompare_configfile( argv[1]);
   vector<float> vsig = trajectory_signature_different(
                                   generate_volume_signature(config.trajectory_ref , 
                                                             config.max_point, 
                                                             config.min_point, 
                                                             GRID_SIZE, 
                                                             VOL_DISCRETE),
                                   generate_volume_signature(config.trajectory_tgt , 
                                                             config.max_point, 
                                                             config.min_point, 
                                                             GRID_SIZE, 
                                                             VOL_DISCRETE ) );
   cout << sd(vsig.begin(), vsig.end()) <<endl;

   return 0;
}
