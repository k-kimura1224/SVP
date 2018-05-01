#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <time.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>

#include "read_data.h"
#include "svpsolver.h"

#define debug 1

using namespace std;

int main( int argc, char** argv){

   cout << fixed << setprecision(3);

   if( argc != 4 )
   {
      cerr << "Error: commandline arguments" << endl;
      return -1;
   }

   int m;
   auto tlimit = atoi( argv[2] );
   auto selection = atoi( argv[3] );

   ReadDim( argv[1], &m);

   assert( m > 0 );

   double *B_ = nullptr;    // [m*m], Colmajor
   B_ = new double[ m * m ];

   ReadData( argv[1], m, B_);

   SVPsolver   svps;
   auto nthreads = 1;

   svps.create_probdata( m, B_ );
   svps.create_sch( m, B_ );
   svps.set_num_thread( nthreads );
   svps.find_min_column();
   svps.compute_bounds();

   bool run = svps.solve_cplex( tlimit, selection );

   if( !run ){
		cout << "could not solve .." << endl;
		return 0;
	}

	svps.disp_bestsol();

   delete[] B_;

   return 0;
}
