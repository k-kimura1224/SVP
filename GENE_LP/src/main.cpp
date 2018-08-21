#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <time.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <vector>

#include "read_data.h"
#include "svpsolver.h"

#define debug 1

using namespace std;



int main( int argc, char** argv){

   cout << fixed << setprecision(3);

   if( argc != 3 )
   {
      cerr << "Error: commandline arguments" << endl;
      return -1;
   }

   int m;
   auto formulationmode = atoi( argv[2] );

   // formulationmode = 1 : normal 01 foumulation
   // formulationmode = 2 : foumulation by using binary digits

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

   bool gen = svps.gene_lpfile( formulationmode, argv[1] );

   if( !gen )
		cout << "could not gene .." << endl;

   delete[] B_;

   return 0;
}
