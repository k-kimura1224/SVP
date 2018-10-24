#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <time.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <unistd.h>

#include "read_data.h"
//#include "vector.h"
//#include "qpsolver.h"
#include "svpsolver.h"


using namespace std;

static
void getProblemName(
   const char* filename,   /* input filename         */
   char*       probname,   /* output problemname     */
   int         maxSize     /* maximum size of p.name */
   )
{
   int   i=0;
   int   j=0;
   int   l;

   /* first find end of string */
   while( filename[i]!=0 )
      ++i;
   l = i;

   /* go back until '.' or '/' or '\' appears */
   while( (i>0) && (filename[i]!='.') && (filename[i]!='/') && (filename[i]!='\\'))
      --i;

   /* if we found '.', search for '/' or '\\' */
   if( filename[i]=='.' ){
      l = i;
      while( (i>0) && (filename[i]!='/') && (filename[i]!='\\') )
         --i;
   }

   /* crrect counter */
   if( (filename[i]=='/') || (filename[i]=='\\') )
      ++i;

   /* copy name */
   while( (i<l) && (filename[i]!=0) ){
      probname[j++] = filename[i++];
      if( j>maxSize-1)
         exit(0);
   }
   probname[j] = 0;

}

int main( int argc, char** argv){

   cout << fixed << setprecision(3);

   int   opt;
   opterr = 0;

   char* filename = nullptr;
   int   nthreads = 1;
   int   timelimit = 5000;
   bool  quiet = false;
   int   memory = 8;
   bool  cutmode = false;
   bool  enumeration = false;

   while ( (opt = getopt(argc, argv, "f:p:t:qhm:ce")) != -1)
   {
      switch ( opt )
      {
         case 'f':
            filename = optarg;
            break;

         case 'p':
            nthreads = atoi( optarg );
            break;

         case 't':
            timelimit = atoi( optarg );
            break;

         case 'q':
            quiet = true;
            break;

         case 'h':
            cout << "Usage: " << argv[0] << " [-f filename] [-p nthreads] [-t timelimit]" << endl;
            break;

         case 'm':
            memory = atoi ( optarg );
            break;

         case 'c':
            cutmode = true;
            break;

         case 'e':
            enumeration = true;
            break;

         default: /* '?' */
            cout << "Usage: " << argv[0] << " [-f filename] [-p nthreads] [-t timelimit]" << endl;
            break;
        }
    }

   if ( filename == nullptr )
   {
      cout << "Error: no input file" << endl;
      cout << "Usage: " << argv[0] << " [-f filename] [-p nthreads] [-t timelimit]" << endl;
      return -1;
   }

   assert( nthreads > 0 );
   assert( memory > 0 );

   memory /= nthreads;

   int      m;

   ReadDim( filename, &m);

   assert( m>0 );

   double   *B_ = NULL;    // [m*m], Colmajor
   B_ =  new double[m*m];

   ReadData( filename, m, B_);

   SVPsolver   svps;
   svps.SVPSsetCutMode( cutmode );
   svps.SVPSsetEnum( enumeration );
   svps.SVPSsetup( m, B_, nthreads, timelimit, memory, quiet,
         false, true, true, true, true, true );

   bool run;
   if( nthreads == 1 )
   {
      run = svps.SVPSsolve();
   }
   else
   {
      run = svps.SVPSparasolve();
   }

   char probname[100];
   getProblemName( filename, probname, 100);

   string com("echo [");
   com += probname;
   com += " @";
   com += to_string(nthreads);
   com += "threads] ";
   if ( svps.SVPSgetListsize() == 0 )
      com += "=== OPTIMAL ===";
   else
      com += "== LIMIT ==";

   com += " TIME: ";
   com += to_string( svps.get_runtime() );

   com += " NOLM: ";
   com += to_string( (int)sqrt(svps.SVPSgetBestval()) );

   com += " NODE: ";
   com += to_string( svps.SVPSgetNnode() );

   if ( svps.SVPSgetListsize() != 0 )
   {
      com += " GAP: ";
      com += to_string( svps.get_gap() );
      com += "%";
   }

   //struct rusage r;
   //if (getrusage(RUSAGE_SELF, &r) != 0) {
   //      /*Failure*/
   //}

   //int mem = floor( (r.ru_maxrss/1024)/1024 );

   //com += " MAXMEM: ";
   //com += to_string( mem );
   //com += "MB";

   com += " >> result.txt";

   //cout << com << endl;
   system(com.c_str());

   svps.disp_bestsol();
   if( !run ){
      cout << endl;
      cout << "could not solve .." << endl;
      return 0;
   }

   delete[] B_;

   return 0;
}
