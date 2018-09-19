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
	const char*	filename,	/*	input filename			  */
	char*			probname,	/*	output problemname	  */
	int			maxSize		/* maximum size of p.name */
	)
{
	int	i=0;
	int	j=0;
	int	l;

	/*	first find end of string */
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
	int	nthreads = 1;
   int   timelimit = 5000;
   bool  quiet = false;

   while ( (opt = getopt(argc, argv, "f:p:t:qh")) != -1)
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
    //for (i = optind; i < argc; i++) {
    //    printf("arg = %s\n", argv[i]);
    //}

	assert( nthreads > 0 );

	int		m;

	ReadDim( filename, &m);

	assert( m>0 );

	double	*B_ = NULL;		// [m*m], Colmajor
	B_	=	new double[m*m];

	ReadData( filename, m, B_);

	SVPsolver	svps;

	svps.create_probdata( m, B_);
	svps.create_sch( m, B_);
	svps.set_num_thread( nthreads );
	svps.find_min_column();
	svps.compute_bounds();
   svps.set_timelimit( timelimit );
   svps.set_quiet( quiet );

	bool run;
	if( nthreads == 1 ){
		run = svps.solve(false, true, 0.0);
	}else{
		run = svps.p_solve();
	}

   char probname[100];
   getProblemName( filename, probname, 100);

   string com("echo [");
   com += probname;
   com += " @";
   com += to_string(nthreads);
   com += "threads] ";
   if ( svps.get_listsize() == 0 )
      com += "=== OPTIMAL ===";
   else
      com += "== TIME OVER ==";

   com += " TIME: ";
   com += to_string( svps.get_runtime() );

   com += " NOLM: ";
   com += to_string( (int)sqrt(svps.get_bestval()) );
   com += " >> result.txt";

   com += " NODE: ";
   com += to_string( svps.get_nnode() );

   if ( svps.get_listsize() != 0 )
   {
      com += " GAP: ";
      com += to_string( svps.get_gap() );
      com += "%";
   }

   //cout << com << endl;
   system(com.c_str());

	svps.disp_bestsol();
	if( !run ){
      cout << endl;
		cout << "could not solve .." << endl;
		return 0;
	}


//
//
//	QPsolver	 qps;
//	qps.set_dim( m );
//	qps.set_obj_quad( m, Q);
//
//	double	*A;
//	A = new double[m];
//	for(int i=0; i<m; i++){
//		if( i < 10 ){
//			A[i] = -1.0;
//		}else{
//			A[i] = 0.0;
//		}
//	}
//	double	*b;
//	b = new double[1];
//	b[0] = -1.0;
//	//qps.set_ineq( m, 1, A, b);
//
//	double	*C;
//	C = new double[m];
//	for(int i=0; i<m; i++){
//		if( (10<=i) && (i<20) ){
//			C[i] = 1.0;
//		}else{
//			C[i] = 0.0;
//		}
//	}
//	double 	*d;
//	d = new double[1];
//	d[0] = 1.0;
//	//qps.set_equ( m, 1, C, d);
//
//	l[0] = 1.0;
//	u[1] = -1.0;
//	qps.set_lb( m, l);
//	qps.set_ub( m, u);
//
//	qps.disp_prob();
//	
//	double	*warm;
//	warm = new double[m];
//
//	warm[0] = 3.0;
//	warm[1] = -3.0;
//	for(int i=2; i<=10; i++){
//		warm[i] = 1.0;
//	}
//	for(int i=11; i<=19; i++){
//		warm[i] = 0.0;
//	}
//	for(int i=20; i<m; i++){
//		warm[i] = 1.0;
//	}
//
//	for(int i=30; i<35; i++){
//		warm[i] = l[i];
//	}
//	for(int i=35; i<m; i++){
//		warm[i] = u[i];
//	}
//
//	qps.set_warm( m, warm);
//	qps.solve();
//
	delete[] B_;
//	delete[] nrm;
//	delete[] l;
//	delete[] u;
//
//	delete[] Q;
//	delete[] A;
//	delete[] b;
//	delete[] C;
//	delete[] d;
//	delete[] warm;

	return 0;
}
