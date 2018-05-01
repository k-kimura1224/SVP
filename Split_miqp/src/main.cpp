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
//#include "vector.h"
//#include "qpsolver.h"
#include "svpsolver.h"

#define debug 1

using namespace std;

int main( int argc, char** argv){

	cout << fixed << setprecision(3);
	if( argc != 3 ){
		cerr << "Error: commandline arguments" << endl;
		return -1;
	}
	int	nthreads = atoi( argv[argc-1] );
	assert( nthreads > 0 );

	int		m;

	ReadDim( argv[1], &m);

	assert( m>0 );

	double	*B_ = NULL;		// [m*m], Colmajor
	B_	=	new double[m*m];

	ReadData( argv[1], m, B_);

	SVPsolver	svps;

	svps.create_probdata( m, B_);
	svps.create_sch( m, B_);
	svps.set_num_thread( nthreads );
	svps.find_min_column();
	svps.compute_bounds();

	bool run;
	if( nthreads == 1 ){
		run = svps.solve(false, true, 0.0, TIMELIMIT);
	}else{
		run = svps.p_solve();
	}

	if( !run ){
		cout << "could not solve .." << endl;
		return 0;
	}

	svps.disp_bestsol();

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
