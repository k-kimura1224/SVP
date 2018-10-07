#ifndef PROBDATA_H__
#define PROBDATA_H__

using namespace std;

class PROB_DATA{
	int			m;			// the number of nodes in L
	double*     B;		// [m*m]
	double*     B_;		// [m*m], Colmajor
	double*     Q;		// [m*m], B'B
	public:

		PROB_DATA();											// default constructor
		PROB_DATA( const PROB_DATA &source );			// copy constructor
		PROB_DATA& operator=( const PROB_DATA& );		// assignment operator
		~PROB_DATA();											// destructor

		void set_data( const int s_m, const double *s_B, const double *s_B_, const double *s_Q );

		int		get_m() const { return m; }
		double*	get_B() const { return B; }
		double*	get_B_() const { return B_; }
		double*	get_Q() const { return Q; }
		double*	get_bvec( const int k ) const { return B_+(k*m); }
};

#endif
