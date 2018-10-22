#ifndef QPSOLVER_H_
#define QPSOLVER_H_

#include <assert.h>
#include "testwatch.h"

using namespace std;


#define  SOLVE_PROBLEM_ACTIVE 1
// 0:dgesv
// 1:dsysv
// 2:dspsv? unimplimented

/*********************************************
   min   x'Qx + p'x
   s.t.  Ax <= b
         Cx = d
         l <= x <= u
         x \in R^n
*********************************************/

class QPsolver{
   const double   *Q;   // [n*n]
   const double   *p;   // [n]
   const double   *A;   // [m*n]
   const double   *b;   // [m]
   const double   *C;   // [k*n]
   const double   *d;   // [k]
   const double   *l;   // [n]
   const double   *u;   // [n]
   int      n;
   int      m;
   int      k;
   const double   *warm;
   double   *bestsol;
   double   bestval;

   double   ep;

   // for solve_activeset_onlybounds
   double*  H;          // max = n*n
   double*  r;          // max = n
   double*  z;          // max = n
   int*     indexes;    // [n], ({nonactive}{active})
   double*  fixed;      // max = n

   // for linear system
   int (*comls)( const double*, const double*, int, double*);

   public:

   QPsolver();                               // default constructor
   QPsolver( const QPsolver &source );       // copy constructor
   QPsolver& operator=( const QPsolver& );   // assignment operator
   ~QPsolver();                              // destructor

   // set
   void  set_dim( const int s_n ) { n = s_n; }
   void  set_obj( const int s_n, const double* s_Q, const double* s_p ){
      assert( n == s_n );
      Q = s_Q;
      p = s_p;
   }
   void  set_obj_quad( const int s_n, const double* s_Q ){
      assert( n == s_n );
      Q = s_Q;
   }
   void  set_obj_line( const int s_n, const double *s_p ){
      assert( n == s_n );
      p = s_p;
   }
   void  set_ineq( const int s_n, const int s_m, const double* s_A, const double *s_b ){
      assert( n == s_n );
      m = s_m;
      A = s_A;
      b = s_b;
   }
   void  set_equ( const int s_n, const int s_k, const double* s_C, const double* s_d ){
      assert( n == s_n );
      k = s_k;
      C = s_C;
      d = s_d;
   }
   void  set_lb( const int s_n, const double* s_l ){
      assert( n == s_n );
      l = s_l;
   }
   void  set_ub( const int s_n, const double* s_u ){
      assert( n == s_n );
      u = s_u;
   }
   void  set_warm( const int s_n, const double* s_warm ){
      assert( n == s_n );
      warm = s_warm;
   }

   void  disp_prob() const;
   void  solve( TESTWATCH* testwatch );
   int   solve_activeset( double *x_new, bool *W_ineq, bool *W_u, bool *W_l);
   int   solve_activeset_onlybounds( double *x_new, bool *W_u, bool *W_l);
   double   compute_stepsize( double *x, double *dx, bool *W_ineq, bool *W_u, bool *W_l);
   bool  update_activeset( double *x_new, bool *W_ineq, bool *W_u, bool *W_l);
   double   compute_objval( double *x );

   double   get_bestval(){ return bestval; }
   double*  get_bestsol(){ return bestsol; }

   void     set_ep( double s_ep ){ ep = s_ep; }
   double   get_ep(){ return ep; }
};

#endif
