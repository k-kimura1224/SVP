#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>
#include <math.h>
#include <list>
#include <iomanip>
#include <omp.h>

#include "svpsolver.h"
#include "probdata.h"
#include "vector.h"
#include "solution.h"
#include "stopwatch.h"
#include "node.h"
#include "Schmidt_manager.h"

//#include <omp.h>

using namespace std;

#define debug_class  0
#define debug  0

SVPsolver::SVPsolver(){ // default constructor
#if debug_class
   cout << "SVPsolver: default constructor" << endl;
#endif
   bestval = -1.0;
   ub = nullptr;
   lb = nullptr;
   GLB = - 1.0e+10;

   listsize = 0;
   nnode = 0;
   nthreads = -1;

   order = nullptr;

   _Appfac = 0.0;
   Appfac = 0.0;

   quiet = false;

   norm = nullptr;
   TIMELIMIT = 86400;

   epsilon = 1.0e-12;
}

SVPsolver::SVPsolver( const SVPsolver &source )
{  // copy constructor
#if debug_class
   cout << "SVPsolver: copy constructor" << endl;
#endif

   probdata = source.probdata;
   GLB = source.GLB;
   bestval = source.bestval;
   bestsol = source.bestsol;
   pool = source.pool;
   _Appfac = source._Appfac;
   Appfac = source.Appfac;
   listsize = source.listsize;
   stopwatch = source.stopwatch;
   testwatch = source.testwatch;
   nnode = source.nnode;
   sch = source.sch;
   quiet = source.quiet;
   nthreads = source.nthreads;
   oa_cpool = source.oa_cpool;
   TIMELIMIT = source.TIMELIMIT;
   epsilon = source.epsilon;

   if ( probdata.get_m() > 0 )
   {
      assert( source.ub != nullptr );
      assert( source.lb != nullptr );

      int m = probdata.get_m();

      ub = new double[m];
      lb = new double[m];

      Copy_vec( source.ub, ub, m );
      Copy_vec( source.lb, lb, m );

      if( source.order != nullptr )
      {
         order = new int[m];
         for (int i = 0; i < m; i++ ){
            order[i] = source.order[i];
         }
      }else{
         order = nullptr;
      }

      if( source.norm != nullptr ){
         norm = new double[m];
         Copy_vec( source.norm, norm, m);
      }else{
         norm = nullptr;
      }
   }else{
      ub = nullptr;
      lb = nullptr;
      order = nullptr;
      norm = nullptr;
   }

   if( listsize > 0 ){
      copy( source.NodeList.begin(), source.NodeList.end(), back_inserter(NodeList));
   }


}

// assignment operator
SVPsolver& SVPsolver::operator=( const SVPsolver& source )
{
#if debug_class
   cout << "SVPsolver: assignment operator" << endl;
#endif

   if( this != &source ){
      probdata = source.probdata;
      GLB = source.GLB;
      bestval = source.bestval;
      bestsol = source.bestsol;
      pool = source.pool;
      _Appfac = source._Appfac;
      Appfac = source.Appfac;
      listsize = source.listsize;
      stopwatch = source.stopwatch;
      testwatch = source.testwatch;
      nnode = source.nnode;
      sch = source.sch;
      quiet = source.quiet;
      nthreads = source.nthreads;
      oa_cpool = source.oa_cpool;
      TIMELIMIT = source.TIMELIMIT;
      epsilon = source.epsilon;

      if( probdata.get_m() > 0 ){
         assert( source.ub != nullptr );
         assert( source.lb != nullptr );

         int m = probdata.get_m();

         delete[] ub;
         delete[] lb;
         ub = new double[m];
         lb = new double[m];

         Copy_vec( source.ub, ub, m);
         Copy_vec( source.lb, lb, m);

         if( source.order != nullptr ){
            delete[] order;
            order = new int[m];
            for(int i=0; i<m; i++){
               order[i] = source.order[i];
            }
         }else{
            order = nullptr;
         }

         if( source.norm != nullptr ){
            delete[] norm;
            norm = new double[m];
            Copy_vec( source.norm, norm, m);
         }else{
            norm = nullptr;
         }
      }else{
         ub = nullptr;
         lb = nullptr;
         order = nullptr;
         norm = nullptr;
      }


      if( listsize > 0 ){
         copy( source.NodeList.begin(), source.NodeList.end(), back_inserter(NodeList));
      }
   }

   return *this;
}

// destructor
SVPsolver::~SVPsolver()
{
#if debug_class
   cout << "SVPsolver: destructor" << endl;
#endif
   delete[] ub;
   delete[] lb;
   list<NODE>().swap(NodeList);
   delete[] order;
   delete[] norm;
   ub = nullptr;
   lb = nullptr;
   order = nullptr;
   norm = nullptr;

}

void SVPsolver::create_probdata(
   int      m,
   double   *B_
   )
{
   assert( m > 0 );
   assert( B_ != nullptr );

   double   *B;      // [m*m],
   double   *Q;      // [m*m],

   B  =  new double[m*m];
   Q  =  new double[m*m];

   TraMat( m, m, B_, B);

   Gen_ZeroVec( m * m, Q);
   Com_mat_AtA( B_, m, m, Q);

   probdata.set_data( m, B, B_, Q);

   delete[] B;
   delete[] Q;

   ub = new double[m];
   lb = new double[m];

   for(int i = 0; i < m; i++ )
   {
      ub[i] = 1.0e+10;
      lb[i] = - 1.0e+10;
   }

   _Appfac = tgamma( ((double)m/2.0) + 1 );
   _Appfac = pow( _Appfac, 1.0/(double)m );

   double absdet = fabs( determinant( B_, m));
   absdet = pow( absdet, 1.0/(double)m );
   _Appfac *= absdet;
   _Appfac /= sqrt( M_PI );
   pool.alloc( 100 );
}

void SVPsolver::create_sch(
   int      m,
   double   *B_
   )
{
   sch.setup( m, B_);
}

void SVPsolver::disp_bestsol()
{
   int   m;
   double *val = nullptr;

   m = probdata.get_m();
   val = bestsol.get_solval();

   assert( val != nullptr );

   cout << endl;
   if( listsize == 0 ){
      cout << "SVPSOLVER found an optimal solution" << endl;
   }else{
      cout << "-- TIMEOVER --" << endl;
   }
   cout << endl;

   cout << "time: " << stopwatch.get_result() << "s" << endl;
   cout << "best value: " << bestval << endl;
   cout << "norm: " << sqrt( bestval ) << endl;
   cout << "AF: " << Appfac << endl;

   if( listsize != 0 ){
      cout << "lower bound: " << GLB << endl;
      cout << "gap: " << 100*(bestval - GLB)/bestval << "%" << endl;
   }
   cout << "Nodes: " << nnode << endl;

   cout << "best solution:" << endl;
   for(int i=0; i<m; i++){
      if( val[i] != 0.0 ){
         cout << " " << i << ": " << val[i] << endl;
      }
   }



}

void  SVPsolver::compute_bounds()
{

   int      m = probdata.get_m();
   double   *B = probdata.get_B();

   assert( bestval > 0.0 );
   assert( m > 0 );

   double *e;
   double *a;
   e = new double[m];
   a = new double[m];
   Gen_ZeroVec( m, e);

   int r;
   auto sqrt_bestval = sqrt( bestval );

   for ( int i = 0; i < m; i++ )
   {
      assert( e[i] == 0.0 );
      e[i] = 1.0;

      r = Com_LS_dgesv( B, e, m, a);
      if ( r != 0 )
      {
         cout << "error: r = " << r << endl;
         exit(-1);
      }
      ub[i] =  floor( sqrt_bestval * Com_nrm( a, m) );
      lb[i] =  -ub[i];

      e[i] = 0.0;
   }

   delete[] e;
   delete[] a;

}


void  SVPsolver::tighten_bounds(
      int      memo,
      double*  set_lb,
      double*  set_ub
      )
{
   int      m = probdata.get_m();

   double*  B_ = probdata.get_B_();
   double*  Q = probdata.get_Q();

   auto     dim = m - ( memo + 1 );
   auto     mdim = m * dim;
   auto     dimdim = dim * dim;

   double*  subB = new double[mdim];
   double*  subBB = new double[dimdim];
   double*  e = new double[dim];
   double*  buf_a = new double[dim];
   double*  buf_b = new double[m];

   Gen_ZeroVec( dim, e );

   auto constant = m * ( memo + 1 );
   for ( auto i = 0; i < mdim; i++ )
      subB[i] = B_[constant + i];

   int ct = 0;

   for ( auto i = memo + 1; i < m; i++ )
   {
      for ( auto j = memo + 1; j < m; j++ )
      {
         subBB[ct] = Q[i*m + j];
         assert( subBB[ct] == Com_dot( &subB[(i-(memo+1))*m], &subB[(j-(memo+1))*m], m));
         ct++;
      }
   }

   int r;
   auto sqrt_bestval = sqrt( bestval );

   double buf;
   constant = memo + 1;

   for ( auto i = 0; i < dim; i++ )
   {
      assert( e[i] == 0.0 );
      e[i] = 1.0;

      r = Com_LS_dgesv( subBB, e, dim, buf_a);
      if ( r != 0 )
      {
         cout << "error: r = " << r << endl;
         exit(-1);
      }

      Com_mat_Ax( subB, m, dim, buf_a, buf_b);

      buf =  floor( sqrt_bestval * Com_nrm( buf_b, m) );
      set_ub[constant+i] =  buf;
      set_lb[constant+i] =  - buf;

      e[i] = 0.0;

   }

   delete[] subB;
   delete[] subBB;
   delete[] e;
   delete[] buf_a;
   delete[] buf_b;
}


void  SVPsolver::disp_log(
   int   k,
   RelaxResult r,
   int   index,
   int   cutoff
   )
{
   assert( k >= 0 );
   assert( k < listsize );

   list<NODE>::iterator it = NodeList.begin();
   advance( it, k);

   cout << stopwatch.get_time() << "s";
   cout << " [" << it->get_index() << "/";
   cout << listsize << "(" << index-1 << ")] ";
   cout << "dpt:" << it->get_dpt() << " ";
   cout << "[ " << it->get_lowerbound() << ", ";
   cout << GLB << ", ";
   cout << bestval << ", ";
   cout << 100*(bestval - GLB)/bestval << "%]";

   cout << " cutoff:" << cutoff;
   cout << " AF:";
   cout << fixed << setprecision(3) << Appfac;

   cout << endl;
}

double SVPsolver::compute_objval(
   double   *x
   )
{
   int      m = probdata.get_m();
   double   *Q = probdata.get_Q();

   assert( m > 0 );
   assert( Q != nullptr );
   assert( x != nullptr );

   double   obj = 0.0;

   double   *Qx; //[m]
   Qx = new double[m];

   Gen_ZeroVec( m, Qx);
   Com_mat_Ax( Q, m, m, x, Qx);

   obj += Com_dot( x, Qx, m);

   delete[] Qx;

   return obj;
}
