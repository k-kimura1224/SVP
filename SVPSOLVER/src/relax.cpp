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

#include "svpsolver.h"
#include "node.h"
#include "probdata.h"
#include "qpdata.h"
#include "qpsolver.h"
#include "solution.h"
#include "vector.h"
#include "Schmidt_manager.h"
#include "cut_pool.h"
#include "cut.h"

using namespace std;

#define debug  0

RelaxResult SVPsolver::solve_relaxation(
   int   selnodeindex
   )
{
   RelaxResult result;
   auto it = NodeList.begin();
   advance( it, selnodeindex );

   if ( bestval - it->get_lowerbound() < 1.0 )
      return INFEASIBLE;

   if ( it->get_zero() == true )
   {
      result = solve_relaxation_BIN_sch( selnodeindex );
   }
   else
   {
      result = solve_relaxation_INT( selnodeindex );
   }

   return result;
}

RelaxResult SVPsolver::solve_relaxation_BIN_sch(
   int   sel
   )
{
   int   m = probdata.get_m();
   auto  it = NodeList.begin();

   advance( it, sel );

   double *u = it->get_ub();
   double *l = it->get_lb();

   int ct = 0;

   for ( int i = 0; i < m; i++ )
   {
      if( u[i] == 0.0 && l[i] == 0.0 )
      {
         sch.set_z_i( i, false );
         ct++;
      }
      else
      {
         sch.set_z_i( i, true );
      }
   }

   if ( ct == m ) return INFEASIBLE;

   assert( sch.get_n() > 0 );

   RelaxResult result;

   sch.compute_GS();

   auto schmin = sch.get_min();

   if( it->get_lowerbound() < schmin )
      it->set_lowerbound( schmin );

   if ( schmin >= bestval )
   {
      result = INFEASIBLE;
#if debug
      cout << "*********************************" << endl;
      cout << "*   Schmidt_manager  #zero: " << ct << endl;
      cout << "*********************************" << endl;
#endif
   }
   else
   {
      result = FEASIBLE;
   }

   return result;
}

RelaxResult SVPsolver::solve_relaxation_INT(
   int   sel
   )
{
   auto it = NodeList.begin();
   advance( it, sel );

   int      m = probdata.get_m();
   double   *Q = probdata.get_Q();
   double   *u = it->get_ub();
   double   *l = it->get_lb();
   double   *warm = it->get_warm();

   int   n = 0;

   for ( int i = 0; i < m; i++ )
   {
      if ( l[i] != u[i] ) n++;
   }

   if ( n == 0 )
      return INFEASIBLE;


   double*  subQ = nullptr;
   double*  subu = nullptr;
   double*  subl = nullptr;
   double*  subw = nullptr;
   double*  p = nullptr;
   double   c_term = 0;

   QPsolver  qps;

   int   ct = 0;

   qps.set_dim( n );

   if( n == m )
   {
      qps.set_obj_quad( m, Q);
      qps.set_lb( m, l);
      qps.set_ub( m, u);
      qps.set_warm( m, warm);
   }
   else
   {
      assert( m - 1 >= n );

      if( qpdata.check_allocation() == false )
         qpdata.alloc( m - 1 );


      subQ = qpdata.get_Qmat();
      assert( subQ != nullptr );

      ct = 0;
      for ( int i = 0; i < m; i++ )
      {
         if ( l[i] != u[i] )
         {
            for(int j=0; j<m; j++)
            {
               if ( l[j] != u[j] )
               {
                  subQ[ct++] = Q[i+(j*m)];
               }
            }
         }
      }

      assert( ct == n*n );
   // printM( n, n, subQ);
   // cout << endl;
   // cout << "[1,0]=" << subQ[0+(1*n)] <<endl;
   // cout << "[0,1]=" << subQ[1+(0*n)] <<endl;
      assert( subQ[0+(1*n)] == subQ[1+(0*n)] );

      if( it->get_sumfixed() != nullptr ){
         p = qpdata.get_pvec();
         assert( p != nullptr );

         ct = 0;
         for(int i=0; i<m; i++){
            if( l[i] != u[i] ){
               p[ct++] = 2*Com_dot( probdata.get_bvec(i), it->get_sumfixed(), m);
            }
         }
         assert( ct == n );

         qps.set_obj( n, subQ, p);
         c_term = Com_dot( it->get_sumfixed(), it->get_sumfixed(), m);
      }else{
         qps.set_obj_quad( n, subQ);
      }

      subu = qpdata.get_uvec();
      subl = qpdata.get_lvec();
      subw = qpdata.get_wvec();

      assert( subu != nullptr );
      assert( subl != nullptr );
      assert( subw != nullptr );

      ct = 0;

      for(int i=0; i<m; i++){
         if( l[i] != u[i] ){
            subu[ct] = u[i];
            subl[ct] = l[i];
            subw[ct] = warm[i];
            ct++;
         }
      }
      assert( ct == n );

      qps.set_ub( n, subu);
      qps.set_lb( n, subl);
      qps.set_warm( n, subw );


   }




   // using oa_cut{{
   double   *A = nullptr;
   double   *b = nullptr;
   if( CUT_OA == true ){
      assert( oa_cpool.get_size() > 0 );
      cout << "stop!!!!!!!" << endl;
      exit(1);

      int   k  = oa_cpool.get_size() * 2;
      A = new double[k*m];
      b = new double[k];

      CUT_PLANE *cut;
      ct = 0;
      for(int i=0; i<oa_cpool.get_size(); i++){
         cut = oa_cpool.get_cut( i );

         assert( cut->get_coef() != nullptr );
         assert( cut->get_lb() != nullptr );
         assert( cut->get_ub() != nullptr );
         assert( cut->get_ub() > cut->get_lb() );

         Copy_vec( cut->get_coef(), &A[ct*m], m);
         for(int j=0; j<m; j++){
            A[(ct*m)+m+j] = - A[(ct*m)+j];
         }
         b[ct] = *(cut->get_ub());
         b[ct+1] = - *(cut->get_lb());

         ct += 2;
      }
      assert( ct == k );

      qps.set_ineq( m, k, A, b);
   }
   // }} using oa_cut

   //qps.disp_prob();

   qps.set_ep( 1e-12 );

   testwatch.start();
   qps.solve();
   testwatch.stop();

   double *relaxvals = new double[m];
   double *qpbestvals = qps.get_bestsol();
   ct = 0;
   for(int i=0; i<m; i++){
      if( l[i] != u[i] ){
         relaxvals[i] = qpbestvals[ct++];
      }else{
         relaxvals[i] = l[i];
      }
   }
   assert( ct == n );

   it->set_relaxsolval( relaxvals );

   double qpbestval = qps.get_bestval() + c_term;

   if( it->get_lowerbound() < qpbestval ){
      it->set_lowerbound( qpbestval );
   }

   RelaxResult result;

   if( bestval - qpbestval  < 1.0 )
   {
      result = INFEASIBLE;
   }else{
      result = FEASIBLE;

      double *vals;

      vals = new double[m];

      for(int i=0; i<m; i++){
         vals[i] = round( relaxvals[i] );
      }

      double objval = compute_objval( vals );

      if( objval == qpbestval ){
         result = GETINTEGER;
      }

      SOLUTION solution;
      solution.set_sol( m, vals, objval);
      pool.add_solution( solution );

      if( bestval > objval ){
         bestval = objval;
         bestsol = solution;
         result = UPDATE;
         Appfac = sqrt( bestval ) / _Appfac;
      }

      delete[] vals;
   }

   delete[] A;
   delete[] b;
   //delete[] subQ;
   //delete[] subu;
   //delete[] subl;
   //delete[] subw;
   //delete[] p;
   delete[] relaxvals;

   return result;
}

