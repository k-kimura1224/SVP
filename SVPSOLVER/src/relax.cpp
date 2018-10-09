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

RelaxResult SVPsolver::SVPSsolveRelaxation(
      NODE&    node
   )
{
   RelaxResult result;

   if ( bestval - node.get_lowerbound() < 1.0 )
      return INFEASIBLE;

   if ( node.get_zero() == true )
   {
      result = SVPSsolveRelaxationBIN( node );
   }
   else
   {
      result = SVPSsolveRelaxationINT( node );
   }

   return result;
}

RelaxResult SVPsolver::SVPSsolveRelaxationBIN(
      NODE&    node
   )
{
   const auto  m = probdata.get_m();

   const auto *u = node.get_ub();
   const auto *l = node.get_lb();

   int ct = 0;

   for ( int i = 0; i < m; ++i )
   {
      if( u[i]  || l[i]  )
      {
         sch.set_z_i( i, true );
      }
      else
      {
         sch.set_z_i( i, false );
         ct++;
      }
   }

   if ( ct == m ) return INFEASIBLE;

   assert( sch.get_n() > 0 );

   RelaxResult result;

   sch.compute_GS();

   auto schmin = sch.get_min();

   if( node.get_lowerbound() < schmin )
      node.set_lowerbound( schmin );

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

RelaxResult SVPsolver::SVPSsolveRelaxationINT(
      NODE&    node
   )
{
   auto  m = probdata.get_m();
   auto  Q = probdata.get_Q();
   auto  u = node.get_ub();
   auto  l = node.get_lb();
   auto  warm = node.get_relaxsolval();

   vector<int> nofixed_index;
   nofixed_index.reserve( m );

   for ( int i = 0; i < m; i++ )
   {
      if ( l[i] != u[i] )
         nofixed_index.push_back( i );
   }

   const int   n = (int) nofixed_index.size();

   if ( n == 0 )
      return INFEASIBLE;


   double*  subQ = nullptr;
   double*  subu = nullptr;
   double*  subl = nullptr;
   double*  subw = nullptr;
   double*  p = nullptr;
   double   c_term = 0.0;
   double*  d_l = nullptr;
   double*  d_u = nullptr;

   QPsolver  qps;

   int   ct = 0;

   qps.set_dim( n );

   if( n == m )
   {
      d_l = new double[m];
      d_u = new double[m];
      for ( int i = 0; i < m; i++ )
      {
         d_l[i] = l[i];
         d_u[i] = u[i];
      }
      qps.set_obj_quad( m, Q);
      qps.set_lb( m, d_l);
      qps.set_ub( m, d_u);
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

      for ( auto i : nofixed_index )
      {
         auto im = i * m;
         for ( auto j : nofixed_index )
         {
            subQ[ct++] = Q[j + im];
         }
      }

      assert( ct == n*n );
      assert( subQ[0+(1*n)] == subQ[1+(0*n)] );

      if( node.get_sumfixed() != nullptr )
      {
         p = qpdata.get_pvec();
         assert( p != nullptr );

         ct = 0;

         auto vec = node.get_sumfixed();
         for ( auto i : nofixed_index )
            p[ct++] = 2 * Com_dot( probdata.get_bvec(i), vec, m );

         assert( ct == n );

         qps.set_obj( n, subQ, p);

         c_term = Com_dot( vec, vec, m);
      }
      else
      {
         qps.set_obj_quad( n, subQ);
      }

      subu = qpdata.get_uvec();
      subl = qpdata.get_lvec();
      subw = qpdata.get_wvec();

      assert( subu != nullptr );
      assert( subl != nullptr );
      assert( subw != nullptr );

      ct = 0;

      for ( auto i : nofixed_index )
      {
         subu[ct] = u[i];
         subl[ct] = l[i];
         subw[ct] = warm[i];
         ct++;
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

   qps.solve( nullptr );
   //qps.solve( &testwatch );

   double *relaxvals = new double[m];
   double *qpbestvals = qps.get_bestsol();
   ct = 0;

   for ( int i = 0; i < m; i++ )
   {
      if( l[i] != u[i] )
         relaxvals[i] = qpbestvals[ct++];
      else
         relaxvals[i] = l[i];
   }

   assert( ct == n );

   node.set_relaxsolval( relaxvals );

   double qpbestval = qps.get_bestval() + c_term;

   if( node.get_lowerbound() < qpbestval )
      node.set_lowerbound( qpbestval );

   RelaxResult result;

   if( bestval - qpbestval  < 1.0 )
   {
      result = INFEASIBLE;
   }
   else
   {
      result = FEASIBLE;

      for( int i = 0; i < m; i++ )
         relaxvals[i] = round( relaxvals[i] );

      double objval = compute_objval( relaxvals );

      if( objval == qpbestval )
         result = GETINTEGER;

      SOLUTION solution;
      solution.set_sol( m, relaxvals, objval);
      pool.add_solution( solution );

      if( bestval > objval ){
         bestval = objval;
         bestsol = solution;
         result = UPDATE;
         Appfac = sqrt( bestval ) / _Appfac;
      }

      // Do not use relaxvals after here
   }

   delete[] A;
   delete[] b;
   delete[] relaxvals;
   delete[] d_l;
   delete[] d_u;

   return result;
}

