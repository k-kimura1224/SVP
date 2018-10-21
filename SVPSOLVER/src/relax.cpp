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
#include "cut.h"

using namespace std;

#define debug  0

RelaxResult SVPsolver::SVPSsolveRelaxation(
      NODE&          node,
      const double*  vars_localub,
      const double*  vars_locallb
   )
{

   if ( bestval - node.get_lowerbound() < 1.0 )
      return INFEASIBLE;

   RelaxResult result = SVPSsolveRelaxationINT( node, vars_localub, vars_locallb );

   return result;
}

RelaxResult SVPsolver::SVPSsolveRelaxationINT(
      NODE&          node,
      const double*  vars_localub,
      const double*  vars_locallb
   )
{
   // copy
   const auto m = probdata.get_m();
   const auto Q = probdata.get_Q();
   const auto ep = epsilon;
   auto& qd = qpdata;

   auto warm = node.get_relaxsolval();

   assert( m > 0 );
   assert( Q != nullptr );
   assert( ep < 1.0 && ep > 0.0 );

   vector<int> nofixed_index;
   nofixed_index.reserve( m );

   double *relaxvals = nullptr;

   for ( auto i = 0; i < m; ++i )
   {
      if ( !Equal( vars_localub[i], vars_locallb[i], ep )  )
         nofixed_index.push_back( i );
   }

   nofixed_index.shrink_to_fit();

   const int n = (int) nofixed_index.size();

   // all values are fixed
   if ( n == 0 )
   {
      relaxvals = new double[m];
      for ( auto i = 0; i < m; ++i )
         relaxvals[i] = vars_localub[i];

      const auto objval = compute_objval( relaxvals );
      SOLUTION solution;

      solution.set_sol( m, relaxvals, objval );
      SVPStrySol( solution, true, true, nullptr );

      delete[] relaxvals;

      return GETINTEGER;
   }

   double* subQ = nullptr;
   double* subu = nullptr;
   double* subl = nullptr;
   double* subw = nullptr;
   double* p = nullptr;
   double c_term = 0.0;

   QPsolver  qps;

   int   ct = 0;

   qps.set_dim( n );

   if( n == m )
   {
      qps.set_obj_quad( m, Q );
      qps.set_lb( m, vars_locallb );
      qps.set_ub( m, vars_localub );
      qps.set_warm( m, warm );
   }
   else
   {
      assert( m - 1 >= n );

      if( qd.check_allocation() == false )
         qd.alloc( m - 1 );

      subQ = qd.get_Qmat();
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
         p = qd.get_pvec();
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
         qps.set_obj_quad( n, subQ );
      }

      subu = qd.get_uvec();
      subl = qd.get_lvec();
      subw = qd.get_wvec();

      assert( subu != nullptr );
      assert( subl != nullptr );
      assert( subw != nullptr );

      ct = 0;

      for ( auto i : nofixed_index )
      {
         subu[ct] = vars_localub[i];
         subl[ct] = vars_locallb[i];
         ct++;
      }

      assert( ct == n );

      ct = 0;

      for ( auto i : nofixed_index )
      {
         subw[ct] = warm[i];
         ct++;
      }

      qps.set_ub( n, subu );
      qps.set_lb( n, subl );
      qps.set_warm( n, subw );
   }

   // using oa_cut{{
   double   *A = nullptr;
   double   *b = nullptr;
   //if( CUT_OA == true ){
   //   assert( oa_cpool.get_size() > 0 );
   //   cout << "stop!!!!!!!" << endl;
   //   exit(1);

   //   int   k  = oa_cpool.get_size() * 2;
   //   A = new double[k*m];
   //   b = new double[k];

   //   CUT_PLANE *cut;
   //   ct = 0;
   //   for(int i=0; i<oa_cpool.get_size(); i++){
   //      cut = oa_cpool.get_cut( i );

   //      assert( cut->get_coef() != nullptr );
   //      assert( cut->get_lb() != nullptr );
   //      assert( cut->get_ub() != nullptr );
   //      assert( cut->get_ub() > cut->get_lb() );

   //      Copy_vec( cut->get_coef(), &A[ct*m], m);
   //      for(int j=0; j<m; j++){
   //         A[(ct*m)+m+j] = - A[(ct*m)+j];
   //      }
   //      b[ct] = *(cut->get_ub());
   //      b[ct+1] = - *(cut->get_lb());

   //      ct += 2;
   //   }
   //   assert( ct == k );

   //   qps.set_ineq( m, k, A, b);
   //}
   // }} using oa_cut

   //qps.disp_prob();

   qps.set_ep( 1e-12 );

   qps.solve( nullptr );
   //qps.solve( &testwatch );

   auto qpbestvals = qps.get_bestsol();
   ct = 0;

   relaxvals = new double[m];
   for ( int i = 0; i < m; ++i )
   {
      if( !Equal( vars_localub[i], vars_locallb[i], ep ) )
         relaxvals[i] = qpbestvals[ct++];
      else
         relaxvals[i] = vars_locallb[i];
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

      if( Equal( objval, qpbestval, ep ) )
         result = GETINTEGER;

      SOLUTION solution;
      solution.set_sol( m, relaxvals, objval);
      SVPStrySol( solution, true, true, nullptr );

      // Do not use relaxvals after here
   }

   delete[] A;
   delete[] b;
   delete[] relaxvals;

   return result;
}

