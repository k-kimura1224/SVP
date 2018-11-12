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

RelaxResult SVPsolver::SVPSsolveTEST(
      NODE&          node
      )
{

   const auto m = probdata.get_m();
   const auto dim = m - node.get_dpt();

   if ( dim > 30 )
      return R_FEASIBLE;

   assert( !ENUM );
   assert( m > 0 );
   assert( dim > 0 && dim <= m );
   assert( dim != m );
   assert( node.get_lowerbound() > 0.1 );
   assert( node.get_dpt() == (int) node.get_fixedvalues().size() );

   //assert( !VsCB.empty() );
   //assert( !LVsCB.empty() );
   //assert( (int) VsCB.size() == m );
   //assert( (int) LVsCB.size() == m );

   //const auto& vscb_ = VsCB[dim - 1];
   //const auto& lvscb_ = LVsCB[dim - 1];
   const auto& dbm_ = DBMs[dim - 1];
   const auto& rs_ = RSNDBM[dim - 1];
   const auto& rsol = node.geti_rsol();
   const auto size = (int)dbm_.size();

   assert( size == (int)rs_.size() );

   //assert( !vscb_.empty() );
   //assert( !lvscb_.empty() );
   //assert( (int) vscb_.size() == dim );
   //assert( (int) lvscb_.size() == dim );
   assert( !rsol.empty() );
   assert( (int) rsol.size() >= dim );

   double tx;
   double M_i;
   double max_M = 0.0;

   #if debug
   cout << "[ index:" << (int) node.get_index() << "]-------------------------------" << endl;
   #endif

   for ( auto i = 0; i < size; ++i )
   {
      tx = 0;
      for ( auto j = 0; j < dim; ++j )
         tx += rsol[j] * dbm_[i][j];

      M_i = round( tx );
      M_i -= tx;
      M_i *= M_i;
      M_i *= rs_[i];

      if ( max_M < M_i )
         max_M = M_i;

      #if debug
      printf("d:%.8f  ", rs_[i]);
      cout << "tx:" <<tx << "  M_i:" << M_i << "  (" << (int) node.get_lowerbound() << ")" << endl;
      #endif
   }

   #if debug || 0
   //cout << "max_M:" << max_M << endl;
   printf("f(bar_x):%d Max_M:%d bestval:%d\n", (int)node.get_lowerbound(), (int)max_M, (int)bestval);
   #endif

   max_M += node.get_lowerbound();

   #if debug
   cout << "f(rsol) + Max:" << max_M << endl;
   #endif

   if ( ceil( max_M ) < bestval )
      return R_FEASIBLE;
   else
   {
      #if debug
      cout << "max_M:infeasible!" << endl;
      #endif
      return R_INFEASIBLE;
   }
}
RelaxResult SVPsolver::SVPSsolveImprovedRelaxation(
      NODE&          node
      )
{

   const auto m = probdata.get_m();
   const auto dim = m - node.get_dpt();

   assert( !ENUM );
   assert( m > 0 );
   assert( dim > 0 && dim <= m );
   assert( dim != m );
   assert( node.get_lowerbound() > 0.1 );
   assert( node.get_dpt() == (int) node.get_fixedvalues().size() );

   assert( !VsCB.empty() );
   assert( !LVsCB.empty() );
   assert( (int) VsCB.size() == m );
   assert( (int) LVsCB.size() == m );

   //const auto& vscb_ = VsCB[dim - 1];
   const auto& lvscb_ = LVsCB[dim - 1];
   const auto& rsol = node.geti_rsol();

   //assert( !vscb_.empty() );
   assert( !lvscb_.empty() );
   //assert( (int) vscb_.size() == dim );
   assert( (int) lvscb_.size() == dim );
   assert( !rsol.empty() );
   assert( (int) rsol.size() >= dim );

   double rsol_i;
   double M_i;
   double max_M = 0.0;

   #if debug
   cout << "[ index:" << (int) node.get_index() << "]-------------------------------" << endl;
   #endif

   //for ( auto i = dim - 1; i > dim - 6; --i )
   //for ( auto i = 0; i < dim; ++i )
   for ( auto i = dim - 1; i && i > dim - 6; --i )
   {

      rsol_i = rsol[i];

      M_i = round( rsol_i );
      M_i -= rsol_i;
      M_i *= M_i;
      M_i *= lvscb_[i];

      if ( max_M < M_i )
         max_M = M_i;

      #if debug
      printf("d:%.8f  ", lvscb_[i]);
      cout << "rsol_i:" << rsol_i << "  M_i:" << M_i << "  (" << (int) node.get_lowerbound() << ")" << endl;
      #endif
   }

   #if debug || 0
   //cout << "max_M:" << max_M << endl;
   printf("f(bar_x):%d Max_M:%d bestval:%d\n", (int)node.get_lowerbound(), (int)max_M, (int)bestval);
   #endif

   max_M += node.get_lowerbound();

   #if debug
   cout << "f(rsol) + Max:" << max_M << endl;
   #endif

   if ( ceil( max_M ) < bestval )
      return R_FEASIBLE;
   else
   {
      #if debug
      cout << "max_M:infeasible!" << endl;
      #endif
      return R_INFEASIBLE;
   }
}

RelaxResult SVPsolver::SVPSsolveRelaxation(
      NODE&          node
   )
{

   assert(0);
   //if ( bestval - node.get_lowerbound() < 1.0 )
      return R_INFEASIBLE;

   //// copy
   //const auto& pd = probdata;
   //const auto m = pd.get_m();
   //const auto Q = pd.get_Q();
   //const auto ep = epsilon;
   //auto& qd = qpdata;

   //const auto& ub = node.get_ub();
   //const auto& lb = node.get_lb();
   //const auto& vec = node.get_sum_fixed();
   //auto& warm = node.geti_rsol();
   //const int dim = (int) ub.size();

   //vector<int> nofixed_index;
   //nofixed_index.reserve( dim );

   //assert( m > 0 );
   //assert( Q != nullptr );
   //assert( ep < 1.0 && ep > 0.0 );
   //assert( dim >= 1 );
   //assert( dim < m );
   //assert( !ub.empty() );
   //assert( !lb.empty() );
   //assert( !warm.empty() );
   //assert( dim == (int) ub.size() );
   //assert( dim == (int) lb.size() );
   //assert( dim == (int) warm.size() );
   //assert( m == (int) vec.size() );

   //for ( auto i = 0; i < dim; ++i )
   //{
   //   if ( ub[i] != lb[i] )
   //      nofixed_index.push_back( i );
   //}

   //nofixed_index.shrink_to_fit();

   //const int n = (int) nofixed_index.size();

   //// all values are fixed
   //if ( n == 0 )
   //{
   //   double* relaxvals = new double[m];

   //   for ( auto i = 0; i < dim; ++i )
   //      relaxvals[i] = ub[i];

   //   auto fixedvalues = node.get_fixedvalues().end();
   //   for ( auto i = dim; i < m; ++i )
   //   {
   //      --fixedvalues;
   //      relaxvals[i] = *fixedvalues;
   //   }

   //   const auto objval = compute_objval( relaxvals );
   //   SOLUTION solution;

   //   solution.set_sol( m, relaxvals, objval );
   //   SVPStrySol( solution, true, true, nullptr );

   //   delete[] relaxvals;

   //   return R_GETINTEGER;
   //}

   //double* subQ = nullptr;
   //double* subu = nullptr;
   //double* subl = nullptr;
   //double* subw = nullptr;
   //double* p = nullptr;
   //double c_term = 0.0;

   //QPsolver  qps;

   //int   ct = 0;

   //qps.set_dim( n );

   //if( n == m )
   //{
   //   assert(0);
   //   //qps.set_obj_quad( m, Q );
   //   //qps.set_lb( m, ub );
   //   //qps.set_ub( m, vars_localub );
   //   //qps.set_warm( m, &warm[0] );
   //}
   //else
   //{
   //   assert( m - 1 >= n );

   //   if( qd.check_allocation() == false )
   //      qd.alloc( m - 1 );

   //   subQ = qd.get_Qmat();
   //   assert( subQ != nullptr );

   //   ct = 0;

   //   for ( auto i : nofixed_index )
   //   {
   //      auto im = i * m;
   //      for ( auto j : nofixed_index )
   //         subQ[ct++] = Q[j + im];
   //   }

   //   assert( ct == n*n );
   //   assert( subQ[0+(1*n)] == subQ[1+(0*n)] );

   //   if( !vec.empty() )
   //   {
   //      assert( (int) vec.size() == m );

   //      p = qd.get_pvec();
   //      assert( p != nullptr );

   //      ct = 0;

   //      assert( !vec.empty() );
   //      assert( m == (int) vec.size() );
   //      for ( auto i : nofixed_index )
   //      {
   //         p[ct] = 0;
   //         auto bvec_i = pd.get_bvec(i);
   //         for ( auto j = 0; j < m; ++j )
   //            p[ct] += 2 * bvec_i[j] * vec[j];
   //         //p[ct++] = 2 * Com_dot( probdata.get_bvec(i), &vec[0], m );
   //         ++ct;
   //      }

   //      assert( ct == n );

   //      qps.set_obj( n, subQ, p);

   //      //c_term = Com_dot( vec, vec, m);
   //      c_term = 0;
   //      for ( auto val : vec )
   //         c_term += val * val;
   //   }
   //   else
   //   {
   //      qps.set_obj_quad( n, subQ );
   //   }

   //   subu = qd.get_uvec();
   //   subl = qd.get_lvec();
   //   subw = qd.get_wvec();

   //   assert( subu != nullptr );
   //   assert( subl != nullptr );
   //   assert( subw != nullptr );

   //   ct = 0;

   //   for ( auto i : nofixed_index )
   //   {
   //      subu[ct] = ub[i];
   //      subl[ct] = lb[i];
   //      ct++;
   //   }

   //   assert( ct == n );

   //   ct = 0;

   //   for ( auto i : nofixed_index )
   //   {
   //      subw[ct] = warm[i];
   //      ct++;
   //   }

   //   qps.set_ub( n, subu );
   //   qps.set_lb( n, subl );
   //   qps.set_warm( n, subw );
   //}

   ////qps.disp_prob();

   //qps.set_ep( 1e-12 );

   //qps.solve( nullptr );
   ////qps.solve( &testwatch );

   //auto qpbestvals = qps.get_bestsol();
   //ct = 0;

   //for ( int i = 0; i < dim; ++i )
   //{
   //   if( ub[i] != lb[i] )
   //      warm[i] = qpbestvals[ct++];
   //   else
   //      warm[i] = ub[i];
   //}

   //assert( ct == n );

   //double qpbestval = qps.get_bestval() + c_term;

   //#if debug
   //cout << c_term << endl;
   //cout << "relax:" << qpbestval << endl;
   //#endif

   //if( node.get_lowerbound() < qpbestval )
   //   node.set_lowerbound( qpbestval );

   //RelaxResult result;

   //if( bestval - qpbestval  < 1.0 )
   //{
   //   result = R_INFEASIBLE;
   //}
   //else
   //{
   //   result = R_FEASIBLE;

   //   vector<double> roundvals( m );
   //   for( int i = 0; i < dim; i++ )
   //      roundvals[i] = round( warm[i] );

   //   auto fixedvalues = node.get_fixedvalues().end();
   //   for ( auto i = dim; i < m; ++i )
   //   {
   //      --fixedvalues;
   //      roundvals[i] = *fixedvalues;
   //   }

   //   auto objval = compute_objval( &roundvals[0] );

   //   if( Equal( objval, qpbestval, ep ) )
   //      result = R_GETINTEGER;

   //   SOLUTION solution;
   //   solution.set_sol( m, &roundvals[0], objval );
   //   SVPStrySol( solution, true, true, nullptr );

   //   // Do not use relaxvals after here
   //}

   //return result;
}

