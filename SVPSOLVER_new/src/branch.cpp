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
#include <utility>

#include "svpsolver.h"
#include "node.h"
#include "probdata.h"
#include "qpsolver.h"
#include "solution.h"
#include "vector.h"

using namespace std;

#define debug  0
#define debug_s  0

//static int find_branchingvariable(
//      const int      m,
//      const double*  relax_solval
//   )
//{
//   double buf;
//   for ( int i = m - 1; 0 <= i; --i )
//   {
//      buf = relax_solval[i];
//      if( !Equal( buf, round( buf ), 1.0e-12 ) )
//         return i;
//   }
//
//   return -1;
//}

BranchResult  SVPsolver::SVPSbranchDefault(
      NODE&       node
   )
{
   assert( !ENUM );

   // copy member variables of SVPsolver
   const auto m = probdata.get_m();
   auto& nodeindex = index;
   const auto& currentbestval = bestval;
   auto& NL = nodelist;

   assert( m > 0 );
   #if debug
   cout << "[Fix";
   cout << " best:" << (int)currentbestval;
   cout << " index:" << (int) node.get_index() << "]-------------------------------" << endl;
   #endif
   // copy member variables of node
   const auto dpt = node.get_dpt();
   const auto dpt1 = dpt + 1;
   const auto& fixedvalues = node.cget_fixedvalues();
   const auto lowerbound = node.get_lowerbound();

   assert( ( fixedvalues.empty() && dpt == 0 ) || ( !fixedvalues.empty() && dpt == (int)fixedvalues.size() ) );

   #if debug
   cout << "fixedvalues: ";
   for ( auto val : fixedvalues )
      cout << val << " ";
   cout << endl;
   #endif

   const auto branchindex = m - 1 - dpt;

   assert( branchindex >= 0 );
   assert( branchindex < m );

   auto& rsol = node.geti_rsol();
   const auto lastvalue = rsol[branchindex];

   assert( branchindex >= 0 && branchindex < m );
   assert( !rsol.empty() );
   assert( m - dpt <= (int) rsol.size() );


   if ( dpt == m - 1 )
   {
      double fixvalue = round( lastvalue );
      double upperbound_new = fixvalue;
      upperbound_new -= lastvalue;
      upperbound_new *= upperbound_new;
      upperbound_new *= *(probdata.get_Q());
      upperbound_new += lowerbound;

      #if debug
      printf("x_%d = %d -> upper: %d (%d)\n", branchindex, (int)fixvalue, (int) upperbound_new, nodeindex);
      #endif
      if ( upperbound_new < currentbestval && (int) upperbound_new )
      {
         SOLUTION setsol;
         double* solvals = new double[m];
         solvals[0] = fixvalue;
         assert( (int) fixedvalues.size() == dpt );
         for ( auto j = 1, k = dpt - 1; j < m; ++j, --k )
            solvals[j] = fixedvalues[k];

         setsol.set_sol( m, solvals, compute_objval(solvals) );
         SVPStrySol( setsol, true, true, nullptr );
         delete[] solvals;
         //assert(0);
      }
      return B_GETINTEGER;
   }

   assert( dpt < m - 1 );

   const auto ceilval = (int) ceil( lastvalue );
   const auto floorval = (int) floor( lastvalue );
   const auto& vsco_i = VsCO[branchindex -1];
   const auto vcov_i = VCOV[branchindex -1];
   const auto m_dpt1 = m - dpt1;
   double init_lowerbound_new;
   int init_fixvalue;
   double sign = 1.0;

   assert( !vsco_i.empty() );
   assert( m - dpt1 == (int) vsco_i.size() );

   if ( ceilval - lastvalue < lastvalue - floorval )
   {
      sign = - 1.0;
      init_fixvalue = ceilval;
   }
   else
      init_fixvalue = floorval;

   assert( init_fixvalue == round( lastvalue ) );

   double init_coef = init_fixvalue;
   init_coef -= lastvalue;

   init_lowerbound_new = init_coef;
   init_lowerbound_new *= init_coef;
   init_lowerbound_new *= vcov_i;
   init_lowerbound_new += lowerbound;

   #if debug
   printf("x_%d = %d -> lower: %f (-)\n", branchindex, (int) init_fixvalue , init_lowerbound_new );
   #endif

   if ( init_lowerbound_new < currentbestval )
   {
      double lowerbound_new;
      double coef;
      auto counter = 2;
      auto fixvalue = init_fixvalue;
      fixvalue += sign;

      function<void(void)> update_fix;
      if ( lowerbound < 1.0e-10 )
         update_fix = [&fixvalue]{ ++fixvalue; };
      else
       update_fix = [&fixvalue, &sign, &counter]{ sign = -sign; fixvalue += sign * counter; ++counter; };

      while ( 1 )
      {
         coef = fixvalue;
         coef -= lastvalue;

         lowerbound_new = coef;
         lowerbound_new *= coef;
         lowerbound_new *= vcov_i;
         lowerbound_new += lowerbound;

         #if debug
         printf("x_%d = %d -> lower: %f (%d)\n", branchindex, (int) fixvalue , lowerbound_new, nodeindex );
         #endif
         if ( lowerbound_new < currentbestval )
         {
            auto& nn = NL.emplace_back_TENDEQUE( dpt1, nodeindex, lowerbound_new, fixedvalues );
            nn.push_back_fixedvalue( fixvalue );

            auto& new_rsol = nn.geti_rsol();
            //new_rsol.assign( &rsol[0], &rsol[m_dpt1] );
            new_rsol.resize( m_dpt1 );

            for ( auto j = 0; j < m_dpt1 ; ++j )
            {
               new_rsol[j] = rsol[j];
               new_rsol[j] -= coef * vsco_i[j];
            }

            ++nodeindex;
         }
         else
            break;
         update_fix();
      }

      node.inc_dpt();
      node.set_lowerbound( init_lowerbound_new );
      node.push_back_fixedvalue( init_fixvalue );

      for ( auto j = 0; j < m_dpt1 ; ++j )
         rsol[j] -= init_coef * vsco_i[j];

      ++nodeindex;
   }
   else
      return B_INFEASIBLE;

   #if debug
   //cout << NL.getListsize() << endl;
   #endif
   assert( !fixedvalues.empty() && dpt+1 == (int)fixedvalues.size() );

   return B_FEASIBLE;
}

BranchResult  SVPsolver::SVPSbranchEnum(
      NODE&       node
   )
{
   assert( ENUM );
   assert( !IMPRELAX );

   // copy member variables of SVPsolver
   const auto m = probdata.get_m();
   const auto& SN = SNOVs;
   const auto& CM = CMGSO;
   const auto ep = epsilon;
   auto& nodeindex = index;
   const auto currentbestval = bestval;
   auto& NL = nodelist;

   assert( m > 0 );
   assert( (int) SN.size() == m );
   assert( (int) CM.size() == m );
   #if debug
   cout << "[Fix";
   cout << " best:" << (int)currentbestval;
   cout << " index:" << (int) node.get_index() << "]-------------------------------" << endl;
   #endif
   // copy member variables of node
   const auto dpt = node.get_dpt();
   const auto dpt1 = dpt + 1;
   const auto& fixedvalues = node.get_fixedvalues();
   const auto& lowerbound = node.get_lowerbound();

   assert( ( fixedvalues.empty() && dpt == 0 )
         || ( !fixedvalues.empty() && dpt == (int)fixedvalues.size() ) );

   #if debug
   cout << "fixedvalues: ";
   for ( auto val : fixedvalues )
      cout << val << " ";
   cout << endl;
   #endif


   const auto branchindex = m - 1 - dpt;
   double sum = 0.0;
   double buf;
   int upper;
   int lower;
   vector<NODE> newnodes;

   assert( branchindex >= 0 );
   assert( branchindex < m );

   // compute bounds on branching variable
   if ( dpt > 0 )
   {
      const auto& vec = CM[dpt];
      assert( (int) vec.size() == dpt + 1 );
      for ( auto i = 0; i < dpt; ++i )
      {
         //cout << vec[i] << " " << fixedvalues[i] << endl;
         sum += vec[i] * fixedvalues[i];
      }
   }

   buf = bestval - lowerbound;
   buf /= SN[branchindex];
   buf = sqrt( buf );

   upper = floor( buf - sum + ep );

   if ( lowerbound < ep )
   {
      lower = 0;
   }
   else
      lower = ceil( - buf - sum - ep );

   #if debug
   printf("fix x_%d : [ %d, %d ] \n", branchindex, lower, upper );
   #endif

   if ( lower > upper )
   {
      #if debug
      printf("x_%d: INFEASIBLE\n", branchindex );
      #endif
      return B_INFEASIBLE;
   }

   if ( dpt == m - 1 )
   {
      double upperbound_new;
      for ( auto i = lower; i <= upper; ++i )
      {
         upperbound_new = sum;
         upperbound_new += (double) i;
         upperbound_new *= upperbound_new;
         upperbound_new *= SN[branchindex];
         upperbound_new += lowerbound;
         upperbound_new += ep;
         upperbound_new = floor( upperbound_new );
         #if debug
         printf("x_%d = %d -> upper: %d (%d)\n", branchindex, i, (int) upperbound_new,
               nodeindex);
         #endif
         if ( upperbound_new < currentbestval && (int) upperbound_new )
         {
            SOLUTION setsol;
            double* solvals = new double[m];
            solvals[0] = i;
            for ( auto j = 1; j < m; ++j )
            {
               solvals[j] = fixedvalues[m-1-j];
            }

         setsol.set_sol( m, solvals, compute_objval(solvals) );
            SVPStrySol( setsol, true, true, nullptr );
            delete[] solvals;
            //assert(0);
         }
      }

      return B_GETINTEGER;
   }

   assert( dpt < m - 1 );

   double lowerbound_new;
   double min_l = currentbestval;
   int found = - 1;
   //const auto m_dpt1 = m - dpt1;
   //double coef;

   newnodes.reserve( upper - lower + 1 );

   //cout << "b:rsol: ";
   //for ( auto r : rsol )
   //   cout << r << " ";
   //cout << endl;
   //printf( "dpt:%d dpt1:%d m_dpt1:%d rsol[m_dpt1]:%f", dpt, dpt1, m_dpt1, rsol[m_dpt1] );
   //cout << endl;


   for ( int i = lower; i <= upper; ++i )
   {
      lowerbound_new = sum;
      lowerbound_new += (double) i;
      lowerbound_new *= lowerbound_new;
      lowerbound_new *= SN[branchindex];
      lowerbound_new += lowerbound;
      #if debug
      printf("x_%d = %d -> lower: %f (%d)\n", branchindex, i, lowerbound_new ,
            nodeindex);
      #endif
      if ( lowerbound_new < currentbestval )
      {
         newnodes.emplace_back( dpt1, nodeindex, lowerbound_new, fixedvalues );
         auto& nn = newnodes.back();
         nn.get_fixedvalues().push_back( i );

         //if ( improvedrelaxation )
         //{
         //   auto& new_rsol = nn.geti_rsol();
         //   new_rsol.assign( &rsol[0], &rsol[m_dpt1] );

         //   coef = i - rsol[m_dpt1];
         //   for ( auto j = 0; j < m_dpt1 ; ++j )
         //      new_rsol[j] -= coef * vsco_i[j];
         //}

         if ( min_l > lowerbound_new )
         {
            found = nodeindex;
            min_l = lowerbound_new;
         }
         ++nodeindex;
      }
   }

   #if debug
   if ( found == - 1 )
      cout << "set: " << "no" << endl;
   else
      cout << "set: " << found << endl;
   #endif

   if ( found == - 1 )
      return B_INFEASIBLE;

   for ( auto& nn : newnodes )
   {
      if ( found != nn.get_index() )
         (NL.*move_back)( nn );
      else
         node = move( nn );
   }

   #if debug
   //cout << NL.getListsize() << endl;
   #endif
   assert( ( fixedvalues.empty() && dpt == 0 )
         || ( !fixedvalues.empty() && dpt+1 == (int)fixedvalues.size() ) );

   return B_FEASIBLE;
}

BranchResult  SVPsolver::SVPSbranchStandard(
      NODE&       node
   )
{
   assert(0);
   return B_FEASIBLE;
   //if ( node.get_presolved() == false )
   //{
   //   #if debug_s
   //   cout << "[presolving.. index:" << (int)node.get_index();
   //   cout << "]-----------------------------" << endl;
   //   #endif
   //   if ( !SVPSpresolveNode( node ) )
   //      return B_INFEASIBLE;
   //}

   //// copy
   //const auto& pd = probdata;
   //const auto m = pd.get_m();
   //auto& nodeindex = index;

   //assert( m > 0 );
   //assert( nodeindex >= 0 );

   //// node
   //NODE C_RIGHT = node;
   //auto& relax_solval_R = C_RIGHT.geti_rsol();
   //const auto branchindex = find_branchingvariable( relax_solval_R.size(), &relax_solval_R[0] );
   //const int ceil_value = ceil( relax_solval_R[branchindex] );
   //const auto set_dpt = node.get_dpt() + 1;

   //assert( !relax_solval_R.empty() );
   //assert( !C_RIGHT.get_sum_fixed().empty() );
   //assert( branchindex >= 0 );
   //assert( branchindex < (int)C_RIGHT.get_sum_fixed().size() );
   //assert( set_dpt > 0 );

   //// ceil <= x_i
   //C_RIGHT.set_lb_i( branchindex, ceil_value );
   //C_RIGHT.set_dpt( set_dpt );
   //C_RIGHT.set_index( nodeindex );
   //++nodeindex;

   //if( ceil_value && ceil_value == C_RIGHT.get_ub()[branchindex] )
   //{
   //   assert( !C_RIGHT.get_sum_fixed().empty() );
   //   assert( (int)C_RIGHT.get_sum_fixed().size() == m );
   //   auto& sumfixed = C_RIGHT.geti_sum_fixed();
   //   const auto bvec = pd.get_bvec( branchindex );
   //   for ( auto i = 0; i < m; ++i )
   //      sumfixed[i] += ceil_value * bvec[i];
   //}

   //relax_solval_R[branchindex] = ceil_value;

   //(nodelist.*move_back)( C_RIGHT );
   //// Do not use C_RIGHT after here

   //// x_i <= floor
   //auto& relax_solval_L = node.geti_rsol();
   //const int floor_value = floor( relax_solval_L[branchindex] );
   //node.set_ub_i( branchindex, floor_value );
   //node.set_dpt( set_dpt );
   //node.set_index( nodeindex );
   //++nodeindex;

   //assert( !relax_solval_L.empty() );

   //if ( floor_value && floor_value == node.get_lb()[branchindex] )
   //{
   //   assert( !node.get_sum_fixed().empty() );
   //   assert( (int)node.get_sum_fixed().size() == m );
   //   auto& sumfixed = node.geti_sum_fixed();
   //   const auto bvec = pd.get_bvec( branchindex );
   //   for ( auto i = 0; i < m; ++i )
   //      sumfixed[i] += floor_value * bvec[i];
   //}

   //relax_solval_L[branchindex] = floor_value;

   //#if debug_s
   //printf("x_%d <= %d and %d <= x_%d\n", branchindex, floor_value, ceil_value, branchindex );
   //#endif

   //return B_FEASIBLE;
}

bool  SVPsolver::SVPSpresolveNode(
      NODE&       node
   )
{
   assert(0);
   return 0;
   //assert( node.get_lowerbound() > 0.0 );
   //assert( !node.get_fixedvalues().empty() );

   //// copy
   //const auto& pd = probdata;
   //const auto m = pd.get_m();
   //const auto B_ = pd.get_B_();

   //const auto sq_best = sq_bestval;

   //assert( m > 0 );
   //assert( Equal( sq_best, sqrt( bestval ), epsilon ) );

   //const auto& fixedvalues = node.get_fixedvalues();
   //const auto numfixed = (int) fixedvalues.size();
   //const auto subdim = m - numfixed;
   //const auto& vscb_ = VsCB[subdim - 1];
   //const auto& lvscb_ = LVsCB[subdim - 1];

   //assert( 0 );
   //// origina lvscb_ is norm, but now it is norm^2

   //assert( numfixed > 0 && numfixed < m );
   //assert( !vscb_.empty() );
   //assert( !lvscb_.empty() );
   //assert( (int) vscb_.size() == subdim );
   //assert( (int) lvscb_.size() == subdim );

   //auto& sumfixed = node.geti_sum_fixed();
   //auto& upper = node.geti_ub();
   //auto& lower = node.geti_lb();
   //auto& rsol = node.geti_rsol();
   //int origindex;
   //int om;
   //double buf_1;
   //double buf_2;

   //assert( sumfixed.empty() );
   //assert( upper.empty() );
   //assert( lower.empty() );
   //assert( rsol.empty() );

   //assert(0);
   //exit(-1);
   //// old sumfixed is empty, but now it is not empty

   //#if debug_s
   //cout << "fixedvalues: ";
   //for ( auto val : node.get_fixedvalues() )
   //   cout << val << " ";
   //cout << endl;
   //#endif

   //auto it = fixedvalues.end();
   //origindex = subdim;

   //sumfixed.assign( m, 0 );
   //upper.resize( subdim );
   //lower.resize( subdim );
   //rsol.resize( subdim );

   //// compute sumfixed
   //do
   //{
   //   --it;

   //   if ( *it )
   //   {
   //      om = origindex * m;
   //      for ( auto i = 0; i < m; ++i )
   //         sumfixed[i] += (*it) * B_[om + i];
   //   }

   //   ++origindex;
   //} while ( it != fixedvalues.begin() );


   //// compute ub, lb and rsol
   //for ( auto i = 0; i < subdim; ++i )
   //{
   //   buf_1 = sq_best;
   //   buf_1 *= lvscb_[i];

   //   assert( !vscb_[i].empty() );
   //   assert( (int) vscb_[i].size() == m );

   //   buf_2 = 0;
   //   for ( auto j = 0; j < m; ++j )
   //      buf_2 += vscb_[i][j] * sumfixed[j];

   //   upper[i] = (int) floor( buf_1 - buf_2 );
   //   lower[i] = (int) ceil( - buf_1 - buf_2 );

   //   if ( lower[i] == upper[i] && lower[i] )
   //   {
   //      om = i * m;
   //      for ( auto j = 0; j < m; ++j )
   //         sumfixed[j] += lower[i] * B_[om + j];
   //   }
   //   else if ( lower[i] > upper[i] )
   //      return false;

   //   rsol[i] = (double) ( lower[i] + upper[i] ) * 0.5;
   //   rsol[i] += 0.1;
   //}

   //#if debug_s
   //for ( auto i = 0; i < subdim; ++i )
   //   printf("%d <= x_%d <= %d \n ", lower[i], i, upper[i]);
   //#endif

   //node.set_presolved( true );

   //return true;
}

void  SVPsolver::SVPSbranch_INT(
      NODE&       node,
      const int   index,
      double*     vars_localub,
      double*     vars_locallb
   )
{
   assert(0);
   //assert( index >= 0 );
   //assert( vars_localub != nullptr );
   //assert( vars_locallb != nullptr );

   //// copy
   //const auto m = probdata.get_m();
   //auto& nodeindex = index;
   //auto& pd = probdata;
   //auto ep = epsilon;

   //assert( m > 0 );
   //assert( nodeindex >= 0 );

   //// node
   //NODE C_RIGHT = node;
   //auto relax_solval_R = C_RIGHT.get_relaxsolval();
   //const auto branchindex = find_branchingvariable( m, relax_solval_R );
   //const int ceil_value = ceil( relax_solval_R[branchindex] );
   //const auto set_dpt = node.get_dpt() + 1;

   //assert( relax_solval_R != nullptr );
   //assert( branchindex >= 0 );
   //assert( branchindex < m );
   //assert( set_dpt > 0 );

   //// ceil <= x_i
   //C_RIGHT.set_branchinfo( branchindex, ceil_value, 'l' );
   //C_RIGHT.set_dpt( set_dpt );
   //C_RIGHT.set_index( nodeindex );

   //if( ceil_value && Equal( (double) ceil_value, vars_localub[branchindex], ep ) )
   //{
   //   if( C_RIGHT.alloc_sumfixed() )
   //      C_RIGHT.set_sumfixed( ceil_value, pd.get_bvec( branchindex ) );
   //   else
   //      C_RIGHT.add_sumfixed( ceil_value, pd.get_bvec( branchindex ) );
   //}

   //relax_solval_R[branchindex] = ceil_value;

   //(nodelist.*move_back)( C_RIGHT );
   //// Do not use C_RIGHT after here

   //// x_i <= floor
   //auto relax_solval_L = node.get_relaxsolval();
   //const int floor_value = floor( relax_solval_L[branchindex] );
   //node.set_branchinfo( branchindex, floor_value, 'u' );
   //node.set_dpt( set_dpt );
   //node.set_index( index + 1 );

   //assert( relax_solval_L != nullptr );

   //if ( floor_value && Equal( vars_locallb[branchindex], floor_value, ep ) )
   //{
   //   if( node.alloc_sumfixed() == true )
   //      node.set_sumfixed( floor_value, pd.get_bvec( branchindex ) );
   //   else
   //      node.add_sumfixed( floor_value, pd.get_bvec( branchindex ) );
   //}

   //relax_solval_L[branchindex] = floor_value;

   //vars_localub[branchindex] = floor_value;
}


